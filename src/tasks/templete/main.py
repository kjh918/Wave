"""
Workflow Task Class Scaffold
- Minimal but production-lean foundation for a registry-based workflow engine
- Includes: Task base class, registry, output validation, run/skip logic with flag files,
  an execution context, a simple DAG runner, and two example tasks (FastQC, PicardDownsampleSam)

You can drop this whole file into your project as `workflow_core.py` or split it
into modules (see section headers). Type hints and docstrings included.

Assumptions pulled from prior chats:
- Use a registry decorator: `@TaskRegistry.register`
- Create a flag file when *all* outputs exist, are non-empty
- Executor should check flag before running and skip completed tasks
- Commands are built with explicit parameters (not a dict of args)
- Singularity containers commonly used with `-B` binds
"""

from __future__ import annotations

import hashlib
import json
import os
import shlex
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, ClassVar, Dict, Iterable, List, Optional, Sequence, Tuple

# =============================================================
# Registry
# =============================================================

class TaskRegistry:
    _REGISTRY: ClassVar[Dict[str, type["Task"]]] = {}

    @classmethod
    def register(cls, klass: type["Task"]) -> type["Task"]:
        """Decorator to register a Task subclass by its TYPE attribute."""
        t = getattr(klass, "TYPE", None)
        if not isinstance(t, str) or not t:
            raise ValueError(f"Task class {klass.__name__} must define TYPE: str")
        if t in cls._REGISTRY:
            raise ValueError(f"Duplicate task TYPE '{t}' for {klass.__name__}")
        cls._REGISTRY[t] = klass
        return klass

    @classmethod
    def get(cls, type_name: str) -> type["Task"]:
        try:
            return cls._REGISTRY[type_name]
        except KeyError:
            raise KeyError(f"Unknown task TYPE '{type_name}'. Registered: {list(cls._REGISTRY)}")

    @classmethod
    def create(cls, type_name: str, **kwargs: Any) -> "Task":
        return cls.get(type_name)(**kwargs)


# =============================================================
# Flags & Output Validation
# =============================================================

@dataclass
class OutputSpec:
    path: Path
    must_exist: bool = True
    non_empty: bool = True

    def is_satisfied(self) -> bool:
        if self.must_exist and not self.path.exists():
            return False
        if self.non_empty and self.path.exists() and self.path.is_file():
            try:
                return self.path.stat().st_size > 0
            except FileNotFoundError:
                return False
        return True


@dataclass
class FlagFiles:
    base: Path

    @property
    def done(self) -> Path:
        return self.base.with_suffix(self.base.suffix + ".done") if self.base.suffix else self.base.with_name(self.base.name + ".done")

    @property
    def fail(self) -> Path:
        return self.base.with_suffix(self.base.suffix + ".fail") if self.base.suffix else self.base.with_name(self.base.name + ".fail")

    @property
    def meta(self) -> Path:
        return self.base.with_suffix(self.base.suffix + ".meta.json") if self.base.suffix else self.base.with_name(self.base.name + ".meta.json")

    def mark_done(self, meta: Optional[Dict[str, Any]] = None) -> None:
        self.fail.unlink(missing_ok=True)
        self.done.parent.mkdir(parents=True, exist_ok=True)
        self.done.write_text("ok\n")
        if meta:
            self.meta.write_text(json.dumps(meta, indent=2, ensure_ascii=False))

    def mark_fail(self, meta: Optional[Dict[str, Any]] = None) -> None:
        self.done.unlink(missing_ok=True)
        self.fail.parent.mkdir(parents=True, exist_ok=True)
        self.fail.write_text("fail\n")
        if meta:
            self.meta.write_text(json.dumps(meta, indent=2, ensure_ascii=False))

    def is_done(self) -> bool:
        return self.done.exists()

    def is_failed(self) -> bool:
        return self.fail.exists()


# =============================================================
# Task Base + Decorator that enforces output validation & flags
# =============================================================

@dataclass
class Task:
    """Base class for workflow tasks.

    Subclasses must define:
      - TYPE: str
      - outputs(): list[OutputSpec]
      - run(): the side-effectful execution (write outputs)

    Optionally override:
      - inputs_signature(): produce a stable hash to scope flags per inputs
      - build_cmd(): return a shell command (for logging/testing)
    """

    TYPE: ClassVar[str]

    work_dir: Path
    name: str
    binds: Sequence[Path] = field(default_factory=lambda: [Path("/storage"), Path("/data")])
    image: Optional[Path] = None  # singularity image
    threads: int = 2

    def outputs(self) -> List[OutputSpec]:
        raise NotImplementedError

    def run(self) -> None:
        """Do the actual work. Should raise on failure."""
        raise NotImplementedError

    # ---- helpers ----
    def inputs_signature(self) -> str:
        """Stable signature used to scope flag paths (e.g., based on primary inputs).
        Default = name. Override if multiple input combinations share a name.
        """
        return hashlib.sha1(self.name.encode()).hexdigest()[:12]

    def flag_files(self) -> FlagFiles:
        flag_base = self.work_dir / f".flag_{self.TYPE}_{self.inputs_signature()}"
        return FlagFiles(flag_base)

    def check_outputs(self) -> bool:
        return all(spec.is_satisfied() for spec in self.outputs())

    # Optional: string command representation for logging/testing
    def build_cmd(self) -> Optional[str]:  # pragma: no cover – override in subclasses
        return None


# Decorator to auto-validate outputs and write flags

def finalize_and_flag(run_method: Callable[[Task], None]) -> Callable[[Task], None]:
    def wrapper(self: Task) -> None:
        flags = self.flag_files()
        try:
            run_method(self)
            if not self.check_outputs():
                raise RuntimeError("Output validation failed: one or more outputs missing or empty")
            flags.mark_done(meta={
                "task": self.TYPE,
                "name": self.name,
                "cmd": self.build_cmd(),
                "threads": self.threads,
            })
        except Exception as e:
            flags.mark_fail(meta={
                "task": self.TYPE,
                "name": self.name,
                "cmd": self.build_cmd(),
                "error": str(e),
            })
            raise
    return wrapper


# =============================================================
# Simple Executor + DAG
# =============================================================

@dataclass
class TaskNode:
    task: Task
    deps: List["TaskNode"] = field(default_factory=list)

    def add_dep(self, other: "TaskNode") -> "TaskNode":
        self.deps.append(other)
        return self


class Executor:
    """Runs TaskNodes in dependency order. Skips tasks with .done flags and valid outputs.
    Fails fast if a dependency fails.
    """

    def __init__(self, dry_run: bool = False):
        self.dry_run = dry_run

    def _run_task(self, node: TaskNode) -> None:
        # deps first
        for d in node.deps:
            self._run_task(d)

        task = node.task
        flags = task.flag_files()

        # Skip if already done and outputs still valid
        if flags.is_done() and task.check_outputs():
            print(f"[SKIP] {task.TYPE}:{task.name} – already completed")
            return

        # If previously failed, we still attempt rerun (fresh)
        if flags.is_failed():
            print(f"[RETRY] {task.TYPE}:{task.name} – previous failure flag found, retrying…")

        cmd = task.build_cmd()
        if self.dry_run:
            print(f"[DRY-RUN] {task.TYPE}:{task.name}\n  CMD: {cmd}")
            return

        print(f"[RUN ] {task.TYPE}:{task.name}")
        if cmd:
            print(f"  CMD: {cmd}")

        task.run()  # decorator will validate & flag
        print(f"[DONE] {task.TYPE}:{task.name}")

    def run(self, roots: Iterable[TaskNode]) -> None:
        seen: set[int] = set()

        def dfs(n: TaskNode) -> None:
            if id(n) in seen:
                return
            seen.add(id(n))
            self._run_task(n)

        for r in roots:
            dfs(r)


# =============================================================
# Shell helpers
# =============================================================

def run_shell(cmd: str, env: Optional[Dict[str, str]] = None) -> None:
    """Run a shell command and raise on non-zero exit."""
    print(f"[SHELL] {cmd}")
    proc = subprocess.run(cmd, shell=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {proc.returncode}")


def singularity_exec(image: Path, binds: Sequence[Path], inner_cmd: str) -> str:
    binds_arg = ",".join(str(b) for b in binds)
    parts = [
        "singularity", "exec",
        "-B", shlex.quote(binds_arg),
        shlex.quote(str(image)),
        inner_cmd,
    ]
    return " ".join(parts)


# =============================================================
# Example Task: FastQC
# =============================================================

@TaskRegistry.register
class FastQCRunner(Task):
    TYPE = "fastqc"

    fastqc_bin: str = "fastqc"
    extract: bool = True

    # inputs
    fastq_r1: Path | None = None
    fastq_r2: Path | None = None

    # outputs
    out_dir: Path | None = None

    def outputs(self) -> List[OutputSpec]:
        assert self.out_dir is not None, "out_dir must be set"
        outs = []
        if self.fastq_r1 is not None:
            outs.append(OutputSpec(self.out_dir / (self.fastq_r1.stem + "_fastqc.zip")))
        if self.fastq_r2 is not None:
            outs.append(OutputSpec(self.out_dir / (self.fastq_r2.stem + "_fastqc.zip")))
        return outs

    def inputs_signature(self) -> str:
        parts = [p for p in [self.fastq_r1, self.fastq_r2] if p]
        key = "|".join(str(p) for p in parts)
        return hashlib.sha1(key.encode()).hexdigest()[:12]

    def build_cmd(self) -> Optional[str]:
        assert self.out_dir is not None, "out_dir must be set"
        inputs = " ".join(shlex.quote(str(p)) for p in [self.fastq_r1, self.fastq_r2] if p)
        flags = "--extract" if self.extract else ""
        inner = f"{self.fastqc_bin} --threads {self.threads} {flags} -o {shlex.quote(str(self.out_dir))} {inputs}"
        if self.image:
            return singularity_exec(self.image, self.binds, inner)
        return inner

    @finalize_and_flag
    def run(self) -> None:
        self.out_dir.mkdir(parents=True, exist_ok=True)
        cmd = self.build_cmd()
        assert cmd is not None
        run_shell(cmd)


# =============================================================
# Example Task: Picard DownsampleSam (explicit args & defaults)
# =============================================================

@TaskRegistry.register
class PicardDownsampleSam(Task):
    TYPE = "picard_downsample"

    picard_jar: Path | None = None  # if using java -jar
    picard_bin: str = "picard"      # or containerized wrapper

    # inputs/outputs
    in_bam: Path | None = None
    out_bam: Path | None = None

    # explicit params with defaults (no dict)
    probability: float = 0.5
    seed: int = 42
    create_index: bool = True
    validation_stringency: str = "LENIENT"  # STRICT|LENIENT|SILENT

    def outputs(self) -> List[OutputSpec]:
        assert self.out_bam is not None, "out_bam must be set"
        outs = [OutputSpec(self.out_bam)]
        if self.create_index:
            outs.append(OutputSpec(self.out_bam.with_suffix(self.out_bam.suffix + ".bai")))
        return outs

    def inputs_signature(self) -> str:
        key = f"{self.in_bam}|{self.probability}|{self.seed}"
        return hashlib.sha1(key.encode()).hexdigest()[:12]

    def _inner_cmd(self) -> str:
        assert self.in_bam and self.out_bam
        args = [
            f"I={shlex.quote(str(self.in_bam))}",
            f"O={shlex.quote(str(self.out_bam))}",
            f"P={self.probability}",
            f"RANDOM_SEED={self.seed}",
            f"VALIDATION_STRINGENCY={self.validation_stringency}",
        ]
        if self.create_index:
            args.append("CREATE_INDEX=true")
        if self.picard_jar:
            return f"java -Xmx4g -jar {shlex.quote(str(self.picard_jar))} DownsampleSam " + " ".join(args)
        else:
            return f"{self.picard_bin} DownsampleSam " + " ".join(args)

    def build_cmd(self) -> Optional[str]:
        inner = self._inner_cmd()
        if self.image:
            return singularity_exec(self.image, self.binds, inner)
        return inner

    @finalize_and_flag
    def run(self) -> None:
        # ensure output dir exists
        assert self.out_bam is not None
        self.out_bam.parent.mkdir(parents=True, exist_ok=True)
        run_shell(self.build_cmd() or "")


# =============================================================
# Example: Wiring a small workflow
# =============================================================

if __name__ == "__main__":
    # Example usage (edit paths to your environment)
    wd = Path("/tmp/work_example")

    fastqc = FastQCRunner(
        work_dir=wd,
        name="sampleA_fastqc",
        threads=8,
        binds=[Path("/storage"), Path("/data")],
        image=Path("/storage/images/fastqc-0.12.1.sif"),
        fastq_r1=Path("/data/reads/sampleA_R1.fastq.gz"),
        fastq_r2=Path("/data/reads/sampleA_R2.fastq.gz"),
        out_dir=wd / "qc",
    )

    downsample = PicardDownsampleSam(
        work_dir=wd,
        name="sampleA_downsample",
        threads=4,
        image=Path("/storage/images/picard-3.2.0.sif"),
        in_bam=Path("/data/bam/sampleA.bam"),
        out_bam=wd / "downsample" / "sampleA.ds0.5.bam",
        probability=0.5,
        seed=777,
        create_index=True,
    )

    # Build a trivial DAG: fastqc has no deps, downsample depends on nothing here
    nodes = [TaskNode(fastqc), TaskNode(downsample)]

    Executor(dry_run=True).run(nodes)  # set dry_run=False to actually execute


# ================================
# src/tasks/task.py (base classes)
# ================================
from __future__ import annotations
import abc
import json
import hashlib
import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any, Iterable, List, Optional

class TaskRegistry:
    _reg: dict[str, type] = {}

    @classmethod
    def register(cls, t: type):
        name = getattr(t, "TYPE", None)
        if not name:
            raise ValueError("Task class must define TYPE")
        if name in cls._reg:
            raise ValueError(f"Duplicate task TYPE: {name}")
        cls._reg[name] = t
        return t

    @classmethod
    def get(cls, name: str):
        if name not in cls._reg:
            raise KeyError(f"Unknown task TYPE: {name}")
        return cls._reg[name]

    @classmethod
    def create(cls, name: str, **kwargs: Any):
        return cls.get(name)(**kwargs)


@dataclass
class FlagFiles:
    base: Path

    @property
    def done(self) -> Path:
        return self.base.with_suffix(self.base.suffix + ".done") if self.base.suffix else self.base.with_name(self.base.name + ".done")

    @property
    def fail(self) -> Path:
        return self.base.with_suffix(self.base.suffix + ".fail") if self.base.suffix else self.base.with_name(self.base.name + ".fail")

    @property
    def meta(self) -> Path:
        return self.base.with_suffix(self.base.suffix + ".meta.json") if self.base.suffix else self.base.with_name(self.base.name + ".meta.json")

    def mark_done(self, meta: Optional[Dict[str, Any]] = None) -> None:
        self.fail.unlink(missing_ok=True)
        self.done.parent.mkdir(parents=True, exist_ok=True)
        self.done.write_text("ok
")
        if meta:
            self.meta.write_text(json.dumps(meta, indent=2, ensure_ascii=False))

    def mark_fail(self, meta: Optional[Dict[str, Any]] = None) -> None:
        self.done.unlink(missing_ok=True)
        self.fail.parent.mkdir(parents=True, exist_ok=True)
        self.fail.write_text("fail
")
        if meta:
            self.meta.write_text(json.dumps(meta, indent=2, ensure_ascii=False))

    def is_done(self) -> bool:
        return self.done.exists()

    def is_failed(self) -> bool:
        return self.fail.exists()


def _normalize_output_paths(outputs: Dict[str, Any]) -> List[Path]:
    paths: List[Path] = []
    for v in outputs.values():
        if v is None:
            continue
        if isinstance(v, (str, Path)):
            paths.append(Path(v))
        elif isinstance(v, Iterable):
            for x in v:
                if x is None:
                    continue
                paths.append(Path(x))
        else:
            # ignore unknown types
            pass
    return paths


def _all_outputs_ok(paths: List[Path]) -> bool:
    for p in paths:
        if not p.exists():
            return False
        try:
            if p.is_file() and p.stat().st_size <= 0:
                return False
        except FileNotFoundError:
            return False
    return True


def run_shell(cmd: str) -> None:
    print(f"[SHELL] {cmd}")
    proc = subprocess.run(cmd, shell=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {proc.returncode}")


class Task(abc.ABC):
    TYPE: str
    DEFAULTS: Dict[str, Any] = {}

    def __init__(self,
                 name: str,
                 tool: str,
                 func: str,
                 threads: int,
                 workdir: Path,
                 inputs: Dict[str, Any] | None = None,
                 outputs: Dict[str, Any] | None = None,
                 params: Dict[str, Any] | None = None):
        self.name = name
        self.tool = tool
        self.func = func
        self.threads = threads  # fixed typo
        self.workdir = Path(workdir)

        self.inputs = inputs or {}
        self.outputs = outputs or {}
        self.params = {**(self.DEFAULTS or {}), **(params or {})}

        self.workdir.mkdir(parents=True, exist_ok=True)

    # ---- Flag helpers ----
    def _sig(self) -> str:
        # generate a stable signature from key inputs + params
        key = json.dumps({
            "name": self.name,
            "inputs": self.inputs,
            "params": self.params,
        }, sort_keys=True, ensure_ascii=False)
        return hashlib.sha1(key.encode()).hexdigest()[:12]

    def _flags(self) -> FlagFiles:
        return FlagFiles(self.workdir / f".flag_{self.TYPE}_{self._sig()}")

    def outputs_ok(self) -> bool:
        return _all_outputs_ok(_normalize_output_paths(self.outputs))

    # ---- Contract ----
    @abc.abstractmethod
    def to_sh(self) -> str:
        """Return a shell command (or a multi-line script) to execute this task."""

    # ---- Execution ----
    def run(self, dry_run: bool = False) -> None:
        flags = self._flags()

        if flags.is_done() and self.outputs_ok():
            print(f"[SKIP] {self.TYPE}:{self.name} already done")
            return

        cmd = self.to_sh()
        if dry_run:
            print(f"[DRY] {self.TYPE}:{self.name}
{cmd}")
            return

        try:
            run_shell(cmd)
            if not self.outputs_ok():
                raise RuntimeError("Output validation failed: missing/empty outputs")
            flags.mark_done(meta={
                "task": self.TYPE,
                "name": self.name,
                "tool": self.tool,
                "func": self.func,
                "threads": self.threads,
                "cmd": cmd,
            })
            print(f"[DONE] {self.TYPE}:{self.name}")
        except Exception as e:
            flags.mark_fail(meta={
                "task": self.TYPE,
                "name": self.name,
                "cmd": cmd,
                "error": str(e),
            })
            print(f"[FAIL] {self.TYPE}:{self.name}: {e}")
            raise


# ============================
# src/tasks/executor.py (mini)
# ============================
from dataclasses import dataclass, field
from typing import List

@dataclass
class Node:
    task: Task
    deps: List["Node"] = field(default_factory=list)

    def add(self, other: "Node") -> "Node":
        self.deps.append(other)
        return self

class Executor:
    def __init__(self, dry_run: bool = False):
        self.dry_run = dry_run
        self._seen: set[int] = set()

    def _run(self, n: Node) -> None:
        if id(n) in self._seen:
            return
        self._seen.add(id(n))
        for d in n.deps:
            self._run(d)
        n.task.run(dry_run=self.dry_run)

    def run(self, roots: List[Node]) -> None:
        for r in roots:
            self._run(r)

# End of files

# ================================
# src/tasks/fastqc/fastqc.py
# ================================
from __future__ import annotations
import shlex
from pathlib import Path
from typing import Dict, Any

from src.tasks.task import Task, TaskRegistry


def _singularity(image: str | Path, binds: list[str | Path], inner_cmd: str) -> str:
    binds_arg = ",".join(str(b) for b in binds)
    return f"singularity exec -B {shlex.quote(binds_arg)} {shlex.quote(str(image))} {inner_cmd}"


@TaskRegistry.register
class FastQCRunner(Task):
    TYPE = "fastqc"
    DEFAULTS: Dict[str, Any] = {
        "extract": True,
        "threads": 2,
        "fastqc_bin": "fastqc",
        "image": None,
        "binds": ["/storage", "/data"],
    }

    def to_sh(self) -> str:
        # inputs
        r1: Path = Path(self.inputs["read1"])  # required
        r2: Path | None = Path(self.inputs["read2"]) if self.inputs.get("read2") else None
        # outputs
        # we don't actually need the outputs to build cmd, but paths are validated post-run
        out_dir: Path = Path(self.workdir)
        out_dir.mkdir(parents=True, exist_ok=True)

        threads: int = int(self.params.get("threads", self.threads or 1))
        extract: bool = bool(self.params.get("extract", True))
        fastqc_bin: str = str(self.params.get("fastqc_bin", "fastqc"))
        image = self.params.get("image")
        binds = self.params.get("binds", ["/storage", "/data"]) or []

        inputs = f"{shlex.quote(str(r1))}"
        if r2:
            inputs += f" {shlex.quote(str(r2))}"

        flags = f"--threads {threads} -o {shlex.quote(str(out_dir))}"
        if extract:
            flags = f"{flags} --extract"

        inner = f"{fastqc_bin} {flags} {inputs}"
        if image:
            return _singularity(image, binds, inner)
        return inner


# ================================
# src/pipeline/loader.py
# ================================
from __future__ import annotations
import re
import yaml
from pathlib import Path
from typing import Any, Dict, List

from src.tasks.task import TaskRegistry
from src.tasks.executor import Node, Executor

PLACEHOLDER_PATTERN = re.compile(r"\{([A-Za-z0-9_]+)\}")


def render_template(s: str, context: Dict[str, Any]) -> str:
    def repl(m: re.Match[str]) -> str:
        key = m.group(1)
        if key not in context:
            raise KeyError(f"Missing placeholder '{key}' for template: {s}")
        return str(context[key])
    return PLACEHOLDER_PATTERN.sub(repl, s)


def _expand_mapping(m: Dict[str, Any], ctx: Dict[str, Any]) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for k, v in m.items():
        if v is None:
            out[k] = None
        elif isinstance(v, str):
            out[k] = render_template(v, ctx)
        elif isinstance(v, list):
            out[k] = [render_template(x, ctx) if isinstance(x, str) else x for x in v]
        elif isinstance(v, dict):
            out[k] = _expand_mapping(v, ctx)
        else:
            out[k] = v
    return out


def load_tasks_from_yaml(yaml_path: str | Path,
                         sample_id: str,
                         work_dir_path: str | Path) -> List[Node]:
    """Parse TASK_LIST YAML and build Node list.

    Supports placeholders: {sample_id}, {work_dir_path}, and {WORK_DIR} (task-local workdir).
    """
    with open(yaml_path, "r") as fh:
        data = yaml.safe_load(fh)

    tasks_cfg = data.get("TASK_LIST", [])
    nodes: List[Node] = []

    for item in tasks_cfg:
        # each item is like {'Rawdata_FastQC': { ... }}
        if not isinstance(item, dict) or len(item) != 1:
            raise ValueError(f"Invalid TASK_LIST item: {item}")
        task_name, cfg = next(iter(item.items()))

        tool = cfg["TOOL"]
        func = cfg.get("FUNC", "run")

        # resolve WORK_DIR first
        task_ctx = {"sample_id": sample_id, "work_dir_path": str(work_dir_path)}
        work_dir_str = render_template(cfg["WORK_DIR"], task_ctx)
        work_dir = Path(work_dir_str)

        # add WORK_DIR to context for I/O expansion
        ctx = {**task_ctx, "WORK_DIR": str(work_dir)}

        inputs = _expand_mapping(cfg.get("INPUT", {}), ctx)
        outputs_raw = _expand_mapping(cfg.get("OUTPUT", {}), ctx)

        # tolerate common typo keys like 'reqd2' -> map to 'read2' if present
        if "reqd2" in outputs_raw and "read2" not in outputs_raw:
            outputs_raw["read2"] = outputs_raw.pop("reqd2")

        params = cfg.get("PARAMS", {})

        TaskClass = TaskRegistry.get(tool)
        task = TaskClass(
            name=task_name,
            tool=tool,
            func=func,
            threads=int(params.get("threads", 1)),
            workdir=work_dir,
            inputs=inputs,
            outputs=outputs_raw,
            params=params,
        )
        nodes.append(Node(task=task))

    return nodes


# ================================
# src/pipeline/run_from_yaml.py (example entrypoint)
# ================================
from __future__ import annotations
import argparse
from pathlib import Path

from src.pipeline.loader import load_tasks_from_yaml
from src.tasks.executor import Executor


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("yaml", type=Path)
    ap.add_argument("sample_id", type=str)
    ap.add_argument("work_dir_path", type=Path)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    nodes = load_tasks_from_yaml(args.yaml, args.sample_id, args.work_dir_path)
    Executor(dry_run=args.dry_run).run(nodes)


if __name__ == "__main__":
    main()

# ================================
# src/tasks/picard/picard.py
# ================================
from __future__ import annotations
import importlib
import shlex
from pathlib import Path
from typing import Any, Dict

from src.tasks.task import Task, TaskRegistry


def _singularity(image: str | Path, binds: list[str | Path], inner_cmd: str) -> str:
    binds_arg = ",".join(str(b) for b in binds)
    return f"singularity exec -B {shlex.quote(binds_arg)} {shlex.quote(str(image))} {inner_cmd}"


_ALIAS = {
    # normalize FUNC values → module basename
    "deduplicate": "_markduplicates",
    "markduplicates": "_markduplicates",
    "rmdup": "_markduplicates",
    "downsampling": "_downsamplesam",
    "downsample": "_downsamplesam",
}


def _load_builder(func: str):
    key = _ALIAS.get(func.lower(), func.lower())
    mod = importlib.import_module(f"src.tasks.picard.{key}")
    # modules must expose a build_cmd(inputs, outputs, params) -> str
    return getattr(mod, "build_cmd")


@TaskRegistry.register
class PicardTask(Task):
    TYPE = "picard"
    DEFAULTS: Dict[str, Any] = {
        "picard_bin": "picard",
        "picard_jar": None,   # if set, uses `java -jar` path instead of picard_bin
        "java_xmx": "4g",
        "validation_stringency": "LENIENT",
        "create_index": True,
        "image": None,
        "binds": ["/storage", "/data"],
        "tmp_dir": None,
        # Function-specific defaults can live in modules, but may also be placed here
    }

    def to_sh(self) -> str:
        builder = _load_builder(self.func)
        inner = builder(self.inputs, self.outputs, {**self.DEFAULTS, **self.params, "threads": self.threads})
        image = self.params.get("image")
        binds = self.params.get("binds", ["/storage", "/data"]) or []
        if image:
            return _singularity(image, binds, inner)
        return inner


# ================================
# src/tasks/picard/_markduplicates.py
# ================================
from __future__ import annotations
import shlex
from pathlib import Path
from typing import Dict, Any


def _picard_prefix(params: Dict[str, Any]) -> str:
    jar = params.get("picard_jar")
    if jar:
        xmx = params.get("java_xmx", "4g")
        return f"java -Xmx{xmx} -jar {shlex.quote(str(jar))}"
    return shlex.quote(str(params.get("picard_bin", "picard")))


def build_cmd(inputs: Dict[str, Any], outputs: Dict[str, Any], params: Dict[str, Any]) -> str:
    """Build MarkDuplicates command.

    Required inputs:
      inputs["in_bam"]
    Required outputs:
      outputs["out_bam"]
    Optional outputs:
      outputs["metrics"] (if absent, default next to out_bam)
    """
    in_bam = Path(inputs["in_bam"]) if "in_bam" in inputs else Path(inputs["bam"])  # tolerate 'bam'
    out_bam = Path(outputs.get("out_bam") or outputs.get("bam") or (in_bam.with_suffix("") .with_suffix(".md.bam")))

    metrics = outputs.get("metrics")
    if not metrics:
        metrics = str(out_bam.with_suffix(out_bam.suffix + ".metrics.txt"))

    create_index = params.get("create_index", True)
    vs = params.get("validation_stringency", "LENIENT")
    tmp_dir = params.get("tmp_dir")
    remove_dup = params.get("remove_duplicates", False)
    assum_sorted = params.get("assume_sorted", False)
    optical_px = params.get("optical_duplicate_pixel_distance")

    args = [
        "MarkDuplicates",
        f"I={shlex.quote(str(in_bam))}",
        f"O={shlex.quote(str(out_bam))}",
        f"M={shlex.quote(str(metrics))}",
        f"VALIDATION_STRINGENCY={vs}",
    ]
    if remove_dup:
        args.append("REMOVE_DUPLICATES=true")
    if create_index:
        args.append("CREATE_INDEX=true")
    if assum_sorted:
        args.append("ASSUME_SORTED=true")
    if optical_px is not None:
        args.append(f"OPTICAL_DUPLICATE_PIXEL_DISTANCE={int(optical_px)}")
    if tmp_dir:
        args.append(f"TMP_DIR={shlex.quote(str(tmp_dir))}")

    return f"{_picard_prefix(params)} {' '.join(args)}"


# ================================
# src/tasks/picard/_downsamplesam.py
# ================================
from __future__ import annotations
import shlex
from pathlib import Path
from typing import Dict, Any

from ._markduplicates import _picard_prefix  # reuse prefix helper


def build_cmd(inputs: Dict[str, Any], outputs: Dict[str, Any], params: Dict[str, Any]) -> str:
    """Build DownsampleSam command.

    Required inputs:
      inputs["in_bam"] or inputs["bam"]
    Required outputs:
      outputs["out_bam"] or outputs["bam"]
    """
    in_bam = Path(inputs.get("in_bam") or inputs["bam"])
    out_bam = Path(outputs.get("out_bam") or outputs.get("bam") or (in_bam.with_suffix("") .with_suffix(".ds.bam")))

    p = params.get("probability")
    if p is None:
        raise ValueError("DownsampleSam requires 'probability' in PARAMS")
    seed = int(params.get("seed", 42))
    create_index = params.get("create_index", True)
    vs = params.get("validation_stringency", "LENIENT")
    tmp_dir = params.get("tmp_dir")

    args = [
        "DownsampleSam",
        f"I={shlex.quote(str(in_bam))}",
        f"O={shlex.quote(str(out_bam))}",
        f"P={p}",
        f"RANDOM_SEED={seed}",
        f"VALIDATION_STRINGENCY={vs}",
    ]
    if create_index:
        args.append("CREATE_INDEX=true")
    if tmp_dir:
        args.append(f"TMP_DIR={shlex.quote(str(tmp_dir))}")

    return f"{_picard_prefix(params)} {' '.join(args)}"


# ================================
# YAML snippets (examples)
# ================================
# - Dedup (MarkDuplicates)
#   - tool: picard, func: deduplicate
# - Downsample (DownsampleSam)
#
# TASK_LIST:
#   - Picard_Dedup:
#       TOOL: picard
#       FUNC: deduplicate
#       WORK_DIR: "{work_dir_path}/{sample_id}/02_Dedup"
#       INPUT:
#         in_bam: "{work_dir_path}/{sample_id}/01_Align/{sample_id}.sorted.bam"
#       OUTPUT:
#         out_bam: "{WORK_DIR}/{sample_id}.dedup.bam"
#         metrics: "{WORK_DIR}/{sample_id}.dedup.metrics.txt"
#       PARAMS:
#         image: /storage/images/picard-3.2.0.sif
#         binds: ["/storage","/data"]
#         create_index: true
#         validation_stringency: LENIENT
#         remove_duplicates: false
#
#   - Picard_Downsample:
#       TOOL: picard
#       FUNC: downsampling
#       WORK_DIR: "{work_dir_path}/{sample_id}/03_Downsample"
#       INPUT:
#         in_bam: "{work_dir_path}/{sample_id}/02_Dedup/{sample_id}.dedup.bam"
#       OUTPUT:
#         out_bam: "{WORK_DIR}/{sample_id}.ds0.5.bam"
#       PARAMS:
#         image: /storage/images/picard-3.2.0.sif
#         binds: ["/storage","/data"]
#         probability: 0.5
#         seed: 777

# ================================
# Directory layout for Picard tasks (modular per-FUNC)
# ================================
# src/tasks/picard/
# ├── __init__.py                 # alias map + importer
# ├── picard.py                   # PicardTask (registry entry)
# ├── markduplicates/
# │   ├── main.py                 # build_cmd() for MarkDuplicates
# │   ├── config.yaml             # defaults/validation hints
# │   └── test/
# │       ├── example_tasklist.yaml
# │       └── run.sh
# ├── deduplicates/               # alias to markduplicates (thin wrapper)
# │   ├── main.py
# │   ├── config.yaml
# │   └── test/
# │       ├── example_tasklist.yaml
# │       └── run.sh
# └── downsampling/
#     ├── main.py                 # build_cmd() for DownsampleSam
#     ├── config.yaml
#     └── test/
#         ├── example_tasklist.yaml
#         └── run.sh


# ================================
# src/tasks/picard/__init__.py
# ================================
from __future__ import annotations
import importlib
from typing import Callable

_ALIAS = {
    "markduplicates": "markduplicates",
    "deduplicates": "deduplicates",  # kept separate dir but can delegate
    "rmdup": "markduplicates",
    "downsampling": "downsampling",
    "downsample": "downsampling",
}


def load_builder(func: str) -> Callable:
    key = _ALIAS.get(func.lower(), func.lower())
    mod = importlib.import_module(f"src.tasks.picard.{key}.main")
    return getattr(mod, "build_cmd")


def load_defaults(func: str) -> dict:
    key = _ALIAS.get(func.lower(), func.lower())
    try:
        mod = importlib.import_module(f"src.tasks.picard.{key}.main")
        return getattr(mod, "DEFAULTS", {})
    except Exception:
        return {}


# ================================
# src/tasks/picard/picard.py (updated to use package importer)
# ================================
from __future__ import annotations
import shlex
from pathlib import Path
from typing import Any, Dict

from src.tasks.task import Task, TaskRegistry
from . import load_builder, load_defaults


def _singularity(image: str | Path, binds: list[str | Path], inner_cmd: str) -> str:
    binds_arg = ",".join(str(b) for b in binds)
    return f"singularity exec -B {shlex.quote(binds_arg)} {shlex.quote(str(image))} {inner_cmd}"


@TaskRegistry.register
class PicardTask(Task):
    TYPE = "picard"
    DEFAULTS: Dict[str, Any] = {
        "picard_bin": "picard",
        "picard_jar": None,
        "java_xmx": "4g",
        "validation_stringency": "LENIENT",
        "create_index": True,
        "image": None,
        "binds": ["/storage", "/data"],
        "tmp_dir": None,
    }

    def to_sh(self) -> str:
        # merge base defaults + per-FUNC defaults + user params
        func_defaults = load_defaults(self.func)
        merged_params = {**self.DEFAULTS, **func_defaults, **(self.params or {}), "threads": self.threads}
        builder = load_builder(self.func)
        inner = builder(self.inputs, self.outputs, merged_params)
        image = merged_params.get("image")
        binds = merged_params.get("binds", ["/storage", "/data"]) or []
        if image:
            return _singularity(image, binds, inner)
        return inner


# ================================
# src/tasks/picard/markduplicates/main.py
# ================================
from __future__ import annotations
import shlex
from pathlib import Path
from typing import Dict, Any

DEFAULTS: Dict[str, Any] = {
    "remove_duplicates": False,
    "assume_sorted": False,
    "optical_duplicate_pixel_distance": None,
}


def _picard_prefix(params: Dict[str, Any]) -> str:
    jar = params.get("picard_jar")
    if jar:
        xmx = params.get("java_xmx", "4g")
        return f"java -Xmx{xmx} -jar {shlex.quote(str(jar))}"
    return shlex.quote(str(params.get("picard_bin", "picard")))


def build_cmd(inputs: Dict[str, Any], outputs: Dict[str, Any], params: Dict[str, Any]) -> str:
    in_bam = Path(inputs.get("in_bam") or inputs["bam"])
    out_bam = Path(outputs.get("out_bam") or outputs.get("bam") or (in_bam.with_suffix("").with_suffix(".md.bam")))
    metrics = outputs.get("metrics") or str(out_bam.with_suffix(out_bam.suffix + ".metrics.txt"))

    args = [
        "MarkDuplicates",
        f"I={shlex.quote(str(in_bam))}",
        f"O={shlex.quote(str(out_bam))}",
        f"M={shlex.quote(str(metrics))}",
        f"VALIDATION_STRINGENCY={params.get('validation_stringency', 'LENIENT')}",
    ]
    if params.get("remove_duplicates", False):
        args.append("REMOVE_DUPLICATES=true")
    if params.get("create_index", True):
        args.append("CREATE_INDEX=true")
    if params.get("assume_sorted", False):
        args.append("ASSUME_SORTED=true")
    if params.get("optical_duplicate_pixel_distance") is not None:
        args.append(f"OPTICAL_DUPLICATE_PIXEL_DISTANCE={int(params['optical_duplicate_pixel_distance'])}")
    if params.get("tmp_dir"):
        args.append(f"TMP_DIR={shlex.quote(str(params['tmp_dir']))}")

    return f"{_picard_prefix(params)} {' '.join(args)}"


# ================================
# src/tasks/picard/markduplicates/config.yaml
# ================================
# defaults and hints for validation (optional)
remove_duplicates: false
assume_sorted: false
optical_duplicate_pixel_distance: null


# ================================
# src/tasks/picard/markduplicates/test/example_tasklist.yaml
# ================================
TASK_LIST:
  - Picard_Dedup:
      TOOL: picard
      FUNC: markduplicates
      WORK_DIR: "{work_dir_path}/{sample_id}/02_Dedup"
      INPUT:
        in_bam: "{work_dir_path}/{sample_id}/01_Align/{sample_id}.sorted.bam"
      OUTPUT:
        out_bam: "{WORK_DIR}/{sample_id}.dedup.bam"
        metrics: "{WORK_DIR}/{sample_id}.dedup.metrics.txt"
      PARAMS:
        image: /storage/images/picard-3.2.0.sif
        binds: ["/storage","/data"]
        create_index: true
        validation_stringency: LENIENT
        remove_duplicates: false


# ================================
# src/tasks/picard/markduplicates/test/run.sh
# ================================
#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR=$(git rev-parse --show-toplevel 2>/dev/null || pwd)
python -m src.pipeline.run_from_yaml \
  "$ROOT_DIR/src/tasks/picard/markduplicates/test/example_tasklist.yaml" \
  SAMPLE01 \
  /tmp/work --dry-run


# ================================
# src/tasks/picard/deduplicates/main.py (delegates to markduplicates)
# ================================
from __future__ import annotations
from src.tasks.picard.markduplicates.main import build_cmd as build_cmd  # re-export
from src.tasks.picard.markduplicates.main import DEFAULTS as DEFAULTS


# ================================
# src/tasks/picard/deduplicates/config.yaml
# ================================
# kept for symmetry; inherits defaults from markduplicates


# ================================
# src/tasks/picard/deduplicates/test/example_tasklist.yaml
# ================================
TASK_LIST:
  - Picard_Dedup:
      TOOL: picard
      FUNC: deduplicates
      WORK_DIR: "{work_dir_path}/{sample_id}/02_Dedup"
      INPUT:
        bam: "{work_dir_path}/{sample_id}/01_Align/{sample_id}.sorted.bam"
      OUTPUT:
        bam: "{WORK_DIR}/{sample_id}.dedup.bam"
      PARAMS:
        image: /storage/images/picard-3.2.0.sif
        binds: ["/storage","/data"]


# ================================
# src/tasks/picard/deduplicates/test/run.sh
# ================================
#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR=$(git rev-parse --show-toplevel 2>/dev/null || pwd)
python -m src.pipeline.run_from_yaml \
  "$ROOT_DIR/src/tasks/picard/deduplicates/test/example_tasklist.yaml" \
  SAMPLE01 \
  /tmp/work --dry-run


# ================================
# src/tasks/picard/downsampling/main.py
# ================================
from __future__ import annotations
import shlex
from pathlib import Path
from typing import Dict, Any

DEFAULTS: Dict[str, Any] = {
    "probability": 0.5,
    "seed": 42,
}


def _picard_prefix(params: Dict[str, Any]) -> str:
    jar = params.get("picard_jar")
    if jar:
        xmx = params.get("java_xmx", "4g")
        return f"java -Xmx{xmx} -jar {shlex.quote(str(jar))}"
    return shlex.quote(str(params.get("picard_bin", "picard")))


def build_cmd(inputs: Dict[str, Any], outputs: Dict[str, Any], params: Dict[str, Any]) -> str:
    in_bam = Path(inputs.get("in_bam") or inputs["bam"])
    out_bam = Path(outputs.get("out_bam") or outputs.get("bam") or (in_bam.with_suffix("").with_suffix(".ds.bam")))

    p = params.get("probability")
    if p is None:
        raise ValueError("DownsampleSam requires 'probability'")

    args = [
        "DownsampleSam",
        f"I={shlex.quote(str(in_bam))}",
        f"O={shlex.quote(str(out_bam))}",
        f"P={p}",
        f"RANDOM_SEED={int(params.get('seed', 42))}",
        f"VALIDATION_STRINGENCY={params.get('validation_stringency', 'LENIENT')}",
    ]
    if params.get("create_index", True):
        args.append("CREATE_INDEX=true")
    if params.get("tmp_dir"):
        args.append(f"TMP_DIR={shlex.quote(str(params['tmp_dir']))}")

    return f"{_picard_prefix(params)} {' '.join(args)}"


# ================================
# src/tasks/picard/downsampling/config.yaml
# ================================
probability: 0.5
seed: 42


# ================================
# src/tasks/picard/downsampling/test/example_tasklist.yaml
# ================================
TASK_LIST:
  - Picard_Downsample:
      TOOL: picard
      FUNC: downsampling
      WORK_DIR: "{work_dir_path}/{sample_id}/03_Downsample"
      INPUT:
        in_bam: "{work_dir_path}/{sample_id}/02_Dedup/{sample_id}.dedup.bam"
      OUTPUT:
        out_bam: "{WORK_DIR}/{sample_id}.ds0.5.bam"
      PARAMS:
        image: /storage/images/picard-3.2.0.sif
        binds: ["/storage","/data"]
        probability: 0.5
        seed: 777


# ================================
# src/tasks/picard/downsampling/test/run.sh
# ================================
#!/usr/bin/env bash
set -euo pipefail
ROOT_DIR=$(git rev-parse --show-toplevel 2>/dev/null || pwd)
python -m src.pipeline.run_from_yaml \
  "$ROOT_DIR/src/tasks/picard/downsampling/test/example_tasklist.yaml" \
  SAMPLE01 \
  /tmp/work --dry-run
