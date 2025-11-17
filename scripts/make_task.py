#!/usr/bin/env python3
from __future__ import annotations
import yaml
from pathlib import Path
import argparse
import os

TEMPLATE = """from __future__ import annotations
from pathlib import Path
from typing import Dict, Any, List, Optional

from wave.core.task import Task
from wave.core.task_registry import register_task
from wave.utils.task_utils import (
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    join_argv_lines,
)

@register_task("{tool}.{func}")
class {class_name}(Task):
    TYPE = "{tool}.{func}"

    INPUTS = {inputs}

    OUTPUTS = {outputs}

    DEFAULTS = {defaults}

    def _build_cmd(self, *, inputs, outputs, params, threads, workdir, sample_id=None):
        # TODO: implement command builder
        raise NotImplementedError("Implement command builder in {class_name}._build_cmd")

    def to_sh(self) -> List[str]:
        p = {**self.DEFAULTS, **(self.params or {})}
        th = int(self.threads or p.get("threads", 1))

        lines = self._build_cmd(
            inputs=self.inputs,
            outputs=self.outputs,
            params=p,
            threads=th,
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id"),
        )
        return join_argv_lines(lines)
"""

def generate_from_meta(meta_path: Path, force: bool = False):
    with open(meta_path, "r") as f:
        meta = yaml.safe_load(f)

    tool = meta.get("tool")
    func = meta.get("func", "main")
    class_name = meta.get("class_name", f"{tool.capitalize()}Task")

    inputs = meta.get("inputs", {})
    outputs = meta.get("outputs", {})
    defaults = meta.get("defaults", {})

    main_py_path = meta_path.parent / "main.py"

    if main_py_path.exists() and not force:
        print(f"[SKIP] {main_py_path} already exists (use --force to overwrite)")
        return

    rendered = TEMPLATE.format(
        tool=tool,
        func=func,
        class_name=class_name,
        inputs=inputs,
        outputs=outputs,
        defaults=defaults,
    )

    with open(main_py_path, "w") as f:
        f.write(rendered)

    print(f"[OK] Generated {main_py_path}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Path to meta.yaml or directory")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    target = Path(args.path)

    if target.is_file() and target.name == "meta.yaml":
        generate_from_meta(target, force=args.force)
        return

    for meta_path in target.rglob("meta.yaml"):
        generate_from_meta(meta_path, force=args.force)

if __name__ == "__main__":
    main()