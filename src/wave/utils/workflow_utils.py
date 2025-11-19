# src/workflow_utils.py
from __future__ import annotations
import re, shlex, importlib, pkgutil
from pathlib import Path
from typing import Any, Dict, Iterable, List

from wave.core.task_registry import TaskRegistry
from wave.core.task import Task


# --------------------------
# task auto-load
# --------------------------
def autoload_tasks(package_root: str = "wave.tasks") -> None:
    """
    wave.tasks 하위에서 *.main 모듈을 자동으로 import.
    각 main.py 안에서 @register_task(...) 가 실행되면서
    TaskRegistry 에 등록된다.
    """
    try:
        pkg = importlib.import_module(package_root)
    except ModuleNotFoundError as e:
        print(f"[WAVE] cannot import package_root='{package_root}': {e}")
        return
    
    for m in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):

        name = m.name

        # 1) main.py만 골라서 import
        if name.endswith(".main"):
            try:
                # print("    import:", name)
                importlib.import_module(name)
            except Exception as e:
                print(f"[WAVE] autoload skip {name}: {e}")
            continue

        # 2) 패키지인 경우, 그 아래의 main들을 위해 패키지만 먼저 로딩
        if m.ispkg:
            try:
                importlib.import_module(name)
            except Exception as e:
                print(f"[WAVE] autoload skip pkg {name}: {e}")
                continue


def skip_finished_tasks(output_dict: Dict[str, Any], done_flag: Path) -> bool:
    """
    output_dict 에서 'dir' key는 제외하고,
    나머지 outputs 의 모든 파일이 존재 & size>0 이면
    done_flag(.done)를 생성하고 True를 반환.
    그렇지 않으면 False.
    """
    outputs = output_dict or {}

    # dir 제외
    file_outputs = {
        k: v for k, v in outputs.items()
        if k != "dir" and v and isinstance(v, str)
    }

    output_count = len(file_outputs)
    ok_count = 0

    for v in file_outputs.values():
        p = Path(v)
        if p.is_file() and p.stat().st_size > 0:
            ok_count += 1

    if output_count > 0 and ok_count == output_count:
        done_flag.write_text("OK\n")
        return True

    return False


# --------------------------
# 템플릿 / 플레이스홀더 관련
# --------------------------
def parse_placeholders(tpl: str) -> List[str]:
    return re.findall(r"\{([a-zA-Z0-9_]+)\}", tpl)


def tpl_to_regex(tpl: str, placeholders: List[str]) -> re.Pattern:
    token_re = re.compile(r"\{(\w+)\}")
    seen = set()
    parts: List[str] = []
    last = 0

    for m in token_re.finditer(tpl):
        parts.append(re.escape(tpl[last:m.start()]))

        name = m.group(1)
        if name in placeholders:
            if name in seen:
                parts.append(rf"(?P={name})")
            else:
                parts.append(rf"(?P<{name}>[^/]+)")
                seen.add(name)
        else:
            parts.append(re.escape(m.group(0)))
        last = m.end()

    parts.append(re.escape(tpl[last:]))
    pattern = "^" + "".join(parts) + "$"
    return re.compile(pattern)


def tpl_to_glob(tpl: str, placeholders: List[str]) -> str:
    pat = tpl
    for ph in placeholders:
        pat = pat.replace(f"{{{ph}}}", "*")
    return pat


def render_value(val: Any, ctx: Dict[str, Any]) -> Any:
    if isinstance(val, str):
        return re.sub(r"\{([a-zA-Z0-9_]+)\}", lambda m: str(ctx.get(m.group(1), "")), val)
    if isinstance(val, list):
        return [render_value(v, ctx) for v in val]
    if isinstance(val, dict):
        return {k: render_value(v, ctx) for k, v in val.items()}
    return val


# --------------------------
# 쉘 라인 조인
# --------------------------
def sh_join(argv: Iterable[str]) -> str:
    SPECIAL_TOKENS = {">", ">>", "<", "|", "2>", "&>", "&&", "||"}
    parts: List[str] = []
    for arg in argv:
        s = str(arg)
        if s in SPECIAL_TOKENS:
            parts.append(s)
        else:
            parts.append(shlex.quote(s))
    return " ".join(parts)


# --------------------------
# Task 인스턴스 생성 헬퍼
# --------------------------
def instantiate_task(
        TaskCls: type[Task],
        *,
        name: str,
        tool: str,
        func: str,
        threads: int,
        workdir: Path,
        inputs: Dict[str, Any],
        outputs: Dict[str, Any],
        params: Dict[str, Any],
    ) -> Task:
    """
    현재 Task 시그니처:
        __init__(self, name, tool, func, threads, workdir, inputs=None, outputs=None, params=None)
    """
    try:
        return TaskCls(
            name=name,
            tool=tool,
            func=func,
            threads=int(threads),
            workdir=workdir,
            inputs=inputs,
            outputs=outputs,
            params=params,
        )
    except TypeError as e:
        raise TypeError(f"[Task instantiate] Cannot create '{name}' ({tool}.{func}) → {e}")