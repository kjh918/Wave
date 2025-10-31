# src/utils/flags.py
from __future__ import annotations
import functools, json, time
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional

# ------------------------------
# 내부 유틸
# ------------------------------
def _iter_output_files(outputs: Any) -> Iterable[Path]:
    if outputs is None:
        return []
    if isinstance(outputs, Mapping):
        for k, v in outputs.items():
            if k == "dir":
                continue
            if isinstance(v, (str, Path)):
                yield Path(v)
    elif isinstance(outputs, (list, tuple)):
        for v in outputs:
            if isinstance(v, (str, Path)):
                yield Path(v)
    elif isinstance(outputs, (str, Path)):
        yield Path(outputs)

def _all_files_ok(paths: Iterable[Path]) -> bool:
    for p in paths:
        try:
            if (not p.is_file()) or p.stat().st_size <= 0:
                return False
        except FileNotFoundError:
            return False
    return True

def _get_task(obj) -> Optional[Any]:
    # kwargs로 전달되었거나, (self, task) 시그니처 등 일반 케이스 지원
    return (
        getattr(obj, "task", None) if hasattr(obj, "task") else None
    )

def _extract_task_from_args_kwargs(args, kwargs) -> Optional[Any]:
    if "task" in kwargs:
        return kwargs["task"]
    # (self, task) or (task,)
    for a in args:
        if hasattr(a, "outputs") and hasattr(a, "workdir"):
            return a
    return None

# ------------------------------
# 실행 전: 이미 완료되었으면 스킵
# ------------------------------
def skip_if_done(flag_name: str = ".done", require_outputs_ok: bool = True, allow_stale_flag: bool = False):
    """
    - <workdir>/<flag_name> 가 존재하고,
    - require_outputs_ok=True 면 outputs 파일들이 모두 존재 & size>0
    인 경우, 실제 실행을 스킵하고 True를 반환.
    allow_stale_flag=True 이면 .done만 있으면 스킵(출력검증 생략)
    """
    def deco(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            task = _extract_task_from_args_kwargs(args, kwargs)
            if task is None:
                return func(*args, **kwargs)

            workdir = Path(getattr(task, "workdir", "."))
            done_flag = workdir / flag_name

            if done_flag.exists():
                if allow_stale_flag or not require_outputs_ok:
                    # 플래그만 있으면 스킵
                    return True
                # 출력 검증
                files = list(_iter_output_files(getattr(task, "outputs", {})))
                if _all_files_ok(files):
                    # 정상 완료 → 스킵
                    return True
                # 플래그는 있는데 출력이 깨진 경우 → 계속 진행(재실행)
            return func(*args, **kwargs)
        return wrapper
    return deco

# ------------------------------
# 실행 후: 완료/실패 플래그 생성
# ------------------------------
def flag_on_complete(flag_name: str = ".done", fail_flag: str = ".failed", write_meta: bool = True):
    def deco(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            result = func(*args, **kwargs)

            task = _extract_task_from_args_kwargs(args, kwargs)
            if task is None:
                return result

            workdir = Path(getattr(task, "workdir", "."))
            outputs = getattr(task, "outputs", {})
            files = list(_iter_output_files(outputs))
            when = time.strftime("%Y-%m-%d %H:%M:%S")

            ok = _all_files_ok(files)
            if ok:
                flag_path = workdir / flag_name
                flag_path.write_text("OK\n")
                if write_meta:
                    (workdir / f"{flag_name}.json").write_text(
                        json.dumps({"status":"OK","timestamp":when,"outputs":[str(p) for p in files]},
                                   indent=2, ensure_ascii=False)
                    )
            else:
                if fail_flag:
                    (workdir / fail_flag).write_text("FAILED\n")
                    if write_meta:
                        (workdir / f"{fail_flag}.json").write_text(
                            json.dumps({"status":"FAILED","timestamp":when,
                                        "missing_or_empty":[str(p) for p in files
                                                            if (not p.exists() or p.stat().st_size<=0)]},
                                       indent=2, ensure_ascii=False)
                        )
            return result
        return wrapper
    return deco

# ------------------------------
# 도우미: 외부에서 직접 체크할 때
# ------------------------------
def is_done(task, flag_name: str = ".done", require_outputs_ok: bool = True) -> bool:
    workdir = Path(getattr(task, "workdir", "."))
    if not (workdir / flag_name).exists():
        return False
    if not require_outputs_ok:
        return True
    files = list(_iter_output_files(getattr(task, "outputs", {})))
    return _all_files_ok(files)
