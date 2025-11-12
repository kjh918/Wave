from __future__ import annotations
import functools, json, time
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, List

# ------------------------------
# 내부 유틸
# ------------------------------
def _iter_output_files(outputs: Any) -> List[Path]:
    paths: List[Path] = []
    if outputs is None:
        return paths
    if isinstance(outputs, Mapping):
        for k, v in outputs.items():
            if k == "dir":
                continue
            if isinstance(v, (str, Path)):
                paths.append(Path(v))
    elif isinstance(outputs, (list, tuple)):
        for v in outputs:
            if isinstance(v, (str, Path)):
                paths.append(Path(v))
    elif isinstance(outputs, (str, Path)):
        paths.append(Path(outputs))
    return paths

def _missing_or_empty(paths: Iterable[Path]) -> List[str]:
    bad: List[str] = []
    for p in paths:
        try:
            if (not p.is_file()) or p.stat().st_size <= 0:
                bad.append(str(p))
        except FileNotFoundError:
            bad.append(str(p))
    return bad

def _all_files_ok(paths: Iterable[Path]) -> bool:
    paths = list(paths)
    if len(paths) == 0:
        # ✅ 출력이 하나도 없으면 "성공"으로 간주하지 않음
        return False
    return len(_missing_or_empty(paths)) == 0

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
    - require_outputs_ok=True 면 outputs 파일들이 모두 존재 & size>0 (그리고 최소 1개 이상)
    인 경우, 실제 실행을 스킵하고 True를 반환.
    - allow_stale_flag=True 이면 .done만 있으면 스킵(출력검증 생략)
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
                    return True
                files = _iter_output_files(getattr(task, "outputs", {}))
                if _all_files_ok(files):
                    return True
            return func(*args, **kwargs)
        return wrapper
    return deco

# ------------------------------
# 실행 후: 완료/실패 플래그 생성
# ------------------------------
def flag_on_complete(flag_name: str = ".done", fail_flag: str = ".failed", write_meta: bool = True):
    """
    성공 기준:
      - outputs에 지정된 파일이 >=1개 존재하고,
      - 전부 존재하며 size > 0
    실패 기준:
      - outputs가 비었거나, 하나라도 없거나 size==0
    """
    def deco(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            result = func(*args, **kwargs)

            task = _extract_task_from_args_kwargs(args, kwargs)
            if task is None:
                return result

            workdir = Path(getattr(task, "workdir", "."))
            outputs = getattr(task, "outputs", {})
            files = _iter_output_files(outputs)
            when = time.strftime("%Y-%m-%d %H:%M:%S")

            ok = _all_files_ok(files)
            if ok:
                if flag_name:
                    (workdir / flag_name).write_text("OK\n")
                    if write_meta:
                        (workdir / f"{flag_name}.json").write_text(
                            json.dumps(
                                {"status": "OK", "timestamp": when, "outputs": [str(p) for p in files]},
                                indent=2, ensure_ascii=False
                            )
                        )
            else:
                bad = _missing_or_empty(files)
                reason = "no_outputs_declared" if len(files) == 0 else "missing_or_empty_outputs"
                if fail_flag:
                    (workdir / fail_flag).write_text("FAILED\n")
                    if write_meta:
                        (workdir / f"{fail_flag}.json").write_text(
                            json.dumps(
                                {
                                    "status": "FAILED",
                                    "timestamp": when,
                                    "reason": reason,
                                    "missing_or_empty": bad,
                                },
                                indent=2, ensure_ascii=False
                            )
                        )
            return result
        return wrapper
    return deco

# ------------------------------
# 외부 체크
# ------------------------------
def is_done(task, flag_name: str = ".done", require_outputs_ok: bool = True) -> bool:
    workdir = Path(getattr(task, "workdir", "."))
    if not (workdir / flag_name).exists():
        return False
    if not require_outputs_ok:
        return True
    files = _iter_output_files(getattr(task, "outputs", {}))
    return _all_files_ok(files)