# wave/utils/flags.py
from __future__ import annotations
import time
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, List, Union, Callable, Optional
import functools

OutputsLike = Union[str, Path, Iterable[str], Mapping[str, Any]]

# -------------------------
# 유틸: outputs -> 파일 경로 리스트로 정규화
# - dict이면 value들 중 'dir' 키는 무시
# - list/tuple이면 요소들을 path로 변환
# - str/path이면 단일 파일로 처리
# -------------------------
def _normalize_outputs(outputs: OutputsLike) -> List[Path]:
    if outputs is None:
        return []
    # dict-like: skip "dir" key
    if isinstance(outputs, Mapping):
        paths = []
        for k, v in outputs.items():
            if k == "dir" or v is None:
                continue
            if isinstance(v, (str, Path)):
                paths.append(Path(v))
        return paths
    # iterable (list/tuple/generator) of strings/paths
    if isinstance(outputs, (list, tuple, set)):
        return [Path(p) for p in outputs if p is not None and p != "dir"]
    # single path string
    return [Path(outputs)]

# -------------------------
# 파일 상태 검사
# -------------------------
def _all_outputs_exist_and_nonzero(paths: Iterable[Path]) -> bool:
    p_list = list(paths)
    if not p_list:
        # outputs가 비어있으면 '검사할 파일이 없음'으로 간주 (caller 정책에 따라 변경 가능)
        return False
    for p in p_list:
        if not p.is_file():
            return False
        try:
            if p.stat().st_size <= 0:
                return False
        except Exception:
            return False
    return True

# -------------------------
# 단순 플래그 생성 함수
# - outputs: dict|list|str 형태로 직접 전달
# - workdir: 플래그를 만들 디렉토리(주로 task workdir)
# - done_name / failed_name: 파일명(기본 ".done", ".failed")
# - write_meta: .json 메타 파일 생성 여부
# 반환: True if marked done, False if marked failed
# -------------------------
def mark_flag_for_outputs(
    outputs: OutputsLike,
    workdir: Union[str, Path],
    done_name: str = ".done",
    failed_name: str = ".failed",
    write_meta: bool = True,
) -> bool:
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)

    paths = _normalize_outputs(outputs)

    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    ok = _all_outputs_exist_and_nonzero(paths)

    if ok:
        done_path = wd / done_name
        done_path.write_text("OK\n")
        if write_meta:
            meta = {
                "status": "OK",
                "timestamp": timestamp,
                "outputs": [str(p) for p in paths],
            }
            (wd / f"{done_name}.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False))
        # remove any stale failed flag
        failed_path = wd / failed_name
        if failed_path.exists():
            try:
                failed_path.unlink()
            except Exception:
                pass
        return True
    else:
        failed_path = wd / failed_name
        failed_path.write_text("FAILED\n")
        if write_meta:
            missing = [str(p) for p in paths if not p.exists() or p.stat().st_size <= 0]
            meta = {
                "status": "FAILED",
                "timestamp": timestamp,
                "missing_or_empty": missing,
                "checked_outputs": [str(p) for p in paths],
            }
            (wd / f"{failed_name}.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False))
        # do not remove .done automatically here
        return False

# -------------------------
# 간단 확인용 함수
# - outputs를 주면 .done 여부만 반환 (workdir에서 .done 파일 직접 확인)
# -------------------------
def is_done(workdir: Union[str, Path], done_name: str = ".done") -> bool:
    return (Path(workdir) / done_name).exists()

# -------------------------
# 편의: 데코레이터 (간단)
# - 함수 실행 후 outputs/workdir 를 찾아 mark_flag_for_outputs 호출
# - 사용법:
#     @auto_flag_on_complete(outputs_arg="outputs", workdir_arg="workdir")
#     def my_runner(..., outputs=..., workdir=...): ...
#
# - 이 데코레이터는 아주 단순하게 kwargs를 살펴봄.
# -------------------------
def auto_flag_on_complete(outputs_arg: str = "outputs", workdir_arg: str = "workdir",
                          done_name: str = ".done", failed_name: str = ".failed",
                          write_meta: bool = True):
    def deco(func: Callable):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            result = func(*args, **kwargs)
            
            # 우선 kwargs에서 직접 찾기
            outputs = kwargs.get(outputs_arg, None)
            workdir = kwargs.get(workdir_arg, None)

            # 다음으로 args에서 (positionally) 찾아보기 (함수 시그니처에 따라 다름)
            if outputs is None or workdir is None:
                # naive scan: dict-like element in args that looks like outputs/workdir
                for a in args:
                    if outputs is None and isinstance(a, (dict, list, tuple, str, Path)):
                        # heuristics: if dict and contains 'outputs' key -> skip (we want direct outputs)
                        if isinstance(a, dict) and ("outputs" in a or "workdir" in a):
                            continue
                        outputs = outputs or a
                    # workdir detection: if arg is string and looks like a path to existing dir
                    if workdir is None and isinstance(a, (str, Path)):
                        ap = Path(a)
                        if ap.exists() and ap.is_dir():
                            workdir = str(ap)

            # if still missing, give up (no flagging)
            if outputs is None or workdir is None:
                return result

            # call simple marker
            try:
                mark_flag_for_outputs(outputs=outputs, workdir=workdir,
                                      done_name=done_name, failed_name=failed_name,
                                      write_meta=write_meta)
            except Exception:
                # 플래깅 실패해도 원래 함수 결과는 그대로 반환
                pass
            return result
        return wrapper
    return deco