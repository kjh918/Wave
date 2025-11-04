# src/tasks/util.py
from __future__ import annotations
from pathlib import Path
from typing import Any, Iterable, Sequence, Optional, List, Union, Callable , Dict
import shlex

def ensure_dir(path: Union[str, Path]) -> str:
    """디렉토리 생성 보장 후 문자열 경로 반환."""
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p.as_posix()

def normalize_binds(binds: Any) -> Optional[List[str]]:
    """
    - None → None
    - 'a,b' → ['a','b']
    - ['a','b'] → ['a','b']
    - 기타 → None
    """
    if binds is None:
        return None
    if isinstance(binds, str):
        vals = [x.strip() for x in binds.split(",") if x.strip()]
        return vals or None
    if isinstance(binds, (list, tuple)):
        return [str(x) for x in binds]
    return None

def singularity_exec_cmd(
        *,
        image: str,
        argv: Sequence[str],
        binds: Optional[Sequence[str]] = None,
        singularity_bin: str = "singularity",
    ) -> List[str]:
    """Singularity exec 커맨드 토큰 배열 생성."""
    cmd: List[str] = [singularity_bin, "exec"]
    for b in (binds or []):
        cmd += ["-B", str(b)]
    cmd.append(str(image))
    cmd += list(map(str, argv))
    return cmd

def join_argv_lines(lines: Iterable[Union[str, Sequence[Any]]]) -> List[str]:
    """
    argv 시퀀스/문자열 혼재를 안전하게 문자열 라인 리스트로 통일.
    """
    out: List[str] = []
    for ln in lines:
        if isinstance(ln, (list, tuple)):
            out.append(shlex.join([str(t) for t in ln]))
        else:
            out.append(str(ln))
    return out

def to_sh_from_builder(
    *,
        builder: Callable[..., Iterable[Union[str, Sequence[str]]]],
        inputs: Dict[str, Any],
        outputs: Dict[str, Any],
        params: Dict[str, Any],
        threads: int,
        workdir: str,
        sample_id: Optional[str] = None,
        ensure_output_dir_key: Optional[str] = None,
    ) -> List[str]:
    # 필요 시 결과 디렉토리 생성
    if ensure_output_dir_key:
        out_dir = outputs.get(ensure_output_dir_key)
        if out_dir:
            outputs[ensure_output_dir_key] = ensure_dir(out_dir)
    # 빌더 호출 → 라인 합치기
    lines = builder(
        inputs=inputs,
        outputs=outputs,
        params=params,
        threads=threads,
        workdir=workdir,
        sample_id=sample_id,
    )
    return join_argv_lines(lines)