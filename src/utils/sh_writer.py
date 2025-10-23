# utils/sh_writer.py (미니멀)
from __future__ import annotations
from pathlib import Path
from typing import Iterable, Sequence
import shlex, datetime

def _line(x) -> str:
    # x가 토큰 리스트면 join, 문자열이면 그대로
    return shlex.join(x) if isinstance(x, (list, tuple)) else str(x)

def write_script_from_cmds(cmds: Iterable[Sequence[str] | str],
                           out_path: str | Path,
                           set_x: bool = False) -> Path:
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "trap 'echo \"[ERR] $(date +%F-%T) $0:$LINENO\" >&2' ERR",
        f"# generated: {datetime.datetime.now().isoformat(timespec='seconds')}",
    ]
    if set_x:
        header.append("set -x")    # 디버그 모드
    body = [_line(c) for c in cmds]
    out.write_text("\n".join([*header, *body, ""]) + "\n")
    out.chmod(0o755)
    return out