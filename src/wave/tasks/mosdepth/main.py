# wave/tasks/mosdepth/mean_depth/main.py
from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional
import os
import shlex

from wave.core.task import Task
from wave.core.task_registry import register_task
from wave.utils.task_utils import (
    ensure_dir,
    join_argv_lines,
)

@register_task("mosdepth.mean_depth")
class MosdepthMeanDepthTask(Task):
    """
    mosdepth로 전체 mean depth 계산

    INPUTS
      bam: 입력 BAM

    OUTPUTS
      dir     : mosdepth 결과 디렉토리
      prefix  : mosdepth prefix
      mean_txt: mean depth를 저장한 txt 파일

    PARAMS
      mosdepth_bin: mosdepth 실행 파일
      threads     : 스레드 수
      extra_args  : ["--no-per-base", "--fast-mode"] 등 추가 옵션
    """

    TYPE = "mosdepth.mean_depth"

    INPUTS: Dict[str, Any] = {
        "bam": {"type": "path", "required": True, "desc": "Input BAM"},
    }
    OUTPUTS: Dict[str, Any] = {
        "dir":     {"type": "dir",  "required": False, "desc": "Output dir"},
        "prefix":  {"type": "path", "required": False, "desc": "Mosdepth prefix (no extension)"},
        "mean_txt": {"type": "path","required": False, "desc": "Mean depth text file"},
    }
    DEFAULTS: Dict[str, Any] = {
        "mosdepth_bin": "/storage/apps/mosdepth-0.3.5/mosdepth",
        "threads": 4,
        "extra_args": ["--no-per-base", "--fast-mode"],
    }

    def _build_cmd(
        self,
        *,
        inputs: Dict[str, Any],
        outputs: Dict[str, Any],
        params: Dict[str, Any],
        threads: int,
        workdir: str,
        sample_id: Optional[str] = None,
    ) -> List[Sequence[str] | str]:

        bam = inputs["bam"]
        out_dir = ensure_dir(outputs.get("dir") or workdir)

        base = sample_id or os.path.splitext(os.path.basename(str(bam)))[0]
        prefix = outputs.get("prefix") or os.path.join(out_dir, base)
        mean_txt = outputs.get("mean_txt") or os.path.join(out_dir, f"{base}.mean_depth.txt")

        mosdepth_bin = str(params.get("mosdepth_bin", self.DEFAULTS["mosdepth_bin"]))
        th = int(params.get("threads", threads or self.DEFAULTS["threads"]))

        extra_args = params.get("extra_args", self.DEFAULTS["extra_args"])
        if isinstance(extra_args, str):
            extra_args = extra_args.split()

        # 1) mosdepth 실행
        argv: List[str] = [mosdepth_bin, "-t", str(th)]
        argv += list(map(str, extra_args))
        argv += [prefix, str(bam)]

        # 2) summary에서 길이 가중 평균 계산해서 mean_txt에 저장
        summary = f"{prefix}.mosdepth.summary.txt"
        awk_cmd = (
            "awk 'NR>1{L+=$2; S+=$2*$4} "
            "END{if(L>0) printf(\"%.6f\\n\", S/L); else print 0}' "
            f"{shlex.quote(summary)} > {shlex.quote(mean_txt)}"
        )

        return [argv, awk_cmd]

    def to_sh(self) -> List[str]:
        p = {**self.DEFAULTS, **(self.params or {})}
        lines = self._build_cmd(
            inputs=self.inputs or {},
            outputs=self.outputs or {"dir": str(self.workdir)},
            params=p,
            threads=int(self.threads or p.get("threads", 4)),
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id"),
        )
        return join_argv_lines(lines)