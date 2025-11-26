# src/tasks/gatk/applybqsr/main.py
from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional
import os
import shlex

from wave.core.task import Task
from wave.core.task_registry import register_task
from wave.utils.task_utils import (
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    to_sh_from_builder,
)

@register_task("glimpse2.phase")
class Glimpse2PhaseTask(Task):
    TYPE = "glimpse2.phase"

    INPUTS = {
        # BAM 파일
        "bam": {"type": "path", "required": True, "desc": "Input BAM"},
        # reference prefix: e.g. reference_panel/split/1000GP.chr22.noNA12878
        "reference": {"type": "path", "required": True, "desc": "Reference prefix (without region suffix)"},
        # chunks 파일: e.g. chunks.chr22.txt
        "chunk": {"type": "path", "required": True, "desc": "Chunk definition file"},
        # 선택: map 파일이 필요한 경우
        "map_file": {"type": "path", "required": False, "desc": "Genetic map file for GLIMPSE2_phase"},
    }

    OUTPUTS = {
        # OUT prefix: e.g. GLIMPSE_impute/NA12878_imputed
        "bcf_prefix": {"type": "path", "required": True, "desc": "Prefix for imputed BCF outputs"},
    }

    DEFAULTS: Dict[str, Any] = {
        "glimpse2_phase": "/bin/GLIMPSE2_phase",
        "image": "/storage/home/jhkim/Apps/GLIMPSE/glimpse2.sif",
        "binds": ["/storage", "/data"],
        "singularity_bin": "singularity",
    }

    def _build_cmd(
        self,
        *,
        inputs,
        outputs,
        params,
        threads,
        workdir,
        sample_id: Optional[str] = None,
    ) -> List[Sequence[str] | str]:

        # ===== INPUT / OUTPUT =====
        bam = inputs["bam"]
        reference_prefix = inputs["reference"]      # REF=...
        out_prefix = outputs["bcf_prefix"]        # OUT=...
        chunks_file = inputs.get("chunk") 
        # ===== PARAMS =====
        glimpse2_phase = str(params.get("glimpse2_phase", "/bin/GLIMPSE2_phase"))
        image = params.get("image")
        binds = normalize_binds(params.get("binds"))
        singularity_bin = str(params.get("singularity_bin", "singularity"))

        # ensure output directory exists
        ensure_dir(os.path.dirname(out_prefix))

        # singularity exec prefix (옵션)
        singularity_cmd: List[str] = []
        if image:
            singularity_cmd = singularity_exec_cmd(
                image=str(image),
                argv=[],
                binds=binds,
                singularity_bin=singularity_bin,
            )

        # 쉘에서 사용할 변수 값은 미리 안전하게 quote
        ref_q = shlex.quote(reference_prefix)
        bam_q = shlex.quote(bam)
        out_q = shlex.quote(out_prefix)
        chunks_q = shlex.quote(chunks_file)
        phase_prog_q = shlex.quote(glimpse2_phase)

        # singularity를 쓰면: singularity exec ... image /bin/GLIMPSE2_phase
        # 안 쓰면: /bin/GLIMPSE2_phase
        if singularity_cmd:
            sing_prefix = " ".join(shlex.quote(x) for x in singularity_cmd)
            phase_call = f"{sing_prefix} {phase_prog_q}"
        else:
            phase_call = phase_prog_q

        # <<< 여기서부터는 bash 스크립트 >>>
        # f-string 안에서 ${VAR} 를 쓰기 위해 {{ }} 로 감싸서 이스케이프
        script = f"""REF={ref_q}
BAM={bam_q}
OUT={out_q}

while read -r ID CHR IRG ORG; do
    # IRG 예: chr22:100000-200000
    REGION_NOCHR="${{IRG#*:}}"
    REGS="${{REGION_NOCHR%-*}}"
    REGE="${{REGION_NOCHR#*-}}"

    REF_BIN="${{REF}}_${{CHR}}_${{REGS}}_${{REGE}}.bin"
    OUT_BCF="${{OUT}}_${{CHR}}_${{REGS}}_${{REGE}}.bcf"

    {phase_call} \\
        --bam-file "${{BAM}}" \\
        --reference "${{REF_BIN}}" \\
        --output "${{OUT_BCF}}" \\
        --threads {int(threads)}
done < {chunks_q}
"""

        # Wave는 이 문자열을 bash로 실행
        return [script]

    def to_sh(self) -> List[str]:
        p = {**self.DEFAULTS, **(self.params or {})}
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs or {},
            outputs=self.outputs or {"bcf_prefix": os.path.join(str(self.workdir), "NA12878_imputed")},
            params=p,
            threads=int(self.threads or 1),
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id") if self.params else None,
            ensure_output_dir_key=None,
        )
