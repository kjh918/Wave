from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional
import os, shlex

from src.tasks.task import Task
from src.tasks.task_registry import register_task
from src.tasks.utils import ensure_dir, to_sh_from_builder, singularity_exec_cmd, normalize_binds

@register_task("picard.downsamplesam_depth")
class PicardDownsampleByDepthTask(Task):
    TYPE = "picard.downsamplesam_depth"

    INPUTS = {
        "bam": {"type": "path", "required": True, "desc": "Input BAM (duplicate-marked)"},
    }
    OUTPUTS = {
        "dir": {"type": "dir", "required": False, "desc": "Output directory"},
    }
    DEFAULTS = {
        # target_mode: "fraction" (원본 대비) 또는 "absolute" (절대 배수, e.g. 0.5x)
        "target_mode": "fraction",
        "target": 0.5,              # fraction이면 0.5=절반, absolute면 0.5=0.5x
        "tolerance": 0.03,          # 허용 상대 오차 (±3%)
        "max_iters": 3,
        "java_bin": "java",
        "picard_jar": "/storage/apps/bin/picard.jar",
        "mosdepth_bin": "/storage/apps/mosdepth-0.3.5/mosdepth",
        "xmx_gb": 8,
        "seed": 42,
        "create_index": True,
        "image": None,
        "binds": None,
        "singularity_bin": "singularity",
    }

    def _mosdepth_block(self, bam: str, prefix: str, threads: int, mosdepth_bin: str) -> List[str]:
        # mean depth 계산 (summary.txt에서 길이 가중 평균)
        lines = [
            f'{mosdepth_bin} -t {threads} --no-per-base --fast-mode {shlex.quote(prefix)} {shlex.quote(bam)}',
            # chrom length bases mean  → 전체 가중 평균
            f"awk 'NR>1{{L+=$2; S+=$2*$4}} END{{if(L>0) printf(\"%.6f\\n\", S/L); else print 0}}' {shlex.quote(prefix)}.mosdepth.summary.txt > {shlex.quote(prefix)}.mean.txt",
        ]
        return lines

    def _read_backtick(self, path: str) -> str:
        # 유틸: 셀에서 $(cat file) 와 같은 효과를 주고 싶으면
        return f'$(cat {shlex.quote(path)};)'

    def _build_cmd(self, *, inputs: Dict[str, Any], outputs: Dict[str, Any],
                   params: Dict[str, Any], threads: int, workdir: str,
                   sample_id: Optional[str] = None) -> List[str]:
        bam = inputs["bam"]


        out_dir = ensure_dir(outputs.get("dir") or workdir)
        base = str(params.get("sample_id"))

        java_bin = str(params.get("java_bin", "java"))
        jar = str(params.get("picard_jar", "/storage/apps/bin/picard.jar"))
        mosdepth_bin = str(params.get("mosdepth_bin", "mosdepth"))
        xmx = int(params.get("xmx_gb", 8))
        seed = int(params.get("seed", 42))
        create_index = str(params.get("create_index", True)).lower()
        target_mode = str(params.get("target_mode", "fraction")).lower()
        target = float(params.get("target", 0.5))
        tol = float(params.get("tolerance", 0.03))
        iters = int(params.get("max_iters", 3))

        # 컨테이너 감싸기 여부
        image = params.get("image")
        binds = normalize_binds(params.get("binds"))
        sng = str(params.get("singularity_bin", "singularity"))

        cmds: List[str] = []

        # 1) baseline depth 측정
        baseline_prefix = os.path.join(out_dir, f"{base}.baseline")
        if os.path.isfile(f'{baseline_prefix}.mean.txt') == False:
            cmds += self._mosdepth_block(bam=bam, prefix=baseline_prefix, threads=threads, mosdepth_bin=mosdepth_bin)
        # bash 변수로 읽어오기
        cmds += [f'BASE_DEPTH=$({{ cat {shlex.quote(baseline_prefix)}.mean.txt;}})',
                 f'echo "[info] Baseline mean depth: $BASE_DEPTH x"']

        # 2) target depth 계산
        target_list = [0.1, 0.2, 0.5, 1, 1.5, 2]
        for target in target_list:
            if target_mode == "fraction":
                cmds += [f'TARGET_DEPTH=$(python3 - <<PY\nprint(round(float("{target}")/float("$BASE_DEPTH"),4))\nPY\n)\n']
            else:
                cmds += [f'TARGET_DEPTH="{target}"']

            cmds += [f'echo "[info] Target mean depth: $TARGET_DEPTH x";']

            loop = f'''
OUT_BAM={shlex.quote(os.path.join(out_dir, f"{base}.downsample.p"))}"$TARGET_DEPTH"_{target}x.bam
{' '.join(map(shlex.quote, [java_bin]))} -Xmx{xmx}g -jar {shlex.quote(jar)} DownsampleSam \
    I={shlex.quote(bam)} \
    O="$OUT_BAM" \
    P="$TARGET_DEPTH" \
    RANDOM_SEED={seed} \
    CREATE_INDEX={create_index}
            '''
            cmds += [loop]
            
        return ['\n'.join(cmds)]

    def to_sh(self) -> List[str]:
        p = {**self.DEFAULTS, **(self.params or {})}
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs or {},
            outputs=self.outputs or {"dir": str(self.workdir)},
            params=p,
            threads=int(self.threads or 4),
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id"),
            ensure_output_dir_key="dir",    
        )