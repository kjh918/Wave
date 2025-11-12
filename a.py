from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional
import os, shlex

from src.tasks.task import Task
from src.tasks.task_registry import register_task
from src.tasks.utils import ensure_dir, to_sh_from_builder, singularity_exec_cmd, normalize_binds

@register_task("picard.downsample_by_depth")
class PicardDownsampleByDepthTask(Task):
    TYPE = "picard.downsample_by_depth"

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
        "mosdepth_bin": "mosdepth",
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
            f'{mosdepth_bin} -t {threads} --no-per-base --fast {shlex.quote(prefix)} {shlex.quote(bam)}',
            # chrom length bases mean  → 전체 가중 평균
            f"awk 'NR>1{{L+=$2; S+=$2*$4}} END{{if(L>0) printf(\"%.6f\\n\", S/L); else print 0}}' {shlex.quote(prefix)}.mosdepth.summary.txt > {shlex.quote(prefix)}.mean.txt",
        ]
        return lines

    def _read_backtick(self, path: str) -> str:
        # 유틸: 셀에서 $(cat file) 와 같은 효과를 주고 싶으면
        return f'$(cat {shlex.quote(path)})'

    def _build_cmd(self, *, inputs: Dict[str, Any], outputs: Dict[str, Any],
                   params: Dict[str, Any], threads: int, workdir: str,
                   sample_id: Optional[str] = None) -> List[str]:
        bam = inputs["bam"]
        out_dir = ensure_dir(outputs.get("dir") or workdir)
        base = sample_id or os.path.splitext(os.path.basename(bam))[0]

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
        cmds += self._mosdepth_block(bam=bam, prefix=baseline_prefix, threads=threads, mosdepth_bin=mosdepth_bin)
        # bash 변수로 읽어오기
        cmds += [f'BASE_DEPTH=$({{ cat {shlex.quote(baseline_prefix)}.mean.txt }})',
                 f'echo "[info] Baseline mean depth: $BASE_DEPTH x"']

        # 2) target depth 계산
        if target_mode == "fraction":
            cmds += [f'TARGET_DEPTH=$(python3 - <<PY\nprint(float("{target}")*float("$BASE_DEPTH"))\nPY)']
        else:
            # absolute (e.g. 0.5x)
            cmds += [f'TARGET_DEPTH="{target}"']

        cmds += [f'echo "[info] Target mean depth: $TARGET_DEPTH x"']

        # 3) 첫 P 설정
        if target_mode == "fraction":
            # 원본 대비 비율 → P = target
            cmds += [f'CUR_P="{target}"']
        else:
            # 절대 배수 → P = min(1, target / base)
            cmds += [f'CUR_P=$(python3 - <<PY\nb=float("$BASE_DEPTH"); t=float("$TARGET_DEPTH"); print(1.0 if b<=0 else min(1.0,max(0.0,t/b)))\nPY)']

        cmds += ['ATTEMPT=1']

        # 4) 반복 보정 루프 (bash)
        #   원본에서 매번 새로 샘플링
        loop = f'''
while true; do
  echo "[info] Attempt $ATTEMPT: P=$CUR_P"
  OUT_BAM={shlex.quote(os.path.join(out_dir, f"{base}.downsample.p"))}$CUR_P.bam

  {' '.join(map(shlex.quote, [java_bin]))} -Xmx{xmx}g -jar {shlex.quote(jar)} DownsampleSam \
      I={shlex.quote(bam)} \
      O="$OUT_BAM" \
      P="$CUR_P" \
      RANDOM_SEED={seed} \
      CREATE_INDEX={create_index}

  # 측정
  VERIFY_PREFIX={shlex.quote(os.path.join(out_dir, f"{base}.p"))}$CUR_P.verify
  {shlex.quote(mosdepth_bin)} -t {threads} --no-per-base --fast "$VERIFY_PREFIX" "$OUT_BAM" >/dev/null 2>&1
  MEAN=$(awk 'NR>1{{L+=$2; S+=$2*$4}} END{{if(L>0) printf("%.6f\\n", S/L); else print 0}}' "$VERIFY_PREFIX.mosdepth.summary.txt")
  echo "[info] Measured mean depth: $MEAN x"

  # 수렴 체크
  RELERR=$(python3 - <<PY
td=float("$TARGET_DEPTH"); md=float("$MEAN")
print(0.0 if td<=0 else abs(md-td)/td)
PY
)
  echo "[info] Relative error: $RELERR"

  python3 - <<'PY'
tol=float("{tol}")
rel=float(os.environ.get("RELERR","1"))
import sys
sys.exit(0 if rel <= tol else 1)
PY
  OK=$?
  if [ "$OK" -eq 0 ]; then
    echo "[ok] Converged within tolerance ±{tol}"
    break
  fi

  # 보정 (원본에서 재샘플링하도록 P 재계산)
  CUR_P=$(python3 - <<PY
cur=float("$CUR_P")
md=float("$MEAN"); td=float("$TARGET_DEPTH")
adj = cur * (td / (md if md>1e-9 else td))
print(min(1.0, max(0.0, adj)))
PY
)
  ATTEMPT=$((ATTEMPT+1))
  if [ "$ATTEMPT" -gt {iters} ]; then
    echo "[warn] Max iterations reached. Keeping last result."
    break
  fi
done
'''
        cmds += [loop]

        # 컨테이너 래핑 필요시 전체를 singularity에서 돌릴 수도 있지만,
        # 여기선 Picard만 컨테이너 감쌌을 때를 고려해 각 라인 그대로 반환.
        return cmds

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