#!/usr/bin/env python3
from __future__ import annotations
from typing import List
from pathlib import Path
import argparse
import os, sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
# Wave 패키지 내부 import (구조에 맞게 조정)
from wave.core.executor import SunGridExecutor
from wave.utils.task_utils import ensure_dir

# atomic task 들 (이미 만들어둔 Task들 가정)
from wave.tasks.picard.downsamplesam.main import PicardDownSampleSamTask
from wave.tasks.mosdepth.main import MosdepthMeanDepthTask

# -----------------------------
# argparse
# -----------------------------
def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Downsample BAM by depth × chromosome using Wave tasks",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--bam", required=True, help="Input BAM (sorted + dedup)")
    p.add_argument("--sample-id", required=True, help="Sample ID")
    p.add_argument(
        "--work-dir",
        required=True,
        help="Workflow root directory (e.g. /path/to/work/sample_root)",
    )
    p.add_argument(
        "--depths",
        default="0.5,1,2",
        help="Comma-separated target depths (e.g. '0.5,1,2')",
    )
    p.add_argument(
        "--chroms",
        default="chr1,chr2,chr3",
        help="Comma-separated chromosomes (e.g. 'chr1,chr2,chrX')",
    )
    p.add_argument(
        "--node",
        default="compute-0-0",
        help="SGE node name (qsub -l h=node)",
    )
    p.add_argument(
        "--threads-depth",
        type=int,
        default=2,
        help="Threads for samtools depth / mosdepth step",
    )
    p.add_argument(
        "--threads-downsample",
        type=int,
        default=4,
        help="Threads for Picard DownsampleSam step",
    )
    return p


# -----------------------------
# main workflow logic
# -----------------------------
def main():
    ap = build_argparser()
    args = ap.parse_args()

    bam = Path(args.bam).resolve()
    sample_id = args.sample_id
    work_root = Path(args.work_dir).resolve()

    depths: List[float] = [float(x) for x in args.depths.split(",") if x.strip()]
    chroms: List[str] = [c.strip() for c in args.chroms.split(",") if c.strip()]

    # 샘플 루트 디렉토리
    sample_root = ensure_dir(work_root / sample_id)
    print(f"[Wave] sample_root = {sample_root}")

    # 공용 executor (로그는 sample_root/log 아래)
    executor = SunGridExecutor(logdir=sample_root / "log")

    # 1) 전체 BAM 기준 baseline depth를 먼저 측정하는 task (선택사항)
    baseline_dir = ensure_dir(sample_root / "baseline_depth")
    baseline_depth_task = SamtoolsDepthTask(
        name=f"{sample_id}.baseline_depth",
        tool="samtools",
        func="depth",                 # 네가 Task에서 어떻게 쓰는지에 맞춰 조정
        threads=args.threads_depth,
        workdir=baseline_dir,
        inputs={
            "bam": str(bam),
            # 필요하다면 regions / bed 파일 추가
        },
        outputs={
            "depth": str(baseline_dir / f"{sample_id}.depth.txt"),
        },
        params={
            # samtools depth 옵션 커스터마이징 가능
        },
    )
    baseline_cmd = baseline_depth_task.to_sh()[0]
    baseline_script = executor.make_script(
        cmd=baseline_cmd,
        job_id=baseline_depth_task.name,
        workdir=str(baseline_depth_task.workdir),
        outputs=baseline_depth_task.outputs,
    )
    baseline_qid = executor.qsub_sh(
        node=args.node,
        script_path=str(baseline_script),
        threads=baseline_depth_task.threads,
        job_id=baseline_depth_task.name,
        random_jobid=False,
        hold_jid=None,
    )
    print(f"[Wave] qsub baseline depth -> {baseline_qid}")

    # 2) chromosome × depth 조합으로 DownsampleSam 태스크를 던짐
    all_qids = []

    for chrom in chroms:
        chrom_root = ensure_dir(sample_root / f"chrom_{chrom}")

        for depth in depths:
            tag = f"{chrom}.{depth}x"
            task_name = f"{sample_id}.downsample.{tag}"

            workdir = ensure_dir(chrom_root / f"downsample_{depth}x")

            # DownsampleSam Task 인스턴스 생성
            down_task = PicardDownsampleSamTask(
                name=task_name,
                tool="picard",
                func="downsamplesam",
                threads=args.threads_downsample,
                workdir=workdir,
                inputs={
                    "bam": str(bam),
                    # chrom별로 bam을 자른 버전을 쓰고 싶다면
                    # 미리 per-chrom bam 생성 Task를 하나 더 두고 여기 input 교체
                },
                outputs={
                    "bam": str(workdir / f"{sample_id}.{chrom}.{depth}x.bam"),
                    # 필요하다면 index, metrics 등도 outputs에 추가
                },
                params={
                    "sample_id": sample_id,
                    "target_depth": depth,
                    "chrom": chrom,
                    # Picard task 내부에서 P 계산을 BAM depth 기준으로 하거나,
                    # 여기서 P 직접 계산을 넣어도 됨
                },
            )

            down_cmd = down_task.to_sh()[0]
            down_script = executor.make_script(
                cmd=down_cmd,
                job_id=down_task.name,
                workdir=str(down_task.workdir),
                outputs=down_task.outputs,
            )

            # baseline depth 끝난 뒤에만 시작하도록 hold_jid=baseline_qid
            qid = executor.qsub_sh(
                node=args.node,
                script_path=str(down_script),
                threads=down_task.threads,
                job_id=down_task.name,
                random_jobid=False,
                hold_jid=baseline_qid,
            )
            all_qids.append(qid)
            print(f"[Wave] qsub {task_name} -> {qid} (hold_jid={baseline_qid})")

    print(f"[Wave] submitted {len(all_qids)} downsample jobs.")
    # 여기서 qid 리스트를 파일로 저장해두거나, 후속 워크플로우에서 읽어와서 후처리하는 것도 가능


if __name__ == "__main__":
    main()