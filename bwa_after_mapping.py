# aftermapping.py
import subprocess
import logging
from pathlib import Path

from ._func import (
    build_cmd_sort_and_index,
    build_cmd_mark_duplicates,
    build_cmd_local_realignment,
    build_cmd_bqsr,
    build_cmd_filter_bam,
    build_cmd_split_bam_by_chrom
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


class AfterMapping:
    """
    Post-alignment steps:
      1) samtools sort + index
      2) GATK4 MarkDuplicates
      3) (optional) GATK3 Local Realignment (RealignerTargetCreator + IndelRealigner)
    실행은 단계별 run_* 메서드 또는 __call__ 로 일괄 수행.
    """

    def __init__(
        self,
        seq_id: str,
        bam_dir: str,
        qc_dir: str,
        reference_fasta: str,
        # 입력 BAM (기본: <bam_dir>/<seq_id>.primary.bam)
        primary_bam: str | None = None,

        # 리소스 & 공통 설정
        threads: int = 8,
        tmp_dir: str | None = None,
        java_xmx: str = "16384m",
        parallel_gc_threads: int = 14,
        bind_paths: list[str] | None = None,

        # 도구 경로 (Singularity 이미지)
        gatk4_image: str = "/storage/images/gatk-4.4.0.0.sif",
        gatk3_image: str = "/storage/images/gatk-3.8-1.sif",

        # Local realignment 옵션
        run_local_realignment: bool = True,
        known_indel1: str | None = None,
        known_indel2: str | None = None,
    ):
        self.seq_id = seq_id
        self.bam_dir = Path(bam_dir)
        self.qc_dir = Path(qc_dir)
        self.reference_fasta = reference_fasta
        self.threads = threads
        self.tmp_dir = tmp_dir
        self.java_xmx = java_xmx
        self.parallel_gc_threads = parallel_gc_threads
        self.bind_paths = bind_paths or ["/storage", "/data"]
        self.gatk4_image = gatk4_image
        self.gatk3_image = gatk3_image
        self.singularity = Path("/storage/apps/singularity/bin/singularity")

        # 입력 기본값
        self.primary_bam = primary_bam or str(self.bam_dir / f"{self.seq_id}.primary.bam")

        # 출력 경로들
        self.sorted_bam = str(self.bam_dir / f"{self.seq_id}.sorted.bam")
        self.sorted_bai = self.sorted_bam + ".bai"
        self.dedup_bam = str(self.bam_dir / f"{self.seq_id}.sorted.dedup.bam")
        self.metrics = str(self.qc_dir / f"{self.seq_id}.mark.duplicates.metrics.txt")
        self.target_intervals = str(self.qc_dir / f"{self.seq_id}.realign.target.intervals")
        self.realigned_bam = str(self.bam_dir / f"{self.seq_id}.sorted.dedup.realign.bam")

        # realignment 제어
        self.run_local_realignment = run_local_realignment
        self.known_snp = '/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.dbsnp138.vcf.gz'
        self.known_indel1 = '/storage/references_and_index/hg38/vcf/Homo_sapiens_assembly38.known_indels.vcf.gz'
        self.known_indel2 = '/storage/references_and_index/hg38/vcf/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'

        # 디렉토리 준비
        self.bam_dir.mkdir(parents=True, exist_ok=True)
        self.qc_dir.mkdir(parents=True, exist_ok=True)

    # --------------------
    # 내부 실행자
    # --------------------
    def _run(self, cmd: list[str]):
        logging.info("CMD: %s", " ".join(cmd))
        subprocess.run(cmd, check=True)

    # --------------------
    # 1) SORT + INDEX
    # --------------------
    def build_sort_index_cmds(self):
        cmd = build_cmd_sort_and_index(
            bam_in=self.primary_bam,
            bam_out=self.sorted_bam,
            threads=self.threads,
        )
        return [
            {'prefix':self.seq_id, 'thread': self.threads, 'cmd':cmd }
            ]

    def run_sort_index(self) -> tuple[str, str]:
        cmds = self.build_sort_index_cmds()
        for c in cmds:
            self._run(c)
        return self.sorted_bam, self.sorted_bai

    # --------------------
    # 2) MARK DUPLICATES
    # --------------------
    def build_mark_duplicates_cmd(self):
        cmd = build_cmd_mark_duplicates(
            gatk_image=self.gatk4_image,
            bam_in=self.sorted_bam,
            bam_out=self.dedup_bam,
            metrics_out=self.metrics,
            tmp_dir=self.tmp_dir,
            java_xmx=self.java_xmx,
            parallel_gc_threads=self.parallel_gc_threads
        )
        return [
            {'prefix':self.seq_id, 'thread': self.parallel_gc_threads, 'cmd':cmd }
            ]

    def run_mark_duplicates(self) -> tuple[str, str]:
        cmd = self.build_mark_duplicates_cmd()
        self._run(cmd)
        # MarkDuplicates에서 --CREATE_INDEX true 사용 → dedup_bam.bai 생성
        return self.dedup_bam, self.dedup_bam + ".bai"

    # --------------------
    # 3) LOCAL REALIGNMENT (옵션)
    # --------------------
    def build_local_realign_cmds(self):
        if not (self.known_indel1 and self.known_indel2):
            raise ValueError("known_indel1/known_indel2 경로가 필요합니다.")
        cmd = build_cmd_local_realignment(
            gatk3_image=self.gatk3_image,
            reference_fasta=self.reference_fasta,
            bam_in=self.dedup_bam,
            known_indel1=self.known_indel1,
            known_indel2=self.known_indel2,
            target_intervals=self.target_intervals,
            bam_out=self.realigned_bam,
            singularity=self.singularity,
            threads=self.threads,
            java_xmx=self.java_xmx,
        )
        return [
            {'prefix':self.seq_id, 'thread': self.parallel_gc_threads, 'cmd':cmd }
            ]


    def run_local_realign(self) -> str | None:
        if not self.run_local_realignment:
            logging.info("Skip local realignment (run_local_realignment=False)")
            return None
        cmds = self.build_local_realign_cmds()
        for c in cmds:
            self._run(c)
        return self.realigned_bam
    
    def build_bqsr_cmds(self):
        """
        Build GATK4 BQSR (BaseRecalibrator + ApplyBQSR)
        Returns [cmd_base_recal, cmd_apply_bqsr]
        """
        if not (self.known_indel1 and self.known_indel2):
            raise ValueError("known indel file paths are required for BQSR.")
        if not hasattr(self, "known_snp") and not getattr(self, "known_snp", None):
            raise ValueError("known_snp file path required for BQSR (set via known_snp=).")

        recal_table = str(self.qc_dir / f"{self.seq_id}.recal.table.txt")
        recal_bam = str(self.bam_dir / f"{self.seq_id}.recal.bam")
        
        cmd_recal = build_cmd_bqsr(
            gatk_image=self.gatk4_image,
            input_bam=self.realigned_bam,
            output_bam=recal_bam,
            reference_fasta=self.reference_fasta,
            output_table=recal_table,
            singularity=self.singularity,
            known_sites=[self.known_snp, self.known_indel1, self.known_indel2],
            threads=self.threads,
            java_xmx=self.java_xmx,
            parallel_gc_threads=self.parallel_gc_threads
        )
        return [
            {'prefix':self.seq_id, 'thread': self.parallel_gc_threads, 'cmd':cmd_recal }
            ]
    
    def build_cmds_filter_bam(self):

        recal_bam = str(self.bam_dir / f"{self.seq_id}.recal.bam")
        filtered_bam = str(self.bam_dir / f"{self.seq_id}.recal.filtered.bam").replace("merged_bam","filtered_bam")
        
        cmd_filter = build_cmd_filter_bam(
            samtools='samtools',
            input_bam = recal_bam,
            output_bam = filtered_bam,
            threads = 8
        )
        return [
            {'prefix':self.seq_id, 'thread': 8, 'cmd':cmd_filter }
            ]

    def build_cmds_split_bam_by_chrom(self):

        filtered_bam = str(self.bam_dir / f"{self.seq_id}.recal.filtered.bam").replace("merged_bam","filtered_bam")
        
        chrom_list =  [f'chr{i}' for i in range(1,23)] + ['chrX','chrY']
        cmd_spliter = build_cmd_split_bam_by_chrom(
            samtools='samtools',
            input_bam = filtered_bam,
            chrom = chrom_list, 
            threads = 4
        )
        return [
            {'prefix':self.seq_id, 'thread': 4, 'cmd':cmd_spliter }
            ]
    # --------------------
    # 일괄 실행
    # --------------------
    def __call__(self) -> dict:
        """
        전체 후처리 실행.
        Returns dict with output paths.
        """
        logging.info("AfterMapping started for %s", self.seq_id)
        sorted_bam, sorted_bai = self.run_sort_index()
        dedup_bam, dedup_bai = self.run_mark_duplicates()
        realigned_bam = self.run_local_realign()

        outputs = {
            "sorted_bam": sorted_bam,
            "sorted_bai": sorted_bai,
            "dedup_bam": dedup_bam,
            "dedup_bai": dedup_bai,
            "realigned_bam": realigned_bam,  # None if skipped
            "target_intervals": self.target_intervals if self.run_local_realignment else None,
            "metrics": self.metrics,
        }
        logging.info("AfterMapping finished for %s", self.seq_id)
        return outputs
