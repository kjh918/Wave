from __future__ import annotations
from typing import Dict, Any, List, Sequence, Optional
import os, shlex

from src.tasks.task import Task
from src.tasks.task_registry import register_task
from src.tasks.utils import (
    ensure_dir,
    normalize_binds,
    singularity_exec_cmd,
    to_sh_from_builder,
)


@register_task("picard.downsamplesam")
class PicardDownSampleSamTask(Task):
    """
    Picard MergeBamAlignment 실행 태스크

    INPUTS:
      unmapped_bam: FastqToSam 결과 BAM
      aligned_bam:  BWA-MEM SAM/BAM (동일 sample)
      reference:    Reference FASTA

    OUTPUTS:
      bam:          결과 primary BAM (default: {workdir}/{sample_id}.primary.bam)
      dir:          결과 디렉토리 (optional)

    PARAMS:
      create_index: bool (default: True)
      max_indels:   int  (default: -1)
      clip_adapters: bool (default: False)
      primary_strategy: str ("MostDistant" 기본)
      retain_attr:  list[str] (default: ["XS"])
      expected_orientations: list[str] (default: ["FR", "RF"])
      java_bin:     str ("java")
      jar:          str ("/storage/apps/bin/picard.jar")
      xmx_gb:       int (16)
      parallel_gc_threads: int (14)
      image:        str | None
      binds:        list[str] | str | None
      singularity_bin: str ("singularity")
    """

    TYPE = "picard.downsamplesam"

    INPUTS: Dict[str, Any] = {
        "unmapped_bam": {"type": "path", "required": True, "desc": "FastqToSam output BAM"},
        "aligned_bam":  {"type": "path", "required": True, "desc": "BWA MEM SAM/BAM"},
        "reference":    {"type": "path", "required": True, "desc": "Reference FASTA"},
    }

    OUTPUTS: Dict[str, Any] = {
        "bam": {"type": "path", "required": False, "desc": "Output primary BAM"},
    }

    DEFAULTS: Dict[str, Any] = {
        "create_index": True,
        "max_indels": -1,
        "clip_adapters": False,
        "primary_strategy": "MostDistant",
        "retain_attr": ["XS"],
        "expected_orientations": ["FR", "RF"],
        "java_bin": "java",
        "jar": "/storage/apps/bin/picard.jar",
        "xmx_gb": 16,
        "parallel_gc_threads": 14,
        "image": "/storage/images/gatk-4.4.0.0.sif",
        "binds": ["/storage", "/data"],
        "singularity_bin": "singularity",
    }

    # ---- 내부 빌더
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

        unmapped_bam = inputs.get("unmapped_bam")
        aligned_bam  = inputs.get("aligned_bam")
        reference    = inputs.get("reference")
        if not all([unmapped_bam, aligned_bam, reference]):
            raise ValueError("[picard.mergebamalignment] Missing required inputs")

        out_dir = ensure_dir(outputs.get("dir") or workdir)
        out_bam = outputs.get("bam") or os.path.join(out_dir, f"{sample_id or 'sample'}.primary.bam")

        java_bin = str(params.get("java_bin", "java"))
        jar_path = str(params.get("jar", "/storage/apps/bin/picard.jar"))
        xmx_gb = int(params.get("xmx_gb", 16))
        pgc = int(params.get("parallel_gc_threads", 14))

        # core Picard arguments
        argv = [
            java_bin,
            f"-XX:ParallelGCThreads={pgc}",
            f"-Xmx{xmx_gb}g",
            "-jar", jar_path,
            "MergeBamAlignment",
            "--UNMAPPED_BAM", unmapped_bam,
            "--ALIGNED_BAM", aligned_bam,
            "--REFERENCE_SEQUENCE", reference,
            "--OUTPUT", out_bam,
            "--CREATE_INDEX", str(params.get("create_index", True)).lower(),
            "--MAX_INSERTIONS_OR_DELETIONS", str(params.get("max_indels", -1)),
            "--CLIP_ADAPTERS", str(params.get("clip_adapters", False)).lower(),
            "--PRIMARY_ALIGNMENT_STRATEGY", str(params.get("primary_strategy", "MostDistant")),
        ]

        # optional list-like parameters
        retain_attr = params.get("retain_attr", ["XS"])
        expected_orient = params.get("expected_orientations", ["FR", "RF"])
        for a in retain_attr:
            argv += ["--ATTRIBUTES_TO_RETAIN", str(a)]
        for o in expected_orient:
            argv += ["--EXPECTED_ORIENTATIONS", str(o)]

        # 컨테이너 래핑
        image = params.get("image")
        if image:
            cmd = singularity_exec_cmd(
                image=str(image),
                argv=argv,
                binds=normalize_binds(params.get("binds")),
                singularity_bin=str(params.get("singularity_bin", "singularity")),
            )
            return [" ".join(map(shlex.quote, cmd))]
        else:
            return [" ".join(map(shlex.quote, argv))]

    # ---- 메인 실행용 to_sh
    def to_sh(self) -> List[str]:
        p = {**self.DEFAULTS, **(self.params or {})}
        return to_sh_from_builder(
            builder=self._build_cmd,
            inputs=self.inputs or {},
            outputs=self.outputs or {"dir": str(self.workdir)},
            params=p,
            threads=int(self.threads or 1),
            workdir=str(self.workdir),
            sample_id=self.params.get("sample_id"),
            ensure_output_dir_key="dir",
        )