# _aligned_bam.py
import subprocess
from pathlib import Path

def _alinged_bam(
        singularity: Path,
        singularity_image: str,
        bind_paths: list[str],
        bwa_index_prefix: str,
        seq_id: str,
        forward_read: str,
        reverse_read: str,
        out_sam_dir: str,
        threads: int,
        read_group_id: str,
        read_group_platform: str,
        read_group_library: str,
        read_group_center: str,
        bwa_extra: list[str] | None = None,
    ) -> tuple[str, str | None]:
    """
    Run bwa mem inside a Singularity container to produce SAM.
    Returns path to <seq_id>.bwa.mem.sam
    """
    out_sam = Path(out_sam_dir) / f"{seq_id}.bwa.mem.sam"
    out_sam.parent.mkdir(parents=True, exist_ok=True)

    # Bind string like: -B /storage,/data
    bind_arg = []
    if bind_paths:
        bind_arg = ["-B", ",".join(bind_paths)]

    rg_str = f"@RG\\tID:{read_group_id}\\tPL:{read_group_platform}\\tLB:{read_group_library}\\tSM:{seq_id}\\tCN:{read_group_center}"

    bwa_cmd = [
        f"{str(singularity)}", "exec", *bind_arg, singularity_image,
        "bwa", "mem",
        "-M",
        "-t", str(threads),
        "-Y",
        "-L", "50,50",
        "-R", f'"{rg_str}"',
        bwa_index_prefix,
        str(forward_read), str(reverse_read),
        ">", str(out_sam)
    ]

    if bwa_extra:
        # insert extra args before inputs if needed; here append for simplicity
        bwa_cmd[ bwa_cmd.index("-R")+2 : bwa_cmd.index(bwa_index_prefix) ] += bwa_extra

    command = ' '.join(bwa_cmd)
    return command

def _unalinged_bam(
        java: Path,  
        seq_id: str,
        picard_jar: str,
        forward_read: str,
        reverse_read: str,
        bam_dir: str,
        read_group_id: str,
        read_group_platform: str,
        read_group_library: str,
        read_group_center: str,
        tmp_dir: str | None = None,
        java_xmx: str = "16384m",
        parallel_gc_threads: int = 14,
    ) -> tuple[str, str | None]:
    """                     
        Picard FastqToSam to create unmapped BAM from paired-end FASTQ.
        Returns [cmd]
    """

    out_bam = Path(bam_dir) / f"{seq_id}.fastqtosam.bam"
    out_bam.parent.mkdir(parents=True, exist_ok=True)

    if tmp_dir:
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    cmd = [
        str(java),
        f"-XX:ParallelGCThreads={parallel_gc_threads}",
        f"-Xmx{java_xmx}",
        "-jar", picard_jar, "FastqToSam",
        "--FASTQ", str(forward_read),
        "--FASTQ2", str(reverse_read),
        "--SAMPLE_NAME", seq_id,
        "--OUTPUT", str(out_bam),
        "--READ_GROUP_NAME", read_group_id,
        "--PLATFORM", read_group_platform,
        "--LIBRARY_NAME", read_group_library,
        "--SEQUENCING_CENTER", read_group_center,
    ]
    if tmp_dir:
        cmd += ["--TMP_DIR", str(tmp_dir)]
    command = ' '.join(cmd)
    return command

def _merge_bam_alignment(
        java: Path,
        picard_jar: str,
        unmapped_bam: str,
        aligned_sam: str,
        reference_fasta: str,
        out_bam_dir: str,
        seq_id: str,
        tmp_dir: str | None = None,
        java_xmx: str = "16384m",
        parallel_gc_threads: int = 14,
    ) -> tuple[str, str | None]:
    """
    Run Picard MergeBamAlignment to create primary BAM (+ BAI).
    Returns (primary_bam, bai_path_or_None)
    """
    primary_bam = Path(out_bam_dir) / f"{seq_id}.primary.bam"
    primary_bam.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(java),
        f"-XX:ParallelGCThreads={parallel_gc_threads}",
        f"-Xmx{java_xmx}",
        "-jar", picard_jar, "MergeBamAlignment",
        "--UNMAPPED_BAM", str(unmapped_bam),
        "--ALIGNED_BAM", str(aligned_sam),
        "--REFERENCE_SEQUENCE", str(reference_fasta),
        "--OUTPUT", str(primary_bam),
        "--CREATE_INDEX", "true",
        "--MAX_INSERTIONS_OR_DELETIONS", "-1",
        "--CLIP_ADAPTERS", "false",
        "--PRIMARY_ALIGNMENT_STRATEGY", "MostDistant",
        "--ATTRIBUTES_TO_RETAIN", "XS",
        "--EXPECTED_ORIENTATIONS", "FR",
        "--EXPECTED_ORIENTATIONS", "RF",
    ]
    if tmp_dir:
        Path(tmp_dir).mkdir(parents=True, exist_ok=True)
        cmd += ["--TMP_DIR", str(tmp_dir)]
    
    command = ' '.join(cmd)
    return command

def _align_bwa_only(
    singularity: Path,
    samtools: Path,
    singularity_image: str,
    bind_paths: list[str],
    bwa_index_prefix: str,
    forward_read: str,
    reverse_read: str,
    out_bam_sorted: str,
    threads: int,
    read_group_id: str,
    read_group_platform: str,
    read_group_library: str,
    read_group_center: str,
    sample_name: str,
    bwa_extra: list[str] | None = None,
) -> list[list[str]]:
    """
    Build a *pipeline*:
      singularity exec bwa mem ... |
      samtools view -b -@ threads |
      samtools sort -@ threads -o out_bam_sorted
    Returns a list of command segments (for Popen pipeline).
    """
    Path(out_bam_sorted).parent.mkdir(parents=True, exist_ok=True)

    bind_arg = []
    if bind_paths:
        bind_arg = ["-B", ",".join(bind_paths)]

    rg_str = (
        f"@RG\\tID:{read_group_id}\\tPL:{read_group_platform}"
        f"\\tLB:{read_group_library}\\tSM:{sample_name}\\tCN:{read_group_center}"
    )

    primary_bam = Path(out_bam_dir) / f"{seq_id}.primary.bam"
    primary_bam.parent.mkdir(parents=True, exist_ok=True)

    bwa_cmd = [
        f"{singularity}", "exec", *bind_arg, singularity_image,
        "bwa", "mem",
        "-M",
        "-t", str(threads),
        "-Y",
        "-L", "50,50",
        "-R", f'"{rg_str}"',
        bwa_index_prefix,
        str(forward_read),
        str(reverse_read),
    ]
    if bwa_extra:
        bwa_cmd += bwa_extra

    cmd += [f"| {samtools}", "view", "-bS", "-o", str(primary_bam)]

    return cmd


def build_cmd_sort_and_index(
        bam_in: str,
        bam_out: str,
        threads: int = 8,
    ) -> list[list[str]]:
    """
    Build samtools sort + index commands
    Returns a list of command lists (for sequential execution)
    """
    Path(bam_out).parent.mkdir(parents=True, exist_ok=True)
    cmd_sort = ["samtools", "sort", "-@", str(threads), "-o", bam_out, bam_in]
    cmd_index = ["samtools", "index", "-b", "-@", str(threads), bam_out]
    
    command = ' '.join(cmd_sort) + '\n' + ' '.join(cmd_index)
    return command

def build_cmd_mark_duplicates(
        gatk_image: str,
        bam_in: str,
        bam_out: str,
        metrics_out: str,
        threads: int = 8,
        tmp_dir: str | None = None,
        java_xmx: str = "16384m",
        parallel_gc_threads: int = 14,
    ) -> list[str]:
    """
    Build GATK4 MarkDuplicates command (via Singularity)
    """
    Path(bam_out).parent.mkdir(parents=True, exist_ok=True)
    Path(metrics_out).parent.mkdir(parents=True, exist_ok=True)
    bind_arg = ["-B", "/storage,/data"]

    cmd = [
        "singularity", "exec", *bind_arg, gatk_image,
        "gatk", "MarkDuplicates",
        "--java-options", f"'-XX:ParallelGCThreads={parallel_gc_threads} -Xmx{java_xmx}'",
        "--INPUT", bam_in,
        "--OUTPUT", bam_out,
        "--METRICS_FILE", metrics_out,
        "--CREATE_INDEX", "true",
        "--REMOVE_SEQUENCING_DUPLICATES", "true",
    ]
    if tmp_dir:
        cmd += ["--TMP_DIR", tmp_dir]
    command = ' '.join(cmd)
    return command


def build_cmd_local_realignment(
        gatk3_image: str,
        reference_fasta: str,
        bam_in: str,
        known_indel1: str,
        known_indel2: str,
        target_intervals: str,
        bam_out: str,
        singularity: str,
        threads: int = 8,
        java_xmx: str = "16384m",
    ) -> list[list[str]]:
    """
    Build GATK3 RealignerTargetCreator + IndelRealigner commands
    Returns [cmd_target_creator, cmd_indel_realigner]
    """
    Path(bam_out).parent.mkdir(parents=True, exist_ok=True)
    Path(target_intervals).parent.mkdir(parents=True, exist_ok=True)
    bind_arg = ["-B", "/storage,/data"]

    cmd_target = [
        f"{str(singularity)}", "exec", *bind_arg, gatk3_image,
        "java", f"-Xmx{java_xmx}", "-jar", "/storage/apps/gatk3/3.8-1/GenomeAnalysisTK.jar",
        "-T", "RealignerTargetCreator",
        "-R", reference_fasta,
        "-I", bam_in,
        "-o", target_intervals,
        "-known", known_indel1,
        "-known", known_indel2,
        "-nt", str(threads),
    ]
    cmd_realign = [
        f"{str(singularity)}", "exec", *bind_arg, gatk3_image,
        "java", f"-Xmx{java_xmx}", "-jar", "/storage/apps/gatk3/3.8-1/GenomeAnalysisTK.jar",
        "-T", "IndelRealigner",
        "-R", reference_fasta,
        "-targetIntervals", target_intervals,
        "-known", known_indel1,
        "-known", known_indel2,
        "-I", bam_in,
        "-o", bam_out,
    ]
    command = ' '.join(cmd_target) + '\n' + ' '.join(cmd_realign)
    return command

    
def build_cmd_bqsr(
        gatk_image: str,
        input_bam: str,
        output_bam: str,
        reference_fasta: str,
        output_table: str,
        known_sites: list[str],
        singularity: str,
        threads: int = 8,
        java_xmx: str = "16384m",
        parallel_gc_threads: int = 14,
        bind_paths:list = ["/storage", "/data"]
    ) -> list[str]:
    """
    Build GATK BaseRecalibrator command (creates recal.table)
    """
    Path(output_table).parent.mkdir(parents=True, exist_ok=True)
    bind_paths = bind_paths or ["/storage", "/data"]
    bind_arg = ["-B", ",".join(bind_paths)]

    cmd = [
        f"{str(singularity)}", "exec", *bind_arg, gatk_image,
        "gatk", "BaseRecalibrator",
        "--java-options", f"'-XX:ParallelGCThreads={parallel_gc_threads} -Xmx{java_xmx}'",
        "--input", input_bam,
        "--reference", reference_fasta,
        "--output", output_table,
    ]
    for site in known_sites:
        cmd += ["--known-sites", site]

    add_cmd = [
        f"{str(singularity)}", "exec", *bind_arg, gatk_image,
        "gatk", "ApplyBQSR",
        "--java-options", f"'-XX:ParallelGCThreads={parallel_gc_threads} -Xmx{java_xmx}'",
        "--input", input_bam,
        "--bqsr-recal-file", output_table,
        "--output", output_bam,
    ]
    command = ' '.join(cmd) + '\n' +' '.join(add_cmd)
    return command


def build_cmd_apply_bqsr(
        gatk_image: str,
        input_bam: str,
        recal_table: str,
        output_bam: str,
        singularity: str,
        threads: int = 8,
        java_xmx: str = "16384m",
        parallel_gc_threads: int = 14,
    ) -> list[str]:
    """
    Build GATK ApplyBQSR command (produces recalibrated BAM)
    """
    Path(output_bam).parent.mkdir(parents=True, exist_ok=True)
    bind_paths = bind_paths or ["/storage", "/data"]
    bind_arg = ["-B", ",".join(bind_paths)]

    cmd = [
        f"{str(singularity)}", "exec", *bind_arg, gatk_image,
        "gatk", "ApplyBQSR",
        "--java-options", f"'-XX:ParallelGCThreads={parallel_gc_threads} -Xmx{java_xmx}'",
        "--input", input_bam,
        "--bqsr-recal-file", recal_table,
        "--output", output_bam,
    ]
    command = ' '.join(cmd)
    return command


def build_cmd_filter_bam(
        samtools: str, 
        input_bam: str,
        output_bam: str,
        threads: int = 8,
    ) -> list[list[str]]:
    """
    Build samtools filter commands:
      1) samtools view with -e expressions
      2) samtools index
    """
    Path(output_bam).parent.mkdir(parents=True, exist_ok=True)

    cmd_filter = [
        f"{samtools}", "view",
        "-b", "-h",
        "-q", "20",               # MAPQ > 20
        "-f", "0x2",              # proper paired only
        "-F", "0x100",            # remove secondary
        "-F", "0x4",              # remove unmapped
        "-e", '"sclen < 20"',       # softclip length < 20
        "-e", '"[NM] < 12"',        # mismatch < 12
        "--threads", str(threads),
        input_bam,
        "-o", output_bam,
    ]
    cmd_index = [
        "samtools", "index",
        "--threads", str(threads),
        output_bam,
    ]
    
    command = ' '.join(cmd_filter) + '\n' + ' '.join(cmd_index)
    return command

def build_cmd_split_bam_by_chrom(
        samtools: str, 
        input_bam: str,
        chrom: str | list[str], 
        threads: int = 2
    ) -> list[list[str]]:

    # Handle both str and list input
    chrom_list = [chrom] if isinstance(chrom, str) else chrom
    cmds = ''
    for c in chrom_list:
        out_bam = input_bam.replace('.bam',f'.{c}.bam')
        cmd = [
            "samtools", "view",
            "-b", "-@", str(threads),
            input_bam, c,
            "-o", str(out_bam),
        ]
        cmds += ' '.join(cmd) + '\n'

    return cmds

    
    def build_cmd_merge_bam(
            samtools: str,
            bam_list: list[str],
            output_bam: str,
            threads: int = 8,
        ) -> list[str]:
        
        Path(output_bam).parent.mkdir(parents=True, exist_ok=True)

        cmd = ["samtools", "merge", "-@", str(threads), output_bam] + bam_list
