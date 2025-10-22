#!/usr/bin/env python3
import subprocess
from pathlib import Path
from typing import List, Optional

def build_picard_downsamplesam_cmd(
        input_bam: str,
        output_bam: str,
        probability: float,
        strategy: str = "ConstantMemory",           # [ConstantMemory | HighAccuracy | Chained]
        accuracy: Optional[float] = None,           # only used for HighAccuracy/Chained
        random_seed: Optional[int] = None,
        create_index: bool = True,
        validation_stringency: str = "LENIENT",     # [STRICT | LENIENT | SILENT]
        java_opts: Optional[List[str]] = None,      # e.g. ["-Xmx8g"]
        picard_jar: str = "picard.jar",             # Path to picard.jar
        run_cmd: bool = False,                      # if True, execute immediately
        singularity_sif: Optional[str] = None,      # path to singularity sif file
        bind_paths: Optional[List[str]] = None,     # singularity bind dirs
        picard_inside: str = "/picard/picard.jar",  # picard.jar inside sif
    ):
    """
    Build (and optionally run) a Picard DownsampleSam command.

    Returns:
        List[str] – the full command as list
    """
    if java_opts is None:
        java_opts = ["-Xmx8g"]
    if bind_paths is None:
        bind_paths = ["/storage", "/data"]

    # Base command selection (Singularity or local)
    if singularity_sif:
        cmd = ["singularity", "exec"]
        for b in bind_paths:
            cmd += ["-B", b]
        cmd += [singularity_sif, "java", *java_opts, "-jar", picard_inside, "DownsampleSam"]
    else:
        cmd = ["java", *java_opts, "-jar", picard_jar, "DownsampleSam"]

    # Add required arguments
    cmd += [
        f"I={input_bam}",
        f"O={output_bam}",
        f"P={probability}",
        f"STRATEGY={strategy}",
        f"CREATE_INDEX={'true' if create_index else 'false'}",
        f"VALIDATION_STRINGENCY={validation_stringency}"
    ]

    # Optional parameters
    if accuracy is not None:
        cmd += [f"ACCURACY={accuracy}"]
    if random_seed is not None:
        cmd += [f"RANDOM_SEED={random_seed}"]

    if run_cmd:
        print("[CMD]", " ".join(cmd))
        subprocess.run(cmd, check=True)

    return cmd


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Picard DownsampleSam wrapper (build/run command)")
    parser.add_argument("-i", "--input", required=True, help="Input BAM/SAM file")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file")
    parser.add_argument("-p", "--probability", required=True, type=float, help="Retain probability (0~1)")
    parser.add_argument("-s", "--strategy", default="ConstantMemory", choices=["ConstantMemory", "HighAccuracy", "Chained"])
    parser.add_argument("-a", "--accuracy", type=float, help="Accuracy (for HighAccuracy/Chained)")
    parser.add_argument("--random-seed", type=int, help="Random seed for reproducibility")
    parser.add_argument("--no-create-index", dest="create_index", action="store_false", help="Disable index creation")
    parser.add_argument("--validation-stringency", default="LENIENT", choices=["STRICT", "LENIENT", "SILENT"])
    parser.add_argument("--picard-jar", default="picard.jar", help="Path to picard.jar (if local)")
    parser.add_argument("--java-opts", nargs="*", default=["-Xmx8g"], help="Java options")
    parser.add_argument("--singularity-sif", help="Path to Singularity SIF for running Picard")
    parser.add_argument("--bind", nargs="*", default=["/storage", "/data"], help="Bind paths for Singularity")
    parser.add_argument("--picard-inside", default="/picard/picard.jar", help="Path to picard.jar inside SIF")
    parser.add_argument("--run", action="store_true", help="Run command immediately")

    args = parser.parse_args()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    cmd = build_picard_downsamplesam_cmd(
        input_bam=args.input,
        output_bam=args.output,
        probability=args.probability,
        strategy=args.strategy,
        accuracy=args.accuracy,
        random_seed=args.random_seed,
        create_index=getattr(args, "create_index", True),
        validation_stringency=args.validation_stringency,
        java_opts=args.java_opts,
        picard_jar=args.picard_jar,
        run_cmd=args.run,
        singularity_sif=args.singularity_sif,
        bind_paths=args.bind,
        picard_inside=args.picard_inside
    )

    print("[CMD]", " ".join(cmd))
