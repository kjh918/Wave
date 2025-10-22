#!/usr/bin/env python3
import subprocess
from pathlib import Path
from typing import List, Optional

def build_samtools_stat_cmd(
        samtools: str,
        input_bam: str,
        output: str ,
        threads: int = 2,  # picard.jar inside sif
    ):
    """
    Build (and optionally run) a Picard DownsampleSam command.

    Returns:
        List[str] – the full command as list
    """
    
    # Add required arguments
    cmd += [
        f"{samtools}", 'stat',
        "-@", f"{thread}",
        f"{input_bam}", '>',
        f"{output}"
    ]

    if run_cmd:
        print("[CMD]", " ".join(cmd))
        subprocess.run(cmd, check=True)

    return cmd


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Samtools stat wrapper (build/run command)")
    parser.add_argument("-s", "--samtools", defaults='samtools', help="Path samtools")
    parser.add_argument("-i", "--input", required=True, help="Input BAM/SAM file")
    parser.add_argument("-o", "--output", required=True, help="Output dir state of BAM file")
    parser.add_argument("-t", "--thread", defaults=2, type=int, help="threads")
    parser.add_argument("-id", "--sample_id", defaults='-', type=str, help="threads")
    parser.add_argument("--run", action="store_true", help="Run command immediately")

    args = parser.parse_args()
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    if args.sample_id == '-':
        sample_id = str(Path(args.input).name)split('.bam','')
    
    output = str(Path(args.output) / f'{sample_id}.stat.txt')

    cmd = build_samtools_stat_cmd(
        samtools=args.samtools,
        input_bam=args.input,
        output=output,
        thread=thread
    )

    print("[CMD]", " ".join(cmd))
