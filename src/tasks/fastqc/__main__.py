# tasks/fastqc/__main__.py
import argparse
from ._func import build_fastqc_cmd
from utils.sh_writer import write_script_from_cmds

ap = argparse.ArgumentParser()
ap.add_argument("--SeqID", required=True)
ap.add_argument("--RawFastqDir", required=True)
ap.add_argument("--qcResDir", required=True)
ap.add_argument("--threads", type=int, default=2)
ap.add_argument("--emit", required=True)  # .sh 저장 경로
ap.add_argument("--run", action="store_true")
args = ap.parse_args()

cmds = build_fastqc_cmd(qcResDir=args.qcResDir, RawFastqDir=args.RawFastqDir,
                        SeqID=args.SeqID, threads=args.threads)
script = write_script_from_cmds(cmds, args.emit, set_x=True)
if args.run:
    import subprocess
    subprocess.run(["bash", str(script)], check=True)