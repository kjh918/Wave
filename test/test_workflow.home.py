
import os
import sys
import argparse
import subprocess
from glob import glob
from typing import Dict, Any, List, Optional
import pyaml 
from pathlib import Path 

ROOT = Path(__file__).resolve().parents[1]    # Wave/
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
from wave.core.workflow import Workflow 


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="run",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--config", dest="config", help="config.yaml")
    return p

def main():
    parser = build_argparser()
    args = parser.parse_args()

    wf = Workflow(args.config)

    # 3️⃣ 샘플 자동 탐색
    samples = wf.discover_samples()
    print(f"\n[WAVE] Found {len(samples)} samples:")

    workflow_dict = wf.build()

    # wf.run()

    




if __name__ == '__main__':

    main()

