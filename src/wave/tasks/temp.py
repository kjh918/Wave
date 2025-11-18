import os
from pathlib import Path
from glob import glob 

ROOT = Path("src/wave/tasks")  # tasks 패키지 경로


def fix_main_imports(main_file: Path):
    """main.py 파일 내용을 읽어서 import 부분 자동 치환"""

    text = main_file.read_text()

    # --- 기존 패턴들 제거/대체
    replacements = {
        # Task
        "from src.wave.core.task import Task": 
        "from wave.core.task import Task",

        # TaskRegistry
        "from src.wave.core.task_registry import register_task":
        "from wave.core.task_registry import register_task",

        # Utils
        "from src.wave.utils.task_utils import":
        "from wave.utils.task_utils import",
    }

    new_text = text
    for old, new in replacements.items():
        new_text = new_text.replace(old, new)

    # 저장
    if new_text != text:
        main_file.write_text(new_text)
        print(f"[FIXED] {main_file}")
    else:
        print(f"[SKIPPED] {main_file}")


def scan_and_fix(root: Path):
    """재귀적으로 tasks 하위의 모든 main.py 탐색"""
    for path in root.rglob("main.py"):
        
        with open(f'{os.path.dirname(path)}/__init__.py','w') as handle:
            pass
        with open(f'{os.path.dirname(os.path.dirname(path))}/__init__.py','w') as handle:
            pass

        fix_main_imports(path)


if __name__ == "__main__":
    print(f"🔍 Scanning for main.py in {ROOT}")
    # scan_and_fix(Path('/storage/home/jhkim/scripts/Task/Wave/src/wave/tasks'))
    scan_and_fix(Path('/Users/kimjihoon/Downloads/GdriveBackup/Projects/Wave/src/wave/tasks'))
    print("\n✨ Import fix completed.")