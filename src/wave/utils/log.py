import time
import sys
from pathlib import Path
from functools import wraps
from contextlib import contextmanager
from io import StringIO


class Logger:
    def __init__(self, run_id: str = None, logdir: Path = Path("./logs")):
        self.run_id = run_id or time.strftime("%Y%m%d_%H%M%S")
        self.logdir = Path(logdir)
        self.session_dir = self.logdir
        self.session_dir.mkdir(parents=True, exist_ok=True)

        self.log_file = self.session_dir / f"log.{self.run_id}.txt"

    # -----------------------------
    #  기본 파일 기록 메서드
    # -----------------------------
    def write(self, text: str):
        with open(self.log_file, "a", encoding="utf-8") as f:
            f.write(text + "\n")

    # -----------------------------
    #  함수 데코레이터 방식
    # -----------------------------
    def log_function(self, func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            self.write(f"[START] {func.__name__} at {time.strftime('%H:%M:%S')}")

            # 함수 내부 print 캡처
            old_stdout = sys.stdout
            sys.stdout = captured = StringIO()

            try:
                result = func(*args, **kwargs)
            finally:
                # 캡처된 출력 복원 + 기록
                sys.stdout = old_stdout
                out = captured.getvalue()
                if out.strip():
                    self.write(f"[PRINT OUTPUT]\n{out.strip()}")

            self.write(f"[END] {func.__name__} at {time.strftime('%H:%M:%S')}")
            self.write("-" * 50)
            return result
        return wrapper

    # -----------------------------
    #  Context Manager 방식
    # -----------------------------
    @contextmanager
    def log_block(self, name="block"):
        self.write(f"[START BLOCK] {name} - {time.strftime('%H:%M:%S')}")

        old_stdout = sys.stdout
        sys.stdout = captured = StringIO()

        try:
            yield
        finally:
            sys.stdout = old_stdout
            out = captured.getvalue()
            if out.strip():
                self.write(f"[BLOCK OUTPUT]\n{out.strip()}")

            self.write(f"[END BLOCK] {name} - {time.strftime('%H:%M:%S')}")
            self.write("-" * 50)
