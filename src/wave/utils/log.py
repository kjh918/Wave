import datetime
import json
import traceback
from functools import wraps

class Logger:
    def __init__(self, func=None, *, task_id=None):
        self.func = func
        self.task_id = task_id
        self.logs = []
        if func is not None:
            wraps(func)(self)

    def __call__(self, *args, **kwargs):
        # 표기: task_id가 있으면 붙이고, 없으면 생략
        task_info = f" ▶ {self.task_id}" if self.task_id is not None else ""

        start_ts = self.timestamp()
        start_time = datetime.datetime.now()

        print(f"[{start_ts}]{task_info} ▶ {self.func.__name__} START")
        try:
            result = self.func(*args, **kwargs)
            end_time = datetime.datetime.now()
            end_ts = end_time.strftime("%Y-%m-%d %H:%M:%S")
            duration = (end_time - start_time).total_seconds()
            print(f"[{end_ts}]{task_info} ▶ {self.func.__name__} END (Process Time : {duration:.4f}s)")
            self.save_log(start_ts, end_ts, duration, args, kwargs, result=result)
            return result
        except Exception as e:
            error_time = datetime.datetime.now()
            error_ts = error_time.strftime("%Y-%m-%d %H:%M:%S")
            duration = (error_time - start_time).total_seconds()

            print(f"[{error_ts}]{task_info} ▶ {self.func.__name__} ERROR (Process Time : {duration:.4f}s Log : {e})")

            tb = traceback.format_exc()
            self.save_log(start_ts, error_ts, duration, args, kwargs, error=str(e), traceback=tb)
            raise  # 에러 재전파 (workflow/executor에서 처리)

    def timestamp(self) -> str:
        """현재 시각을 'YYYY-MM-DD HH:MM:SS' 형식 문자열로 반환"""
        return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def save_log(self, start_ts, end_ts, duration, args, kwargs, result=None, error=None, traceback=None):
        entry = {
            "task_id": self.task_id,
            "function": self.func.__name__ if self.func else None,
            "start_time": start_ts,
            "end_time": end_ts,
            "duration_sec": duration,
            "args": args,
            "kwargs": kwargs,
            "result": result,
            "error": error,
            "traceback": traceback,
        }
        self.logs.append(entry)

    def get_logs(self):
        return self.logs

    def save_logs_to_file(self, path: str, mode: str = "a"):
        with open(path, mode, encoding="utf-8") as f:
            for rec in self.logs:
                f.write(json.dumps(rec, ensure_ascii=False) + "")

# 기본 로거
def logger(func):
    return Logger(func)

# task_id 버전 로거
def logger_with_task(task_id: str):
    def decorator(func):
        return Logger(func, task_id=task_id)
    return decorator