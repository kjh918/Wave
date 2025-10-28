# src/tasks/loader.py
from __future__ import annotations
import importlib
from typing import Type
from src.tasks.task import Task

class TaskLoadError(RuntimeError): ...

def load_task_class(tool: str) -> Type[Task]:
    """
    버전 없이 tool 폴더의 __main__에서 TASK_CLASS를 가져온다.
    예: src/tasks/fastqc/__main__.py → TASK_CLASS
    """
    mod_name = f"src.tasks.{tool}.__main__"
    try:
        mod = importlib.import_module(mod_name)
    except Exception as e:
        raise TaskLoadError(f"Cannot import {mod_name}: {e}")
    cls = getattr(mod, "TASK_CLASS", None)
    # print(issubclass(cls, Task))
    # print(Task)
    # print(cls)
    # exit()
    # if not isinstance(cls, type) or not issubclass(cls, Task):
    #     raise TaskLoadError(f"TASK_CLASS missing or invalid in {mod_name}")
    return cls
