from __future__ import annotations
import abc
from pathlib import Path
from typing import Any, Dict, Iterable, List, Callable

# ---- Registry ----
class TaskRegistry:
    _reg: dict[str, type] = {}

    @classmethod
    def register(cls, t: type):
        name = getattr(t, "TYPE", None)
        if not name:
            raise ValueError("Task class must define TYPE")
        if name in cls._reg:
            raise ValueError(f"Task TYPE already registered: {name}")
        cls._reg[name] = t
        return t

    @classmethod
    def get(cls, name: str):
        if name not in cls._reg:
            raise KeyError(f"Unknown task TYPE: {name}")
        return cls._reg[name]

    @classmethod
    def all(cls) -> Dict[str, type]:
        return dict(cls._reg)

# ---- Base Task ----
class Task(abc.ABC):
    """
    모든 Task 클래스의 기본 부모.
    각 Runner(FastQCRunner, FastPRunner 등)가 상속받아서 사용.
    """

    def __init__(
        self,
        workdir: str,
        inputs: Optional[Dict[str, Any]] = None,
        outputs: Optional[Dict[str, Any]] = None,
        params: Optional[Dict[str, Any]] = None,
    ):
        self.workdir = Path(workdir)
        self.inputs = inputs or {}
        self.outputs = outputs or {}
        self.params = params or {}

        # Task 실행 전 출력 디렉토리 보장
        self.workdir.mkdir(parents=True, exist_ok=True)
