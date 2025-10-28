# src/tasks/task.py
from __future__ import annotations
import abc
from typing import Dict, Any
from pathlib import Path

class TaskRegistry:
    _reg: dict[str, type] = {}

    @classmethod
    def register(cls, t: type):
        name = getattr(t, "TYPE", None)
        if not name:
            raise ValueError("Task class must define TYPE")
        cls._reg[name] = t
        return t

    @classmethod
    def get(cls, name: str):
        if name not in cls._reg:
            raise KeyError(f"Unknown task TYPE: {name}")
        return cls._reg[name]

class Task(abc.ABC):
    TYPE: str
    DEFAULTS: Dict[str, Any] = {}

    def __init__(self, name: str, workdir: Path,
                 inputs: Dict[str, Any] = None,
                 outputs: Dict[str, Any] = None,
                 params: Dict[str, Any] = None):
        self.name = name
        self.workdir = Path(workdir)
        self.inputs = inputs or {}
        self.outputs = outputs or {}
        self.params = {**(self.DEFAULTS or {}), **(params or {})}
        self.workdir.mkdir(parents=True, exist_ok=True)

    @abc.abstractmethod
    def to_sh(self):
        ...