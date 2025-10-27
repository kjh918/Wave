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
    각 Task는 아래 스키마를 갖고, `to_sh()` 로 실행할 쉘 라인을 돌려줍니다.
    - INPUTS:   필수/옵션 인풋 키들(설명/템플릿 허용)
    - OUTPUTS:  산출물 경로나 템플릿
    - DEFAULTS: 기본 파라미터 값
    - OPTIONAL: 옵션 파라미터 키 집합
    """
    TYPE: str  # subclass 필수 지정
    INPUTS: Dict[str, Any] = {}
    OUTPUTS: Dict[str, Any] = {}
    DEFAULTS: Dict[str, Any] = {}
    OPTIONAL: set[str] | list[str] = set()

    def __init__(self, name: str, version: str, params: Dict[str, Any], workdir: Path):
        # 최소 보관만; 검증/템플릿 전개는 다음 단계에서 추가
        self.name = name
        self.version = version
        self.params = {**self.DEFAULTS, **(params or {})}
        self.workdir = Path(workdir)
