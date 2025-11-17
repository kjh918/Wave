# src/wave/core/task_registry.py
from __future__ import annotations
from typing import Dict, Type, List, TYPE_CHECKING

if TYPE_CHECKING:
    from .task import Task  # 순환 import 방지용 타입 힌트

class TaskRegistry:
    """아주 심플한 전역 태스크 레지스트리."""
    _REG: Dict[str, Type["Task"]] = {}

    @classmethod
    def register(cls, key: str, task_cls: Type["Task"]) -> None:
        """
        key: "fastqc", "gatk4.baserecalibrator", "picard.fastqtosam" 처럼
             TOOL[.FUNC] 형태로 통일
        """
        if key in cls._REG:
            # 덮어써도 되지만, 디버깅 편하게 경고만 찍을 수도 있음
            # print(f"[WAVE] Task '{key}' is already registered, overriding.")
            pass
        cls._REG[key] = task_cls

    @classmethod
    def get(cls, key: str) -> Type["Task"]:
        print(key)
        try:
            return cls._REG[key]
        except KeyError:
            raise KeyError(
                f"Unknown task TYPE: {key}. "
                f"Registered: {list(cls._REG.keys())}"
            )

    @classmethod
    def list_keys(cls) -> List[str]:
        return list(cls._REG.keys())


def register_task(key: str):
    """
    사용 예:
        @register_task("fastqc")
        class FastQCRunner(Task): ...
        
        @register_task("gatk4.baserecalibrator")
        class GatkBaseRecalibratorTask(Task): ...
    """
    def deco(cls: Type["Task"]):
        TaskRegistry.register(key, cls)
        return cls
    return deco