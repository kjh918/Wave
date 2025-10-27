# src/tasks/loader.py
from __future__ import annotations
import importlib
from typing import Type, Optional
from src.tasks.task import Task

class TaskLoadError(RuntimeError): ...

    def load_task_class(tool: str, func: str, version: Optional[str]) -> Type[Task]:
        """
        tool/func/version 기반으로 모듈을 import하고, TASK_IMPLS에서 버전에 해당하는 Task 클래스를 반환.
        """
        mod_name = f"src.tasks.{tool}.{func}"
        try:
            mod = importlib.import_module(mod_name)
        except Exception as e:
            raise TaskLoadError(f"Failed to import module '{mod_name}': {e}")

        # 1) 명시 맵 우선
        impls = getattr(mod, "TASK_IMPLS", None)
        if isinstance(impls, dict):
            if version and version in impls:
                return impls[version]
            # fallback: "", "default" 키, 또는 가장 최신으로 보이는 키
            for k in ("", "default"):
                if k in impls:
                    return impls[k]
            # 버전 문자열이 존재하면 가장 큰 버전(숫자 비교) 선택 시도
            try:
                def _key(v):
                    # "v1.2.3" → (1,2,3)
                    return tuple(int(x) for x in v.lstrip("v").split("."))
                best = sorted((k for k in impls.keys() if k), key=_key, reverse=True)[0]
                return impls[best]
            except Exception:
                pass

        # 2) 단일 클래스 노출 (TaskClass)
        cls = getattr(mod, "TaskClass", None)
        if cls:
            return cls

        # 3) 모듈 내에서 Task 서브클래스 한 개만 존재하는 경우 자동 선택
        candidates = []
        for v in mod.__dict__.values():
            try:
                if isinstance(v, type) and issubclass(v, Task) and v is not Task:
                    candidates.append(v)
            except Exception:
                pass
        if len(candidates) == 1:
            return candidates[0]

        raise TaskLoadError(f"No usable Task class found in '{mod_name}'. "
                            f"Expose TASK_IMPLS or TaskClass.")
