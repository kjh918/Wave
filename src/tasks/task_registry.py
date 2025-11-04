# src/tasks/task_registry.py
from __future__ import annotations
from typing import Dict, Type, Optional

class TaskRegistry:
    _registry: Dict[str, Type] = {}

    @classmethod
    def register(cls, cls_obj: Type):
        key = getattr(cls_obj, "TYPE", None)
        if not key:
            raise ValueError(f"Task class {cls_obj.__name__} missing TYPE")
        key = key.lower()
        if key in cls._registry and cls._registry[key] is not cls_obj:
            raise KeyError(f"Task TYPE already registered: {key}")
        cls._registry[key] = cls_obj
        return cls_obj

    @classmethod
    def get(cls, key: str) -> Type:
        key = key.lower()
        if key not in cls._registry:
            raise KeyError(f"Unknown task TYPE: {key}")
        return cls._registry[key]

def register_task(type_name: Optional[str] = None):
    def deco(cls_obj: Type):
        if type_name:
            cls_obj.TYPE = type_name
        return TaskRegistry.register(cls_obj)
    return deco