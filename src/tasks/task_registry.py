# task_registry.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Type, Optional, List

@dataclass(frozen=True)
class TaskKey:
    name: str         # 안정 ID (예: "fastqc")
    version: str = "" # "v1", "v2" 등. 빈 문자열이면 기본

class TaskBase: ...
# ↑ 공통 인터페이스 (build_command, output_files 등)

class TaskRegistry:
    def __init__(self):
        # name -> version -> class
        self._impls: Dict[str, Dict[str, Type[TaskBase]]] = {}
        # alias -> TaskKey(name, version)
        self._aliases: Dict[str, TaskKey] = {}

    def register(self, name: str, cls: Type[TaskBase], version: str = "", aliases: Optional[List[str]] = None):
        self._impls.setdefault(name, {})
        self._impls[name][version] = cls
        for al in aliases or []:
            self._aliases[al] = TaskKey(name, version)

    def resolve(self, id_or_alias: str, version: Optional[str] = None) -> Type[TaskBase]:
        # 1) alias면 본 이름으로
        if id_or_alias in self._aliases:
            tk = self._aliases[id_or_alias]
            name, v = tk.name, (version or tk.version)
        else:
            name, v = id_or_alias, (version or "")
        if name not in self._impls:
            raise KeyError(f"Unknown task '{name}'")
        impls = self._impls[name]
        if v in impls:                     # 지정 버전
            return impls[v]
        if "" in impls and not version:    # 기본 구현
            return impls[""]
        # 최신 버전 선택(규칙: vN 숫자 최대)
        if not version:
            cand = sorted((k for k in impls if k.startswith("v")), key=lambda s: int(s[1:]), reverse=True)
            if cand: return impls[cand[0]]
        raise KeyError(f"No implementation for '{name}' (version='{version}')")

TASKS = TaskRegistry()

def register_task(name: str, version: str = "", aliases: Optional[List[str]] = None):
    def deco(cls: Type[TaskBase]):
        TASKS.register(name=name, cls=cls, version=version, aliases=aliases)
        return cls
    return deco
