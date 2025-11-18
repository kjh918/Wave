# src/wave/core/loader.py
from __future__ import annotations
import importlib
import pkgutil
from typing import Optional

from wave.core.task_registry import TaskRegistry  # 디버깅 용도


def autoload_tasks(package_root: str = "wave.tasks") -> None:
    """
    wave.tasks 하위에서 *.main 모듈을 자동으로 import.
    각 main.py 안에서 @register_task(...) 가 실행되면서
    TaskRegistry 에 등록된다.
    """
    try:
        pkg = importlib.import_module(package_root)
    except ModuleNotFoundError as e:
        print(f"[WAVE] cannot import package_root='{package_root}': {e}")
        return

    print(f"[WAVE] autoload from {package_root}, path={pkg.__path__}")

    for m in pkgutil.walk_packages(pkg.__path__, pkg.__name__ + "."):
        # m.name 예시:
        #   wave.tasks.fastqc
        #   wave.tasks.fastp
        #   wave.tasks.gatk4
        #   wave.tasks.gatk4.haplotypecaller
        #   wave.tasks.gatk4.haplotypecaller.main
        name = m.name
        print("  - found:", name, "ispkg=", m.ispkg)

        # 1) main.py만 골라서 import
        if name.endswith(".main"):
            try:
                print("    import:", name)
                importlib.import_module(name)
            except Exception as e:
                print(f"[WAVE] autoload skip {name}: {e}")
            continue

        # 2) 패키지인 경우, 그 아래의 main들을 위해 패키지만 먼저 로딩
        if m.ispkg:
            try:
                importlib.import_module(name)
            except Exception as e:
                print(f"[WAVE] autoload skip pkg {name}: {e}")
                continue

    print(f"[WAVE] registered tasks: {list(TaskRegistry._reg.keys())}")