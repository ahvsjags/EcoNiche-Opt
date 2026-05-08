from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

from econiche.model import EcoNicheOpt
from econiche.module import DEFAULT_STATES, EcoNicheConfig, EcoNicheModule, EcoNicheResult
from econiche_opt.api import EcoNicheOptClassifier
from econiche_opt.datasets import load_demo_multicohort

try:
    __version__ = version("econiche-opt")
except PackageNotFoundError:
    __version__ = "0.2.0"

__all__ = [
    "DEFAULT_STATES",
    "EcoNicheConfig",
    "EcoNicheModule",
    "EcoNicheOptClassifier",
    "EcoNicheOpt",
    "EcoNicheResult",
    "__version__",
    "load_demo_multicohort",
]
