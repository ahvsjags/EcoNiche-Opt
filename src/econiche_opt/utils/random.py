from __future__ import annotations

import random

import numpy as np


def set_global_seed(seed: int = 42) -> int:
    random.seed(seed)
    np.random.seed(seed)
    return seed
