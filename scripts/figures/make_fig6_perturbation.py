from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.plotting import write_pending_figure

write_pending_figure(ROOT / "figures/fig6_perturbation.pdf", "Perturbation reversal", "candidate axes only; no treatment recommendation")
print("Wrote figures/fig6_perturbation.pdf")
