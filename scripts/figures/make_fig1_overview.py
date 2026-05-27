from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.plotting import write_pending_figure

write_pending_figure(ROOT / "figures/fig1_overview.pdf", "Study design and EcoNiche-Opt overview", "registry -> priors -> model -> benchmark -> hypotheses")
print("Wrote figures/fig1_overview.pdf")
