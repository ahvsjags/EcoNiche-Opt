from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.plotting import write_pending_figure

write_pending_figure(ROOT / "figures/fig3_module_network.pdf", "Identified multicellular resistance niche", "Run EcoNiche-Opt and network priors for full graph")
print("Wrote figures/fig3_module_network.pdf")
