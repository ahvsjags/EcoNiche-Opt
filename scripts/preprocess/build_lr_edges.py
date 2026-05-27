from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    out = ROOT / "data/priors/ligand_receptor_edges.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    edges = pd.DataFrame(
        [
            ("CXCL12", "CXCR4", "caf_ecm_exclusion", "tnk_effector"),
            ("TGFB1", "TGFBR1", "caf_ecm_exclusion", "tnk_effector"),
            ("IL10", "IL10RA", "myeloid_suppression", "tnk_effector"),
        ],
        columns=["ligand", "receptor", "source_state", "target_state"],
    )
    edges.to_csv(out, sep="\t", index=False)
    print(f"Wrote ligand-receptor prior to {out}")


if __name__ == "__main__":
    main()
