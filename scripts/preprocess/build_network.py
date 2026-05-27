from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    out = ROOT / "data/priors/string_edges.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    edges = pd.DataFrame(
        [
            ("AXL", "VIM", 0.9),
            ("HLA-A", "B2M", 0.95),
            ("GZMB", "PRF1", 0.9),
            ("PDCD1", "LAG3", 0.8),
            ("COL1A1", "FN1", 0.85),
            ("S100A8", "S100A9", 0.95),
        ],
        columns=["gene_a", "gene_b", "confidence"],
    )
    edges.to_csv(out, sep="\t", index=False)
    edges.to_parquet(ROOT / "data/priors/network_adjacency.parquet", index=False)
    print(f"Wrote network prior to {out}")


if __name__ == "__main__":
    main()
