from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def table_status(path: Path, label: str) -> str:
    if not path.exists():
        return f"[RESULT_PENDING: {label}]"
    frame = pd.read_csv(path, sep="\t")
    return f"{label}: {len(frame)} rows generated."


def main() -> None:
    summary = "\n".join(
        [
            "# Result Summaries",
            "",
            table_status(ROOT / "results/demo/lodo_metrics.tsv", "Demo LODO metrics"),
            table_status(ROOT / "results/real/lodo_metrics.tsv", "Real-data LODO metrics"),
            table_status(ROOT / "results/real/model_comparison_bootstrap.tsv", "Model comparison"),
            table_status(ROOT / "results/perturbation/prioritized_perturbations.tsv", "Perturbation prioritization"),
            "",
        ]
    )
    out = ROOT / "paper/result_summaries.md"
    out.write_text(summary, encoding="utf-8")
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
