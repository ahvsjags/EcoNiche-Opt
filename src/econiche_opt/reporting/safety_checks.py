from __future__ import annotations

from pathlib import Path

import pandas as pd

from econiche_opt.reporting.claim_gate import gate_claim


def validate_manuscript_safety(path: str | Path, evidence_path: str | Path | None = None) -> pd.DataFrame:
    manuscript = Path(path)
    evidence = pd.read_csv(evidence_path, sep="\t") if evidence_path and Path(evidence_path).exists() else None
    rows = []
    for line_no, line in enumerate(manuscript.read_text(encoding="utf-8").splitlines(), start=1):
        result = gate_claim(line, evidence=evidence)
        if not result.allowed:
            rows.append(
                {
                    "file": str(manuscript),
                    "line": line_no,
                    "status": result.status,
                    "reason": result.reason,
                    "text": line,
                }
            )
    return pd.DataFrame(rows)
