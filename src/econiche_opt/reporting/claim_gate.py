from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import pandas as pd

SUPPORTED_STATUSES = {"SUPPORTED", "NOT_SIGNIFICANT", "RESULT_PENDING"}
RESTRICTED_WORDS = re.compile(
    r"\b(outperform(?:s|ed|ing)?|superior|best|beats?|significant(?:ly)? better|state-of-the-art)\b",
    flags=re.IGNORECASE,
)


@dataclass(frozen=True)
class ClaimGateResult:
    allowed: bool
    status: str
    reason: str
    evidence_rows: int = 0


def _has_supported_evidence(evidence: pd.DataFrame | None, alpha: float) -> bool:
    if evidence is None or evidence.empty:
        return False
    q_cols = [col for col in evidence.columns if col.lower() in {"q_value", "fdr", "q"}]
    p_cols = [col for col in evidence.columns if col.lower() in {"p_value", "p"}]
    status_ok = True
    if "status" in evidence.columns:
        status_ok = evidence["status"].astype(str).str.upper().isin(SUPPORTED_STATUSES).any()
    q_ok = any(pd.to_numeric(evidence[col], errors="coerce").le(alpha).any() for col in q_cols)
    p_ok = any(pd.to_numeric(evidence[col], errors="coerce").le(alpha).any() for col in p_cols)
    paired_ok = True
    if "test" in evidence.columns:
        paired_ok = evidence["test"].astype(str).str.contains("paired|delong|bootstrap", case=False, regex=True).any()
    return bool(status_ok and paired_ok and (q_ok or p_ok))


def gate_claim(claim: str, evidence: pd.DataFrame | None = None, alpha: float = 0.05) -> ClaimGateResult:
    if not RESTRICTED_WORDS.search(claim):
        return ClaimGateResult(True, "ALLOWED_DESCRIPTIVE", "Claim does not assert superiority.", 0)
    if _has_supported_evidence(evidence, alpha):
        return ClaimGateResult(True, "SUPPORTED", "Superiority language is backed by statistical evidence.", len(evidence))
    return ClaimGateResult(
        False,
        "RESULT_PENDING",
        "Superiority or best-model language requires paired bootstrap/DeLong evidence with FDR support.",
        0 if evidence is None else len(evidence),
    )


def sanitize_claim(claim: str, evidence: pd.DataFrame | None = None) -> str:
    result = gate_claim(claim, evidence=evidence)
    if result.allowed:
        return claim
    return RESTRICTED_WORDS.sub("does not yet demonstrate superiority over", claim)


def validate_claims(claims: Iterable[str], evidence: pd.DataFrame | None = None) -> list[ClaimGateResult]:
    return [gate_claim(claim, evidence=evidence) for claim in claims]


def validate_claim_file(path: str | Path, evidence_path: str | Path | None = None) -> ClaimGateResult:
    text = Path(path).read_text(encoding="utf-8")
    evidence = None
    if evidence_path and Path(evidence_path).exists():
        evidence = pd.read_csv(evidence_path, sep="\t")
    return gate_claim(text, evidence=evidence)
