import pandas as pd

from econiche_opt.reporting.claim_gate import gate_claim, sanitize_claim


def test_claim_gate_blocks_unsupported_superiority():
    result = gate_claim("EcoNiche-Opt outperforms all baselines.")
    assert not result.allowed
    assert result.status == "RESULT_PENDING"


def test_claim_gate_allows_supported_superiority():
    evidence = pd.DataFrame([{"test": "paired_bootstrap", "q_value": 0.01, "status": "SUPPORTED"}])
    result = gate_claim("EcoNiche-Opt is superior to the baseline.", evidence=evidence)
    assert result.allowed
    assert result.status == "SUPPORTED"


def test_sanitize_claim_rewrites_unsupported_language():
    sanitized = sanitize_claim("EcoNiche-Opt is the best model.")
    assert "best" not in sanitized.lower()
