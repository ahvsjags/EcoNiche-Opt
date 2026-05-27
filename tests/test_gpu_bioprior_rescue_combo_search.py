from __future__ import annotations

from scripts.analysis.run_gpu_bioprior_rescue_combo_search import candidate_specs


def test_lipid_pi3k_prior_specs_are_constrained_and_include_robust_candidate():
    specs = candidate_specs("lipid_pi3k")
    candidates = {spec["candidate"] for spec in specs}
    genes = {spec["gene"] for spec in specs if spec["gene"]}

    assert genes == {"PLA2G2D", "PIK3CD"}
    assert "0.80*base+0.20*rz__PLA2G2D" in candidates
    assert "base_rescue_robust" in candidates


def test_robust_policy_excludes_percentile_candidates():
    specs = candidate_specs("lipid_pi3k", transform_policy="robust_only")
    candidates = {spec["candidate"] for spec in specs}

    assert "0.80*base+0.20*rz__PLA2G2D" in candidates
    assert "0.80*base+0.20*pct__PLA2G2D" not in candidates


def test_immune_prior_specs_include_checkpoint_and_tcell_axes():
    specs = candidate_specs("immune")
    genes = {spec["gene"] for spec in specs if spec["gene"]}

    assert {"CD247", "TIGIT", "SLAMF7", "MAP4K1"}.issubset(genes)
