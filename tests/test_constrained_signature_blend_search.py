from __future__ import annotations

import pandas as pd

from scripts.analysis.run_constrained_signature_blend_search import blend_score, candidate_weight_specs


def test_constrained_blend_candidates_are_normalized_and_include_module_anchor():
    specs = candidate_weight_specs()
    assert any(spec["candidate"] == "single__EcoNiche-Opt-ModulePriorFixed" for spec in specs)
    assert any(spec["candidate"].startswith("module_plus_immune_family") for spec in specs)
    for spec in specs:
        assert abs(sum(spec["weights"].values()) - 1.0) < 1e-9
        assert all(weight > 0 for weight in spec["weights"].values())


def test_blend_score_uses_available_weighted_features():
    table = pd.DataFrame(
        {
            "EcoNiche-Opt-ModulePriorFixed": [1.0, 2.0],
            "CXCL9": [3.0, 5.0],
        },
        index=["S1", "S2"],
    )
    score = blend_score(table, {"EcoNiche-Opt-ModulePriorFixed": 0.75, "CXCL9": 0.25, "MISSING": 1.0})
    assert score.round(6).tolist() == [1.5, 2.75]
