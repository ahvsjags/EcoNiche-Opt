from __future__ import annotations

import json
from pathlib import Path

from scripts.analysis.run_gpu_lipid_pair_rescue_search import _candidate_specs
from scripts.reporting.make_gpu_lipid_pair_rescue_package import PATTERN


ROOT = Path(__file__).resolve().parents[1]


def test_component_dominant_specs_constrain_base_weight():
    specs = _candidate_specs("component_dominant")
    assert specs
    assert all(float(spec["weight_base"]) <= 0.35 for spec in specs)
    assert any("pct__PLA2G2D" in str(spec["candidate"]) and "rz__PIK3CD" in str(spec["candidate"]) for spec in specs)


def test_lipid_pair_package_freezes_selected_candidate():
    spec_path = ROOT / "deliverables/gpu_lipid_pair_rescue_package_20260528/gpu_lipid_pair_rescue_scoring_spec.json"
    assert spec_path.exists()
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    assert spec["model_id"] == "EcoNiche-Opt-GPU-LipidPI3K-PairRescue"
    assert PATTERN.fullmatch(spec["candidate"])
    assert spec["claim_boundary"]["external_or_holdout_labels_used_for_training"] is False
    assert spec["claim_boundary"]["external_or_holdout_labels_used_for_feature_selection"] is False
    assert spec["locked_score"]["locked_threshold"] > 0
