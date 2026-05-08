from econiche.module import EcoNicheConfig, EcoNicheModule
from econiche.objective import module_size_penalty, robust_objective


def test_large_module_has_larger_size_penalty():
    cfg = EcoNicheConfig(min_genes_per_state=2, max_genes_per_state=5)
    small = EcoNicheModule({"s": {"g1", "g2", "g3"}})
    large = EcoNicheModule({"s": {f"g{i}" for i in range(20)}})

    assert module_size_penalty(small, cfg) < module_size_penalty(large, cfg)


def test_biological_terms_raise_objective_and_penalties_lower_it():
    cfg = EcoNicheConfig()
    good = robust_objective(
        cfg,
        auc_values=[0.8, 0.82],
        auprc_values=[0.7, 0.72],
        ece_mean=0.05,
        biological_terms={"cell_specificity": 0.8, "pathway": 0.6, "network": 0.5, "lr": 0.4, "stability": 0.7},
        penalties={"size": 0.1, "batch": 0.0, "leakage": 0.0, "redundancy": 0.0, "therapy_confounding": 0.0},
    )["score"]
    penalized = robust_objective(
        cfg,
        auc_values=[0.8, 0.82],
        auprc_values=[0.7, 0.72],
        ece_mean=0.05,
        biological_terms={"cell_specificity": 0.0, "pathway": 0.0, "network": 0.0, "lr": 0.0, "stability": 0.0},
        penalties={"size": 3.0, "batch": 2.0, "leakage": 1.0, "redundancy": 2.0, "therapy_confounding": 2.0},
    )["score"]

    assert good > penalized
