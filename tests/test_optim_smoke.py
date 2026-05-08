from econiche.demo import make_synthetic_data
from econiche.model import EcoNicheOpt
from econiche.module import EcoNicheConfig
from econiche.priors import make_default_cell_state_priors


def test_optimizer_recovers_planted_signal_on_synthetic_data():
    demo = make_synthetic_data(n_cohorts=3, n_per_cohort=36, random_state=7)
    priors = make_default_cell_state_priors(demo["genes"])
    model = EcoNicheOpt(EcoNicheConfig(population_size=30, generations=12, random_state=7), priors=priors)

    result = model.fit(demo["X_by_cohort"], demo["y_by_cohort"], demo["metadata_by_cohort"])
    recovered = set(result.module_table()["gene"]) & set(demo["planted_genes"])

    assert len(recovered) >= 6
    assert result.lodo_metrics["AUROC"].mean() > 0.75
    assert not result.history.empty
    assert {"pred_prob", "EcoNicheScore", "sample_id", "cohort", "true_label"}.issubset(result.predictions.columns)
