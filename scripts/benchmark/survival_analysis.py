from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche.survival import cox_placeholder
from econiche.expression import clean_expression_matrix, load_entrez_symbol_map, read_table_matrix
from econiche.scoring import compute_state_scores
from econiche.module import EcoNicheModule


def main() -> None:
    out = ROOT / "results/real/survival_cox.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    survival = run_gse183924_survival()
    if survival is None:
        pred_path = ROOT / "results/real/lodo_predictions.tsv"
        metadata = pd.read_csv(pred_path, sep="\t") if pred_path.exists() else pd.DataFrame()
        cox_placeholder(metadata).to_csv(out, sep="\t", index=False)
        (ROOT / "results/real/survival_km.tsv").write_text("status\treason\nRESULT_PENDING\tSurvival columns not materialized\n", encoding="utf-8")
    else:
        cox, km = survival
        cox.to_csv(out, sep="\t", index=False)
        km.to_csv(ROOT / "results/real/survival_km.tsv", sep="\t", index=False)
    print(f"Wrote survival status to {out}")


def run_gse183924_survival() -> tuple[pd.DataFrame, pd.DataFrame] | None:
    module_path = ROOT / "results/real/econiche_module.tsv"
    coef_path = ROOT / "results/real/coefficients.tsv"
    expr_path = ROOT / "data/raw/GSE183924/suppl/GSE183924_FPKM.txt.gz"
    clinical_path = ROOT / "data/raw/GSE183924/suppl/GSE183924_Clinical_Data_with_CIBERSORT_analysis_for_GEO.xlsx"
    if not (module_path.exists() and coef_path.exists() and expr_path.exists() and clinical_path.exists()):
        return None
    module_table = pd.read_csv(module_path, sep="\t")
    module = EcoNicheModule({state: set(frame["gene"].astype(str)) for state, frame in module_table.groupby("state")})
    directions = dict(zip(module_table["gene"].astype(str), module_table["direction"].astype(int)))
    expr = clean_expression_matrix(read_table_matrix(expr_path, sep="\t"), entrez_map=load_entrez_symbol_map(ROOT / "data/priors/Homo_sapiens.gene_info.gz"))
    state_scores = compute_state_scores(expr, module, directions, normalize=True)
    coefficients = pd.read_csv(coef_path, sep="\t").set_index("feature")["coefficient"].to_dict()
    score = pd.Series(float(coefficients.get("intercept", 0.0)), index=state_scores.index)
    for state in state_scores.columns:
        score = score + state_scores[state] * float(coefficients.get(state, 0.0))
    clinical = pd.read_excel(clinical_path)
    clinical = clinical.rename(columns={"RNA-Seq ID": "sample_id", "RFS (months)": "rfs_months", "RFS Censor": "rfs_censor"})
    clinical = clinical.set_index("sample_id", drop=False)
    common = sorted(set(score.index) & set(clinical.index))
    if len(common) < 5:
        return None
    df = clinical.loc[common, ["sample_id", "rfs_months", "rfs_censor"]].copy()
    df["EcoNicheScore"] = score.loc[common].values
    df["rfs_event"] = df["rfs_censor"].astype(float)
    from lifelines import CoxPHFitter

    cox_input = df[["rfs_months", "rfs_event", "EcoNicheScore"]].dropna()
    if len(cox_input) < 5 or cox_input["rfs_event"].nunique() < 2:
        cox = pd.DataFrame([{"cohort": "GSE183924", "analysis": "cox", "status": "RESULT_PENDING", "reason": "too_few_events"}])
    else:
        fitter = CoxPHFitter().fit(cox_input, duration_col="rfs_months", event_col="rfs_event")
        cox = fitter.summary.reset_index().rename(columns={"covariate": "term"})
        cox["cohort"] = "GSE183924"
        cox["endpoint"] = "relapse_free_survival"
        cox["status"] = "available"
    median = df["EcoNicheScore"].median()
    df["risk_group"] = ["high" if value >= median else "low" for value in df["EcoNicheScore"]]
    km = df.groupby("risk_group").agg(n=("sample_id", "count"), median_rfs_months=("rfs_months", "median"), events=("rfs_event", "sum")).reset_index()
    km["cohort"] = "GSE183924"
    return cox, km


if __name__ == "__main__":
    main()
