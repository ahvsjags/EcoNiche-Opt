from __future__ import annotations

import pandas as pd


def cox_placeholder(metadata: pd.DataFrame, score_column: str = "EcoNicheScore") -> pd.DataFrame:
    required = {score_column, "os", "os_event"}
    missing = sorted(required - set(metadata.columns))
    if missing:
        return pd.DataFrame(
            [{"analysis": "cox_proportional_hazards", "status": "RESULT_PENDING", "reason": f"missing columns: {','.join(missing)}"}]
        )
    try:
        from lifelines import CoxPHFitter
    except Exception:
        return pd.DataFrame([{"analysis": "cox_proportional_hazards", "status": "unavailable_with_reason", "reason": "lifelines_not_installed"}])
    df = metadata[[score_column, "os", "os_event"]].dropna()
    if len(df) < 5:
        return pd.DataFrame([{"analysis": "cox_proportional_hazards", "status": "RESULT_PENDING", "reason": "too_few_samples"}])
    fitter = CoxPHFitter().fit(df, duration_col="os", event_col="os_event")
    summary = fitter.summary.reset_index().rename(columns={"covariate": "term"})
    summary["status"] = "available"
    return summary
