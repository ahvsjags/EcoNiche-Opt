from __future__ import annotations

from typing import Any

import pandas as pd

from econiche.demo import make_synthetic_data


def load_demo_multicohort(
    n_cohorts: int = 3,
    n_per_cohort: int = 24,
    random_state: int = 42,
    response_positive: bool = True,
) -> dict[str, Any]:
    """Return a small synthetic multicohort dataset for package examples.

    The demo data are not used as real evidence. They exist only so users can
    test installation, APIs, and serialization without downloading public ICB
    cohorts.
    """

    demo = make_synthetic_data(n_cohorts=n_cohorts, n_per_cohort=n_per_cohort, random_state=random_state)
    if response_positive:
        demo["y_response_by_cohort"] = {
            cohort: pd.Series(1 - labels.astype(int), index=labels.index, name="response_label")
            for cohort, labels in demo["y_by_cohort"].items()
        }
    else:
        demo["y_response_by_cohort"] = {
            cohort: pd.Series(labels.astype(int), index=labels.index, name="nonresponse_label")
            for cohort, labels in demo["y_by_cohort"].items()
        }
    return demo
