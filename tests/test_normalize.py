import numpy as np
import pandas as pd

from econiche.normalize import rank_gaussian_normalize


def test_rank_gaussian_normalize_preserves_rowwise_order_and_finite_values():
    x = pd.DataFrame(
        [[1.0, 10.0, 5.0], [7.0, 3.0, 9.0]],
        index=["s1", "s2"],
        columns=["g1", "g2", "g3"],
    )

    out = rank_gaussian_normalize(x)

    assert out.shape == x.shape
    assert np.isfinite(out.to_numpy()).all()
    assert out.loc["s1", "g2"] > out.loc["s1", "g3"] > out.loc["s1", "g1"]
    assert out.loc["s2", "g3"] > out.loc["s2", "g1"] > out.loc["s2", "g2"]
