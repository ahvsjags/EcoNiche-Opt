# No Fabrication Policy

EcoNiche-Opt separates executable demo results, public real-data results, and restricted or unavailable resources.

- Real cohort outputs must be produced by registered download and preprocessing commands.
- Controlled-access datasets are marked `ACCESS_RESTRICTED` and are never replaced by synthetic or manually invented values.
- Missing benchmark, survival, single-cell, perturbation, or baseline results are emitted as `RESULT_PENDING` or `unavailable_with_reason`.
- Superiority language is blocked unless paired bootstrap or DeLong comparisons with FDR support are present.
- Perturbation outputs are hypothesis-generation artifacts and are not treatment recommendations.
