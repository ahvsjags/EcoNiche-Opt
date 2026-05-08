from __future__ import annotations


def estimate_matrix_megabytes(n_samples: int, n_genes: int, bytes_per_value: int = 8) -> float:
    return float(n_samples * n_genes * bytes_per_value / (1024**2))


def assert_matrix_within_memory(n_samples: int, n_genes: int, max_megabytes: float = 2048.0) -> None:
    estimated = estimate_matrix_megabytes(n_samples, n_genes)
    if estimated > max_megabytes:
        raise MemoryError(
            f"Estimated expression matrix size {estimated:.1f} MB exceeds limit {max_megabytes:.1f} MB."
        )
