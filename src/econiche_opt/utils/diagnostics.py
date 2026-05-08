from __future__ import annotations


class EcoNicheDiagnosticError(RuntimeError):
    """Raised for user-actionable pipeline errors."""


def require_file(path, message: str | None = None) -> None:
    from pathlib import Path

    if not Path(path).exists():
        raise EcoNicheDiagnosticError(message or f"Required file not found: {path}")
