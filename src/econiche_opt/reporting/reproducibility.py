from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path


def generate_reproducibility_report(out: str | Path, mode: str = "demo") -> Path:
    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# EcoNiche-Opt Reproducibility Report",
        "",
        f"- mode: {mode}",
        f"- generated_utc: {datetime.now(timezone.utc).isoformat()}",
        "- data_policy: real outputs must come from registered download/preprocess pipelines",
        "- restricted_data_policy: ACCESS_RESTRICTED datasets are documented, not imputed",
        "- validation_entrypoint: python -m econiche_opt.cli validate-project --mode demo",
        "",
        "## Core commands",
        "",
        "```bash",
        "python -m econiche_opt.cli make-demo",
        "python -m econiche_opt.cli run-benchmark --demo",
        "python -m econiche_opt.cli make-figures --demo",
        "python -m econiche_opt.cli make-manuscript --demo",
        "python -m econiche_opt.cli validate-project --mode demo",
        "```",
    ]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path
