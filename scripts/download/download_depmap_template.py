from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
out = ROOT / "results/perturbation/depmap_download_status.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("status\treason\nRESULT_PENDING\tConfigure DepMap release source\n", encoding="utf-8")
print(f"Wrote {out}")
