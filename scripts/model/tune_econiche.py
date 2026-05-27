from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
out = ROOT / "results/real/tuning_status.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("status\treason\nRESULT_PENDING\tNested CV tuning hooks are available through LODO scaffold; configure real cohorts first\n", encoding="utf-8")
print(f"Wrote {out}")
