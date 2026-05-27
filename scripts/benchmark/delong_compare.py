from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
out = ROOT / "results/real/delong_compare.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("status\treason\nRESULT_PENDING\tDeLong implementation not run; use paired bootstrap output as primary comparison scaffold\n", encoding="utf-8")
print(f"Wrote {out}")
