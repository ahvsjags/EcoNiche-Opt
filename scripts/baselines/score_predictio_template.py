from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
out = ROOT / "results/real/predictio_status.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("model_name\tstatus\treason\nPredictIO\tunavailable_with_reason\tConfigure PredictIO signature before scoring\n", encoding="utf-8")
print(f"Wrote {out}")
