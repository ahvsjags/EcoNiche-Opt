from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
out = ROOT / "results/real/tide_status.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("model_name\tstatus\treason\nTIDE\tunavailable_with_reason\tExternal TIDE service/code not bundled; local template signatures are in run_baselines.py\n", encoding="utf-8")
print(f"Wrote {out}")
