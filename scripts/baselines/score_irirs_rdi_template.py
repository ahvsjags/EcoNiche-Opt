from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
out = ROOT / "results/real/iris_rdi_status.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("model_name\tstatus\treason\nIRIS_RDI\tunavailable_with_reason\tRequires published interaction model artifacts not bundled\n", encoding="utf-8")
print(f"Wrote {out}")
