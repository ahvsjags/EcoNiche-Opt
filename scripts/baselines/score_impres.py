from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
out = ROOT / "results/real/impres_status.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
out.write_text("model_name\tstatus\treason\nIMPRES\tunavailable_with_reason\tOriginal pairwise parser requires curated checkpoint gene panel metadata\n", encoding="utf-8")
print(f"Wrote {out}")
