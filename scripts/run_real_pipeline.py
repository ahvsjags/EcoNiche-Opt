from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.pipeline.real import RealPipelineConfig, assert_real_pipeline_ok, run_real_pipeline


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", default="config/data_registry.yml")
    parser.add_argument("--results-dir", default="results/real")
    parser.add_argument("--out-dir", default="results/real_pipeline")
    parser.add_argument("--execute-download", action="store_true")
    parser.add_argument("--execute-preprocess", action="store_true")
    parser.add_argument("--execute-training", action="store_true")
    parser.add_argument("--skip-secondary", action="store_true")
    parser.add_argument("--allow-missing-results", action="store_true")
    args = parser.parse_args()

    cfg = RealPipelineConfig(
        root=ROOT,
        registry=ROOT / args.registry,
        results_dir=ROOT / args.results_dir,
        out_dir=ROOT / args.out_dir,
        execute_download=args.execute_download,
        execute_preprocess=args.execute_preprocess,
        execute_training=args.execute_training,
        execute_secondary=not args.skip_secondary,
        strict_existing_results=not args.allow_missing_results,
    )
    manifest = run_real_pipeline(cfg)
    print(manifest[["stage", "status", "reason"]].to_string(index=False))
    assert_real_pipeline_ok(manifest)


if __name__ == "__main__":
    main()
