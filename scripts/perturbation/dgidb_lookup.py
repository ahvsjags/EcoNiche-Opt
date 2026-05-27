from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from econiche_opt.perturbation.reversal import build_reversal_signature  # noqa: E402

DGIDB_GRAPHQL = "https://dgidb.org/api/graphql"


def _chunks(values: list[str], size: int) -> list[list[str]]:
    return [values[idx : idx + size] for idx in range(0, len(values), size)]


def query_dgidb(genes: list[str], timeout: int = 45) -> pd.DataFrame:
    query = """
    query($genes: [String!]) {
      interactions(geneNames: $genes, first: 1000) {
        nodes {
          id
          interactionTypes { type directionality }
          drug { name conceptId }
          gene { name conceptId }
          sources { sourceDbName }
          publications { pmid }
        }
      }
    }
    """
    rows = []
    for chunk in _chunks(genes, 40):
        response = requests.post(DGIDB_GRAPHQL, json={"query": query, "variables": {"genes": chunk}}, timeout=timeout)
        response.raise_for_status()
        payload = response.json()
        if payload.get("errors"):
            raise RuntimeError(str(payload["errors"]))
        for node in payload.get("data", {}).get("interactions", {}).get("nodes", []):
            interactions = node.get("interactionTypes") or []
            sources = node.get("sources") or []
            publications = node.get("publications") or []
            rows.append(
                {
                    "target_gene": (node.get("gene") or {}).get("name", ""),
                    "drug_name": (node.get("drug") or {}).get("name", ""),
                    "drug_concept_id": (node.get("drug") or {}).get("conceptId", ""),
                    "interaction_types": ";".join(
                        sorted(
                            {
                                f"{item.get('type', '')}:{item.get('directionality', '')}".strip(":")
                                for item in interactions
                                if item
                            }
                        )
                    ),
                    "sources": ";".join(sorted({str(item.get("sourceDbName", "")) for item in sources if item})),
                    "pmids": ";".join(sorted({str(item.get("pmid", "")) for item in publications if item.get("pmid")})),
                    "status": "DGIdb_evidence",
                }
            )
    return pd.DataFrame(rows).drop_duplicates() if rows else pd.DataFrame()


def _pending(reason: str) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "target_gene": "RESULT_PENDING",
                "drug_name": "RESULT_PENDING",
                "drug_concept_id": "",
                "interaction_types": "",
                "sources": "",
                "pmids": "",
                "status": "RESULT_PENDING",
                "reason": reason,
            }
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Annotate EcoNiche module genes with DGIdb drug-gene interactions.")
    parser.add_argument("--module", default="results/real/econiche_module.tsv")
    parser.add_argument("--out", default="results/perturbation/dgidb_hits.tsv")
    parser.add_argument("--max-genes", type=int, default=80)
    parser.add_argument("--timeout", type=int, default=45)
    args = parser.parse_args()

    module_path = ROOT / args.module if not Path(args.module).is_absolute() else Path(args.module)
    out = ROOT / args.out if not Path(args.out).is_absolute() else Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    if not module_path.exists():
        result = _pending(f"module table missing: {module_path}")
    else:
        module = pd.read_csv(module_path, sep="\t")
        signature = build_reversal_signature(module, top_genes=args.max_genes)
        genes = list(dict.fromkeys(signature["resistance_up"] + signature["resistance_down"]))[: args.max_genes]
        if not genes:
            result = _pending("no module genes available for DGIdb query")
        else:
            try:
                result = query_dgidb(genes, timeout=args.timeout)
                if result.empty:
                    result = _pending("DGIdb returned no interactions for selected module genes")
            except (requests.RequestException, RuntimeError) as exc:
                result = _pending(f"DGIdb query failed: {exc}")
    result.to_csv(out, sep="\t", index=False)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
