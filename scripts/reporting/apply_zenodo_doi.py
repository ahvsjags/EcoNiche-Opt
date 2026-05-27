from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DOI_RE = re.compile(r"^10\.5281/zenodo\.\d+$")
DOI_URL_PREFIX = "https://doi.org/"


DEFAULT_TEXT_FILES = [
    Path("CITATION.cff"),
    Path("README.md"),
    Path("DATA_RESULTS_FIGURES_UPLOAD_NOTES.md"),
    Path("paper/econiche_opt_manuscript_en_v1_20260509.md"),
    Path("paper/Journal of Translational Medicine投稿/EcoNiche-Opt_JTM_Main_Manuscript.md"),
    Path("paper/communications_medicine_submission/communications_medicine_submission_readiness.md"),
    Path("deliverables/zenodo_manual_publication_steps_20260528.md"),
]

DEFAULT_MANIFEST = Path("deliverables/zenodo_release_metadata_20260527/zenodo_release_manifest.json")


def normalize_doi(value: str) -> str:
    doi = value.strip()
    for prefix in ("https://doi.org/", "http://doi.org/", "https://dx.doi.org/", "http://dx.doi.org/"):
        if doi.lower().startswith(prefix):
            doi = doi[len(prefix) :]
            break
    doi = doi.strip().rstrip("/")
    if not DOI_RE.fullmatch(doi):
        raise ValueError(f"Not a Zenodo DOI: {value!r}. Expected pattern 10.5281/zenodo.<record_id>.")
    return doi


def _replace_between_markers(text: str, start: str, end: str, replacement: str) -> str:
    pattern = re.compile(re.escape(start) + r".*?" + re.escape(end), flags=re.S)
    block = f"{start}\n{replacement.rstrip()}\n{end}"
    if pattern.search(text):
        return pattern.sub(block, text)
    suffix = "\n" if text.endswith("\n") else "\n\n"
    return text + suffix + block + "\n"


def update_citation_cff(path: Path, doi: str) -> bool:
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    changed = False
    for idx, line in enumerate(lines):
        if line.startswith("doi:"):
            new_line = f"doi: {doi}"
            if line != new_line:
                lines[idx] = new_line
                changed = True
            break
    else:
        lines.append(f"doi: {doi}")
        changed = True
    if changed:
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return changed


def update_release_manifest(path: Path, doi: str) -> bool:
    data = json.loads(path.read_text(encoding="utf-8"))
    before = json.dumps(data, sort_keys=True)
    data["zenodo_doi"] = doi
    data["zenodo_url"] = DOI_URL_PREFIX + doi
    data["doi_status"] = "doi_minted"
    data["required_before_citation"] = [
        "Zenodo DOI minted and recorded.",
        "Use the DOI URL in manuscript Code availability, CITATION.cff, README, and release notes.",
    ]
    after = json.dumps(data, sort_keys=True)
    if before != after:
        path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        return True
    return False


def update_data_upload_notes(path: Path, doi: str) -> bool:
    text = path.read_text(encoding="utf-8")
    new = re.sub(r"`zenodo_doi=RESULT_PENDING`", f"`zenodo_doi={doi}`", text)
    new = re.sub(
        r"Zenodo metadata are prepared with\s+`zenodo_doi=[^`]+`; no DOI should be cited until Zenodo or an\s+institutional archive mints a real identifier\.",
        f"Zenodo metadata are archived with DOI `{doi}` ({DOI_URL_PREFIX}{doi}).",
        new,
        flags=re.S,
    )
    if new != text:
        path.write_text(new, encoding="utf-8")
        return True
    return False


def update_readme(path: Path, doi: str) -> bool:
    text = path.read_text(encoding="utf-8")
    body = (
        f"Frozen release DOI: `{doi}` ({DOI_URL_PREFIX}{doi}).\n\n"
        "This DOI refers to the archived public software and shareable reproducibility materials. "
        "Controlled-access or licensed source datasets remain governed by their original repositories."
    )
    new = _replace_between_markers(text, "<!-- ZENODO_DOI_START -->", "<!-- ZENODO_DOI_END -->", body)
    if new != text:
        path.write_text(new, encoding="utf-8")
        return True
    return False


def update_manual_steps(path: Path, doi: str) -> bool:
    text = path.read_text(encoding="utf-8")
    new = text.replace("Status: RESULT_PENDING until Zenodo mints a real DOI.", f"Status: DOI minted as `{doi}`.")
    new = re.sub(r"After Zenodo publishes the record, copy the minted DOI\. It must match the pattern `10\.5281/zenodo\.<record_id>`\.", f"Zenodo DOI: `{doi}` ({DOI_URL_PREFIX}{doi}).", new)
    new = new.replace("The Zenodo DOI is not complete until the record is published by Zenodo. Until then, manuscript and citation files must retain `RESULT_PENDING` or omit the DOI.", f"The Zenodo DOI has been minted as `{doi}` and can be cited as {DOI_URL_PREFIX}{doi}.")
    if new != text:
        path.write_text(new, encoding="utf-8")
        return True
    return False


def update_code_availability(path: Path, doi: str) -> bool:
    text = path.read_text(encoding="utf-8")
    sentence = f"The archived release DOI is `{doi}` ({DOI_URL_PREFIX}{doi})."
    if "archived release DOI" in text:
        new = re.sub(r"The archived release DOI is `10\.5281/zenodo\.\d+` \(https://doi\.org/10\.5281/zenodo\.\d+\)\.", sentence, text)
    else:
        archive_pat = re.compile(r"(release-specific source archive is available at `https://github\.com/ahvsjags/EcoNiche-Opt/archive/refs/tags/v0\.3\.4-gpu-lipid-pair-rescue-20260528\.zip`\.)")
        new = archive_pat.sub(r"\1 " + sentence, text, count=1)
        if new == text:
            new = text.rstrip() + "\n\n" + sentence + "\n"
    if new != text:
        path.write_text(new, encoding="utf-8")
        return True
    return False


def apply_doi(root: Path, doi_value: str, dry_run: bool = False) -> list[Path]:
    doi = normalize_doi(doi_value)
    changed: list[Path] = []

    updates = [
        (Path("CITATION.cff"), update_citation_cff),
        (DEFAULT_MANIFEST, update_release_manifest),
        (Path("DATA_RESULTS_FIGURES_UPLOAD_NOTES.md"), update_data_upload_notes),
        (Path("README.md"), update_readme),
        (Path("deliverables/zenodo_manual_publication_steps_20260528.md"), update_manual_steps),
        (Path("paper/econiche_opt_manuscript_en_v1_20260509.md"), update_code_availability),
        (Path("paper/Journal of Translational Medicine投稿/EcoNiche-Opt_JTM_Main_Manuscript.md"), update_code_availability),
        (Path("paper/communications_medicine_submission/communications_medicine_submission_readiness.md"), update_code_availability),
    ]

    if dry_run:
        return [rel for rel, _ in updates if (root / rel).exists()]

    for rel, updater in updates:
        path = root / rel
        if path.exists() and updater(path, doi):
            changed.append(rel)
    return changed


def main() -> None:
    parser = argparse.ArgumentParser(description="Apply a minted Zenodo DOI to release citation files.")
    parser.add_argument("--doi", required=True, help="Zenodo DOI, e.g. 10.5281/zenodo.12345678")
    parser.add_argument("--root", default=str(ROOT))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    try:
        doi = normalize_doi(args.doi)
        changed = apply_doi(Path(args.root), doi, dry_run=args.dry_run)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    print(json.dumps({"doi": doi, "dry_run": args.dry_run, "targets": [str(p) for p in changed]}, ensure_ascii=False))


if __name__ == "__main__":
    main()
