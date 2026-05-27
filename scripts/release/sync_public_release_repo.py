from __future__ import annotations

import argparse
import fnmatch
import hashlib
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]

ROOT_FILES = [
    ".gitignore",
    "AGENTS.md",
    "CITATION.cff",
    "DATA_RESULTS_FIGURES_UPLOAD_NOTES.md",
    "Dockerfile",
    "environment.yml",
    "LICENSE",
    "Makefile",
    "MANIFEST.in",
    "pyproject.toml",
    "pytest.ini",
    "README.md",
    "requirements.txt",
]

PAPER_ROOT_FILES = [
    "article_figure_captions_20260508.md",
    "article_suite_audit_20260508.md",
    "econiche_opt_cover_letter_en_v1_20260509.docx",
    "econiche_opt_cover_letter_en_v1_20260509.md",
    "econiche_opt_manuscript_en_v1_20260509.docx",
    "econiche_opt_manuscript_en_v1_20260509.md",
    "econiche_opt_supporting_en_v1_20260509.docx",
    "econiche_opt_supporting_en_v1_20260509.md",
    "manuscript_econiche_opt_article_draft_20260508.md",
    "references.bib",
]

FULL_DIRS = [
    ".github",
    "config",
    "docs",
    "r-package",
    "schemas",
    "scripts",
    "src",
    "tests",
    "workflow",
]

SELECTED_DIRS = [
    "data/metadata",
    "data/priors",
    "data/interim",
    "data/processed",
    "deliverables",
    "figures/article",
    "paper/communications_medicine_submission",
    "paper/Journal of Translational Medicine投稿",
    "results",
    "tables",
]

EXCLUDE_NAMES = {
    "__pycache__",
    ".pytest_cache",
    ".snakemake",
    ".git",
    "tmp",
    "build",
    "dist",
}

EXCLUDE_PATTERNS = [
    "*.pyc",
    "*.part",
    "*.tiff",
    "*.tif",
    "*.zip",
    "*.expr.tsv",
    "*expression*.tsv",
    "goal_status.yml",
    "github_release_status_*.json",
    "top_tier_target_audit_*.md",
    "top_tier_target_audit_*.tsv",
    "zenodo_api_publication_status_*.json",
    "zenodo_upload_archive_*/*",
    "deliverables/github_release_status_*.json",
    "deliverables/top_tier_target_audit_*.md",
    "deliverables/top_tier_target_audit_*.tsv",
    "deliverables/zenodo_api_publication_status_*.json",
    "deliverables/zenodo_upload_archive_*/*",
    "zenodo_release_metadata_*/ZENODO_RELEASE_CHECKLIST.md",
    "zenodo_release_metadata_*/zenodo_release_manifest.json",
    "deliverables/zenodo_release_metadata_*/ZENODO_RELEASE_CHECKLIST.md",
    "deliverables/zenodo_release_metadata_*/zenodo_release_manifest.json",
    "*/processed_inputs/*",
    "*/tmp/*",
]

TARGET_CLEANUP_PATTERNS = [
    "deliverables/github_release_status_*.json",
    "deliverables/top_tier_target_audit_*.md",
    "deliverables/top_tier_target_audit_*.tsv",
    "deliverables/zenodo_release_metadata_*/ZENODO_RELEASE_CHECKLIST.md",
    "deliverables/zenodo_release_metadata_*/zenodo_release_manifest.json",
    "docs/goal_status.yml",
]


def should_skip(path: Path, source_root: Path) -> bool:
    rel = path.relative_to(source_root).as_posix()
    if any(part in EXCLUDE_NAMES for part in path.parts):
        return True
    return any(fnmatch.fnmatch(rel, pattern) or fnmatch.fnmatch(path.name, pattern) for pattern in EXCLUDE_PATTERNS)


def copy_file(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def cleanup_target_only_files(target: Path) -> list[Path]:
    removed: list[Path] = []
    for pattern in TARGET_CLEANUP_PATTERNS:
        for path in target.glob(pattern):
            if path.is_file():
                path.unlink()
                removed.append(path)
    return removed


def copy_tree_filtered(source_root: Path, target_root: Path) -> list[Path]:
    copied: list[Path] = []
    for src in source_root.rglob("*"):
        if src.is_dir() or should_skip(src, source_root):
            continue
        rel = src.relative_to(source_root)
        dst = target_root / rel
        copy_file(src, dst)
        copied.append(dst)
    return copied


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", default="github_release/econiche-opt-package-20260508")
    parser.add_argument("--manifest", default="public_release_sync_manifest.tsv")
    args = parser.parse_args()

    target = (ROOT / args.target).resolve()
    if not target.exists() or not (target / ".git").exists():
        raise SystemExit(f"Target is not a Git release repository: {target}")
    if not str(target).startswith(str(ROOT.resolve())):
        raise SystemExit(f"Refusing to sync outside project root: {target}")

    removed = cleanup_target_only_files(target)

    copied: list[Path] = []
    for rel in ROOT_FILES:
        src = ROOT / rel
        if src.exists() and src.is_file():
            copy_file(src, target / rel)
            copied.append(target / rel)

    for rel in FULL_DIRS:
        src = ROOT / rel
        if src.exists() and src.is_dir():
            copied.extend(copy_tree_filtered(src, target / rel))

    data_readme = ROOT / "data" / "README.md"
    if data_readme.exists():
        copy_file(data_readme, target / "data" / "README.md")
        copied.append(target / "data" / "README.md")

    external_target = target / "data" / "external"
    external_target.mkdir(parents=True, exist_ok=True)
    for src in (ROOT / "data" / "external").glob("*.tsv"):
        copy_file(src, external_target / src.name)
        copied.append(external_target / src.name)
    external_readme = ROOT / "data" / "external" / "README.md"
    if external_readme.exists():
        copy_file(external_readme, external_target / "README.md")
        copied.append(external_target / "README.md")

    for rel in SELECTED_DIRS:
        src = ROOT / rel
        if src.exists() and src.is_dir():
            copied.extend(copy_tree_filtered(src, target / rel))

    for filename in PAPER_ROOT_FILES:
        src = ROOT / "paper" / filename
        if src.exists() and src.is_file():
            dst = target / "paper" / filename
            copy_file(src, dst)
            copied.append(dst)

    rows = ["relative_path\tsize_bytes\tsha256"]
    for path in sorted(set(copied)):
        rel = path.relative_to(target).as_posix()
        rows.append(f"{rel}\t{path.stat().st_size}\t{sha256(path)}")
    manifest_path = target / args.manifest
    manifest_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"Copied {len(set(copied))} public files into {target}")
    if removed:
        print(f"Removed {len(set(removed))} post-release verification files from {target}")
    print(f"Wrote {manifest_path}")


if __name__ == "__main__":
    main()
