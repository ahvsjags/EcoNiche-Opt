# Zenodo publication steps for EcoNiche-Opt v0.3.4

Status: RESULT_PENDING until Zenodo mints a real DOI.

## Frozen GitHub release

- Repository: https://github.com/ahvsjags/EcoNiche-Opt
- Release tag: v0.3.4-gpu-lipid-pair-rescue-20260528
- Release URL: https://github.com/ahvsjags/EcoNiche-Opt/releases/tag/v0.3.4-gpu-lipid-pair-rescue-20260528
- Verified tag commit: f7d8dc0ca3fbd42224e5d21b7fdec8f182247883
- Version: 0.3.4

## Recommended Zenodo route

1. Log in to https://zenodo.org with the account that will own the archive.
2. Open GitHub integration: https://zenodo.org/account/settings/github/
3. Enable archiving for `ahvsjags/EcoNiche-Opt`.
4. Confirm that the release `v0.3.4-gpu-lipid-pair-rescue-20260528` is archived.
5. If Zenodo does not archive the existing release automatically, create a new GitHub release using the same tag after integration is enabled, or use the manual upload route below.
6. After Zenodo publishes the record, copy the minted DOI. It must match the pattern `10.5281/zenodo.<record_id>`.
7. Replace `RESULT_PENDING` in `CITATION.cff`, manuscript Code availability, README/release notes, and `deliverables/zenodo_release_metadata_20260527/zenodo_release_manifest.json`.
8. Rerun:

```powershell
python scripts\validation\audit_top_tier_targets.py
python scripts\validation\audit_submission_package.py
python -m econiche_opt.cli validate-goals --goal-file docs\goal_status.yml
```

## API publication route

If you prefer a scripted route, create a Zenodo access token with `deposit:write` and `deposit:actions` scopes and set it only as an environment variable:

```powershell
$env:ZENODO_TOKEN='PASTE_ZENODO_TOKEN_HERE'
python scripts\reporting\publish_zenodo_release.py --execute --publish --apply-doi
```

The script builds the archive from the frozen GitHub release tag, creates a Zenodo deposition, uploads the zip archive, publishes the record, extracts the minted DOI, and then calls `scripts\reporting\apply_zenodo_doi.py` to update citation/manuscript files. Without `--execute`, the same script runs a dry-run and writes `deliverables/zenodo_api_publication_status_20260528.json`; without `--publish`, it only creates an unpublished draft and must not be cited as a DOI-bearing record.

To test the flow without publishing to production Zenodo:

```powershell
$env:ZENODO_TOKEN='PASTE_SANDBOX_ZENODO_TOKEN_HERE'
python scripts\reporting\publish_zenodo_release.py --sandbox --execute
```

## Manual upload fallback

If GitHub integration is unavailable, create a new Zenodo software upload and upload the source archive from:

https://github.com/ahvsjags/EcoNiche-Opt/archive/refs/tags/v0.3.4-gpu-lipid-pair-rescue-20260528.zip

Use metadata from:

- `deliverables/zenodo_release_metadata_20260527/.zenodo.json`

Do not upload controlled-access raw data, licensed data, local caches, or expression matrices excluded by `DATA_RESULTS_FIGURES_UPLOAD_NOTES.md`.

## Claim boundary

The GitHub release is complete and public. The Zenodo DOI is not complete until the record is published by Zenodo. Until then, manuscript and citation files must retain `RESULT_PENDING` or omit the DOI.
