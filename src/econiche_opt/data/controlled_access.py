from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


CONTROLLED_TOKENS = ("dbgap", "ega", "controlled", "restricted", "access_restricted")


@dataclass(frozen=True)
class ControlledAccessDecision:
    status: str
    reason: str
    instructions: str


def is_controlled_access(access_text: str | None) -> bool:
    value = (access_text or "").lower()
    return any(token in value for token in CONTROLLED_TOKENS)


def controlled_access_decision(accession: str, access_text: str | None) -> ControlledAccessDecision:
    if is_controlled_access(access_text):
        return ControlledAccessDecision(
            status="ACCESS_RESTRICTED",
            reason=f"{accession} is marked as controlled or restricted access.",
            instructions="Do not fabricate data. Obtain approved access, then place real files under data/raw/controlled/.",
        )
    return ControlledAccessDecision(
        status="PUBLIC_OR_MANUAL",
        reason=f"{accession} is not marked as controlled in the registry.",
        instructions="Proceed only through the registered downloader or documented manual import path.",
    )


def assert_no_controlled_plaintext(path: str | Path) -> None:
    text = Path(path).read_text(encoding="utf-8") if Path(path).exists() else ""
    forbidden = ["password=", "token=", "dbgap_key", "ega_password"]
    hits = [item for item in forbidden if item.lower() in text.lower()]
    if hits:
        raise ValueError(f"Controlled-access credential material found in {path}: {', '.join(hits)}")
