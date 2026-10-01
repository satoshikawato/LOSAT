#!/usr/bin/env python3
"""Known frozen raw-output mismatches (frozen_mismatch_allowlist.json).

The frozen expectations are never rewritten. A listed case passes only with the
listed actual output; a listed case that matches its frozen hash again (or that a
run covering its program never executes) is reported as stale so the entry is
removed.
"""
from __future__ import annotations

import json
from pathlib import Path

ALLOWLIST = Path(__file__).with_name("frozen_mismatch_allowlist.json")
SCHEMA = "losat-frozen-mismatch-allowlist-v1"


class Allowlist:
    def __init__(self, path: Path = ALLOWLIST):
        document = json.loads(Path(path).read_text())
        if document.get("schema") != SCHEMA:
            raise ValueError(f"{path}: unexpected schema {document.get('schema')!r}")
        self.path = Path(path)
        self.entries = {(entry["program"], entry["case_id"]): entry for entry in document["entries"]}

    # NCBI reference: c++/src/objtools/align_format/tabular.cpp:1098-1108
    # x_PrintField(*iter); m_Ostream << "\n";
    # Raw output bytes are compared by SHA-256 without normalization.
    def classify(self, program: str, case_id: str, expected: str, actual: str | None, executed: bool = True) -> str:
        """Returns 'match', 'allowed', 'stale' or 'mismatch'."""
        entry = self.entries.get((program, case_id))
        equal = executed and actual == expected
        if entry is None or entry["frozen_sha256"] != expected:
            return "match" if equal else "mismatch"
        if equal:
            return "stale"
        if executed and actual == entry["allowed_actual_sha256"]:
            return "allowed"
        return "mismatch"

    def unexecuted(self, executed: set[tuple[str, str]], programs: set[str]) -> list[tuple[str, str]]:
        """Entries of the given programs that were not executed."""
        return sorted(key for key in self.entries if key[0] in programs and key not in executed)
