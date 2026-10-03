#!/usr/bin/env python3
"""Build INVENTORY.tsv: the SD inventory of NCBI's dc-megablast and blastn-short path and the port.

inventory/<R>.tsv are the rows of range R (A application and arguments, B the C++ options
layer and the option checks, C core setup and parameters, D lookup/scan/ungapped
extension, E gapped extension/traceback/hit saving, F formatting), written before the port
by read-only agents (conventions in inventory/COMMON.md); `status` is LOSAT before the
port (a92fa902f). inventory/result_<R>.tsv give, per row, what the code of the port does
(`e2i_result`, conventions in inventory/RESULT_COMMON.md), found by a second set of
read-only agents on 90c5f0181; their extra rows (X1, X2, ...) are NCBI functions the
inventory missed. RESOLUTION records how a GAP row was then fixed (`e2i_final`); every
other row keeps its result.

Usage: build_inventory.py  (writes INVENTORY.tsv next to this script and prints the counts)
"""
from __future__ import annotations

import collections
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
RANGES = "ABCDEF"
FIELDS = ["range", "row", "ncbi_file", "ncbi_lines", "ncbi_function", "branch", "status_before", "e2i_result",
          "e2i_final", "losat_location", "evidence", "notes"]
TEMPLATE_LENGTH = ("ported: -template_length is checked against 16/18/21 with a base-10 conversion, as "
                   "CArgAllowIntegerSet (blast_input_aux.hpp:214-222,239), so 0x12 is a parser error "
                   "(blastinput/value_parsers.rs blastn_template_length; tests/cli_v2.rs), commit {fix}")
# (range, row) -> what became of a GAP row after the result agents' report.
RESOLUTION = {
    ("A", "27"): TEMPLATE_LENGTH,
    ("A", "X1"): TEMPLATE_LENGTH,
}
FIX = "c45e57d85"


def read(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))


def main() -> int:
    rows = []
    for name in RANGES:
        results = read(HERE / "inventory" / f"result_{name}.tsv")
        inventory = read(HERE / "inventory" / f"{name}.tsv")
        assert len([r for r in results if not r["row"].startswith("X")]) == len(inventory), name
        for result in results:
            key = (name, result["row"])
            final = result["e2i_result"]
            if final == "GAP":
                final = RESOLUTION[key].format(fix=FIX)
            rows.append({
                "range": name, "row": result["row"], "ncbi_file": result["ncbi_file"],
                "ncbi_lines": result["ncbi_lines"], "ncbi_function": result["ncbi_function"],
                "branch": result["branch"], "status_before": result["status_before"],
                "e2i_result": result["e2i_result"], "e2i_final": final,
                "losat_location": result["losat_location"], "evidence": result["evidence"],
                "notes": result["notes"],
            })
    with open(HERE / "INVENTORY.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_NONE,
                                escapechar="\\")
        writer.writeheader()
        writer.writerows(rows)
    before = collections.Counter(row["status_before"] for row in rows)
    after = collections.Counter(row["e2i_result"] for row in rows)
    final = collections.Counter(row["e2i_final"].split(":")[0] for row in rows)
    print(f"{len(rows)} rows; before {dict(before)}; result {dict(after)}; final {dict(final)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
