#!/usr/bin/env python3
"""Builds INVENTORY.tsv from the stage-2 range tables and the orchestrator's decisions.

Every row of stage2/{A,B,C,D,E1,E2,F}.tsv gets an id (<range>-<n>) and an `e2g_action`:
`-` for faithful and n/a rows; for the other rows the action of the first matching rule
below (see stage2/REVIEW.md for the checks behind each decision). `e2g_result` records
what the transpile did (filled from RESULTS as items are finished). The build fails if a
divergent or unported row has no rule.

Usage: build_inventory.py   (writes INVENTORY.tsv next to this script)
"""
from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
RANGES = ("A", "B", "C", "D", "E1", "E2", "F")
FIELDS = ["id", "range", "ncbi_file", "ncbi_lines", "ncbi_function", "branch", "losat_location",
          "status", "impact", "notes", "e2g_action", "e2g_result"]

# Transpile items (T), new explicit rejections (R), re-examinations (V), kept rejections
# (K) and approved exceptions (X, PD-LOSAT-CLI-NONSEARCH-DIFFERENCES, 2026-10-02).
ACTIONS = {
    "T1": "transpile: initial hit order (concatenated query offset, stable sort)",
    "T2": "transpile: BlastGetStartForGappedAlignmentNucl in Int4 arithmetic",
    "T3": "transpile: preliminary e-value reap comparison as in C",
    "T4": "transpile: traceback in the stored HSP list order",
    "T5": "transpile: ascending query offsets in small and standard lookup cells",
    "T6": "transpile: diagonal array whenever the query block is at most 8000",
    "T7": "transpile: BATCH_SIZE, CHUNK_SIZE, OVERLAP_CHUNK_SIZE",
    "T8": "transpile: Karlin-Altschul failure after a first batch of invalid queries",
    "T9": "transpile: warnings per query batch",
    "T10": "transpile: missing -subject error",
    "T11": "transpile: CTOOLKIT_COMPATIBLE (bits)",
    "T12": "transpile: PRE_FETCH_SEQS_LIMIT error",
    "T13": "transpile: exact %#8.3g for Lambda, K, H",
    "T14": "transpile: outfmt 0 write failure (BLAST failed to write output, exit 6); outfmt 6/7 is X-write",
    "R1": "new explicit rejection: BL2SEQ_LEGACY",
    "R2": "new explicit rejection: NCBI diagnostics/registry environment variables and .ncbirc keys that change output",
    "R3": "new explicit rejection: file names that are not UTF-8 (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES)",
    "V1": "re-examine: reward - penalty above 3000",
    "K-fasta": "keep-rejected until session SF ports CFastaReader before S17 (TD-12, plan §10, DW-13)",
    "K-layer": "keep-rejected: needs an NCBI C++ layer or feature outside the app's option scope",
    "K-crash": "keep-rejected: NCBI crashes or misbehaves",
    "K-mode": "keep-rejected: another program mode (rmblastn scoring, other tasks, -ungapped)",
    "K-limit": "keep-rejected: LOSAT limit where NCBI's 32-bit arithmetic overflows",
    "X-cli": "approved exception: argument syntax errors and -help (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES 1)",
    "X-threads": "approved exception: -num_threads with -subject, no NCBI thread warnings (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES 2)",
    "X-oom": "approved exception: allocation failure aborts (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES 4)",
}

# (range regex, function regex, branch/notes regex, status regex) -> action
RULES = [
    (r"D|E1", r"score_compare_match|Blast_InitHitListSortByScore|BLAST_GetGappedScore", r"", r"divergent", "T1"),
    (r"E1", r"BlastGetStartForGappedAlignmentNucl", r"", r"divergent", "T2"),
    (r"E1", r"s_Blast_HSPListReapByPrelimEvalue", r"", r"divergent", "T3"),
    (r"E2", r"Blast_TracebackFromHSPList", r"", r"divergent", "T4"),
    (r"D", r"s_BlastSmallNaLookupFinalize|BlastSmallNaLookupTableNew|s_BlastNaLookupFinalize|BlastNaLookupTableNew", r"", r"divergent", "T5"),
    (r"D", r"BlastExtendWordNew", r"", r"divergent", "T6"),
    (r".*", r"SplitQuery_GetOverlapChunkSize|SplitQuery_GetChunkSize|GetQueryBatchSize", r"", r".*", "T7"),
    (r"B|C", r"CreateScoreBlock|Blast_ScoreBlkKbpGappedCalc", r"invalid|later batch", r"rejected", "T8"),
    (r"A", r"CreateWarningsForSeqDataInTitle", r"timing", r"divergent", "T9"),
    (r"A", r"CBlastDatabaseArgs::ExtractAlgorithmOptions", r"neither -db nor -subject", r"divergent", "T10"),
    (r"A|F", r"kBits", r"", r".*", "T11"),
    (r"A", r"s_PreFetchSeqs", r"", r".*", "T12"),
    (r"F", r"PrintKAParameters", r"", r"divergent", "T13"),
    (r".*", r".*", r"BL2SEQ_LEGACY", r"unported", "R1"),
    (r"A", r"CNcbiApplication environment", r"", r"unported", "R2"),
    (r"C", r"Blast_KarlinBlkUngappedCalc", r"3000", r"rejected", "V1"),
    (r"A", r"CBlastnApp::Init|CNcbiApplicationAPI::AppMain", r"", r"divergent", "X-cli"),
    (r"A", r"CMTArgs", r"", r"unported", "X-threads"),
    (r"A", r"GetSubjectFile", r"", r"divergent", "R3"),
    (r"A", r"CATCH_ALL", r"BLAST failed to write output", r"divergent", "T14"),
    (r"A", r"CATCH_ALL", r"Out of memory", r"unported", "X-oom"),
    (r"A", r"CFasta|FastaDefline|IsIStreamEmpty|x_FastaToSeqLoc|ReadOneSeq|GetNextSeqBatch|CATCH_ALL", r"", r"rejected", "K-fasta"),
    (r"A", r"CArg_Double", r"", r"rejected", "K-layer"),
    (r"F", r"GenerateDefline|x_CleanAndCompress", r"", r"rejected", "K-crash"),
    (r".*", r"s_MultiSeqGetTotLen|s_BlastGreedyAlignMemAlloc", r"", r"rejected", "K-limit"),
    (r".*", r".*", r"rmblastn|reward == 0|reward <= 0|matrix_only|m_DisableKAStats|-ungapped|is_gapped false|"
                 r"m_IsUngappedSearch|blastn-short|dc-megablast", r"rejected", "K-mode"),
    (r".*", r"CreateTask|GetTasks|CTaskCmdLineArgs", r"", r"rejected", "K-mode"),
    (r"B|C", r"SetupSubjects_OMF|BLAST_GapAlignSetUp", r"", r"rejected", "K-fasta"),
    (r".*", r".*", r"", r"rejected", "K-layer"),
]

# Filled as the transpile items are finished: action -> result text.
RESULTS: dict[str, str] = {}


def main() -> int:
    rows = []
    for name in RANGES:
        with open(HERE / "stage2" / f"{name}.tsv", newline="") as handle:
            for index, row in enumerate(csv.DictReader(handle, delimiter="\t"), 1):
                row = {key: (value or "").strip() for key, value in row.items() if key}
                row["id"] = f"{name}-{index:03d}"
                action = "-"
                if row["status"] not in ("faithful", "n/a"):
                    text = row["branch"] + " " + row["notes"]
                    for rng, function, branch, status, act in RULES:
                        if (re.fullmatch(rng, row["range"]) and re.search(function, row["ncbi_function"])
                                and re.search(branch, text) and re.fullmatch(status, row["status"])):
                            action = act
                            break
                    else:
                        print(f"no rule: {row['id']} {row['status']} {row['ncbi_function']}", file=sys.stderr)
                        return 1
                row["e2g_action"] = action
                row["e2g_result"] = RESULTS.get(action, "")
                rows.append(row)
    with open(HERE / "INVENTORY.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} rows")
    print("status:", dict(Counter(row["status"] for row in rows)))
    counts = Counter(row["e2g_action"] for row in rows if row["e2g_action"] != "-")
    for action, count in sorted(counts.items()):
        print(f"{action}\t{count}\t{ACTIONS[action]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
