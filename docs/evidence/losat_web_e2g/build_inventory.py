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
    "K-limit": "keep-rejected: LOSAT limit where NCBI's 32-bit arithmetic overflows (the greedy gap-cost limit 32767 is conservative: NCBI fails from about 600000, audit c)",
    "X-cli": "approved exception: argument syntax errors and -help (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES 1)",
    "X-threads": "approved exception: -num_threads with -subject, no NCBI thread warnings (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES 2)",
    "X-oom": "approved exception: allocation failure aborts (PD-LOSAT-CLI-NONSEARCH-DIFFERENCES 4)",
    "XD1": "approved exception: a query chunk that NCBI would split again is searched once (PD-LOSAT-NCBI-DEFECTS 1)",
    "XD2": "approved exception: outfmt 0 titles of punctuation stop at the end of the string (PD-LOSAT-NCBI-DEFECTS 2)",
    "P1": "port: -evalue with a signed infinity or NaN, which NCBI searches with (PD-LOSAT-NCBI-DEFECTS rule 4)",
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
    (r"B", r"SplitQuery_CreateChunkData", r"", r"divergent", "XD1"),
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
    # Maintainer decision 2026-10-02 (PD-LOSAT-NCBI-DEFECTS): the infinite and NaN e-values
    # that NCBI accepts are ported.
    (r"A", r"CArg_Double", r"NaN", r"rejected", "P1"),
    (r"A", r"CArg_Double", r"", r"rejected", "K-layer"),
    # Audit (c), 2026-10-02: NCBI decodes HTML character references (NStr::HtmlDecode) in the
    # titles, it does not crash; only the titles of punctuation crash it (x_CleanAndCompress).
    (r"F", r"GenerateDefline", r"HTML|Html|html|&", r"rejected", "K-layer"),
    # Maintainer decision 2026-10-02: the titles of punctuation are an approved exception.
    (r"F", r"GenerateDefline|x_CleanAndCompress", r"", r"rejected", "XD2"),
    (r".*", r"s_MultiSeqGetTotLen|s_BlastGreedyAlignMemAlloc", r"", r"rejected", "K-limit"),
    (r".*", r".*", r"rmblastn|reward == 0|reward <= 0|matrix_only|m_DisableKAStats|-ungapped|is_gapped false|"
                 r"m_IsUngappedSearch|blastn-short|dc-megablast", r"rejected", "K-mode"),
    (r".*", r"CreateTask|GetTasks|CTaskCmdLineArgs", r"", r"rejected", "K-mode"),
    (r"B|C", r"SetupSubjects_OMF|BLAST_GapAlignSetUp", r"", r"rejected", "K-fasta"),
    (r".*", r".*", r"", r"rejected", "K-layer"),
]

# Filled as the transpile items are finished: action -> result text (S07+++b, 2026-10-02).
RESULTS: dict[str, str] = {
    "T1": "faithful after 019ba8ea4; concatenated q_start, stable sort; fixtures pal.* (40 X+revcomp(X) queries) match NCBI before and after",
    "T2": "faithful after b5940fba9 (+ fixture de42604ff); Int4 offset; overflow hunt (2762 commands, overflow-checked build): the only overflow site, 188 commands, all equal NCBI before and after; 926 sweep combinations equal",
    "T3": "faithful after 35e3e6d7f; !(evalue > cutoff) in both reaps; unit test; reachable with the NaN e-values ported in 30713884f (P1)",
    "T4": "faithful after 7eb018f70; traced in the stored order (heap lists in e-value order); for BLASTN equal to score order; one interval tree equivalent to per-query trees; fixtures prelim.gaps10_*, ties.gaps10_1_1",
    "T5": "faithful after 435f97afd; ascending cells for the small and standard tables (TaskConfig::mb_lookup); unit test; fixtures rep.*, rep2.*",
    "T6": "faithful after 0b3b851c7; diagonal array whenever the block is at most 8000; fixtures sq.* (12 queries, block 5155)",
    "T7": "faithful after 8ccf079ba, 5dac71f72 and ad5fa9c85; NCBI's StringToInt and size_t/TSeqPos arithmetic; CHUNK_SIZE=1000 without BATCH_SIZE gives NCBI's empty-batch error after the outfmt 0 prolog (exit 3, reproduced by maintainer decision, PD-LOSAT-NCBI-DEFECTS); a chunk that NCBI would split again (CCoreException, exit 3) is searched once, approved exception 1 of PD-LOSAT-NCBI-DEFECTS after d846e9bbe (fixtures env.resplit_*, validation resplit/); explicit rejection: non-integers (CStringException text names the build's files, exit 255) and a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE where it splits the batch (NCBI's size_t chunk ranges wrap: a CCoreException, or chunks with gaps); where such a pair does not split, as NCBI (fixtures env.negative_pair_*); 28 fixtures env.*; 22 batch_sweep runs x 60 cases and 4 split_check runs with the variables: 0 differ",
    "T8": "faithful after baf180fbb; reports of the batches before, no epilog, then BLAST engine error (exit 3); a failing batch with an invalid query (NCBI crashes) stays rejected; fixtures kaerror.later_batch.fmt{0,6,7}",
    "T9": "faithful after baf180fbb and e37099f44; title warnings when a batch is read, invalid-query warnings with its report, written between the query reports as NCBI posts them on cerr (tied to cout: before the report of the batch's first query and before the query's preamble); fixtures warnings.batches.fmt{0,6}, warnings.batch1000.fmt7 and *.merged (2>&1); order/: 456 runs (19 inputs, outfmt 0/6/7, merged, separate, /dev/full; BATCH_SIZE unset, 1, 200) equal NCBI",
    "T10": "faithful after 1a0fd98c1; -subject optional for the parser, NCBI's error (exit 1) before -query and -out are opened; fixtures nosubject.*",
    "T11": "faithful after e4b5c4a63; kBits from CTOOLKIT_COMPATIBLE (any value); fixtures ctoolkit.*; 55 outfmt 0 cases of every program x 3 environments: 51 same, 4 TBLASTX (outfmt 0 not implemented, S08)",
    "T12": "faithful after 650b02771; integers accepted (no output change), non-integers rejected explicitly (CStringException text names the build's files, exit 255); fixtures env.prefetch*",
    "T13": "faithful after a99527f20; exact %#8.3g; unit test with C's strings",
    "T14": "faithful after 5cd9cc3cf and 48a9ee0c2; outfmt 0 write failure: BLAST failed to write output, exit 6 (oracle -out /dev/full); a closed pipe in outfmt 0 (NCBI: ended by SIGPIPE) gives the same message and exit 6, approved exception 5 of PD-LOSAT-CLI-NONSEARCH-DIFFERENCES 1.1 (DW-14); after e37099f44 a failed outfmt 0 write stops at the flush before the first warning, as NCBI's first flush after the version line (no warning is printed); fixtures write.devfull_fmt0, write.devfull_warnings_fmt0",
    "XD1": "approved exception after d846e9bbe (PD-LOSAT-NCBI-DEFECTS 1): each chunk searched once; validation resplit/: 42 configurations, NCBI exit 3 in all, LOSAT equals NCBI's unsplit output in 35 and NCBI's output at the largest overlap without a second split in 74 of 84 runs, the others differ only by chunk-boundary HSPs as NCBI's own splits; fixtures env.resplit_* (oracle_env)",
    "XD2": "approved exception after c72452236 (PD-LOSAT-NCBI-DEFECTS 2): the cleanup stops at the end of the string (', ,' gives ', '); NCBI crashes (SIGSEGV) in all 4 validation runs with hits, LOSAT's report equals NCBI's with placeholder deflines once the placeholders are replaced (punct_defline/); fixtures punct.nohit_* (no hits on such subjects: NCBI runs) match NCBI",
    "P1": "ported after 30713884f: the limit is removed, ncbi_double also reads nan(n-char-sequence), the web API reads -evalue with the CLI's parser; the cutoff stays 1 and no HSP is reaped (unit test); fixtures evalue.*; sweep e2g-sweep-AC5: 440 cases (+inf, -nan, +nan(1), 1e999; megablast, blastn, word 7, -subject_besthit, IUPAC pool), 0 differ",
    "R1": "explicit rejection after 7b63980b2 (any value; an empty query still ends with Query is Empty!)",
    "R2": "explicit rejection after 219c2c49e: DIAG_*, NCBI_CONFIG_* (except entries that change no output), ABORT_ON_THROW, stack-trace and LOG_* parameters, non-Boolean BLAST_USAGE_REPORT; blastn.ini and .ncbirc on NCBI's search path accepted only with entries that change no output (fixture ncbirc.harmless_fmt0); oracle research e2g-r2r3",
    "R3": "explicit rejection after 6182cef73 (-subject, -query, -out not UTF-8, in NCBI's open order)",
    "V1": "rejection removed after 82e8c7593 except reward 32767 / penalty -32768 (BlastScoreBlkMaxScoreSet leaves BLAST_SCORE_MAX/MIN out of the range: a reward of 32767 is counted outside the frequency array and can crash NCBI, a penalty of -32768 makes every query invalid); sweep e2g-v1-sweep: 512 compared, 448 identical, the 64 others reward 32767 or megablast gap costs above 32767 (K-limit); fixtures scores.*",
}


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
