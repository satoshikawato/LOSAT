#!/usr/bin/env python3
"""Build INVENTORY.tsv: the S08 inventory of NCBI's tblastx outfmt 0/6/7 path and the port.

inventory/<R>.tsv are the rows of range R (A application and orchestration, B alignments,
C description table, D outfmt 7/6, E what the report receives from the search, F residues
and genetic codes, G the hit list), written before the port by read-only agents
(conventions in inventory/COMMON.md); `status` is LOSAT's TBLASTX before S08.
inventory/result_<R>.tsv give, per row, what the final code of the port did
(`s08_result`, conventions in inventory/RESULT_COMMON.md), found by a second set of
read-only agents on the commit 0d533ba76. The rows they marked GAP were then fixed or
rejected; RESOLUTION below records how (`s08_final`). Rows E X1 and X2 were added by the
result agents. The independent audit of 968fa98c8 (README, "独立監査") changed the final
state of some rows (AUDIT_FINAL) and found NCBI paths without a row, which range H adds
(AUDIT_ROWS; status before S08 `-`). Its fixes are the commits 7fbbfad96 (round 1),
aa4b6f8c3, 6017f0ea2 and 5fa53b3f0 (round 2), and bc521f450 (round 3).

Usage: build_inventory.py  (writes INVENTORY.tsv next to this script and prints the counts)
"""
from __future__ import annotations

import collections
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent
RANGES = "ABCDEFG"
FIELDS = ["range", "row", "ncbi_file", "ncbi_lines", "ncbi_function", "branch", "status_before", "s08_result",
          "s08_final", "losat_location", "evidence", "notes"]
INPUTS = ("rejected: query and subject read as BLASTN reads them (blastn/input.rs check_deflines_of, "
          "check_residues_of, check_records_have_residues_of, check_utf8_file_name), 06953c68d")
# (range, row) -> what became of a GAP row after the result agents' report.
RESOLUTION = {
    ("A", "6"): INPUTS + " (empty deflines and deflines starting with white space, in every format)",
    ("A", "26"): INPUTS + " (residues that are not IUPAC nucleotide letters)",
    ("A", "41"): INPUTS + " (a file name that is not UTF-8)",
    ("A", "52"): "ported: the query is opened, -out created, then the query read (run_impl.rs run, "
                 "search_cli), 06953c68d",
    ("A", "116"): "ported: the threshold written as a C++ stream writes a double (report/pairwise.rs "
                  "cpp_default_double), 06953c68d",
    ("A", "126"): INPUTS + " (control characters, non-ASCII bytes, found in the bytes of the file)",
    ("D", "2"): "rejected: an empty query from a stream without a position (run_impl.rs search_cli), "
                "06953c68d; -query - stays an argument error (deferred)",
    ("D", "17"): INPUTS + " (the deflines NCBI reads differently; Subject_N ids not reached)",
    ("D", "37"): INPUTS,
    ("E", "29"): "ported: the init hit list sorted by score_compare_match per subject chunk "
                 "(blast_gapalign.rs sort_init_hsps_by_score_ncbi), 06953c68d; fixture tie.init_seg_no",
    ("E", "X1"): "rejected: -culling_limit above 0 (run_impl.rs search), 06953c68d",
    ("E", "X2"): "ported: BLAST_LargeGapSumE in NCBI's order of evaluation (stats/sum_statistics.rs "
                 "ncbi_large_gap_sum_e), 06953c68d; fixture sume.large_gap_seg_no",
    ("F", "13"): INPUTS + " (residues that are not IUPAC nucleotide letters)",
    ("F", "27"): "ported: U read as T before the search (blastn/input.rs with_u_as_t), 06953c68d; "
                 "fixtures input.rna_*",
    ("G", "32"): "rejected: -culling_limit above 0 (run_impl.rs search), 06953c68d",
}
AUDIT = "7fbbfad96"
# (range, row) -> the final state after the audit's fixes.
AUDIT_FINAL = {
    ("A", "24"): f"ported: -evalue 0 fails with NCBI's hit saving check (run_impl.rs check_hit_saving_options), "
                 f"{AUDIT}; deferred to S08+ (TD-13): NCBI's texts and exit 1 of the other option checks after "
                 "parsing (-threshold 0, malformed -seg, -word_size), which LOSAT's parser stops with exit 2",
    ("E", "50"): f"ported: a SEG window, locut or hicut that is not above 0 keeps NCBI's default "
                 f"(blast_filter.c:1147-1154; value_parsers.rs SegSpec::params), {AUDIT}; fixtures options.seg_*",
    ("A", "69"): "rejected: a record without residues (blastn/input.rs check_records_have_residues_of), 06953c68d",
    ("D", "51"): "rejected: a record without residues (blastn/input.rs check_records_have_residues_of), 06953c68d",
}
# Rows of range H: (ncbi_file, ncbi_lines, ncbi_function, branch, s08_final, losat_location, evidence).
AUDIT_ROWS = [
    ("src/algo/blast/core/blast_options.c", "1518-1523", "BlastHitSavingOptionsValidate", "-evalue <= 0",
     f"ported: {AUDIT}", "tblastx/blast_engine/run_impl.rs check_hit_saving_options",
     "fixtures options.evalue0_*; adapter validate"),
    ("src/algo/blast/core/blast_options.c", "1303-1311 and the other validators", "BLAST_ValidateOptions",
     "-threshold 0, -word_size", "deferred: NCBI's texts and exit 1 (S08+, TD-13); LOSAT's parser stops (exit 2)",
     "blastinput/value_parsers.rs", "audit (a) F5"),
    ("src/algo/blast/blastinput/blast_args.cpp", "375-384, 396-406",
     "CFilteringArgs::x_TokenizeFilteringArgs / ExtractAlgorithmOptions", "-seg WINDOW LOCUT HICUT",
     f"ported: single spaces, an int window, {AUDIT}; deferred: NCBI's texts and exit 1 of the errors (S08+)",
     "blastinput/value_parsers.rs parse_seg_filtering", "tests/cli_v2.rs; audit (a) F2"),
    ("src/algo/blast/core/blast_filter.c", "1147-1154", "BlastSetUp_Filter", "SEG parameters not above 0",
     f"ported: {AUDIT}", "blastinput/value_parsers.rs SegSpec::params", "fixtures options.seg_*; audit (a) F3"),
    ("src/algo/blast/api/local_blast.cpp", "54-103", "SplitQuery_GetChunkSize", "CHUNK_SIZE",
     f"ported: the error for a size not divisible by 3, the rejection of a value NCBI cannot convert, {AUDIT}",
     "tblastx/blast_engine/run_impl.rs check_query_split_environment", "fixtures env.chunk*"),
    ("src/algo/blast/api/split_query_aux_priv.cpp", "53-60, 100-111",
     "SplitQuery_GetOverlapChunkSize / SplitQuery_CalculateNumChunks", "OVERLAP_CHUNK_SIZE",
     f"ported: no effect on an ungapped search, the rejection of a value NCBI cannot convert, {AUDIT}",
     "tblastx/blast_engine/run_impl.rs check_query_split_environment", "fixture env.overlap_negative_fmt6"),
    ("src/app/blast/blast_app_util.cpp", "206-210", "InitializeSubject", "BL2SEQ_LEGACY",
     f"rejected: the legacy bl2seq report, {AUDIT}", "tblastx/blast_engine/run_impl.rs check_unsupported_environment",
     "audit (a) F4"),
    ("src/algo/blast/format/blast_format.cpp", "1540, 1550-1551", "CBlastFormat::PrintOneResultSet",
     "the titles of the shown subjects",
     f"ported: the outfmt 0 title checks of the subjects with hits only, after the search, {AUDIT}",
     "tblastx/blast_engine/run_impl.rs check_shown_subject_titles; tblastn/args.rs run_local",
     "fixtures title.*; punct_defline.py; audit (c) F1"),
    ("src/corelib/ncbistr.cpp", "4523-4590", "NStr::HtmlDecode", "character references in a title",
     f"rejected: exactly where NCBI decodes (its scan and entity table), {AUDIT}",
     "report/defline.rs ncbi_nucleotide_title_is_decoded", "html_titles.py; fixture title.kept_fmt0"),
    ("src/objtools/readers/fasta_reader_utils.cpp", "215-225", "CFastaReader::ParseDefLine", "a DEL byte",
     f"ported: the title ends at the first byte below a space only, {AUDIT}",
     "tblastx/blast_engine/run_impl.rs check_report_titles", "fixtures input.del_deflines_*"),
    ("src/app/blast/blast_app_util.hpp", "252-255", "CATCH_ALL: std::ios::failure", "standard output closed",
     "exception: approved exception 6 of PD-LOSAT-CLI-NONSEARCH-DIFFERENCES (maintainer, 2026-10-03, DW-17); "
     "Rust's runtime opens /dev/null on a closed standard output before main, so LOSAT writes the report there "
     "and succeeds (the check of 7fbbfad96 failed callers with a /dev/null opened read and write; removed in "
     "1117e8c17)",
     "cli.rs report_standard_output", "tests/run_local_tblastx.rs a_standard_output_on_dev_null_is_written"),
    ("src/app/blast/tblastx_app.cpp", "132-137", "CTblastxApp::Run", "Query is Empty!, then BATCH_SIZE",
     f"ported: BATCH_SIZE before the queries and the deferred subject checks; a subject bio cannot read "
     f"rejected where the search starts, {AUDIT}", "tblastx/blast_engine/run_impl.rs run, search_cli",
     "fixture input.subject_blank_first_empty_query; tests a_batch_size_that_is_not_an_integer_is_rejected"),
    ("src/algo/blast/core/*.c", "(audit (a) COVERAGE.tsv section S7)", "114 functions of the search core",
     "default and accepted options",
     "faithful: cited by LOSAT's references; no difference in the audit's differential runs, except "
     "BLAST_Cutoffs (the next row)",
     "LOSAT/src/algorithm/tblastx", "audit (a): about 985 cases, round 2 1600 random runs; investigations/"),
    ("src/algo/blast/core/blast_stat.c", "4090-4135", "BLAST_Cutoffs", "a large -evalue against a small search space",
     "ported: at least the caller's 1, aa4b6f8c3", "tblastx/ncbi_cutoffs.rs blast_cutoffs_from_one",
     "fixture cutoff.floor_evalue_1e10; audit (b) F-1, (a) round 2 N3"),
    ("src/algo/blast/core/blast_stat.c", "4040-4063", "BlastKarlinEtoS_simple", "-",
     "faithful: the score from E", "tblastx/ncbi_cutoffs.rs cutoff_score_from_evalue", "audit (a) round 2"),
    ("src/algo/blast/core/blast_stat.c", "4157-", "BLAST_KarlinStoE_simple", "-", "faithful: the e-value of a score",
     "stats/", "audit (a) round 2: 1600 random runs equal"),
    ("src/algo/blast/core/blast_stat.c", "4418-, 4491-", "BLAST_SmallGapSumE / BLAST_UnevenGapSumE", "-",
     "faithful", "stats/sum_statistics.rs", "audit (a) round 2"),
    ("src/algo/blast/core/blast_stat.c", "5041-", "BLAST_ComputeLengthAdjustment", "-", "faithful",
     "tblastx/ncbi_cutoffs.rs compute_eff_lengths_tblastx", "equal footers in every outfmt 0 run"),
    ("src/algo/blast/core/blast_engine.c", "221-262", "s_GetNextSubjectChunk", "subjects over 5,000,000 nt",
     "faithful", "tblastx/blast_engine/run_impl.rs (subject chunks)",
     "audit (a) round 2: 28 runs with subjects of 5.53 and 11.03 Mnt"),
    ("src/algo/blast/core/lookup_wrap.c", "99-100", "LookupTableWrapInit", "the bone type of the lookup table",
     "faithful", "tblastx/lookup/backbone.rs", "audit (a) round 2: 28 runs"),
    ("src/app/blast/blast_app_util.cpp", "728-750", "s_PreFetchSeqs", "PRE_FETCH_SEQS_LIMIT",
     "rejected: a value that NCBI cannot convert, aa4b6f8c3",
     "tblastx/blast_engine/run_impl.rs check_unsupported_environment",
     "tests/run_local_tblastx.rs; audit (b) F-4b, (a) round 2 N2"),
    ("src/objtools/readers/fasta.cpp", "375-384", "CFastaReader::ReadOneSeq", "lines before the first defline",
     "ported: a subject bio cannot read is deferred only after lines NCBI skips; other text rejected at "
     "read, 5fa53b3f0; the records after those lines are checked and warned about at read, and a file of "
     "such lines only has no subject (exit 3), bc521f450", "blastn/input.rs "
     "only_skipped_lines_before_first_defline, from_first_defline; tblastx run_impl.rs run",
     "audit (a) round 2 N1, round 3 N1-ii and N1-iii, (b) N-b1, (d) L1; fixtures "
     "input.subject_comment_only_{empty_query,fmt0}, input.subject_comment_first_title_empty_query"),
    ("src/algo/blast/api/seqsrc_multiseq.cpp", "-", "the subject sequence source, SplitQuery_SetEffectiveSearchSpace",
     "-", "faithful: rows E43, E44", "LOSAT/src/algorithm/tblastx", "audit (a): equal footers in every outfmt 0 run"),
]


def read(path: Path) -> list[dict[str, str]]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))


def main() -> int:
    rows = []
    for name in RANGES:
        inventory = read(HERE / "inventory" / f"{name}.tsv")
        results = {row["row"]: row for row in read(HERE / "inventory" / f"result_{name}.tsv")}
        numbered = [(str(index), row) for index, row in enumerate(inventory, 1)]
        numbered += [(key, None) for key in results if not key.isdigit()]
        for key, before in numbered:
            result = results[key]
            final = RESOLUTION.get((name, key), result["s08_result"])
            if result["s08_result"] == "GAP" and (name, key) not in RESOLUTION:
                raise SystemExit(f"GAP row {name} {key} has no resolution")
            rows.append({
                "range": name, "row": key,
                "ncbi_file": result["ncbi_file"], "ncbi_lines": result["ncbi_lines"],
                "ncbi_function": result["ncbi_function"], "branch": result["branch"],
                "status_before": before["status"] if before else "-",
                "s08_result": result["s08_result"], "s08_final": final,
                "losat_location": result["losat_location"], "evidence": result["evidence"],
                "notes": result["notes"],
            })
    rows = [{**row, "s08_final": AUDIT_FINAL.get((row["range"], row["row"]), row["s08_final"])} for row in rows]
    for number, (ncbi_file, lines, function, branch, final, location, evidence) in enumerate(AUDIT_ROWS, 1):
        rows.append({"range": "H", "row": str(number), "ncbi_file": ncbi_file, "ncbi_lines": lines,
                     "ncbi_function": function, "branch": branch, "status_before": "-", "s08_result": "-",
                     "s08_final": final, "losat_location": location, "evidence": evidence,
                     "notes": "added by the independent audit"})
    with open(HERE / "INVENTORY.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n",
                                quoting=csv.QUOTE_NONE, escapechar="\\")
        writer.writeheader()
        writer.writerows(rows)
    before = collections.Counter(row["status_before"] for row in rows)
    result = collections.Counter(row["s08_result"] for row in rows)
    final = collections.Counter(row["s08_final"].split(":")[0] for row in rows)
    print(f"rows {len(rows)}")
    print("status_before", dict(before))
    print("s08_result", dict(result))
    print("s08_final", dict(final))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
