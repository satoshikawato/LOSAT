#!/usr/bin/env python3
"""Comparison-only Stage B source, CLI, and genetic-code checks.

NCBI: c++/src/objects/seqfeat/gc.prt:105-357; blastinput/blast_args.cpp:997-1056;
blastinput/tblastn_args.cpp:45-152; api/blast_aux.cpp:588-613.
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import re
import subprocess

COMMIT = "598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4"
IDS = (1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16,
       21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33)
ROOT = Path(__file__).resolve().parents[3]
EVIDENCE = Path(__file__).resolve().parent


def run(cmd: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(cmd, text=True, capture_output=True, check=False)


def check_source(source_repo: Path, out: Path) -> None:
    result = run(["git", "-C", str(source_repo), "show", f"{COMMIT}:c++/src/objects/seqfeat/gc.prt"])
    assert result.returncode == 0, result.stderr
    gc = result.stdout
    source_rows = [(int(code), aa) for code, aa in re.findall(
        r"\bid\s+(\d+)\s*,\s*ncbieaa\s+\"([A-Z*]+)\"", gc)]
    fixture_rows = [(int(parts[0]), parts[1]) for parts in (
        line.split("\t") for line in (EVIDENCE / "gc_prt_27.tsv").read_text().splitlines()[1:])]
    assert source_rows == fixture_rows
    assert tuple(code for code, _ in source_rows) == IDS
    assert all(len(aa) == 64 for _, aa in source_rows)
    api = EVIDENCE.parent / "tlosan_stage_a" / "api_20260923_verified"
    code32 = (api / "code32_fmt6.out").read_text()
    code1 = (api / "subject32_g1_fmt6.out").read_text()
    assert code32.startswith("q1\ts_code32\t100.000\t120\t0\t0\t1\t120\t1\t360\t3.39e-93\t251\n")
    assert code1.startswith("q1\ts_code32\t95.833\t120\t5\t0\t1\t120\t1\t360\t1.75e-85\t231\n")
    (out / "genetic_code_comparison.tsv").write_text(
        "id\tsource_codons\tstop_codons\tsource\n" + "".join(
            f"{code}\t64\t{aa.count('*')}\tgc.prt@{COMMIT}\n" for code, aa in source_rows))
    (out / "source_and_api.txt").write_text(
        f"gc.prt commit: {COMMIT}\n"
        f"gc.prt sha256: {hashlib.sha256(gc.encode()).hexdigest()}\n"
        "source entries: 27; source codons: 1728; TSV entries: 27; TSV codons: 1728\n"
        "code32 Stage A C++ API format-6: identity 100.000, bits 251, E 3.39e-93\n"
        "same subject with code1 C++ API: identity 95.833, bits 231, E 1.75e-85\n")


def check_scoring(source_repo: Path, out: Path) -> None:
    """Compare all Rust scoring CLI pairs with pinned blast_stat.c."""
    source = run(["git", "-C", str(source_repo), "show",
                  f"{COMMIT}:c++/src/algo/blast/core/blast_stat.c"])
    assert source.returncode == 0, source.stderr
    rust = (ROOT / "LOSAT/src/algorithm/tblastn/scoring.rs").read_text()
    rows = []
    for name in ("blosum45", "blosum50", "blosum62", "blosum80", "blosum90",
                 "pam250", "pam30", "pam70", "prot_idenity"):
        label = "IDENTITY" if name == "prot_idenity" else name.upper()
        c_table = re.search(rf"static array_of_8 {name}_values\[[^]]+\] = \{{(.*?)\n\}};",
                            source.stdout, re.S)
        c_prefs = re.search(rf"static Int4 {'prot_identity' if name == 'prot_idenity' else name}_prefs\[[^]]+\] = \{{(.*?)\n\}};",
                            source.stdout, re.S)
        rust_table = re.search(rf"const {label}: &\[\(i32, i32\)\] = &\[(.*?)\];", rust, re.S)
        rust_pref = re.search(rf'"{label}" => \(\((\d+), (\d+)\), {label}\)', rust)
        assert c_table and c_prefs and rust_table and rust_pref, label
        c_pairs = [(int(a), int(b)) for a, b in re.findall(r"\{(\d+),\s*(\d+),", c_table.group(1))]
        rust_pairs = [(int(a), int(b)) for a, b in re.findall(r"\((\d+),\s*(\d+)\)", rust_table.group(1))]
        prefs = re.findall(r"BLAST_MATRIX_(?:BEST|NOMINAL)", c_prefs.group(1))
        assert c_pairs == rust_pairs and len(prefs) == len(c_pairs) + 1, label
        best = c_pairs[prefs.index("BLAST_MATRIX_BEST") - 1]
        assert best == tuple(map(int, rust_pref.groups())), label
        rows.append(f"{label}\t{len(c_pairs)}\t{best[0]}/{best[1]}\tPASS\n")
    (out / "matrix_gap_comparison.tsv").write_text(
        "matrix\tallowed_pairs\tpreferred_gap\tsource_match\n" + "".join(rows))


def check_cli(binary: Path, out: Path) -> None:
    query = ROOT / "docs/evidence/tlosan_stage_a/fixtures/query.faa"
    subject = ROOT / "docs/evidence/tlosan_stage_a/fixtures/subject_code1.fna"
    base = [str(binary), "tblastn", "-query", str(query), "-subject", str(subject)]
    rows = []
    def case(name: str, args: list[str], expected: str) -> None:
        completed = run(base + args)
        diagnostic = (completed.stderr + completed.stdout).strip().splitlines()
        message = diagnostic[0] if diagnostic else ""
        assert completed.returncode != 0, name
        assert expected in completed.stderr, (name, completed.stderr)
        rows.append((name, completed.returncode, message))
    for code in IDS:
        case(f"db_gencode={code}", ["-db_gencode", str(code)], "local search is unimplemented")
    for code in (0, 7, 8, 17, 20, 34, 255, 256, -1, "abc"):
        case(f"invalid_gencode={code}", ["-db_gencode", str(code)], "genetic code")
    case("default", [], "local search is unimplemented")
    case("outfmt=6", ["-outfmt", "6"], "local search is unimplemented")
    case("outfmt=7", ["-outfmt", "7"], "local search is unimplemented")
    case("explicit_bool_values", ["-soft_masking", "true", "-sum_stats", "false"], "local search is unimplemented")
    case("tblastn-fast", ["-task", "tblastn-fast"], "unsupported TBLASTN task")
    case("db_search", ["-db", "db"], "cannot be used with")
    case("psi_checkpoint", ["-in_pssm", "pssm.asn"], "cannot be used with")
    case("remote", ["-remote"], "unsupported TBLASTN remote")
    case("ungapped_composition", ["-ungapped"], "Composition-adjusted searched are not supported")
    case("ungapped_comp_off", ["-ungapped", "-comp_based_stats", "F"], "local search is unimplemented")
    case("compo_suffix", ["-comp_based_stats", "2u"], "local search is unimplemented")
    case("compo_unknown_ncbi", ["-comp_based_stats", "bogus"], "local search is unimplemented")
    case("word_size=7", ["-word_size", "7"], "local search is unimplemented")
    case("word_size=8", ["-word_size", "8"], "Word-size must be less than 8")
    case("threshold=0", ["-threshold", "0"], "Non-zero threshold required")
    case("evalue=0", ["-evalue", "0"], "expect value or cutoff score must be greater than zero")
    case("matrix=PAM30", ["-matrix", "PAM30"], "local search is unimplemented")
    case("matrix=IDENTITY", ["-matrix", "identity"], "local search is unimplemented")
    case("matrix=BAD", ["-matrix", "BAD"], "not a supported matrix")
    case("gap=1/1", ["-gapopen", "1", "-gapextend", "1"], "not supported for BLOSUM62")
    case("matrix=PAM30_gap=9/1", ["-matrix", "PAM30", "-gapopen", "9", "-gapextend", "1"], "local search is unimplemented")
    case("matrix=IDENTITY_word=6", ["-matrix", "IDENTITY", "-word_size", "6"], "Word size larger than 5")
    case("outfmt=5", ["-outfmt", "5"], "unsupported TBLASTN outfmt")
    case("unported_max_hsps", ["-max_hsps", "2"], "unsupported TBLASTN option")
    case("subject_loc", ["-subject_loc", "1-30"], "unsupported TBLASTN -subject_loc")
    case("subject_loc_remote", ["-subject_loc", "1-30", "-remote"], "cannot be used with")
    case("remote_num_threads", ["-remote", "-num_threads", "2"], "cannot be used with")
    direct_psi_remote = run([str(binary), "tblastn", "-subject", str(subject), "-in_pssm", "p.asn", "-remote"])
    assert direct_psi_remote.returncode != 0 and "cannot be used with" in direct_psi_remote.stderr
    rows.append(("remote_in_pssm", direct_psi_remote.returncode, direct_psi_remote.stderr.strip().splitlines()[0]))
    direct_db_range = run([str(binary), "tblastn", "-query", str(query), "-db", "db", "-subject_loc", "1-30"])
    assert direct_db_range.returncode != 0 and "cannot be used with" in direct_db_range.stderr
    rows.append(("db_subject_loc", direct_db_range.returncode, direct_db_range.stderr.strip().splitlines()[0]))
    for name, argv, expected in (
        ("db_only", [str(binary), "tblastn", "-query", str(query), "-db", "db"], "unsupported TBLASTN database"),
        ("psi_only", [str(binary), "tblastn", "-subject", str(subject), "-in_pssm", "p.asn"], "unsupported PSI-TBLASTN"),
    ):
        completed = run(argv)
        assert completed.returncode != 0 and expected in completed.stderr, (name, completed.stderr)
        rows.append((name, completed.returncode, completed.stderr.strip().splitlines()[0]))
    tblastx_code32 = run([str(binary), "tblastx", "-query", str(subject),
                          "-subject", str(subject), "-db_gencode", "32"])
    assert tblastx_code32.returncode != 0 and "unsupported genetic code" in tblastx_code32.stderr
    rows.append(("tblastx_db_gencode=32_rejected", tblastx_code32.returncode,
                 tblastx_code32.stderr.strip().splitlines()[0]))
    (out / "cli_accept_reject.tsv").write_text(
        "case\texit_code\tobserved_diagnostic\n" + "".join(
            f"{name}\t{exit_code}\t{message}\n" for name, exit_code, message in rows))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=ROOT / "LOSAT/target/release/LOSAT")
    parser.add_argument("--source-repo", type=Path,
                        default=Path("/mnt/c/Users/genom/GitHub/ncbi-blast"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    check_source(args.source_repo, args.output)
    check_scoring(args.source_repo, args.output)
    check_cli(args.binary.resolve(), args.output)
    print("PASS: 27 source tables, 1728 codons, 9 matrix gap tables, code32 API control, CLI acceptance/rejection")


if __name__ == "__main__":
    main()
