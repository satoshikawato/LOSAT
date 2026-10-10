#!/usr/bin/env python3
"""Focused tests for the NCBI unit-test case ledger.

NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/blastfilter_unit_test.cpp:2014-2016
```c++
#if SEQLOC_MIX_QUERY_OK
// Test the masking with a Seq-loc of type mix
BOOST_AUTO_TEST_CASE(DustSeqlocMix) {
```

The tests freeze how NCBI's Boost.Test cases are listed (live, fixture and
compiled-out cases) and how the ledger checker ties citations, ledger rows
and LOSAT tests together; they do not test BLAST behavior.
"""

from __future__ import annotations

import collections
import importlib.util
import os
import sys
import tempfile
import unittest
from pathlib import Path


CHECKER_PATH = Path(__file__).with_name("ncbi_unit_case_ledger.py")
SPEC = importlib.util.spec_from_file_location("ncbi_unit_case_ledger", CHECKER_PATH)
assert SPEC is not None and SPEC.loader is not None
ledger = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ledger
SPEC.loader.exec_module(ledger)

REPO_ROOT = Path(__file__).resolve().parents[2]
COMMIT = "598d8ae6a72b923127ba2fbfaffd48e4c83bfbf4"
API = "c++/src/algo/blast/unit_tests/api"

NCBI_FIXTURE = """\
#include <corelib/test_boost.hpp>
BOOST_AUTO_TEST_SUITE(fixture)

BOOST_AUTO_TEST_CASE(LiveCase)
{
    const char* text = "}{ BOOST_AUTO_TEST_CASE(NotACase)";
    BOOST_REQUIRE(text);
}

BOOST_FIXTURE_TEST_CASE(FixtureCase, CFixture) {
    // a closing brace in a comment: }
    BOOST_REQUIRE(true);
}

BOOST_AUTO_TEST_CASE_TIMEOUT(TimeoutCase, 3);
BOOST_AUTO_TEST_CASE(TimeoutCase)
{
}

/* BOOST_AUTO_TEST_CASE(CommentedOut) { } */

#if 0
BOOST_AUTO_TEST_CASE(DeadCase) {
}
#endif

#if SEQLOC_MIX_QUERY_OK
BOOST_AUTO_TEST_CASE(SeqlocMixCase) {
#if 0
    int a = 1;
#else
    int a = 2;
#endif
}
#endif

#if 1
BOOST_AUTO_TEST_CASE(IfOneCase) {
}
#else
BOOST_AUTO_TEST_CASE(IfOneElseCase) {
}
#endif

#ifdef NCBI_THREADS
BOOST_AUTO_TEST_CASE(ThreadsCase) {
}
#endif

#define DECLARE_TEST(name, wordsize)                    \\
BOOST_AUTO_TEST_CASE( name##ScanOffsetSize##wordsize ) { \\
    Run(wordsize);                                      \\
}

DECLARE_TEST(Tiny, 4);
DECLARE_TEST(Disco_Coding_16_,
             11)
#if 0
DECLARE_TEST(Dead, 8);
#endif

BOOST_AUTO_TEST_SUITE_END()
"""


def write_repo(root: Path, rust: dict[str, str], ledger_rows: list[str], cases: list[tuple] | None = None) -> None:
    if cases is None:
        cases = [
            ("blasthits", f"{API}/blasthits_unit_test.cpp", "testHSPListSort", 1037, 1102, "AUTO"),
            ("blasthits", f"{API}/blasthits_unit_test.cpp", "testSubjectBestHit", 1300, 1340, "AUTO"),
            ("aalookup", f"{API}/aalookup_unit_test.cpp", "testDebruijnPSSM", 310, 433, "DISABLED"),
        ]
    ledger.write_cases(root / ledger.DEFAULT_CASES, COMMIT, [ledger.Case(*row) for row in cases])
    (root / ledger.DEFAULT_LEDGER).write_text(
        "\n".join(["\t".join(ledger.LEDGER_HEADER), *ledger_rows]) + "\n", encoding="utf-8"
    )
    for rel, text in rust.items():
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")


CITED_TEST = """\
pub fn sort_hsps() {}

#[cfg(test)]
mod tests {
    // NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/blasthits_unit_test.cpp:1094-1099 testHSPListSort
    // ```c++
    // Blast_HSPListSortByScore(m_HspList);
    // ```
    #[test]
    fn hsp_list_sort_matches_ncbi() {}

    #[test]
    fn subject_best_hit() {}
}
"""


def row(module: str, case: str, status: str, test: str = "", note: str = "", since: str = "2026-10-10") -> str:
    return "\t".join((module, case, status, test, note, since))


class ExtractTests(unittest.TestCase):
    def test_live_fixture_and_compiled_out_cases(self) -> None:
        cases = ledger.extract_cases_from_text("m", "f.cpp", NCBI_FIXTURE)
        kinds = {case.case: case.kind for case in cases}
        self.assertEqual(
            kinds,
            {
                "LiveCase": "AUTO",
                "FixtureCase": "FIXTURE",
                "TimeoutCase": "AUTO",
                "DeadCase": "DISABLED",
                "SeqlocMixCase": "DISABLED",
                "IfOneCase": "AUTO",
                "IfOneElseCase": "DISABLED",
                "ThreadsCase": "AUTO",
                "DECLARE_TEST": "MACRO",
            },
        )
        # The timeout line is not a case; the macro is one row as written, from
        # the #define to its last invocation, with the runtime count in the note.
        self.assertEqual(len(cases), 9)
        spans = {case.case: (case.line_start, case.line_end) for case in cases}
        notes = {case.case: case.note for case in cases}
        self.assertEqual(spans["DECLARE_TEST"], (50, 59))
        self.assertEqual(
            notes["DECLARE_TEST"],
            "2 cases at run time from DECLARE_TEST(...) at lines 55-59, named <name>ScanOffsetSize<wordsize>; "
            "1 more compiled out",
        )
        self.assertEqual(notes["LiveCase"], "")
        self.assertEqual(spans["LiveCase"], (4, 8))
        self.assertEqual(spans["FixtureCase"], (10, 13))
        self.assertEqual(spans["SeqlocMixCase"], (28, 34))

    def test_extract_reads_every_source_and_rejects_duplicate_names(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            src = Path(directory) / "c++" / "src"
            for _, folder, pattern in ledger.SOURCES:
                name = pattern.replace("*", "fake_unit_test")
                path = src / folder / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(f"BOOST_AUTO_TEST_CASE(Case_{path.stem}) {{\n}}\n", encoding="utf-8")
            cases = ledger.extract_cases(Path(directory))
            self.assertEqual(len(cases), len(ledger.SOURCES))
            self.assertIn(f"{API}/blasthits_unit_test.cpp", {case.file for case in cases})
            extra = src / "algo/blast/unit_tests/blast_format/second_unit_test.cpp"
            extra.write_text("BOOST_AUTO_TEST_CASE(Case_fake_unit_test) {\n}\n", encoding="utf-8")
            with self.assertRaises(ledger.LedgerError):
                ledger.extract_cases(Path(directory))

    def test_committed_case_list_counts(self) -> None:
        cases = ledger.load_cases(REPO_ROOT / ledger.DEFAULT_CASES)
        self.assertEqual(cases.commit, COMMIT)
        total = collections.Counter(case.module for case in cases.cases)
        disabled = collections.Counter(case.module for case in cases.cases if case.kind == "DISABLED")
        fixture = collections.Counter(case.module for case in cases.cases if case.kind == "FIXTURE")
        expected = {
            "optionshandle": 75, "blastfilter": 67, "bl2seq": 98, "blastinput": 96,
            "blasthits": 25, "blastsetup": 59, "split_query": 33,
        }
        for module, count in expected.items():
            self.assertEqual(total[module], count, module)
        self.assertEqual(total["blast_format"] + total["align_format"], 38)
        self.assertEqual(disabled["blastfilter"], 1)
        self.assertEqual(disabled["bl2seq"], 3)
        self.assertEqual(fixture["optionshandle"], 61)
        self.assertEqual(len(cases.cases), 679)
        self.assertEqual(total["ntscan"], 2)  # DiscontigTwoSubjects + the DECLARE_TEST macro
        macro = cases.by_module_case[("ntscan", "DECLARE_TEST")]
        self.assertEqual((macro.line_start, macro.line_end, macro.kind), (846, 905, "MACRO"))
        self.assertTrue(macro.note.startswith("44 cases at run time"), macro.note)
        self.assertEqual([case.case for case in cases.cases if case.kind == "MACRO"], ["DECLARE_TEST"])
        self.assertNotIn(("ntscan", "TinyScanOffsetSize4"), cases.by_module_case)

    @unittest.skipUnless(os.environ.get("NCBI_SRC"), "NCBI_SRC is not set")
    def test_committed_case_list_matches_a_fresh_extract(self) -> None:
        ncbi_src = Path(os.environ["NCBI_SRC"])
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / "NCBI_CASES.tsv"
            ledger.write_cases(out, ledger.ncbi_commit(ncbi_src), ledger.extract_cases(ncbi_src))
            self.assertEqual(
                out.read_text(encoding="utf-8"),
                (REPO_ROOT / ledger.DEFAULT_CASES).read_text(encoding="utf-8"),
            )


class CheckTests(unittest.TestCase):
    def run_check(self, rust: dict[str, str], rows: list[str]) -> list[str]:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_repo(root, rust, rows)
            failures, _ = ledger.check(
                root, ledger.load_cases(root / ledger.DEFAULT_CASES), ledger.load_ledger(root / ledger.DEFAULT_LEDGER)
            )
        return [f"{f.module}\t{f.case}\t{f.reason}" for f in failures]

    def test_consistent_repository_passes(self) -> None:
        failures = self.run_check(
            {"LOSAT/src/hits.rs": CITED_TEST},
            [
                row("blasthits", "testHSPListSort", "ported", "LOSAT/src/hits.rs::hsp_list_sort_matches_ncbi"),
                row("blasthits", "testSubjectBestHit", "partial", "LOSAT/src/hits.rs::subject_best_hit"),
                row("aalookup", "testDebruijnPSSM", "n-a", note="NCBI compiles it out"),
            ],
        )
        self.assertEqual(failures, [])

    def test_citation_must_name_a_case_within_its_lines_above_a_test(self) -> None:
        unknown = CITED_TEST.replace("testHSPListSort", "testNoSuchCase")
        failures = self.run_check({"LOSAT/src/hits.rs": unknown}, [])
        self.assertTrue(any("cited case not in NCBI_CASES.tsv" in f for f in failures), failures)

        outside = CITED_TEST.replace(":1094-1099 ", ":1094-1200 ")
        failures = self.run_check({"LOSAT/src/hits.rs": outside}, [row("blasthits", "testHSPListSort", "ported")])
        self.assertTrue(any("outside the case (1037-1102)" in f for f in failures), failures)

        production = CITED_TEST.replace("    #[test]\n    fn hsp_list_sort", "    fn hsp_list_sort", 1)
        failures = self.run_check({"LOSAT/src/hits.rs": production}, [row("blasthits", "testHSPListSort", "ported")])
        self.assertTrue(any("not a #[test]" in f for f in failures), failures)

        malformed = CITED_TEST.replace("testHSPListSort\n", "testHSPListSort trailing words\n")
        failures = self.run_check({"LOSAT/src/hits.rs": malformed}, [])
        self.assertTrue(any("malformed citation" in f for f in failures), failures)

        other_commit = CITED_TEST.replace("(598d8ae6)", "(1234567a)")
        failures = self.run_check({"LOSAT/src/hits.rs": other_commit}, [row("blasthits", "testHSPListSort", "ported")])
        self.assertTrue(any("citation commit 1234567a is not" in f for f in failures), failures)

        detached = CITED_TEST.replace("    // ```\n    #[test]\n    fn hsp_list_sort", "    // ```\n\n    #[test]\n    fn hsp_list_sort")
        failures = self.run_check({"LOSAT/src/hits.rs": detached}, [row("blasthits", "testHSPListSort", "ported")])
        self.assertTrue(any("not directly above a fn" in f for f in failures), failures)

        ignored = CITED_TEST.replace("    #[test]\n    fn hsp_list_sort", "    #[test]\n    #[ignore]\n    fn hsp_list_sort")
        failures = self.run_check({"LOSAT/src/hits.rs": ignored}, [row("blasthits", "testHSPListSort", "ported")])
        self.assertTrue(any("which is #[ignore]d" in f for f in failures), failures)

    def test_macro_case_is_cited_by_its_body_and_invocation_lines(self) -> None:
        cases = [
            ("ntscan", f"{API}/ntscan_unit_test.cpp", "DECLARE_TEST", 846, 905, "MACRO",
             "44 cases at run time from DECLARE_TEST(...) at lines 855-905, named <name>ScanOffsetSize<wordsize>"),
        ]
        source = """\
#[cfg(test)]
mod tests {
    // NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/ntscan_unit_test.cpp:847-853,855 DECLARE_TEST
    #[test]
    fn tiny_scan_offset_size_4() {}
}
"""
        for lines, expected in (("847-853,855", []), ("847-853,906", ["outside the case (846-905)"])):
            with tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                write_repo(
                    root,
                    {"LOSAT/src/scan.rs": source.replace("847-853,855", lines)},
                    [row("ntscan", "DECLARE_TEST", "partial", "LOSAT/src/scan.rs::tiny_scan_offset_size_4")],
                    cases,
                )
                failures, _ = ledger.check(
                    root, ledger.load_cases(root / ledger.DEFAULT_CASES), ledger.load_ledger(root / ledger.DEFAULT_LEDGER)
                )
            reasons = [f.reason for f in failures]
            self.assertEqual(len(reasons), len(expected), reasons)
            for reason, text in zip(reasons, expected):
                self.assertIn(text, reason)

    def test_multi_line_attributes_are_stepped_over(self) -> None:
        source = CITED_TEST.replace(
            "    #[test]\n    fn hsp_list_sort",
            "    #[cfg(all(\n        target_arch = \"x86_64\",\n        feature = \"simd\"\n    ))]\n"
            "    #[test]\n    #[cfg_attr(\n        miri,\n        allow(dead_code)\n    )]\n    fn hsp_list_sort",
        )
        failures = self.run_check(
            {"LOSAT/src/hits.rs": source},
            [row("blasthits", "testHSPListSort", "ported", "LOSAT/src/hits.rs::hsp_list_sort_matches_ncbi")],
        )
        self.assertEqual(failures, [])

    def test_cited_case_needs_a_ported_or_partial_row(self) -> None:
        failures = self.run_check({"LOSAT/src/hits.rs": CITED_TEST}, [])
        self.assertIn("blasthits\ttestHSPListSort\tcited in code but has no ledger row (LOSAT/src/hits.rs:5)", failures)
        failures = self.run_check({"LOSAT/src/hits.rs": CITED_TEST}, [row("blasthits", "testHSPListSort", "to-port")])
        self.assertTrue(any("cited in code but ledger status is to-port" in f for f in failures), failures)

    def test_ledger_rows_are_well_formed(self) -> None:
        failures = self.run_check(
            {},
            [
                row("blasthits", "testHSPListSort", "done"),
                row("blasthits", "testSubjectBestHit", "superseded"),
                row("blasthits", "testNoSuchCase", "n-a", note="x"),
                row("aalookup", "testDebruijnPSSM", "n-a", since="today"),
                row("aalookup", "testDebruijnPSSM", "n-a"),
            ],
        )
        joined = "\n".join(failures)
        self.assertIn("status 'done' is not one of", joined)
        self.assertIn("superseded row needs the reason in note", joined)
        self.assertIn("names a case not in NCBI_CASES.tsv", joined)
        self.assertIn("since 'today' is not YYYY-MM-DD", joined)
        self.assertIn("duplicate ledger row", joined)

    def test_removed_test_fails_until_the_ledger_is_fixed(self) -> None:
        rows = [row("blasthits", "testSubjectBestHit", "partial", "LOSAT/src/hits.rs::subject_best_hit")]
        removed = CITED_TEST.replace("    #[test]\n    fn subject_best_hit() {}\n", "")
        failures = self.run_check({"LOSAT/src/hits.rs": removed}, rows + [row("blasthits", "testHSPListSort", "ported")])
        self.assertIn("blasthits\ttestSubjectBestHit\tlosat_test fn not found: LOSAT/src/hits.rs::subject_best_hit", failures)
        self.assertIn("blasthits\t-\tledger counts 2 linked ported/partial cases, code backs 1", failures)

        moved = self.run_check({"LOSAT/tests/hits.rs": CITED_TEST}, rows)
        self.assertTrue(any("losat_test file not found: LOSAT/src/hits.rs" in f for f in moved), moved)

    def test_ported_row_must_be_linked_but_partial_may_carry_evidence(self) -> None:
        failures = self.run_check(
            {},
            [
                row("blasthits", "testHSPListSort", "ported"),
                row("blasthits", "testSubjectBestHit", "partial", note="evidence (b28e42e1): LOSAT/src/x.rs:10"),
                row("aalookup", "testDebruijnPSSM", "partial"),
            ],
        )
        self.assertEqual(
            failures,
            [
                "aalookup\ttestDebruijnPSSM\tpartial row has no test, citation or evidence note",
                "blasthits\ttestHSPListSort\tported row has no citation line and no losat_test",
            ],
        )

    def test_fixture_reference_must_exist(self) -> None:
        failures = self.run_check(
            {"LOSAT/tests/fixtures/case.tsv": "row_7\n"},
            [
                row("blasthits", "testHSPListSort", "e2e", "LOSAT/tests/fixtures/case.tsv#row_7"),
                row("blasthits", "testSubjectBestHit", "e2e", "LOSAT/tests/fixtures/missing.tsv;LOSAT/tests/fixtures/case.tsv#row_8;LOSAT/tests"),
            ],
        )
        self.assertEqual(
            failures,
            [
                "blasthits\ttestSubjectBestHit\tlosat_test fixture file not found: LOSAT/tests/fixtures/missing.tsv",
                "blasthits\ttestSubjectBestHit\tlosat_test fixture id not found: LOSAT/tests/fixtures/case.tsv#row_8",
                "blasthits\ttestSubjectBestHit\tlosat_test fixture file not found: LOSAT/tests",
            ],
        )

    def test_ported_and_partial_rows_need_a_test_fn_not_a_file(self) -> None:
        failures = self.run_check(
            {"LOSAT/src/hits.rs": CITED_TEST.replace("// NCBI unit test", "// NCBI reference"), "LOSAT/tests/fixtures/case.tsv": "x\n"},
            [
                row("blasthits", "testHSPListSort", "ported", "LOSAT/src/hits.rs"),
                row("blasthits", "testSubjectBestHit", "partial", "LOSAT/tests/fixtures/case.tsv"),
            ],
        )
        self.assertEqual(
            failures,
            [
                "blasthits\ttestHSPListSort\tlosat_test names a Rust file without ::fn: LOSAT/src/hits.rs",
                "blasthits\ttestHSPListSort\tported row names no #[test] fn and no citation (fixtures alone are e2e)",
                "blasthits\ttestSubjectBestHit\tpartial row names no #[test] fn and no citation (fixtures alone are e2e)",
            ],
        )

    def test_ignored_test_does_not_back_a_row(self) -> None:
        source = CITED_TEST.replace("    #[test]\n    fn subject_best_hit", "    #[test]\n    #[ignore]\n    fn subject_best_hit")
        failures = self.run_check(
            {"LOSAT/src/hits.rs": source},
            [
                row("blasthits", "testHSPListSort", "ported"),
                row("blasthits", "testSubjectBestHit", "partial", "LOSAT/src/hits.rs::subject_best_hit"),
            ],
        )
        self.assertIn("blasthits\ttestSubjectBestHit\tlosat_test fn is #[ignore]d: LOSAT/src/hits.rs::subject_best_hit", failures)


class ReportTests(unittest.TestCase):
    def test_report_counts_statuses_and_unlisted_cases(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_repo(
                root,
                {},
                [
                    row("blasthits", "testHSPListSort", "ported", "LOSAT/src/hits.rs::t"),
                    row("aalookup", "testDebruijnPSSM", "partial", note="evidence"),
                ],
            )
            text = ledger.report(
                ledger.load_cases(root / ledger.DEFAULT_CASES), ledger.load_ledger(root / ledger.DEFAULT_LEDGER)
            )
        lines = text.splitlines()
        self.assertEqual(lines[0], "| NCBI module | cases | ported | partial | e2e | to-port | n-a | superseded | unlisted |")
        self.assertIn("| blasthits | 2 | 1 | 0 | 0 | 0 | 0 | 0 | 1 |", lines)
        self.assertIn("| aalookup | 1 | 0 | 1 | 0 | 0 | 0 | 0 | 0 |", lines)
        self.assertIn("| **total** | **3** | **1** | **1** | **0** | **0** | **0** | **0** | **1** |", lines)
        self.assertTrue(lines[-1].endswith("(no citation, no losat_test): 1."))


if __name__ == "__main__":
    unittest.main()
