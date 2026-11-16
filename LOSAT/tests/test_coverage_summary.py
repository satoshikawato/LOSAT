#!/usr/bin/env python3
"""Tests for the merge of the nightly coverage shards (`coverage_summary.py`).

No NCBI counterpart: CI tooling over cargo-llvm-cov's lcov export.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT_PATH = Path(__file__).with_name("coverage_summary.py")
SPEC = importlib.util.spec_from_file_location("coverage_summary", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
cs = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = cs
SPEC.loader.exec_module(cs)


def record(path: str, functions: list[tuple[int, str, int]], lines: dict[int, int]) -> str:
    """One lcov record as llvm-cov writes it: FN lines, FNDA lines, DA lines."""
    out = [f"SF:{path}"]
    out += [f"FN:{line},{name}" for line, name, _ in functions]
    out += [f"FNDA:{count},{name}" for _, name, count in functions]
    out += [f"FNF:{len(functions)}", "FNH:0"]
    out += [f"DA:{line},{count}" for line, count in lines.items()]
    out += [f"LF:{len(lines)}", "LH:0", "end_of_record"]
    return "\n".join(out) + "\n"


class MergeTests(unittest.TestCase):
    def test_a_line_is_covered_when_one_shard_ran_it(self) -> None:
        files: dict[str, cs.FileCoverage] = {}
        cs.read_lcov(record("/w/src/a.rs", [], {1: 0, 2: 3, 3: 0}), files)
        cs.read_lcov(record("/w/src/a.rs", [], {1: 5, 2: 1, 3: 0}), files)
        self.assertEqual(files["/w/src/a.rs"].lines, {1: 5, 2: 4, 3: 0})
        self.assertEqual(files["/w/src/a.rs"].line_counts(), (2, 3))

    def test_a_line_only_one_shard_lists_is_counted(self) -> None:
        files: dict[str, cs.FileCoverage] = {}
        cs.read_lcov(record("/w/src/a.rs", [], {1: 0}), files)
        cs.read_lcov(record("/w/src/a.rs", [], {2: 0}) + record("/w/src/b.rs", [], {7: 2}), files)
        self.assertEqual(files["/w/src/a.rs"].line_counts(), (0, 2))
        self.assertEqual(files["/w/src/b.rs"].line_counts(), (1, 1))

    def test_records_that_start_on_one_line_are_one_function(self) -> None:
        # Two instantiations of one generic function (line 10), a closure (line 11) and a copy
        # of the function in a second test binary (the same name listed again).
        files: dict[str, cs.FileCoverage] = {}
        text = record("/w/src/a.rs", [(10, "_RINvf_u8", 0), (10, "_RINvf_u16", 0),
                                      (11, "_RNCf0", 0), (10, "_RINvf_u8", 0)], {10: 0, 11: 0})
        cs.read_lcov(text, files)
        self.assertEqual(files["/w/src/a.rs"].function_counts(), (0, 2))
        cs.read_lcov(record("/w/src/a.rs", [(10, "_RINvf_u16", 4)], {10: 4}), files)
        self.assertEqual(files["/w/src/a.rs"].function_counts(), (1, 2))

    def test_the_merged_lcov_reads_back_to_the_same_counts(self) -> None:
        files: dict[str, cs.FileCoverage] = {}
        cs.read_lcov(record("/w/src/a.rs", [(1, "f", 2), (1, "g", 0), (5, "h", 0)],
                            {1: 2, 2: 2, 5: 0}), files)
        cs.read_lcov(record("/w/src/b.rs", [(3, "k", 1)], {3: 1}), files)
        text = cs.write_lcov(files)
        self.assertIn("SF:/w/src/a.rs\n", text)
        self.assertIn("FNF:2\nFNH:1\n", text)
        self.assertIn("LF:3\nLH:2\nend_of_record\n", text)
        again: dict[str, cs.FileCoverage] = {}
        cs.read_lcov(text, again)
        self.assertEqual(again, files)

    def test_a_count_without_its_fn_line_is_skipped(self) -> None:
        files: dict[str, cs.FileCoverage] = {}
        text = "SF:/w/src/a.rs\nFNDA:3,orphan\nFN:4,f\nFNDA:2,f\nDA:4,2\nend_of_record\n"
        cs.read_lcov(text, files)
        self.assertEqual(files["/w/src/a.rs"].functions, {"f": [4, 2]})
        # A name that an earlier record listed is still not counted without its FN line here.
        cs.read_lcov("SF:/w/src/a.rs\nFNDA:5,f\nend_of_record\n", files)
        self.assertEqual(files["/w/src/a.rs"].functions, {"f": [4, 2]})

    def test_lines_outside_a_record_are_skipped(self) -> None:
        files: dict[str, cs.FileCoverage] = {}
        cs.read_lcov("TN:\nDA:1,1\n" + record("/w/src/a.rs", [], {1: 1}) + "DA:9,9\n", files)
        self.assertEqual(files["/w/src/a.rs"].lines, {1: 1})


class SummaryTests(unittest.TestCase):
    def test_modules_are_programs_under_algorithm_and_directories_elsewhere(self) -> None:
        src = Path("/w/LOSAT/src")
        self.assertEqual(cs.module_of("/w/LOSAT/src/algorithm/blastn/lookup/x.rs", src),
                         "algorithm/blastn")
        self.assertEqual(cs.module_of("/w/LOSAT/src/algorithm/mod.rs", src), "algorithm")
        self.assertEqual(cs.module_of("/w/LOSAT/src/seqdb/isam.rs", src), "seqdb")
        self.assertEqual(cs.module_of("/w/LOSAT/src/lib.rs", src), "src/*.rs")
        self.assertEqual(cs.module_of("/w/LOSAT/build.rs", src), "(outside src)")
        # A directory whose name only starts like `src` is not under it.
        self.assertEqual(cs.module_of("/w/LOSAT/src2/a.rs", src), "(outside src)")

    def test_most_uncovered_lines_come_first(self) -> None:
        files: dict[str, cs.FileCoverage] = {}
        cs.read_lcov(record("/w/src/seqdb/a.rs", [(1, "f", 1)], {1: 1, 2: 0}), files)
        cs.read_lcov(record("/w/src/core/b.rs", [(1, "g", 0)], {1: 0, 2: 0, 3: 0}), files)
        text = cs.summary(files, Path("/w/src"), "0123456789abcdef", ["coverage-shard-1: success"])
        self.assertIn("Commit 012345678; shards: coverage-shard-1: success.", text)
        self.assertIn("| lines | 1 | 5 | 20.0 |", text)
        self.assertIn("| functions | 1 | 2 | 50.0 |", text)
        self.assertLess(text.index("| `core` | 1 | 0.0 | 3 | 0.0 |"),
                        text.index("| `seqdb` | 1 | 50.0 | 1 | 100.0 |"))
        self.assertNotIn("did not pass", text)
        self.assertIn("inline `#[cfg(test)]` modules of `src` files are counted", text)

    def test_a_shard_that_failed_or_is_missing_is_named(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            ok = root / "coverage-shard-1"
            ok.mkdir()
            (ok / "lcov.info").write_text(record(str(root / "src/a.rs"), [(1, "f", 1)], {1: 1}))
            (ok / "outcome").write_text("success\n")
            failed = root / "coverage-shard-2"
            failed.mkdir()
            (failed / "lcov.info").write_text(record(str(root / "src/a.rs"), [], {2: 0}))
            (failed / "outcome").write_text("failure\n")
            empty = root / "coverage-shard-3"
            empty.mkdir()
            (empty / "outcome").write_text("failure\n")
            out = root / "out"
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                rc = cs.main(["--src", str(root / "src"), "--shards", "4", "--out", str(out),
                              str(ok), str(failed), str(empty)])
            self.assertEqual(rc, 0)
            text = stdout.getvalue()
            self.assertIn("shards: coverage-shard-1: success, coverage-shard-2: failure, "
                          "coverage-shard-3: no report, coverage-shard-4: no report.", text)
            self.assertIn("did not pass its tests or gave no report", text)
            self.assertIn("| lines | 1 | 2 | 50.0 |", text)
            self.assertEqual((out / "files.tsv").read_text().splitlines(),
                             ["file\tlines_covered\tlines\tfunctions_covered\tfunctions",
                              "a.rs\t1\t2\t1\t1"])
            self.assertIn("LF:2\nLH:1\n", (out / "lcov.info").read_text())

    def test_no_report_at_all_fails_and_names_every_shard(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stdout = io.StringIO()
            # An unmatched glob reaches the script as the pattern itself.
            with contextlib.redirect_stdout(stdout):
                rc = cs.main(["--src", str(root), "--shards", "2", "--out", str(root / "out"),
                              str(root / "coverage-shard-*")])
            self.assertEqual(rc, 1)
            self.assertIn("shards: coverage-shard-1: no report, coverage-shard-2: no report.",
                          stdout.getvalue())

    def test_shards_are_named_in_number_order(self) -> None:
        names = ["coverage-shard-10", "coverage-shard-2", "coverage-shard-1"]
        ordered = sorted(((n, "success") for n in names), key=cs.shard_order)
        self.assertEqual([n for n, _ in ordered],
                         ["coverage-shard-1", "coverage-shard-2", "coverage-shard-10"])


if __name__ == "__main__":
    unittest.main()
