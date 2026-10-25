#!/usr/bin/env python3
"""Tests for ci_fast_regressions.py and frozen_allowlist.py (no LOSAT binary is run)."""
from __future__ import annotations

import contextlib
import io
import json
import re
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import blastn_regression_fixtures as fixtures  # noqa: E402
import ci_fast_regressions as fast  # noqa: E402
import db_search_fixtures as db_fixtures  # noqa: E402
import tblastx_regression_fixtures as tblastx_fixtures  # noqa: E402
import xinclude_fixtures as xinclude  # noqa: E402
from frozen_allowlist import Allowlist  # noqa: E402

sys.path.insert(0, str(fast.FIXTURE_DIR))
import run_oracle  # noqa: E402  (the outfmt 0 manifest's reader and writer)

FROZEN, KNOWN, OTHER = "f" * 64, "a" * 64, "0" * 64


def allowlist() -> Allowlist:
    document = {"schema": "losat-frozen-mismatch-allowlist-v1", "entries": [
        {"program": "blastn", "case_id": "known", "frozen_sha256": FROZEN, "allowed_actual_sha256": KNOWN}]}
    handle = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
    json.dump(document, handle)
    handle.close()
    return Allowlist(Path(handle.name))


def row(case_id: str, output: str, expected: str = FROZEN, program: str = "blastn") -> dict[str, str]:
    return {"program": program, "case_id": case_id, "command": "[]", "exit": "0", "output_sha256": output,
            "stderr_sha256": "e", "expected_sha256": expected, "expected_source": "gate_a"}


class SamClRuleTests(unittest.TestCase):
    """The one rule of the SAM (outfmt 17) rows: the CL tag of the @PG line is not compared."""

    WITH = b"@HD\tVN:1.2\n@SQ\tSN:Query_1\tLN:3\n@PG\tID:0\tVN:2.17.0+\tCL:/x/blastn -a \tPN:blastn\nr\t0\n"
    WITHOUT = b"@HD\tVN:1.2\n@SQ\tSN:Query_1\tLN:3\n@PG\tID:0\tVN:2.17.0+\tPN:blastn\nr\t0\n"

    def test_a_row_that_differs_in_the_command_line_only_is_equal(self):
        other = self.WITH.replace(b"/x/blastn -a ", b"losat blastn -b -num_threads 4 ")
        self.assertEqual(db_fixtures.strip_sam_cl(self.WITH), self.WITHOUT)
        self.assertEqual(db_fixtures.strip_sam_cl(other), self.WITHOUT)
        self.assertEqual(db_fixtures.strip_sam_cl(self.WITHOUT), self.WITHOUT)

    def test_any_other_byte_is_compared(self):
        for old, new in ((b"r\t0", b"r\t1"), (b"LN:3", b"LN:4"), (b"PN:blastn", b"PN:blastx"),
                         (b"VN:2.17.0+", b"VN:2.17.1+"), (b"@HD", b"@XX")):
            changed = self.WITH.replace(old, new)
            self.assertNotEqual(db_fixtures.strip_sam_cl(changed), self.WITHOUT, (old, new))

    def test_the_rule_applies_to_outfmt_17_only(self):
        for outfmt in ("17", "17 SQ", " 17 SR ", "17 SQ SR"):
            self.assertTrue(db_fixtures.is_sam_outfmt(outfmt), outfmt)
        for outfmt in ("0", "6", "5", "16", "7 qseqid", "170", ""):
            self.assertFalse(db_fixtures.is_sam_outfmt(outfmt), outfmt)


class AllowlistTests(unittest.TestCase):
    def test_classify(self):
        listed = allowlist()
        self.assertEqual(listed.classify("blastn", "known", FROZEN, KNOWN), "allowed")
        self.assertEqual(listed.classify("blastn", "known", FROZEN, FROZEN), "stale")
        self.assertEqual(listed.classify("blastn", "known", FROZEN, OTHER), "mismatch")
        self.assertEqual(listed.classify("blastn", "known", FROZEN, KNOWN, executed=False), "mismatch")
        self.assertEqual(listed.classify("blastn", "known", OTHER, OTHER), "match")
        self.assertEqual(listed.classify("blastn", "other", FROZEN, FROZEN), "match")
        self.assertEqual(listed.classify("blastn", "other", FROZEN, KNOWN), "mismatch")
        self.assertEqual(listed.unexecuted({("blastn", "other")}, {"blastn"}), [("blastn", "known")])
        self.assertEqual(listed.unexecuted(set(), {"blastp"}), [])

    def test_repository_allowlist_is_empty_and_gate_a_reads_the_v011_manifest(self):
        self.assertEqual(Allowlist().entries, {})
        # PD-NCBI-PLATFORM-VARIANCE 1.3: the approved Sakai output is the Gate A
        # expectation itself, read from the certification driver's manifest.
        expected = fast.capture_outputs.gate_a()
        self.assertEqual(len(expected), 43)
        self.assertEqual(
            expected[("blastn", "Sakai.MG1655.megablast")],
            "ac8177764f35b2fb92777f8b1d01ea646cc55ae21b05f1f486ee231348a34a45",
        )
        # The S02 baseline's outputs are the v0.1.1 expectations in all 43 cases.
        baseline = fast.read_tsv(fast.BASELINE)
        gate_a_rows = {key: row for key, row in baseline.items() if row["expected_source"] == "gate_a"}
        self.assertEqual(set(gate_a_rows), set(expected))
        self.assertEqual({key: row["output_sha256"] for key, row in gate_a_rows.items()}, expected)


class CheckRowsTests(unittest.TestCase):
    def test_baseline_and_frozen_checks(self):
        baseline = {("blastn", "known"): row("known", KNOWN), ("blastn", "same"): row("same", FROZEN)}
        failures, allowed = fast.check_rows([row("known", KNOWN), row("same", FROZEN)], baseline, allowlist(), {"blastn"}, True)
        self.assertEqual((failures, allowed), ([], ["blastn/known"]))
        failures, _ = fast.check_rows([row("known", FROZEN), row("same", FROZEN)], baseline, allowlist(), {"blastn"}, False)
        self.assertTrue(any("output_sha256 differs" in line for line in failures))
        self.assertTrue(any("remove it" in line for line in failures))
        failures, _ = fast.check_rows([row("same", FROZEN)], baseline, allowlist(), {"blastn"}, True)
        self.assertTrue(any("not executed" in line for line in failures))
        failures, _ = fast.check_rows([row("same", FROZEN)], baseline, allowlist(), {"blastn"}, False)
        self.assertEqual(failures, [])


class BlastnFixtureTests(unittest.TestCase):
    def test_manifest_matches_cases_and_files(self):
        for module in (fixtures, tblastx_fixtures):
            rows = module.read_manifest()
            self.assertEqual([(r["case_id"], r["argv"], r["losat_extra"], r.get("env") or "", r.get("oracle_env") or "")
                              for r in rows], module.CASES, module.__name__)
            for r in rows:
                data = (module.FIXTURES / f"{r['case_id']}.out").read_bytes()
                self.assertEqual(module.sha256(data), r["stdout_sha256"], r["case_id"])
                err = module.FIXTURES / f"{r['case_id']}.err"
                self.assertEqual(module.sha256(err.read_bytes()) if err.exists() else "", r["stderr_sha256"], r["case_id"])


class DbSearchFixtureTests(unittest.TestCase):
    def test_manifest_rows_match_frozen_files(self):
        rows = db_fixtures.read_manifest()[1]
        self.assertEqual({r["program"] for r in rows}, set(db_fixtures.PROGRAMS))
        # ci_fast runs the -db rows of every program the manifest has, at 1 and 4 threads.
        self.assertEqual(set(fast.DB_SEARCH_PROGRAMS), set(db_fixtures.PROGRAMS))
        self.assertTrue(set(fast.DB_SEARCH_PROGRAMS) <= set(fast.PROGRAMS))
        self.assertEqual(len({r["row_id"] for r in rows}), len(rows))
        for r in rows:
            if r["kind"] != "ncbi":
                continue
            data = (db_fixtures.FROZEN / f"{r['row_id']}.out").read_bytes()
            self.assertEqual(db_fixtures.sha256(data), r["stdout_sha256"], r["row_id"])
            self.assertEqual(str(len(data)), r["stdout_bytes"], r["row_id"])
            err = db_fixtures.FROZEN / f"{r['row_id']}.err"
            self.assertEqual(db_fixtures.sha256(err.read_bytes()) if err.exists() else "", r["stderr_sha256"], r["row_id"])


    def test_env_tokens_are_environment_not_arguments(self):
        row = {"extra_args": "BLASTDB=tests/seqdb/taxdb/good BATCH_SIZE=100 -evalue 1e-5 -word_size 3",
               "db": "d", "query": "q", "outfmt": "16"}
        self.assertEqual(db_fixtures.search_argv(row),
                         ["-db", "d", "-query", "q", "-evalue", "1e-5", "-word_size", "3", "-outfmt", "16"])
        env = db_fixtures.row_env(row, {"PATH": "/bin"})
        self.assertEqual((env["BLASTDB"], env["BATCH_SIZE"], env["PATH"]),
                         ("tests/seqdb/taxdb/good", "100", "/bin"))
        self.assertNotIn("CHUNK_SIZE", env)
        # a token of another name is an argument (and is not silently dropped)
        self.assertFalse(db_fixtures.is_env_token("-evalue"))
        self.assertFalse(db_fixtures.is_env_token("NCBI=1"))

    def test_clean_env_drops_what_changes_a_result_and_seals_home(self):
        dirty = {"PATH": "/bin", "HOME": "/home/someone", "USERPROFILE": "C:/Users/someone",
                 "NCBI_CONFIG__BLAST__BLASTDB": "/x", "NCBI_CONFIG_PATH": "/y", "NCBI": "/z",
                 "BLASTDB": "/d", "BLAST_USE_NCBI_ENV": "yes", "LOSAT_X_SOMETHING": "1",
                 "BATCH_SIZE": "1", "CHUNK_SIZE": "2", "OLD_FSC": "1", "BL2SEQ_LEGACY": "1", "LANG": "C"}
        clean = db_fixtures.clean_env(dirty, "/empty/home")
        self.assertEqual(clean, {"PATH": "/bin", "HOME": "/empty/home", "USERPROFILE": "/empty/home",
                                 "LANG": "C"})
        # nothing of the caller's environment is changed
        self.assertEqual(dirty["HOME"], "/home/someone")


def fake_search(directory: Path, master: bytes, sides: list[bytes], stdout: bytes = b"", extra: str = "") -> list[str]:
    """A stand-in for a program that writes the files of outfmt 13/14: the command to run it."""
    script = directory / "fake.py"
    lines = [
        "import sys",
        "out = sys.argv[sys.argv.index('-out') + 1]",
        f"open(out, 'wb').write({master!r})",
        f"for n, side in enumerate({sides!r}, 1):",
        "    open(f'{out}_{n}.xml', 'wb').write(side)",
    ]
    if extra:
        lines.append(f"open(out + {extra!r}, 'wb').write(b'x')")
    lines.append(f"sys.stdout.buffer.write({stdout!r})")
    script.write_text("\n".join(lines) + "\n")
    return [sys.executable, str(script), "blastp", "-outfmt", "14"]


class XincludeFixtureTests(unittest.TestCase):
    """The `files` columns of the manifests (outfmt 13 and 14: a master file and side files)."""

    def manifest_lines(self, path: Path) -> tuple[list[str], list[list[str]]]:
        lines = [line for line in path.read_text().split("\n") if line and not line.startswith("#")]
        return lines[0].split("\t"), [line.split("\t") for line in lines[1:]]

    def test_rows_without_the_columns_keep_their_text(self):
        # the outfmt 0 manifest: a freeze writes every row as it was read (no trailing empty cells added)
        header, rows = self.manifest_lines(run_oracle.MANIFEST)
        self.assertEqual(header[run_oracle.BASE_COLUMNS:], xinclude.FILES_COLUMNS)
        for cells in rows:
            row = dict(zip(header, cells))
            self.assertEqual(run_oracle.row_cells(header, row), cells, cells[0])
        # the -db manifest likewise, and its header is the old one only while no row has the columns
        header, rows = self.manifest_lines(db_fixtures.MANIFEST)
        parsed = db_fixtures.read_manifest()[1]
        self.assertEqual(header, db_fixtures.header_cells(parsed))
        for cells, row in zip(rows, parsed):
            self.assertEqual(db_fixtures.row_cells(row), cells, cells[0])
        old = {column: "" for column in db_fixtures.COLUMNS}
        self.assertEqual(db_fixtures.header_cells([old]), db_fixtures.COLUMNS)
        self.assertEqual(len(db_fixtures.row_cells(old)), len(db_fixtures.COLUMNS))
        self.assertEqual(db_fixtures.header_cells([{**old, "files": "xinclude"}]),
                         db_fixtures.COLUMNS + xinclude.FILES_COLUMNS)

    def test_files_rows_match_the_frozen_master_and_side_files(self):
        found = 0
        for manifest_rows, frozen, key in (
                (run_oracle.read_manifest()[2], run_oracle.ENGINE / "tests/fixtures/outfmt0", "fixture_id"),
                (db_fixtures.read_manifest()[1], db_fixtures.FROZEN, "row_id")):
            for row in manifest_rows:
                if not xinclude.wants_files(row):
                    self.assertEqual(row.get("side_sha256", ""), "", row[key])
                    continue
                found += 1
                row_id, outfmt = row[key], row["outfmt"]
                master = (frozen / f"{row_id}.out").read_bytes()
                self.assertEqual(xinclude.sha256(master), row["stdout_sha256"], row_id)
                self.assertEqual(str(len(master)), row["stdout_bytes"], row_id)
                hashes = row["side_sha256"].split(",") if row["side_sha256"] else []
                names = [xinclude.side_name(row_id, n, outfmt) for n in range(1, len(hashes) + 1)]
                self.assertEqual(sorted(path.name for path in frozen.glob(f"{row_id}.out_*")), sorted(names), row_id)
                self.assertEqual([xinclude.sha256((frozen / name).read_bytes()) for name in names], hashes, row_id)
                # the master lists the side files by those names (13: no end of line after the last)
                listed = re.findall(rb'href="([^"]*)"' if outfmt == "14" else rb'"File": "([^"]*)"', master)
                self.assertEqual([name.decode() for name in listed], names, row_id)
                self.assertEqual(master.endswith(b"\n"), outfmt == "14", row_id)
        self.assertGreaterEqual(found, 9)

    def test_side_files_are_collected_in_the_order_of_their_numbers(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            for n in (10, 2, 1):
                (directory / f"r.out_{n}.xml").write_bytes(str(n).encode())
            (directory / "r.out").write_bytes(b"master")
            master, sides, stray = xinclude.collect(directory, "r", "14")
            self.assertEqual((master, stray), (b"master", ["side file numbers [1, 2, 10] are not 1..3"]))
            self.assertEqual(sides, [b"1", b"2", b"10"])
            (directory / "r.out_10.xml").unlink()
            (directory / "other.txt").write_bytes(b"")
            (directory / "r.out_3.json").write_bytes(b"")
            _, sides, stray = xinclude.collect(directory, "r", "14")
            self.assertEqual((sides, stray), ([b"1", b"2"], ["other.txt", "r.out_3.json"]))

    def test_a_run_is_compared_with_the_frozen_master_side_files_and_standard_output(self):
        def first_difference(expected: bytes, actual: bytes) -> str:
            return f"expected {expected!r}, got {actual!r}"

        with tempfile.TemporaryDirectory() as tmp:
            directory = Path(tmp)
            frozen = directory / "frozen"
            frozen.mkdir()
            row = {"row_id": "r", "outfmt": "14", "stdout_sha256": xinclude.sha256(b"M"), "stdout_bytes": "1",
                   "side_sha256": ",".join(xinclude.sha256(side) for side in (b"A", b"B"))}
            written = xinclude.run_with_files(fake_search(directory, b"M", [b"A", b"B"]), "r", "14")
            self.assertEqual(xinclude.observed(written)["side_sha256"], row["side_sha256"])
            xinclude.freeze(written, frozen, "r", "14")
            self.assertEqual(sorted(path.name for path in frozen.iterdir()), ["r.out", "r.out_1.xml", "r.out_2.xml"])
            self.assertEqual(xinclude.problems(written, row, frozen, first_difference), [])
            # a different master, side file, number of side files, a standard output, a stray file
            for master, sides, stdout, extra, expect in (
                    (b"N", [b"A", b"B"], b"", "", "master"),
                    (b"M", [b"A", b"C"], b"", "", "side file 2"),
                    (b"M", [b"A"], b"", "", "1 side files, expected 2"),
                    (b"M", [b"A", b"B"], b"x", "", "stdout is not empty"),
                    (b"M", [b"A", b"B"], b"", "_stray", "unexpected files")):
                other = xinclude.run_with_files(fake_search(directory, master, sides, stdout, extra), "r", "14")
                found_problems = xinclude.problems(other, row, frozen, first_difference)
                self.assertTrue(any(expect in line for line in found_problems), (expect, found_problems))
            # a freeze refuses standard output
            noisy = xinclude.run_with_files(fake_search(directory, b"M", [b"A"], b"x"), "r", "14")
            with self.assertRaises(SystemExit):
                xinclude.freeze(noisy, frozen, "r", "14")

    def test_a_comparison_counts_a_stray_file_as_a_difference_and_writes_nothing(self):
        # run_oracle.py in compare mode (no --freeze) with a stand-in for the NCBI program: a stray file
        # or standard output is a difference of that fixture (the comparison goes on), --out stays empty
        # and the manifest is not written.
        def compare(stray: bool, stdout: bool) -> tuple[int, str, list[str]]:
            with tempfile.TemporaryDirectory() as tmp:
                bin_dir, out = Path(tmp) / "bin", Path(tmp) / "out"
                bin_dir.mkdir()
                out.mkdir()
                program = bin_dir / "blastp"
                program.write_text("\n".join([
                    f"#!{sys.executable}", "import sys",
                    "if sys.argv[1] == '-version':", "    print('blastp: 2.17.0+'); sys.exit(0)",
                    "out = sys.argv[sys.argv.index('-out') + 1]",
                    "open(out, 'wb').write(b'M')", "open(out + '_1.xml', 'wb').write(b'A')",
                    "open(out + '_stray', 'wb').write(b'x')" if stray else "",
                    "sys.stdout.write('x')" if stdout else ""]) + "\n")
                program.chmod(0o755)
                good = {"fixture_id": "r", "program": "blastp", "outfmt": "14", "files": "xinclude",
                        "stdout_sha256": xinclude.sha256(b"M"), "stdout_bytes": "1", "stderr_sha256": "",
                        "side_sha256": xinclude.sha256(b"A")}
                rows = [good, {**good, "fixture_id": "s"}]
                saved = (sys.argv, run_oracle.read_manifest, run_oracle.search_argv)
                sys.argv = ["run_oracle.py", "--bin-dir", str(bin_dir), "--out", str(out)]
                run_oracle.read_manifest = lambda: ([], list(good), rows)
                run_oracle.search_argv = lambda row: ["-outfmt", "14"]
                captured = io.StringIO()
                try:
                    with contextlib.redirect_stdout(captured):
                        status = run_oracle.main()
                finally:
                    sys.argv, run_oracle.read_manifest, run_oracle.search_argv = saved
                return status, captured.getvalue(), [path.name for path in out.iterdir()]

        status, printed, written = compare(stray=False, stdout=False)
        self.assertEqual((status, written), (0, []), printed)
        for stray, stdout in ((True, False), (False, True)):
            status, printed, written = compare(stray, stdout)
            self.assertEqual((status, written), (1, []), printed)
            self.assertIn("differing=['r', 's']", printed)

    def test_the_files_column_and_the_format_are_validated(self):
        self.assertFalse(xinclude.wants_files({}))
        self.assertFalse(xinclude.wants_files({"files": ""}))
        self.assertTrue(xinclude.wants_files({"files": "xinclude"}))
        with self.assertRaises(SystemExit):
            xinclude.wants_files({"files": "yes"})
        self.assertEqual((xinclude.extension("13"), xinclude.extension("14")), ("json", "xml"))
        with self.assertRaises(SystemExit):
            xinclude.extension("16")
        self.assertEqual(xinclude.side_name("id", 3, "13"), "id.out_3.json")


class SelectionTests(unittest.TestCase):
    def test_paths_select_programs(self):
        everything = set(fast.PROGRAMS)
        self.assertEqual(fast.select_programs(["docs/README.md", "web/app/src/main.ts"]), set())
        self.assertEqual(fast.select_programs(["LOSAT/src/core/blast_util.rs"]), everything)
        self.assertEqual(fast.select_programs(["LOSAT/Cargo.lock"]), everything)
        self.assertEqual(fast.select_programs([".github/workflows/ci.yml"]), everything)
        # A program's file selects the programs whose code imports its directory, directly or through
        # another program, not every program: blastn imports no other program, so a tblastn file does
        # not select it.
        tblastn = fast.select_programs(["LOSAT/src/algorithm/tblastn/args.rs"])
        self.assertEqual(tblastn, fast.program_dependents()["tblastn"])
        self.assertIn("tblastn", tblastn)
        self.assertNotIn("blastn", tblastn)
        blastn = fast.select_programs(["LOSAT/src/algorithm/blastn/hsp.rs"])
        self.assertIn("blastn", blastn)
        # TBLASTX uses BLASTN's hit list (blastn/hsp.rs HitList) and blank-input check.
        self.assertIn("tblastx", blastn)
        for program in fast.PROGRAMS:
            self.assertIn(program, fast.program_dependents()[program])


class StartupTraceNoticeTests(unittest.TestCase):
    """`check_startup_trace_notice` against stub executables standing in for `losat`."""

    def check(self, body: str) -> list[str]:
        with tempfile.TemporaryDirectory() as tmp:
            program = Path(tmp) / "losat"
            program.write_text(f"#!{sys.executable}\nimport os, sys\n{body}\n")
            program.chmod(0o755)
            return fast.check_startup_trace_notice(program)

    def test_a_default_build_passes(self):
        self.assertEqual(self.check(
            "print('losat 1.0')\n"
            "if os.environ.get('LOSAT_STARTUP_TRACE'):\n"
            f"    sys.stderr.write({fast.STARTUP_TRACE_NOTICE!r})"), [])

    def test_a_build_that_reads_the_variable_fails(self):
        failures = self.check(
            "print('losat 1.0')\n"
            "if os.environ.get('LOSAT_STARTUP_TRACE'):\n"
            "    sys.stderr.write('[startup] enter main\\n')")
        self.assertEqual(len(failures), 1)
        self.assertIn("default build", failures[0])

    def test_stderr_without_the_variable_and_a_failing_exit_fail(self):
        failures = self.check("print('losat 1.0')\nsys.stderr.write('warning\\n')\nsys.exit(1)")
        self.assertTrue(any("exited" in failure for failure in failures))
        self.assertTrue(any("stderr without" in failure for failure in failures))


if __name__ == "__main__":
    unittest.main()
