#!/usr/bin/env python3
"""Focused tests for the NCBI-only platform characterization driver.

NCBI references:
- ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:38-45
- ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1098-1108
These tests preserve the default twelve-field order and emitted row order.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tempfile
import unittest

import yaml


TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parents[1]
SCRIPT_PATH = TESTS_DIR / "characterize_ncbi_platform_variance_v010.py"
WORKFLOW_PATH = REPO_ROOT / ".github/workflows/ncbi-only-platform-characterization.yml"


def load_driver():
    spec = importlib.util.spec_from_file_location(
        "ncbi_characterization_tested", SCRIPT_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load characterization driver")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


driver = load_driver()


class InventoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.authority, cls.catalog = driver.load_platform_authority(REPO_ROOT)
        cls.plan = driver.build_command_plan(
            REPO_ROOT,
            cls.authority,
            cls.catalog,
            {
                program: Path(f"${{NCBI_BIN}}/{program}")
                for program in ("blastn", "blastp", "tblastx")
            },
            Path("${EVIDENCE_ROOT}"),
        )

    def test_exact_six_representatives_and_order(self) -> None:
        self.assertEqual(
            [(pair["program"], pair["case_id"]) for pair in self.plan],
            [
                ("blastn", "PesePMNV.MjPMNV.task_blastn"),
                ("blastn", "Sakai.MG1655.megablast"),
                ("blastp", "pairwise_default_serial"),
                ("tblastx", "p03_mela_pemojnva"),
                ("tblastx", "p14_ap027131_ap027133_query4"),
                ("tblastx", "d06_ap027131_ap027133_db4"),
            ],
        )
        self.assertEqual(len(self.plan), driver.EXPECTED_REPRESENTATIVE_COUNT)
        self.assertEqual(
            len(self.plan) * 2, driver.EXPECTED_PLATFORM_EXECUTIONS
        )
        self.assertEqual(
            len(self.plan) * 2 * 3, driver.EXPECTED_TOTAL_EXECUTIONS
        )
        self.assertEqual(driver.EXPECTED_TOTAL_AUTHORITATIVE_EXECUTIONS, 18)
        self.assertEqual(driver.EXPECTED_TOTAL_DIAGNOSTIC_EXECUTIONS, 18)
        self.assertEqual(driver.EXPECTED_TOTAL_EXECUTIONS, 36)

    def test_commands_come_from_existing_authorities(self) -> None:
        by_case = {step["case_id"]: step for step in self.plan}
        pese = by_case["PesePMNV.MjPMNV.task_blastn"]
        self.assertEqual(
            pese["authoritative"]["semantics"],
            {
                "task": "blastn",
                "outfmt": "7",
                "threads": "1",
                "query_gencode": "1",
                "db_gencode": "1",
            },
        )
        sakai = by_case["Sakai.MG1655.megablast"]
        self.assertEqual(sakai["authoritative"]["semantics"]["task"], "megablast")
        self.assertEqual(
            by_case["p14_ap027131_ap027133_query4"]["authoritative"]
            ["semantics"]["query_gencode"],
            "4",
        )
        self.assertEqual(
            by_case["d06_ap027131_ap027133_db4"]["authoritative"]
            ["semantics"]["db_gencode"],
            "4",
        )
        for pair in self.plan:
            for role in ("authoritative", "diagnostic"):
                search = pair[role]
                self.assertNotIn("LOSAT", search["argv"][0])
                self.assertIn(pair["program"], search["argv"][0])
                self.assertEqual(search["semantics"]["threads"], "1")
            self.assertEqual(
                pair["diagnostic"]["semantics"]["outfmt"],
                driver.DIAGNOSTIC_OUTFMT,
            )

    def test_exact_ten_unique_fixture_paths(self) -> None:
        observed = {
            details["repository_relative"]
            for pair in self.plan
            for details in pair["inputs"].values()
        }
        self.assertEqual(
            observed,
            {
                "LOSAT/tests/fasta/AP027131.fasta",
                "LOSAT/tests/fasta/AP027133.fasta",
                "LOSAT/tests/fasta/AP027152.fasta",
                "LOSAT/tests/fasta/AP027202.fasta",
                "LOSAT/tests/fasta/MG1655.fna",
                "LOSAT/tests/fasta/MelaMJNV.fasta",
                "LOSAT/tests/fasta/PajaWSV.faa",
                "LOSAT/tests/fasta/PemoMJNVA.fasta",
                "LOSAT/tests/fasta/Sakai.fna",
                "LOSAT/tests/fasta/SicyWSV.faa",
            },
        )

    def test_each_pair_differs_only_in_outfmt_and_reuses_input_paths(self) -> None:
        for pair in self.plan:
            comparison = driver.compare_pair_commands(
                pair["authoritative"]["argv"], pair["diagnostic"]["argv"]
            )
            self.assertTrue(comparison["identical_non_output_arguments"])
            self.assertEqual(
                comparison["diagnostic_outfmt"], driver.DIAGNOSTIC_OUTFMT
            )
            self.assertEqual(len(comparison["argument_differences"]), 1)
            self.assertEqual(
                driver.input_arguments(
                    pair["authoritative"]["argv"], pair["program"]
                ),
                driver.input_arguments(
                    pair["diagnostic"]["argv"], pair["program"]
                ),
            )

    def test_pair_identity_records_identical_input_hashes(self) -> None:
        for pair in self.plan:
            for details in pair["inputs"].values():
                details.update(
                    {
                        "sha256": "a" * 64,
                        "byte_length": 12,
                        "newline_profile": {
                            "lf_count": 2,
                            "bare_lf_count": 2,
                            "crlf_count": 0,
                            "bare_cr_count": 0,
                        },
                        "lf_only": True,
                    }
                )
            identity = driver._pair_identity_record(pair, "b" * 64)
            self.assertTrue(identity["identical_non_output_arguments"])
            self.assertTrue(identity["identical_input_sha256"])
            self.assertFalse(identity["diagnostic_output_is_authoritative"])


class FixtureByteTests(unittest.TestCase):
    def test_git_blob_read_ignores_changed_checkout_newlines(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            subprocess.run(["git", "init", "-q"], cwd=root, check=True)
            subprocess.run(
                ["git", "config", "user.email", "test@example.invalid"],
                cwd=root,
                check=True,
            )
            subprocess.run(
                ["git", "config", "user.name", "Fixture Test"], cwd=root, check=True
            )
            fixture = root / "fixture.fa"
            fixture.write_bytes(b">q\nACGT\n")
            subprocess.run(["git", "add", "fixture.fa"], cwd=root, check=True)
            subprocess.run(
                ["git", "commit", "-q", "-m", "fixture"], cwd=root, check=True
            )
            source_sha = (
                subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root)
                .decode("ascii")
                .strip()
            )
            fixture.write_bytes(b">q\r\nACGT\r\n")
            _, observed = driver._git_blob(
                root, PurePosixPath("fixture.fa"), source_sha=source_sha
            )
            self.assertEqual(observed, b">q\nACGT\n")
            self.assertNotEqual(observed, fixture.read_bytes())

    def test_lf_crlf_and_bare_cr_profiles_are_separate(self) -> None:
        data = b"a\r\nb\nc\rd"
        self.assertEqual(
            driver.newline_profile(data),
            {
                "lf_count": 2,
                "bare_lf_count": 1,
                "crlf_count": 1,
                "bare_cr_count": 1,
            },
        )
        self.assertEqual(driver.lf_diagnostic_bytes(data), b"a\nb\nc\nd")
        self.assertEqual(
            driver.require_lf_fixture(b">q\nACGT\n", "fixture.fa")["crlf_count"],
            0,
        )
        with self.assertRaises(driver.CharacterizationFailure):
            driver.require_lf_fixture(b">q\r\nACGT\r\n", "fixture.fa")

    def test_three_manifest_gate_detects_a_platform_difference(self) -> None:
        common = {
            "repository_relative": "LOSAT/tests/fasta/q.fa",
            "git_blob_oid": "a" * 40,
            "expected_git_blob_raw_sha256": "b" * 64,
            "staged_relative": "fixtures/LOSAT/tests/fasta/q.fa",
            "staged_sha256": "b" * 64,
            "byte_length": 4,
            "newline_profile": {
                "lf_count": 1,
                "bare_lf_count": 1,
                "crlf_count": 0,
                "bare_cr_count": 0,
            },
            "lf_only": True,
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = []
            for platform_id in ("windows-x64", "macos-arm64", "macos-x64"):
                path = root / f"{platform_id}.json"
                path.write_text(
                    json.dumps(
                        {
                            "platform_id": platform_id,
                            "candidate_sha": "d" * 40,
                            "fixture_source_sha": driver.PR6_CERT_SHA_V3,
                            "inputs": [dict(common)],
                        }
                    ),
                    encoding="utf-8",
                )
                paths.append(path)
            output = root / "gate.json"
            driver.verify_fixture_invariance(paths, output)
            self.assertEqual(
                json.loads(output.read_text(encoding="utf-8"))["status"],
                "CROSS_PLATFORM_FIXTURES_BYTE_IDENTICAL",
            )
            changed = json.loads(paths[-1].read_text(encoding="utf-8"))
            changed["inputs"][0]["staged_sha256"] = "c" * 64
            paths[-1].write_text(json.dumps(changed), encoding="utf-8")
            with self.assertRaises(driver.CharacterizationFailure):
                driver.verify_fixture_invariance(paths, output)


class StructuredComparisonTests(unittest.TestCase):
    def test_reports_exact_row_field_values_and_key_order(self) -> None:
        linux = b"q\ts\t73.413\t3020\t101\t4\t1\t3020\t1\t3020\t0.0\t900\n"
        platform = b"q\ts\t73.444\t3016\t105\t4\t1\t3020\t1\t3020\t0.0\t900\r\n"
        linux_rows = driver.parse_structured_rows(linux)
        platform_rows = driver.parse_structured_rows(platform)
        comparison = driver.structured_comparison(
            platform_rows, linux_rows, platform, linux
        )
        self.assertTrue(comparison["row_key_order_equal"])
        self.assertEqual(comparison["differing_row_count"], 1)
        self.assertEqual(
            comparison["differing_rows"][0]["field_differences"],
            [
                {"field": "pident", "linux": "73.413", "platform": "73.444"},
                {"field": "length", "linux": "3020", "platform": "3016"},
                {"field": "mismatch", "linux": "101", "platform": "105"},
            ],
        )
        self.assertNotEqual(
            comparison["platform_raw_sha256"], comparison["linux_raw_sha256"]
        )

    def test_sakai_and_d06_contexts_remain_distinct(self) -> None:
        representatives = {item.case_id: item for item in driver.REPRESENTATIVES}
        sakai = driver._case_authority_context(
            representatives["Sakai.MG1655.megablast"]
        )
        d06 = driver._case_authority_context(
            representatives["d06_ap027131_ap027133_db4"]
        )
        self.assertEqual(sakai["existing_authority"], "SOURCE_UNDETERMINED_ACCEPTED")
        self.assertEqual(d06["existing_authority"], "approved_db_gencode_deviation")
        self.assertNotEqual(sakai["existing_authority"], d06["existing_authority"])


class DiagnosticSeparationTests(unittest.TestCase):
    def test_pair_capture_preserves_raw_bytes_and_authoritative_runs_first(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            authority, catalog = driver.load_platform_authority(REPO_ROOT)
            pair = driver.build_command_plan(
                REPO_ROOT,
                authority,
                catalog,
                {
                    program: Path(f"/official/{program}")
                    for program in ("blastn", "blastp", "tblastx")
                },
                root,
            )[0]
            calls = []
            original_run_bytes = driver.run_bytes

            def fake_run_bytes(command, cwd):
                outfmt = driver._option_value(command, "-outfmt")
                calls.append(outfmt)
                output = Path(driver._option_value(command, "-out"))
                payload = (
                    b"authoritative\r\n"
                    if outfmt != driver.DIAGNOSTIC_OUTFMT
                    else b"diagnostic\r\n"
                )
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_bytes(payload)
                return subprocess.CompletedProcess(command, 0, b"", b"")

            driver.run_bytes = fake_run_bytes
            try:
                authoritative, _ = driver._execute_paired_search(
                    REPO_ROOT, pair, "authoritative"
                )
                diagnostic, _ = driver._execute_paired_search(
                    REPO_ROOT, pair, "diagnostic"
                )
            finally:
                driver.run_bytes = original_run_bytes

            self.assertEqual(calls, ["7", driver.DIAGNOSTIC_OUTFMT])
            self.assertEqual(authoritative, b"authoritative\r\n")
            self.assertEqual(diagnostic, b"diagnostic\r\n")
            case_dir = root / "blastn" / "PesePMNV.MjPMNV.task_blastn"
            self.assertEqual(
                (case_dir / "authoritative" / "ncbi.raw.out").read_bytes(),
                b"authoritative\r\n",
            )
            self.assertFalse((case_dir / "ncbi.execution.raw.out").exists())

    def test_rich_diagnostic_field_order_and_standard_projection(self) -> None:
        values = [
            "q",
            "s",
            "73.413",
            "3020",
            "101",
            "4",
            "1",
            "3020",
            "1",
            "3020",
            "0.0",
            "900",
            "2210",
            "10AC5",
            "AC-GT",
            "ACTGT",
        ]
        diagnostic_rows = driver.parse_diagnostic_rows(
            ("\t".join(values) + "\n").encode("utf-8")
        )
        authoritative_rows = driver.parse_structured_rows(
            ("\t".join(values[:12]) + "\n").encode("utf-8")
        )
        self.assertEqual(driver.standard_projection(diagnostic_rows), authoritative_rows)
        self.assertEqual(diagnostic_rows[0]["score"], "2210")
        self.assertEqual(diagnostic_rows[0]["btop"], "10AC5")
        self.assertEqual(diagnostic_rows[0]["qseq"], "AC-GT")
        self.assertEqual(diagnostic_rows[0]["sseq"], "ACTGT")

    def test_diagnostic_path_cannot_satisfy_authoritative_loader(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            authoritative = root / "authoritative" / "ncbi.raw.out"
            diagnostic = root / "diagnostic" / "ncbi.raw.out"
            authoritative.parent.mkdir()
            diagnostic.parent.mkdir()
            authoritative.write_bytes(b"authoritative\n")
            diagnostic.write_bytes(b"diagnostic\n")
            self.assertEqual(
                driver.load_platform_authoritative_raw(authoritative),
                b"authoritative\n",
            )
            with self.assertRaises(driver.CharacterizationFailure):
                driver.load_platform_authoritative_raw(diagnostic)

    def test_changed_hsp_record_keeps_authority_and_diagnostic_separate(self) -> None:
        linux = b"q\ts\t90\t10\t1\t0\t1\t10\t1\t10\t1e-4\t20\n"
        platform = b"q\ts\t80\t10\t2\t0\t1\t10\t1\t10\t1e-4\t20\n"
        diagnostic_values = platform.decode("utf-8").rstrip("\n").split("\t")
        diagnostic_values.extend(["40", "5AC4", "ACGT", "ATGT"])
        diagnostic_rows = driver.parse_diagnostic_rows(
            ("\t".join(diagnostic_values) + "\n").encode("utf-8")
        )
        comparison = driver.structured_comparison(
            driver.parse_structured_rows(platform),
            driver.parse_structured_rows(linux),
            platform,
            linux,
        )
        records = driver.changed_hsp_diagnostic_records(
            comparison, diagnostic_rows
        )
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["platform_diagnostic"]["score"], "40")
        self.assertFalse(records[0]["diagnostic_used_for_authoritative_pass_fail"])

    def test_alternative_representation_requires_distinct_rich_evidence(self) -> None:
        expected = {"windows-x64", "macos-arm64", "macos-x64"}
        common = {"btop": "10", "qseq": "ACGT", "sseq": "ACGT"}
        rows = [
            {"platform_id": platform_id, **common}
            for platform_id in sorted(expected)
        ]
        finding, distinct_count, matched = driver.rich_representation_finding(
            rows, expected
        )
        self.assertNotIn("ALTERNATIVE_EQUAL_SCORING", finding)
        self.assertEqual(distinct_count, 1)
        self.assertEqual(matched, expected)

        rows[-1] = {
            "platform_id": "windows-x64",
            "btop": "4AC5",
            "qseq": "AC-GT",
            "sseq": "ACTGT",
        }
        finding, distinct_count, _ = driver.rich_representation_finding(
            rows, expected
        )
        self.assertIn("ALTERNATIVE_EQUAL_SCORING", finding)
        self.assertEqual(distinct_count, 2)


class WorkflowContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.authority, _ = driver.load_platform_authority(REPO_ROOT)

    def test_yaml_parses_and_has_exact_topology(self) -> None:
        document = yaml.load(
            WORKFLOW_PATH.read_text(encoding="utf-8"), Loader=yaml.BaseLoader
        )
        jobs = document["jobs"]
        for job_name in ("stage-fixtures", "characterize"):
            matrix = jobs[job_name]["strategy"]["matrix"]["include"]
            self.assertEqual(
                [(item["platform_id"], item["runner"]) for item in matrix],
                [
                    ("windows-x64", "windows-2025"),
                    ("macos-arm64", "macos-15"),
                    ("macos-x64", "macos-15-intel"),
                ],
            )

    def test_workflow_has_exact_thirty_six_searches_and_no_losat_run(self) -> None:
        text = WORKFLOW_PATH.read_text(encoding="utf-8")
        self.assertEqual(
            text.count(
                "Run exactly six authoritative and six diagnostic NCBI searches"
            ),
            1,
        )
        self.assertEqual(driver.EXPECTED_PLATFORM_AUTHORITATIVE_EXECUTIONS, 6)
        self.assertEqual(driver.EXPECTED_PLATFORM_DIAGNOSTIC_EXECUTIONS, 6)
        self.assertEqual(driver.EXPECTED_PLATFORM_EXECUTIONS, 12)
        self.assertEqual(driver.EXPECTED_TOTAL_AUTHORITATIVE_EXECUTIONS, 18)
        self.assertEqual(driver.EXPECTED_TOTAL_DIAGNOSTIC_EXECUTIONS, 18)
        self.assertEqual(driver.EXPECTED_TOTAL_EXECUTIONS, 36)
        self.assertNotIn("cargo build", text)
        self.assertNotIn("target/release/LOSAT", text)
        self.assertNotIn("--native-bin", text)
        self.assertNotIn("certify --", text)
        self.assertNotIn("-version", text)

    def test_workflow_aggregates_three_platforms_without_linux_rerun(self) -> None:
        document = yaml.load(
            WORKFLOW_PATH.read_text(encoding="utf-8"), Loader=yaml.BaseLoader
        )
        analysis = document["jobs"]["analyze-rich-diagnostics"]
        self.assertEqual(
            analysis["needs"], ["characterize", "retained-linux-evidence"]
        )
        text = WORKFLOW_PATH.read_text(encoding="utf-8")
        self.assertIn("analyze-platforms", text)
        self.assertEqual(text.count("--platform-evidence"), 3)
        self.assertIn("18 authoritative plus 18 diagnostic", text)
        self.assertNotIn("linux-rerun", text.lower())

    def test_official_archive_inventory_matches_accepted_platform_authority(
        self,
    ) -> None:
        text = WORKFLOW_PATH.read_text(encoding="utf-8")
        self.assertEqual(
            {
                platform_id: (spec.archive_name, spec.archive_md5)
                for platform_id, spec in self.authority.PLATFORMS.items()
            },
            {
                "windows-x64": (
                    "ncbi-blast-2.17.0+-x64-win64.tar.gz",
                    "dcd973097407a2910061ff4fb51b09fb",
                ),
                "macos-arm64": (
                    "ncbi-blast-2.17.0+-aarch64-macosx.tar.gz",
                    "8dc685eb284713db76de41a4dabf96ad",
                ),
                "macos-x64": (
                    "ncbi-blast-2.17.0+-x64-macosx.tar.gz",
                    "cde7979c0edca21da0567ab172b68b0e",
                ),
            },
        )
        for spec in self.authority.PLATFORMS.values():
            self.assertIn(spec.archive_name, text)
            self.assertIn(spec.runner_label, text)
        self.assertNotIn("brew install", text)
        self.assertNotIn("choco install", text)


if __name__ == "__main__":
    unittest.main()
