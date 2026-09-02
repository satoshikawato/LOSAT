#!/usr/bin/env python3
"""Focused tests for the NCBI-only platform characterization driver.

NCBI references:
- ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:38-45
- ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1098-1108
These tests preserve the default twelve-field order and emitted row order.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

import yaml


TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parents[1]
SCRIPT_PATH = TESTS_DIR / "characterize_ncbi_platform_variance_v010.py"
WORKFLOW_PATH = REPO_ROOT / ".github/workflows/ncbi-only-platform-characterization.yml"
REFERENCE_ROOT = TESTS_DIR / "retained_linux_ncbi_reference_v010"
REFERENCE_MANIFEST_PATH = REFERENCE_ROOT / "linux_reference_manifest.json"


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
        self.assertEqual(len(self.plan) * 2, driver.EXPECTED_PLATFORM_EXECUTIONS)
        self.assertEqual(len(self.plan) * 2 * 3, driver.EXPECTED_TOTAL_EXECUTIONS)
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
            by_case["p14_ap027131_ap027133_query4"]["authoritative"]["semantics"][
                "query_gencode"
            ],
            "4",
        )
        self.assertEqual(
            by_case["d06_ap027131_ap027133_db4"]["authoritative"]["semantics"][
                "db_gencode"
            ],
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
            self.assertEqual(comparison["diagnostic_outfmt"], driver.DIAGNOSTIC_OUTFMT)
            self.assertEqual(len(comparison["argument_differences"]), 1)
            self.assertEqual(
                driver.input_arguments(pair["authoritative"]["argv"], pair["program"]),
                driver.input_arguments(pair["diagnostic"]["argv"], pair["program"]),
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


class RetainedLinuxReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.manifest_bytes = REFERENCE_MANIFEST_PATH.read_bytes()
        cls.manifest = json.loads(cls.manifest_bytes.decode("utf-8"))
        cls.payloads = {
            Path(record["repository_relative"])
            .relative_to(driver.RETAINED_LINUX_REFERENCE_ROOT.as_posix())
            .as_posix(): REPO_ROOT.joinpath(
                *PurePosixPath(record["repository_relative"]).parts
            ).read_bytes()
            for record in cls.manifest["files"]
        }
        driver._validate_linux_reference_payloads(cls.manifest, cls.payloads)

    def test_exact_six_cases_and_twelve_bound_payloads(self) -> None:
        self.assertEqual(
            [
                (record["program"], record["case_id"])
                for record in self.manifest["cases"]
            ],
            [
                ("blastn", "PesePMNV.MjPMNV.task_blastn"),
                ("blastn", "Sakai.MG1655.megablast"),
                ("blastp", "pairwise_default_serial"),
                ("tblastx", "p03_mela_pemojnva"),
                ("tblastx", "p14_ap027131_ap027133_query4"),
                ("tblastx", "d06_ap027131_ap027133_db4"),
            ],
        )
        self.assertEqual(self.manifest["case_count"], 6)
        self.assertEqual(self.manifest["tracked_payload_file_count"], 12)
        observed = {
            path.relative_to(REFERENCE_ROOT).as_posix()
            for path in REFERENCE_ROOT.rglob("*")
            if path.is_file() and path != REFERENCE_MANIFEST_PATH
        }
        self.assertEqual(observed, set(self.payloads))

    def test_every_payload_matches_manifest_without_transformation(self) -> None:
        for record in self.manifest["files"]:
            artifact_relative = (
                Path(record["repository_relative"])
                .relative_to(driver.RETAINED_LINUX_REFERENCE_ROOT.as_posix())
                .as_posix()
            )
            data = self.payloads[artifact_relative]
            self.assertEqual(len(data), record["byte_length"])
            self.assertEqual(hashlib.sha256(data).hexdigest(), record["sha256"])
            self.assertEqual(driver.newline_profile(data), record["newline_profile"])

    def test_exact_command_metadata_is_retained_for_each_case(self) -> None:
        for representative in driver.REPRESENTATIVES:
            relative = PurePosixPath(
                representative.program, representative.case_id, "command.json"
            ).as_posix()
            record = json.loads(self.payloads[relative].decode("utf-8"))
            self.assertEqual(
                PurePosixPath(record["argv"][0]).name, representative.program
            )
            self.assertEqual(driver._option_value(record["argv"], "-num_threads"), "1")
            driver._validate_retained_command(representative, self.payloads[relative])

    def test_provenance_binds_pr5_runtime_manifest_and_direct_oracle_record(
        self,
    ) -> None:
        provenance = self.manifest["provenance"]
        self.assertEqual(
            provenance["retained_evidence_manifest"]["sha256"],
            driver.PR5_EVIDENCE_MANIFEST_SHA256,
        )
        self.assertEqual(
            provenance["runtime_evidence_identity"]["runtime_cert_sha"],
            "5845d22ed9842449628a647f29b8c6762511ca59",
        )
        oracle = provenance["oracle_identity"]
        self.assertEqual(
            oracle["source_retained_evidence_path"], "oracle_identities.json"
        )
        self.assertIn("directly", oracle["provenance_statement"])
        self.assertIn("no companion installation", oracle["provenance_statement"])

    def test_git_blob_export_ignores_changed_checkout_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            copied = root.joinpath(*driver.RETAINED_LINUX_REFERENCE_ROOT.parts)
            copied.parent.mkdir(parents=True)
            shutil.copytree(REFERENCE_ROOT, copied)
            subprocess.run(["git", "init", "-q"], cwd=root, check=True)
            subprocess.run(
                ["git", "config", "user.email", "test@example.invalid"],
                cwd=root,
                check=True,
            )
            subprocess.run(
                ["git", "config", "user.name", "Reference Test"],
                cwd=root,
                check=True,
            )
            subprocess.run(
                ["git", "config", "core.autocrlf", "false"],
                cwd=root,
                check=True,
            )
            subprocess.run(
                ["git", "add", driver.RETAINED_LINUX_REFERENCE_ROOT.as_posix()],
                cwd=root,
                check=True,
            )
            subprocess.run(
                ["git", "commit", "-q", "-m", "retained reference"],
                cwd=root,
                check=True,
            )
            source_sha = (
                subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root)
                .decode("ascii")
                .strip()
            )
            relative = "blastn/PesePMNV.MjPMNV.task_blastn/ncbi.raw.out"
            original = self.payloads[relative]
            changed = copied.joinpath(*PurePosixPath(relative).parts)
            changed.write_bytes(original.replace(b"\n", b"\r\n"))
            manifest_bytes, _, payloads = driver._tracked_linux_reference_from_git(
                root, source_sha
            )
            self.assertEqual(manifest_bytes, self.manifest_bytes)
            self.assertEqual(payloads[relative], original)
            self.assertNotEqual(payloads[relative], changed.read_bytes())


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
    def test_pair_capture_preserves_raw_bytes_and_authoritative_runs_first(
        self,
    ) -> None:
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
        self.assertEqual(
            driver.standard_projection(diagnostic_rows), authoritative_rows
        )
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
        records = driver.changed_hsp_diagnostic_records(comparison, diagnostic_rows)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["platform_diagnostic"]["score"], "40")
        self.assertFalse(records[0]["diagnostic_used_for_authoritative_pass_fail"])

    def test_alternative_representation_requires_distinct_rich_evidence(self) -> None:
        expected = {"windows-x64", "macos-arm64", "macos-x64"}
        common = {"btop": "10", "qseq": "ACGT", "sseq": "ACGT"}
        rows = [
            {"platform_id": platform_id, **common} for platform_id in sorted(expected)
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
        finding, distinct_count, _ = driver.rich_representation_finding(rows, expected)
        self.assertIn("ALTERNATIVE_EQUAL_SCORING", finding)
        self.assertEqual(distinct_count, 2)


# NCBI references:
# - ncbi-blast/c++/src/objtools/align_format/format_flags.cpp:38-45,
#   57-68, 86-101, 114-119, 144-146 defines the invariant and rich fields.
# - ncbi-blast/c++/src/objtools/align_format/tabular.cpp:762-773 obtains score,
#   E-value, and bit score; 778-840 obtains IDs; 900-1027 derives alignment
#   length, qseq, sseq, and BTOP; 1031-1069 derives endpoints.
class ExhaustiveRichComparisonTests(unittest.TestCase):
    def diagnostic_row(self, row_index=1, **overrides):
        row = {
            "row_index": row_index,
            "query_id": "q",
            "subject_id": "s",
            "pident": "100.000",
            "length": "10",
            "mismatch": "0",
            "gapopen": "0",
            "qstart": "1",
            "qend": "10",
            "sstart": "1",
            "send": "10",
            "evalue": "1e-5",
            "bitscore": "20",
            "score": "40",
            "btop": "10",
            "qseq": "ACGT",
            "sseq": "ACGT",
        }
        row.update(overrides)
        return row

    def groups(self, rows_by_platform, triggers=None):
        return driver.exhaustive_rich_representation_groups(
            "blastn", "case", rows_by_platform, triggers or {}
        )

    def three_platform_rich_difference(self, **identity_overrides):
        common = self.diagnostic_row(**identity_overrides)
        alternative = self.diagnostic_row(
            btop="4AC5",
            qseq="AC-GT",
            sseq="ACTGT",
            **identity_overrides,
        )
        return {
            "platform-a": [alternative],
            "platform-b": [common],
            "platform-c": [common],
        }

    def test_rich_only_difference_needs_no_standard_field_seed(self) -> None:
        groups = self.groups(self.three_platform_rich_difference())
        self.assertEqual(len(groups), 1)
        group = groups[0]
        self.assertFalse(group["standard_authoritative_fields_changed"])
        self.assertFalse(group["authoritative_seeded"])
        self.assertTrue(group["rich_only"])
        self.assertFalse(group["diagnostic_used_for_authoritative_pass_fail"])
        self.assertEqual(
            set(group["alignment_identity"]),
            {
                "program",
                "case_id",
                *driver.RICH_INVARIANT_FIELDS,
                "raw_score",
            },
        )

    def test_authoritative_seeded_rich_difference_remains_distinct(self) -> None:
        rows = self.three_platform_rich_difference()
        identity = driver.alignment_representation_identity(
            "blastn", "case", rows["platform-a"][0]
        )
        groups = self.groups(rows, {identity: [{"field": "pident"}]})
        self.assertEqual(len(groups), 1)
        self.assertTrue(groups[0]["standard_authoritative_fields_changed"])
        self.assertTrue(groups[0]["authoritative_seeded"])
        self.assertFalse(groups[0]["rich_only"])
        self.assertEqual(
            groups[0]["triggering_authoritative_differences"],
            [{"field": "pident"}],
        )

    def test_identical_rich_rows_are_not_reported(self) -> None:
        self.assertEqual(
            self.groups(
                {
                    platform_id: [self.diagnostic_row()]
                    for platform_id in ("platform-a", "platform-b", "platform-c")
                }
            ),
            [],
        )

    def test_duplicate_invariant_identities_preserve_multiplicity(self) -> None:
        duplicate_rows = [self.diagnostic_row(1), self.diagnostic_row(2)]
        groups = self.groups(
            {
                "platform-a": duplicate_rows,
                "platform-b": [self.diagnostic_row()],
                "platform-c": [self.diagnostic_row()],
            }
        )
        self.assertEqual(len(groups), 1)
        evidence = {
            record["platform_id"]: record
            for record in groups[0]["platform_representation_multisets"]
        }
        self.assertEqual(evidence["platform-a"]["row_count"], 2)
        self.assertEqual(
            evidence["platform-a"]["representations"][0]["multiplicity"], 2
        )
        self.assertEqual(
            evidence["platform-b"]["representations"][0]["multiplicity"], 1
        )

    def test_row_order_only_differences_do_not_create_groups(self) -> None:
        first = self.diagnostic_row(btop="10", qseq="ACGT", sseq="ACGT")
        second = self.diagnostic_row(
            row_index=2, btop="4AC5", qseq="AC-GT", sseq="ACTGT"
        )
        self.assertEqual(
            self.groups(
                {
                    "platform-a": [first, second],
                    "platform-b": [second, first],
                    "platform-c": [first, second],
                }
            ),
            [],
        )

    def test_differing_score_prevents_rich_grouping(self) -> None:
        rows = self.three_platform_rich_difference()
        rows["platform-a"][0]["score"] = "41"
        self.assertEqual(self.groups(rows), [])

    def test_differing_endpoint_prevents_rich_grouping(self) -> None:
        for field, value in (
            ("qstart", "2"),
            ("qend", "9"),
            ("sstart", "2"),
            ("send", "9"),
        ):
            with self.subTest(field=field):
                rows = self.three_platform_rich_difference()
                rows["platform-a"][0][field] = value
                self.assertEqual(self.groups(rows), [])

    def test_differing_length_evalue_or_bitscore_prevents_grouping(self) -> None:
        for field, value in (
            ("length", "11"),
            ("evalue", "2e-5"),
            ("bitscore", "21"),
        ):
            with self.subTest(field=field):
                rows = self.three_platform_rich_difference()
                rows["platform-a"][0][field] = value
                self.assertEqual(self.groups(rows), [])

    def test_platform_relationship_is_derived_without_named_divergent_host(
        self,
    ) -> None:
        group = self.groups(self.three_platform_rich_difference())[0]
        classes = [
            item["platform_ids"] for item in group["platform_equivalence_classes"]
        ]
        self.assertEqual(classes, [["platform-a"], ["platform-b", "platform-c"]])
        platform_hashes = {
            item["platform_id"]: item["representation_multiset_sha256"]
            for item in group["platform_representation_multisets"]
        }
        for equivalence_class in group["platform_equivalence_classes"]:
            self.assertEqual(
                {
                    platform_hashes[platform_id]
                    for platform_id in equivalence_class["platform_ids"]
                },
                {equivalence_class["representation_multiset_sha256"]},
            )
        self.assertNotIn("windows", json.dumps(group).lower())

    def test_seeded_rich_only_and_exhaustive_counts_are_distinct(self) -> None:
        first = self.three_platform_rich_difference(qstart="1", qend="10")
        second = self.three_platform_rich_difference(qstart="21", qend="30")
        rows = {
            platform_id: [first[platform_id][0], second[platform_id][0]]
            for platform_id in first
        }
        seeded_identity = driver.alignment_representation_identity(
            "blastn", "case", first["platform-a"][0]
        )
        counts = driver.rich_representation_group_counts(
            self.groups(rows, {seeded_identity: [{"field": "pident"}]})
        )
        self.assertEqual(
            counts,
            {
                "authoritative_seeded_equal_score_representation_group_count": 1,
                "rich_only_equal_score_representation_group_count": 1,
                "exhaustive_equal_score_representation_group_count": 2,
            },
        )

    def test_linux_rich_evidence_is_explicitly_unavailable(self) -> None:
        linux = b"q\ts\t100\t10\t0\t0\t1\t10\t1\t10\t1e-5\t20\n"
        evidence = driver.retained_linux_rich_evidence(linux)
        self.assertFalse(evidence["rich_comparison_available"])
        self.assertFalse(evidence["rich_fields_available"])
        self.assertEqual(
            evidence["unavailable_fields"], ["score", "btop", "qseq", "sseq"]
        )
        self.assertFalse(evidence["inference_or_synthesis_performed"])
        self.assertFalse(evidence["rerun_performed"])


class WorkflowContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.authority, _ = driver.load_platform_authority(REPO_ROOT)

    def _job_run_text(self, job_name: str) -> str:
        document = yaml.load(
            WORKFLOW_PATH.read_text(encoding="utf-8"), Loader=yaml.BaseLoader
        )
        return "\n".join(
            step.get("run", "") for step in document["jobs"][job_name]["steps"]
        )

    def assert_no_blast_or_losat_invocation(self, run_text: str) -> None:
        self.assertIsNone(
            re.search(
                r"(?:^|\s|[;&|])(?:[^\s]*/)?"
                r"(?:blastn|blastp|tblastx|blast_formatter)(?=\s|$)",
                run_text,
            )
        )
        self.assertIsNone(
            re.search(
                r"(?:^|\s|[;&|])(?:[^\s]*/)?LOSAT(?:\.exe|\.wasm)?(?=\s|$)",
                run_text,
            )
        )

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
        self.assertEqual(analysis["needs"], ["characterize", "retained-linux-evidence"])
        text = WORKFLOW_PATH.read_text(encoding="utf-8")
        self.assertIn("analyze-platforms", text)
        self.assertEqual(text.count("--platform-evidence"), 3)
        self.assertIn("18 authoritative plus 18 diagnostic", text)
        self.assertNotIn("linux-rerun", text.lower())

    def test_retained_linux_export_is_hosted_git_blob_transport_only(self) -> None:
        document = yaml.load(
            WORKFLOW_PATH.read_text(encoding="utf-8"), Loader=yaml.BaseLoader
        )
        job = document["jobs"]["retained-linux-evidence"]
        self.assertEqual(job["runs-on"], "ubuntu-24.04")
        serialized = json.dumps(job, sort_keys=True)
        self.assertNotIn("self-hosted", serialized)
        self.assertNotIn("/mnt/c/", serialized)
        self.assertNotIn("LINUX_EVIDENCE_ROOT", serialized)
        self.assertIn("export-linux-reference", serialized)
        self.assertIn("--candidate-sha", serialized)
        self.assertIn("ncbi-characterization-retained-linux", serialized)
        run_text = self._job_run_text("retained-linux-evidence")
        self.assertNotIn("curl ", run_text)
        self.assertNotIn("cargo ", run_text)
        self.assertNotIn("run_losat", run_text)
        self.assertNotIn("--native-bin", run_text)
        self.assert_no_blast_or_losat_invocation(run_text)

    def test_analysis_job_cannot_execute_blast_or_losat(self) -> None:
        run_text = self._job_run_text("analyze-rich-diagnostics")
        self.assertNotIn("cargo ", run_text)
        self.assertNotIn("run_losat", run_text)
        self.assertNotIn("--native-bin", run_text)
        self.assert_no_blast_or_losat_invocation(run_text)

    def test_main_bootstrap_and_protected_contracts_are_unmodified(self) -> None:
        changed = set(
            subprocess.check_output(
                ["git", "diff", "--name-only", driver.PR6_CERT_SHA_V3],
                cwd=REPO_ROOT,
            )
            .decode("utf-8")
            .splitlines()
        )
        protected = {
            "LOSAT/src/main.rs",
            "LOSAT/tests/platform_native_v010_canonical.tsv",
            ".github/workflows/platform-native-certification-v010.yml",
            "LOSAT/tests/certify_platform_native_v010.py",
            "LOSAT/tests/test_certify_platform_native_v010.py",
        }
        self.assertTrue(changed.isdisjoint(protected))
        self.assertFalse(
            any(path.startswith("docs/product_decisions/") for path in changed)
        )
        self.assertFalse(any(path.startswith(".agents/plans/") for path in changed))

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
