from __future__ import annotations

# NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:828-832
# ```c++
# ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) {
#     tabinfo.SetFields(**itr, *m_Scope, &m_ScoringMatrix);
#     tabinfo.Print();
# }
# ```
# These tests protect the raw output and existing-authority orchestration boundary.

from collections import Counter
from dataclasses import replace
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest import mock


SCRIPT = Path(__file__).with_name("certify_platform_native_v010.py")
SPEC = importlib.util.spec_from_file_location("certify_platform_native_v010", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {SCRIPT}")
platform_cert = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = platform_cert
SPEC.loader.exec_module(platform_cert)


class PlatformCertificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo_root = Path(__file__).resolve().parents[2]
        cls.catalog = platform_cert.load_catalog(cls.repo_root)

    def test_frozen_manifest_is_exactly_the_pr5_43_row_authority(self) -> None:
        self.assertEqual(len(self.catalog.canonical_rows), 43)
        self.assertEqual(
            Counter(row["program"] for row in self.catalog.canonical_rows),
            Counter({"blastn": 14, "blastp": 9, "tblastx": 20}),
        )
        self.assertEqual(
            Counter(
                row["classification"]
                for row in self.catalog.canonical_rows
                if row["program"] == "blastn"
            ),
            Counter({"EXACT_TEXT": 13, "SOURCE_UNDETERMINED_ACCEPTED": 1}),
        )
        sakai = self.catalog.canonical[("blastn", "Sakai.MG1655.megablast")]
        self.assertEqual(sakai["classification"], "SOURCE_UNDETERMINED_ACCEPTED")
        self.assertEqual(sakai["contract"], "SOURCE_UNDETERMINED_ACCEPTED")

    def test_representatives_are_exactly_the_six_pr6_cases(self) -> None:
        self.assertEqual(
            platform_cert.REPRESENTATIVES,
            (
                ("blastn", "PesePMNV.MjPMNV.task_blastn", "ordinary_exact"),
                ("blastn", "Sakai.MG1655.megablast", "source_undetermined"),
                ("blastp", "pairwise_default_serial", "representative"),
                ("tblastx", "p03_mela_pemojnva", "ordinary_exact"),
                (
                    "tblastx",
                    "p14_ap027131_ap027133_query4",
                    "query_gencode_exact",
                ),
                (
                    "tblastx",
                    "d06_ap027131_ap027133_db4",
                    "db_gencode_deviation",
                ),
            ),
        )
        self.assertNotIn(
            ("tblastx", "p11_avclpv_psclpv"),
            {(program, case_id) for program, case_id, _ in platform_cert.REPRESENTATIVES},
        )

    def test_plan_is_exactly_61_executions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            steps = platform_cert.build_steps(
                self.repo_root,
                root,
                self.catalog,
                root / "LOSAT",
                {
                    "blastn": root / "blastn",
                    "blastp": root / "blastp",
                    "tblastx": root / "tblastx",
                },
            )
        self.assertEqual(len(steps), 61)
        self.assertEqual(
            Counter(step.kind for step in steps),
            Counter({"matrix": 43, "oracle": 6, "repeatability": 12}),
        )
        self.assertEqual(
            [step.execution_index for step in steps], list(range(1, 62))
        )
        matrix_keys = {
            (step.program, step.case_id) for step in steps if step.kind == "matrix"
        }
        self.assertEqual(matrix_keys, set(self.catalog.canonical))
        oracle_keys = {
            (step.program, step.case_id) for step in steps if step.kind == "oracle"
        }
        self.assertEqual(
            oracle_keys,
            {
                (program, case_id)
                for program, case_id, _ in platform_cert.REPRESENTATIVES
            },
        )

    def test_planned_fixture_paths_use_pr5_root_without_resolved_path_leaks(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            steps = platform_cert.build_steps(
                self.repo_root,
                root,
                self.catalog,
                root / "LOSAT",
                {
                    "blastn": root / "blastn",
                    "blastp": root / "blastp",
                    "tblastx": root / "tblastx",
                },
            )
        self.assertEqual(
            platform_cert.HISTORICAL_LEXICAL_ROOT,
            "/tmp/losat-pr5-runtime-cert-5845d22/LOSAT",
        )
        matrix = [step for step in steps if step.kind == "matrix"]
        self.assertEqual(len(matrix), 43)
        blastn = [step for step in matrix if step.program == "blastn"]
        path_sensitive = {
            step.case_id
            for step in blastn
            if platform_cert._option_value(step.command, "--outfmt") == "7"
        }
        self.assertEqual(
            path_sensitive,
            {
                "PesePMNV.MjPMNV.task_blastn",
                "PmeNMV.MjPMNV.task_blastn",
                "PmeNMV.PesePMNV.task_blastn",
                "PeseMJNV.PemoMJNVB.task_blastn",
                "MelaMJNV.PemoMJNVA.task_blastn",
                "PemoMJNVA.PeseMJNV.task_blastn",
                "MjeNMV.MelaMJNV.task_blastn",
                "MjPMNV.MlPMNV.task_blastn",
                "NZ_CP006932.NZ_CP006932.task_blastn",
                "EDL933.Sakai.megablast",
                "Sakai.MG1655.megablast",
                "compact.multi_query.outfmt7",
                "compact.multi_query.no_hit.outfmt7",
            },
        )
        headerless = [step for step in blastn if step.case_id not in path_sensitive]
        self.assertEqual(
            [step.case_id for step in headerless], ["compact.multi_query.outfmt6"]
        )
        prefix = f"{platform_cert.HISTORICAL_LEXICAL_ROOT}/"
        for step in steps:
            for _, _, value in platform_cert._input_arguments(
                step.command, step.program
            ):
                self.assertTrue(value.startswith(prefix))
                self.assertNotIn(str(self.repo_root), value)
                self.assertNotIn("\\", value)
        for step in matrix:
            if step.program in {"blastp", "tblastx"}:
                self.assertEqual(
                    platform_cert._option_value(step.command, "--outfmt"), "6"
                )
                self.assertEqual(
                    step.expected_losat_sha256,
                    self.catalog.canonical[(step.program, step.case_id)][
                        "losat_sha256"
                    ],
                )
        self.assertEqual(
            len(platform_cert.required_fixture_relatives(steps)),
            platform_cert.EXPECTED_STAGED_FIXTURE_COUNT,
        )

    def test_fixture_staging_preserves_source_bytes_and_records_both_paths(
        self,
    ) -> None:
        relative = Path("tests/fasta/platform_cert_stage_unit_fixture.fa")
        lexical = platform_cert.historical_lexical_path(relative.as_posix())
        physical = platform_cert.controlled_physical_path(lexical)
        payload = b">fixture\r\nACGTN\r\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "LOSAT" / relative
            source.parent.mkdir(parents=True)
            source.write_bytes(payload)
            output = root / "evidence"
            output.mkdir()
            step = platform_cert.SearchStep(
                execution_index=1,
                step_id="matrix:blastn:fixture",
                kind="matrix",
                program="blastn",
                case_id="fixture",
                semantic_class="EXACT_TEXT",
                run_index=1,
                command=(
                    "LOSAT",
                    "blastn",
                    "-q",
                    lexical,
                    "-s",
                    lexical,
                    "--outfmt",
                    "7",
                    "-o",
                    "unused.out",
                ),
                environment=(),
                output_rel="unused.out.partial",
                final_output_rel="unused.out",
                expected_losat_sha256="0" * 64,
            )
            try:
                records = platform_cert.stage_required_fixtures(root, [step], output)
                self.assertEqual(len(records), 1)
                self.assertEqual(physical.read_bytes(), payload)
                self.assertEqual(records[0].lexical_target_path, lexical)
                self.assertEqual(records[0].physical_target_path, str(physical))
                self.assertEqual(records[0].sha256, hashlib.sha256(payload).hexdigest())
                evidence = json.loads(
                    (output / "fixture-staging.json").read_text(encoding="utf-8")
                )
                self.assertEqual(evidence["status"], "VERIFIED")
                self.assertEqual(evidence["staged_input_count"], 1)
            finally:
                physical.unlink(missing_ok=True)

    def test_physical_and_lexical_controlled_paths_are_not_conflated(self) -> None:
        lexical = platform_cert.historical_lexical_path("tests/fasta/fixture.fa")
        with mock.patch.object(
            platform_cert.os.path,
            "abspath",
            return_value=r"D:\tmp\losat-pr5-runtime-cert-5845d22\LOSAT\tests\fasta\fixture.fa",
        ):
            physical = platform_cert.controlled_physical_path(lexical)
        self.assertEqual(
            lexical,
            "/tmp/losat-pr5-runtime-cert-5845d22/LOSAT/tests/fasta/fixture.fa",
        )
        self.assertNotEqual(str(physical), lexical)

    def test_controlled_path_preflight_proves_literal_shell_free_argv(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            record = platform_cert.preflight_controlled_fixture_path(output)
            evidence = json.loads(
                (output / "controlled-path-preflight.json").read_text(encoding="utf-8")
            )
        self.assertEqual(record["status"], "VERIFIED")
        self.assertEqual(
            evidence["child_observation"]["argv"], record["lexical_sentinel_path"]
        )
        self.assertFalse(evidence["subprocess_contract"]["shell"])
        self.assertFalse(evidence["subprocess_contract"]["msys_or_git_bash_involved"])

    def test_controlled_path_preflight_is_fail_closed_before_search(self) -> None:
        lexical = platform_cert.historical_lexical_path(
            ".platform-native-certification-path-preflight"
        )
        altered = json.dumps(
            {
                "argv": r"D:\tmp\rewritten-by-a-shell",
                "sha256": hashlib.sha256(
                    b"LOSAT_PLATFORM_CERTIFICATION_PATH_PREFLIGHT\n"
                ).hexdigest(),
            }
        ).encode()
        completed = platform_cert.subprocess.CompletedProcess(
            args=[sys.executable, "-c", "<probe>", lexical],
            returncode=0,
            stdout=altered,
            stderr=b"",
        )
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            with mock.patch.object(
                platform_cert.subprocess, "run", return_value=completed
            ):
                with self.assertRaises(platform_cert.CertificationFailure):
                    platform_cert.preflight_controlled_fixture_path(output)
            evidence = json.loads(
                (output / "controlled-path-preflight.json").read_text(encoding="utf-8")
            )
        self.assertEqual(evidence["status"], "FAILED")
        self.assertIn("did not preserve", evidence["error"])

    def test_matrix_hashing_preserves_raw_output_without_rewriting(self) -> None:
        payload = (
            b"# BLASTN 2.17.0+\r\n"
            b"# Database: /tmp/losat-pr5-runtime-cert-5845d22/LOSAT/tests/fasta/s.fa\r\n"
            b"q\ts\t100.000\r\n"
        )
        raw_hash = hashlib.sha256(payload).hexdigest()
        self.assertNotEqual(
            raw_hash,
            hashlib.sha256(platform_cert.canonical_newlines(payload)).hexdigest(),
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output_rel = "matrix/blastn/fixture/losat.out.partial"
            partial = root / output_rel
            writer = (
                "from pathlib import Path;import sys;"
                "Path(sys.argv[1]).parent.mkdir(parents=True,exist_ok=True);"
                "Path(sys.argv[1]).write_bytes(bytes.fromhex(sys.argv[2]))"
            )
            step = platform_cert.SearchStep(
                execution_index=1,
                step_id="matrix:blastn:fixture",
                kind="matrix",
                program="blastn",
                case_id="fixture",
                semantic_class="EXACT_TEXT",
                run_index=1,
                command=(sys.executable, "-c", writer, str(partial), payload.hex()),
                environment=(),
                output_rel=output_rel,
                final_output_rel=output_rel.removesuffix(".partial"),
                expected_losat_sha256=raw_hash,
            )
            identity = {
                "source_sha": "a" * 40,
                "platform_id": "test",
                "losat": {"sha256": "b" * 64},
            }
            record = platform_cert._execute_step(
                root, root, identity, self.catalog, step
            )
            final_output = root / step.final_output_rel
            self.assertEqual(final_output.read_bytes(), payload)
            self.assertIn(b"# Database:", final_output.read_bytes())
            self.assertEqual(record["output_sha256"], raw_hash)
            self.assertEqual(
                record["verification"]["classification"],
                "CANONICAL_PR5_RAW_BYTES",
            )

    def test_python_tarfile_extracts_only_after_pinned_verification(self) -> None:
        payload = b"official oracle fixture\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / "fixture.tar.gz"
            with tarfile.open(archive, mode="w:gz") as bundle:
                member = tarfile.TarInfo("ncbi-blast-2.17.0+/bin/blastn")
                member.size = len(payload)
                bundle.addfile(member, io.BytesIO(payload))
            digest = hashlib.md5(
                archive.read_bytes(), usedforsecurity=False
            ).hexdigest()
            spec = replace(
                platform_cert.PLATFORMS["macos-arm64"],
                archive_name=archive.name,
                archive_md5=digest,
            )
            checksum = root / f"{archive.name}.md5"
            destination = root / "oracle"
            checksum.write_text(f"{'0' * 32}  {archive.name}\n", encoding="ascii")
            with self.assertRaises(platform_cert.CertificationFailure):
                platform_cert.extract_verified_archive(
                    spec, archive, checksum, destination
                )
            self.assertFalse(destination.exists())
            checksum.write_text(f"{digest}  {archive.name}\n", encoding="ascii")
            record = platform_cert.extract_verified_archive(
                spec, archive, checksum, destination
            )
            self.assertEqual(
                (destination / "ncbi-blast-2.17.0+/bin/blastn").read_bytes(),
                payload,
            )
            self.assertEqual(record["status"], "EXTRACTED_VERIFIED_ARCHIVE")
            self.assertEqual(record["member_count"], 1)

    def test_archive_requires_pinned_published_record_and_content(self) -> None:
        content = b"official archive fixture"
        digest = hashlib.md5(content, usedforsecurity=False).hexdigest()
        spec = replace(
            platform_cert.PLATFORMS["macos-arm64"],
            archive_name="fixture.tar.gz",
            archive_md5=digest,
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            archive = root / spec.archive_name
            checksum = root / f"{spec.archive_name}.md5"
            archive.write_bytes(content)
            checksum.write_text(
                f"{digest}  {spec.archive_name}\n", encoding="ascii"
            )
            record = platform_cert.verify_archive(spec, archive, checksum)
            self.assertEqual(record["archive_md5"], digest)
            checksum.write_text(
                f"{'0' * 32}  {spec.archive_name}\n", encoding="ascii"
            )
            with self.assertRaises(platform_cert.CertificationFailure):
                platform_cert.verify_archive(spec, archive, checksum)

    def test_line_ending_diagnostic_is_narrow(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            crlf = root / "crlf.out"
            lf = root / "lf.out"
            biological_difference = root / "different.out"
            crlf.write_bytes(b"a\tb\r\nc\td\r\n")
            lf.write_bytes(b"a\tb\nc\td\n")
            biological_difference.write_bytes(b"a\tb\nc\te\n")
            self.assertTrue(platform_cert.line_endings_only(crlf, lf))
            self.assertFalse(
                platform_cert.line_endings_only(crlf, biological_difference)
            )
            self.assertFalse(platform_cert.line_endings_only(lf, lf))

    def test_blastn_direct_check_accepts_only_crlf_diagnostic(self) -> None:
        row = "q\ts\t100.000\t4\t0\t0\t1\t4\t1\t4\t1e-10\t20.0"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ncbi = root / "ncbi.out"
            losat = root / "losat.out"
            ncbi.write_bytes((row + "\r\n").encode())
            losat.write_bytes((row + "\n").encode())
            step = platform_cert.SearchStep(
                execution_index=44,
                step_id="oracle:blastn:PesePMNV.MjPMNV.task_blastn",
                kind="oracle",
                program="blastn",
                case_id="PesePMNV.MjPMNV.task_blastn",
                semantic_class="ordinary_exact",
                run_index=1,
                command=("blastn",),
                environment=(),
                output_rel="oracle.out.partial",
                final_output_rel="oracle.out",
                expected_losat_sha256=hashlib.sha256(losat.read_bytes()).hexdigest(),
            )
            result = platform_cert.classify_oracle(
                self.catalog, step, ncbi, losat, root / "diagnostic"
            )
            self.assertEqual(result["authority_classification"], "EXACT_DATA")
            self.assertEqual(result["platform_diagnostic"], "LINE_ENDINGS_ONLY")

            ncbi.write_bytes((row.replace("20.0", "21.0") + "\r\n").encode())
            with self.assertRaises(platform_cert.CertificationFailure):
                platform_cert.classify_oracle(
                    self.catalog, step, ncbi, losat, root / "rejected"
                )

    def test_binary_architecture_reads_pe_and_macho_headers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pe = bytearray(128)
            pe[:2] = b"MZ"
            pe[60:64] = (64).to_bytes(4, "little")
            pe[64:68] = b"PE\0\0"
            pe[68:70] = (0x8664).to_bytes(2, "little")
            pe_path = root / "LOSAT.exe"
            pe_path.write_bytes(pe)
            self.assertEqual(platform_cert.binary_architecture(pe_path), "x86_64")

            macho_path = root / "LOSAT"
            macho_path.write_bytes(
                b"\xcf\xfa\xed\xfe" + (0x0100000C).to_bytes(4, "little") + b"\0" * 64
            )
            self.assertEqual(platform_cert.binary_architecture(macho_path), "arm64")

    def test_resume_identity_binds_source_binary_oracles_and_manifest(self) -> None:
        identity = {
            "evidence_schema": 1,
            "source_sha": "a" * 40,
            "platform_id": "macos-arm64",
            "runner_label": "macos-15",
            "toolchain": {"cargo_V": "cargo 1.92.0"},
            "losat": {"sha256": "b" * 64},
            "oracle_archive": {"archive_md5": "c" * 32},
            "oracles": {
                "blastn": {"sha256": "d" * 64},
                "blastp": {"sha256": "e" * 64},
                "tblastx": {"sha256": "f" * 64},
            },
            "canonical_manifest": {"sha256": "1" * 64},
            "controlled_fixtures": {
                "lexical_root": platform_cert.HISTORICAL_LEXICAL_ROOT,
                "staged_input_count": 34,
                "inputs": [],
            },
        }
        key = platform_cert._identity_resume_key(identity)
        changed = json.loads(json.dumps(identity))
        changed["losat"]["sha256"] = "2" * 64
        self.assertNotEqual(key, platform_cert._identity_resume_key(changed))
        changed = json.loads(json.dumps(identity))
        changed["controlled_fixtures"]["lexical_root"] = "/tmp/not-authoritative"
        self.assertNotEqual(key, platform_cert._identity_resume_key(changed))

    def test_manual_workflow_has_only_the_three_pr6_platform_jobs(self) -> None:
        workflow = (
            self.repo_root / ".github" / "workflows" / "platform-native-certification.yml"
        ).read_text(encoding="utf-8")
        self.assertIn("workflow_dispatch:", workflow)
        self.assertNotIn("pull_request:", workflow)
        self.assertNotIn("schedule:", workflow)
        self.assertIn("max-parallel: 3", workflow)
        self.assertEqual(workflow.count("platform_id:"), 3)
        for runner in ("windows-2025", "macos-15", "macos-15-intel"):
            self.assertIn(f"runner: {runner}", workflow)
        for excluded in ("ubuntu-latest", "wasm32", "benchmark", "soak"):
            self.assertNotIn(excluded, workflow.lower())
        self.assertIn("--expected-sha \"$CANDIDATE_SHA\"", workflow)
        self.assertIn("preflight-path", workflow)
        self.assertIn("extract-verified-archive", workflow)
        self.assertNotIn("tar -xzf", workflow)
        self.assertIn("actions/upload-artifact@v4", workflow)
        self.assertIn("if: ${{ always() }}", workflow)

    def test_resume_import_revalidates_verified_output_hash(self) -> None:
        payload = b"canonical output\n"
        output_hash = hashlib.sha256(payload).hexdigest()
        identity = {
            "evidence_schema": 1,
            "source_sha": "a" * 40,
            "platform_id": "macos-arm64",
            "runner_label": "macos-15",
            "toolchain": {"cargo_V": "cargo 1.92.0"},
            "losat": {"sha256": "b" * 64},
            "oracle_archive": {"archive_md5": "c" * 32},
            "oracles": {
                "blastn": {"sha256": "d" * 64},
                "blastp": {"sha256": "e" * 64},
                "tblastx": {"sha256": "f" * 64},
            },
            "canonical_manifest": {"sha256": "1" * 64},
            "controlled_fixtures": {
                "lexical_root": platform_cert.HISTORICAL_LEXICAL_ROOT,
                "staged_input_count": 34,
                "inputs": [],
            },
        }
        step = platform_cert.SearchStep(
            execution_index=1,
            step_id="matrix:blastn:fixture",
            kind="matrix",
            program="blastn",
            case_id="fixture",
            semantic_class="EXACT_TEXT",
            run_index=1,
            command=("LOSAT",),
            environment=(),
            output_rel="matrix/blastn/fixture/losat.out.partial",
            final_output_rel="matrix/blastn/fixture/losat.out",
            expected_losat_sha256=output_hash,
        )
        record = {
            "status": "VERIFIED",
            "step_id": step.step_id,
            "execution_index": 1,
            "kind": "matrix",
            "program": "blastn",
            "case_id": "fixture",
            "run_index": 1,
            "expected_losat_sha256": output_hash,
            "source_sha": identity["source_sha"],
            "platform_id": identity["platform_id"],
            "losat_binary_sha256": identity["losat"]["sha256"],
            "output_rel": step.final_output_rel,
            "output_sha256": output_hash,
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prior = root / "prior"
            current = root / "current"
            prior.mkdir()
            current.mkdir()
            (prior / "identity.json").write_text(
                json.dumps(identity), encoding="utf-8"
            )
            prior_output = prior / step.final_output_rel
            prior_output.parent.mkdir(parents=True)
            prior_output.write_bytes(payload)
            completion = platform_cert._step_record_path(prior, step)
            completion.parent.mkdir(parents=True)
            completion.write_text(json.dumps(record), encoding="utf-8")
            imported = platform_cert.import_resume(
                prior, current, identity, [step], self.catalog
            )
            self.assertEqual(set(imported), {step.step_id})
            self.assertEqual((current / step.final_output_rel).read_bytes(), payload)

            prior_output.write_bytes(b"corrupt\n")
            second = root / "second"
            second.mkdir()
            with self.assertRaises(platform_cert.CertificationFailure):
                platform_cert.import_resume(
                    prior, second, identity, [step], self.catalog
                )


if __name__ == "__main__":
    unittest.main()
