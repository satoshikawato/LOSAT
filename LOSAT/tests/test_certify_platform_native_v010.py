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
import copy
import contextlib
from dataclasses import replace
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
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
        cls.candidate_sha = platform_cert.resolve_git_head_sha(cls.repo_root)
        cls.catalog = platform_cert.load_catalog(cls.repo_root, cls.candidate_sha)
        cls.authority_path = SCRIPT.with_name("ncbi_platform_variance_v010.json")
        cls.authority = platform_cert.load_native_authority(
            cls.repo_root, cls.candidate_sha, cls.authority_path
        )
        cls.authority_document = copy.deepcopy(cls.authority.document)

    def assert_authority_rejected(self, mutate) -> None:
        document = copy.deepcopy(self.authority_document)
        mutate(document)
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "authority.json"
            path.write_text(json.dumps(document), encoding="utf-8")
            with self.assertRaises(platform_cert.CertificationFailure):
                platform_cert._load_native_authority_bytes(path.read_bytes(), path)

    def commit_files(self, root: Path, files: dict[str, bytes]) -> str:
        subprocess.run(["git", "init", "--quiet"], cwd=root, check=True)
        for relative, payload in files.items():
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
        subprocess.run(["git", "add", "--all"], cwd=root, check=True)
        subprocess.run(
            [
                "git",
                "-c",
                "user.name=LOSAT Certifier Test",
                "-c",
                "user.email=certifier-test@example.invalid",
                "commit",
                "--quiet",
                "-m",
                "Create immutable-byte fixture",
            ],
            cwd=root,
            check=True,
        )
        return platform_cert.resolve_git_head_sha(root)

    def authority_for_payload(
        self,
        platform_id: str,
        program: str,
        case_id: str,
        payload: bytes,
        **changes,
    ):
        outputs = copy.deepcopy(self.authority.outputs)
        key = (platform_id, program, case_id)
        outputs[key] = {
            "program": program,
            "case_id": case_id,
            "raw_sha256": hashlib.sha256(payload).hexdigest(),
            "byte_length": len(payload),
            "newline_profile": platform_cert.newline_profile(payload),
            "data_row_count": platform_cert.output_data_row_count(payload),
            **changes,
        }
        return replace(self.authority, outputs=outputs)

    def native_invocation_identity(self, step, platform_id: str) -> dict[str, object]:
        platform_record = self.authority.platforms[platform_id]
        archive = platform_record["archive"]
        checksum = archive["checksum"]
        executable = next(
            item
            for item in platform_record["executables"]
            if item["program"] == step.program
        )
        return {
            "source_sha": self.candidate_sha,
            "platform_id": platform_id,
            "authority_version": self.authority.authority_version,
            "authority_file_sha256": self.authority.file_sha256,
            "oracle_archive": {
                "archive_name": archive["filename"],
                "archive_checksum_algorithm": checksum["algorithm"],
                "archive_md5": checksum["value"],
            },
            "oracles": {
                step.program: {
                    "path": step.command[0],
                    "sha256": executable["sha256"],
                    "architecture": executable["architecture"],
                }
            },
        }

    def create_aggregate_fixture(self, root: Path):
        canonical_rows = copy.deepcopy(self.catalog.canonical_rows)
        for row in canonical_rows:
            payload = f"matrix:{row['program']}:{row['case_id']}\n".encode()
            row["losat_sha256"] = hashlib.sha256(payload).hexdigest()
        canonical = {(row["program"], row["case_id"]): row for row in canonical_rows}
        catalog = replace(
            self.catalog, canonical_rows=canonical_rows, canonical=canonical
        )
        outputs = copy.deepcopy(self.authority.outputs)
        for platform_id in sorted(self.authority.platforms):
            for program, case_id, _ in platform_cert.REPRESENTATIVES:
                payload = f"{platform_id}:oracle:{program}:{case_id}\n".encode()
                key = (platform_id, program, case_id)
                outputs[key] = {
                    "program": program,
                    "case_id": case_id,
                    "raw_sha256": hashlib.sha256(payload).hexdigest(),
                    "byte_length": len(payload),
                    "newline_profile": platform_cert.newline_profile(payload),
                    "data_row_count": 1,
                }
        authority = replace(self.authority, outputs=outputs)
        topology = platform_cert.expected_completion_topology(catalog)
        for platform_id in sorted(authority.platforms):
            artifact = root / f"platform-native-{platform_id}"
            (artifact / "completions").mkdir(parents=True)
            losat_sha256 = hashlib.sha256(f"{platform_id}:LOSAT".encode()).hexdigest()
            authority_platform = authority.platforms[platform_id]
            archive = authority_platform["archive"]
            checksum = archive["checksum"]
            oracle_identities = {
                program: {
                    "sha256": next(
                        executable["sha256"]
                        for executable in authority_platform["executables"]
                        if executable["program"] == program
                    ),
                    "architecture": next(
                        executable["architecture"]
                        for executable in authority_platform["executables"]
                        if executable["program"] == program
                    ),
                }
                for program in ("blastn", "blastp", "tblastx")
            }
            platform_cert.atomic_write_json(
                artifact / "identity.json",
                {
                    "source_sha": "a" * 40,
                    "platform_id": platform_id,
                    "authority_version": authority.authority_version,
                    "authority_file_sha256": authority.file_sha256,
                    "os": {
                        "system": authority_platform["platform_system"],
                        "machine": authority_platform["normalized_machine"],
                    },
                    "losat": {
                        "sha256": losat_sha256,
                        "architecture": platform_cert.PLATFORMS[
                            platform_id
                        ].binary_arch,
                    },
                    "oracle_archive": {
                        "archive_name": archive["filename"],
                        "archive_checksum_algorithm": checksum["algorithm"],
                        "archive_md5": checksum["value"],
                    },
                    "oracles": oracle_identities,
                    "execution_contract": {
                        "matrix": 43,
                        "oracle": 6,
                        "repeatability": 12,
                        "per_platform": 61,
                        "gate_a_losat_canonical": 43,
                        "gate_b_native_ncbi_reference": 6,
                    },
                },
            )
            planned = []
            for expected in topology.values():
                planned.append(
                    {
                        **expected,
                        "command": [],
                        "environment": (
                            [["RAYON_NUM_THREADS", "1"]]
                            if expected["kind"] in {"matrix", "repeatability"}
                            and expected["program"] == "tblastx"
                            else []
                        ),
                        "output_rel": f"{expected['output_rel']}.partial",
                        "final_output_rel": expected["output_rel"],
                    }
                )
            platform_cert.atomic_write_json(
                artifact / "command_plan.json",
                {
                    "source_sha": "a" * 40,
                    "platform_id": platform_id,
                    "authority_version": authority.authority_version,
                    "authority_file_sha256": authority.file_sha256,
                    "execution_count": 61,
                    "controlled_fixture_root": platform_cert.HISTORICAL_LEXICAL_ROOT,
                    "executions": planned,
                },
            )
            platform_cert.atomic_write_json(
                artifact / "summary.json",
                {
                    "decision": "PLATFORM_NATIVE_CERTIFIED",
                    "source_sha": "a" * 40,
                    "platform_id": platform_id,
                    "authority_version": authority.authority_version,
                    "authority_file_sha256": authority.file_sha256,
                    "search_executions": 61,
                    "gate_a_losat_canonical": {
                        "total": 43,
                        "passed": 43,
                        "all_exact_pr5_raw_bytes": True,
                    },
                    "gate_b_platform_native_ncbi_reference": {
                        "total": 6,
                        "passed": 6,
                        "all_exact_authority_fingerprints": True,
                    },
                    "repeatability": {
                        "representatives": 6,
                        "runs_per_representative": 3,
                        "additional_executions": 12,
                        "all_repeatable": True,
                    },
                    "sakai_ratchet": {
                        "expected_losat_parity_class": "SOURCE_UNDETERMINED_ACCEPTED",
                        "observed_losat_parity_class": "SOURCE_UNDETERMINED_ACCEPTED",
                        "passed": True,
                    },
                    "tblastx_deviation_ratchet": {
                        "accepted_deviation_ceiling": 6,
                        "d06_expected_losat_contract": "approved_db_gencode_deviation",
                        "d06_expected_losat_parity_class": "HSP_SET_DIFF",
                        "d06_observed_losat_parity_class": "HSP_SET_DIFF",
                        "passed": True,
                    },
                },
            )
            platform_cert.atomic_write_json(
                artifact / "state.json",
                {
                    "platform_id": platform_id,
                    "status": "CERTIFIED",
                    "completed_execution_count": 61,
                    "completed_by_kind": {
                        "matrix": 43,
                        "oracle": 6,
                        "repeatability": 12,
                    },
                    "completed_step_ids": sorted(topology),
                    "expected_execution_count": 61,
                },
            )
            for step_id, expected in topology.items():
                if expected["kind"] == "oracle":
                    payload = (
                        f"{platform_id}:oracle:{expected['program']}:"
                        f"{expected['case_id']}\n"
                    ).encode()
                    verification = {
                        "gate": "PLATFORM_NATIVE_NCBI_REFERENCE",
                        "invocation_identity": (
                            platform_cert.expected_native_invocation_evidence(
                                authority,
                                platform_id,
                                expected["program"],
                                expected["case_id"],
                            )
                        ),
                        "native_fingerprint": (
                            platform_cert.expected_native_fingerprint_evidence(
                                authority,
                                platform_id,
                                expected["program"],
                                expected["case_id"],
                            )
                        ),
                        "native_vs_losat_diagnostic": {
                            "runtime_enforced": False,
                            "classification": "DIAGNOSTIC_ERROR",
                        },
                    }
                else:
                    payload = (
                        f"matrix:{expected['program']}:{expected['case_id']}\n"
                    ).encode()
                    verification = {
                        "gate": (
                            "LOSAT_CANONICAL"
                            if expected["kind"] == "matrix"
                            else "LOSAT_REPEATABILITY"
                        ),
                        "classification": "CANONICAL_PR5_RAW_BYTES",
                        "canonical_sha256": expected["expected_losat_sha256"],
                    }
                output_rel = expected["output_rel"]
                output_path = artifact / output_rel
                output_path.parent.mkdir(parents=True, exist_ok=True)
                output_path.write_bytes(payload)
                digest = hashlib.sha256(payload).hexdigest()
                platform_cert.atomic_write_json(
                    artifact / "completions" / f"{step_id.replace(':', '__')}.json",
                    {
                        "status": "VERIFIED",
                        "source_sha": "a" * 40,
                        "platform_id": platform_id,
                        "authority_version": authority.authority_version,
                        "authority_file_sha256": authority.file_sha256,
                        "losat_binary_sha256": losat_sha256,
                        **expected,
                        "output_sha256": digest,
                        "verification": verification,
                    },
                )
            for required in (
                "cross_platform_equalities.tsv",
                "direct_ncbi_checks.tsv",
                "repeatability.tsv",
                "fixture-staging.json",
                "controlled-path-preflight.json",
            ):
                platform_cert.atomic_write_bytes(artifact / required, b"unit fixture\n")
            rows = []
            for path in sorted(artifact.rglob("*")):
                if path.is_file() and path.name != "FINAL_EVIDENCE_MANIFEST.sha256":
                    rows.append(
                        f"{platform_cert.sha256_path(path)}  "
                        f"{path.relative_to(artifact).as_posix()}\n"
                    )
            platform_cert.atomic_write_bytes(
                artifact / "FINAL_EVIDENCE_MANIFEST.sha256",
                "".join(rows).encode(),
            )
        return authority, catalog

    def refresh_evidence_manifest(self, artifact: Path) -> None:
        rows = []
        for path in sorted(artifact.rglob("*")):
            if path.is_file() and path.name != "FINAL_EVIDENCE_MANIFEST.sha256":
                rows.append(
                    f"{platform_cert.sha256_path(path)}  "
                    f"{path.relative_to(artifact).as_posix()}\n"
                )
        platform_cert.atomic_write_bytes(
            artifact / "FINAL_EVIDENCE_MANIFEST.sha256",
            "".join(rows).encode(),
        )

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

    def test_native_authority_has_exact_bounded_cardinalities(self) -> None:
        self.assertEqual(
            self.authority.authority_version, platform_cert.NATIVE_AUTHORITY_VERSION
        )
        self.assertEqual(len(self.authority.representatives), 6)
        self.assertEqual(len(self.authority.platforms), 3)
        self.assertEqual(len(self.authority.outputs), 18)
        self.assertEqual(len(self.authority.diagnostics), 18)
        provenance = self.authority.document["characterization_provenance"]
        self.assertEqual(len(provenance["rich_groups"]), 11)
        self.assertEqual(
            provenance["rich_group_counts"],
            {"exhaustive": 11, "authoritative_seeded": 7, "rich_only": 4},
        )
        self.assertEqual(
            self.authority.file_sha256,
            hashlib.sha256(
                platform_cert.read_git_blob_bytes(
                    self.repo_root,
                    self.candidate_sha,
                    platform_cert.NATIVE_AUTHORITY_REPO_RELATIVE,
                )
            ).hexdigest(),
        )

    def test_strict_json_rejects_duplicate_keys(self) -> None:
        raw = self.authority_path.read_text(encoding="utf-8")
        duplicate = raw.replace(
            '  "schema_version": 1',
            '  "schema_version": 1,\n  "schema_version": 1',
            1,
        )
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "duplicate.json"
            path.write_text(duplicate, encoding="utf-8")
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure, "duplicate authority JSON key"
            ):
                platform_cert._load_native_authority_bytes(path.read_bytes(), path)

    def test_runtime_authority_hash_uses_git_blob_when_worktree_is_crlf(self) -> None:
        raw = platform_cert.read_git_blob_bytes(
            self.repo_root,
            self.candidate_sha,
            platform_cert.NATIVE_AUTHORITY_REPO_RELATIVE,
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            relative = platform_cert.NATIVE_AUTHORITY_REPO_RELATIVE.as_posix()
            candidate_sha = self.commit_files(root, {relative: raw})
            authority_path = root / relative
            authority_path.write_bytes(raw.replace(b"\n", b"\r\n"))
            authority = platform_cert.load_native_authority(
                root, candidate_sha, authority_path
            )
        self.assertEqual(
            authority.file_sha256,
            "fb156568b854e64196cb0bcc369cd32fd1b1d298f654c5d230ce64ed62f10c2f",
        )
        self.assertEqual(authority.document, self.authority.document)

    def test_git_blob_authority_retains_strict_duplicate_key_validation(self) -> None:
        raw = platform_cert.read_git_blob_bytes(
            self.repo_root,
            self.candidate_sha,
            platform_cert.NATIVE_AUTHORITY_REPO_RELATIVE,
        )
        duplicate = raw.replace(
            b'  "schema_version": 1',
            b'  "schema_version": 1,\n  "schema_version": 1',
            1,
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            relative = platform_cert.NATIVE_AUTHORITY_REPO_RELATIVE.as_posix()
            candidate_sha = self.commit_files(root, {relative: duplicate})
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure, "duplicate authority JSON key"
            ):
                platform_cert.load_native_authority(
                    root, candidate_sha, root / relative
                )

    def test_git_blob_reader_and_candidate_identity_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate_sha = self.commit_files(root, {"tracked.txt": b"tracked\n"})
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure, "candidate Git blob is missing"
            ):
                platform_cert.read_git_blob_bytes(root, candidate_sha, "missing.txt")
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure,
                "candidate Git object is missing or is not a commit",
            ):
                platform_cert.read_git_blob_bytes(root, "0" * 40, "tracked.txt")
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure, "certification SHA mismatch"
            ):
                platform_cert.validate_git_identity(root, "0" * 40)

    def test_canonical_manifest_is_parsed_and_identified_from_git_bytes(self) -> None:
        raw = platform_cert.read_git_blob_bytes(
            self.repo_root,
            self.candidate_sha,
            platform_cert.CANONICAL_MANIFEST_REPO_RELATIVE,
        )
        self.assertEqual(
            hashlib.sha256(raw).hexdigest(),
            "7e65cfb4de0dc575c47d059f1e29ddc9999fc50887ac769fb8f299a9ec97760c",
        )
        self.assertEqual(
            hashlib.sha256(raw.replace(b"\n", b"\r\n")).hexdigest(),
            "199e123e7980b01a68c6ad9605b7c7dec7310e06e7d339b476def52b52c06e4e",
        )
        with mock.patch.object(
            platform_cert,
            "read_git_blob_bytes",
            wraps=platform_cert.read_git_blob_bytes,
        ) as read_blob:
            catalog = platform_cert.load_catalog(self.repo_root, self.candidate_sha)
        read_blob.assert_any_call(
            self.repo_root,
            self.candidate_sha,
            platform_cert.CANONICAL_MANIFEST_REPO_RELATIVE,
        )
        self.assertEqual(catalog.canonical_rows, self.catalog.canonical_rows)

    def test_strict_authority_rejects_unknown_fields_versions_and_platforms(
        self,
    ) -> None:
        mutations = {
            "unknown field": lambda d: d.update({"unexpected": True}),
            "authority version": lambda d: d.update({"authority_version": "v-next"}),
            "schema version": lambda d: d.update({"schema_version": 2}),
            "unknown platform": lambda d: d["normative_acceptance"]["platforms"][
                0
            ].update({"platform_id": "linux-nearest"}),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                self.assert_authority_rejected(mutate)

    def test_strict_authority_rejects_cardinality_duplicates_and_missing_tuples(
        self,
    ) -> None:
        def duplicate_output(document) -> None:
            outputs = document["normative_acceptance"]["platforms"][0][
                "expected_outputs"
            ]
            outputs.append(copy.deepcopy(outputs[0]))

        mutations = {
            "five representatives": lambda d: d["normative_acceptance"][
                "representatives"
            ].pop(),
            "two platforms": lambda d: d["normative_acceptance"]["platforms"].pop(),
            "missing output": lambda d: d["normative_acceptance"]["platforms"][0][
                "expected_outputs"
            ].pop(),
            "duplicate output": duplicate_output,
            "missing diagnostic": lambda d: d["diagnostic_metadata"].pop(),
            "ten rich groups": lambda d: d["characterization_provenance"][
                "rich_groups"
            ].pop(),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                self.assert_authority_rejected(mutate)

    def test_strict_authority_rejects_normative_identity_mutations(self) -> None:
        def mutate_nested(path, value):
            def apply(document) -> None:
                target = document
                for item in path[:-1]:
                    target = target[item]
                target[path[-1]] = value

            return apply

        mutations = {
            "archive filename": mutate_nested(
                ("normative_acceptance", "platforms", 0, "archive", "filename"),
                "other.tar.gz",
            ),
            "checksum algorithm": mutate_nested(
                (
                    "normative_acceptance",
                    "platforms",
                    0,
                    "archive",
                    "checksum",
                    "algorithm",
                ),
                "sha1",
            ),
            "checksum value": mutate_nested(
                (
                    "normative_acceptance",
                    "platforms",
                    0,
                    "archive",
                    "checksum",
                    "value",
                ),
                "0" * 32,
            ),
            "query length": mutate_nested(
                (
                    "normative_acceptance",
                    "representatives",
                    0,
                    "query",
                    "byte_length",
                ),
                0,
            ),
            "query newline": mutate_nested(
                (
                    "normative_acceptance",
                    "representatives",
                    0,
                    "query",
                    "newline_profile",
                    "bare_cr_count",
                ),
                1,
            ),
            "subject length": mutate_nested(
                (
                    "normative_acceptance",
                    "representatives",
                    0,
                    "subject",
                    "byte_length",
                ),
                0,
            ),
            "subject newline": mutate_nested(
                (
                    "normative_acceptance",
                    "representatives",
                    0,
                    "subject",
                    "newline_profile",
                    "bare_cr_count",
                ),
                1,
            ),
            "outfmt": mutate_nested(
                (
                    "normative_acceptance",
                    "representatives",
                    0,
                    "authoritative_outfmt",
                ),
                "6",
            ),
            "environment": mutate_nested(
                (
                    "normative_acceptance",
                    "representatives",
                    0,
                    "environment_overrides",
                ),
                {"RAYON_NUM_THREADS": "1"},
            ),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                self.assert_authority_rejected(mutate)

    def test_platform_sanity_is_exact_after_narrow_machine_normalization(self) -> None:
        self.assertEqual(
            platform_cert.validate_native_platform_sanity(
                self.authority,
                "windows-x64",
                observed_system="Windows",
                observed_machine="AMD64",
            )["platform_id"],
            "windows-x64",
        )
        self.assertEqual(
            platform_cert.validate_native_platform_sanity(
                self.authority,
                "macos-arm64",
                observed_system="Darwin",
                observed_machine="aarch64",
            )["platform_id"],
            "macos-arm64",
        )
        for platform_id, system, machine in (
            ("windows-x64", "Darwin", "x86_64"),
            ("macos-arm64", "Darwin", "x86_64"),
            ("nearest-platform", "Darwin", "x86_64"),
        ):
            with self.subTest(platform_id=platform_id):
                with self.assertRaisesRegex(
                    platform_cert.CertificationFailure,
                    "NATIVE_NCBI_AUTHORITY_MISS",
                ):
                    platform_cert.validate_native_platform_sanity(
                        self.authority,
                        platform_id,
                        observed_system=system,
                        observed_machine=machine,
                    )

    def test_native_invocation_rejects_archive_executable_input_command_and_environment_mutations(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary)
            oracles = {
                "blastn": output / "blastn",
                "blastp": output / "blastp",
                "tblastx": output / "tblastx",
            }
            steps = platform_cert.build_steps(
                self.repo_root,
                output,
                self.catalog,
                output / "LOSAT",
                oracles,
            )
            step = next(
                item
                for item in steps
                if item.step_id == "oracle:blastn:PesePMNV.MjPMNV.task_blastn"
            )
            platform_record = self.authority.platforms["macos-x64"]
            archive = platform_record["archive"]
            checksum = archive["checksum"]
            executable = next(
                item
                for item in platform_record["executables"]
                if item["program"] == "blastn"
            )
            identity = {
                "source_sha": self.candidate_sha,
                "platform_id": "macos-x64",
                "authority_version": self.authority.authority_version,
                "authority_file_sha256": self.authority.file_sha256,
                "oracle_archive": {
                    "archive_name": archive["filename"],
                    "archive_checksum_algorithm": checksum["algorithm"],
                    "archive_md5": checksum["value"],
                },
                "oracles": {
                    "blastn": {
                        "path": step.command[0],
                        "sha256": executable["sha256"],
                        "architecture": executable["architecture"],
                    }
                },
            }
            representative = self.authority.representatives[
                ("blastn", "PesePMNV.MjPMNV.task_blastn")
            ]
            controlled = {
                representative[role]["lexical_path"]: self.repo_root
                / representative[role]["repository_relative"]
                for role in ("query", "subject")
            }

            def verify(authority, observed_identity, observed_step) -> None:
                with mock.patch.object(
                    platform_cert,
                    "controlled_physical_path",
                    side_effect=lambda lexical: controlled[lexical],
                ):
                    platform_cert.verify_native_invocation_identity(
                        self.repo_root,
                        output,
                        authority,
                        observed_identity,
                        observed_step,
                    )

            verify(self.authority, identity, step)
            mutations = []
            changed_identity = copy.deepcopy(identity)
            changed_identity["oracle_archive"]["archive_name"] = "other.tar.gz"
            mutations.append(("archive", self.authority, changed_identity, step))
            changed_identity = copy.deepcopy(identity)
            changed_identity["oracles"]["blastn"]["sha256"] = "0" * 64
            mutations.append(("executable", self.authority, changed_identity, step))
            for role, field, value in (
                ("query", "sha256", "0" * 64),
                ("query", "byte_length", 1),
                (
                    "query",
                    "newline_profile",
                    {
                        "lf_count": 1,
                        "crlf_count": 1,
                        "bare_lf_count": 0,
                        "bare_cr_count": 0,
                    },
                ),
                ("subject", "sha256", "0" * 64),
                ("subject", "byte_length", 1),
                (
                    "subject",
                    "newline_profile",
                    {
                        "lf_count": 1,
                        "crlf_count": 1,
                        "bare_lf_count": 0,
                        "bare_cr_count": 0,
                    },
                ),
            ):
                representatives = copy.deepcopy(self.authority.representatives)
                representatives[("blastn", "PesePMNV.MjPMNV.task_blastn")][role][
                    field
                ] = value
                mutations.append(
                    (
                        f"{role} {field}",
                        replace(self.authority, representatives=representatives),
                        identity,
                        step,
                    )
                )
            command = list(step.command)
            command[command.index("-task")] = "-not-task"
            mutations.append(
                (
                    "command",
                    self.authority,
                    identity,
                    replace(step, command=tuple(command)),
                )
            )
            command = list(step.command)
            command[command.index("-outfmt") + 1] = "6"
            mutations.append(
                (
                    "outfmt",
                    self.authority,
                    identity,
                    replace(step, command=tuple(command)),
                )
            )
            command = list(step.command)
            command[command.index("-out") + 1] = str(output / "other.out")
            mutations.append(
                (
                    "output path",
                    self.authority,
                    identity,
                    replace(step, command=tuple(command)),
                )
            )
            mutations.append(
                (
                    "environment",
                    self.authority,
                    identity,
                    replace(step, environment=(("RAYON_NUM_THREADS", "1"),)),
                )
            )
            for label, authority, observed_identity, observed_step in mutations:
                with self.subTest(label=label):
                    with self.assertRaisesRegex(
                        platform_cert.CertificationFailure,
                        "NATIVE_NCBI_AUTHORITY_MISS",
                    ):
                        verify(authority, observed_identity, observed_step)

    def test_changed_committed_gate_b_fixture_hard_fails_without_normalization(
        self,
    ) -> None:
        representative = self.authority.representatives[
            ("blastn", "PesePMNV.MjPMNV.task_blastn")
        ]
        files = {}
        for role in ("query", "subject"):
            relative = representative[role]["repository_relative"]
            files[relative] = platform_cert.read_git_blob_bytes(
                self.repo_root, self.candidate_sha, relative
            )
        files[representative["query"]["repository_relative"]] += b"CHANGED\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            changed_sha = self.commit_files(root, files)
            output = root / "evidence"
            output.mkdir()
            steps = platform_cert.build_steps(
                self.repo_root,
                output,
                self.catalog,
                root / "LOSAT",
                {
                    "blastn": root / "blastn",
                    "blastp": root / "blastp",
                    "tblastx": root / "tblastx",
                },
            )
            step = next(
                item
                for item in steps
                if item.step_id == "oracle:blastn:PesePMNV.MjPMNV.task_blastn"
            )
            identity = self.native_invocation_identity(step, "macos-x64")
            identity["source_sha"] = changed_sha
            controlled = {}
            for role in ("query", "subject"):
                expected = representative[role]
                destination = root / "controlled" / f"{role}.fasta"
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_bytes(
                    platform_cert.read_git_blob_bytes(
                        root, changed_sha, expected["repository_relative"]
                    )
                )
                controlled[expected["lexical_path"]] = destination
            with (
                mock.patch.object(
                    platform_cert,
                    "controlled_physical_path",
                    side_effect=lambda lexical: controlled[lexical],
                ),
                self.assertRaisesRegex(
                    platform_cert.CertificationFailure,
                    "native invocation query input identity changed",
                ),
            ):
                platform_cert.verify_native_invocation_identity(
                    root,
                    output,
                    self.authority,
                    identity,
                    step,
                )

    def test_native_raw_fingerprint_is_decisive_for_every_field(self) -> None:
        payload = b"q\ts\t100.000\t4\t0\t0\t1\t4\t1\t4\t1e-10\t20.0\r\n"
        key = ("windows-x64", "blastn", "PesePMNV.MjPMNV.task_blastn")
        authority = self.authority_for_payload(*key, payload)
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "native.out"
            output.write_bytes(payload)
            result = platform_cert.verify_native_raw_fingerprint(
                authority, *key, output
            )
            self.assertEqual(result["classification"], "NATIVE_NCBI_FINGERPRINT_EXACT")
            mutations = {
                "raw hash": {"raw_sha256": "0" * 64},
                "byte length": {"byte_length": len(payload) + 1},
                "newline": {
                    "newline_profile": platform_cert.newline_profile(
                        payload.replace(b"\r\n", b"\n")
                    )
                },
                "row count": {"data_row_count": 2},
            }
            for label, change in mutations.items():
                with self.subTest(label=label):
                    rejected = self.authority_for_payload(*key, payload, **change)
                    with self.assertRaisesRegex(
                        platform_cert.CertificationFailure,
                        "NATIVE_NCBI_AUTHORITY_MISS",
                    ):
                        platform_cert.verify_native_raw_fingerprint(
                            rejected, *key, output
                        )

    def test_windows_lf_equal_but_raw_different_still_fails(self) -> None:
        expected = b"q\ts\r\n"
        observed = b"q\ts\n"
        key = ("windows-x64", "blastp", "pairwise_default_serial")
        authority = self.authority_for_payload(*key, expected)
        self.assertEqual(
            platform_cert.canonical_newlines(expected),
            platform_cert.canonical_newlines(observed),
        )
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "native.out"
            output.write_bytes(observed)
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure, "NATIVE_NCBI_AUTHORITY_MISS"
            ):
                platform_cert.verify_native_raw_fingerprint(authority, *key, output)

    def test_tblastx_reordering_and_semantic_equivalence_cannot_rescue_gate_b(
        self,
    ) -> None:
        cases = (
            (
                ("macos-x64", "tblastx", "p03_mela_pemojnva"),
                b"q1\ts1\t1\nq2\ts2\t2\n",
                b"q2\ts2\t2\nq1\ts1\t1\n",
            ),
            (
                ("macos-x64", "blastp", "pairwise_default_serial"),
                b"q\ts\t100.000\t4\t0\t0\t1\t4\t1\t4\t1e-10\t20.0\n",
                b"q\ts\t100.000\t4\t0\t0\t1\t4\t1\t4\t1e-10\t21.0\n",
            ),
        )
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "native.out"
            for key, expected, observed in cases:
                with self.subTest(key=key):
                    authority = self.authority_for_payload(*key, expected)
                    before = copy.deepcopy(authority.outputs)
                    output.write_bytes(observed)
                    for _ in range(2):
                        with self.assertRaisesRegex(
                            platform_cert.CertificationFailure,
                            "NATIVE_NCBI_AUTHORITY_MISS",
                        ):
                            platform_cert.verify_native_raw_fingerprint(
                                authority, *key, output
                            )
                    self.assertEqual(authority.outputs, before)

    def test_selector_has_no_linux_nearest_or_case_fallback(self) -> None:
        for selector in (
            ("linux-x64", "blastp", "pairwise_default_serial"),
            ("macos-arm64", "blastp", "nearest-case"),
        ):
            with self.subTest(selector=selector):
                with self.assertRaisesRegex(
                    platform_cert.CertificationFailure,
                    "NATIVE_NCBI_AUTHORITY_MISS",
                ):
                    platform_cert.resolve_native_authority(self.authority, *selector)

    def test_losat_and_native_axes_remain_orthogonal(self) -> None:
        platform_cert.validate_native_authority_catalog(self.authority, self.catalog)
        sakai = self.authority.representatives[("blastn", "Sakai.MG1655.megablast")]
        d06 = self.authority.representatives[("tblastx", "d06_ap027131_ap027133_db4")]
        self.assertEqual(sakai["losat_parity_class"], "SOURCE_UNDETERMINED_ACCEPTED")
        self.assertEqual(d06["losat_contract"], "approved_db_gencode_deviation")
        self.assertEqual(d06["losat_parity_class"], "HSP_SET_DIFF")
        self.assertEqual(
            sum(
                row["contract"] == "approved_db_gencode_deviation"
                for row in self.catalog.canonical_rows
                if row["program"] == "tblastx"
            ),
            6,
        )
        fingerprint = self.authority.outputs[
            ("windows-x64", "blastn", "Sakai.MG1655.megablast")
        ]
        diagnostic = self.authority.diagnostics[
            ("windows-x64", "blastn", "Sakai.MG1655.megablast")
        ]
        self.assertNotIn("losat_parity_class", fingerprint)
        self.assertFalse(diagnostic["runtime_enforced"])

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
            {
                (program, case_id)
                for program, case_id, _ in platform_cert.REPRESENTATIVES
            },
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
        self.assertEqual([step.execution_index for step in steps], list(range(1, 62)))
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
        self.assertTrue(
            all(step.environment == () for step in steps if step.kind == "oracle")
        )
        self.assertTrue(
            all(
                step.environment == (("RAYON_NUM_THREADS", "1"),)
                for step in steps
                if step.kind in {"matrix", "repeatability"}
                and step.program == "tblastx"
            )
        )
        self.assertTrue(
            all(
                step.environment == ()
                for step in steps
                if step.kind in {"matrix", "repeatability"}
                and step.program != "tblastx"
            )
        )
        for step in (item for item in steps if item.kind == "oracle"):
            observed = list(step.command)
            observed[0] = "{oracle}"
            observed[observed.index("-out") + 1] = "{output}"
            expected = self.authority.representatives[(step.program, step.case_id)][
                "ordered_argv_template"
            ]
            self.assertEqual(observed, expected)

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

    def test_fixture_staging_uses_lf_git_blob_not_crlf_worktree_bytes(
        self,
    ) -> None:
        relative = Path("tests/fasta/platform_cert_stage_unit_fixture.fa")
        lexical = platform_cert.historical_lexical_path(relative.as_posix())
        git_payload = b">fixture\nACGTN\n"
        worktree_payload = git_payload.replace(b"\n", b"\r\n")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "LOSAT" / relative
            candidate_sha = self.commit_files(
                root, {f"LOSAT/{relative.as_posix()}": git_payload}
            )
            source.write_bytes(worktree_payload)
            output = root / "evidence"
            output.mkdir()
            physical = root / "controlled" / relative
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
            with mock.patch.object(
                platform_cert,
                "controlled_physical_path",
                return_value=physical,
            ):
                records = platform_cert.stage_required_fixtures(
                    root, candidate_sha, [step], output
                )
                self.assertEqual(len(records), 1)
                self.assertEqual(physical.read_bytes(), git_payload)
                self.assertEqual(records[0].lexical_target_path, lexical)
                self.assertEqual(records[0].physical_target_path, str(physical))
                self.assertEqual(
                    records[0].source_git_object,
                    f"{candidate_sha}:LOSAT/{relative.as_posix()}",
                )
                self.assertEqual(
                    records[0].sha256, hashlib.sha256(git_payload).hexdigest()
                )
                self.assertEqual(records[0].byte_length, len(git_payload))
                self.assertEqual(
                    records[0].newline_profile,
                    platform_cert.newline_profile(git_payload),
                )
                evidence = json.loads(
                    (output / "fixture-staging.json").read_text(encoding="utf-8")
                )
                self.assertEqual(evidence["status"], "VERIFIED")
                self.assertEqual(evidence["staged_input_count"], 1)
                self.assertEqual(
                    evidence["inputs"][0]["source_git_object"],
                    records[0].source_git_object,
                )

    def test_all_controlled_fixtures_and_pese_gate_b_use_candidate_git_bytes(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "evidence"
            output.mkdir()
            oracles = {
                "blastn": root / "blastn",
                "blastp": root / "blastp",
                "tblastx": root / "tblastx",
            }
            steps = platform_cert.build_steps(
                self.repo_root,
                output,
                self.catalog,
                root / "LOSAT",
                oracles,
            )
            physical_root = root / "controlled"

            def physical_for(lexical: str) -> Path:
                relative = platform_cert._fixture_relative_from_lexical(lexical)
                return physical_root.joinpath(*relative.parts)

            with mock.patch.object(
                platform_cert,
                "controlled_physical_path",
                side_effect=physical_for,
            ):
                records = platform_cert.stage_required_fixtures(
                    self.repo_root,
                    self.candidate_sha,
                    steps,
                    output,
                )
                self.assertEqual(
                    len(records), platform_cert.EXPECTED_STAGED_FIXTURE_COUNT
                )
                for record in records:
                    blob = platform_cert.read_git_blob_bytes(
                        self.repo_root,
                        self.candidate_sha,
                        record.source_repo_relative,
                    )
                    self.assertEqual(
                        Path(record.physical_target_path).read_bytes(), blob
                    )
                    self.assertEqual(record.sha256, hashlib.sha256(blob).hexdigest())

                representative = self.authority.representatives[
                    ("blastn", "PesePMNV.MjPMNV.task_blastn")
                ]
                by_relative = {
                    record.source_repo_relative: record for record in records
                }
                for role in ("query", "subject"):
                    expected = representative[role]
                    record = by_relative[expected["repository_relative"]]
                    self.assertEqual(record.sha256, expected["sha256"])
                    self.assertEqual(record.byte_length, expected["byte_length"])
                    self.assertEqual(
                        record.newline_profile, expected["newline_profile"]
                    )

                step = next(
                    item
                    for item in steps
                    if item.step_id == "oracle:blastn:PesePMNV.MjPMNV.task_blastn"
                )
                platform_record = self.authority.platforms["macos-x64"]
                archive = platform_record["archive"]
                checksum = archive["checksum"]
                executable = next(
                    item
                    for item in platform_record["executables"]
                    if item["program"] == "blastn"
                )
                identity = {
                    "source_sha": self.candidate_sha,
                    "platform_id": "macos-x64",
                    "authority_version": self.authority.authority_version,
                    "authority_file_sha256": self.authority.file_sha256,
                    "oracle_archive": {
                        "archive_name": archive["filename"],
                        "archive_checksum_algorithm": checksum["algorithm"],
                        "archive_md5": checksum["value"],
                    },
                    "oracles": {
                        "blastn": {
                            "path": step.command[0],
                            "sha256": executable["sha256"],
                            "architecture": executable["architecture"],
                        }
                    },
                }
                invocation = platform_cert.verify_native_invocation_identity(
                    self.repo_root,
                    output,
                    self.authority,
                    identity,
                    step,
                )
                self.assertEqual(
                    invocation["inputs"]["query"]["sha256"],
                    "d78996c7897146934aee27db0df03c5655b67e75d4216e2a0b2dcf1c94f27093",
                )
                self.assertEqual(
                    invocation["inputs"]["subject"]["sha256"],
                    "0b9f7dfc0bc4720aea9c196d9d6b5fc6addc4058cb9c64ebba2ea6657f6e8fb6",
                )

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

    def test_gate_b_fingerprint_cannot_compensate_for_gate_a_failure(self) -> None:
        payload = b"native-fingerprint-matches-this\n"
        gate_b_authority = self.authority_for_payload(
            "macos-x64", "blastp", "pairwise_default_serial", payload
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output_rel = "matrix/blastp/pairwise_default_serial/losat.out.partial"
            partial = root / output_rel
            writer = (
                "from pathlib import Path;import sys;"
                "Path(sys.argv[1]).parent.mkdir(parents=True,exist_ok=True);"
                "Path(sys.argv[1]).write_bytes(bytes.fromhex(sys.argv[2]))"
            )
            step = platform_cert.SearchStep(
                execution_index=1,
                step_id="matrix:blastp:pairwise_default_serial",
                kind="matrix",
                program="blastp",
                case_id="pairwise_default_serial",
                semantic_class="source_defined_exact",
                run_index=1,
                command=(sys.executable, "-c", writer, str(partial), payload.hex()),
                environment=(),
                output_rel=output_rel,
                final_output_rel=output_rel.removesuffix(".partial"),
                expected_losat_sha256="0" * 64,
            )
            identity = {
                "source_sha": "a" * 40,
                "platform_id": "macos-x64",
                "losat": {"sha256": "b" * 64},
                "authority_version": gate_b_authority.authority_version,
                "authority_file_sha256": gate_b_authority.file_sha256,
            }
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure,
                "CROSS_PLATFORM_LOSAT_DIVERGENCE",
            ):
                platform_cert._execute_step(
                    root,
                    root,
                    identity,
                    self.catalog,
                    step,
                    gate_b_authority,
                )

    def test_gate_a_does_not_consume_mutated_gate_b_fingerprints(self) -> None:
        payload = b"canonical-losat-output\n"
        mutated_authority = self.authority_for_payload(
            "macos-x64",
            "blastp",
            "pairwise_default_serial",
            b"unrelated-native-fingerprint\n",
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output_rel = "matrix/blastp/pairwise_default_serial/losat.out.partial"
            partial = root / output_rel
            writer = (
                "from pathlib import Path;import sys;"
                "Path(sys.argv[1]).parent.mkdir(parents=True,exist_ok=True);"
                "Path(sys.argv[1]).write_bytes(bytes.fromhex(sys.argv[2]))"
            )
            step = platform_cert.SearchStep(
                execution_index=1,
                step_id="matrix:blastp:pairwise_default_serial",
                kind="matrix",
                program="blastp",
                case_id="pairwise_default_serial",
                semantic_class="source_defined_exact",
                run_index=1,
                command=(sys.executable, "-c", writer, str(partial), payload.hex()),
                environment=(),
                output_rel=output_rel,
                final_output_rel=output_rel.removesuffix(".partial"),
                expected_losat_sha256=hashlib.sha256(payload).hexdigest(),
            )
            identity = {
                "source_sha": "a" * 40,
                "platform_id": "macos-x64",
                "losat": {"sha256": "b" * 64},
                "authority_version": mutated_authority.authority_version,
                "authority_file_sha256": mutated_authority.file_sha256,
            }
            record = platform_cert._execute_step(
                root, root, identity, self.catalog, step, mutated_authority
            )
            self.assertEqual(record["output_sha256"], step.expected_losat_sha256)
            self.assertEqual(record["verification"]["gate"], "LOSAT_CANONICAL")

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
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure,
                "NATIVE_NCBI_AUTHORITY_MISS",
            ):
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
            checksum.write_text(f"{digest}  {spec.archive_name}\n", encoding="ascii")
            record = platform_cert.verify_archive(spec, archive, checksum)
            self.assertEqual(record["archive_md5"], digest)
            checksum.write_text(f"{'0' * 32}  {spec.archive_name}\n", encoding="ascii")
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure,
                "NATIVE_NCBI_AUTHORITY_MISS",
            ):
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

    def test_blastn_semantic_comparison_is_diagnostic_only(self) -> None:
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
            result = platform_cert.report_native_vs_losat_diagnostic(
                self.catalog, step, ncbi, losat, root / "diagnostic"
            )
            self.assertEqual(result["classification"], "EXACT_DATA")
            self.assertEqual(result["platform_diagnostic"], "LINE_ENDINGS_ONLY")
            self.assertFalse(result["runtime_enforced"])

            ncbi.write_bytes((row.replace("20.0", "21.0") + "\r\n").encode())
            result = platform_cert.report_native_vs_losat_diagnostic(
                self.catalog, step, ncbi, losat, root / "value-difference"
            )
            self.assertEqual(result["classification"], "VALUE_DIFF")
            self.assertFalse(result["runtime_enforced"])

    def test_exact_pese_fingerprint_passes_before_value_difference_diagnostic(
        self,
    ) -> None:
        losat_payload = b"q\ts\t100.000\t4\t0\t0\t1\t4\t1\t4\t1e-10\t20.0\n"
        native_payload = b"q\ts\t100.000\t4\t0\t0\t1\t4\t1\t4\t1e-10\t21.0\n"
        authority = self.authority_for_payload(
            "macos-x64",
            "blastn",
            "PesePMNV.MjPMNV.task_blastn",
            native_payload,
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            losat = root / "matrix/blastn/PesePMNV.MjPMNV.task_blastn/losat.out"
            losat.parent.mkdir(parents=True)
            losat.write_bytes(losat_payload)
            output_rel = "oracle/blastn/PesePMNV.MjPMNV.task_blastn/ncbi.out.partial"
            partial = root / output_rel
            writer = (
                "from pathlib import Path;import sys;"
                "Path(sys.argv[1]).parent.mkdir(parents=True,exist_ok=True);"
                "Path(sys.argv[1]).write_bytes(bytes.fromhex(sys.argv[2]))"
            )
            step = platform_cert.SearchStep(
                execution_index=44,
                step_id="oracle:blastn:PesePMNV.MjPMNV.task_blastn",
                kind="oracle",
                program="blastn",
                case_id="PesePMNV.MjPMNV.task_blastn",
                semantic_class="ordinary_exact",
                run_index=1,
                command=(
                    sys.executable,
                    "-c",
                    writer,
                    str(partial),
                    native_payload.hex(),
                ),
                environment=(),
                output_rel=output_rel,
                final_output_rel=output_rel.removesuffix(".partial"),
                expected_losat_sha256=hashlib.sha256(losat_payload).hexdigest(),
            )
            identity = {
                "source_sha": "a" * 40,
                "platform_id": "macos-x64",
                "authority_version": authority.authority_version,
                "authority_file_sha256": authority.file_sha256,
                "losat": {"sha256": "b" * 64},
            }
            invocation = {"unit_test": "exact identity verified first"}
            with mock.patch.object(
                platform_cert,
                "verify_native_invocation_identity",
                return_value=invocation,
            ):
                record = platform_cert._execute_step(
                    root, root, identity, self.catalog, step, authority
                )
            verification = record["verification"]
            self.assertEqual(
                verification["native_fingerprint"]["classification"],
                "NATIVE_NCBI_FINGERPRINT_EXACT",
            )
            self.assertEqual(
                verification["native_vs_losat_diagnostic"]["classification"],
                "VALUE_DIFF",
            )
            self.assertFalse(
                verification["native_vs_losat_diagnostic"]["runtime_enforced"]
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
            "authority_version": self.authority.authority_version,
            "authority_file_sha256": self.authority.file_sha256,
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
        changed = json.loads(json.dumps(identity))
        changed["authority_file_sha256"] = "3" * 64
        self.assertNotEqual(key, platform_cert._identity_resume_key(changed))
        changed = json.loads(json.dumps(identity))
        changed["authority_version"] = "other"
        self.assertNotEqual(key, platform_cert._identity_resume_key(changed))

    def test_manual_workflow_has_only_the_three_pr6_platform_jobs(self) -> None:
        workflow = (
            self.repo_root
            / ".github"
            / "workflows"
            / "platform-native-certification.yml"
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
        self.assertIn('--expected-sha "$CANDIDATE_SHA"', workflow)
        self.assertIn("preflight-path", workflow)
        self.assertIn("extract-verified-archive", workflow)
        self.assertNotIn("tar -xzf", workflow)
        self.assertIn("actions/upload-artifact@v4", workflow)
        self.assertIn("if: ${{ always() }}", workflow)
        self.assertIn("cross-platform-aggregation:", workflow)
        self.assertIn("needs: [platform-certification]", workflow)
        self.assertEqual(workflow.count("CROSS_PLATFORM_NATIVE_CERTIFIED"), 0)
        self.assertEqual(
            Path(platform_cert.__file__)
            .read_text(encoding="utf-8")
            .count('print("CROSS_PLATFORM_NATIVE_CERTIFIED")'),
            1,
        )
        self.assertEqual(workflow.count("actions/download-artifact@v4"), 5)

    def test_resume_import_revalidates_verified_output_hash(self) -> None:
        payload = b"canonical output\n"
        output_hash = hashlib.sha256(payload).hexdigest()
        identity = {
            "evidence_schema": 1,
            "authority_version": self.authority.authority_version,
            "authority_file_sha256": self.authority.file_sha256,
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
            "authority_version": identity["authority_version"],
            "authority_file_sha256": identity["authority_file_sha256"],
            "output_rel": step.final_output_rel,
            "output_sha256": output_hash,
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prior = root / "prior"
            current = root / "current"
            prior.mkdir()
            current.mkdir()
            (prior / "identity.json").write_text(json.dumps(identity), encoding="utf-8")
            prior_output = prior / step.final_output_rel
            prior_output.parent.mkdir(parents=True)
            prior_output.write_bytes(payload)
            completion = platform_cert._step_record_path(prior, step)
            completion.parent.mkdir(parents=True)
            completion.write_text(json.dumps(record), encoding="utf-8")
            imported = platform_cert.import_resume(
                prior,
                current,
                self.repo_root,
                identity,
                [step],
                self.catalog,
                self.authority,
            )
            self.assertEqual(set(imported), {step.step_id})
            self.assertEqual((current / step.final_output_rel).read_bytes(), payload)

            prior_output.write_bytes(b"corrupt\n")
            second = root / "second"
            second.mkdir()
            with self.assertRaises(platform_cert.CertificationFailure):
                platform_cert.import_resume(
                    prior,
                    second,
                    self.repo_root,
                    identity,
                    [step],
                    self.catalog,
                    self.authority,
                )

    def test_cross_platform_success_token_is_suppressed_on_failed_subgate(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            captured_stdout = io.StringIO()
            captured_stderr = io.StringIO()
            with (
                contextlib.redirect_stdout(captured_stdout),
                contextlib.redirect_stderr(captured_stderr),
            ):
                status = platform_cert.main(
                    [
                        "aggregate",
                        "--input-root",
                        str(root),
                        "--output",
                        str(root / "aggregate.json"),
                        "--expected-sha",
                        "a" * 40,
                        "--authority",
                        str(self.authority_path),
                    ]
                )
            self.assertEqual(status, 1)
            self.assertNotIn(
                "CROSS_PLATFORM_NATIVE_CERTIFIED", captured_stdout.getvalue()
            )
            self.assertIn("PLATFORM_CERTIFICATION_FAILED", captured_stderr.getvalue())

    def test_aggregate_requires_three_complete_authority_bound_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_root = root / "input"
            input_root.mkdir()
            authority, catalog = self.create_aggregate_fixture(input_root)
            output = root / "aggregate.json"
            result = platform_cert.aggregate_platform_evidence(
                input_root, output, "a" * 40, authority, catalog
            )
            self.assertEqual(result["decision"], "CROSS_PLATFORM_NATIVE_CERTIFIED")
            self.assertEqual(result["platform_count"], 3)
            self.assertEqual(result["total_search_executions"], 183)
            self.assertEqual(
                {item["platform_id"] for item in result["platforms"]},
                {"windows-x64", "macos-arm64", "macos-x64"},
            )
            diagnostic_completion = (
                input_root
                / "platform-native-windows-x64"
                / "completions"
                / "oracle__blastn__PesePMNV.MjPMNV.task_blastn.json"
            )
            diagnostic = json.loads(diagnostic_completion.read_text(encoding="utf-8"))[
                "verification"
            ]["native_vs_losat_diagnostic"]
            self.assertEqual(diagnostic["classification"], "DIAGNOSTIC_ERROR")
            self.assertFalse(diagnostic["runtime_enforced"])
            failed = input_root / "platform-native-windows-x64" / "summary.json"
            summary = json.loads(failed.read_text(encoding="utf-8"))
            summary["gate_a_losat_canonical"]["passed"] = 42
            failed.write_text(json.dumps(summary), encoding="utf-8")
            self.refresh_evidence_manifest(failed.parent)
            with self.assertRaises(platform_cert.CertificationFailure):
                platform_cert.aggregate_platform_evidence(
                    input_root,
                    root / "rejected.json",
                    "a" * 40,
                    authority,
                    catalog,
                )

    def test_aggregate_rejects_every_required_evidence_subgate(self) -> None:
        def update_json(path: Path, mutate) -> None:
            value = json.loads(path.read_text(encoding="utf-8"))
            mutate(value)
            platform_cert.atomic_write_json(path, value)

        cases = (
            (
                "topology",
                "aggregate completion identity failed",
                "completions/matrix__blastn__PesePMNV.MjPMNV.task_blastn.json",
                lambda value: value.__setitem__("execution_index", 99),
                True,
            ),
            (
                "gate_a",
                "aggregate Gate A completion failed",
                "completions/matrix__blastn__PesePMNV.MjPMNV.task_blastn.json",
                lambda value: value["verification"].__setitem__(
                    "classification", "NOT_CANONICAL"
                ),
                True,
            ),
            (
                "gate_b_invocation",
                "aggregate Gate B completion failed",
                "completions/oracle__blastn__PesePMNV.MjPMNV.task_blastn.json",
                lambda value: value["verification"]["invocation_identity"].__setitem__(
                    "archive_name", "other.tar.gz"
                ),
                True,
            ),
            (
                "gate_b_fingerprint",
                "aggregate Gate B completion failed",
                "completions/oracle__blastn__PesePMNV.MjPMNV.task_blastn.json",
                lambda value: value["verification"]["native_fingerprint"].__setitem__(
                    "byte_length",
                    value["verification"]["native_fingerprint"]["byte_length"] + 1,
                ),
                True,
            ),
            (
                "repeatability",
                "aggregate repeatability completion failed",
                "completions/repeatability__blastn__PesePMNV.MjPMNV.task_blastn__run2.json",
                lambda value: value["verification"].__setitem__(
                    "classification", "DIFFER"
                ),
                True,
            ),
            (
                "sakai",
                "aggregate Sakai ratchet failed",
                "summary.json",
                lambda value: value["sakai_ratchet"].__setitem__(
                    "observed_losat_parity_class", "PLATFORM_EXCEPTION"
                ),
                True,
            ),
            (
                "d06",
                "aggregate TBLASTX deviation ratchet failed",
                "summary.json",
                lambda value: value["tblastx_deviation_ratchet"].__setitem__(
                    "accepted_deviation_ceiling", 7
                ),
                True,
            ),
            (
                "authority",
                "aggregate platform identity mismatch",
                "identity.json",
                lambda value: value.__setitem__("authority_file_sha256", "0" * 64),
                True,
            ),
            (
                "manifest",
                "final evidence integrity mismatch",
                "direct_ncbi_checks.tsv",
                None,
                False,
            ),
        )
        for label, message, relative, mutate, refresh in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                input_root = root / "input"
                input_root.mkdir()
                authority, catalog = self.create_aggregate_fixture(input_root)
                artifact = input_root / "platform-native-windows-x64"
                target = artifact / relative
                if mutate is None:
                    target.write_bytes(target.read_bytes() + b"tampered\n")
                else:
                    update_json(target, mutate)
                if refresh:
                    self.refresh_evidence_manifest(artifact)
                with self.assertRaisesRegex(
                    platform_cert.CertificationFailure, message
                ):
                    platform_cert.aggregate_platform_evidence(
                        input_root,
                        root / "aggregate.json",
                        "a" * 40,
                        authority,
                        catalog,
                    )


if __name__ == "__main__":
    unittest.main()
