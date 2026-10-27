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
import re
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


SAKAI = ("blastn", "Sakai.MG1655.megablast")
SAKAI_GATE_A_STEPS = [11, 52, 53]
# The aggregate fixture's frozen bytes of the case and a known different output.
SAKAI_FIXTURE_FROZEN = b"matrix:blastn:Sakai.MG1655.megablast\n"
SAKAI_FIXTURE_ALLOWED = b"allowed:blastn:Sakai.MG1655.megablast\n"
SAKAI_FROZEN_PR5 = "3f27c1f1396b59ecb7e78827b5da390e9c97a276c05cbc7a6c82f56bd0a03262"
SAKAI_APPROVED = "ac8177764f35b2fb92777f8b1d01ea646cc55ae21b05f1f486ee231348a34a45"
# The allowlist entry that PD-NCBI-PLATFORM-VARIANCE 1.1 to 1.2 carried; 1.3 moved
# its output into the v0.1.1 manifest and emptied the list.
HISTORICAL_SAKAI_ENTRY = {
    "program": SAKAI[0],
    "case_id": SAKAI[1],
    "frozen_sha256": SAKAI_FROZEN_PR5,
    "allowed_actual_sha256": SAKAI_APPROVED,
    "frozen_gate_a_steps": SAKAI_GATE_A_STEPS,
}


class PlatformCertificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.repo_root = Path(__file__).resolve().parents[2]
        cls.candidate_sha = platform_cert.resolve_git_head_sha(cls.repo_root)
        cls.catalog = platform_cert.load_catalog(cls.repo_root, cls.candidate_sha)
        cls.authority_path = SCRIPT.with_name(
            platform_cert.NATIVE_AUTHORITY_REPO_RELATIVE.name
        )
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

    def fixture_allowlist(self, *extra_entries, **changes):
        entry = {
            "program": SAKAI[0],
            "case_id": SAKAI[1],
            "frozen_sha256": hashlib.sha256(SAKAI_FIXTURE_FROZEN).hexdigest(),
            "allowed_actual_sha256": hashlib.sha256(SAKAI_FIXTURE_ALLOWED).hexdigest(),
            "frozen_gate_a_steps": SAKAI_GATE_A_STEPS,
            **changes,
        }
        document = {
            "schema": "losat-frozen-mismatch-allowlist-v1",
            "entries": [entry, *extra_entries],
        }
        return platform_cert.Allowlist(
            Path("frozen_mismatch_allowlist.json"),
            raw=json.dumps(document).encode(),
        )

    def fixture_allowed_cases(self) -> list[dict[str, object]]:
        return [
            {
                "program": SAKAI[0],
                "case_id": SAKAI[1],
                "frozen_sha256": hashlib.sha256(SAKAI_FIXTURE_FROZEN).hexdigest(),
                "allowed_actual_sha256": hashlib.sha256(
                    SAKAI_FIXTURE_ALLOWED
                ).hexdigest(),
                "frozen_gate_a_steps": SAKAI_GATE_A_STEPS,
            }
        ]

    def create_aggregate_fixture(self, root: Path, allowlist=None):
        allowed_cases = [] if allowlist is None else self.fixture_allowed_cases()
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
            allowlist_identity = (
                {}
                if allowlist is None
                else {
                    "frozen_mismatch_allowlist": {
                        "git_object": (
                            f"{'a' * 40}:LOSAT/tests/frozen_mismatch_allowlist.json"
                        ),
                        "sha256": allowlist.file_sha256,
                        "gate_a_cases": allowed_cases,
                    }
                }
            )
            platform_cert.atomic_write_json(
                artifact / "identity.json",
                {
                    "source_sha": "a" * 40,
                    "platform_id": platform_id,
                    "authority_version": authority.authority_version,
                    "authority_file_sha256": authority.file_sha256,
                    **allowlist_identity,
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
                        "all_exact_canonical_raw_bytes": not allowed_cases,
                        "allowed_known_mismatches": allowed_cases,
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
                        "classification": "CANONICAL_RAW_BYTES",
                        "canonical_sha256": expected["expected_losat_sha256"],
                    }
                    if (
                        allowed_cases
                        and (expected["program"], expected["case_id"]) == SAKAI
                    ):
                        payload = SAKAI_FIXTURE_ALLOWED
                        verification.update(
                            classification="ALLOWED_KNOWN_MISMATCH",
                            allowed_actual_sha256=hashlib.sha256(payload).hexdigest(),
                            allowlist_file_sha256=allowlist.file_sha256,
                        )
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
        self.assertEqual(
            set(self.authority.platforms),
            {"windows-x64", "macos-arm64", "macos-x64", "linux-aarch64"},
        )
        self.assertEqual(set(self.authority.platforms), set(platform_cert.PLATFORMS))
        self.assertEqual(len(self.authority.outputs), 24)
        self.assertEqual(len(self.authority.diagnostics), 24)
        provenance = self.authority.document["characterization_provenance"]
        self.assertEqual(provenance["characterization_run_id"], 33582928402)
        self.assertEqual(len(provenance["rich_groups"]), 11)
        self.assertEqual(
            provenance["rich_group_counts"],
            {"exhaustive": 11, "authoritative_seeded": 7, "rich_only": 4},
        )
        added = self.authority.document["added_platform_characterizations"]
        self.assertEqual(
            [
                (
                    item["characterization_run_id"],
                    item["candidate_sha"],
                    [run["platform_id"] for run in item["platform_runs"]],
                )
                for item in added
            ],
            [
                (
                    *platform_cert.ADDED_PLATFORM_CHARACTERIZATIONS["linux-aarch64"],
                    ["linux-aarch64"],
                )
            ],
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

    def test_authority_v011_carries_v010_and_the_characterized_record(self) -> None:
        def git_bytes(relative: str) -> bytes:
            return platform_cert.read_git_blob_bytes(
                self.repo_root, self.candidate_sha, relative
            )

        v010_raw = git_bytes("LOSAT/tests/ncbi_platform_variance_v010.json")
        # The v0.1.0 file is historical and is never edited.
        self.assertEqual(
            hashlib.sha256(v010_raw).hexdigest(),
            "fb156568b854e64196cb0bcc369cd32fd1b1d298f654c5d230ce64ed62f10c2f",
        )
        v010 = json.loads(v010_raw)
        v011 = self.authority.document
        self.assertEqual(
            v011["normative_acceptance"]["representatives"],
            v010["normative_acceptance"]["representatives"],
        )
        self.assertEqual(
            v011["normative_acceptance"]["platforms"][:3],
            v010["normative_acceptance"]["platforms"],
        )
        self.assertEqual(v011["diagnostic_metadata"][:18], v010["diagnostic_metadata"])
        self.assertEqual(
            v011["characterization_provenance"], v010["characterization_provenance"]
        )
        record_raw = git_bytes(
            "docs/evidence/ncbi_platform_variance_v011/"
            "linux-aarch64-authority-record.json"
        )
        record = json.loads(record_raw)
        self.assertEqual(
            v011["normative_acceptance"]["platforms"][3],
            record["normative_platform_record"],
        )
        self.assertEqual(v011["diagnostic_metadata"][18:], record["diagnostic_metadata"])
        (added,) = v011["added_platform_characterizations"]
        self.assertEqual(
            added["derived_record_sha256"], hashlib.sha256(record_raw).hexdigest()
        )
        self.assertEqual(added["candidate_sha"], record["candidate_sha"])
        self.assertEqual(added["platform_runs"], [record["platform_run"]])
        self.assertEqual(added["execution_counts"], record["execution_counts"])

    def test_strict_authority_binds_one_characterization_source_per_platform(
        self,
    ) -> None:
        def added(document):
            return document["added_platform_characterizations"][0]

        def v010_platform_in_added(document) -> None:
            run = copy.deepcopy(
                document["characterization_provenance"]["platform_runs"][0]
            )
            added(document)["platform_runs"].append(run)

        def added_platform_in_rich_group(document) -> None:
            group = document["characterization_provenance"]["rich_groups"][0]
            multiset = copy.deepcopy(group["platform_representation_multisets"][0])
            multiset["platform_id"] = "linux-aarch64"
            group["platform_representation_multisets"].append(multiset)

        def added_platform_in_v010_runs(document) -> None:
            document["characterization_provenance"]["platform_runs"].append(
                copy.deepcopy(added(document)["platform_runs"][0])
            )

        mutations = {
            "run id": lambda d: added(d).update(
                {"characterization_run_id": 33582928402}
            ),
            "candidate": lambda d: added(d).update({"candidate_sha": "0" * 40}),
            "normative": lambda d: added(d).update({"runtime_enforced": True}),
            "counts": lambda d: added(d)["execution_counts"].update({"losat": 1}),
            "runner label": lambda d: added(d)["platform_runs"][0].update(
                {"runner_label": "ubuntu-24.04"}
            ),
            "unknown field": lambda d: added(d).update({"unexpected": True}),
            "v0.1.0 platform claimed twice": v010_platform_in_added,
            "added platform in rich groups": added_platform_in_rich_group,
            "added platform in v0.1.0 runs": added_platform_in_v010_runs,
            "two added entries": lambda d: d["added_platform_characterizations"].append(
                copy.deepcopy(added(d))
            ),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                self.assert_authority_rejected(mutate)

    def test_strict_json_rejects_duplicate_keys(self) -> None:
        raw = self.authority_path.read_text(encoding="utf-8")
        duplicate = raw.replace(
            '  "schema_version": 2',
            '  "schema_version": 2,\n  "schema_version": 2',
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
            "4f01d4adfa59fd912574098a17a4783fd25ffc6705c5045cf8e99685d004184a",
        )
        self.assertEqual(authority.document, self.authority.document)

    def test_git_blob_authority_retains_strict_duplicate_key_validation(self) -> None:
        raw = platform_cert.read_git_blob_bytes(
            self.repo_root,
            self.candidate_sha,
            platform_cert.NATIVE_AUTHORITY_REPO_RELATIVE,
        )
        duplicate = raw.replace(
            b'  "schema_version": 2',
            b'  "schema_version": 2,\n  "schema_version": 2',
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
            "9bfac81fff33cd770729540494e9c8df1e76b798374d120c26dd1945e5bd1611",
        )
        self.assertEqual(
            hashlib.sha256(raw.replace(b"\n", b"\r\n")).hexdigest(),
            "d3cec365ca3ee9ef5aceeda3bedc3341c6e5b92ed9a2210712f1f4945349fd56",
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
        read_blob.assert_any_call(
            self.repo_root,
            self.candidate_sha,
            platform_cert.CANONICAL_V010_REPO_RELATIVE,
        )
        self.assertEqual(catalog.canonical_rows, self.catalog.canonical_rows)

    def test_canonical_v011_is_v010_with_only_the_approved_sakai_output(self) -> None:
        v010 = platform_cert.read_git_blob_bytes(
            self.repo_root,
            self.candidate_sha,
            platform_cert.CANONICAL_V010_REPO_RELATIVE,
        )
        self.assertEqual(
            hashlib.sha256(v010).hexdigest(), platform_cert.CANONICAL_V010_SHA256
        )
        _, v010_lines = platform_cert._canonical_lines(v010)
        v010_rows = list(platform_cert.csv.DictReader(v010_lines, delimiter="\t"))
        changed = [
            (old, new)
            for old, new in zip(v010_rows, self.catalog.canonical_rows)
            if old != new
        ]
        self.assertEqual(len(v010_rows), 43)
        self.assertEqual(len(self.catalog.canonical_rows), 43)
        self.assertEqual(len(changed), 1)
        old, new = changed[0]
        self.assertEqual((old["program"], old["case_id"]), SAKAI)
        self.assertEqual(old["losat_sha256"], SAKAI_FROZEN_PR5)
        self.assertEqual(new, {**old, "losat_sha256": SAKAI_APPROVED})
        # The contract and class stay those of the authority's representative.
        sakai = self.catalog.canonical[SAKAI]
        self.assertEqual(sakai["contract"], "SOURCE_UNDETERMINED_ACCEPTED")
        self.assertEqual(sakai["classification"], "SOURCE_UNDETERMINED_ACCEPTED")
        # The approved output is the retained Linux NCBI BLAST+ 2.17.0 output.
        retained = {
            entry["retained_linux_raw_sha256"]
            for entry in self.authority.document["diagnostic_metadata"]
            if (entry["program"], entry["case_id"]) == SAKAI
        }
        self.assertEqual(retained, {SAKAI_APPROVED})

    def test_canonical_manifest_must_derive_from_the_v010_file(self) -> None:
        real = platform_cert.read_git_blob_bytes
        canonical = platform_cert.CANONICAL_MANIFEST_REPO_RELATIVE
        v010 = platform_cert.CANONICAL_V010_REPO_RELATIVE
        v011_raw = real(self.repo_root, self.candidate_sha, canonical)
        v010_raw = real(self.repo_root, self.candidate_sha, v010)
        other_row = b"PesePMNV.MjPMNV.task_blastn\tEXACT_TEXT\tEXACT_TEXT\tfcd9236c"
        self.assertEqual(v011_raw.count(other_row), 1)
        cases = {
            "v0.1.0 canonical manifest changed": {
                v010: v010_raw.replace(b"PR 5", b"PR5", 1)
            },
            "canonical manifest is not v0.1.0 with the approved outputs": {
                canonical: v011_raw.replace(other_row, other_row[:-1] + b"d"),
            },
            "canonical manifest provenance comments changed": {
                canonical: v011_raw.replace(b"# DERIVED_FROM_SHA256", b"# DERIVED_SHA256"),
            },
        }
        for message, replaced in cases.items():
            def fake(root, sha, relative, replaced=replaced):
                key = platform_cert.PurePosixPath(str(relative))
                if key in replaced:
                    return replaced[key]
                return real(root, sha, relative)

            with self.subTest(message=message), mock.patch.object(
                platform_cert, "read_git_blob_bytes", side_effect=fake
            ):
                with self.assertRaisesRegex(platform_cert.CertificationFailure, message):
                    platform_cert.load_catalog(self.repo_root, self.candidate_sha)
        # Behind the v0.1.0 file hash: the approved output must replace exactly
        # the frozen PR 5 hash of a case that v0.1.0 lists.
        sakai_row = next(
            line
            for line in v010_raw.splitlines(keepends=True)
            if line.startswith(b"blastn\tSakai.MG1655.megablast\t")
        )
        for message, altered in (
            (
                "approved output replaces another hash",
                v010_raw.replace(SAKAI_FROZEN_PR5.encode(), b"0" * 64),
            ),
            ("approved output case missing", v010_raw.replace(sakai_row, b"")),
        ):
            with self.subTest(message=message), mock.patch.object(
                platform_cert,
                "CANONICAL_V010_SHA256",
                hashlib.sha256(altered).hexdigest(),
            ):
                with self.assertRaisesRegex(platform_cert.CertificationFailure, message):
                    platform_cert.derive_canonical_v011_rows(altered)

    def test_strict_authority_rejects_unknown_fields_versions_and_platforms(
        self,
    ) -> None:
        mutations = {
            "unknown field": lambda d: d.update({"unexpected": True}),
            "authority version": lambda d: d.update({"authority_version": "v-next"}),
            "schema version": lambda d: d.update({"schema_version": 1}),
            "v0.1.0 version": lambda d: d.update(
                {"authority_version": "ncbi-platform-variance-v0.1.0"}
            ),
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
            "no added characterization": lambda d: d[
                "added_platform_characterizations"
            ].clear(),
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
            if platform_cert._option_value(step.command, "-outfmt") == "7"
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
                    platform_cert._option_value(step.command, "-outfmt"), "6"
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
                    "-outfmt",
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
                "CANONICAL_RAW_BYTES",
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

    def test_binary_architecture_reads_pe_macho_and_elf_headers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pe = bytearray(128)
            pe[:2] = b"MZ"
            pe[60:64] = (64).to_bytes(4, "little")
            pe[64:68] = b"PE\0\0"
            pe[68:70] = (0x8664).to_bytes(2, "little")
            pe_path = root / "losat.exe"
            pe_path.write_bytes(pe)
            self.assertEqual(platform_cert.binary_architecture(pe_path), "x86_64")

            macho_path = root / "LOSAT"
            macho_path.write_bytes(
                b"\xcf\xfa\xed\xfe" + (0x0100000C).to_bytes(4, "little") + b"\0" * 64
            )
            self.assertEqual(platform_cert.binary_architecture(macho_path), "arm64")

            # ELF64: e_ident[EI_DATA] at 5, e_machine at 18 (EM_AARCH64 0xB7,
            # EM_X86_64 0x3E), in the byte order EI_DATA names.
            for byte_order, ei_data in (("little", 1), ("big", 2)):
                for machine, expected in ((0xB7, "arm64"), (0x3E, "x86_64")):
                    elf = bytearray(64)
                    elf[:4] = b"\x7fELF"
                    elf[4] = 2
                    elf[5] = ei_data
                    elf[18:20] = machine.to_bytes(2, byte_order)
                    elf_path = root / f"LOSAT-elf-{byte_order}-{machine:x}"
                    elf_path.write_bytes(bytes(elf))
                    self.assertEqual(
                        platform_cert.binary_architecture(elf_path), expected
                    )
            elf = bytearray(64)
            elf[:4] = b"\x7fELF"
            elf[5] = 1
            elf[18:20] = (0xF3).to_bytes(2, "little")
            other_path = root / "LOSAT-riscv"
            other_path.write_bytes(bytes(elf))
            self.assertEqual(platform_cert.binary_architecture(other_path), "unknown")

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
            "frozen_mismatch_allowlist": {"sha256": "4" * 64, "gate_a_cases": []},
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
        changed = json.loads(json.dumps(identity))
        changed["frozen_mismatch_allowlist"]["sha256"] = "5" * 64
        self.assertNotEqual(key, platform_cert._identity_resume_key(changed))
        changed = json.loads(json.dumps(identity))
        del changed["frozen_mismatch_allowlist"]
        with self.assertRaisesRegex(
            platform_cert.CertificationFailure, "frozen mismatch allowlist identity"
        ):
            platform_cert._identity_resume_key(changed)

    def test_manual_workflow_has_only_the_authority_platform_jobs(self) -> None:
        workflow = (
            self.repo_root
            / ".github"
            / "workflows"
            / "platform-native-certification.yml"
        ).read_text(encoding="utf-8")
        self.assertIn("workflow_dispatch:", workflow)
        self.assertNotIn("pull_request:", workflow)
        self.assertNotIn("schedule:", workflow)
        self.assertIn("max-parallel: 4", workflow)
        self.assertEqual(workflow.count("platform_id:"), 4)
        for spec in platform_cert.PLATFORMS.values():
            self.assertIn(f"runner: {spec.runner_label}", workflow)
            self.assertIn(f"target_triple: {spec.target_triple}", workflow)
            self.assertIn(f"oracle_archive: {spec.archive_name}", workflow)
            self.assertIn(f"name: platform-native-{spec.platform_id}", workflow)
        for excluded in ("ubuntu-latest", "wasm32", "benchmark", "soak"):
            self.assertNotIn(excluded, workflow.lower())
        self.assertIn('--expected-sha "$CANDIDATE_SHA"', workflow)
        self.assertIn("preflight-path", workflow)
        self.assertIn("extract-verified-archive", workflow)
        self.assertNotIn("tar -xzf", workflow)
        # Actions are pinned to a commit with the release tag as a comment.
        self.assertRegex(workflow, r"actions/upload-artifact@[0-9a-f]{40} # v4\.")
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
        self.assertEqual(
            len(re.findall(r"actions/download-artifact@[0-9a-f]{40} # v4\.", workflow)), 6
        )

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
            "frozen_mismatch_allowlist": {"sha256": "4" * 64, "gate_a_cases": []},
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
            "verification": {
                "gate": "LOSAT_CANONICAL",
                "classification": "CANONICAL_RAW_BYTES",
                "canonical_sha256": output_hash,
            },
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

    def test_aggregate_requires_every_complete_authority_bound_artifact(self) -> None:
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
            self.assertEqual(result["platform_count"], 4)
            self.assertEqual(result["total_search_executions"], 244)
            self.assertEqual(
                {item["platform_id"] for item in result["platforms"]},
                {"windows-x64", "macos-arm64", "macos-x64", "linux-aarch64"},
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

    # NCBI reference: ncbi-blast/c++/src/objtools/align_format/tabular.cpp:1098-1108
    #   x_PrintField(*iter); m_Ostream << "\n";
    # The tests below cover frozen_mismatch_allowlist.json in Gate A: raw bytes
    # are still compared by SHA-256, and a listed case passes only with its
    # recorded actual output.
    def plan_rows(self, catalog=None):
        topology = platform_cert.expected_completion_topology(catalog or self.catalog)
        return [
            (
                step["execution_index"],
                step["kind"],
                step["program"],
                step["case_id"],
                step["expected_losat_sha256"],
            )
            for step in topology.values()
        ]

    def fixture_catalog(self):
        canonical_rows = copy.deepcopy(self.catalog.canonical_rows)
        for row in canonical_rows:
            payload = f"matrix:{row['program']}:{row['case_id']}\n".encode()
            row["losat_sha256"] = hashlib.sha256(payload).hexdigest()
        canonical = {(row["program"], row["case_id"]): row for row in canonical_rows}
        return replace(self.catalog, canonical_rows=canonical_rows, canonical=canonical)

    def gate_a_step(self, root: Path, kind: str, payload: bytes, expected: bytes):
        writer = (
            "from pathlib import Path;import sys;"
            "Path(sys.argv[1]).parent.mkdir(parents=True,exist_ok=True);"
            "Path(sys.argv[1]).write_bytes(bytes.fromhex(sys.argv[2]))"
        )
        output_rel = f"{kind}/{SAKAI[0]}/{SAKAI[1]}/losat.out.partial"
        return platform_cert.SearchStep(
            execution_index=11,
            step_id=f"{kind}:{SAKAI[0]}:{SAKAI[1]}",
            kind=kind,
            program=SAKAI[0],
            case_id=SAKAI[1],
            semantic_class="SOURCE_UNDETERMINED_ACCEPTED",
            run_index=1,
            command=(
                sys.executable,
                "-c",
                writer,
                str(root / output_rel),
                payload.hex(),
            ),
            environment=(),
            output_rel=output_rel,
            final_output_rel=output_rel.removesuffix(".partial"),
            expected_losat_sha256=hashlib.sha256(expected).hexdigest(),
        )

    def test_repository_allowlist_is_empty_and_sakai_is_exact_in_v011(self) -> None:
        allowlist = platform_cert.load_frozen_allowlist(
            self.repo_root, self.candidate_sha
        )
        raw = platform_cert.read_git_blob_bytes(
            self.repo_root,
            self.candidate_sha,
            platform_cert.FROZEN_ALLOWLIST_REPO_RELATIVE,
        )
        self.assertEqual(allowlist.file_sha256, hashlib.sha256(raw).hexdigest())
        self.assertEqual(allowlist.entries, {})
        cases = platform_cert.bind_frozen_allowlist(allowlist, self.plan_rows())
        self.assertEqual(cases, [])
        self.assertEqual(platform_cert.allowed_known_mismatch_lines(cases), [])
        self.assertEqual(
            platform_cert.frozen_allowlist_identity("a" * 40, allowlist, cases),
            {
                "git_object": f"{'a' * 40}:LOSAT/tests/frozen_mismatch_allowlist.json",
                "sha256": allowlist.file_sha256,
                "gate_a_cases": [],
            },
        )
        topology = platform_cert.expected_completion_topology(self.catalog)
        self.assertEqual(
            [
                (step_id, step["expected_losat_sha256"])
                for step_id, step in topology.items()
                if step["execution_index"] in SAKAI_GATE_A_STEPS
            ],
            [
                ("matrix:blastn:Sakai.MG1655.megablast", SAKAI_APPROVED),
                ("repeatability:blastn:Sakai.MG1655.megablast:run2", SAKAI_APPROVED),
                ("repeatability:blastn:Sakai.MG1655.megablast:run3", SAKAI_APPROVED),
            ],
        )
        self.assertEqual(
            platform_cert.gate_a_verification(
                "matrix", *SAKAI, SAKAI_APPROVED, SAKAI_APPROVED, allowlist
            )["classification"],
            "CANONICAL_RAW_BYTES",
        )
        # The old frozen PR 5 bytes now fail Gate A, and the entry that allowed
        # the approved output against them no longer binds to the plan.
        with self.assertRaisesRegex(
            platform_cert.CertificationFailure, "CROSS_PLATFORM_LOSAT_DIVERGENCE"
        ):
            platform_cert.gate_a_verification(
                "matrix", *SAKAI, SAKAI_APPROVED, SAKAI_FROZEN_PR5, allowlist
            )
        historical = platform_cert.Allowlist(
            Path("frozen_mismatch_allowlist.json"),
            raw=json.dumps(
                {
                    "schema": "losat-frozen-mismatch-allowlist-v1",
                    "entries": [HISTORICAL_SAKAI_ENTRY],
                }
            ).encode(),
        )
        with self.assertRaisesRegex(
            platform_cert.CertificationFailure, "frozen hash outside the plan"
        ):
            platform_cert.bind_frozen_allowlist(historical, self.plan_rows())

    def test_frozen_allowlist_is_read_from_the_candidate_git_blob(self) -> None:
        relative = platform_cert.FROZEN_ALLOWLIST_REPO_RELATIVE.as_posix()
        raw = platform_cert.read_git_blob_bytes(
            self.repo_root, self.candidate_sha, relative
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate_sha = self.commit_files(root, {relative: raw})
            (root / relative).write_bytes(raw.replace(b"\n", b"\r\n") + b"edited")
            allowlist = platform_cert.load_frozen_allowlist(root, candidate_sha)
            self.assertEqual(allowlist.file_sha256, hashlib.sha256(raw).hexdigest())
            self.assertEqual(allowlist.entries, {})
        entry = dict(HISTORICAL_SAKAI_ENTRY)
        document = {**json.loads(raw), "entries": [entry]}
        listed = json.dumps(document).encode()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate_sha = self.commit_files(root, {relative: listed})
            allowlist = platform_cert.load_frozen_allowlist(root, candidate_sha)
            self.assertEqual(set(allowlist.entries), {SAKAI})
        rejected = {
            "duplicate case": json.dumps(
                {**document, "entries": [entry, entry]}
            ).encode(),
            "duplicate key": listed.replace(
                b'"schema":', b'"schema": "other", "schema":', 1
            ),
            "schema": json.dumps({**document, "schema": "other"}).encode(),
            "no entries": json.dumps({"schema": document["schema"]}).encode(),
            "not an object": b"[]",
            "not JSON": b"{",
        }
        other = {**entry, "case_id": "another_case"}
        other.pop("frozen_gate_a_steps")
        for label, entries in (
            ("program is not a string", [entry, {**other, "program": 1}]),
            ("case ID is null", [entry, {**other, "case_id": None}]),
            ("entry is not an object", [entry, "blastn"]),
            (
                "no frozen hash",
                [{k: v for k, v in entry.items() if k != "frozen_sha256"}],
            ),
            ("allowed is not a SHA-256", [{**entry, "allowed_actual_sha256": "AC81"}]),
            ("steps are not integers", [{**entry, "frozen_gate_a_steps": [11.0]}]),
            ("steps are not a list", [{**entry, "frozen_gate_a_steps": "11"}]),
        ):
            rejected[label] = json.dumps({**document, "entries": entries}).encode()
        for label, payload in rejected.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                candidate_sha = self.commit_files(root, {relative: payload})
                with self.assertRaises(platform_cert.CertificationFailure):
                    platform_cert.load_frozen_allowlist(root, candidate_sha)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            candidate_sha = self.commit_files(root, {"tracked.txt": b"tracked\n"})
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure, "candidate Git blob is missing"
            ):
                platform_cert.load_frozen_allowlist(root, candidate_sha)

    def test_allowlist_binding_requires_the_exact_gate_a_plan(self) -> None:
        plan = self.plan_rows(self.fixture_catalog())
        self.assertEqual(
            platform_cert.bind_frozen_allowlist(self.fixture_allowlist(), plan),
            self.fixture_allowed_cases(),
        )
        other_suite = {
            "program": "tblastn",
            "case_id": "not_in_gate_a",
            "frozen_sha256": "1" * 64,
            "allowed_actual_sha256": "2" * 64,
        }
        self.assertEqual(
            platform_cert.bind_frozen_allowlist(
                self.fixture_allowlist(other_suite), plan
            ),
            self.fixture_allowed_cases(),
        )
        representatives = {
            (program, case_id) for program, case_id, _ in platform_cert.REPRESENTATIVES
        }
        step_index, row = next(
            (index, row)
            for index, row in enumerate(self.catalog.canonical_rows, 1)
            if (row["program"], row["case_id"]) not in representatives
        )
        matrix_only = {
            "program": row["program"],
            "case_id": row["case_id"],
            "frozen_sha256": hashlib.sha256(
                f"matrix:{row['program']}:{row['case_id']}\n".encode()
            ).hexdigest(),
            "allowed_actual_sha256": "2" * 64,
            "frozen_gate_a_steps": [step_index],
        }
        self.assertEqual(
            sorted(
                case["frozen_gate_a_steps"]
                for case in platform_cert.bind_frozen_allowlist(
                    self.fixture_allowlist(matrix_only), plan
                )
            ),
            sorted([[step_index], SAKAI_GATE_A_STEPS]),
        )
        rejected = {
            "a step is missing": self.fixture_allowlist(frozen_gate_a_steps=[11, 52]),
            "another step": self.fixture_allowlist(frozen_gate_a_steps=[11, 52, 54]),
            "steps not declared": self.fixture_allowlist(frozen_gate_a_steps=None),
            "allowed equals frozen": self.fixture_allowlist(
                allowed_actual_sha256=hashlib.sha256(SAKAI_FIXTURE_FROZEN).hexdigest()
            ),
            "allowed is not a SHA-256": self.fixture_allowlist(
                allowed_actual_sha256="AC81"
            ),
            "steps for another frozen hash": self.fixture_allowlist(
                frozen_sha256="1" * 64
            ),
            "steps for a case outside the plan": self.fixture_allowlist(
                {**other_suite, "frozen_gate_a_steps": [1]}
            ),
        }
        for label, allowlist in rejected.items():
            with (
                self.subTest(label=label),
                self.assertRaisesRegex(
                    platform_cert.CertificationFailure, "frozen mismatch allowlist"
                ),
            ):
                platform_cert.bind_frozen_allowlist(allowlist, plan)

    def test_gate_a_allowlist_passes_only_the_recorded_actual_output(self) -> None:
        allowlist = self.fixture_allowlist()
        frozen = hashlib.sha256(SAKAI_FIXTURE_FROZEN).hexdigest()
        allowed = hashlib.sha256(SAKAI_FIXTURE_ALLOWED).hexdigest()
        self.assertEqual(
            platform_cert.gate_a_verification(
                "matrix", *SAKAI, frozen, allowed, allowlist
            ),
            {
                "gate": "LOSAT_CANONICAL",
                "classification": "ALLOWED_KNOWN_MISMATCH",
                "canonical_sha256": frozen,
                "allowed_actual_sha256": allowed,
                "allowlist_file_sha256": allowlist.file_sha256,
            },
        )
        self.assertEqual(
            platform_cert.gate_a_verification(
                "repeatability", *SAKAI, frozen, allowed, allowlist
            )["gate"],
            "LOSAT_REPEATABILITY",
        )
        with self.assertRaisesRegex(
            platform_cert.CertificationFailure, "STALE_FROZEN_MISMATCH_ALLOWLIST"
        ):
            platform_cert.gate_a_verification(
                "matrix", *SAKAI, frozen, frozen, allowlist
            )
        with self.assertRaisesRegex(
            platform_cert.CertificationFailure,
            f"CROSS_PLATFORM_LOSAT_DIVERGENCE.*allowed known mismatch {allowed}",
        ):
            platform_cert.gate_a_verification(
                "matrix", *SAKAI, frozen, "0" * 64, allowlist
            )
        # Without the allowlist, and for every case it does not list, only the
        # frozen hash passes.
        with self.assertRaisesRegex(
            platform_cert.CertificationFailure, "CROSS_PLATFORM_LOSAT_DIVERGENCE"
        ):
            platform_cert.gate_a_verification("matrix", *SAKAI, frozen, allowed, None)
        with self.assertRaisesRegex(
            platform_cert.CertificationFailure, "CROSS_PLATFORM_LOSAT_DIVERGENCE"
        ):
            platform_cert.gate_a_verification(
                "matrix", "blastn", "other_case", frozen, allowed, allowlist
            )
        self.assertEqual(
            platform_cert.gate_a_verification(
                "matrix", "blastn", "other_case", frozen, frozen, allowlist
            ),
            {
                "gate": "LOSAT_CANONICAL",
                "classification": "CANONICAL_RAW_BYTES",
                "canonical_sha256": frozen,
            },
        )

    def test_execution_reports_an_allowed_known_mismatch_and_fails_a_stale_entry(
        self,
    ) -> None:
        allowlist = self.fixture_allowlist()
        identity = {
            "source_sha": "a" * 40,
            "platform_id": "test",
            "losat": {"sha256": "b" * 64},
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            step = self.gate_a_step(
                root, "matrix", SAKAI_FIXTURE_ALLOWED, SAKAI_FIXTURE_FROZEN
            )
            record = platform_cert._execute_step(
                root, root, identity, self.catalog, step, None, allowlist
            )
            self.assertEqual(
                (root / step.final_output_rel).read_bytes(), SAKAI_FIXTURE_ALLOWED
            )
            self.assertEqual(
                record["expected_losat_sha256"],
                hashlib.sha256(SAKAI_FIXTURE_FROZEN).hexdigest(),
            )
            self.assertEqual(
                record["output_sha256"],
                hashlib.sha256(SAKAI_FIXTURE_ALLOWED).hexdigest(),
            )
            self.assertEqual(
                record["verification"]["classification"], "ALLOWED_KNOWN_MISMATCH"
            )
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure, "CROSS_PLATFORM_LOSAT_DIVERGENCE"
            ):
                platform_cert._execute_step(root, root, identity, self.catalog, step)
            stale = self.gate_a_step(
                root, "repeatability", SAKAI_FIXTURE_FROZEN, SAKAI_FIXTURE_FROZEN
            )
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure, "STALE_FROZEN_MISMATCH_ALLOWLIST"
            ):
                platform_cert._execute_step(
                    root, root, identity, self.catalog, stale, None, allowlist
                )
            other = self.gate_a_step(
                root, "matrix", b"another output\n", SAKAI_FIXTURE_FROZEN
            )
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure, "CROSS_PLATFORM_LOSAT_DIVERGENCE"
            ):
                platform_cert._execute_step(
                    root, root, identity, self.catalog, other, None, allowlist
                )

    def test_resume_accepts_an_allowed_output_only_with_the_allowlist(self) -> None:
        allowlist = self.fixture_allowlist()
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
            "frozen_mismatch_allowlist": platform_cert.frozen_allowlist_identity(
                "a" * 40, allowlist, self.fixture_allowed_cases()
            ),
            "controlled_fixtures": {
                "lexical_root": platform_cert.HISTORICAL_LEXICAL_ROOT,
                "staged_input_count": 34,
                "inputs": [],
            },
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prior = root / "prior"
            prior.mkdir()
            step = self.gate_a_step(
                prior, "matrix", SAKAI_FIXTURE_ALLOWED, SAKAI_FIXTURE_FROZEN
            )
            record = platform_cert._execute_step(
                prior, prior, identity, self.catalog, step, None, allowlist
            )
            (prior / "identity.json").write_text(json.dumps(identity), encoding="utf-8")
            current = root / "current"
            current.mkdir()
            imported = platform_cert.import_resume(
                prior,
                current,
                self.repo_root,
                identity,
                [step],
                self.catalog,
                self.authority,
                allowlist,
            )
            self.assertEqual(
                imported[step.step_id]["verification"], record["verification"]
            )
            self.assertEqual(
                (current / step.final_output_rel).read_bytes(), SAKAI_FIXTURE_ALLOWED
            )
            second = root / "second"
            second.mkdir()
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure, "resume LOSAT output changed"
            ):
                platform_cert.import_resume(
                    prior,
                    second,
                    self.repo_root,
                    identity,
                    [step],
                    self.catalog,
                    self.authority,
                )
            # The prior record's verification must be the one this allowlist gives.
            completion = platform_cert._step_record_path(prior, step)
            for label, mutate in (
                (
                    "another allowed hash",
                    lambda value: value["verification"].__setitem__(
                        "allowed_actual_sha256", "0" * 64
                    ),
                ),
                (
                    "another allowlist file",
                    lambda value: value["verification"].__setitem__(
                        "allowlist_file_sha256", "0" * 64
                    ),
                ),
                ("no verification", lambda value: value.pop("verification")),
                (
                    "verification is null",
                    lambda value: value.__setitem__("verification", None),
                ),
            ):
                with self.subTest(label=label):
                    value = json.loads(json.dumps(record))
                    mutate(value)
                    completion.write_text(json.dumps(value), encoding="utf-8")
                    target = root / f"rejected-{label.replace(' ', '-')}"
                    target.mkdir()
                    with self.assertRaisesRegex(
                        platform_cert.CertificationFailure,
                        "resume Gate A verification changed",
                    ):
                        platform_cert.import_resume(
                            prior,
                            target,
                            self.repo_root,
                            identity,
                            [step],
                            self.catalog,
                            self.authority,
                            allowlist,
                        )
            # A prior record with the frozen bytes of a listed case is stale.
            stale = json.loads(json.dumps(record))
            (prior / step.final_output_rel).write_bytes(SAKAI_FIXTURE_FROZEN)
            stale["output_sha256"] = hashlib.sha256(SAKAI_FIXTURE_FROZEN).hexdigest()
            stale["verification"] = {
                "gate": "LOSAT_CANONICAL",
                "classification": "CANONICAL_RAW_BYTES",
                "canonical_sha256": stale["output_sha256"],
            }
            completion.write_text(json.dumps(stale), encoding="utf-8")
            third = root / "third"
            third.mkdir()
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure, "resume LOSAT output changed"
            ):
                platform_cert.import_resume(
                    prior,
                    third,
                    self.repo_root,
                    identity,
                    [step],
                    self.catalog,
                    self.authority,
                    allowlist,
                )

    def test_final_evidence_names_the_allowed_known_mismatch(self) -> None:
        platform_id = "macos-arm64"
        catalog = self.fixture_catalog()
        allowed_cases = self.fixture_allowed_cases()
        allowed_sha256 = hashlib.sha256(SAKAI_FIXTURE_ALLOWED).hexdigest()
        identity = {
            "source_sha": "a" * 40,
            "platform_id": platform_id,
            "authority_version": self.authority.authority_version,
            "authority_file_sha256": self.authority.file_sha256,
        }

        def completions(steps, allowed_keys):
            completed = {}
            for step in steps:
                if step.kind == "oracle":
                    output_sha256 = "9" * 64
                    verification = {
                        "native_fingerprint": (
                            platform_cert.expected_native_fingerprint_evidence(
                                self.authority, platform_id, step.program, step.case_id
                            )
                        ),
                        "native_vs_losat_diagnostic": {
                            "classification": "DIAGNOSTIC_ERROR",
                            "platform_diagnostic": "none",
                        },
                    }
                else:
                    allowed = (step.program, step.case_id) in allowed_keys
                    output_sha256 = (
                        allowed_sha256 if allowed else step.expected_losat_sha256
                    )
                    verification = {
                        "classification": (
                            "ALLOWED_KNOWN_MISMATCH"
                            if allowed
                            else "CANONICAL_RAW_BYTES"
                        )
                    }
                completed[step.step_id] = {
                    "kind": step.kind,
                    "output_sha256": output_sha256,
                    "verification": verification,
                }
            return completed

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            steps = platform_cert.build_steps(
                self.repo_root,
                root,
                catalog,
                root / "LOSAT",
                {
                    "blastn": root / "blastn",
                    "blastp": root / "blastp",
                    "tblastx": root / "tblastx",
                },
            )
            output = root / "allowed"
            platform_cert.finalize_evidence(
                output, identity, steps, completions(steps, {SAKAI}), 61, allowed_cases
            )
            summary = json.loads((output / "summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["decision"], "PLATFORM_NATIVE_CERTIFIED")
            self.assertEqual(
                summary["gate_a_losat_canonical"],
                {
                    "total": 43,
                    "passed": 43,
                    "all_exact_canonical_raw_bytes": False,
                    "allowed_known_mismatches": allowed_cases,
                },
            )
            self.assertEqual(
                summary["cross_platform_losat"],
                {"total": 43, "all_exact_bytes": False, "allowed_known_mismatches": 1},
            )
            equalities = (output / "cross_platform_equalities.tsv").read_text(
                encoding="utf-8"
            )
            self.assertEqual(equalities.count("\tEXACT_BYTES\n"), 42)
            self.assertIn(
                f"blastn\tSakai.MG1655.megablast\t"
                f"{hashlib.sha256(SAKAI_FIXTURE_FROZEN).hexdigest()}\t"
                f"{allowed_sha256}\tALLOWED_KNOWN_MISMATCH\n",
                equalities,
            )
            direct = (output / "direct_ncbi_checks.tsv").read_text(encoding="utf-8")
            sakai_row = next(
                line for line in direct.splitlines() if "Sakai.MG1655.megablast" in line
            )
            self.assertIn(f"\t{allowed_sha256}\t", sakai_row)

            exact = root / "exact"
            platform_cert.finalize_evidence(
                exact, identity, steps, completions(steps, set()), 61
            )
            exact_summary = json.loads(
                (exact / "summary.json").read_text(encoding="utf-8")
            )
            self.assertEqual(
                exact_summary["gate_a_losat_canonical"],
                {
                    "total": 43,
                    "passed": 43,
                    "all_exact_canonical_raw_bytes": True,
                    "allowed_known_mismatches": [],
                },
            )
            self.assertEqual(
                exact_summary["cross_platform_losat"],
                {"total": 43, "all_exact_bytes": True, "allowed_known_mismatches": 0},
            )
            for label, allowed_keys, cases in (
                ("allowed step outside the allowlist", {SAKAI}, []),
                ("listed case recorded as exact", set(), allowed_cases),
            ):
                with (
                    self.subTest(label=label),
                    self.assertRaisesRegex(
                        platform_cert.CertificationFailure,
                        "Gate A allowlist coverage changed",
                    ),
                ):
                    platform_cert.finalize_evidence(
                        root / "rejected",
                        identity,
                        steps,
                        completions(steps, allowed_keys),
                        61,
                        cases,
                    )

    def test_aggregate_accepts_only_the_bound_allowed_known_mismatch(self) -> None:
        def update_json(path: Path, mutate) -> None:
            value = json.loads(path.read_text(encoding="utf-8"))
            mutate(value)
            platform_cert.atomic_write_json(path, value)

        allowlist = self.fixture_allowlist()
        frozen_sha256 = hashlib.sha256(SAKAI_FIXTURE_FROZEN).hexdigest()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_root = root / "input"
            input_root.mkdir()
            authority, catalog = self.create_aggregate_fixture(input_root, allowlist)
            result = platform_cert.aggregate_platform_evidence(
                input_root,
                root / "aggregate.json",
                "a" * 40,
                authority,
                catalog,
                allowlist,
            )
            self.assertEqual(result["decision"], "CROSS_PLATFORM_NATIVE_CERTIFIED")
            self.assertEqual(
                result["frozen_mismatch_allowlist"]["gate_a_cases"],
                self.fixture_allowed_cases(),
            )
            self.assertEqual(
                [
                    (item["gate_a_passed"], item["gate_a_allowed_known_mismatches"])
                    for item in result["platforms"]
                ],
                [(43, 1)] * 4,
            )
            # The same evidence is not certified when the aggregator's allowlist
            # does not list the case, or lists another output: the platforms'
            # allowlist identity is then not the aggregator's.
            for label, other in (
                ("no allowlist", None),
                (
                    "another allowed output",
                    self.fixture_allowlist(allowed_actual_sha256="7" * 64),
                ),
            ):
                with (
                    self.subTest(label=label),
                    self.assertRaisesRegex(
                        platform_cert.CertificationFailure,
                        "aggregate frozen mismatch allowlist identity failed",
                    ),
                ):
                    platform_cert.aggregate_platform_evidence(
                        input_root,
                        root / "rejected.json",
                        "a" * 40,
                        authority,
                        catalog,
                        other,
                    )

        def stale(artifact: Path, name: str) -> None:
            completion = artifact / "completions" / name
            value = json.loads(completion.read_text(encoding="utf-8"))
            (artifact / value["output_rel"]).write_bytes(SAKAI_FIXTURE_FROZEN)
            value["output_sha256"] = frozen_sha256
            value["verification"] = {
                "gate": value["verification"]["gate"],
                "classification": "CANONICAL_RAW_BYTES",
                "canonical_sha256": frozen_sha256,
            }
            platform_cert.atomic_write_json(completion, value)

        def third_output(artifact: Path) -> None:
            completion = (
                artifact / "completions" / "matrix__blastn__Sakai.MG1655.megablast.json"
            )
            value = json.loads(completion.read_text(encoding="utf-8"))
            payload = b"neither frozen nor allowed\n"
            (artifact / value["output_rel"]).write_bytes(payload)
            value["output_sha256"] = hashlib.sha256(payload).hexdigest()
            value["verification"]["allowed_actual_sha256"] = value["output_sha256"]
            platform_cert.atomic_write_json(completion, value)

        cases = (
            (
                "listed case has a third output and claims the allowed class",
                "aggregate Gate A completion failed",
                third_output,
            ),
            (
                "summary claims exact PR 5 bytes",
                "aggregate Gate A failed",
                lambda artifact: update_json(
                    artifact / "summary.json",
                    lambda value: value["gate_a_losat_canonical"].update(
                        all_exact_canonical_raw_bytes=True, allowed_known_mismatches=[]
                    ),
                ),
            ),
            (
                "another allowlist file",
                "aggregate frozen mismatch allowlist identity failed",
                lambda artifact: update_json(
                    artifact / "identity.json",
                    lambda value: value["frozen_mismatch_allowlist"].__setitem__(
                        "sha256", "0" * 64
                    ),
                ),
            ),
            (
                "listed matrix step is the frozen output again",
                "aggregate Gate A completion failed",
                lambda artifact: stale(
                    artifact, "matrix__blastn__Sakai.MG1655.megablast.json"
                ),
            ),
            (
                "listed repeatability step is the frozen output again",
                "aggregate repeatability completion failed",
                lambda artifact: stale(
                    artifact, "repeatability__blastn__Sakai.MG1655.megablast__run3.json"
                ),
            ),
            (
                "allowed classification on a case that is not listed",
                "aggregate Gate A completion failed",
                lambda artifact: update_json(
                    artifact
                    / "completions"
                    / "matrix__blastn__PesePMNV.MjPMNV.task_blastn.json",
                    lambda value: value["verification"].update(
                        classification="ALLOWED_KNOWN_MISMATCH",
                        allowed_actual_sha256=value["output_sha256"],
                        allowlist_file_sha256=allowlist.file_sha256,
                    ),
                ),
            ),
        )
        for label, message, mutate in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                input_root = root / "input"
                input_root.mkdir()
                authority, catalog = self.create_aggregate_fixture(
                    input_root, allowlist
                )
                artifact = input_root / "platform-native-macos-x64"
                mutate(artifact)
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
                        allowlist,
                    )
        # Evidence in which the listed case has its frozen bytes again on every
        # platform is not certified while the allowlist still lists it (a stale
        # entry); the platforms' allowlist identity is the aggregator's here.
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_root = root / "input"
            input_root.mkdir()
            authority, catalog = self.create_aggregate_fixture(input_root, allowlist)
            for artifact in sorted(input_root.iterdir()):
                for name in (
                    "matrix__blastn__Sakai.MG1655.megablast.json",
                    "repeatability__blastn__Sakai.MG1655.megablast__run2.json",
                    "repeatability__blastn__Sakai.MG1655.megablast__run3.json",
                ):
                    stale(artifact, name)
                self.refresh_evidence_manifest(artifact)
            with self.assertRaisesRegex(
                platform_cert.CertificationFailure,
                "aggregate Gate A completion failed",
            ):
                platform_cert.aggregate_platform_evidence(
                    input_root,
                    root / "aggregate.json",
                    "a" * 40,
                    authority,
                    catalog,
                    allowlist,
                )


if __name__ == "__main__":
    unittest.main()
