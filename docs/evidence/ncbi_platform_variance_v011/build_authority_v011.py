#!/usr/bin/env python3
"""Assemble LOSAT/tests/ncbi_platform_variance_v011.json from v0.1.0 and one added-platform record.

v0.1.1 = v0.1.0 with authority_version "ncbi-platform-variance-v0.1.1" and schema_version 2;
its representatives, three platform records, 18 diagnostics and characterization_provenance are
copied unchanged. The added platform's normative record is appended to
normative_acceptance.platforms, its six diagnostics to diagnostic_metadata, and its
characterization source to added_platform_characterizations. The record is the
`linux-aarch64-authority-record.json` that `classify-added-platform` of the candidate's
characterize_ncbi_platform_variance_v010.py wrote in the characterization run; its SHA-256 is
bound as derived_record_sha256.

usage: build_authority_v011.py --v010 V010.json --record RECORD.json --run-id N --output OUT.json
"""

import argparse
import hashlib
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v010", type=Path, required=True)
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--run-id", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    authority = json.loads(args.v010.read_bytes())
    record_bytes = args.record.read_bytes()
    record = json.loads(record_bytes)
    if authority["authority_version"] != "ncbi-platform-variance-v0.1.0":
        raise SystemExit("input is not authority v0.1.0")
    if record["status"] != "ADDED_PLATFORM_RECORD_DERIVED_NOT_ACCEPTED":
        raise SystemExit("input is not an added-platform record")
    platform_id = record["platform_id"]
    if any(
        item["platform_id"] == platform_id
        for item in authority["normative_acceptance"]["platforms"]
    ):
        raise SystemExit(f"{platform_id} is already in the authority")

    authority["authority_version"] = "ncbi-platform-variance-v0.1.1"
    authority["schema_version"] = 2
    authority["normative_acceptance"]["platforms"].append(
        record["normative_platform_record"]
    )
    authority["diagnostic_metadata"].extend(record["diagnostic_metadata"])
    authority["added_platform_characterizations"] = [
        {
            "candidate_sha": record["candidate_sha"],
            "characterization_run_id": args.run_id,
            "derived_record_sha256": hashlib.sha256(record_bytes).hexdigest(),
            "execution_counts": record["execution_counts"],
            "platform_runs": [record["platform_run"]],
            "runtime_enforced": False,
        }
    ]
    args.output.write_bytes(
        (json.dumps(authority, indent=2, sort_keys=True) + "\n").encode("utf-8")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
