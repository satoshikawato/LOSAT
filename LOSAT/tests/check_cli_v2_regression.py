#!/usr/bin/env python3
"""Compare CLI v2 with frozen pre-migration output and the NCBI oracle."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess


# NCBI c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-94:
# kArgQuery("query"), kArgSubject("subject"), kArgNumThreads("num_threads").
# Frozen before CLI v2 from e2abf5848b09309c761c44a373d3bd188daafd80
# plus the existing working tree. All four outputs also matched NCBI 2.17.0+.
CASES = [
    ("blastn", "blastn_parity_compact.fasta", "blastn_parity_compact.fasta",
     ["-task", "blastn", "-word_size", "11"], "6",
     "ac4a20b64c8da3eaea725ed09dc89880294f58e971b7207c943489faacdbe03e"),
    ("blastn", "blastn_parity_compact.fasta", "blastn_parity_compact.fasta",
     ["-task", "blastn", "-word_size", "11"], "7",
     "e8b6e54f9f42f7f73f1d91c0cf8dc465942b8234f22eb5826f7253f084c9baf2"),
    ("blastp", "SicyWSV.faa", "PajaWSV.faa", ["-max_hsps", "1"], "6",
     "fd4b010800e32ce6c823cb38b42a10b7845f3342edae892acccc8f554f9edf34"),
    ("tblastx", "LC738874.fasta", "LC738875.fasta", [], "6",
     "86c05a04efb50e4026720e2d44fe2db2e6446f9594e174f3fde56931d09d5b49"),
]


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--losat-bin", type=Path, default=root / "LOSAT/target/release/LOSAT")
    parser.add_argument("--ncbi-dir", type=Path, help="Otherwise use NCBI programs on PATH")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    environment = {key: value for key, value in os.environ.items() if not key.startswith("LOSAT_")}
    records = []
    for program, query, subject, extra, outfmt, expected in CASES:
        oracle = str(args.ncbi_dir / program) if args.ncbi_dir else program
        version = subprocess.check_output([oracle, "-version"], text=True).splitlines()[0]
        for threads in (1, 2):
            command = [program, "-query", f"LOSAT/tests/fasta/{query}",
                       "-subject", f"LOSAT/tests/fasta/{subject}", "-outfmt", outfmt,
                       "-num_threads", str(threads), "-evalue", "10", *extra]
            case = f"{program}-n{threads}-fmt{outfmt}"
            result = subprocess.run([str(args.losat_bin.resolve()), *command], cwd=root,
                                    env=environment, capture_output=True, check=True)
            (args.output_dir / f"{case}.raw").write_bytes(result.stdout)
            digest = hashlib.sha256(result.stdout).hexdigest()
            record = {"case": case, "argv": command, "sha256": digest,
                      "baseline_equal": digest == expected, "ncbi_version": version}
            if threads == 1:
                reference = subprocess.run([oracle, *command[1:]], cwd=root,
                                           capture_output=True, check=True)
                (args.output_dir / f"{case}.ncbi.raw").write_bytes(reference.stdout)
                record["ncbi_equal"] = result.stdout == reference.stdout
            records.append(record)
            (args.output_dir / "results.json").write_text(json.dumps(records, indent=2) + "\n")
            if not record["baseline_equal"] or not record.get("ncbi_equal", True):
                raise SystemExit(f"raw-byte regression: {case}; see {args.output_dir}")
            print(f"{case}: exact", flush=True)


    common = ["blastp", "-query", "LOSAT/tests/fasta/SicyWSV.faa",
              "-subject", "LOSAT/tests/fasta/PajaWSV.faa", "-outfmt", "6",
              "-num_threads", "1", "-max_hsps", "1"]
    # NCBI api/blast_options_handle.cpp:381-400: Create(eBlastp, locality).
    # v0.1.0 exposes only the ordinary, comparison-supported task.
    command = common + ["-task", "blastp"]
    result = subprocess.run([str(args.losat_bin.resolve()), *command], cwd=root,
                            env=environment, capture_output=True, check=True)
    digest = hashlib.sha256(result.stdout).hexdigest()
    expected = "fd4b010800e32ce6c823cb38b42a10b7845f3342edae892acccc8f554f9edf34"
    (args.output_dir / "blastp-explicit-task.raw").write_bytes(result.stdout)
    records.append({"case": "blastp-explicit-task", "argv": command,
                    "sha256": digest, "baseline_equal": digest == expected})
    if digest != expected:
        raise SystemExit("explicit ordinary blastp task changed output")
    for name, extra in [("short-hidden", ["-task", "blastp-short"]),
                        ("fast-hidden", ["-task", "blastp-fast"]),
                        ("sw-removed", ["-use_sw_tback"])]:
        command = common + extra
        result = subprocess.run([str(args.losat_bin.resolve()), *command], cwd=root,
                                env=environment, capture_output=True)
        (args.output_dir / f"{name}.stderr").write_bytes(result.stderr)
        rejected = result.returncode == 2 and result.stdout == b"" and extra[-1].encode() in result.stderr
        records.append({"case": name, "argv": command, "rejected_at_parser": rejected,
                        "exit_code": result.returncode})
        if not rejected:
            raise SystemExit(f"public capability was not rejected at parsing: {name}")
        print(f"{name}: rejected at parser", flush=True)
    (args.output_dir / "results.json").write_text(json.dumps(records, indent=2) + "\n")


if __name__ == "__main__":
    main()
