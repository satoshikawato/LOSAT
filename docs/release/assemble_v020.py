#!/usr/bin/env python3
"""Assemble and execute the three exact LOSAT v0.2.0 command artifacts."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = Path(__file__).with_name("v0.2.0_rc_contract.json")
HOST_FILES = {
    "serial_command_wasi": ("run_losat_wasi.js", "wasi_artifact.js"),
    "threaded_command_wasi": (
        "run_losat_wasi_threads.js",
        "wasi_thread_host.js",
        "wasi_shared_memory.js",
        "wasi_artifact.js",
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_hash(path: Path, expected: str) -> None:
    if not path.is_file() or sha256(path) != expected:
        raise RuntimeError(f"SHA-256 mismatch or missing input: {path}")


def run(command: list[str], *, cwd: Path | None = None) -> subprocess.CompletedProcess[bytes]:
    result = subprocess.run(command, cwd=cwd, capture_output=True, check=False)
    if result.returncode:
        raise RuntimeError(
            f"command failed ({result.returncode}): {command!r}\n"
            f"stderr: {result.stderr.decode(errors='replace')}"
        )
    return result


def archive_name(record: dict[str, str]) -> str:
    return record["filename"].removesuffix(".tar.gz")


def add_bytes(tar: tarfile.TarFile, name: str, data: bytes, executable: bool) -> None:
    info = tarfile.TarInfo(name)
    info.size = len(data)
    info.mode = 0o755 if executable else 0o644
    info.mtime = 0
    info.uid = info.gid = 0
    info.uname = info.gname = ""
    tar.addfile(info, io.BytesIO(data))


def assemble(record: dict[str, str], output: Path) -> None:
    ident = record["id"]
    prefix = archive_name(record)
    binary = ROOT / record["binary"]
    require_hash(binary, record["binary_sha256"])
    members = [(binary, "LOSAT" if ident == "native_linux_x64" else "LOSAT.wasm", True)]
    members += [(ROOT / "LOSAT/tests" / name, name, False) for name in HOST_FILES.get(ident, ())]
    members.append((ROOT / "LOSAT/LICENSE", "LICENSE", False))
    with output.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as zipped:
        with tarfile.open(fileobj=zipped, mode="w") as tar:
            for source, name, executable in members:
                add_bytes(tar, f"{prefix}/{name}", source.read_bytes(), executable)


def extract(archive: Path, destination: Path, prefix: str) -> Path:
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            if not member.isfile() or not member.name.startswith(prefix + "/"):
                raise RuntimeError(f"unsafe archive member: {member.name}")
            if ".." in Path(member.name).parts:
                raise RuntimeError(f"unsafe archive member: {member.name}")
        tar.extractall(destination, filter="data")
    return destination / prefix


# NCBI reference: c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
# const string kArgQuery("query"), kArgSubject("subject");
# const string kArgDbGeneticCode("db_gencode"), kArgNumThreads("num_threads");
# This is a package smoke invocation of the already certified local CLI contract.
def check_extracted(record: dict[str, str], archive: Path, smoke: dict[str, str]) -> dict[str, object]:
    ident = record["id"]
    with tempfile.TemporaryDirectory(prefix="losat-v020-extract-") as temp:
        temp_path = Path(temp)
        extracted = extract(archive, temp_path, archive_name(record))
        query = temp_path / "query.faa"
        subject = temp_path / "subject.fna"
        shutil.copyfile(ROOT / smoke["query"], query)
        shutil.copyfile(ROOT / smoke["subject"], subject)
        require_hash(query, smoke["query_sha256"])
        require_hash(subject, smoke["subject_sha256"])
        binary = extracted / ("LOSAT" if ident == "native_linux_x64" else "LOSAT.wasm")
        require_hash(binary, record["binary_sha256"])
        if ident == "native_linux_x64":
            prefix = [str(binary)]
            threads = "1"
        else:
            runner = "run_losat_wasi.js" if ident == "serial_command_wasi" else "run_losat_wasi_threads.js"
            prefix = ["node", "--no-warnings", "--experimental-wasi-unstable-preview1", str(extracted / runner), str(binary)]
            threads = "1" if ident == "serial_command_wasi" else "4"
        version = run(prefix + ["--version"]).stdout.decode().strip()
        if version != "losat 0.2.0":
            raise RuntimeError(f"unexpected extracted version: {version!r}")
        command = prefix + [
            "tblastn", "-task", "tblastn", "-query", str(query), "-subject", str(subject),
            "-db_gencode", "32", "-outfmt", "6", "-num_threads", threads,
        ]
        result = run(command)
        output_sha = hashlib.sha256(result.stdout).hexdigest()
        if output_sha != smoke["outfmt6_sha256"]:
            raise RuntimeError(f"extracted {ident} output mismatch: {output_sha}")
        return {
            "version": version,
            "command": [part.replace(str(temp_path), "{extract_root}") for part in command],
            "stdout_bytes": len(result.stdout),
            "stdout_sha256": output_sha,
            "stderr_sha256": hashlib.sha256(result.stderr).hexdigest(),
        }


def main() -> None:
    global ROOT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--crate", type=Path, required=True)
    args = parser.parse_args()
    ROOT = args.repo_root.resolve()
    contract = json.loads(CONTRACT.read_text())
    candidate = contract["candidate_sha"]
    run(["git", "cat-file", "-e", f"{candidate}^{{commit}}"], cwd=ROOT)
    # The script and evidence live outside the Cargo package; every package input
    # must still be byte-identical to the exact source candidate.
    run(["git", "diff", "--quiet", candidate, "--", "LOSAT", "README.md"], cwd=ROOT)
    if not str(args.output_dir.resolve()).startswith("/tmp/"):
        raise RuntimeError("output directory must be under /tmp")
    output_dir = args.output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise RuntimeError(f"output directory is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    crate = args.crate.resolve()
    require_hash(crate, contract["source_package"]["sha256"])
    crate_out = output_dir / contract["source_package"]["filename"]
    shutil.copyfile(crate, crate_out)
    records = []
    for artifact in contract["artifacts"]:
        archive = output_dir / artifact["filename"]
        assemble(artifact, archive)
        require_hash(archive, artifact["archive_sha256"])
        smoke = check_extracted(artifact, archive, contract["smoke"])
        records.append({
            "id": artifact["id"], "target": artifact["target"],
            "filename": artifact["filename"], "archive_bytes": archive.stat().st_size,
            "archive_sha256": sha256(archive), "binary_sha256": artifact["binary_sha256"],
            "extracted_smoke": smoke,
        })
    handoff = {
        "candidate_sha": candidate,
        "contract_sha256": sha256(CONTRACT),
        "source_package": {
            "filename": crate_out.name, "bytes": crate_out.stat().st_size,
            "sha256": sha256(crate_out),
        },
        "artifacts": records,
        "status": "LOCAL_ARTIFACTS_VERIFIED",
    }
    (output_dir / "handoff.json").write_text(json.dumps(handoff, indent=2) + "\n")
    names = [record["filename"] for record in records] + [crate_out.name, "handoff.json"]
    (output_dir / "SHA256SUMS").write_text(
        "".join(f"{sha256(output_dir / name)}  {name}\n" for name in names)
    )
    print(json.dumps({"status": handoff["status"], "candidate_sha": candidate,
                      "archive_hashes": {row["id"]: row["archive_sha256"] for row in records}}, indent=2))


if __name__ == "__main__":
    main()
