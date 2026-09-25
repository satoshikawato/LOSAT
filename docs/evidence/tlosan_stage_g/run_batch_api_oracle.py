#!/usr/bin/env python3
"""Run the pinned comparison-only NCBI API oracle at CLI TBLASTN query boundaries."""
from __future__ import annotations

import hashlib
import subprocess
import sys
import tempfile
from pathlib import Path

# NCBI c++/src/app/blast/tblastn_app.cpp:250-252,275-301:
# CBlastInput(..., GetQueryBatchSize()); GetNextSeqBatch; CLocalBlast.Run()
# NCBI c++/src/algo/blast/blastinput/blast_input_aux.cpp:105-127:
# case eTblastn: retval = 20000;
# NCBI c++/src/algo/blast/blastinput/blast_input.cpp:135-165:
# the sequence that makes size_read >= batch_size stays in that batch.
BATCH_RESIDUES = 20_000


def batches(path: Path) -> list[bytes]:
    records: list[tuple[bytes, int]] = []
    header = b""
    sequence: list[bytes] = []
    for line in path.read_bytes().splitlines(keepends=True):
        if line.startswith(b">"):
            if header:
                records.append((header + b"".join(sequence), sum(len(part.strip()) for part in sequence)))
            header, sequence = line, []
        else:
            sequence.append(line)
    if header:
        records.append((header + b"".join(sequence), sum(len(part.strip()) for part in sequence)))
    if not records:
        raise ValueError("no FASTA records")
    output: list[bytes] = []
    index = 0
    while index < len(records):
        chunk: list[bytes] = []
        residues = 0
        while index < len(records) and residues < BATCH_RESIDUES:
            record, size = records[index]
            chunk.append(record if record.endswith(b"\n") else record + b"\n")
            residues += size
            index += 1
        output.append(b"".join(chunk))
    return output


def compose(parts: list[bytes], fmt: int, query_count: int) -> bytes:
    if len(parts) == 1:
        return parts[0]
    if fmt == 6:
        return b"".join(parts)
    if fmt == 7:
        marker = b"# BLAST processed "
        bodies = []
        for part in parts:
            end = part.rfind(marker)
            if end < 0:
                raise ValueError("API outfmt7 processed-query footer missing")
            bodies.append(part[:end])
        return b"".join(bodies) + f"# BLAST processed {query_count} queries\n".encode()
    if fmt == 0:
        query_marker = b"\n\nQuery= "
        footer_marker = b"\n\n  Database: "
        intro = parts[0][:parts[0].find(query_marker)]
        sections = []
        for part in parts:
            start, end = part.find(query_marker), part.rfind(footer_marker)
            if start < 0 or end < start:
                raise ValueError("API outfmt0 query or final database section missing")
            sections.append(part[start:end])
        footer = parts[-1][parts[-1].rfind(footer_marker):]
        return intro + b"".join(sections) + footer
    raise ValueError(f"unsupported outfmt {fmt}")


def main() -> int:
    if len(sys.argv) != 6:
        raise SystemExit("usage: run_batch_api_oracle.py API QUERY SUBJECT CODE OUTFMT")
    api, query, subject = map(Path, sys.argv[1:4])
    code, fmt = int(sys.argv[4]), int(sys.argv[5])
    chunks = batches(query)
    parts = []
    with tempfile.TemporaryDirectory(prefix="tlosan-stageg-api-batches-") as directory:
        for index, chunk in enumerate(chunks):
            fixture = Path(directory) / f"query_batch_{index}.faa"
            fixture.write_bytes(chunk)
            command = [str(api), str(fixture), str(subject), str(code), str(fmt)]
            result = subprocess.run(command, capture_output=True)
            if result.returncode:
                sys.stderr.buffer.write(result.stderr)
                return result.returncode
            parts.append(result.stdout)
            print(f"API_BATCH\t{index}\t{hashlib.sha256(chunk).hexdigest()}\t{len(chunk)}", file=sys.stderr)
    sys.stdout.buffer.write(compose(parts, fmt, sum(chunk.count(b">", 0, len(chunk)) for chunk in chunks)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
