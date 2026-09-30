#!/usr/bin/env python3
"""Compare BLASTN input handling and query batches with NCBI (Session S07++; comparison only).

This is docs/evidence/losat_web_e2c/check_inputs.py with the expectations of S07++. LOSAT
now searches NCBI's adaptive query batches, so the cases that S07+ rejected because they
depended on the batches after the first are compared byte for byte ("same"). These are the
reports of invalid queries after the first batch, and the gapped X-drop of gap costs beyond
the tables with queries of different compositions. A Karlin-Altschul table error after a
first batch of invalid queries stays an explicit rejection. Two cases more check the batch
size from which NCBI splits a batch into query chunks, which LOSAT rejects explicitly.

Usage: check_inputs.py --bin-dir DIR --losat LOSAT --work DIR
"""
import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("e2c_check_inputs", HERE.parent / "losat_web_e2c" / "check_inputs.py")
e2c = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(e2c)

BATCHED = {"invalid_at_end.fmt7", "long_invalid_run.fmt0"}
e2c_cases = e2c.cases
e2c_make_inputs = e2c.make_inputs
# NCBI splits a query batch of at least 2 x (1000000 - 100) residues with -task blastn into
# query chunks (CQuerySplitter), which LOSAT rejects; one residue less is not split.
SPLIT = 2 * (1_000_000 - 100)


def make_inputs(work: Path) -> None:
    e2c_make_inputs(work)
    edl933 = "".join(line.strip() for line in (e2c.ENGINE / "tests/fasta/EDL933.fna").read_text().splitlines()[1:])
    (work / "split_query.fa").write_text(">split\n" + edl933[:SPLIT] + "\n")
    (work / "unsplit_query.fa").write_text(">unsplit\n" + edl933[:SPLIT - 1] + "\n")
    (work / "split_subject.fa").write_text(">window\n" + edl933[1_105_000:1_125_000] + "\n")


def cases(work: Path) -> list[tuple]:
    rows = []
    for name, argv, expect, *extra in e2c_cases(work):
        if name in BATCHED or name.startswith("batches.many_iupac."):
            expect = "same"
        rows.append((name, argv, expect, *extra))
    split = ["-subject", f"{work}/split_subject.fa", "-task", "blastn", "-outfmt", "6"]
    rows += [("s07pp.query_split", ["-query", f"{work}/split_query.fa", *split], "losat-rejects"),
             ("s07pp.query_unsplit", ["-query", f"{work}/unsplit_query.fa", *split], "same")]
    return rows


e2c.cases = cases
e2c.make_inputs = make_inputs

if __name__ == "__main__":
    sys.exit(e2c.main())
