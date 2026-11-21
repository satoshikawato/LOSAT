#!/usr/bin/env python3
"""Freezes NCBI BLAST+ 2.17.0 windowmasker's output for each row of cases.tsv (run in this
directory, as tests/masker_fixtures.rs runs LOSAT): ncbi/<name>.out, .err, .rc, with this
directory's path in the standard error written as @CWD@. First makes db/wm and db/wm_ps with
NCBI makeblastdb from in/multi.fa (the second with -parse_seqids), db/multi_ids and db/wm_b
from in/multi_ids.fa (the first with -parse_seqids), and db/wm_tx from in/tx.fa (-parse_seqids
-taxid_map in/tx_taxid.txt). NCBI is run here only, never by the tests. A row's first words of
the form NAME=value set that variable for its run (BLASTDB: where the taxonomy database is).

Usage: python3 -I make_fixtures.py --ncbi-bin <dir of windowmasker and makeblastdb> [--only ROW...]
(--only: those rows, and only the databases that are missing; the other rows' files stay.)"""
import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# Rows whose standard output and error NCBI also writes into one file (`2>&1`, a run of its
# own): ncbi/<name>.merged. The reader's warnings come between records, so where they fall among
# the output's bytes shows NCBI's stdout buffer (`cerr` is tied to `cout`).
MERGED = {"r_letters", "r_letters_fasta", "r_letters_slb", "r_gap_fasta_ps", "r_gap_mi",
          "d_letters_fasta"}


def environment():
    env = {k: v for k, v in os.environ.items()
           if not k.startswith(("LOSAT_", "NCBI_", "BLAST"))}
    env["NCBI_CONFIG_PATH"] = HERE
    env["BLAST_USAGE_REPORT"] = "false"
    return env


def cases():
    rows = []
    for line in open(os.path.join(HERE, "cases.tsv")):
        line = line.rstrip("\n")
        if not line or line.startswith("#"):
            continue
        name, arguments = line.split("\t", 1)
        words = shlex.split(arguments)
        variables = {}
        while words and re.fullmatch(r"[A-Z_]+=.*", words[0], re.S):
            name_, value = words.pop(0).split("=", 1)
            variables[name_] = value
        stdin = None
        if "<" in words:
            at = words.index("<")
            stdin = words[at + 1]
            del words[at:at + 2]
        rows.append((name, words, stdin, variables))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi-bin", required=True)
    ap.add_argument("--only", nargs="+")
    a = ap.parse_args()
    windowmasker = os.path.join(a.ncbi_bin, "windowmasker")
    makeblastdb = os.path.join(a.ncbi_bin, "makeblastdb")
    env = environment()
    db = os.path.join(HERE, "db")
    if not a.only:
        shutil.rmtree(db, ignore_errors=True)
    os.makedirs(db, exist_ok=True)
    for name, source, extra in (("wm", "multi.fa", []), ("wm_ps", "multi.fa", ["-parse_seqids"]),
                                ("multi_ids", "multi_ids.fa", ["-parse_seqids"]), ("wm_b", "multi_ids.fa", []),
                                ("wm_tx", "tx.fa", ["-parse_seqids", "-taxid_map", "in/tx_taxid.txt"])):
        if os.path.exists(os.path.join(db, f"{name}.nsq")):
            continue
        subprocess.run([makeblastdb, "-in", f"in/{source}", "-dbtype", "nucl", "-out", f"db/{name}"] + extra,
                       cwd=HERE, env=env, check=True, stdout=subprocess.DEVNULL, stdin=subprocess.DEVNULL)
    out_dir = os.path.join(HERE, "ncbi")
    if not a.only:
        shutil.rmtree(out_dir, ignore_errors=True)
    os.makedirs(out_dir, exist_ok=True)
    for name, words, stdin, variables in cases():
        if a.only and name not in a.only:
            continue
        source = open(os.path.join(HERE, stdin), "rb") if stdin else subprocess.DEVNULL
        r = subprocess.run([windowmasker] + words, cwd=HERE, env={**env, **variables}, stdin=source,
                           capture_output=True, timeout=600)
        status = r.returncode if r.returncode >= 0 else 128 - r.returncode
        err = r.stderr.replace(HERE.encode(), b"@CWD@")
        for ext, data in (("out", r.stdout), ("err", err), ("rc", f"{status}\n".encode())):
            with open(os.path.join(out_dir, f"{name}.{ext}"), "wb") as f:
                f.write(data)
        print(f"{name}\t{status}\t{len(r.stdout)}\t{len(err)}")
        if name in MERGED:
            source = open(os.path.join(HERE, stdin), "rb") if stdin else subprocess.DEVNULL
            with open(os.path.join(out_dir, f"{name}.merged"), "wb") as one:
                subprocess.run([windowmasker] + words, cwd=HERE, env={**env, **variables}, stdin=source,
                               stdout=one, stderr=one, timeout=600)
            with open(os.path.join(out_dir, f"{name}.merged"), "rb") as f:
                merged = f.read().replace(HERE.encode(), b"@CWD@")
            with open(os.path.join(out_dir, f"{name}.merged"), "wb") as f:
                f.write(merged)
    return 0


if __name__ == "__main__":
    sys.exit(main())
