#!/usr/bin/env python3
"""Freezes NCBI BLAST+ 2.17.0 windowmasker's output for each row of cases.tsv (run in this
directory, as tests/masker_fixtures.rs runs LOSAT): ncbi/<name>.out, .err, .rc, with this
directory's path in the standard error written as @CWD@. A row whose words name @OUT@ runs with
a new empty directory there; the files the run leaves in it follow its standard output, each as
a line "@OUT@/<name> <length>" after a newline and the file's bytes (tests/masker_fixtures.rs
reads LOSAT's runs the same way), and the directory's path in the standard error is written
@OUT@. A row's first words of the form NAME=value set that variable for its run (BLASTDB: where
the taxonomy database is). First makes db/wm and db/wm_ps with NCBI makeblastdb from in/multi.fa
(the second with -parse_seqids), db/multi_ids and db/wm_b from in/multi_ids.fa (the first with
-parse_seqids), db/cd from in/cd.fa, db/wm_tx from in/tx.fa (-parse_seqids -taxid_map
in/tx_taxid.txt), and the counts files of the -convert rows in counts/ with NCBI windowmasker
(units of 4 letters of dm2L_60k.fa in each format; the ascii file without its ##pct lines, with
a unit size of 0 and of 256, and with ##pct lines whose values NCBI fails to read; the binary
file with a cut ##pct line; the optimized binary file of -t_high 1000, whose counts above 511
overflow; the files of -t_high_pct 0; the binary counts written to the standard output). NCBI is
run here only, never by the tests.

Usage: python3 -I make_fixtures.py --ncbi-bin <dir of windowmasker and makeblastdb> [--only ROW...]
(--only: those rows, and only the databases that are missing; the other rows' files stay.)"""
import argparse
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

# Rows whose standard output and error NCBI also writes into one file (`2>&1`, a run of its
# own): ncbi/<name>.merged. The reader's warnings come between records, so where they fall among
# the output's bytes shows NCBI's stdout buffer (`cerr` is tied to `cout`).
MERGED = {"r_letters", "r_letters_fasta", "r_letters_slb", "r_gap_fasta_ps", "r_gap_mi",
          "d_letters_fasta", "k_binary", "cd_two", "v_binary_stdout"}


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


def make_counts(windowmasker, env):
    """The counts files of the -convert rows, in counts/."""
    counts = os.path.join(HERE, "counts")
    shutil.rmtree(counts, ignore_errors=True)
    os.makedirs(counts)
    base = ["-mk_counts", "-in", "../winmask/fasta/dm2L_60k.fa", "-unit", "4"]
    runs = [(f"u4.{f}", ["-sformat", f, "-meta", "dm2L 60k"]) for f in ("ascii", "binary", "oascii", "obinary")]
    runs += [("u4_ovf.obinary", ["-sformat", "obinary", "-t_high", "1000"]),
             ("u4_thpct0.ascii", ["-t_high_pct", "0"]),
             ("u4_thpct0.obinary", ["-sformat", "obinary", "-t_high_pct", "0"])]
    for name, extra in runs:
        subprocess.run([windowmasker] + base + extra + ["-out", f"counts/{name}"], cwd=HERE, env=env,
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, stdin=subprocess.DEVNULL)
    r = subprocess.run([windowmasker] + base + ["-sformat", "binary"], cwd=HERE, env=env, check=True,
                       capture_output=True, stdin=subprocess.DEVNULL)
    with open(os.path.join(counts, "u4_stdout.binary"), "wb") as f:
        f.write(r.stdout)
    with open(os.path.join(counts, "u4.ascii"), "rb") as f:
        lines = f.read().split(b"\n")
    with open(os.path.join(counts, "u4_nopct.ascii"), "wb") as f:
        f.write(b"\n".join(line for line in lines if not line.startswith(b"##pct:")))
    # The unit-size line (the first line after the metadata) 0, and 256 (0 as NCBI's Uint1).
    first = next(i for i, line in enumerate(lines) if not line.startswith(b"##"))
    for name, size in (("u0.ascii", b"0"), ("u256.ascii", b"256")):
        with open(os.path.join(counts, name), "wb") as f:
            f.write(b"\n".join(lines[:first] + [size] + lines[first + 1:]))
    # ##pct: lines that NCBI's istringstream reads with failed or missing values (each 0).
    words = [b"##pct: 0 0", b"##pct: 1", b"##pct: abc 50", b"##pct: 3 xyz", b"##pct: -4 12.5",
             b"##pct: 99999999999 13", b"##pct: 6 99.5 extra words", b"##pct:7 100"]
    kept = [line for line in lines if not line.startswith(b"##pct:")]
    at = next(i for i, line in enumerate(kept) if not line.startswith(b"##"))
    with open(os.path.join(counts, "u4_pctwords.ascii"), "wb") as f:
        f.write(b"\n".join(kept[:at] + words + kept[at:]))
    # The binary file's ##pct: line cut after two values, a word that is not a number between.
    with open(os.path.join(counts, "u4.binary"), "rb") as f:
        data = f.read()
    length = int.from_bytes(data[4:8], "little")
    metadata = b"\0".join(b"##pct: 500 0 0.1 abc 0.3" if field.startswith(b"##pct:") else field
                          for field in data[8:8 + length].split(b"\0"))
    with open(os.path.join(counts, "u4_pctcut.binary"), "wb") as f:
        f.write(data[:4] + len(metadata).to_bytes(4, "little") + metadata + data[8 + length:])


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
                                ("cd", "cd.fa", []),
                                ("wm_tx", "tx.fa", ["-parse_seqids", "-taxid_map", "in/tx_taxid.txt"])):
        if os.path.exists(os.path.join(db, f"{name}.nsq")):
            continue
        subprocess.run([makeblastdb, "-in", f"in/{source}", "-dbtype", "nucl", "-out", f"db/{name}"] + extra,
                       cwd=HERE, env=env, check=True, stdout=subprocess.DEVNULL, stdin=subprocess.DEVNULL)
    make_counts(windowmasker, env)
    out_dir = os.path.join(HERE, "ncbi")
    if not a.only:
        shutil.rmtree(out_dir, ignore_errors=True)
    os.makedirs(out_dir, exist_ok=True)
    for name, words, stdin, variables in cases():
        if a.only and name not in a.only:
            continue
        source = open(os.path.join(HERE, stdin), "rb") if stdin else subprocess.DEVNULL
        out = tempfile.mkdtemp(prefix=f"wm-{name}-") if any("@OUT@" in w for w in words) else None
        if out:
            words = [w.replace("@OUT@", out) for w in words]
        r = subprocess.run([windowmasker] + words, cwd=HERE, env={**env, **variables}, stdin=source,
                           capture_output=True, timeout=600)
        status = r.returncode if r.returncode >= 0 else 128 - r.returncode
        stdout, err = r.stdout, r.stderr
        if out:
            for file in sorted(os.listdir(out)):
                data = open(os.path.join(out, file), "rb").read()
                stdout += f"\n@OUT@/{file} {len(data)}\n".encode() + data
            shutil.rmtree(out)
            err = err.replace(out.encode(), b"@OUT@")
        err = err.replace(HERE.encode(), b"@CWD@")
        for ext, data in (("out", stdout), ("err", err), ("rc", f"{status}\n".encode())):
            with open(os.path.join(out_dir, f"{name}.{ext}"), "wb") as f:
                f.write(data)
        print(f"{name}\t{status}\t{len(stdout)}\t{len(err)}")
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
