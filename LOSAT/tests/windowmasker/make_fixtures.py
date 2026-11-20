#!/usr/bin/env python3
"""Freezes NCBI BLAST+ 2.17.0 windowmasker's output for each row of cases.tsv (run in this
directory, as tests/masker_fixtures.rs runs LOSAT): ncbi/<name>.out, .err, .rc, with this
directory's path in the standard error written as @CWD@. A row whose words name @OUT@ runs with
a new empty directory there; the files the run leaves in it follow its standard output, each as
a line "@OUT@/<name> <length>" after a newline and the file's bytes (tests/masker_fixtures.rs
reads LOSAT's runs the same way), and the directory's path in the standard error is written
@OUT@. First makes db/wm and db/wm_ps with NCBI makeblastdb from in/multi.fa (the second with
-parse_seqids), db/cd from in/cd.fa, and the counts files of the -convert rows in counts/ with
NCBI windowmasker (units of 4 letters of dm2L_60k.fa in each format; the ascii file without its
##pct lines; the optimized binary file of -t_high 1000, whose counts above 511 overflow; the
files of -t_high_pct 0; the binary counts written to the standard output). NCBI is run here
only, never by the tests.

Usage: python3 -I make_fixtures.py --ncbi-bin <dir of windowmasker and makeblastdb>"""
import argparse
import os
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
        stdin = None
        if "<" in words:
            at = words.index("<")
            stdin = words[at + 1]
            del words[at:at + 2]
        rows.append((name, words, stdin))
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi-bin", required=True)
    a = ap.parse_args()
    windowmasker = os.path.join(a.ncbi_bin, "windowmasker")
    makeblastdb = os.path.join(a.ncbi_bin, "makeblastdb")
    env = environment()
    db = os.path.join(HERE, "db")
    shutil.rmtree(db, ignore_errors=True)
    os.makedirs(db)
    for name, fasta, extra in (("wm", "multi", []), ("wm_ps", "multi", ["-parse_seqids"]), ("cd", "cd", [])):
        subprocess.run([makeblastdb, "-in", f"in/{fasta}.fa", "-dbtype", "nucl", "-out", f"db/{name}"] + extra,
                       cwd=HERE, env=env, check=True, stdout=subprocess.DEVNULL, stdin=subprocess.DEVNULL)
    make_counts(windowmasker, env)
    out_dir = os.path.join(HERE, "ncbi")
    shutil.rmtree(out_dir, ignore_errors=True)
    os.makedirs(out_dir)
    for name, words, stdin in cases():
        source = open(os.path.join(HERE, stdin), "rb") if stdin else subprocess.DEVNULL
        out = tempfile.mkdtemp(prefix=f"wm-{name}-") if any("@OUT@" in w for w in words) else None
        if out:
            words = [w.replace("@OUT@", out) for w in words]
        r = subprocess.run([windowmasker] + words, cwd=HERE, env=env, stdin=source,
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
                subprocess.run([windowmasker] + words, cwd=HERE, env=env, stdin=source,
                               stdout=one, stderr=one, timeout=600)
            with open(os.path.join(out_dir, f"{name}.merged"), "rb") as f:
                merged = f.read().replace(HERE.encode(), b"@CWD@")
            with open(os.path.join(out_dir, f"{name}.merged"), "wb") as f:
                f.write(merged)
    return 0


if __name__ == "__main__":
    sys.exit(main())
