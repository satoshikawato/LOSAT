#!/usr/bin/env python3
"""Build the database of the `db.xs.blastn.bat.fetch.*` rows (line RG item X-short-bat-2, 2026-11-24).

NCBI BLAST+ 2.17.0 makeblastdb is the comparison oracle only. `dmg` holds three random subjects of 60 bases with
parsed ids (s1, s2, s3); the header of s2 (OID 1) starts with 8 bytes 0xFF, so a read of it fails. Under
`-negative_seqidlist neg_s3.txt` (s3 is removed) NCBI reads the header of each sequence its search aligns again
(`check_oid_exclusion`): the query batch of qb (= s2) fails there, after the report of qa (= s1), whose search
aligns s1 only. The database is written with the creation date in its index file, so the frozen outputs need the
committed files, not a new run of this script.

usage: make_xsb2n.py --ncbi-bin DIR [--out LOSAT/tests/db_search/xsb2n]"""
import argparse, os, random, subprocess, tempfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ncbi-bin", required=True)
    parser.add_argument("--out", default=os.path.dirname(os.path.abspath(__file__)))
    args = parser.parse_args()
    out = args.out
    rng = random.Random(20261124)
    seqs = ["".join(rng.choice("ACGT") for _ in range(60)) for _ in range(3)]
    with open(f"{out}/dmg.fna", "w") as f:
        f.write("".join(f">s{k + 1}\n{s}\n" for k, s in enumerate(seqs)))
    with open(f"{out}/dq.fna", "w") as f:
        f.write(f">qa\n{seqs[0]}\n>qb\n{seqs[1]}\n")
    with open(f"{out}/neg_s3.txt", "w") as f:
        f.write("s3\n")
    makeblastdb = os.path.join(args.ncbi_bin, "makeblastdb")
    env = dict(os.environ, BLAST_USAGE_REPORT="false")
    subprocess.run([makeblastdb, "-in", "dmg.fna", "-dbtype", "nucl", "-parse_seqids", "-out", "dmg"],
                   cwd=out, env=env, check=True, stdout=subprocess.DEVNULL)
    # The header of s2 starts where s1's ends: the length of the headers of a database of s1 alone.
    with tempfile.TemporaryDirectory() as tmp:
        with open(f"{tmp}/one.fna", "w") as f:
            f.write(f">s1\n{seqs[0]}\n")
        subprocess.run([makeblastdb, "-in", "one.fna", "-dbtype", "nucl", "-parse_seqids", "-out", "one"],
                       cwd=tmp, env=env, check=True, stdout=subprocess.DEVNULL)
        at = os.path.getsize(f"{tmp}/one.nhr")
    with open(f"{out}/dmg.nhr", "r+b") as f:
        f.seek(at)
        f.write(b"\xff" * 8)


if __name__ == "__main__":
    main()
