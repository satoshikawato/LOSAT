#!/usr/bin/env python3
"""Build the masked database fixtures of LOSAT/tests/seqdb (mask data of makeblastdb -mask_data).

NCBI BLAST+ 2.17.0 is the comparison oracle only: dustmasker and segmasker write mask
files, makeblastdb writes the databases with their mask column (`.naa .nab .nac`, `.paa .pab .pac`),
and blastdbcmd writes what NCBI reads back (expected/mask_*.tsv, expected/mask_*.algorithms.txt).
The Rust tests read the committed files and never run an NCBI program.

usage: make_mask_fixtures.py --ncbi-bin DIR --work DIR [--out LOSAT/tests/seqdb]

No database here uses -parse_seqids: the masks name the sequences by their number in the
FASTA file (`local id N`), which is what both dustmasker and makeblastdb do without the flag.
Rerunning changes the dates stored in the index and column files.
"""
import argparse
import random
import re
import shutil
import subprocess
from pathlib import Path

NUCL_EXTS = [".nin", ".nsq", ".nhr", ".naa", ".nab", ".nac"]
PROT_EXTS = [".pin", ".psq", ".phr", ".paa", ".pab", ".pac"]


def write_fasta(path, records, width=60):
    with open(path, "w") as f:
        for title, seq in records:
            f.write(f">{title}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + "\n")


def nucleotide_records():
    rng = random.Random(4004)

    def rand(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    repeat = rand(120)
    return [
        ("m01 microsatellite and a repeat", rand(200) + "AT" * 40 + rand(150) + repeat + rand(90)),
        ("m02 poly-A at the start", "A" * 90 + rand(260)),
        ("m03 short, no mask", rand(40)),
        ("m04 CAG repeat at the end", rand(300) + "CAG" * 30),
        ("m05 one base", "A"),
        ("m06 two low-complexity runs", rand(150) + "GA" * 35 + rand(400) + "TTTTC" * 14 + rand(120)),
        ("m07 the repeat twice", repeat + rand(80) + repeat),
        ("m08 ambiguity in a run", rand(100) + "ACAC" * 10 + "N" * 5 + "ACAC" * 10 + rand(100)),
    ]


def protein_records():
    rng = random.Random(5005)
    aa = "ACDEFGHIKLMNPQRSTVWY"

    def rand(n):
        return "".join(rng.choice(aa) for _ in range(n))

    return [
        ("p01 a low-complexity stretch", rand(80) + "QQQQQQQQQQQQQQQQPPPPQQQQ" + rand(60)),
        ("p02 no mask", rand(70)),
        ("p03 two stretches", "SSSSGSSSSGSSSSGSSSS" + rand(90) + "EEEEKKKEEEEKKKEEEE" + rand(20)),
        ("p04 short", rand(12)),
        ("p05 stretch at the end", rand(100) + "AAAAAAAAAAAAAAAAAAAAAAA"),
    ]


def hand_mask(path, program, options, entries):
    """Blast-db-mask-info in ASN.1 text: `entries` are (local id, [(from, to)]) with inclusive ends."""
    lines = ["Blast-db-mask-info ::= {", "  algo-id 0,", f"  algo-program {program},",
             f"  algo-options \"{options}\",", "  masks {", "    masks {"]
    locs = []
    for local_id, ranges in entries:
        ints = ",\n".join(f"          {{\n            from {a},\n            to {b},\n"
                          f"            id local id {local_id}\n          }}" for a, b in ranges)
        locs.append("      packed-int {\n" + ints + "\n      }")
    lines.append(",\n".join(locs))
    lines += ["    },", "    more FALSE", "  }", "}", ""]
    Path(path).write_text("\n".join(lines))


def run(cmd, cwd, check=True):
    result = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True)
    if check and result.returncode != 0:
        raise SystemExit(f"{' '.join(map(str, cmd))}\n{result.stdout!r}{result.stderr!r}")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi-bin", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parent))
    args = ap.parse_args()
    bin_dir = Path(args.ncbi_bin)
    work = Path(args.work)
    out = Path(args.out)
    if work.exists():
        raise SystemExit(f"{work} exists; give a new work directory")
    work.mkdir(parents=True)

    write_fasta(work / "mask_nucl.fna", nucleotide_records())
    write_fasta(work / "mask_prot.faa", protein_records())
    run([bin_dir / "dustmasker", "-in", "mask_nucl.fna", "-outfmt", "maskinfo_asn1_bin",
         "-out", "dust.asnb"], work)
    run([bin_dir / "segmasker", "-in", "mask_prot.faa", "-outfmt", "maskinfo_asn1_bin",
         "-out", "seg.asnb"], work)
    # Masks written by hand under the windowmasker program: a range at the start, one that ends
    # at the last base, a one-base range, two ranges in one sequence that touch, and the one-base
    # sequence m05 masked whole.
    hand_mask(work / "wm.asnt", "windowmasker", "", [
        (1, [(0, 9), (300, 309), (500, 529)]),
        (2, [(349, 349)]),
        (4, [(10, 19), (20, 29), (380, 389)]),
        (5, [(0, 0)]),
        (7, [(0, 119), (200, 319)]),
    ])
    # The same masks for a second database, with other options (another algorithm id).
    hand_mask(work / "wm_opts.asnt", "windowmasker", "-t_extend 3", [(3, [(5, 14)]), (6, [(0, 99)])])

    # A ':' in the options (stored as 0x01 in the description) for the two ways a description is built.
    hand_mask(work / "wm_colon.asnt", "windowmasker", "-x a:b", [(1, [(0, 9), (300, 309)]), (4, [(10, 19)])])

    def makedb(name, fasta, dbtype, version, *extra):
        result = run([bin_dir / "makeblastdb", "-in", fasta, "-dbtype", dbtype, "-out", name,
                      "-blastdb_version", str(version), *extra], work)
        if result.stderr:
            raise SystemExit(f"makeblastdb {name}: {result.stderr!r}")

    makedb("mask_nucl_v4", "mask_nucl.fna", "nucl", 4, "-mask_data", "dust.asnb,wm.asnt",
           "-title", "LOSAT mask fixture, nucleotide, version 4")
    makedb("mask_nucl_v5", "mask_nucl.fna", "nucl", 5, "-mask_data", "dust.asnb,wm.asnt")
    makedb("mask_nucl_mv", "mask_nucl.fna", "nucl", 4, "-mask_data", "dust.asnb,wm.asnt",
           "-max_file_sz", "400B", "-title", "masked nucleotide in volumes")
    # -mask_id: user algorithms (ids 100, 101); the ':' of a description is stored as 0x01.
    makedb("mask_nucl_ids", "mask_nucl.fna", "nucl", 4, "-mask_data", "dust.asnb,wm.asnt",
           "-mask_id", "dustX,wmY", "-mask_desc", "dust: from dustmasker,hand written")
    # A second database whose only algorithm, "wmZ", also has the id 100 (remapped in an alias),
    # and one whose "wmY" has the id 100 (the id of "wmY" in mask_nucl_ids is 101).
    makedb("mask_nucl_z", "mask_nucl.fna", "nucl", 4, "-mask_data", "wm_opts.asnt", "-mask_id", "wmZ")
    makedb("mask_nucl_y", "mask_nucl.fna", "nucl", 4, "-mask_data", "wm_opts.asnt", "-mask_id", "wmY")
    makedb("mask_nucl_opts", "mask_nucl.fna", "nucl", 4, "-mask_data", "wm_opts.asnt,wm.asnt")
    makedb("mask_nucl_colon", "mask_nucl.fna", "nucl", 4, "-mask_data", "wm_colon.asnt")
    makedb("mask_nucl_colon_id", "mask_nucl.fna", "nucl", 4, "-mask_data", "wm_colon.asnt",
           "-mask_id", "c:id")
    makedb("mask_prot_v4", "mask_prot.faa", "prot", 4, "-mask_data", "seg.asnb")
    makedb("mask_nucl_plain", "mask_nucl.fna", "nucl", 4)
    (work / "mask_remap.nal").write_bytes(
        b"TITLE algorithms of four volumes\n"
        b"DBLIST mask_nucl_ids mask_nucl_plain mask_nucl_z mask_nucl_y mask_nucl_v4\n")

    dbs = [("mask_nucl_v4", "nucl"), ("mask_nucl_v5", "nucl"), ("mask_nucl_mv", "nucl"),
           ("mask_nucl_ids", "nucl"), ("mask_nucl_z", "nucl"), ("mask_nucl_y", "nucl"),
           ("mask_nucl_opts", "nucl"), ("mask_prot_v4", "prot"), ("mask_remap", "nucl"),
           ("mask_nucl_colon", "nucl"), ("mask_nucl_colon_id", "nucl")]
    expected = out / "expected"
    db_out = out / "db"
    for name, dbtype in dbs:
        info = run([bin_dir / "blastdbcmd", "-db", name, "-dbtype", dbtype, "-info"], work).stdout.decode()
        m = re.search(r"\nAvailable filtering algorithms applied to database sequences:\n\n.*?\n(?=\nVolumes:)",
                      info, re.S)
        if not m:
            raise SystemExit(f"no algorithm table for {name}:\n{info}")
        # GetAvailableMaskAlgorithmDescriptions(): from the leading newline to the last row's newline.
        (expected / f"{name}.algorithms.txt").write_text(m.group(0))
        ids = [int(line.split()[0]) for line in m.group(0).splitlines()[4:] if line.strip()]
        rows = []
        for algo in ids:
            masks = run([bin_dir / "blastdbcmd", "-db", name, "-dbtype", dbtype, "-entry", "all",
                         "-outfmt", f"%o\t%m{algo}"], work).stdout.decode()
            rows += [f"{algo}\t{line}\n" for line in masks.splitlines()]
        (expected / f"{name}.masks.tsv").write_text("".join(rows))

    for path in sorted(work.iterdir()):
        if not path.name.startswith("mask_") or path.suffix in (".fna", ".faa"):
            continue
        if path.suffix in NUCL_EXTS + PROT_EXTS + [".nal", ".ndb"]:
            shutil.copyfile(path, db_out / path.name)
    shutil.copyfile(work / "mask_nucl.fna", out / "mask_nucl.fna")
    shutil.copyfile(work / "mask_prot.faa", out / "mask_prot.faa")
    for f in ("wm.asnt", "wm_opts.asnt", "wm_colon.asnt", "dust.asnb", "seg.asnb"):
        shutil.copyfile(work / f, out / f"mask_{f}")


if __name__ == "__main__":
    main()
