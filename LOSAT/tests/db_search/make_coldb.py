#!/usr/bin/env python3
"""Builds tests/db_search/coldb/: the databases with column files of the rows db.col.* (line SW, SW8 follow-up 2).

Usage: make_coldb.py NCBI_BIN [OUT]   (NCBI_BIN holds NCBI BLAST+ 2.17.0's makeblastdb, segmasker, dustmasker and
windowmasker; OUT defaults to coldb/ next to this script)

A column of a volume is the files `.?Xa` (index), `.?Xb` (big-endian data) and `.?Xc` (little-endian data) for a letter X;
`makeblastdb -mask_data` writes the masks as the column `BlastDb/MaskData` (`.paa .pab .pac`, `.naa .nab .nac`). NCBI's
search programs open a volume's columns only through `CSeqDBVol::x_OpenAllColumns`, which runs for a masking algorithm
of `-db_soft_mask` or `-db_hard_mask` alone, so a search without these options reads none of these files.

The records and taxids are those of taxfilter/tfp.faa (23 protein globins), tfp.map, taxfilter/tfn.fna (16 globin cDNAs)
and tfn.map:

- `p5m`: version 5 protein, `-parse_seqids`, `-taxid_map`, the masks of `segmasker -parse_seqids` (algorithm 21, seg).
- `p4m`: version 4 protein without `-parse_seqids`, the masks of `segmasker` (which, like makeblastdb without
  `-parse_seqids`, names the records by their number).
- `n5m`: version 5 nucleotide, `-parse_seqids`, `-taxid_map`, two volumes (`-max_file_sz 4KB`: `n5m.00` and `n5m.01`,
  each with its column), the masks of `dustmasker -parse_seqids` and `windowmasker -parse_seqids` (two algorithms).
- `n4m`: version 4 nucleotide without `-parse_seqids`, the masks of `dustmasker`.
- `dmg/CASE/p4m` (every case of DAMAGE) and `dmg/CASE/n4m` (the cases of NUCL_DAMAGE): copies of `p4m` and `n4m` whose
  column is damaged as DAMAGE says; NCBI searches each as the undamaged database.

`.pjs` and `.njs` (the build's paths and dates; no search reads them) are dropped. A rebuild changes the dates stored in
the index and column files.
"""
import os
import shutil
import subprocess
import sys

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "taxfilter")

# Each case: what is done to every column of the volume (index `.?Xa`, data `.?Xb` and `.?Xc`), and the error that
# NCBI's column reader would raise if a search opened the column (CSeqDBVol::x_OpenAllColumns, CSeqDBColumn).
DAMAGE = {
    "idx_cut": "the index cut to half its length",
    "idx_junk": "every byte of the index replaced by 0xA5",
    "idx_empty": "the index empty",
    "idx_ver": "the index's format version (its first 4 bytes, big-endian) 2 instead of 1",
    "data_cut": "both data files cut to half their length",
    "data_junk": "every byte of both data files replaced by 0xA5",
    "be_only": "the little-endian data file removed (the big-endian one alone)",
    "le_only": "the big-endian data file removed (the little-endian one alone)",
    "no_data": "both data files removed (the index alone)",
    "no_idx": "the index removed (the data files alone)",
    "oids": "the column of a database of the first 6 records (6 OIDs, not the volume's)",
    "dup": "the column copied under the next letter (two columns of one title)",
}
NUCL_DAMAGE = ("idx_junk", "data_cut", "be_only")


def run(cwd, binary, *args):
    subprocess.run([binary, *args], check=True, cwd=cwd, stdout=subprocess.DEVNULL)


def read_fasta(path):
    records = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                records.append([line])
            else:
                records[-1].append(line)
    return records


def columns(directory, base, kind):
    """The letters X of the columns of the volume `base` (`.?Xa` present)."""
    return [x for x in "abcdefghijklmnopqrstuvwxyz" if os.path.exists(os.path.join(directory, f"{base}.{kind}{x}a"))]


def volumes(directory, base, kind):
    """The volume names of `base` (`base`, or `base.00`, `base.01`, ... when the database has several volumes)."""
    names = sorted(n[: -len(f".{kind}in")] for n in os.listdir(directory) if n.endswith(f".{kind}in"))
    return [n for n in names if n == base or n.startswith(base + ".")]


def damage(directory, base, kind, case, oids_column=None):
    """Applies DAMAGE[case] to every column of every volume of `directory/base` (kind `p` or `n`). `oids_column` is
    the path without extension of a volume whose columns `oids` copies in."""
    for volume in volumes(directory, base, kind):
        stem = os.path.join(directory, volume)
        letters = columns(directory, volume, kind)
        for x in letters:
            index, big, little = (f"{stem}.{kind}{x}{e}" for e in "abc")
            if case == "idx_cut":
                truncate(index)
            elif case == "idx_junk":
                junk(index)
            elif case == "idx_empty":
                open(index, "wb").close()
            elif case == "idx_ver":
                with open(index, "r+b") as f:
                    f.write((2).to_bytes(4, "big"))
            elif case == "data_cut":
                truncate(big)
                truncate(little)
            elif case == "data_junk":
                junk(big)
                junk(little)
            elif case == "be_only":
                os.remove(little)
            elif case == "le_only":
                os.remove(big)
            elif case == "no_data":
                os.remove(big)
                os.remove(little)
            elif case == "no_idx":
                os.remove(index)
            elif case == "oids":
                for e in "abc":
                    shutil.copyfile(f"{oids_column}.{kind}{x}{e}", f"{stem}.{kind}{x}{e}")
            elif case == "dup":
                free = next(y for y in "abcdefghijklmnopqrstuvwxyz" if y not in letters)
                for e in "abc":
                    if os.path.exists(f"{stem}.{kind}{x}{e}"):
                        shutil.copyfile(f"{stem}.{kind}{x}{e}", f"{stem}.{kind}{free}{e}")
                break
            else:
                raise SystemExit(f"unknown case {case}")


def truncate(path):
    size = os.path.getsize(path)
    with open(path, "r+b") as f:
        f.truncate(size // 2)


def junk(path):
    size = os.path.getsize(path)
    with open(path, "wb") as f:
        f.write(b"\xa5" * size)


def copy_db(src_dir, base, dst_dir):
    os.makedirs(dst_dir, exist_ok=True)
    for name in os.listdir(src_dir):
        if name.startswith(base + "."):
            shutil.copyfile(os.path.join(src_dir, name), os.path.join(dst_dir, name))


def build(bin_dir, out):
    tool = lambda name: os.path.join(bin_dir, name)  # noqa: E731
    os.makedirs(out, exist_ok=True)
    work = os.path.join(out, "work")
    os.makedirs(work, exist_ok=True)
    tfp, tfn = os.path.join(SRC, "tfp.faa"), os.path.join(SRC, "tfn.fna")
    pmap, nmap = os.path.join(SRC, "tfp.map"), os.path.join(SRC, "tfn.map")

    # The masks.
    run(work, tool("segmasker"), "-in", tfp, "-parse_seqids", "-outfmt", "maskinfo_asn1_bin", "-out", "seg_ids.asnb")
    run(work, tool("segmasker"), "-in", tfp, "-outfmt", "maskinfo_asn1_bin", "-out", "seg.asnb")
    run(work, tool("dustmasker"), "-in", tfn, "-parse_seqids", "-outfmt", "maskinfo_asn1_bin", "-out", "dust_ids.asnb")
    run(work, tool("dustmasker"), "-in", tfn, "-outfmt", "maskinfo_asn1_bin", "-out", "dust.asnb")
    run(work, tool("windowmasker"), "-in", tfn, "-mk_counts", "-parse_seqids", "-out", "wm.counts")
    run(work, tool("windowmasker"), "-in", tfn, "-ustat", "wm.counts", "-parse_seqids", "-outfmt", "maskinfo_asn1_bin",
        "-out", "wm_ids.asnb")
    w = lambda name: os.path.join(work, name)  # noqa: E731

    # The databases.
    run(out, tool("makeblastdb"), "-in", tfp, "-dbtype", "prot", "-out", "p5m", "-parse_seqids", "-taxid_map", pmap,
        "-blastdb_version", "5", "-mask_data", w("seg_ids.asnb"), "-title", "coldb p5m")
    run(out, tool("makeblastdb"), "-in", tfp, "-dbtype", "prot", "-out", "p4m", "-blastdb_version", "4",
        "-mask_data", w("seg.asnb"), "-title", "coldb p4m")
    run(out, tool("makeblastdb"), "-in", tfn, "-dbtype", "nucl", "-out", "n5m", "-parse_seqids", "-taxid_map", nmap,
        "-blastdb_version", "5", "-max_file_sz", "4KB", "-mask_data", w("dust_ids.asnb") + "," + w("wm_ids.asnb"),
        "-title", "coldb n5m")
    run(out, tool("makeblastdb"), "-in", tfn, "-dbtype", "nucl", "-out", "n4m", "-blastdb_version", "4",
        "-mask_data", w("dust.asnb"), "-title", "coldb n4m")

    # The columns of 6 OIDs for the case `oids` (record 6 of tfp.faa has seg masks).
    for name, records in (("p6.faa", read_fasta(tfp)[:6]), ("n6.fna", read_fasta(tfn)[:6])):
        with open(w(name), "w") as f:
            f.writelines(line for record in records for line in record)
    run(work, tool("segmasker"), "-in", "p6.faa", "-outfmt", "maskinfo_asn1_bin", "-out", "seg6.asnb")
    run(work, tool("dustmasker"), "-in", "n6.fna", "-outfmt", "maskinfo_asn1_bin", "-out", "dust6.asnb")
    run(work, tool("makeblastdb"), "-in", "p6.faa", "-dbtype", "prot", "-out", "p6", "-blastdb_version", "4",
        "-mask_data", "seg6.asnb")
    run(work, tool("makeblastdb"), "-in", "n6.fna", "-dbtype", "nucl", "-out", "n6", "-blastdb_version", "4",
        "-mask_data", "dust6.asnb")

    for small, kind in (("p6", "p"), ("n6", "n")):
        if not columns(work, small, kind):
            raise SystemExit(f"{small} has no column")
    for base, kind, cases, small in (("p4m", "p", DAMAGE, "p6"), ("n4m", "n", NUCL_DAMAGE, "n6")):
        for case in cases:
            target = os.path.join(out, "dmg", case)
            copy_db(out, base, target)
            damage(target, base, kind, case, oids_column=w(small))

    for directory, _, files in os.walk(out):
        for name in files:
            if name.endswith((".pjs", ".njs")):
                os.remove(os.path.join(directory, name))
    shutil.rmtree(work)


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        raise SystemExit(__doc__)
    default = os.path.join(os.path.dirname(os.path.abspath(__file__)), "coldb")
    build(sys.argv[1], sys.argv[2] if len(sys.argv) == 3 else default)
