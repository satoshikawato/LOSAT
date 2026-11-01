#!/usr/bin/env python3
"""Builds tests/db_search/multidb/: the databases of the rows db.mdb.* (DI-8f: -db with several names).

Usage: make_multidb.py MAKEBLASTDB   (writes multidb/*.fa, multidb/*.map and the databases below)

The records are those of taxfilter/tfp.faa (23 protein globins) and taxfilter/tfn.fna (16 globin cDNAs); the definition
lines (some records carry a second one) stay as they are, the taxids come from taxfilter/tfp.map and tfn.map.

- `pa` / `pb`: version 5 protein, `-parse_seqids`, `-taxid_map`; `pa` = records 1-12 of tfp.faa, `pb` = records 13-23 and
  `Q00001.1`, an exact copy of record 2's sequence (taxid 9606), so that two databases hold an identical sequence.
- `na` / `nb`: version 5 nucleotide, `-parse_seqids`, `-taxid_map`; `na` = records 1-8 of tfn.fna, `nb` = records 9-16 and
  `ENSTTEST0001.1`, a copy of record 1's sequence (taxid 9606).
- `nnp`: version 5 nucleotide without `-parse_seqids` (records of nb).
- `n4a` / `n4b`: version 4 nucleotide without `-parse_seqids` (records of na / nb).
- `p4a` / `p4b`: version 4 protein without `-parse_seqids` (records of pa / pb).
- `d1/g` and `d2/g`: version 5 `-parse_seqids` nucleotide of the records of na and nb, one base name in two directories.
- `et`: version 5 `-parse_seqids` nucleotide of the records of na with an empty title.
- `al.nal`: an alias file (TITLE multidb alias of na, DBLIST na).
- `sl.txt`: a seqid list of the first accession of na, nb, pa and pb.

`.pjs` and `.njs` (the metadata with the build's paths and dates; no search reads them) are dropped.
"""
import os
import re
import subprocess
import sys

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "multidb")
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "taxfilter")


def makeblastdb(binary, *args):
    subprocess.run([binary, *args], check=True, cwd=HERE, stdout=subprocess.DEVNULL)


def read_fasta(name):
    records = []
    with open(os.path.join(SRC, name)) as f:
        for line in f:
            if line.startswith(">"):
                records.append([line, []])
            else:
                records[-1][1].append(line)
    return records


def read_map(name):
    with open(os.path.join(SRC, name)) as f:
        return dict(line.split() for line in f if line.strip())


def accessions(header):
    """The first token of each definition line (control-A separated); `sp|P69905.2|` gives P69905.2."""
    found = []
    for part in header[1:].rstrip("\n").split("\x01"):
        token = part.split(" ", 1)[0]
        found.append(token.split("|")[1] if "|" in token else token)
    return found


def write_set(base, records, taxids):
    """Writes base.fa and, when `taxids` is not None, base.map (the map lines of the accessions in the records)."""
    with open(os.path.join(HERE, base + ".fa"), "w") as f:
        for header, lines in records:
            f.write(header)
            f.writelines(lines)
    if taxids is not None:
        with open(os.path.join(HERE, base + ".map"), "w") as f:
            for header, _ in records:
                for acc in accessions(header):
                    if acc in taxids:
                        f.write(f"{acc} {taxids[acc]}\n")


def main(binary):
    os.makedirs(os.path.join(HERE, "d1"), exist_ok=True)
    os.makedirs(os.path.join(HERE, "d2"), exist_ok=True)
    tfp, tfn = read_fasta("tfp.faa"), read_fasta("tfn.fna")
    pmap, nmap = read_map("tfp.map"), read_map("tfn.map")
    pmap["Q00001.1"] = "9606"
    nmap["ENSTTEST0001.1"] = "9606"
    pcopy = [">sp|Q00001.1| copy of record 2 [tie]\n", tfp[1][1]]
    ncopy = [">ENSTTEST0001.1 copy of record 1\n", tfn[0][1]]
    sets = {
        "pa": (tfp[:12], pmap), "pb": (tfp[12:] + [pcopy], pmap),
        "na": (tfn[:8], nmap), "nb": (tfn[8:] + [ncopy], nmap),
    }
    for base, (records, taxids) in sets.items():
        write_set(base, records, taxids)
    write_set("nnp", sets["nb"][0], None)
    for base, src in (("n4a", "na"), ("n4b", "nb"), ("p4a", "pa"), ("p4b", "pb"), ("d1/g", "na"), ("d2/g", "nb"),
                      ("et", "na")):
        write_set(base, sets[src][0], None)
    for base, kind in (("pa", "prot"), ("pb", "prot"), ("na", "nucl"), ("nb", "nucl")):
        makeblastdb(binary, "-in", base + ".fa", "-dbtype", kind, "-out", base, "-parse_seqids", "-blastdb_version", "5",
                    "-taxid_map", base + ".map", "-title", "multidb " + base)
    makeblastdb(binary, "-in", "nnp.fa", "-dbtype", "nucl", "-out", "nnp", "-blastdb_version", "5", "-title", "multidb nnp")
    for base, kind in (("n4a", "nucl"), ("n4b", "nucl"), ("p4a", "prot"), ("p4b", "prot")):
        makeblastdb(binary, "-in", base + ".fa", "-dbtype", kind, "-out", base, "-blastdb_version", "4",
                    "-title", "multidb " + base)
    for d in ("d1", "d2"):
        makeblastdb(binary, "-in", f"{d}/g.fa", "-dbtype", "nucl", "-out", f"{d}/g", "-parse_seqids",
                    "-blastdb_version", "5", "-title", f"multidb {d} g")
    makeblastdb(binary, "-in", "et.fa", "-dbtype", "nucl", "-out", "et", "-parse_seqids", "-blastdb_version", "5",
                "-title", "")
    for name in ("pa.pjs", "pb.pjs", "p4a.pjs", "p4b.pjs", "na.njs", "nb.njs", "nnp.njs", "n4a.njs", "n4b.njs",
                 "d1/g.njs", "d2/g.njs", "et.njs"):
        path = os.path.join(HERE, name)
        if os.path.exists(path):
            os.remove(path)
    with open(os.path.join(HERE, "al.nal"), "w") as f:
        f.write("TITLE multidb alias of na\nDBLIST na\n")
    firsts = [accessions(sets[b][0][0][0])[0] for b in ("na", "nb", "pa", "pb")]
    with open(os.path.join(HERE, "sl.txt"), "w") as f:
        f.write("".join(acc + "\n" for acc in firsts))


if __name__ == "__main__":
    main(sys.argv[1])
