#!/usr/bin/env python3
"""Write the inputs of the `db.sw.*` rows of line PS-SW's fix round 1 (PSIBLAST `-use_sw_tback`, 2026-11-27).

NCBI BLAST+ 2.17.0 makeblastdb is the comparison oracle only. From tks_W (`tests/fasta/rps/tks_w.faa`), tks_M
(`tests/fasta/rps/tks_m.faa`) and sfail (`tests/fasta/x4/sw_subjects.faa`: tks_M with D first, twice), and the tk PSSM
`tests/fixtures/pssm/ps4_tk_paja.asn.gz`, it writes:
- `tests/fasta/pssw/sfail.faa`: sfail alone, a subject whose Smith-Waterman redo fails (two HSPs, and the query's first
  letter scores D at 0 or less), so that a search whose last redone match is sfail's ends in `Traceback failed`.
- `tests/fasta/pssw/smm.faa`: tks_M twice, a subject whose redo never ends in NCBI (two HSPs, M against the query's first
  letter above 0).
- the version 5 databases (`-parse_seqids`) `tests/db_search/pssw/psw` (sx_frag, whose only line has a ti, then sfail)
  and `pswk` (sexact, tks_W itself, then sx_frag and sfail), searched with `tests/db_search/x4excl/neg.txt`: the negative
  list keeps sx_frag without a line, and NCBI's fetch of the redo excludes it, which sets the redo's status back to 0.
  The redo runs the matches in OID order at `-comp_based_stats 0` and in e-value order at 1 to 3, so the excluded
  subject comes before sfail's failing match in the one and after it in the other.
- `tests/fixtures/pssm/ps_sw_tk_scores.asn.gz`: the tk PSSM without the Pssm member `intermediateData` (scores only, no
  frequency ratios), gzipped without a name or time.
The databases are written with the creation date in their index files, so the frozen outputs need the committed files,
not a new run of this script.

With `--fix2-only` (fix round 2, 2026-11-28) it writes only the inputs of that round, from the proteins of
`tests/db_search/t5/hd.faa` (hd0, which `tests/db_search/t5/hd_q0.faa` is; hd1; hd2, a chance match of it), tks_W, tks_M
and sfail, and makes the version 5 databases with `-parse_seqids`:
- `tests/fixtures/pssm/ps_sw_hd0.asn.gz`: the PSSM of hd_q0 from NCBI psiblast (`-subject hd.faa -num_iterations 1
  -save_pssm_after_last_round`), gzipped without a name or time.
- `tests/db_search/pssw/pfu8` (REVIEW_fix1 MF-1's input): OID 0 hd2 (`NSW1.1 weak`) with its title a UTF8String, OIDs 1
  to 10 ten copies of hd0 (`NSS<k>.1`), OID 11 hd1 (`NSW2.1`, which `tests/db_search/t5/fx_neg2.txt` removes).
- `tests/db_search/pssw/pfu2t`: OID 0 hd2 (`NSW1.1 weak`) with its title's tag an OctetString, OID 1 hd0 (`NSS1.1`) with
  its title's tag an INTEGER, OID 2 hd1 (`NSW2.1`): NCBI's exception names the tag of the header that its redo reads
  first.
- `tests/db_search/pssw/pende`: OID 0 smm (tks_M twice, whose redo never ends in NCBI), OID 1 stk (YP_009220601.1 of
  `tests/fasta/WSSV.faa`, the query of the tk PSSM, whose match comes before smm's by e-value) with its title's tag an
  OctetString.
- `tests/fasta/pssw/sflip.faa`: sfail, then sb (tks_W with H8N and V63N): in round 1 of tks_W sb's match is redone last,
  in round 2 (the PSSM of that round) sfail's, which fails after sb's was kept.

usage: make_pssw.py --ncbi-bin DIR [--fix2-only]"""
import gzip
import importlib.util
import os
import subprocess
import sys

args = sys.argv[1:]
if len(args) not in (2, 3) or args[0] != "--ncbi-bin" or args[2:] not in ([], ["--fix2-only"]):
    sys.exit(__doc__)
NCBI = args[1]
FIX2 = args[2:] == ["--fix2-only"]
D = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(os.path.dirname(D))


def fasta(path):
    """The records of `path` as {first word of the title: residues}."""
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line
    return seqs


tks_w = fasta(f"{TESTS}/fasta/rps/tks_w.faa")["tks_W"]
tks_m = fasta(f"{TESTS}/fasta/rps/tks_m.faa")["tks_M"]
sfail = fasta(f"{TESTS}/fasta/x4/sw_subjects.faa")["sfail"]
assert sfail == "D" + tks_m[1:] + "D" + tks_m[1:], "sfail is tks_M with D first, twice"
assert "MAGRVE" in tks_w
ENV = {"PATH": "/usr/bin:/bin", "HOME": D, "BLAST_USAGE_REPORT": "false"}


def fix2():
    """The inputs of fix round 2 (`--fix2-only`)."""
    spec = importlib.util.spec_from_file_location("make_t5", f"{TESTS}/db_search/t5/make_t5.py")
    t5 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t5)
    hd = fasta(f"{TESTS}/db_search/t5/hd.faa")
    assert hd["hd0"] == fasta(f"{TESTS}/db_search/t5/hd_q0.faa")["q0"]
    stk = fasta(f"{TESTS}/fasta/WSSV.faa")["YP_009220601.1"]
    assert stk[1:].startswith(tks_w[1:]) and len(stk) == 398, "stk is the tk PSSM's query, tks_W its start"

    def tag(new):
        return lambda record: t5.with_first_tag(record, b"\xa0\x80", 0x1A, new)

    dbs = {
        "pfu8": ([("NSW1.1", "weak", hd["hd2"])] + [(f"NSS{k}.1", f"strong copy {k}", hd["hd0"]) for k in range(1, 11)]
                 + [("NSW2.1", "other", hd["hd1"])], {0: tag(0x0C)}),
        "pfu2t": ([("NSW1.1", "weak", hd["hd2"]), ("NSS1.1", "strong", hd["hd0"]), ("NSW2.1", "other", hd["hd1"])],
                  {0: tag(0x04), 1: tag(0x02)}),
        "pende": ([("dbj|BAA200501.1|", "smm tks_M twice", tks_m + tks_m),
                   ("dbj|BAA200502.1|", "stk the tk PSSM's query", stk)], {1: tag(0x04)}),
    }
    for name, (records, edits) in dbs.items():
        fa = f"{D}/{name}.fa"
        with open(fa, "w") as f:
            f.write("".join(f">{i} {t}\n{r}\n" for i, t, r in records))
        subprocess.run([f"{NCBI}/makeblastdb", "-in", fa, "-dbtype", "prot", "-parse_seqids", "-out", f"{D}/{name}",
                        "-title", name, "-blastdb_version", "5"], check=True, stdout=subprocess.DEVNULL, env=ENV)
        os.remove(fa)
        for oid, edit in edits.items():
            t5.patch_record_any(f"{D}/{name}", "p", oid, edit)

    pssm = f"{D}/hd0.asn"
    subprocess.run([f"{NCBI}/psiblast", "-query", f"{TESTS}/db_search/t5/hd_q0.faa", "-subject",
                    f"{TESTS}/db_search/t5/hd.faa", "-num_iterations", "1", "-out_pssm", pssm,
                    "-save_pssm_after_last_round", "-outfmt", "6"], check=True, stdout=subprocess.DEVNULL, env=ENV)
    with open(f"{TESTS}/fixtures/pssm/ps_sw_hd0.asn.gz", "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as f:
            f.write(open(pssm, "rb").read())
    os.remove(pssm)

    sb = tks_w[:7] + "N" + tks_w[8:62] + "N" + tks_w[63:]
    assert (tks_w[7], tks_w[62]) == ("H", "V") and len(sb) == len(tks_w)
    with open(f"{TESTS}/fasta/pssw/sflip.faa", "w") as f:
        f.write(f">sfail tks_M with D first, twice\n{sfail}\n>sb tks_W with H8N and V63N\n{sb}\n")


if FIX2:
    fix2()
    sys.exit(0)

os.makedirs(f"{TESTS}/fasta/pssw", exist_ok=True)
with open(f"{TESTS}/fasta/pssw/sfail.faa", "w") as f:
    f.write(f">sfail tks_M with D first, twice\n{sfail}\n")
with open(f"{TESTS}/fasta/pssw/smm.faa", "w") as f:
    f.write(f">smm tks_M twice\n{tks_m}{tks_m}\n")

dbs = {
    "psw": [
        ("gnl|ti|7204", "sx_frag six letters of tks_W", "MAGRVD"),
        ("dbj|BAA200302.1|", "sfail tks_M with D first, twice", sfail),
    ],
    "pswk": [
        ("dbj|BAA200401.1|", "sexact tks_W itself", tks_w),
        ("gnl|ti|7304", "sx_frag six letters of tks_W", "MAGRVD"),
        ("dbj|BAA200403.1|", "sfail tks_M with D first, twice", sfail),
    ],
}
for name, records in dbs.items():
    fa = f"{D}/{name}.fa"
    with open(fa, "w") as f:
        f.write("".join(f">{i} {t}\n{r}\n" for i, t, r in records))
    subprocess.run([f"{NCBI}/makeblastdb", "-in", fa, "-dbtype", "prot", "-parse_seqids", "-out", f"{D}/{name}",
                    "-title", name, "-blastdb_version", "5"], check=True, stdout=subprocess.DEVNULL)
    os.remove(fa)

tk = gzip.open(f"{TESTS}/fixtures/pssm/ps4_tk_paja.asn.gz", "rt").read()
start, end = tk.index("    intermediateData {\n"), tk.index("    finalData {\n")
with open(f"{TESTS}/fixtures/pssm/ps_sw_tk_scores.asn.gz", "wb") as raw:
    with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as f:
        f.write((tk[:start] + tk[end:]).encode())
