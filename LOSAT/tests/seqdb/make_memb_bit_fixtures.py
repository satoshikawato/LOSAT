#!/usr/bin/env python3
"""Build and freeze the MEMB_BIT and MASKLIST interaction fixtures (LOSAT/tests/seqdb/memb_bit).

Uses official NCBI BLAST+ 2.17.0 binaries as comparison oracle only: the Rust integration test
(`LOSAT/tests/memb_bit.rs`) reads the frozen files and never runs any NCBI binary.

Usage:
    python3 make_memb_bit_fixtures.py --ncbi-bin DIR --work DIR [--out LOSAT/tests/seqdb/memb_bit]
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


# --- Minimal ASN.1 BER encoder for Blast-def-line-set (ASN1_BlastDefLine user object) ---
def ber_len(n: int) -> bytes:
    if n < 0x80:
        return bytes([n])
    b = n.to_bytes((n.bit_length() + 7) // 8, "big")
    return bytes([0x80 | len(b)]) + b


def ber_int(val: int) -> bytes:
    if val == 0:
        b = b"\x00"
    else:
        nbytes = (val.bit_length() + 8) // 8
        b = val.to_bytes(nbytes, "big", signed=True)
        if val < 0:
            while len(b) > 1 and b[0] == 0xFF and (b[1] & 0x80):
                b = b[1:]
    return b"\x02" + ber_len(len(b)) + b


def ber_vis(s: str) -> bytes:
    b = s.encode("ascii")
    return b"\x1a" + ber_len(len(b)) + b


def ber_seq(items: list[bytes]) -> bytes:
    return b"\x30\x80" + b"".join(items) + b"\x00\x00"


def ber_ctx(tag: int, inner: bytes) -> bytes:
    return bytes([0xA0 | tag]) + b"\x80" + inner + b"\x00\x00"


def ber_seqid_gi(gi: int) -> bytes:
    return ber_ctx(11, ber_int(gi))


def ber_seqid_text(choice_tag: int, acc: str, ver: int = 1, name: str | None = None) -> bytes:
    parts = []
    if name is not None:
        parts.append(ber_ctx(0, ber_vis(name)))
    parts.append(ber_ctx(1, ber_vis(acc)))
    parts.append(ber_ctx(3, ber_int(ver)))
    return ber_ctx(choice_tag, ber_seq(parts))


def ber_seqid_pdb(mol: str, chain: str) -> bytes:
    return ber_ctx(14, ber_seq([ber_ctx(0, ber_vis(mol)), ber_ctx(3, ber_vis(chain))]))


def ber_defline(title: str, seqids: list[bytes], memberships: list[int] | None = None) -> bytes:
    parts = [ber_ctx(0, ber_vis(title)), ber_ctx(1, ber_seq(seqids))]
    if memberships is not None:
        parts.append(ber_ctx(3, ber_seq([ber_int(m) for m in memberships])))
    return ber_seq(parts)


def make_asn_bioseq(seq_id_text: str, mol: str, seq_str: str, deflines: list[bytes]) -> str:
    bdls_hex = ber_seq(deflines).hex().upper()
    mol_kw, code_kw = ("aa", "ncbieaa") if mol == "aa" else ("dna", "iupacna")
    return f"""seq {{
  id {{ {seq_id_text} }},
  descr {{
    user {{
      type str "ASN1_BlastDefLine",
      data {{
        {{
          label str "ASN1_BlastDefLine",
          data oss {{ '{bdls_hex}'H }}
        }}
      }}
    }}
  }},
  inst {{
    repr raw,
    mol {mol_kw},
    length {len(seq_str)},
    seq-data {code_kw} "{seq_str}"
  }}
}}"""


def wrap_bioseq_set(entries: list[str]) -> str:
    return "Bioseq-set ::= {\n  seq-set {\n" + ",\n".join(entries) + "\n  }\n}\n"


PSEQ0 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVK"
PSEQ1 = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY"
PSEQ2 = "GVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG"
PSEQ3 = "MVLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASE"
PSEQ4 = "MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITK"
PSEQ5 = "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASE"

# 80-bp high-complexity nucleotide sequences whose frame +1 translates to clean proteins
NSEQ0 = "ATGAAGACTGCATACATTGCTAAGCAAAGACAAATCTCCTTTGTGAAAAGTCACTTTAGTCGTCAATTGGAAGAACGTCT"
NSEQ1 = "GCTTGTGACGAATTTGGACACATCAAGCTCATGAACCCTCAACGTTCTACCGTTTGGTATGCTTGCGATGAATTCGGTCA"
NSEQ2 = "GGTGTGCTGAGCCCTGCTGACAAGACCAACGTGAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAGTATGGTGC"


def run_ncbi(ncbi_bin: Path, lock: Path | None, cmd: list[str], cwd: Path, env_extra: dict | None = None):
    env = {
        "PATH": f"{ncbi_bin}:/usr/bin:/bin",
        "BLAST_USAGE_REPORT": "false",
        "LC_ALL": "C",
        "HOME": str(cwd),
    }
    if env_extra:
        env.update(env_extra)
    full_cmd = ([f"flock", str(lock)] if lock else []) + [str(ncbi_bin / cmd[0])] + cmd[1:]
    proc = subprocess.run(full_cmd, cwd=cwd, env=env, capture_output=True, timeout=60)
    return proc.returncode, proc.stdout, proc.stderr


def build_databases(ncbi_bin: Path, lock: Path | None, work: Path):
    # 1. Protein ASN.1 with memberships covering bits 1, 2, 3, 31, 32, multi-int front(), empty, absent, zero
    prot_asn = wrap_bioseq_set([
        # OID 0: dl0 bit 1 (gi 1001, sp P00001.1), dl1 bit 2 (gi 1002, pdb 1ABC|A)
        make_asn_bioseq(
            'gi 1001, swissprot { name "M1_A", accession "P00001", version 1 }',
            "aa",
            PSEQ0,
            [
                ber_defline("seq0_dl0 bit1", [ber_seqid_gi(1001), ber_seqid_text(7, "P00001", 1, "M1_A")], [1]),
                ber_defline("seq0_dl1 bit2", [ber_seqid_gi(1002), ber_seqid_pdb("1ABC", "A")], [2]),
            ],
        ),
        # OID 1: dl0 bits 1&3 (5), dl1 multi-int {1, 2} (front() is 1, NOT 2!) + predicted XP_ for OID_MASK_TYPE 1
        make_asn_bioseq(
            'gi 1003, other { accession "NP_000001", version 1 }',
            "aa",
            PSEQ1,
            [
                ber_defline("seq1_dl0 bit1and3", [ber_seqid_gi(1003), ber_seqid_text(9, "NP_000001", 1)], [5]),
                ber_defline("seq1_dl1 bit1model_twoints", [ber_seqid_gi(1004), ber_seqid_text(9, "XP_000001", 1)], [1, 2]),
            ],
        ),
        # OID 2: no memberships field at all
        make_asn_bioseq(
            'gi 1005, other { accession "NP_000002", version 1 }',
            "aa",
            PSEQ2,
            [
                ber_defline("seq2_dl0 nomemb", [ber_seqid_gi(1005), ber_seqid_text(9, "NP_000002", 1)], None),
            ],
        ),
        # OID 3: empty memberships { }
        make_asn_bioseq(
            'gi 1006, other { accession "NP_000003", version 1 }',
            "aa",
            PSEQ3,
            [
                ber_defline("seq3_dl0 emptymemb", [ber_seqid_gi(1006), ber_seqid_text(9, "NP_000003", 1)], []),
            ],
        ),
        # OID 4: dl0 bit 32 (INT_MIN = -2147483648), dl1 bit 31 (0x40000000 = 1073741824)
        make_asn_bioseq(
            'gi 1007, other { accession "NP_000004", version 1 }',
            "aa",
            PSEQ4,
            [
                ber_defline("seq4_dl0 bit32", [ber_seqid_gi(1007), ber_seqid_text(9, "NP_000004", 1)], [-2147483648]),
                ber_defline("seq4_dl1 bit31", [ber_seqid_gi(1008), ber_seqid_text(9, "NP_000005", 1)], [1073741824]),
            ],
        ),
        # OID 5: memberships { 0 }
        make_asn_bioseq(
            'gi 1009, other { accession "NP_000006", version 1 }',
            "aa",
            PSEQ5,
            [
                ber_defline("seq5_dl0 membzero", [ber_seqid_gi(1009), ber_seqid_text(9, "NP_000006", 1)], [0]),
            ],
        ),
    ])

    prot_b_asn = wrap_bioseq_set([
        make_asn_bioseq(
            'gi 2010, other { accession "NP_000010", version 1 }',
            "aa",
            PSEQ0,
            [ber_defline("volb_0 bit2", [ber_seqid_gi(2010), ber_seqid_text(9, "NP_000010", 1)], [2])],
        ),
        make_asn_bioseq(
            'gi 2011, other { accession "NP_000011", version 1 }',
            "aa",
            PSEQ3,
            [ber_defline("volb_1 bit1", [ber_seqid_gi(2011), ber_seqid_text(9, "NP_000011", 1)], [1])],
        ),
    ])

    nucl_asn = wrap_bioseq_set([
        # OID 0: dl0 GI 2001 bit 1 (masked 5-26), dl1 GI 2002 bit 2 (unmasked)
        make_asn_bioseq(
            'gi 2001, other { accession "NM_000001", version 1 }',
            "dna",
            NSEQ0,
            [
                ber_defline("nseq0_dl0 bit1", [ber_seqid_gi(2001), ber_seqid_text(9, "NM_000001", 1)], [1]),
                ber_defline("nseq0_dl1 bit2", [ber_seqid_gi(2002), ber_seqid_text(9, "NM_000002", 1)], [2]),
            ],
        ),
        # OID 1: GI 2003 no memberships
        make_asn_bioseq(
            'gi 2003, other { accession "NM_000003", version 1 }',
            "dna",
            NSEQ1,
            [
                ber_defline("nseq1_dl0 nomemb", [ber_seqid_gi(2003), ber_seqid_text(9, "NM_000003", 1)], None),
            ],
        ),
        # OID 2: GI 2004 bit 2 (masked 30-51)
        make_asn_bioseq(
            'gi 2004, other { accession "NM_000004", version 1 }',
            "dna",
            NSEQ2,
            [
                ber_defline("nseq2_dl0 bit2", [ber_seqid_gi(2004), ber_seqid_text(9, "NM_000004", 1)], [2]),
            ],
        ),
    ])

    gi_mask_asnt = """Blast-db-mask-info ::= {
  algo-id 11,
  algo-program dust,
  algo-options "",
  masks {
    masks {
      int { from 5, to 25, id gi 2001 },
      int { from 30, to 50, id gi 2004 }
    },
    more FALSE
  }
}
"""

    db_dir = work / "db"
    db_dir.mkdir(parents=True)
    (db_dir / "prot.asnt").write_text(prot_asn)
    (db_dir / "prot_b.asnt").write_text(prot_b_asn)
    (db_dir / "nucl.asnt").write_text(nucl_asn)
    (db_dir / "gmask.asnt").write_text(gi_mask_asnt)

    for cmd in [
        ["makeblastdb", "-in", "prot.asnt", "-input_type", "asn1_txt", "-dbtype", "prot", "-parse_seqids", "-out", "pm5", "-title", "pm5"],
        ["makeblastdb", "-in", "prot.asnt", "-input_type", "asn1_txt", "-dbtype", "prot", "-parse_seqids", "-blastdb_version", "4", "-out", "pm4", "-title", "pm4"],
        ["makeblastdb", "-in", "prot_b.asnt", "-input_type", "asn1_txt", "-dbtype", "prot", "-parse_seqids", "-out", "pm5_b", "-title", "pm5_b"],
        ["makeblastdb", "-in", "nucl.asnt", "-input_type", "asn1_txt", "-dbtype", "nucl", "-parse_seqids", "-out", "nm5_gm", "-title", "nm5_gm", "-mask_data", "gmask.asnt", "-gi_mask", "-gi_mask_name", "gm_test"],
    ]:
        rc, out, err = run_ncbi(ncbi_bin, lock, cmd, db_dir)
        # Single-volume -gi_mask in NCBI 2.17.0 fails at metadata (mutually recursive alias); files are written.
        if rc != 0 and "-gi_mask" not in cmd:
            raise RuntimeError(f"{cmd} failed (rc={rc}): {out!r} {err!r}")

    # Fix single-volume nm5_gm.nal (PD-LOSAT-NCBI-DEFECTS Q6=B)
    (db_dir / "nm5_gm.nal").write_text("TITLE nm5_gm\nDBLIST nm5_gm\nMASKLIST gm_test\n")

    # Build a small v4 DB with empty 0-byte .phr record for OID 0 (`pz4`) to test filters_header_sets
    (db_dir / "pz4.fa").write_text(f">p0\n{PSEQ0}\n")
    rc, out, err = run_ncbi(ncbi_bin, lock, ["makeblastdb", "-in", "pz4.fa", "-dbtype", "prot", "-blastdb_version", "4", "-out", "pz4", "-title", "pz4"], db_dir)
    if rc != 0:
        raise RuntimeError(f"pz4 failed: {err!r}")
    # Zero out .phr and patch .pin header offsets so OID 0 has [0, 0) header bytes
    (db_dir / "pz4.phr").write_bytes(b"")
    pin = bytearray((db_dir / "pz4.pin").read_bytes())
    # v4 .pin ends with (num_oids + 1) hdr offsets + (num_oids + 1) seq offsets (each 4B big-endian)
    # 1 OID -> last 16 bytes are hdr[0], hdr[1], seq[0], seq[1]; zero hdr[0] and hdr[1]
    pin[-16:-8] = b"\x00" * 8
    (db_dir / "pz4.pin").write_bytes(bytes(pin))

    # Build OIDLIST file for pm4 (keeps all 6 OIDs: 0..5)
    # Format: 4B BE max_oid (5) + 1 byte bits (0b11111100 = 0xFC)
    (db_dir / "all6.msk").write_bytes((5).to_bytes(4, "big") + bytes([0xFC]))

    # Remove input staging files from db_dir
    for tmp in ["prot.asnt", "prot_b.asnt", "nucl.asnt", "gmask.asnt", "pz4.fa"]:
        (db_dir / tmp).unlink()

    # Create nogless subdirectory (nm5_gm without .nog, with gm_test.*)
    nog = db_dir / "nogless"
    nog.mkdir()
    for p in db_dir.iterdir():
        if p.name.startswith("nm5_gm.") and p.suffix not in (".nog", ".nal"):
            shutil.copy2(p, nog / p.name)
        elif p.name.startswith("gm_test."):
            shutil.copy2(p, nog / p.name)

    # Create helper list files in db_dir
    (db_dir / "gi1001.gil").write_text("1001\n1003\n1004\n")
    (db_dir / "si_p00001.txt").write_text("P00001.1\n")

    # Alias files in db_dir
    aliases = {
        "mb_0.pal": "DBLIST pm5\nMEMB_BIT 0\n",
        "mb_1.pal": "DBLIST pm5\nMEMB_BIT 1\n",
        "mb_2.pal": "DBLIST pm5\nMEMB_BIT 2\n",
        "mb_3.pal": "DBLIST pm5\nMEMB_BIT 3\n",
        "mb_4.pal": "DBLIST pm5\nMEMB_BIT 4\n",
        "mb_31.pal": "DBLIST pm5\nMEMB_BIT 31\n",
        "mb_32.pal": "DBLIST pm5\nMEMB_BIT 32\n",
        "mb_plus1.pal": "DBLIST pm5\nMEMB_BIT +1\n",
        "mb_abc.pal": "DBLIST pm5\nMEMB_BIT abc\n",
        "mb_neg1.pal": "DBLIST pm5\nMEMB_BIT -1\n",
        "mb_empty.pal": "DBLIST pm5\nMEMB_BIT \n",
        "mb_over.pal": "DBLIST pm5\nMEMB_BIT 4294967296\n",
        "ord_si_mb.pal": "DBLIST pm5\nSEQIDLIST a b\nMEMB_BIT abc\n",
        "ord_mb_tx.pal": "DBLIST pm5\nMEMB_BIT abc\nTAXIDLIST a b\n",
        "nest_1_1.pal": "DBLIST mb_1\nMEMB_BIT 1\n",
        "nest_1_0.pal": "DBLIST mb_0\nMEMB_BIT 1\n",
        "nest_0_1.pal": "DBLIST mb_1\nMEMB_BIT 0\n",
        "conflict_plain.pal": "DBLIST mb_1 mb_2\n",
        "conflict_nseq.pal": "DBLIST mb_1 mb_2\nNSEQ 6\nLENGTH 360\n",
        "conflict_gi.pal": "DBLIST mb_1 mb_2\nGILIST gi1001.gil\n",
        "direct_and_sub.pal": "DBLIST pm5 mb_1\n",
        "mb_b_2.pal": "DBLIST pm5_b\nMEMB_BIT 2\n",
        "toprange_mb0.pal": "DBLIST pm5\nFIRST_OID 2\nLAST_OID 4\nMEMB_BIT 0\n",
        "mb1_range.pal": "DBLIST pm5\nFIRST_OID 1\nLAST_OID 3\nMEMB_BIT 1\n",
        "mb1_om1.pal": "DBLIST pm4\nOIDLIST all6.msk\nOID_MASK_TYPE 1\nMEMB_BIT 1\n",
        "mb1_gi.pal": "DBLIST pm4\nGILIST gi1001.gil\nMEMB_BIT 1\n",
        "mb2_si_v4.pal": "DBLIST pm4\nSEQIDLIST si_p00001.txt\nMEMB_BIT 2\n",
        "pz_mb0.pal": "DBLIST pz4\nMEMB_BIT 0\n",
        "pz_mb1.pal": "DBLIST pz4\nMEMB_BIT 1\n",
        "nmb_gm_1.nal": "TITLE nmb_gm_1\nDBLIST nm5_gm\nMASKLIST gm_test\nMEMB_BIT 1\n",
        "nmb_gm_2.nal": "TITLE nmb_gm_2\nDBLIST nm5_gm\nMASKLIST gm_test\nMEMB_BIT 2\n",
        "nmb_mask_empty.nal": "TITLE empty\nDBLIST nm5_gm\nMASKLIST \n",
        "nc1.nal": "DBLIST nm5_gm\nMEMB_BIT 1\n",
        "nc2.nal": "DBLIST nm5_gm\nMEMB_BIT 2\n",
        "nc_nseq.nal": "DBLIST nc1 nc2\nNSEQ 3\nLENGTH 240\n",
    }
    for name, content in aliases.items():
        (db_dir / name).write_text(content)

    # Alias files in nogless
    (nog / "a1.nal").write_text("DBLIST nm5_gm\nMASKLIST gm_test\nMEMB_BIT 1\n")
    (nog / "a2.nal").write_text("DBLIST nm5_gm\nMASKLIST gm_test\nMEMB_BIT 2\n")
    (nog / "c1.nal").write_text("DBLIST nm5_gm\nMEMB_BIT 1\n")
    (nog / "c2.nal").write_text("DBLIST nm5_gm\nMEMB_BIT 2\n")
    (nog / "cc.nal").write_text("DBLIST c1 c2\nMASKLIST gm_test\nNSEQ 3\nLENGTH 240\n")

    # Queries in work/in/
    in_dir = work / "in"
    in_dir.mkdir(parents=True)
    (in_dir / "qp.fa").write_text(
        f">q0\n{PSEQ0}\n>q1\n{PSEQ1}\n>q2\n{PSEQ2}\n>q4\n{PSEQ4}\n"
    )
    (in_dir / "qn.fa").write_text(
        f">qn0\n{NSEQ0}\n>qn1\n{NSEQ1}\n>qn2\n{NSEQ2}\n"
    )


CASES = [
    # (name, cwd_rel, words)
    # Basic MEMB_BIT filtering in blastdbcmd
    ("bdc_mb0", "db", ["blastdbcmd", "-db", "mb_0", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_mb1", "db", ["blastdbcmd", "-db", "mb_1", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_mb2", "db", ["blastdbcmd", "-db", "mb_2", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_mb3", "db", ["blastdbcmd", "-db", "mb_3", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_mb4", "db", ["blastdbcmd", "-db", "mb_4", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_mb31", "db", ["blastdbcmd", "-db", "mb_31", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_mb32", "db", ["blastdbcmd", "-db", "mb_32", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_mb_plus1", "db", ["blastdbcmd", "-db", "mb_plus1", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_mb1_fasta", "db", ["blastdbcmd", "-db", "mb_1", "-dbtype", "prot", "-entry", "all"]),
    ("bdc_mb1_ctrla", "db", ["blastdbcmd", "-db", "mb_1", "-dbtype", "prot", "-entry", "all", "-ctrl_a"]),
    ("bdc_mb1_target", "db", ["blastdbcmd", "-db", "mb_1", "-dbtype", "prot", "-entry", "P00001.1", "-target_only", "-outfmt", "%o|%a|%t"]),
    ("bdc_mb1_info", "db", ["blastdbcmd", "-db", "mb_1", "-dbtype", "prot", "-info"]),
    # Conversion errors and exact ComputeMasks parse order
    ("bdc_mb_abc", "db", ["blastdbcmd", "-db", "mb_abc", "-dbtype", "prot", "-info"]),
    ("bdc_mb_neg1", "db", ["blastdbcmd", "-db", "mb_neg1", "-dbtype", "prot", "-info"]),
    ("bdc_mb_empty", "db", ["blastdbcmd", "-db", "mb_empty", "-dbtype", "prot", "-info"]),
    ("bdc_mb_over", "db", ["blastdbcmd", "-db", "mb_over", "-dbtype", "prot", "-info"]),
    ("bdc_ord_si_mb", "db", ["blastdbcmd", "-db", "ord_si_mb", "-dbtype", "prot", "-info"]),
    ("bdc_ord_mb_tx", "db", ["blastdbcmd", "-db", "ord_mb_tx", "-dbtype", "prot", "-info"]),
    # Nesting, conflicts, direct volume, multi-volume, ranges, OID_MASK_TYPE, GILIST, SEQIDLIST
    ("bdc_nest_1_1", "db", ["blastdbcmd", "-db", "nest_1_1", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_nest_1_0", "db", ["blastdbcmd", "-db", "nest_1_0", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_nest_0_1", "db", ["blastdbcmd", "-db", "nest_0_1", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_conflict_plain_info", "db", ["blastdbcmd", "-db", "conflict_plain", "-dbtype", "prot", "-info"]),
    ("bdc_conflict_plain_all", "db", ["blastdbcmd", "-db", "conflict_plain", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a"]),
    ("bdc_conflict_gi_info", "db", ["blastdbcmd", "-db", "conflict_gi", "-dbtype", "prot", "-info"]),
    ("bdc_direct_and_sub", "db", ["blastdbcmd", "-db", "direct_and_sub", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_multivol", "db", ["blastdbcmd", "-db", "mb_1 mb_b_2", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_toprange_mb0", "db", ["blastdbcmd", "-db", "toprange_mb0", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a"]),
    ("bdc_mb1_range", "db", ["blastdbcmd", "-db", "mb1_range", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_mb1_om1", "db", ["blastdbcmd", "-db", "mb1_om1", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_mb1_gi", "db", ["blastdbcmd", "-db", "mb1_gi", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    ("bdc_mb2_si_v4", "db", ["blastdbcmd", "-db", "mb2_si_v4", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a|%t"]),
    # Empty .phr header record (filters_header_sets)
    ("bdc_pz_mb0", "db", ["blastdbcmd", "-db", "pz_mb0", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a"]),
    ("bdc_pz_mb1", "db", ["blastdbcmd", "-db", "pz_mb1", "-dbtype", "prot", "-entry", "all", "-outfmt", "%o|%a"]),
    # MASKLIST + MEMB_BIT + nogless + conflict
    ("bdc_nmb_gm1", "db", ["blastdbcmd", "-db", "nmb_gm_1", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o|%a|%m0"]),
    ("bdc_nmb_gm2", "db", ["blastdbcmd", "-db", "nmb_gm_2", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o|%a|%m0"]),
    ("bdc_mask_empty_info", "db", ["blastdbcmd", "-db", "nmb_mask_empty", "-dbtype", "nucl", "-info"]),
    ("bdc_nogless_a1", "db/nogless", ["blastdbcmd", "-db", "a1", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o|%a|%m0"]),
    ("bdc_nogless_a2", "db/nogless", ["blastdbcmd", "-db", "a2", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o|%a|%m0"]),
    ("bdc_nogless_cc", "db/nogless", ["blastdbcmd", "-db", "cc", "-dbtype", "nucl", "-entry", "all", "-outfmt", "%o|%a|%m0"]),
    # blastdbcheck (verbosity 3, -no_isam, deterministic -stride 1)
    ("chk_mb1", "db", ["blastdbcheck", "-db", "mb_1", "-dbtype", "prot", "-verbosity", "3", "-no_isam", "-stride", "1"]),
    ("chk_conflict_plain", "db", ["blastdbcheck", "-db", "conflict_plain", "-dbtype", "prot", "-verbosity", "3", "-no_isam", "-stride", "1"]),
    ("chk_conflict_gi", "db", ["blastdbcheck", "-db", "conflict_gi", "-dbtype", "prot", "-verbosity", "3", "-no_isam", "-stride", "1"]),
    # Search programs (blastp, blastx, psiblast, blastn, tblastn, tblastx)
    ("bp_mb1_f6", "db", ["blastp", "-query", "../in/qp.fa", "-db", "mb_1", "-outfmt", "6"]),
    ("bp_mb2_f6", "db", ["blastp", "-query", "../in/qp.fa", "-db", "mb_2", "-outfmt", "6"]),
    ("bp_mb4_lineless_f6", "db", ["blastp", "-query", "../in/qp.fa", "-db", "mb_4", "-outfmt", "6"]),
    ("bp_mb32_f6", "db", ["blastp", "-query", "../in/qp.fa", "-db", "mb_32", "-outfmt", "6"]),
    ("bp_mb1_f0", "db", ["blastp", "-query", "../in/qp.fa", "-db", "mb_1", "-outfmt", "0"]),
    ("bp_mb2_user_si", "db", ["blastp", "-query", "../in/qp.fa", "-db", "mb_2", "-seqidlist", "si_p00001.txt", "-outfmt", "6"]),
    ("bp_mb1_om1", "db", ["blastp", "-query", "../in/qp.fa", "-db", "mb1_om1", "-outfmt", "6"]),
    ("bp_conflict_nseq", "db", ["blastp", "-query", "../in/qp.fa", "-db", "conflict_nseq", "-outfmt", "6"]),
    ("bp_conflict_gi", "db", ["blastp", "-query", "../in/qp.fa", "-db", "conflict_gi", "-outfmt", "6"]),
    ("psi_mb1_f6", "db", ["psiblast", "-query", "../in/qp.fa", "-db", "mb_1", "-outfmt", "6"]),
    ("psi_conflict_nseq", "db", ["psiblast", "-query", "../in/qp.fa", "-db", "conflict_nseq", "-outfmt", "6"]),
    ("bx_mb1_f6", "db", ["blastx", "-query", "../in/qn.fa", "-db", "mb_1", "-outfmt", "6"]),
    ("bx_mb2_f6", "db", ["blastx", "-query", "../in/qn.fa", "-db", "mb_2", "-outfmt", "6"]),
    ("bn_nmb_gm1_hm", "db", ["blastn", "-query", "../in/qn.fa", "-db", "nmb_gm_1", "-dust", "no", "-db_hard_mask", "0", "-outfmt", "6"]),
    ("bn_nmb_gm2_hm", "db", ["blastn", "-query", "../in/qn.fa", "-db", "nmb_gm_2", "-dust", "no", "-db_hard_mask", "0", "-outfmt", "6"]),
    ("bn_nogless_a1_hm", "db/nogless", ["blastn", "-query", "../../in/qn.fa", "-db", "a1", "-dust", "no", "-db_hard_mask", "0", "-outfmt", "6"]),
    ("bn_nogless_a2_hm", "db/nogless", ["blastn", "-query", "../../in/qn.fa", "-db", "a2", "-dust", "no", "-db_hard_mask", "0", "-outfmt", "6"]),
    ("bn_nc_nseq", "db", ["blastn", "-query", "../in/qn.fa", "-db", "nc_nseq", "-outfmt", "6"]),
    ("tbn_nmb_gm1", "db", ["tblastn", "-query", "../in/qp.fa", "-db", "nmb_gm_1", "-outfmt", "6"]),
    ("tbx_nmb_gm2", "db", ["tblastx", "-query", "../in/qn.fa", "-db", "nmb_gm_2", "-outfmt", "6"]),
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ncbi-bin", required=True, type=Path)
    ap.add_argument("--work", required=True, type=Path)
    ap.add_argument("--out", default=Path(__file__).resolve().parent / "memb_bit", type=Path)
    ap.add_argument("--lock", default=None, type=Path)
    args = ap.parse_args()

    work = args.work.resolve()
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    build_databases(args.ncbi_bin.resolve(), args.lock, work)

    out = args.out.resolve()
    if out.exists():
        shutil.rmtree(out)
    shutil.copytree(work / "db", out / "db")
    shutil.copytree(work / "in", out / "in")

    exp_dir = out / "expected"
    exp_dir.mkdir(parents=True)

    root_bytes = str(work).encode("utf-8")
    case_lines = []
    for name, cwd_rel, words in CASES:
        rc, stdout, stderr = run_ncbi(args.ncbi_bin.resolve(), args.lock, words, work / cwd_rel)
        (exp_dir / f"{name}.rc").write_text(f"{rc}\n")
        (exp_dir / f"{name}.stdout").write_bytes(stdout.replace(root_bytes, b"<R>"))
        (exp_dir / f"{name}.stderr").write_bytes(stderr.replace(root_bytes, b"<R>"))
        case_lines.append("\t".join([name, cwd_rel] + words))

    (out / "cases.tsv").write_text("\n".join(case_lines) + "\n")
    print(f"Froze {len(CASES)} cases into {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
