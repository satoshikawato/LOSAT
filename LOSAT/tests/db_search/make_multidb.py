#!/usr/bin/env python3
"""Builds tests/db_search/multidb/: the databases of the rows db.mdb.* (DI-8f: -db with several names).

Usage: make_multidb.py [--ordsl-only | --mixdb2-only | --mixtax-only | --mixsfx-only | --loadord-only | --loadord2-only | --loadord3-only | --loadord4-only | --loadord5-only | --loadord6-only | --loadord7-only | --v4late-only | --pfoff-only | --pfoff2-only | --pfoff3-only | --rf1-only | --rf2-only] MAKEBLASTDB   (writes multidb/*.fa, multidb/*.map and the databases below)

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
- `pgim`: version 5 protein, `-parse_seqids`, records 1-2 of tfp.faa whose only Seq-id is a giim (`gim|70001`, `gim|70002`;
  line DBID mixdb, review 1 H1: NCBI's prefetch looks a giim up in every LMDB file of a list).
- `al.nal`: an alias file (TITLE multidb alias of na, DBLIST na).
- `mix_u.nal`: an alias file of `na nnp` whose FIRST_OID 9 and LAST_OID 17 keep the sequences of nnp alone (line DBID mixdb).
- `sl.txt`: a seqid list of the first accession of na, nb, pa and pb.
- `ov_p` / `ov_n` (line DBID ordsl): `makeblastdb -input_type blastdb` of ../taxfilter/tfu_p and tfu_n (version 5, no
  parsed Seq-ids, four records) in two volumes of two records (`-max_file_sz` 400 and 500 bytes) under one LMDB file:
  makeblastdb stores `gnl|BL_ORD_ID|n` with n the OID in the database and keys the accession table with it, and NCBI's
  header reads add the volume's first OID (2) to the stored n of the second volume.
- `sl_ord0.txt`, `sl_ord1.txt`, `sl_ord2.txt`, `sl_ord9.txt`, `sl_ord13.txt`: the seqid list `gnl|BL_ORD_ID|n`;
  `sl_ord1.bsl` the binary list of `sl_ord1.txt` (`blastdb_aliastool -seqid_file_out`).
- `pxm` / `nxm` (line DBID ordsl, review 1): version 5, `-parse_seqids`, two records named to sort after `pa` / `na`: the
  letters of record 13 of tfp.faa (record 9 of tfn.fna) under two definition lines, `ref|NP_900003.1|` (`NM_`) and
  `gnl|BL_ORD_ID|7`, and a random sequence (seed 20261115) under `ref|NP_900005.1|gnl|BL_ORD_ID|0`, whose tag makeblastdb
  stores as the string "0" (a tag with a leading 0 is not made a number). `xq.faa`: that random protein as a query.
  Lists `sl_pxm7.txt` and `sl_pxm19.txt` (`NP_900003.1` and `gnl|BL_ORD_ID|7` / `19`: 7 plus the 12 sequences of pa),
  `sl_nxm7.txt`, `sl_ord7.txt`, `sl_ord15.txt` (7 plus the 8 of na), `sl_ord19.txt`, `sl_np5.txt` (`NP_900005.1`).
  `make_multidb.py --ordsl-only MAKEBLASTDB` builds only these and the ones above.
- `put` (line DBID mixdb2, review 1 H1): version 5 protein without `-parse_seqids`, the three proteins of
  ../ids/ids_q.faa (those of ../ids/noacc_v5) with taxid 9606 (`-taxid`), so that its loads look the taxonomy up in a
  list with `pa`.
- `mix2_q.faa`, `mix2_qx.fna`, `mix2_q.fna`, `mix2_qt.faa` (line DBID mixdb2): two queries each, the first of which
  hits only a database without the accession table (`-evalue 1e-10`), the second the parsed one (../taxfilter/q.faa or
  q.fna): `qp1` of ../ids/ids_q.faa (noacc_v5, put), the same back-translated (one codon per amino acid), the
  `FBtr0083245` record of nnp.fa (glob1), and the longest ORF of that record (153 aa).
  `make_multidb.py --mixdb2-only MAKEBLASTDB` builds only `put` and these.
- `nut` (line DBID mixtax): version 5 nucleotide without `-parse_seqids`, the records of nnp.fa with taxid 9606
  (`-taxid`): `put`'s counterpart for the nucleotide programs, whose loads look the taxonomy up in a list with `na`.
  `mix2_qc.faa` and `mix2_qcx.fna` (line DBID mixtax): the first query of mix2_q.faa and mix2_qx.fna, then one whose
  letters are those of ../taxfilter/q.faa followed by `qp1` (back-translated in the second), which hits both `pa` and
  `put`. `putl` and `nutl` (line DBID mixtax, review 1 L2): `qp2` of ../ids/ids_q.faa (back-translated in `nutl`) with
  taxid 99990, the last entry of ../../seqdb/taxdb/<variant>, without the table; `mix2_ql.faa`, `mix2_qlx.fna`,
  `mix2_ql.fna` and `mix2_qlt.faa`: the first query of mix2_q.faa, mix2_qx.fna, mix2_q.fna and mix2_qt.faa, then the
  letters of ../taxfilter/q.faa or q.fna followed by `qp2` (back-translated in the nucleotide files), which hits `pa` or
  `na`, `put` or `nut` and the volume `putl` or `nutl` that the first query does not.
  `make_multidb.py --mixtax-only MAKEBLASTDB` builds only `nut`, `putl`, `nutl` and these.
- `putb` (line RG item T5c, merge 3): `put` with protein titles ending in `]` (`alpha, partial [Homo sapiens]`,
  `chi, partial [ ]`, `MAG: (plastid) [ ]` for the three proteins of ../ids/ids_q.faa, taxid 9606): in a list with
  `pa`, the first batch of mix2_qc.faa / mix2_qcx.fna reports `alpha, ...`, whose end NCBI's
  `x_AdjustProteinTitleSuffix` rewrites in the tabular ids but not in the XML2 and JSON `Hit` titles (the definition
  line's), before the second batch's prefetch stops. `make_multidb.py --mixsfx-only MAKEBLASTDB` builds only `putb`.
- `ldp` and `ldpn` (line DBID loadord): version 5, `-parse_seqids`, taxid 9606: `ref|NP_000002.1|` (`NM_`) with the letters
  of ../taxfilter/q.faa (q.fna), then `ref|NP_000900.1|` (`NM_`) with random letters (seed 20261116). Beside a database
  without the accession table, NCBI's prefetch loads the batch's Seq-align ids in the order of their Seq-id handles,
  `gnl|BL_ORD_ID|n` by n (0 last of the numbers) and the accession by its number 2, and stops at the accession's lookup:
  `lda` (`ldan`) holds `qp1` of ../ids/ids_q.faa (the glob1 record `FBtr0083245` of nnp.fa) alone, at OID 0 in a list with
  `ldp` (`ldpn`), which loads after the stop; `ldb` (`ldbn`) a random sequence and then that one, at OID 1, which loads
  before it; `ldu` (`ldun`) that one alone, sorted after `ldp` (`ldpn`), at OID 2: the same packed number 2, which NCBI
  orders by the address of the handles. Without `-parse_seqids`, taxid 9606 (`-taxid`), so that their loads open the
  taxonomy database. `ldg` (`ldgn`): `gi|1|ref|NP_000002.1|` (`NM_`), the letters of `ldp`'s first record under a gi whose
  load (handle 1, before the accession's 2) registers the accession in NCBI's data loader, so that NCBI never looks it up
  (review 1 of loadord, H1); `ldq` (`ldqn`): `lcl|Query_1` with those letters, which NCBI's scope finds among the queries
  (H2); both `-parse_seqids`. `ldo_qc.faa`: one query, the letters of ../taxfilter/q.faa followed by `qp1`, which hits both
  (psiblast's batch is one query). `make_multidb.py --loadord-only MAKEBLASTDB` builds only these.
- `ldc` (`ldcn`) and `ldh` (`ldhn`) (line DBID loadord, DECISIONS 60 (6)): version 5, `-parse_seqids`, taxid 9606, the
  letters of `ldp`'s first record: `ldc` a copy of it (`ref|NP_000002.1|`), the same Seq-align id, whose load is the
  stopping lookup itself (beside `ldp`, NCBI stops in a later batch too); `ldh` `gi|9000|ref|NP_000002.1|`, which holds
  the accession but whose gi (handle 9000) loads after it in one batch, so NCBI stops. `make_multidb.py
  --loadord2-only MAKEBLASTDB` builds only these (the databases of `--loadord-only` keep their dates).
- Review 2 of line DBID loadord (DECISIONS 60 (7), (8)): version 5, `-parse_seqids`, taxid 9606. `ldq2` (`ldq2n`):
  `lcl|Query_2` with the letters of `ldp`'s first record, the generated id of a second query; `ldx` (`ldxn`):
  `gi|9001|ref|NP_000002.1|` (`NM_`) with random letters (seed 20261116 + 3), a holder of the stopping accession that no
  query hits; `ldd` (`lddn`): `gi|9000|ref|NP_000002.1|` (`NM_`) with the letters of `lda` (`ldan`), a holder that the
  first query of mix2_q.faa (mix2_qx.fna, mix2_q.fna, mix2_qt.faa) hits in a batch before the stop; `ldmy`:
  `lcl|myprot` with the letters of `ldp`'s first record. `ldo_short.faa`: the first 40 letters of `qp1`, then the
  letters of ../taxfilter/q.faa, so that `-query_loc 50-140` drops the first record (which keeps `Query_1`);
  `ldo_short.fna` and `ldo_shortx.fna`: the first 120 letters of glob1 (of qp1 back-translated), then
  ../taxfilter/q.fna, for `-query_loc 150-420`. `ldo_my.asn`: ../../fixtures/pssm/q1_tfp.asn.gz with the query's
  local id `Query_1` changed to `myprot` (the PSSM of the globin of ../taxfilter/q.faa). `make_multidb.py
  --loadord3-only MAKEBLASTDB` builds only these (the databases before keep their dates).
- Review 3 of line DBID loadord (DECISIONS 60 (10)): version 5, `-parse_seqids`, taxid 9606, the letters of `ldp`'s
  first record, the stopping accessions beside it in one batch: `ldr` (`ldrn`) `ref|NP_000001.1|` (`NM_`, packed 1,
  before `ldp`'s 2), `lds` (`ldsn`) `ref|NP_000003.1|` (packed 3, after it), and `ldi` (`ldin`)
  `gi|1|ref|NP_000003.1|`, a holder of `lds`'s accession whose gi (handle 1) loads before both. `make_multidb.py
  --loadord4-only MAKEBLASTDB` builds only these (the databases before keep their dates).
- Review 4 of line DBID loadord (DECISIONS 60 (12)): version 4, `-parse_seqids`, taxid 9606, one record each,
  `ref|NP_000002.1|` (`NM_`) with the letters of `ldp`'s first record (copies of one record): `l5c` (`l5cn`), `l5d`
  (`l5dn`) and `l5e` (`l5en`); `l5d`'s string ISAM index (`.psi`, `.nsi`) is a directory holding a placeholder file,
  which NCBI's `CSeqDBIsam::IndexExists` (`CFile::Exists`, a regular file) takes for no ISAM: its lookups find nothing
  in that volume and go on to the next (H1). `make_multidb.py --loadord5-only MAKEBLASTDB` builds only these (the
  databases before keep their dates).
- Review 5 of line DBID loadord (DECISIONS 60 (15), line RG K4h's E2 in loadord's shapes): `ldl` (`ldln`), version 5
  without `-parse_seqids`, taxid 99990 (the last entry of ../../seqdb/taxdb/offset_last_past_end, whose lookup writes
  the offset error), two records of `lda`'s letters (`qp1`, glob1): beside `ldp`, OID 0 (`gnl|BL_ORD_ID|0`) loads after
  the stop and OID 1 before it, in one volume. The queries `ldo_e2.faa` (BLASTP), `ldo_e2x.fna` (BLASTX), `ldo_e2.fna`
  (BLASTN, TBLASTX) and `ldo_e2t.faa` (TBLASTN): `q0` without hits, `q1` titled `q1 a<U+00E9>b z` (in UTF-8, at which NCBI's
  XML stream cuts its outfmt 5 `Iteration`, so its report makes none of its loads) with `lda`'s letters repeated (qp1
  16 times, back-translated in `ldo_e2x.fna`; glob1 6 times in `ldo_e2.fna`; glob1's ORF 16 times in `ldo_e2t.faa`), so
  that its HSPs turn its batch's prefetch off under `PRE_FETCH_SEQS_LIMIT`, and `q2` with the letters of
  ../taxfilter/q.faa (q.fna) and then those of `lda` (`ldan`) once, which hits the stopping accession and the unparsed
  record. TBLASTX counts a Seq-align per subject (an ungapped search), so its `q1` must hit more subjects than `q2`:
  `ldyn`, version 5 without `-parse_seqids`, taxid 9606, four records of random letters (seed 20261116 + 5), and
  `ldo_e2y.fna`, whose `q1` holds glob1 and then those letters. `make_multidb.py --loadord6-only MAKEBLASTDB` builds
  only these (the databases before keep their dates).
- Review 6 of line DBID loadord (DECISIONS 60 (17), H1): `ldv`, version 5, `-parse_seqids`, taxid 9606,
  `gi|9002|ref|NP_000002.1|` (a holder of `ldp`'s stopping accession) with a variant of `lda`'s letters (qp1, one letter
  in ten replaced, seed 20261116 + 7): the first query of mix2_q.faa hits it after `lda`'s exact copy, so psiblast's
  traceback holds it beyond `-max_target_seqs 1`, which no report shows. `make_multidb.py --loadord7-only MAKEBLASTDB`
  builds only it (the databases before keep their dates).
- Line DBID v4late (review 3 of line DBID loadord, H3): version 4 protein, `-parse_seqids`, taxid 9606, one record
  each, beside `l5c` in a list whose later volume's string ISAM a test makes unreadable (mode 000, which git does not
  keep): `l6o` holds `l5c`'s `ref|NP_000002.1|` with other letters (qp1 of ../ids/ids_q.faa), so NCBI's lookup of its
  Seq-align id gives `l5c`'s record first and never opens `l6o`'s ISAM; `l6u` holds `ref|NP_000004.1|`, which only its
  ISAM can hold, so NCBI's lookup opens it; `l5b` (a name before `l5c`'s) holds `ref|NP_000999.1|ref|NP_000002.1|`
  with qp1's letters, so that NCBI's lookup of `l5e`'s id gives it first, whose definition line `GetBioseq` rejects
  for that id (`s_SeqDB_SeqIdIn` stops at the first id, `oid headers do not contain target gi/seq_id`; review 1 of
  the item, H1). `make_multidb.py --v4late-only MAKEBLASTDB` builds only these (the databases before keep their
  dates).
- `pfo_q.faa`, `pfo_qx.fna`, `pfo_q.fna`, `pfo_qt.faa` (line DBID pfoff): the two queries of mix2_q.faa, mix2_qx.fna,
  mix2_q.fna and mix2_qt.faa in the other order, so that with `PRE_FETCH_SEQS_LIMIT=1 BATCH_SIZE=1` (the prefetch
  off) the first query's report reads a parsed subject whose lookup throws (`DB contains no accession info.`) and the
  second's does not. `pfo_qf.fna`, `pfo_qft.faa` (rows `.taxfirst.` of the nucleotide databases, review 1 of item
  pfoff-1, L1): `na`'s HBB record (`ENST00000335295.4`; `pfo_qft.faa` its longest ORF, translated), then the query of
  mix2_q.fna (mix2_qt.faa) that hits `nnp`/`nut` alone, so that the first query's first `Hit` is a parsed subject of
  `na` whose lookup throws and no subject of `nut` (with a taxid) comes before it, where pfo_q.fna's first query
  shows `nut`'s copy of HBA1 (`ENSTTEST0001.1`) first. `pfo_qs.faa`, `pfo_qsx.fna`, `pfo_qs.fna`, `pfo_qst.faa`,
  `pfo_qss.fna` (rows `.taxstop.`, the same review): a first query whose report throws at its first `Hit` and shows a
  subject of `put`/`nut` (with a taxid) after it, which NCBI never loads, with more Seq-aligns than the second, whose
  batch's prefetch (on below the rows' `PRE_FETCH_SEQS_LIMIT`) loads such a subject and then stops at one whose lookup
  throws: BLASTP's and BLASTX's first query the letters of mix2_q's (mix2_qx's) parsed query, then half of its other
  query's, twice; the nucleotide programs' `na`'s HBB twice (its ORF for TBLASTN), then half of the query that hits
  `nnp`; the second query that query's letters, then the parsed one's (TBLASTX: a third of them, then the first 200
  of the parsed one's). `make_multidb.py --pfoff-only MAKEBLASTDB` writes only these (no database; MAKEBLASTDB is not
  run).
- `pfz`, `pfzq`, `pfo_qz.faa`, `pfo_qh.fna` (line DBID pfoff-2, review 2 of the item, M1, and review 1, H1): `pfz` (version
  5, `-parse_seqids`) holds `ref|XP_900001.1|`, whose alignment to `pfo_qz.faa` (`-threshold 9 -comp_based_stats 0`)
  has no identities, then `ref|XP_900002.1|`, whose alignment has some; `pfzq` the first alone. In a list with
  `ids/noacc_v5` both throw where the batch's prefetch is off, so a row of `nident` alone reads nothing of the first
  (it prints 0 and the run goes on) and ends the run at the second. `pfo_qh.fna`: `na`'s HBB record, whose rows all
  throw in a list with `nut`, then twice the letters of the query of mix2_q.fna that hits `nnp`/`nut` alone, whose two
  HSPs turn its batch's prefetch off at `PRE_FETCH_SEQS_LIMIT=1`, so that its rows load the subject and look its name
  up. `make_multidb.py --pfoff2-only MAKEBLASTDB` builds only these.
- `pf3p`, `pf3u`, `pfo_qr.faa` (line DBID pfoff-3, review 1 of the item, H1): `pf3p` (version 5, `-parse_seqids`)
  holds `ref|XP_810001.1|`, two separated copies of segments of `pfo_qr.faa`, so its alignment has two HSPs; `pf3u`
  (version 5, no accession table, taxid 9606) holds `bsub`, a mutated copy of one segment, weak enough that psiblast's
  first round does not report it and its second round does. In the list `pf3p pf3u` with the prefetch off, psiblast's
  second round shows the two-part description table of the first `-num_descriptions` HSPs: with 2, the two HSPs of
  `XP_810001.1` (whose lookup throws) fill it, and `bsub` is neither shown nor loaded (no text of a damaged taxonomy
  database); with 3 it is. `pfo_qm.faa` (review 2 of the item, L7) holds `mix2_qc.faa`'s `mix2_c` and then
  `pfo_q.faa`'s `P69905.2`, whose reports on `pa noacc_v5` have 13 and 12 Seq-aligns: with `BATCH_SIZE=1
  PRE_FETCH_SEQS_LIMIT=12` the first batch's prefetch is off and the second's on, and that prefetch throws.
  `pfo_qe.faa` (review 3 of the item, N1) holds `pfo_qr.faa`'s query, whose report on `pf3p pf3u` ends the run at
  `XP_810001.1` where its batch's prefetch is off (`BATCH_SIZE=1 PRE_FETCH_SEQS_LIMIT=1`), then `bsub`'s letters, whose
  batch has one Seq-align, so its prefetch is on and loads `bsub` (with `BLASTDB=btd_empty` NCBI crashes there): NCBI
  never gets there. `make_multidb.py --pfoff3-only MAKEBLASTDB` builds only these.
- `n3p`, `n3u` and the queries `rf1_*.faa`, `rf1_*.fna` (line DBID rf-1, 2026-11-27; the rows db.rf1.*): `n3p` and `n3u`
  are `pf3p` and `pf3u` back-translated (one codon for each residue, `RF1_CODONS`), nucleotide, with the same options;
  `rf1_q0.faa` / `rf1_q0.fna` hold `q0_nohit`, which hits nothing; `rf1_q0e` adds `bsub`'s letters (back-translated in
  the `.fna`), whose batch's prefetch loads `bsub` (taxid 9606), and `rf1_e.faa` holds those letters alone (psiblast
  searches each query alone: its first query); `rf1_r` holds `pfo_qe.faa`'s `pfo_r` (back-translated
  in the `.fna`), whose report shows `pf3p`'s hit, and `rf1_re` adds `bsub`'s letters. With `BLASTDB=btd_empty` (an
  empty taxdb.btd) NCBI crashes at the load of `bsub`, after the reports of the batches before (with `BATCH_SIZE=1`,
  those of `q0_nohit` or `pfo_r`). `rf1_q0t.faa` holds `q0_nohit`, then `taxfilter/q.faa`'s P69905, whose report shows
  taxfilter's `tfe_nr` OID 3 (an empty stored title, taxid 9606: NCBI makes the title from the name of the taxid, and
  crashes at the sequence's load). `make_multidb.py --rf1-only MAKEBLASTDB` builds only these.
- `rf2p`, `rf2n` and the query `rf2_q0t.fna` (line DBID rf-2, 2026-11-28; the rows db.rf2.*): one sequence each with parsed
  Seq-ids and an empty title (`ref|XP_820003.1|`, `ref|XM_820003.1|`; the letters of taxfilter's `tfe.faa` and `tfe.fna`
  record 3), `-parse_seqids -taxid 9606`, so that NCBI makes the title from the name of the taxid (with
  `BLASTDB=btd_empty` it crashes at the sequence's load); `rf2_q0t.fna` holds `rf1_q0.fna`'s `q0_nohit`, then
  `taxfilter/q.fna`'s HBA1 cDNA, which hits both (as `rf1_q0t.faa`'s P69905 does). `make_multidb.py --rf2-only
  MAKEBLASTDB` builds only these, and the queries `rf1_rps.faa` and `rf1_rps.fna` (line DBID rf-1 review 1 L3; rpsblast's
  and rpstblastn's rows db.rf1.rps.*, db.rf1.rpst.*): `q0_nohit`, then `../../fasta/rps/q649.faa`'s YP_009220649.1
  (back-translated in the `.fna`), which hits a profile of each volume of the taxid copy of `rps_mini3v`.

`.pjs` and `.njs` (the metadata with the build's paths and dates; no search reads them) are dropped.
"""
import os
import random
import re
import subprocess
import sys

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "multidb")
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "taxfilter")


def makeblastdb(binary, *args):
    subprocess.run([binary, *args], check=True, cwd=HERE, stdout=subprocess.DEVNULL)


def read_fasta(name):
    return read_fasta_path(os.path.join(SRC, name))


def read_fasta_path(path):
    records = []
    with open(path) as f:
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


def make_pgim(binary, tfp):
    """Writes pgim.fa (records 1-2 of tfp.faa, each under a giim id before its title) and builds pgim."""
    records = [[f">gim|{70001 + i} " + header[1:].split(" ", 1)[1].split("\x01")[0].rstrip("\n") + "\n", lines]
               for i, (header, lines) in enumerate(tfp[:2])]
    write_set("pgim", records, None)
    makeblastdb(binary, "-in", "pgim.fa", "-dbtype", "prot", "-out", "pgim", "-parse_seqids", "-blastdb_version", "5",
                "-title", "multidb pgim")


def make_ordsl(binary):
    """Builds ov_p and ov_n and writes the ordinal-id seqid lists (line DBID ordsl)."""
    for kind, source, base, size in (("prot", "tfu_p", "ov_p", "400"), ("nucl", "tfu_n", "ov_n", "500")):
        makeblastdb(binary, "-in", os.path.join(SRC, source), "-input_type", "blastdb", "-dbtype", kind, "-out", base,
                    "-blastdb_version", "5", "-max_file_sz", size, "-title", "multidb " + base)
        os.remove(os.path.join(HERE, f"{base}.{kind[0]}js"))
    for n in (0, 1, 2, 9, 13):
        with open(os.path.join(HERE, f"sl_ord{n}.txt"), "w") as f:
            f.write(f"gnl|BL_ORD_ID|{n}\n")
    aliastool = os.path.join(os.path.dirname(binary), "blastdb_aliastool")
    subprocess.run([aliastool, "-seqid_file_in", "sl_ord1.txt", "-seqid_file_out", "sl_ord1.bsl"], check=True, cwd=HERE,
                   stdout=subprocess.DEVNULL)
    rng = random.Random(20261115)
    tfp, tfn = read_fasta("tfp.faa"), read_fasta("tfn.fna")
    for base, kind, source, prefix, alphabet, length in (("pxm", "prot", tfp[12], "NP", "ACDEFGHIKLMNPQRSTVWY", 120),
                                                          ("nxm", "nucl", tfn[8], "NM", "ACGT", 300)):
        letters = "".join(rng.choice(alphabet) for _ in range(length))
        random_lines = [letters[at:at + 60] + "\n" for at in range(0, length, 60)]
        records = [[f">ref|{prefix}_900003.1| two-line t3a\x01gnl|BL_ORD_ID|7 two-line t3b\n", source[1]],
                   [f">ref|{prefix}_900005.1|gnl|BL_ORD_ID|0 string tag t5\n", random_lines]]
        write_set(base, records, None)
        makeblastdb(binary, "-in", base + ".fa", "-dbtype", kind, "-out", base, "-parse_seqids", "-blastdb_version", "5",
                    "-title", "multidb " + base)
        os.remove(os.path.join(HERE, f"{base}.{kind[0]}js"))
        if kind == "prot":
            with open(os.path.join(HERE, "xq.faa"), "w") as f:
                f.write(">xq the random protein of pxm\n" + "".join(random_lines))
    lists = {"sl_pxm7.txt": "NP_900003.1\ngnl|BL_ORD_ID|7\n", "sl_pxm19.txt": "NP_900003.1\ngnl|BL_ORD_ID|19\n",
             "sl_nxm7.txt": "NM_900003.1\ngnl|BL_ORD_ID|7\n", "sl_ord7.txt": "gnl|BL_ORD_ID|7\n",
             "sl_ord15.txt": "gnl|BL_ORD_ID|15\n", "sl_ord19.txt": "gnl|BL_ORD_ID|19\n", "sl_np5.txt": "NP_900005.1\n"}
    for name, text in lists.items():
        with open(os.path.join(HERE, name), "w") as f:
            f.write(text)


def make_mixdb2(binary):
    """Builds put and writes the queries mix2_q*.f* (line DBID mixdb2)."""
    ids_q = os.path.join(os.path.dirname(HERE), "ids", "ids_q.faa")
    makeblastdb(binary, "-in", ids_q, "-dbtype", "prot", "-out", "put", "-blastdb_version", "5", "-taxid", "9606",
                "-title", "multidb put")
    os.remove(os.path.join(HERE, "put.pjs"))

    def letters(path, index=0):
        header, lines = read_fasta_path(path)[index]
        return header.rstrip("\n"), "".join(line.strip() for line in lines)

    qp1 = letters(ids_q)[1]
    glob1 = next("".join(line.strip() for line in lines) for header, lines in read_fasta_path(os.path.join(HERE, "nnp.fa"))
                 if header.split()[0] == ">FBtr0083245")
    bases = "TCAG"
    amino = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    code = {a + b + c: amino[16 * i + 4 * j + k] for i, a in enumerate(bases) for j, b in enumerate(bases)
            for k, c in enumerate(bases)}
    orf = ""
    for frame in range(3):
        protein = "".join(code.get(glob1[at:at + 3], "X") for at in range(frame, len(glob1) - 2, 3))
        for piece in protein.split("*"):
            start = piece.find("M")
            if start >= 0 and len(piece) - start > len(orf):
                orf = piece[start:]
    codon = {"A": "GCT", "R": "CGT", "N": "AAC", "D": "GAC", "C": "TGC", "Q": "CAG", "E": "GAA", "G": "GGT", "H": "CAC",
             "I": "ATC", "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTC", "P": "CCG", "S": "TCT", "T": "ACC", "W": "TGG",
             "Y": "TAC", "V": "GTT"}
    protein_q, nucleotide_q = letters(os.path.join(SRC, "q.faa")), letters(os.path.join(SRC, "q.fna"))
    for name, records in (
            ("mix2_q.faa", [(">mix2_u qp1 of ids/noacc_v5 (no accession table)", qp1), protein_q]),
            ("mix2_qx.fna", [(">mix2_u qp1 of ids/noacc_v5 back-translated", "".join(codon[a] for a in qp1)),
                             nucleotide_q]),
            ("mix2_q.fna", [(">mix2_u glob1 of multidb/nnp (no accession table)", glob1), nucleotide_q]),
            ("mix2_qt.faa", [(">mix2_u the glob1 ORF of multidb/nnp", orf), protein_q])):
        with open(os.path.join(HERE, name), "w") as f:
            for header, seq in records:
                f.write(header + "\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60)) + "\n")


def make_mixtax(binary):
    """Builds nut, putl and nutl and writes the queries mix2_qc*.f* and mix2_ql*.f* (line DBID mixtax)."""
    makeblastdb(binary, "-in", "nnp.fa", "-dbtype", "nucl", "-out", "nut", "-blastdb_version", "5", "-taxid", "9606",
                "-title", "multidb nut")
    os.remove(os.path.join(HERE, "nut.njs"))

    def records(name):
        out = []
        for header, lines in read_fasta_path(os.path.join(HERE, name)):
            out.append((header.rstrip("\n"), "".join(line.strip() for line in lines)))
        return out

    (u_header, qp1), _ = records("mix2_q.faa")
    (ux_header, qp1_x), _ = records("mix2_qx.fna")
    (un_header, glob1), _ = records("mix2_q.fna")
    (ut_header, orf), _ = records("mix2_qt.faa")
    globin = records(os.path.join(SRC, "q.faa"))[0][1]
    globin_x = records(os.path.join(SRC, "q.fna"))[0][1]
    qp2 = records(os.path.join(os.path.dirname(HERE), "ids", "ids_q.faa"))[1][1]
    codon = {"A": "GCT", "R": "CGT", "N": "AAC", "D": "GAC", "C": "TGC", "Q": "CAG", "E": "GAA", "G": "GGT", "H": "CAC",
             "I": "ATC", "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTC", "P": "CCG", "S": "TCT", "T": "ACC", "W": "TGG",
             "Y": "TAC", "V": "GTT"}
    qp2_x = "".join(codon[a] for a in qp2)
    for name, kind, header, seq in (("putl", "prot", ">qp2 unrelated protein", qp2),
                                    ("nutl", "nucl", ">qp2_x unrelated protein back-translated", qp2_x)):
        with open(os.path.join(HERE, name + ".fa"), "w") as f:
            f.write(header + "\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60)) + "\n")
        makeblastdb(binary, "-in", name + ".fa", "-dbtype", kind, "-out", name, "-blastdb_version", "5", "-taxid",
                    "99990", "-title", "multidb " + name)
        os.remove(os.path.join(HERE, name + (".pjs" if kind == "prot" else ".njs")))
    for name, pair in (
            ("mix2_qc.faa", [(u_header, qp1), (">mix2_c the letters of taxfilter/q.faa, then qp1", globin + qp1)]),
            ("mix2_qcx.fna", [(ux_header, qp1_x),
                              (">mix2_c the letters of taxfilter/q.fna, then qp1 back-translated", globin_x + qp1_x)]),
            ("mix2_ql.faa", [(u_header, qp1), (">mix2_l the letters of taxfilter/q.faa, then qp2", globin + qp2)]),
            ("mix2_qlx.fna", [(ux_header, qp1_x),
                              (">mix2_l the letters of taxfilter/q.fna, then qp2 back-translated", globin_x + qp2_x)]),
            ("mix2_ql.fna", [(un_header, glob1),
                             (">mix2_l the letters of taxfilter/q.fna, then qp2 back-translated", globin_x + qp2_x)]),
            ("mix2_qlt.faa", [(ut_header, orf), (">mix2_l the letters of taxfilter/q.faa, then qp2", globin + qp2)])):
        with open(os.path.join(HERE, name), "w") as f:
            for header, seq in pair:
                f.write(header + "\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60)) + "\n")


def make_mixsfx(binary):
    """Builds putb (line RG item T5c, merge 3): `put` with protein titles ending in `]`."""
    titles = ("alpha, partial [Homo sapiens]", "chi, partial [ ]", "MAG: (plastid) [ ]")
    records = read_fasta_path(os.path.join(os.path.dirname(HERE), "ids", "ids_q.faa"))
    with open(os.path.join(HERE, "putb.fa"), "w") as f:
        for title, (_, lines) in zip(titles, records, strict=True):
            seq = "".join(line.strip() for line in lines)
            f.write(">" + title + "\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60)) + "\n")
    makeblastdb(binary, "-in", "putb.fa", "-dbtype", "prot", "-out", "putb", "-blastdb_version", "5", "-taxid",
                "9606", "-title", "multidb putb")
    os.remove(os.path.join(HERE, "putb.pjs"))


def make_loadord(binary):
    """Builds ldp, lda, ldb, ldu and their nucleotide counterparts ldpn, ldan, ldbn, ldun (line DBID loadord)."""
    rng = random.Random(20261116)

    def letters(path, pick=lambda header: True):
        return next("".join(line.strip() for line in lines) for header, lines in read_fasta_path(path) if pick(header))

    def random_letters(alphabet, length):
        return "".join(rng.choice(alphabet) for _ in range(length))

    qp1 = letters(os.path.join(os.path.dirname(HERE), "ids", "ids_q.faa"))
    glob1 = letters(os.path.join(HERE, "nnp.fa"), lambda header: header.split()[0] == ">FBtr0083245")
    for kind, suffix, prefix, globin, unparsed, alphabet, length in (
            ("prot", "", "NP", letters(os.path.join(SRC, "q.faa")), qp1, "ACDEFGHIKLMNPQRSTVWY", 120),
            ("nucl", "n", "NM", letters(os.path.join(SRC, "q.fna")), glob1, "ACGT", 360)):
        sets = {
            "ldp": [(f">ref|{prefix}_000002.1| loadord packed number 2", globin),
                    (f">ref|{prefix}_000900.1| loadord random letters", random_letters(alphabet, length))],
            "lda": [(">loadord the unparsed hit at OID 0", unparsed)],
            "ldb": [(">loadord random letters", random_letters(alphabet, length)),
                    (">loadord the unparsed hit at OID 1", unparsed)],
            "ldu": [(">loadord the unparsed hit at OID 2", unparsed)],
            "ldg": [(f">gi|1|ref|{prefix}_000002.1| loadord the accession of ldp under a gi", globin)],
            "ldq": [(">lcl|Query_1 loadord a local id of the form of a generated query id", globin)],
        }
        for base, records in sets.items():
            name = base + suffix
            with open(os.path.join(HERE, name + ".fa"), "w") as f:
                for header, seq in records:
                    f.write(header + "\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60)) + "\n")
            parse = ["-parse_seqids"] if base in ("ldp", "ldg", "ldq") else []
            makeblastdb(binary, "-in", name + ".fa", "-dbtype", kind, "-out", name, *parse, "-blastdb_version", "5",
                        "-taxid", "9606", "-title", "multidb " + name)
            os.remove(os.path.join(HERE, name + (".pjs" if kind == "prot" else ".njs")))
    globin = letters(os.path.join(SRC, "q.faa"))
    with open(os.path.join(HERE, "ldo_qc.faa"), "w") as f:
        seq = globin + qp1
        f.write(">ldo_c the letters of taxfilter/q.faa, then qp1\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60))
                + "\n")


def make_loadord2(binary):
    """Builds ldc, ldh and their nucleotide counterparts ldcn, ldhn (line DBID loadord, DECISIONS 60 (6))."""

    def letters(path):
        return next("".join(line.strip() for line in lines) for header, lines in read_fasta_path(path))

    for kind, suffix, prefix, globin in (("prot", "", "NP", letters(os.path.join(SRC, "q.faa"))),
                                         ("nucl", "n", "NM", letters(os.path.join(SRC, "q.fna")))):
        sets = {
            "ldc": [(f">ref|{prefix}_000002.1| loadord a copy of ldp's first record", globin)],
            "ldh": [(f">gi|9000|ref|{prefix}_000002.1| loadord the accession of ldp under a gi after it", globin)],
        }
        for base, records in sets.items():
            name = base + suffix
            with open(os.path.join(HERE, name + ".fa"), "w") as f:
                for header, seq in records:
                    f.write(header + "\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60)) + "\n")
            makeblastdb(binary, "-in", name + ".fa", "-dbtype", kind, "-out", name, "-parse_seqids", "-blastdb_version",
                        "5", "-taxid", "9606", "-title", "multidb " + name)
            os.remove(os.path.join(HERE, name + (".pjs" if kind == "prot" else ".njs")))


def make_loadord3(binary):
    """Builds ldq2, ldx, ldd, ldmy and the nucleotide ldq2n, ldxn, lddn, the queries ldo_short.faa, ldo_short.fna and
    ldo_shortx.fna, and the PSSM ldo_my.asn (review 2 of line DBID loadord)."""
    import gzip
    rng = random.Random(20261116 + 3)

    def letters(path, pick=lambda header: True):
        return next("".join(line.strip() for line in lines) for header, lines in read_fasta_path(path) if pick(header))

    def write(name, records):
        with open(os.path.join(HERE, name), "w") as f:
            for header, seq in records:
                f.write(header + "\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60)) + "\n")

    qp1 = letters(os.path.join(os.path.dirname(HERE), "ids", "ids_q.faa"))
    glob1 = letters(os.path.join(HERE, "nnp.fa"), lambda header: header.split()[0] == ">FBtr0083245")
    qp1_x = letters(os.path.join(HERE, "mix2_qx.fna"))
    globin, globin_x = letters(os.path.join(SRC, "q.faa")), letters(os.path.join(SRC, "q.fna"))
    for kind, suffix, prefix, own, unparsed, alphabet, length in (
            ("prot", "", "NP", globin, qp1, "ACDEFGHIKLMNPQRSTVWY", 120),
            ("nucl", "n", "NM", globin_x, glob1, "ACGT", 360)):
        sets = {
            "ldq2": [(">lcl|Query_2 loadord a local id of the form of a second generated query id", own)],
            "ldx": [(f">gi|9001|ref|{prefix}_000002.1| loadord the accession of ldp under a gi, letters no query hits",
                     "".join(rng.choice(alphabet) for _ in range(length)))],
            "ldd": [(f">gi|9000|ref|{prefix}_000002.1| loadord the accession of ldp under a gi, the letters of lda", unparsed)],
        }
        if kind == "prot":
            sets["ldmy"] = [(">lcl|myprot loadord the local id of the PSSM query of ldo_my.asn", own)]
        for base, records in sets.items():
            name = base + suffix
            write(name + ".fa", records)
            makeblastdb(binary, "-in", name + ".fa", "-dbtype", kind, "-out", name, "-parse_seqids", "-blastdb_version",
                        "5", "-taxid", "9606", "-title", "multidb " + name)
            os.remove(os.path.join(HERE, name + (".pjs" if kind == "prot" else ".njs")))
    write("ldo_short.faa", [(">ldo_s the first 40 letters of qp1", qp1[:40]),
                            (">ldo_g the letters of taxfilter/q.faa", globin)])
    write("ldo_short.fna", [(">ldo_s the first 120 letters of glob1", glob1[:120]),
                            (">ldo_g the letters of taxfilter/q.fna", globin_x)])
    write("ldo_shortx.fna", [(">ldo_s the first 120 letters of qp1 back-translated", qp1_x[:120]),
                             (">ldo_g the letters of taxfilter/q.fna", globin_x)])
    with gzip.open(os.path.join(os.path.dirname(os.path.dirname(HERE)), "fixtures", "pssm", "q1_tfp.asn.gz"), "rt") as f:
        pssm = f.read()
    assert pssm.count('local str "Query_1"') == 1
    with open(os.path.join(HERE, "ldo_my.asn"), "w") as f:
        f.write(pssm.replace('local str "Query_1"', 'local str "myprot"'))


def make_loadord4(binary):
    """Builds ldr, lds, ldi and their nucleotide counterparts ldrn, ldsn, ldin (review 3 of line DBID loadord)."""

    def letters(path):
        return next("".join(line.strip() for line in lines) for header, lines in read_fasta_path(path))

    for kind, suffix, prefix, globin in (("prot", "", "NP", letters(os.path.join(SRC, "q.faa"))),
                                         ("nucl", "n", "NM", letters(os.path.join(SRC, "q.fna")))):
        sets = {
            "ldr": [(f">ref|{prefix}_000001.1| loadord an accession of packed number 1, before ldp's", globin)],
            "lds": [(f">ref|{prefix}_000003.1| loadord an accession of packed number 3, after ldp's", globin)],
            "ldi": [(f">gi|1|ref|{prefix}_000003.1| loadord the accession of lds under a gi that loads first", globin)],
        }
        for base, records in sets.items():
            name = base + suffix
            with open(os.path.join(HERE, name + ".fa"), "w") as f:
                for header, seq in records:
                    f.write(header + "\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60)) + "\n")
            makeblastdb(binary, "-in", name + ".fa", "-dbtype", kind, "-out", name, "-parse_seqids", "-blastdb_version",
                        "5", "-taxid", "9606", "-title", "multidb " + name)
            os.remove(os.path.join(HERE, name + (".pjs" if kind == "prot" else ".njs")))


def make_loadord5(binary):
    """Builds l5c, l5d, l5e and their nucleotide counterparts l5cn, l5dn, l5en (review 4 of line DBID loadord, H1)."""

    def letters(path):
        return next("".join(line.strip() for line in lines) for header, lines in read_fasta_path(path))

    for kind, suffix, prefix, letter, globin in (("prot", "", "NP", "p", letters(os.path.join(SRC, "q.faa"))),
                                                 ("nucl", "n", "NM", "n", letters(os.path.join(SRC, "q.fna")))):
        for base in ("l5c", "l5d", "l5e"):
            name = base + suffix
            with open(os.path.join(HERE, name + ".fa"), "w") as f:
                f.write(f">ref|{prefix}_000002.1| loadord a copy of ldp's first record, version 4\n"
                        + "\n".join(globin[at:at + 60] for at in range(0, len(globin), 60)) + "\n")
            makeblastdb(binary, "-in", name + ".fa", "-dbtype", kind, "-out", name, "-parse_seqids", "-blastdb_version",
                        "4", "-taxid", "9606", "-title", "multidb " + name)
            for extra in (".pjs", ".njs"):
                path = os.path.join(HERE, name + extra)
                if os.path.exists(path):
                    os.remove(path)
        # l5d's string ISAM index: a directory, which CSeqDBIsam::IndexExists (CFile::Exists, a regular file) does not
        # take for one, so that NCBI's lookups find nothing in that volume.
        index = os.path.join(HERE, f"l5d{suffix}.{letter}si")
        os.remove(index)
        os.mkdir(index)
        with open(os.path.join(index, "placeholder"), "w") as f:
            f.write("A directory in place of the string ISAM index (make_multidb.py --loadord5-only).\n")


def make_loadord6(binary):
    """Builds ldl, ldln, ldyn and the queries ldo_e2.faa, ldo_e2x.fna, ldo_e2.fna, ldo_e2t.faa and ldo_e2y.fna (review 5
    of line DBID loadord)."""

    def letters(path, pick=lambda header: True):
        return next("".join(line.strip() for line in lines) for header, lines in read_fasta_path(path) if pick(header))

    def write(name, records):
        with open(os.path.join(HERE, name), "w", encoding="utf-8") as f:
            for header, seq in records:
                f.write(header + "\n" + "\n".join(seq[at:at + 60] for at in range(0, len(seq), 60)) + "\n")

    qp1 = letters(os.path.join(os.path.dirname(HERE), "ids", "ids_q.faa"))
    glob1 = letters(os.path.join(HERE, "nnp.fa"), lambda header: header.split()[0] == ">FBtr0083245")
    qp1_x = letters(os.path.join(HERE, "mix2_qx.fna"))
    orf = letters(os.path.join(HERE, "mix2_qt.faa"))
    globin, globin_x = letters(os.path.join(SRC, "q.faa")), letters(os.path.join(SRC, "q.fna"))
    for kind, suffix, unparsed in (("prot", "", qp1), ("nucl", "n", glob1)):
        name = "ldl" + suffix
        write(name + ".fa", [(">loadord the unparsed hit at OID 0, taxid 99990", unparsed),
                             (">loadord the unparsed hit at OID 1, taxid 99990", unparsed)])
        makeblastdb(binary, "-in", name + ".fa", "-dbtype", kind, "-out", name, "-blastdb_version", "5", "-taxid",
                    "99990", "-title", "multidb " + name)
        os.remove(os.path.join(HERE, name + (".pjs" if kind == "prot" else ".njs")))
    rng = random.Random(20261116 + 5)
    other = "".join(rng.choice("ACGT") for _ in range(360))
    write("ldyn.fa", [(f">loadord random letters that only q1 of ldo_e2y.fna hits, copy {n}", other) for n in range(4)])
    makeblastdb(binary, "-in", "ldyn.fa", "-dbtype", "nucl", "-out", "ldyn", "-blastdb_version", "5", "-taxid", "9606",
                "-title", "multidb ldyn")
    os.remove(os.path.join(HERE, "ldyn.njs"))
    none_p, none_n = "WWWWWWWWWWHHHHHHHHHHWWWWWWWWWW", "ACGTTGCAAGTCCTAGGATCCGTAGCATGCAATTGGCCAATCGATCGTTAACCGGTACG"
    cut = ">q1 a\u00e9b z"
    for name, nohit, repeated, two in (("ldo_e2.faa", none_p, qp1 * 16, globin + qp1),
                                       ("ldo_e2x.fna", none_n, qp1_x * 16, globin_x + qp1_x),
                                       ("ldo_e2.fna", none_n, glob1 * 6, globin_x + glob1),
                                       ("ldo_e2t.faa", none_p, orf * 16, globin + orf),
                                       ("ldo_e2y.fna", none_n, glob1 + other, globin_x + glob1)):
        write(name, [(">q0 none", nohit), (cut, repeated), (">q2 two", two)])


def make_loadord7(binary):
    """Builds ldv (review 6 of line DBID loadord, H1)."""
    qp1 = next("".join(line.strip() for line in lines)
               for header, lines in read_fasta_path(os.path.join(os.path.dirname(HERE), "ids", "ids_q.faa")))
    rng = random.Random(20261116 + 7)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    variant = "".join(rng.choice(alphabet.replace(letter, "")) if rng.random() < 0.1 else letter for letter in qp1)
    with open(os.path.join(HERE, "ldv.fa"), "w") as f:
        f.write(">gi|9002|ref|NP_000002.1| loadord the accession of ldp under a gi, a variant of lda's letters\n"
                + "\n".join(variant[at:at + 60] for at in range(0, len(variant), 60)) + "\n")
    makeblastdb(binary, "-in", "ldv.fa", "-dbtype", "prot", "-out", "ldv", "-parse_seqids", "-blastdb_version", "5",
                "-taxid", "9606", "-title", "multidb ldv")
    os.remove(os.path.join(HERE, "ldv.pjs"))


def make_v4late(binary):
    """Builds l6o, l6u and l5b (line DBID v4late)."""
    qp1 = next("".join(line.strip() for line in lines)
               for header, lines in read_fasta_path(os.path.join(os.path.dirname(HERE), "ids", "ids_q.faa")))
    for name, header in (("l6o", ">ref|NP_000002.1| v4late the accession of l5c with other letters, version 4"),
                         ("l6u", ">ref|NP_000004.1| v4late an accession that only this volume holds, version 4"),
                         ("l5b", ">ref|NP_000999.1|ref|NP_000002.1| v4late the accession of l5c behind another one, "
                                 "version 4")):
        with open(os.path.join(HERE, name + ".fa"), "w") as f:
            f.write(header + "\n" + "\n".join(qp1[at:at + 60] for at in range(0, len(qp1), 60)) + "\n")
        makeblastdb(binary, "-in", name + ".fa", "-dbtype", "prot", "-out", name, "-parse_seqids", "-blastdb_version",
                    "4", "-taxid", "9606", "-title", "multidb " + name)
        os.remove(os.path.join(HERE, name + ".pjs"))


def longest_orf(letters):
    """The longest open reading frame (from an ATG to the stop) of the plus strand of `letters`, translated with the
    standard genetic code, without the stop."""
    bases = "TCAG"
    amino = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    code = {a + b + c: amino[16 * i + 4 * j + k] for i, a in enumerate(bases) for j, b in enumerate(bases)
            for k, c in enumerate(bases)}
    best = ""
    for frame in range(3):
        protein = "".join(code.get(letters[at:at + 3], "X") for at in range(frame, len(letters) - 2, 3))
        for piece in protein.split("*"):
            start = piece.find("M")
            if start >= 0 and len(piece) - start > len(best):
                best = piece[start:]
    return best


def make_pfoff(_binary):
    """Writes the queries pfo_q*.f* (line DBID pfoff): those of mix2_q*.f* in the other order, pfo_qf.fna and
    pfo_qft.faa (na's HBB record, its ORF, then the query that hits nnp alone), and pfo_qs*.f* (the rows .taxstop.)."""
    for source, name in (("mix2_q.faa", "pfo_q.faa"), ("mix2_qx.fna", "pfo_qx.fna"), ("mix2_q.fna", "pfo_q.fna"),
                         ("mix2_qt.faa", "pfo_qt.faa")):
        records = read_fasta_path(os.path.join(HERE, source))
        with open(os.path.join(HERE, name), "w") as f:
            for header, lines in reversed(records):
                f.write(header + "".join(lines))
    hbb = next(lines for header, lines in read_fasta_path(os.path.join(HERE, "na.fa"))
               if header.startswith(">ENST00000335295.4 "))
    letters = "".join(line.strip() for line in hbb).upper()
    for source, name, header, query in (
            ("mix2_q.fna", "pfo_qf.fna", ">pfo_f the letters of multidb/na's ENST00000335295.4 (HBB)", letters),
            ("mix2_qt.faa", "pfo_qft.faa", ">pfo_ft the longest ORF of multidb/na's ENST00000335295.4 (HBB)",
             longest_orf(letters))):
        unparsed = next((h, lines) for h, lines in read_fasta_path(os.path.join(HERE, source))
                        if h.startswith(">mix2_u "))
        with open(os.path.join(HERE, name), "w") as f:
            f.write(header + "\n")
            for at in range(0, len(query), 60):
                f.write(query[at:at + 60] + "\n")
            f.write(unparsed[0] + "".join(unparsed[1]))

    def letters_of(source, unparsed):
        """The letters of the query of `source` that hits nnp/noacc_v5 alone (`unparsed`) or of the other one."""
        for header, lines in read_fasta_path(os.path.join(HERE, source)):
            if header.startswith(">mix2_u ") == unparsed:
                return "".join(line.strip() for line in lines)

    def write_queries(name, queries):
        with open(os.path.join(HERE, name), "w") as f:
            for header, query in queries:
                f.write(header + "\n")
                for at in range(0, len(query), 60):
                    f.write(query[at:at + 60] + "\n")

    for source, name in (("mix2_q.faa", "pfo_qs.faa"), ("mix2_qx.fna", "pfo_qsx.fna")):
        parsed, unparsed = letters_of(source, False), letters_of(source, True)
        first = parsed + unparsed[:len(unparsed) // 2]
        write_queries(name, ((">pfo_s1 twice the letters of the parsed query of " + source + ", then half of the other",
                              first + first),
                             (">pfo_s2 the letters of the other query of " + source + ", then those of the parsed one",
                              unparsed + parsed)))
    for source, name, hbb in (("mix2_q.fna", "pfo_qs.fna", letters), ("mix2_qt.faa", "pfo_qst.faa", longest_orf(letters))):
        parsed, unparsed = letters_of(source, False), letters_of(source, True)
        write_queries(name, ((">pfo_s1 twice multidb/na's ENST00000335295.4 (HBB" + (", its ORF" if hbb != letters else "")
                              + "), then half of the query of " + source + " that hits nnp", hbb + hbb +
                              unparsed[:len(unparsed) // 2]),
                             (">pfo_s2 the letters of the query of " + source + " that hits nnp, then those of the other",
                              unparsed + parsed)))
    parsed, unparsed = letters_of("mix2_q.fna", False), letters_of("mix2_q.fna", True)
    write_queries("pfo_qss.fna", ((">pfo_s1 twice multidb/na's ENST00000335295.4 (HBB), then half of the query of "
                                   "mix2_q.fna that hits nnp", letters + letters + unparsed[:len(unparsed) // 2]),
                                  (">pfo_s3 a third of the query of mix2_q.fna that hits nnp, then 200 letters of the "
                                   "other", unparsed[:len(unparsed) // 3] + parsed[:200])))


def make_pfoff2(binary):
    """Builds pfz and pfzq and writes pfo_qz.faa and pfo_qh.fna (line DBID pfoff-2)."""
    zero = "WWWWWWWWWWHHKYEVIMIKLVMFVIEEIFIMEVLIFLVLLEVFVMYREYMILRMYILLFKIMILVLFDMEKHHWWWWWWWWWW"
    some = "CCCCCCCCCCRFDIVLVRMILYIVDCCCCCCCCCC"
    for name, records in (("pfz", ((">ref|XP_900001.1| pfoff-2 an alignment of no identities", zero),
                                   (">ref|XP_900002.1| pfoff-2 an alignment of identities", some))),
                          ("pfzq", ((">ref|XP_900001.1| pfoff-2 an alignment of no identities", zero),))):
        with open(os.path.join(HERE, name + ".fa"), "w") as f:
            for header, letters in records:
                f.write(header + "\n" + letters + "\n")
        makeblastdb(binary, "-in", name + ".fa", "-dbtype", "prot", "-out", name, "-parse_seqids", "-blastdb_version",
                    "5", "-title", "multidb " + name)
        os.remove(os.path.join(HERE, name + ".pjs"))
    with open(os.path.join(HERE, "pfo_qz.faa"), "w") as f:
        f.write(">pfo_z a query whose alignment to pfz's XP_900001.1 has no identities\n"
                "PPPPPPPPPPGGRFDIVLVRMILYIVDDVYVLDIMVYMIMMDIYILFKDFLVMKLFVMMYRVLVMIMYELDRGGPPPPPPPPPP\n")
    hbb = next(lines for header, lines in read_fasta_path(os.path.join(HERE, "na.fa"))
               if header.startswith(">ENST00000335295.4 "))
    letters = "".join(line.strip() for line in hbb).upper()
    unparsed = next("".join(line.strip() for line in lines)
                    for header, lines in read_fasta_path(os.path.join(HERE, "mix2_q.fna"))
                    if header.startswith(">mix2_u "))
    with open(os.path.join(HERE, "pfo_qh.fna"), "w") as f:
        for header, query in ((">pfo_h1 the letters of multidb/na's ENST00000335295.4 (HBB)", letters),
                              (">pfo_h2 twice the letters of the query of mix2_q.fna that hits nnp", unparsed + unparsed)):
            f.write(header + "\n")
            for at in range(0, len(query), 60):
                f.write(query[at:at + 60] + "\n")


def make_pfoff3(binary):
    """Builds pf3p and pf3u and writes pfo_qr.faa, pfo_qe.faa and pfo_qm.faa (line DBID pfoff-3, reviews 1 to 3 of the
    item)."""
    query = ("IWVCGTEYVGHAEKMQHNPTFNQPWRLSAHCDYAQSMWQKPMTSECEELMWDMVKHMFIVHQEPSSGAHKPKKKQQLQGPWWWINPNWRWMRIPNIQTLTVIHD"
             "WWLNCHASYYTWFRYNVTKVDQQACCEEEYPPIRRFPWWHEELGPY")
    copies = ("MRIPWIQTLTVIHDWWLNCHASYYTWFRYNVTKVDQQACCEEENPPIRRFPWWHEELGPYHGAMVYYLNPTPLFSCGQYPERIDYRRPDTIWVCGTEYVVHAEKMQ"
              "HNPTFNQPWRLSAHCDYAQSMWQRPMTSECEELMWDTVKHMFIV")
    weak = "EECSRDAVFYDGYKRSAFHRVWLECKSNDYVGCHMWEFSFTCFFPQSFGIAPHPATFVEGMFCSEVNCQWLDWCGREEGCAQLKRIHVHKWAKTM"
    for name, header, letters, options in (
            ("pf3p", ">ref|XP_810001.1| pfoff-3 two copies of segments of pfo_qr", copies, ["-parse_seqids"]),
            ("pf3u", ">bsub pfoff-3 a subject with a taxid", weak, ["-taxid", "9606"])):
        with open(os.path.join(HERE, name + ".fa"), "w") as f:
            f.write(header + "\n")
            for at in range(0, len(letters), 60):
                f.write(letters[at:at + 60] + "\n")
        makeblastdb(binary, "-in", name + ".fa", "-dbtype", "prot", "-out", name, *options, "-blastdb_version", "5",
                    "-title", "multidb " + name)
        os.remove(os.path.join(HERE, name + ".pjs"))
    with open(os.path.join(HERE, "pfo_qr.faa"), "w") as f:
        f.write(">pfo_r a query whose second psiblast round shows two HSPs of pf3p before pf3u's bsub\n")
        for at in range(0, len(query), 60):
            f.write(query[at:at + 60] + "\n")
    # Review 3 of the item, N1: pfo_r, whose report ends the run where a batch's prefetch is off, then bsub's letters,
    # whose batch's prefetch is on and loads bsub (its taxid's lookup), which NCBI never makes.
    with open(os.path.join(HERE, "pfo_qe.faa"), "w") as f:
        for header, letters in ((">pfo_r a query whose report ends the run", query),
                                (">pfo_e bsub's letters, a later batch whose prefetch is on", weak)):
            f.write(header + "\n")
            for at in range(0, len(letters), 60):
                f.write(letters[at:at + 60] + "\n")
    records = []
    for name, which in (("mix2_qc.faa", 1), ("pfo_q.faa", 0)):
        with open(os.path.join(HERE, name)) as f:
            parts = [[]]
            for line in f:
                if line.startswith(">"):
                    parts.append([])
                parts[-1].append(line)
        records.append("".join(parts[which + 1]))
    with open(os.path.join(HERE, "pfo_qm.faa"), "w") as f:
        f.write("".join(records))


# One codon for each residue, of the back-translations of `make_rf1`.
RF1_CODONS = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA",
              "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACC", "V": "GTT",
              "W": "TGG", "Y": "TAT"}


def make_rf1(binary):
    """Builds n3p and n3u from pf3p.fa and pf3u.fa (`make_pfoff3`) and writes the rf1 queries (line DBID rf-1)."""
    def letters(path, index=0):
        return "".join(line.strip() for line in read_fasta_path(os.path.join(HERE, path))[index][1])

    def back(protein):
        return "".join(RF1_CODONS[residue] for residue in protein)

    def write(name, records):
        with open(os.path.join(HERE, name), "w") as f:
            for header, sequence in records:
                f.write(header + "\n")
                for at in range(0, len(sequence), 60):
                    f.write(sequence[at:at + 60] + "\n")

    copies, weak, query = letters("pf3p.fa"), letters("pf3u.fa"), letters("pfo_qe.faa")
    for name, header, protein, options in (
            ("n3p", ">ref|XM_810001.1| rf-1 pf3p's letters back-translated", copies, ["-parse_seqids"]),
            ("n3u", ">bsubn rf-1 pf3u's bsub back-translated, a subject with a taxid", weak, ["-taxid", "9606"])):
        write(name + ".fa", [(header, back(protein))])
        makeblastdb(binary, "-in", name + ".fa", "-dbtype", "nucl", "-out", name, *options, "-blastdb_version", "5",
                    "-title", "multidb " + name)
        os.remove(os.path.join(HERE, name + ".njs"))
    q0_protein = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQQIAAALE"
    q0_nucleotide = ("ATGCGTACGTAGCTAGCTAGCTAGCTGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAG"
                     "CTAGCTAGCTAGCTAGCTAGCTAG")
    q0 = ">q0_nohit a query without hits"
    e = ">pfo_e bsub's letters, a later batch whose prefetch loads bsub"
    r = ">pfo_r a query whose report shows pf3p's hit"
    write("rf1_q0.faa", [(q0, q0_protein)])
    write("rf1_q0e.faa", [(q0, q0_protein), (e, weak)])
    write("rf1_e.faa", [(e, weak)])
    write("rf1_r.faa", [(r, query)])
    write("rf1_re.faa", [(r, query), (e, weak)])
    write("rf1_q0.fna", [(q0, q0_nucleotide)])
    write("rf1_q0e.fna", [(q0, q0_nucleotide), (e, back(weak))])
    write("rf1_r.fna", [(r, back(query))])
    write("rf1_re.fna", [(r, back(query)), (e, back(weak))])
    write_rf1_title(write, (q0, q0_protein))


def write_rf1_title(write, q0):
    """Writes rf1_q0t.faa (rf-1 R1): `q0_nohit`, then taxfilter/q.faa's P69905, whose report shows taxfilter's `tfe_nr`
    OID 3, a sequence with an empty stored title and taxid 9606, whose title NCBI makes from the name of the taxid."""
    header, lines = read_fasta_path(os.path.join(HERE, "..", "taxfilter", "q.faa"))[0]
    write("rf1_q0t.faa", [q0, (header.rstrip("\n"), "".join(line.strip() for line in lines))])


def make_rf2(binary):
    """Builds rf2p and rf2n and writes rf2_q0t.fna (line DBID rf-2)."""
    for name, dbtype, path, header in (("rf2p", "prot", "tfe.faa", ">ref|XP_820003.1|"),
                                       ("rf2n", "nucl", "tfe.fna", ">ref|XM_820003.1|")):
        letters = "".join(line.strip() for line in read_fasta(path)[3][1])
        with open(os.path.join(HERE, name + ".fa"), "w") as f:
            f.write(header + "\n")
            for at in range(0, len(letters), 60):
                f.write(letters[at:at + 60] + "\n")
        makeblastdb(binary, "-in", name + ".fa", "-dbtype", dbtype, "-out", name, "-parse_seqids", "-taxid", "9606",
                    "-blastdb_version", "5", "-title", "multidb " + name)
        os.remove(os.path.join(HERE, name + (".pjs" if dbtype == "prot" else ".njs")))
    with open(os.path.join(HERE, "rf2_q0t.fna"), "w") as f:
        for path in (os.path.join(HERE, "rf1_q0.fna"), os.path.join(SRC, "q.fna")):
            with open(path) as part:
                f.write(part.read())
    write_rf1_rps()


def write_rf1_rps():
    """Writes rf1_rps.faa and rf1_rps.fna (line DBID rf-1, review 1 L3): `q0_nohit` (rf1_q0.faa, rf1_q0.fna), then
    ../../fasta/rps/q649.faa's YP_009220649.1 (in rf1_rps.fna back-translated with RF1_CODONS), which hits a profile of
    each volume of rpsblast's rps_mini3v, whose taxids (tx3v00 to tx3v02: 9606, 10090, 562) the loads look up."""
    q0 = read_fasta_path(os.path.join(HERE, "rf1_q0.faa"))[0]
    q0n = read_fasta_path(os.path.join(HERE, "rf1_q0.fna"))[0]
    header, lines = read_fasta_path(os.path.join(HERE, "..", "..", "fasta", "rps", "q649.faa"))[0]
    protein = "".join(line.strip() for line in lines)
    for name, first, letters in (("rf1_rps.faa", q0, protein),
                                 ("rf1_rps.fna", q0n, "".join(RF1_CODONS[residue] for residue in protein))):
        with open(os.path.join(HERE, name), "w") as f:
            f.write(first[0] + "".join(first[1]))
            f.write(header)
            for at in range(0, len(letters), 60):
                f.write(letters[at:at + 60] + "\n")


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
    make_pgim(binary, tfp)
    for base, kind in (("n4a", "nucl"), ("n4b", "nucl"), ("p4a", "prot"), ("p4b", "prot")):
        makeblastdb(binary, "-in", base + ".fa", "-dbtype", kind, "-out", base, "-blastdb_version", "4",
                    "-title", "multidb " + base)
    for d in ("d1", "d2"):
        makeblastdb(binary, "-in", f"{d}/g.fa", "-dbtype", "nucl", "-out", f"{d}/g", "-parse_seqids",
                    "-blastdb_version", "5", "-title", f"multidb {d} g")
    makeblastdb(binary, "-in", "et.fa", "-dbtype", "nucl", "-out", "et", "-parse_seqids", "-blastdb_version", "5",
                "-title", "")
    for name in ("pa.pjs", "pb.pjs", "p4a.pjs", "p4b.pjs", "na.njs", "nb.njs", "nnp.njs", "n4a.njs", "n4b.njs",
                 "d1/g.njs", "d2/g.njs", "et.njs", "pgim.pjs"):
        path = os.path.join(HERE, name)
        if os.path.exists(path):
            os.remove(path)
    with open(os.path.join(HERE, "al.nal"), "w") as f:
        f.write("TITLE multidb alias of na\nDBLIST na\n")
    with open(os.path.join(HERE, "mix_u.nal"), "w") as f:
        f.write("TITLE multidb mix_u: the sequences of nnp in a list with na\nDBLIST na nnp\nFIRST_OID 9\nLAST_OID 17\n")
    firsts = [accessions(sets[b][0][0][0])[0] for b in ("na", "nb", "pa", "pb")]
    with open(os.path.join(HERE, "sl.txt"), "w") as f:
        f.write("".join(acc + "\n" for acc in firsts))


if __name__ == "__main__":
    if sys.argv[1] == "--ordsl-only":
        make_ordsl(sys.argv[2])
    elif sys.argv[1] == "--mixdb2-only":
        make_mixdb2(sys.argv[2])
    elif sys.argv[1] == "--mixtax-only":
        make_mixtax(sys.argv[2])
    elif sys.argv[1] == "--mixsfx-only":
        make_mixsfx(sys.argv[2])
    elif sys.argv[1] == "--loadord-only":
        make_loadord(sys.argv[2])
    elif sys.argv[1] == "--loadord2-only":
        make_loadord2(sys.argv[2])
    elif sys.argv[1] == "--loadord3-only":
        make_loadord3(sys.argv[2])
    elif sys.argv[1] == "--loadord4-only":
        make_loadord4(sys.argv[2])
    elif sys.argv[1] == "--loadord5-only":
        make_loadord5(sys.argv[2])
    elif sys.argv[1] == "--loadord6-only":
        make_loadord6(sys.argv[2])
    elif sys.argv[1] == "--loadord7-only":
        make_loadord7(sys.argv[2])
    elif sys.argv[1] == "--v4late-only":
        make_v4late(sys.argv[2])
    elif sys.argv[1] == "--pfoff-only":
        make_pfoff(sys.argv[2])
    elif sys.argv[1] == "--pfoff2-only":
        make_pfoff2(sys.argv[2])
    elif sys.argv[1] == "--pfoff3-only":
        make_pfoff3(sys.argv[2])
    elif sys.argv[1] == "--rf1-only":
        make_rf1(sys.argv[2])
    elif sys.argv[1] == "--rf2-only":
        make_rf2(sys.argv[2])
    else:
        main(sys.argv[1])
        make_ordsl(sys.argv[1])
        make_mixdb2(sys.argv[1])
        make_mixtax(sys.argv[1])
        make_mixsfx(sys.argv[1])
        make_loadord(sys.argv[1])
        make_loadord2(sys.argv[1])
        make_loadord3(sys.argv[1])
        make_loadord4(sys.argv[1])
        make_loadord5(sys.argv[1])
        make_loadord6(sys.argv[1])
        make_loadord7(sys.argv[1])
        make_v4late(sys.argv[1])
        make_pfoff(sys.argv[1])
        make_pfoff2(sys.argv[1])
        make_pfoff3(sys.argv[1])
        make_rf1(sys.argv[1])
        make_rf2(sys.argv[1])
