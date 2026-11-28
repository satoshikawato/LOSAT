# `db_search_long/`: frozen rows of the long-sequence searches

`dbmask_long.<task>.6.out`: NCBI BLAST+ 2.17.0's rows of `tests/blastn_db_soft_long.rs` (series NP, M4a-2a, review 1
T4): `-db_soft_mask` on a 5.6 Mbp database sequence that the preliminary stage cuts into two chunks, with unsorted,
overlapping and nested masks and masks across the chunk boundary. The test writes its inputs (`long.fna`,
`long_masks.asnt`, `lq.fna`) with a fixed PRNG; `DBMASK_LONG_DIR=<dir>` keeps them in `<dir>`. Frozen from those
inputs with NCBI's own makeblastdb and blastn (`BLAST_USAGE_REPORT=false`):

```
makeblastdb -in long.fna -dbtype nucl -out long -mask_data long_masks.asnt
blastn -query lq.fna -db long -db_soft_mask 11 -task megablast -outfmt 6 -num_threads 1 > dbmask_long.megablast.6.out
blastn -query lq.fna -db long -db_soft_mask 11 -task dc-megablast -outfmt 6 -num_threads 1 > dbmask_long.dc-megablast.6.out
blastn -query lq.fna -db long -db_soft_mask 11 -task blastn -outfmt 6 -num_threads 1 > dbmask_long.blastn.6.out
blastn -query lq.fna -db long -db_soft_mask 11 -task blastn-short -evalue 1e-5 -outfmt 6 -num_threads 1 > dbmask_long.blastn-short.6.out
```

`dbmask_vols.<task>.6.out`: NCBI BLAST+ 2.17.0's rows of the second test of `tests/blastn_db_soft_long.rs` (review 2 T1):
`-db_soft_mask` on a database of two volumes (`-max_file_sz 64000`) of 400 sequences of about 1,000 letters, most
with up to three masks (some overlapping the one before), searched in many blocks at 1, 4 and 8 threads. The test
writes `vols.fna`, `vols_masks.asnt` and `vq.fna`; `DBMASK_VOLS_DIR=<dir>` keeps them. NCBI's own makeblastdb gives the
same `.nsq`, `.nhr`, `.nab` and `.nac` files (the dates and the title differ). Frozen with:

```
makeblastdb -in vols.fna -dbtype nucl -out vols -max_file_sz 64000 -mask_data vols_masks.asnt
blastn -query vq.fna -db vols -db_soft_mask 11 -task megablast -outfmt 6 -num_threads 1 > dbmask_vols.megablast.6.out
blastn -query vq.fna -db vols -db_soft_mask 11 -task dc-megablast -outfmt 6 -num_threads 1 > dbmask_vols.dc-megablast.6.out
blastn -query vq.fna -db vols -db_soft_mask 11 -task blastn -outfmt 6 -num_threads 1 > dbmask_vols.blastn.6.out
blastn -query vq.fna -db vols -db_soft_mask 11 -task blastn-short -evalue 1e-5 -outfmt 6 -num_threads 1 > dbmask_vols.blastn-short.6.out
```

NCBI's rows at `-num_threads 8` are the same. Without `-db_soft_mask` NCBI finds 24, 24, 115 and 24 rows (masked: 14, 14,
94 and 24): the queries over a first mask with 10 letters on each side lose their seeds.

`dbmask_hard_long.<task>.6.out`: NCBI BLAST+ 2.17.0's rows (`-outfmt "6 std nident"`) of `tests/blastn_db_hard_long.rs`
(series NP, M4a-2b, review T3 of ad-M4a-6): `-db_hard_mask` on two sequences longer than NCBI's chunk, `hard1` (5,300,000
letters, the one-letter mask 1001: its last stretch is cut into a full chunk and a final chunk, review F3) and `hard2`
(5,400,000 letters, masks 1,000-1,099 and 5,100,002-5,100,049), and `hard3` (3,000 letters, adjacent, one-letter and
three-letter masks and a mask that runs to its end after them), searched at 1 and 4 threads. The test writes `hl.fna`,
`hl_masks.asnt` and `hlq.fna`; `DBMASK_HARD_LONG_DIR=<dir>` keeps them. Frozen with:

```
makeblastdb -in hl.fna -dbtype nucl -out hl -mask_data hl_masks.asnt
blastn -query hlq.fna -db hl -db_hard_mask 11 -task megablast -outfmt "6 std nident" -num_threads 1 > dbmask_hard_long.megablast.6.out
blastn -query hlq.fna -db hl -db_hard_mask 11 -task dc-megablast -outfmt "6 std nident" -num_threads 1 > dbmask_hard_long.dc-megablast.6.out
blastn -query hlq.fna -db hl -db_hard_mask 11 -task blastn -outfmt "6 std nident" -num_threads 1 > dbmask_hard_long.blastn.6.out
blastn -query hlq.fna -db hl -db_hard_mask 11 -task blastn-short -evalue 1e-5 -outfmt "6 std nident" -num_threads 1 > dbmask_hard_long.blastn-short.6.out
```

`dbmask_prot_hard_long.<program>.<mode>.6.out`: NCBI BLAST+ 2.17.0's rows (`-outfmt "6 std nident positive btop"`) of `tests/protein_db_hard_long.rs` (series NP, M4b-3d): `-db_hard_mask 21` on five protein sequences: `phl_lr_mid` (5,000,300 aa, masks `[10, 20)` and `[5_000_250, 5_000_280)` — a > 5,000,000 hard range in the middle), `phl_lr_po` (5,000,300 aa, masks `[620, 650)`, `[50, 80)`, `[5_000_250, 5_000_280)` — out-of-order masks combined with a > 5,000,000 hard range), `phl_lr_plain` (5,000,300 aa, no masks — unmasked > 5,000,000 sibling), `phl_po` (400 aa, out-of-order masks `[250, 280)`, `[100, 130)`, `[320, 340)`), and `phl_ow` (240 aa, overlapping masks with a one-word short chunk whose HSP merges with a later chunk). The test writes `phl.faa`, `phl_masks.asnt`, `phlq.faa`, and `phlq.fna`; `DBMASK_PROT_HARD_LONG_DIR=<dir>` keeps them. Frozen from official NCBI BLAST+ 2.17.0 (`-num_threads 1`, matching `-num_threads 4`). A second test in the same file checks decision B-d on `overlap_past_list`: `po_last` (static `LastRangeOverlap` refusal), `po_ow_last` (dynamic `last_range_overlap_merge_refusal`), and `lr_last` (a > 5,000,000-residue final hard range merging within its own range with `DBSEQ_CHUNK_OVERLAP = 100`).

`phi_long.<case>.6.out`: NCBI BLAST+ 2.17.0's rows of `tests/psiblast_phi_long.rs` (line PS, item PS8a-2): psiblast
`-phi_pattern` (`W-C-x(2)-H-W`) on subjects longer than NCBI's part of 5,000,000 residues, searched in overlapping
parts: a pattern where the second part starts (`phl_5m1`), in the overlap of two parts (`phl_ov`, its hit counted in
both parts, which every E-value reads), across 5,000,000 (`phl_c2`), in both overlaps of three parts (`phl_3`) and the
two places of the second query in an overlap (`phl_2`), with two short subjects, and `-sum_stats true` with
`-evalue 1e-138` (the linked E-values free HSPs that the run without sum statistics keeps). The test writes `phl_s.faa`,
`phl_q.faa`, `phl_q2.faa` and `phl.pat` with a fixed generator; `PHI_LONG_DIR=<dir>` keeps them. Frozen with NCBI's own
makeblastdb and psiblast (`F="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
score nident"`, `P="-phi_pattern phl.pat -num_iterations 1 -outfmt $F -num_threads 1"`):

```
makeblastdb -in phl_s.faa -dbtype prot -out phl
psiblast -query phl_q.faa -subject phl_s.faa $P > phi_long.sub.6.out
psiblast -query phl_q.faa -db phl $P > phi_long.db.6.out
psiblast -query phl_q2.faa -subject phl_s.faa $P > phi_long.q2.sub.6.out
psiblast -query phl_q2.faa -db phl $P > phi_long.q2.db.6.out
psiblast -query phl_q.faa -subject phl_s.faa $P -evalue 1e-138 > phi_long.e138.sub.6.out
psiblast -query phl_q.faa -subject phl_s.faa $P -evalue 1e-138 -sum_stats true > phi_long.ss.e138.sub.6.out
psiblast -query phl_q2.faa -db phl $P -evalue 1e-138 -sum_stats true > phi_long.ss.e138.q2.db.6.out
```

NCBI's rows at `-num_threads 8` are the same (with `-subject`, NCBI warns that it searches with one thread).

`phi_merge.<case>.out`: NCBI BLAST+ 2.17.0's bytes of the second test of `tests/psiblast_phi_long.rs` (line PS, PS8a-2
fix round 1): psiblast `-phi_pattern` with the zinc-finger pattern `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H` (eight
places of 21 to 24 residues in TF3A_YEAST) on a subject of 5,000,300 residues (`pz_s`, letters without `C` and `H`)
with TF3A at 4,999,771, across the end of the first part. HSPs of the first and second parts that belong to different
places merge, and the merged HSP's Seq-align (its start plus its edit script) ends after its own ends: past the end of
the query `pz_q1` (TF3A, 429 residues: `qend` 442 and 439, `-` in outfmt 0's query row past 429) and inside the query
`pz_q2` (TF3A and 40 residues that the subject does not have; `qend` 258 where the HSP's own end is 248). `qcovhsp`
is the HSP's own span, `qcovs` the Seq-aligns'. The test writes `pz_q1.faa`, `pz_q2.faa`, `pz_s.faa` and `zf.pat`
with a fixed generator; `PHI_MERGE_DIR=<dir>` keeps them. Frozen in that folder with NCBI's own makeblastdb and
psiblast (`F="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore score nident
positive gaps qcovhsp qcovs"`, `P="-phi_pattern zf.pat -num_iterations 1 -num_threads 1"`), for `n` in 1 and 2:

```
makeblastdb -in pz_s.faa -dbtype prot -out pz
psiblast -query pz_qn.faa -subject pz_s.faa $P -outfmt "$F" > phi_merge.qn.sub.6.out
psiblast -query pz_qn.faa -db pz $P -outfmt "$F" > phi_merge.qn.db.6.out
psiblast -query pz_qn.faa -subject pz_s.faa $P -outfmt "$F" -sum_stats true > phi_merge.qn.ss.sub.6.out
psiblast -query pz_qn.faa -db pz $P -outfmt "$F" -sum_stats true > phi_merge.qn.ss.db.6.out
psiblast -query pz_qn.faa -subject pz_s.faa $P -outfmt 0 -evalue 1e-180 > phi_merge.qn.sub.0.out
psiblast -query pz_qn.faa -db pz $P -outfmt 0 -evalue 1e-100 -sum_stats true > phi_merge.qn.ss.db.0.out
```

The third query `pz_q3` (TF3A's first 81 residues, with its place at 50 of 24 residues, then its place at 164 of 22,
which ends the query: 103 residues) is put at 4,999,909 of `pz_r` (the same letters as `pz_s` otherwise). The merged
HSP aligns its pattern region from the last place's start with the first place's length, two letters past the query's
end: NCBI reads the query's sentinel and the zero byte after it there (`qend` 115, `-` in outfmt 0's query row):

```
makeblastdb -in pz_r.faa -dbtype prot -out pr
psiblast -query pz_q3.faa -subject pz_r.faa $P -outfmt "$F" > phi_merge.q3.sub.6.out
psiblast -query pz_q3.faa -db pr $P -outfmt "$F" > phi_merge.q3.db.6.out
psiblast -query pz_q3.faa -subject pz_r.faa $P -outfmt "$F" -sum_stats true > phi_merge.q3.ss.sub.6.out
psiblast -query pz_q3.faa -subject pz_r.faa $P -outfmt 0 > phi_merge.q3.sub.0.out
```

The test compares outfmt 0 in LOSAT's words (`wording::losat_output_bytes`) and without the line of the database's
`Posted date:`.
