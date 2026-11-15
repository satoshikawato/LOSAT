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
