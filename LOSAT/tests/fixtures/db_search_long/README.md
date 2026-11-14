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
