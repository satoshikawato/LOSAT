# `tblastx_db_soft_lead/`: frozen rows of searches whose chunks have ranges that start before them

`<name>.6.out`: NCBI BLAST+ 2.17.0's rows of `tests/tblastx_db_soft_lead.rs` (series NP, M4b-T2): TBLASTX and TBLASTN
`-db_soft_mask` on 15.3 Mbp database sequences whose masks are stored out of order, so that a later range of the second
chunk of frame +1 starts before the chunk (in `lead2`, after a first range raised to one word on the frame's last
letter: the one-hit word finder of `-window_size 0` scans none of that chunk, as its first range is shorter than a
word, and the two-hit one goes on to the later ranges). The test writes its inputs (`q.faa`, `qx.fna`, `lead.fna`,
`lead_masks.asnt`, `before.fna`, `before_masks.asnt`); `TBLASTX_DB_SOFT_LEAD_DIR=<dir>` keeps them in `<dir>`. Frozen
from those inputs with NCBI's own makeblastdb and programs (`BLAST_USAGE_REPORT=false`; the rows at `-num_threads 2`
are the same):

```
makeblastdb -in lead.fna -dbtype nucl -out lead -mask_data lead_masks.asnt
makeblastdb -in before.fna -dbtype nucl -out before -mask_data before_masks.asnt
tblastx -query qx.fna -db lead -db_soft_mask 11 -seg no -window_size 0 -outfmt 6 -num_threads 1 > lead.tblastx_w0.6.out
tblastn -query q.faa -db lead -db_soft_mask 11 -seg no -window_size 0 -outfmt 6 -num_threads 1 > lead.tblastn_w0.6.out
tblastn -query q.faa -db lead -db_soft_mask 11 -seg no -outfmt 6 -num_threads 1 > lead.tblastn.6.out
tblastx -query qx.fna -db before -db_soft_mask 11 -seg no -window_size 0 -outfmt 6 -num_threads 1 > before.tblastx_w0.6.out
tblastx -query qx.fna -db before -db_soft_mask 11 -seg no -outfmt 6 -num_threads 1 > before.tblastx.6.out
tblastn -query q.faa -db before -db_soft_mask 11 -seg no -outfmt 6 -num_threads 1 > before.tblastn.6.out
tblastn -query q.faa -db before -db_soft_mask 11 -seg no -ungapped -comp_based_stats 0 -window_size 0 -outfmt 6 -num_threads 1 > before.tblastn_ungapped_w0.6.out
```

`tblastn -query q.faa -db before -db_soft_mask 11 -seg no -window_size 0` ends in a segmentation fault (exit 139, at 1
and 2 threads), with no rows; the test checks LOSAT's refusal there.
