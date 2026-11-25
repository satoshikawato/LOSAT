# `tblastn_db_soft_lead_start/`: frozen rows for pre-chunk short alignments and one-thread diagonal-table clears

`<name>.6.out`: NCBI BLAST+ 2.17.0 `-num_threads 1` outputs for `tests/tblastn_db_soft_lead_start.rs`
(series NP, M4b-T3): `tblastn` and `tblastx` with `-db_soft_mask` on nucleotide database sequences
longer than 15,000,000 bases whose masks are stored out of order, so that a later range of the
second chunk of frame +1 (starting at residue 4,999,900) scans letters before the chunk.

## Layout

- Query protein `q.faa` (118 residues) and its coding DNA `qx.fna` (354 bases, one codon per
  residue): seven unique cores (`a1`, `a2`, `a3`, `e1`, `b1`, `c1`, `d1`) separated by `GWG`
  flanks.
- `q1a1`, `q1a2`, `q1c1`, `q1d1`: one 15,300,000-base poly-A sequence with out-of-order dust
  masks `[(14_999_400, 14_999_849), (14_999_880, 15_179_999), (3_000_000, 3_299_999)]` (stored
  order), placing the case core at chunk-2 offset `-4` (`q1a1`, 11 aa, gapped start `+1`), `-5`
  (`q1a2`, 11 aa, gapped start `0`), `-6` (`q1c1`, 11 aa, gapped start `-1`), or `-4` (`q1d1`,
  12 aa, window read), plus a control copy at residue 5,019,900.
- `C0`, `C2`: poly-A filler sequences (`t0` 10 bases, `f1` 40,000,000 bases, `f2`..`f16` 15 x
  15,000,000 bases, `flast` 3,434,890 bases in `C0` or 3,434,889 bases in `C2`), followed by `L`
  (OID 18, 15,300,000 bases with the same out-of-order masks, `a1` at `-4` and control `+20_000`,
  `e1` at `-1` and control `+30_000`, and `b1` at `-30`). In NCBI's one-thread `-window_size 0`
  search, the diagonal offset before `L` is `536_870_912` (`INT4_MAX / 4 + 1`) in `C0` (clearing
  the diagonal table after chunk 1 of frame +1 so pre-chunk words in chunk 2 are skipped) and
  `536_870_910` in `C2` (no clear after chunk 1).

The test writes the FASTA and `.asnt` files itself; setting `TBLASTN_DB_SOFT_LEAD_START_DIR=<dir>`
keeps them in `<dir>`.

## Frozen commands (NCBI BLAST+ 2.17.0, `BLAST_USAGE_REPORT=false`)

```
RICH="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qframe sframe evalue bitscore score nident positive btop"

for db in q1a1 q1a2 q1c1 q1d1 C0 C2; do
  makeblastdb -in ${db}.fna -dbtype nucl -out ${db} -mask_data ${db}.asnt
done

tblastn -query q.faa -db q1a1 -db_soft_mask 11 -seg no -window_size 0 -outfmt 6 -num_threads 1 > q1a1.tblastn_w0.6.out
tblastn -query q.faa -db q1a1 -db_soft_mask 11 -seg no -window_size 0 -outfmt "$RICH" -num_threads 1 > q1a1.tblastn_w0rich.6.out
tblastn -query q.faa -db q1a2 -db_soft_mask 11 -seg no -window_size 0 -outfmt 6 -num_threads 1 > q1a2.tblastn_w0.6.out

tblastn -query q.faa -db C0 -db_soft_mask 11 -seg no -ungapped -comp_based_stats 0 -window_size 0 -outfmt 6 -num_threads 1 > C0.tblastn_ungw0.6.out
tblastn -query q.faa -db C0 -db_soft_mask 11 -seg no -window_size 0 -comp_based_stats 0 -outfmt "$RICH" -num_threads 1 > C0.tblastn_w0cbs0rich.6.out
tblastn -query q.faa -db C0 -db_soft_mask 11 -seg no -window_size 0 -outfmt 6 -num_threads 1 > C0.tblastn_w0.6.out
tblastx -query qx.fna -db C0 -db_soft_mask 11 -seg no -window_size 0 -outfmt 6 -num_threads 1 > C0.tblastx_w0.6.out

tblastn -query q.faa -db C2 -db_soft_mask 11 -seg no -ungapped -comp_based_stats 0 -window_size 0 -outfmt 6 -num_threads 1 > C2.tblastn_ungw0.6.out
tblastx -query qx.fna -db C2 -db_soft_mask 11 -seg no -window_size 0 -outfmt 6 -num_threads 1 > C2.tblastx_w0.6.out
```

`tblastn -query q.faa -db q1c1|q1d1|C2 -db_soft_mask 11 -seg no -window_size 0 -num_threads 1`
each end in a segmentation fault (exit 139) in NCBI BLAST+ 2.17.0; the test checks LOSAT's
refusal (exit 1, empty stdout, naming the database sequence, frame, and chunk start) at 1 and 4
threads.
