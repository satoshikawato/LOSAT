# TBLASTN and TBLASTX `-db_hard_mask` overlapping-mask fixtures

Frozen inputs and NCBI BLAST+ 2.17.0 (`-num_threads 1`) outputs for
`LOSAT/tests/t_db_hard_overlap.rs` (series NP, M4b-3c part 2).

## Tool and command lines

Official NCBI BLAST+ 2.17.0 release binaries (`makeblastdb`, `tblastn`, `tblastx`).

Database build (for each `<db>` in `s2`, `nest`, `g3`, `uns`, `lk`, `fb`):

```sh
makeblastdb -in <db>.fna -dbtype nucl -out <db> -mask_data <db>.asnt
```

Searches (every search uses `-db_hard_mask 11 -seg no -outfmt "6 std sframe" -num_threads 1`,
with query `q.faa` for `tblastn` and `qx.fna` for `tblastx`):

- `tblastn.s2.def.1a.6.out`: `tblastn` on `s2` (default options)
- `tblastn.s2.w0.1a.6.out`: `tblastn` on `s2` (`-window_size 0 -evalue 1000`)
- `tblastn.s2.ung.1a.6.out`: `tblastn` on `s2` (`-ungapped -comp_based_stats 0 -evalue 1000`)
- `tblastx.s2.xdef.1a.6.out`: `tblastx` on `s2` (`-evalue 1000`)
- `tblastx.s2.xw0.1a.6.out`: `tblastx` on `s2` (`-evalue 1000 -window_size 0`)
- `tblastn.nest.def.1a.6.out`: `tblastn` on `nest` (default options)
- `tblastn.nest.w0.1a.6.out`: `tblastn` on `nest` (`-window_size 0 -evalue 1000`)
- `tblastx.nest.xw0.1a.6.out`: `tblastx` on `nest` (`-evalue 1000 -window_size 0`)
- `tblastn.fb.w0.1a.6.out`: `tblastn` on `fb` (`-window_size 0 -evalue 1000`)
- `tblastx.fb.xw0.1a.6.out`: `tblastx` on `fb` (`-evalue 1000 -window_size 0`)

## Layouts

Each nucleotide sequence is 720 bases (240 residues in frame +1), except `fb`'s first sequence `F`
which is 722 bases (240 residues in every frame). Mask intervals below are in frame +1 residues
unless noted as bases:

- `s2`: masks `[40, 60)`, `[55, 70)`, `[72, 90)`; `WWW` at residues 60 and 70 (an inverted range
  at 60 followed by a 2-residue range at 70).
- `nest`: nested masks `[40, 200)`, `[60, 100)`; the query's first 33 residues at 120 (inside the
  outer mask, searched by `[100, 240)`) and at 205.
- `g3`: masks `[40, 60)`, `[55, 80)`, `[75, 90)`; `WWW` at 60 and 80 (two consecutive inverted
  ranges; NCBI `tblastn` two-hit gapped search splits its interval tree without end).
- `uns`: masks out of start order `[100, 120)`, `[40, 60)` (refused before the search).
- `lk`: `A` with masks `[40, 200)`, `[60, 238)` (chunks of 40, -140, and 2 residues per plus
  frame) followed by `B` with mask `[230, 235)` (`B`'s hits reach `A`'s entries on NCBI's table
  of diagonals; refused before the search).
- `fb`: `F` (722 bases) with masks `[3, 717)`, `[6, 717)` in bases (chunks of 1, -237, and 1
  residues per frame, storing no entries at `-window_size 0` while moving NCBI's one-thread offset
  to -1410) followed by `B` (`B` meets only entries that no sequence wrote).
