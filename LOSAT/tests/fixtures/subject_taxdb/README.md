# subject_taxdb

NCBI BLAST+ 2.17.0 (`env_blast` on ws, the official release's bytes on the cases checked) for the taxonomy
database of a `-subject` search's tabular reports (line RG item E1w-2b, ad-RG-E-15, 2026-11-14;
`remaining-gaps-20261112/E/e1w2b/probe/probe_g3.sh`). Inputs `tests/fasta/outfmt11/arch_*` (blastn and tblastx:
`arch_tbx_q.fna` / `arch_tbx_s.fna`), run with `env -i PATH HOME=/nonexistent BLAST_USAGE_REPORT=false`.

`<program>.<format>.<environment>.{out,err}`: the standard output and error.
- format: `sing` = `6 qseqid sseqid ssciname scomname sblastname sskingdom`, `plur` = `6 qseqid sseqid sscinames`,
  `csv20` = `20 qseqid sscinames`.
- environment: `none` = no `BLASTDB`, `taxdb` = `BLASTDB=<tests/seqdb/taxdb/good>`,
  `badenv` = `NCBI_CONFIG__BLAST__BLASTDB` set to the bytes `x\xffy`.
All runs exit 0. Used by `tests/library_out_taxdb.rs`.
