- Web ABI v1: BLASTN, BLASTP and TBLASTX read their FASTA inputs with NCBI's reader, as the CLI and ABI v2 do
  (plan TD-1 revised, DW-48). The inputs that v1 refused because the `bio` crate read them differently from
  NCBI now behave as on the command line: non-ASCII and control bytes in deflines, residues that the reader
  removes with a warning, and records without residues. The arguments, formats and request errors stay frozen.
- LOSAT no longer depends on the `bio` crate: 75 packages left `Cargo.lock`, among them the three unmaintained
  crates that the nightly `cargo audit` ignored, so that job has no `--ignore` list any more.
