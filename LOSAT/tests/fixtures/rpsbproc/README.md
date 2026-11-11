# rpsbproc fixtures (PS16a)

Inputs and frozen NCBI outputs for `losat rpsbproc` (NCBI rpsbproc 0.5.1, the 2025-05-13 package).

## Provenance

- Archives (`arch/*.asn`): NCBI BLAST+ 2.17.0 `rpsblast -outfmt 11`, one thread, against CDD v3.22
  (`Cdd_LE.tar.gz` of 2026-10-27, run in its directory as `-db Cdd`), queries from E. coli K-12 MG1655
  (`GCF_000005845.2_ASM584v2_protein.faa`). The commands are in `ARCH_COMMANDS.tsv`. `lcl.asn` searches a profile
  database made by `makeprofiledb` whose ids are `lcl|` (LOSAT's rpsblast fixture `prof_rps`). `bin.asn` is the binary
  archive (`ARCHIVE_FORMAT=bin`), `one.xml` is `-outfmt 5`, `garbage.txt` is not an archive.
- `arch/bad/*.asn`: `one.asn` edited by hand (truncated, an unknown member, a duplicate member, a trailing comma, no
  `results`, no `alignments`, comments, CRLF line ends, two objects without a line end between them, and so on).
  Keep their bytes as they are (CRLF included).
- `data/`: the CDD v3.22 data files of NCBI's FTP site (taken 2026-11-11; `family_superfamily_links` as rewritten on
  2026-11-04), cut to the rows that the archives reach. `cdtrack.txt` is empty, as NCBI publishes it for v3.22.
  For every case, NCBI's output with this cut data equals its output with the full files.
- `oracle/<case>.{out,err,rc}` (and `.outfile` for `-o`): NCBI's rpsbproc 0.5.1 Linux x64 binary of the 2025-05-13
  package, run in a directory that holds `arch/` and `data/`, with the arguments of `CASES.tsv` and stdin from
  `/dev/null` unless the case names a file. Bioconda's `rpsbproc 0.5.1` (built on the BLAST+ 2.17.0 toolkit) gave the
  same bytes for every case.

## CASES.tsv

`case`, `use` (`compare`: LOSAT must give the same stdout, stderr and exit code; `refuse`: LOSAT refuses until a later
item: `-t feats` and `-t both` need sites and motifs, XML input needs the XML reader), `stdin` (`-` = empty), `setup`
(`data-without:<file>`: a directory `data_missing/` holding the data files except that one), `arguments` (split on
spaces).
