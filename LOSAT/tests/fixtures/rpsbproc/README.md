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
- The cases of the PS16a review (2026-11-11; the same two NCBI binaries, same bytes), run as above:
  - `root/*.asn`, which the tests copy to the run directory itself (the report's header names the `-i` file as given):
    `bracket_member.asn`, `neg0_length.asn`, `org_descr.asn` and `nolength.asn` are `arch/one.asn`, an edited copy of
    it and `arch/one.asn` again (a member id in brackets, `length -0`, an `org` descriptor, no `inst.length`; oracle
    requests O3-2, O3-1, O3-6, O3-9 of the review's angle 3). `q_tab_utf8.asn` is NCBI rpsblast 2.17.0's own archive
    (`-db Cdd -evalue 0.01 -outfmt 11`, CDD v3.22) of two copies of `NP_414783.1` whose definition lines hold a TAB and
    a UTF-8 letter (O3-10).
  - `data_with/<dir>/`: the data files of `data/` that the case `d_<name>` edits; the tests make `<dir>/` from `data/`
    and these files (oracle request A4-1 of the review's angle 4): `d_acc0` adds a `cddid.tbl` row with an empty
    accession, `d_wsline` a line of one space, `d_annot10` a `cddannot.dat` line whose coordinates are its last
    field, `d_bs12` writes the threshold of `smart00421` as `0x10`, `d_bsnul` writes the bit score of `CHL00113` as
    `1x`, a 0 byte and `y` (NCBI quotes the string up to the 0 byte, exit 2), and `d_ct3`'s `cdtrack.txt` is the line
    `smart00421` alone.
  - `xml_missing_input`: `-i` names a missing file and the standard input is XML (oracle request M4 of the review).

## CASES.tsv

`case`, `use` (`compare`: LOSAT must give the same stdout, stderr and exit code; `refuse`: LOSAT refuses until a later
item: `-t feats` and `-t both` need sites and motifs, XML input needs the XML reader), `stdin` (`-` = empty), `setup`
(`data-without:<file>`: a directory `data_missing/` holding the data files except that one; `data-with:<dir>`: a
directory `<dir>/` holding the data files with those of `data_with/<dir>/` in their place), `arguments` (split on
spaces).
