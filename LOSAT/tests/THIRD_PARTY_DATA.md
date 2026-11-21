# Third-party data in the tests

LOSAT's own work in this directory (the tests, scripts and fixtures that LOSAT wrote) is under the MIT
License (the file `LICENSE`). The sequences, annotation texts and database rows in the test inputs and in the
frozen outputs come from other projects. This file says whose they are, under what terms, and where they are. The
README of the directory that holds an input says how it was made; this file is the one place for the notices.

These data are in the repository and in the source archives that GitHub makes from it. They are not in the
crate, the pylosat wheel or the release archives of the binaries.

The file lists below are not exhaustive: a frozen output that prints a title or a sequence of one of these records
carries the data too. To find them, search the directory for the record ids named here.

## Data that ask for attribution

### UniProtKB/Swiss-Prot (CC BY 4.0)

- What: sequences and definition lines (names, the text of `RecName: Full=`) of UniProtKB/Swiss-Prot records: 23
  globins (`db_search/taxfilter/tfp.faa`, and the databases and inputs made from it), `P69905.2` and other
  Swiss-Prot records in the query and subject files of `db_search/` and `fasta/psiblast/`, and five records of
  UniProtKB/Swiss-Prot release 2026_03 (`O15067`, `P68501`, `P12042`, `B1YBM5`, `P68502`; `fasta/outfmt0/`). The
  globins and `P69905.2` come from the `swissprot` BLAST database that NCBI distributes; the release of that copy is
  not recorded.
- Where: about 300 files: `db_search/` (inputs, `multidb/`, `coldb/`, `taxfilter/` and others), `fasta/outfmt0/`,
  `fasta/psiblast/`, `blastdbcmd/`, and the frozen outputs under `fixtures/` that print them.
- Changes: the records are written into test FASTA files and BLAST databases. Some records carry an added second
  definition line, and the titles of the five 2026_03 records in `fasta/outfmt0/` are shortened or rewritten (for
  example `q_mtb MTB_ONCMY (UniProtKB P68501) metallothionein B, Oncorhynchus mykiss`).
- License: Creative Commons Attribution 4.0 International (CC BY 4.0), <https://creativecommons.org/licenses/by/4.0/>.
  Source: UniProt Consortium, <https://www.uniprot.org/>; terms: <https://www.uniprot.org/help/license>. The statement
  in the LICENSE file of UniProt's FTP directory reads: "We have chosen to apply the Creative Commons Attribution 4.0
  International (CC BY 4.0) License (https://creativecommons.org/licenses/by/4.0/) to all copyrightable parts of our
  databases. All databases and documents in the UniProt FTP directory may be copied and redistributed freely, without
  advance permission, provided that this copyright statement is reproduced with each copy."

### FlyBase, through Ensembl BDGP6.54 (CC BY 4.0)

- What: five *Drosophila melanogaster* cDNA records (`FBtr0078009`, `FBtr0100310`, `FBtr0083245`, `FBtr0301286`,
  `FBtr0078640`) with their gene ids (`FBgn...`), taken from Ensembl's BDGP6.54 `cdna` file. Several have the Ensembl
  definition line cut to the identifier.
- Where: 11 FASTA files (`blastx_diag_bitfield/query.fna`, `fasta/outfmt0/bx_best_hsp_query.fna`,
  `db_search/taxfilter/tfn.fna` and 8 in `db_search/multidb/`), the BLAST databases and id maps made from them
  (`db_search/`, `blastdbcmd/db/`), and the frozen outputs that print the ids; about 230 files in all.
- License: FlyBase states that all annotations and data produced by FlyBase are distributed under CC BY 4.0
  (<https://wiki.flybase.org/wiki/FlyBase:About>, "FlyBase Licenses and Availability"), and that different terms may
  apply to data from other resources that FlyBase does not distribute. Source: FlyBase, <https://flybase.org/>.

### SGD, through Ensembl R64-1-1 (CC BY 4.0)

- What: three *Saccharomyces cerevisiae* protein records, `YLR106C` (MDN1), `YBR142W` (MAK5) and `YOL075C`, from
  Ensembl's R64-1-1 `pep` file, with the title cut to the identifier in two of the files.
- Where: `blastx_diag_bitfield/subjects.faa`, `fasta/outfmt0/bx_best_hsp_subject.faa`, and about 30 files that print
  the ids (`blastx_diag_bitfield/`, `fixtures/outfmt0/`, `fixtures/db_search/`, `db_search/bxbest.*`).
- License: Saccharomyces Genome Database (SGD), <https://www.yeastgenome.org/>, Creative Commons Attribution 4.0
  International, <https://creativecommons.org/licenses/by/4.0/>. This is as SGD's help pages state it; SGD's own terms
  page could not be read directly when this notice was written (2026-11-21).

### TIGRFAMs, through NCBI's Conserved Domain Database (CDD) v3.22 (CC BY-SA 4.0)

- What: the description text of 29 TIGRFAMs families (rows `TIGR...` of CDD's `cddid.tbl`: accession, short name and
  description). The accessions and short names of other TIGRFAMs families also appear in the outputs of `rpsbproc`
  and the other rpsblast-family tests.
- Where: the 10 `cddid.tbl` files of `fixtures/rpsbproc/` that hold a TIGR row with a description (`data/` and
  `data_with/{badmotif,ctprev,d_acc0,d_wsline,nt,ntprev,sites,xnew,xnewprev}/`). They are excerpts of the file
  `cddid.tbl` of CDD v3.22 (NCBI's FTP site, taken 2026-11-11), cut to the rows that the test archives reach. The
  TIGR rows are otherwise unchanged, except that in `data_with/{nt,ntprev,xnew,xnewprev}/` the description is cut to
  24 characters.
- License: Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0),
  <https://creativecommons.org/licenses/by-sa/4.0/>, full text
  <https://creativecommons.org/licenses/by-sa/4.0/legalcode>. NCBI (National Library of Medicine, NIH) holds the
  rights since April 2018 and publishes the data under this license
  (<https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/license_and_availability.txt>). The TIGRFAMs rows in the files named
  above stay under CC BY-SA 4.0; the rest of this directory is not part of them.
- NCBI asks that use of the TIGRFAMs database cite the four articles that describe it:
  - Haft DH, Loftus BJ, Richardson DL, Yang F, Eisen JA, Paulsen IT, White O. TIGRFAMs: a protein family resource for
    the functional identification of proteins. Nucleic Acids Res. 2001 Jan 1;29(1):41-3. PMID: 11125044
  - Haft DH, Selengut JD, White O. The TIGRFAMs database of protein families. Nucleic Acids Res. 2003 Jan
    1;31(1):371-3. PMID: 12520025
  - Selengut JD, Haft DH, Davidsen T, Ganapathy A, Gwinn-Giglio M, Nelson WC, Richter AR, White O. TIGRFAMs and Genome
    Properties: tools for the assignment of molecular function and biological process in prokaryotic genomes. Nucleic
    Acids Res. 2007 Jan;35(Database issue):D260-4. Epub 2006 Dec 6. PMID: 17151080
  - Haft DH, Selengut JD, Richter RA, Harkins D, Basu MK, Beck E. TIGRFAMs and Genome Properties in 2013. Nucleic
    Acids Res. 2013 Jan;41(Database issue):D387-95. doi: 10.1093/nar/gks1234. Epub 2012 Nov 28. PMID: 23197656

## Data with terms that are partly stated

### SMART, through CDD v3.22

- What: the description text of 12 SMART families (rows `smart...` of `cddid.tbl`), and the `smart` rows of CDD's
  `bitscore_specific.txt` (bit-score thresholds, as NCBI publishes them in CDD). No SMART hidden Markov models or
  alignments are in the tests.
- Where: the 8 `cddid.tbl` files of `fixtures/rpsbproc/` that hold a `smart` row with a description (`data/` and
  `data_with/{ctprev,d_acc0,d_wsline,nt,ntprev,xnew,xnewprev}/`) and the `bitscore_specific.txt` files next to them.
  Same excerpt rule as for TIGRFAMs above.
- Terms: SMART (<https://smart.embl.de/about.cgi>) states that its data (hidden Markov models, alignments and
  thresholds) are available via a license that is free to academics and not to commercial users. It does not state
  terms for the description text of a family. The source is SMART, <https://smart.embl.de/>, as distributed in CDD.

## Data that need no attribution

These are listed so that the list of the third-party data is complete.

- Pfam, through CDD: 75 `pfam` rows of `cddid.tbl`. CC0 1.0 (<https://interpro-documentation.readthedocs.io/en/latest/license.html>).
- NCBI's own CDD rows (`cd`, `cl`, `sd`, `COG`, `PRK`, `PLN`, `CHL`, `PTZ`, `PHA` and `NF` accessions), the tables and
  the rpsblast archives of CDD v3.22 under `fixtures/rpsbproc/`; NCBI RefSeq sequences, NCBI Taxonomy rows
  and the NCBI unit-test data of `seqdb/db/`: NCBI places no restriction on their use or distribution
  (<https://www.ncbi.nlm.nih.gov/home/about/policies/>), and the policies note that submitters or third parties may hold rights in
  some records.
- INSDC (DDBJ, ENA, GenBank) records: the INSDC policy places no restriction on use
  (<https://www.insdc.org/policy/>).
- Ensembl: the *D. melanogaster* BDGP6.54 genome pieces (`winmask/`), the GRCh38 and GRCz11 cDNAs (`taxfilter/tfn.fna`
  and the files made from it). Ensembl's terms place no restriction on the data
  (<https://www.ensembl.org/info/about/legal/disclaimer.html>); the READMEs of these directories name the source.
- Protein Data Bank: only made-up or reused ids (`pdb|1ABC|A` and similar) are in the tests. The PDB data are CC0 1.0
  (<https://www.wwpdb.org/about/usage-policies>).
