# rpsbproc fixtures (PS16a-PS16d, PS16-Q1)

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

`case`, `use` (`compare`: LOSAT must give the same stdout, stderr and exit code; `linux`: the same, on Linux only; `refuse`: LOSAT refuses, exit 1, with
the bytes of `refused/<case>.{out,err,rc}`: the refusals of PS16b, PS16c and PS16d below), `stdin` (`-` = empty), `setup`
(`data-without:<file>`: a directory `data_missing/` holding the data files except that one; `data-with:<dir>`: a
directory `<dir>/` holding the data files with those of `data_with/<dir>/` in their place; `data-with:<a>+<dir>`: the
same with those of `data_with/<a>/` first, then those of `data_with/<dir>/`), `arguments` (split on
spaces). PS16a's cases of `-t feats`, `-t both` and the default `-t` (named `refuse_*` until PS16b) are `t_feats`,
`t_both` and `t_default`.

## Sites and motifs (PS16b)

- Archives added for PS16b (the same NCBI rpsblast 2.17.0 run, `-evalue 0.01`, commands in `ARCH_COMMANDS.tsv`):
  - `mstate.asn` (`BATCH_SIZE=1`: three sessions): NP_414781.1 and NP_416506.1 (cd08071, "Zn binding site", motif
    `HHD`), then NP_416542.1 (cd00438, "catalytic site", motif `HKY`). NCBI's `MotifCheck` keeps the positions of the
    first motif site it checks in the process (a function-static list) and checks every later motif site against them:
    the second query's `HHD` site prints, the third query's `HKY` site does not, across sessions.
  - `mfirst.asn`: NP_416542.1 alone; its `HKY` site prints.
  - `gen2.asn`: NP_416215.1 and NP_415978.1, whose hits have generic sites of hierarchy roots (with `ctprev`).
  - `motdom.asn`: NP_416435.1, a `MOTIFS` line of a curated domain (a source PSSM-ID other than 0).
  - `sdsite.asn`: NP_416698.1, sites and motifs of a structural domain (`sd`, source 0).
  - `sd2.asn` (PS16a): its second query NP_418053.1 has structural-domain hits and prints no `SITES` block.
- `data/` stays PS16a's. The PS16b archives read data sets of their own under `data_with/`, cut by the same rule
  (the rows the named archives reach), so that no file passes 64 KB:
  - `sites/`: the PS16b archives; `cdtrack.txt` empty (from `data/`), as NCBI publishes it for v3.22.
  - `ctprev/`: the PS16b archives and `three.asn`, with `cdtrack.txt` = their rows of the CDD v3.23 preview
    `cdtrack.txt` (`pub/mmdb/cdd/preview/cdtrack.txt`, 2026-11-04) and the `cddid.tbl` rows of the parents and roots
    those rows name: the `Generic` sites.
  - `rootpid/`: `gen2.asn` as `ctprev/`, but the parents and roots that `gen2.asn` does not hit (cd04410) are
    renumbered to PSSM-ID 100,000,000 + their own: a Generic site whose root has a PSSM-ID of 100,000,000 or more
    prints the root's accession instead of the number.
  - `badmotif/`: `mfirst.asn` as `sites/`, with the motif of cd00438 "catalytic site" written `K#R`, which NCBI's
    PROSITE parser refuses (exit 3, the output ends inside `SITES`).
- For every PS16b case, NCBI's output with the cut data equals its output with the full files (`ctprev/`: the full
  v3.22 files with the preview `cdtrack.txt`), and Bioconda's rpsbproc 0.5.1 gave the same bytes as NCBI's binary.
- The cases of the PS16b review (2026-11-12; oracle requests R1, R2, R5, R6, R20, R21; NCBI's binary and Bioconda's
  gave the same bytes), run as above. `data_with/rd1/`, `rd5/`, `rd4/`, `rd13/` and `rd14/` hold only `cddannot.dat`,
  that of `sites/` with the motif of cd00438 "catalytic site" (`HKY`) replaced; the cases lay them over `sites/`
  (`data-with:sites+rd1`):
  - `rd1` `Hx(0,10)K` and `rd5` `H{P}(0,10)K`: NCBI's match reads the residues `HKY` of the site beyond their end
    (the query has 185), but every byte there gives the same answer: the site prints.
  - `rd4` `K` 0xE9 `R`: a parse error whose message quotes the byte as it is (exit 3).
  - `rd13` `Hx(99999999999999999999)K` (`atol` gives `LONG_MAX`) and `rd14` `Hx(200000000000000)K`: the minimal map
    passes the string's `max_size()` (`basic_string::_M_replace_aux`) or needs 2e14 bytes (`std::bad_alloc`); exit 3,
    the report cut as for a parse error (NCBI run under `prlimit --as=4294967296`).
  - `arch/mfirst_delta.asn`: `mfirst.asn` with its query as `repr delta`, the residues in two `ext delta` literals
    (NCBIeaa, 100 and 85 residues); case `both.mfirst_delta` with `sites/`.
- The cases of the second PS16b review (2026-11-12, oracle O2 on ws: NCBI's binary and Bioconda's 0.5.1, the same
  bytes), as above: `rd15` `Hx(4611686018427387903)K` and `rd18` `Hx(9223372036854775806)K`, minimal maps of 2^62 and
  2^63 - 1 letters, which NCBI's C++20 builds allocate (`max_size()` is 2^63 - 1) and fail with `std::bad_alloc`
  (exit 3, the report cut as for a parse error).
- `refused/<case>.{out,err,rc}`: LOSAT's own bytes (written from LOSAT's output and checked by hand: standard output is
  NCBI's up to the line before the refused site or query, standard error ends with the refusal). The refusals of
  PS16b (PD-LOSAT-NCBI-DEFECTS 1.15, explicit rejections):
  - `refuse_rd2.mfirst`: `data_with/rd2/` `Hx(0,10)W` (review oracle R4): NCBI's match reads beyond the residues
    `HKY` and its answer depends on the bytes there (NCBI's stack: the review's runs dropped the site, but no line of
    NCBI's code decides it).
  - `refuse_rd19.mfirst`: `data_with/rd19/` `H-{<}(0,67108863)-{<}(0,67108863)-K` (review 2 oracle O3): NCBI holds
    1.5 GiB of match records at once (exit 0 under a 4 GiB address space, `std::bad_alloc` under 1.5 GiB).
  - `refuse_rd20.mfirst`: `data_with/rd20/` `Hx(536870910)K` (O3): the minimal map's string of 512 MB grows to
    1 GiB while the old buffer lives (as rd19 on NCBI).
  - `refuse_nodata`: `arch/q_nodata.asn`, `mfirst.asn` without the query's `seq-data` (review oracle R19), with
    `sites/`: the site's positions are outside the (empty) residues, where NCBI reads memory outside the string
    (its two builds print different bytes).

## Nucleotide queries (PS16c)

- Archives added for PS16c: NCBI BLAST+ 2.17.0 `rpstblastn -outfmt 11`, one thread, against CDD v3.22 (run in its
  directory as `-db Cdd`), queries cut from the E. coli K-12 MG1655 genome (`NC_000913.3`) as 5 kb pieces named
  `mgNNNNN_<from>-<to>` (1-based genome positions); commands in `ARCH_COMMANDS.tsv`. They are stored gzip-compressed
  (`arch/*.asn.gz`, `gzip -9n`); the tests unpack each to `arch/<name>.asn` in the run directory before running, so the
  report's header names `arch/<name>.asn`.
  - `nt67.asn` (mg00067), `nt63.asn` (mg00063), `nt49.asn` (mg00049): hits in plus and minus frames, specific sites
    with codon ranges (minus strand descending), superfamilies, and Generic sites of hierarchy roots with `ntprev`.
  - `nt68.asn` (mg00068): `MOTIFS` lines of a structural domain (ANK repeats, source 0) on the plus strand;
    `nt68rc.asn`: the reverse complement of the same piece (title `..._rc`), the same motifs on the minus strand.
  - `nt60.asn`, `nt60t1.asn`, `nt60t2.asn`: mg00060 and the same piece 1 and 2 bases shorter at its end; its hits are
    on the minus strand, whose frame numbers follow the query length modulo 3.
  - `nt04amb.asn`: mg00004 with `N` at every 97th base from base 51, `R` at every 211th from 81 and `Y` at every 307th
    from 131 (`ncbi4na` Seq-data); a site codon with an ambiguous base prints `X`.
  - `nt67seg.asn` (`-seg yes`: `results.masks` with frames), `nt67e10.asn` (the default E-value).
  - `ntpair.asn`: mg00063 and mg00067 in one batch (one session, two queries).
  - `ntnohit.asn`: mg00001 to mg00010 with `-strand minus`, where NCBI's rpstblastn finds no hit: every query prints
    `DOMAINS` and `ENDDOMAINS` with nothing between.
  - `nt94had.asn`: bases 2,301-3,360 of mg00094 (genome 467,301-468,360): a HAD-domain hit whose "HAD signature
    motif I" site carries a motif.
  - `ntmix_p_n.asn` and `ntmix_n_p.asn`: PS16b's `mfirst.asn` (a protein query whose `HKY` catalytic site carries a
    motif) and `nt94had.asn` one after the other, in both orders. NCBI's motif state (`st_motif_pos`) is one per
    process, so it crosses from the protein session to the nucleotide one and back: the second session's motif site
    is dropped.
- `data_with/nt/`: the six data files cut to the rows these archives and `mfirst.asn` reach (the rule of `data/`);
  the description column of `cddid.tbl` (never printed) is cut to 24 characters to keep the file small.
  `data_with/ntprev/` (used as `data-with:nt+ntprev`): the rows of the v3.23 preview `cdtrack.txt` for those domains
  and `nt`'s `cddid.tbl` with the parents and roots those rows name. For every PS16c case, NCBI's output with these
  data equals its output with the full v3.22 files (with the preview `cdtrack.txt` for `ntprev`), and bioconda's
  rpsbproc 0.5.1 gives the same bytes.
- The cases of the PS16c review (2026-11-13; NCBI's binary and bioconda's 0.5.1 gave the same bytes, with these data
  and with the full v3.22 files), archives edited by hand as the review's oracles were (`gzip -9n`):
  - `ntnull1.asn`: `nt67.asn` with a second `Std-seg` in its first alignment whose only location is `null`
    (`loc { null NULL }`). NCBI's `seqlocs[0]->IsNull() || seqlocs[1]->IsNull()` skips it without reading a second
    location: the output is that of `nt67.asn`.
  - `ntshort2.asn`: `nt67.asn` with the Seq-data `iupacna "AC"` (its `length` stays 5000). With `-t doms` NCBI prints
    the hits (`ntshort2.doms`); where it translates the query for sites, see `refuse_ntshort2` below.
- `refused/` holds the refusals of PS16c (PD-LOSAT-NCBI-DEFECTS 1.16) too, each LOSAT's output cut where NCBI's goes on
  (checked against NCBI's output of the same case):
  - `refuse_ntshort2`: `ntshort2.asn` with `-t both`: NCBI's `TranslateAll` erases from an empty string after
    `ENDDOMAINS` (segmentation fault in both builds).
  - `refuse_natoprot`: `arch/natoprot.asn`, `three.asn` with `mol na` (the protein queries read as nucleotide ones:
    their Dense-seg hits in frame 1, their sites on the translation of their letters), `-t both -m std`: a site at
    residue 117 of a translation of 58 residues, where NCBI reads outside the string (NCBI's binary, bioconda's and
    NCBI's binary with the full data files print three different lines).
  - `refuse_ntnoids1`: `nt67.asn` whose first alignment has no ids in its segments (its query id comes from the
    first location): NCBI reads the hit id past the ids (segmentation fault after the `SESSION` line).
  (`refuse_ntplusmax` of PS16c is the `widened` row `ntplusmax.doms.std` since PS16-Q1, below.)

## The NCBI C++ Toolkit options (line RG, item E1c)

- The rows `e1c_*` (2026-11-13): `-dryrun` (rpsbproc's own `DryRun` after its `Init`), `-logfile log.txt` and
  `-conffile c.ini`. `oracle/<case>.{out,err,rc}` and `.outfile` / `.log` (the `out.txt` and `log.txt` the run
  made): NCBI's rpsbproc 0.5.1, the binary of the 2025-05-13 package and bioconda's (on the BLAST+ 2.17.0 toolkit),
  which gave the same bytes, run as the test runs LOSAT (the run directory above, no `LOSAT_` or `NCBI_` variable,
  `NCBI_CONFIG_PATH` = the run directory, stdin from `/dev/null` or the named file).
- `conf/`: the registry files of the setups `conf:<file>` (copied as `c.ini`), `ini:<file>` (as `rpsbproc.ini`) and
  `ncbirc:<file>` (as `.ncbirc`): `blast.ini` is `[BLAST]` alone; `garbage.ini` a line that is not an entry,
  `sec.ini` an unclosed section, `quote.ini` a value with a quote inside, `cont.ini` a continued value then a line
  that is not an entry (NCBI's report, exit 2); `bom2.ini` a UTF-8 byte order mark and one byte more (NCBI's
  `Error: (110.4) Error reading the registry after line 1: `, the run goes on); `bom1.ini` the first byte of one (as a
  `.ncbirc`: the same post; as a `.ini` NCBI posts nothing).
- `arch/bad/subj_sequences.asn`, `matrix_integer.asn`, `evalue_rawscore.asn` (the review of item E1c): `arch/one.asn`
  edited by hand: `subject sequences { }`, `MatrixName` `value integer 5`, `EvalueThreshold` `value cutoff raw-score
  5`. NCBI writes the `SESSION` line's tab before each getter throws; with `-o log.txt -logfile log.txt` its stream
  keeps the part of that line until the program ends, and writes it over the start of the exception's report.
  `matrix_integer_longdb.asn` is `matrix_integer.asn` with a database name of 5,000 bytes, which libstdc++'s buffer of
  8,192 bytes keeps (an insertion is written directly only when it does not fit in what is left);
  `matrix_integer_db8165.asn` one of 8,165 bytes, after which the buffer holds 8,190 bytes and takes the `char` of the
  next tab (`put`) without writing (the review 2 of item E1c).
- `use` `linux` (`linux_cout_longdb`): run only on Linux (`rpsbproc_gives_ncbis_bytes_through_a_pipe_on_linux`):
  `-logfile /dev/stdout` with the standard output a pipe, where glibc's buffer of `cout` (4,096 bytes) decides the
  order of the report and the diagnostics.
- With `-dryrun` NCBI reads no input and passes the XML input, and LOSAT gives its bytes (`e1c_dry_xml`).

## XML input (PS16d)

- XML added for PS16d: NCBI BLAST+ 2.17.0 `rpsblast -outfmt 5` and `rpstblastn -outfmt 5`, one thread, against CDD
  v3.22 (run in its directory as `-db Cdd`), on the queries of the PS16a-c archives (same options, so the same hits)
  and four new nucleotide queries; commands in `ARCH_COMMANDS.tsv`. They are stored gzip-compressed
  (`arch/x_*.xml.gz`, `gzip -9n`); the tests unpack each to `arch/x_<name>.xml` in the run directory before running,
  so the report's header names `arch/x_<name>.xml`. NCBI's rpsbproc posts `Warning: XML format support is deprecated!`
  and reads each `BlastOutput` document as one session (`SESSION 1 RPSBLAST 2.17.0+ Cdd BLOSUM62 0.01`: the version,
  the database, the matrix and `-evalue` of the document); unlike the archive path, `-e` drops the HSPs whose E-value
  is above it, and a query without a kept HSP prints as `Peptide` with no `DOMAINS` block, even a nucleotide one.
  - Protein (`data/`, `data-with:sites`): `x_three`, `x_one`, `x_tiny3` (`BATCH_SIZE=1`: still one document), `x_sd2`,
    `x_src2`, `x_qgap`, `x_dctr_e10` (the default E-value: the `-e` cases), `x_tie1`, `x_tie2`, `x_notitle`, `x_parsed`
    (`-parse_deflines`), `x_mfirst`, `x_gen2`, `x_motdom`, `x_sdsite`, `x_mstate` (the PS16b sites and motifs, the motif
    state across the queries of one document).
  - Nucleotide (`data-with:nt`): `x_ntr67`, `x_ntr67e10`, `x_ntr67seg`, `x_ntr68rc`, `x_ntr94had`, `x_ntr60`,
    `x_ntr60t1`, `x_ntr63`, `x_ntr04amb`, `x_ntpair`, `x_ntnohit` (the PS16c queries; NCBI builds a frame's residues
    from the HSP's query string, `-` where no HSP covers it, instead of translating the query).
  - New queries (`data_with/xnew/`, the rows their hits reach, the rule of `data_with/nt/`; `data_with/xnewprev/`
    the v3.23 preview `cdtrack.txt` rows, as `data-with:xnew+xnewprev`): `x_lcasem` (`three` with residues 21-60 lower
    case and `-lcase_masking`), `x_special` (deflines with `&`, `<`, `>`, quotes, `]]>`, `-->`; the second one holds
    `&amp;` and `&lt;x&gt;` as text, which NCBI's writer leaves unescaped and its reader decodes), `x_mix3` (no-hit
    queries around `one` and `three`), `x_nohit`, `x_lcl` (a profile DB with `lcl|` ids: `Hit_accession` `dom_rr1`
    gives PSSM-ID 0), `x_ntspec1` (mg00017 alone: its minus-strand PRK00915 hit is Specific here and Non-specific in
    the archive path), `x_ntspec2` (mg00006 with bases 1,001-1,600 lower case and `-lcase_masking`: the same for
    COG1047), `x_ntmix2` (a 600 nt polyA, mg00067 and a 6 nt query), `x_ntp` (`-parse_deflines`, an id without a
    title: `No definition line`).
  - Streams made by hand: `x_cat2` (`x_one` then `x_three`: two sessions), `x_cattrunc` (`x_one` then the first half
    of `x_three`: the cut document prints nothing), `x_lead` (white space before the document), `x_crlf`, `x_garbage`
    and `x_thenasn` (text or an ASN.1 archive after the document end the run), `x_emptyit` (no `Iteration`),
    `x_xml2` (`-outfmt 16`, BLAST XML2: no session), `x_mixpn` / `x_mixnp` (`x_mfirst` and `x_ntr94had` in both
    orders: the motif state crosses the two documents).
- For every PS16d case, NCBI's output with these data equals its output with the full v3.22 files (with the preview
  `cdtrack.txt` for the `*prev` sets), and bioconda's rpsbproc 0.5.1 gives the same bytes.
- Made by hand (PS16d's implementer, 2026-11-13): `arch/x_h_<name>.xml.gz`, edits of `x_one`, `x_three` and `x_ntr67`
  for the edges of NCBI's XML reader and of `ParseBlastOutput` (generator `PR/tmp/impl16d/scripts/gen_h1.py` and
  `gen_h2.py` of the task folder); `oracle/x.h_*` are NCBI's bytes in the fixture layout (both 0.5.1 builds the same,
  and the same with the full v3.22 files). The getters of OPTIONAL members that NCBI's code calls unguarded
  (`x.h_nomatrix`: the session line cut before the matrix; `x.h_noalen`, `x.h_noframe`; `x.h_noalen.e`: not reached for
  an HSP above `-e`), libstdc++'s `basic_string::replace` out of range (`x.h_from0`: a hit from residue 0; `x.h_qlen2`,
  `x.h_ntqlen2`: queries shorter than their hits), gaps and leading or trailing gaps in the alignment strings
  (`x.h_qgap`, `x.h_gaps2`, `x.h_lead`; `x.h_hshort1`: `Hsp_hseq` one letter short, its terminating 0 read as a
  residue), a query string that the alignments make longer (`x.h_qlen50`), entity and character references
  (`x.h_defent`: UTF-8 of `&#233;` and `&#x1F600;`, `&#0;` dropped; `x.h_defbigref`: a reference of 2^32 - 1 ends the
  text and a CDATA section follows), CDATA (`x.h_cdata`), NCBI's end-of-line handling (`x.h_defcr`), the declared
  encoding (`x.h_latin1`, `x.h_cp1252`: bytes from 0x80 converted to UTF-8; `x.h_ebcdic`: unsupported, no session),
  a comment between members (`x.h_comment`) and inside text (`x.h_defcomment`: the document fails), an unknown
  member, a missing mandatory member and two members out of order after a good document (`x.h_unknown`,
  `x.h_missmand`, `x.h_order`: one session), `<Iteration_hits/>` (`x.h_selfhits`), a byte order mark right after the
  first document's last `>` (`x.h_bom2`: skipped, two sessions; the PS16d review's `r_bom2tight`) and after a line end
  (`x.h_bom2nl`: not skipped, one session), `-2147483648` and `2147483648` as
  INTEGER values (`x.h_intmin`: read as INT_MIN, unused; `x.h_intover`: overflow, no session), an empty REAL
  (`x.h_evalempty`: 0), `Hit_accession` beyond `long` and with leading white space (`x.h_acchuge`, `x.h_acclead`),
  a score beyond `int` (`x.h_scorehuge`: NCBI's `(int)` of it is never read), frames on a protein query
  (`x.h_frame1`, `x.h_frame4`, `x.h_framem1`; `x.h_tframe2nd`: frame 1 on the first two HSPs of `x_three`) and frame
  0 and 4 on a nucleotide one (`x.h_ntframe0`, `x.h_ntframe4`), a query length below the hits' ends
  (`x.h_ntqlen1500`).
- NaN values (PS16d review S1-1 and S3-3; NCBI's `HUGE_VAL/HUGE_VAL`, printed `-nan`): `Parameters_expect` `nan`
  (`x.h_nanexpect`, the review's `x_nan`), an archive's E-value threshold `NOT-A-NUMBER` (`nan_threshold`,
  `arch/bad/a_nan.asn`, the review's `a_nan`) and a hit's `e_value` `NOT-A-NUMBER` (`nan_hit`,
  `arch/bad/a_nanhit.asn`, both edits of `arch/one.asn`), a kept `Hsp_evalue` or `Hsp_bit-score` `nan` (`x.h_evalnan`,
  `x.h_bitnan`), every third `Hsp_evalue` of `x_three` (`x.h_tnan3`), and in `x_ntr67e10` (47 HSPs, 18 and 22 in frames
  +1 and -1: NCBI's introsort partitions with NaN keys) every third `Hsp_evalue` (`x.h_ntnanthird`), every one
  (`x.h_ntnanall`) and every fourth `Hsp_bit-score` (`x.h_ntbitnan`). NCBI's bytes in the fixture layout, both builds
  the same, cut data = full data.
- `refused/refuse_x*`: LOSAT's bytes for the hand-made XML that it refuses (PD-LOSAT-NCBI-DEFECTS 1.17): an
  attribute (`refuse_xattr`), a namespace prefix (`refuse_xns`), a control character (`refuse_xctrl`; NCBI replaces it
  and posts a warning), an `Hsp_hseq` two letters short (`refuse_xhshort2`),
  a protein query of length -3 (`refuse_xqlenneg`: NCBI's string of 2^32 - 3 bytes, `std::bad_alloc` by the address
  space), a nucleotide query of 1 base (`refuse_xntqlen1`: frames of 1,431,655,765 residues) and of 2 bases
  (`refuse_xntqlen2`: empty frames, whose residues NCBI writes over the strings' terminator, undefined). Each is NCBI's
  standard output cut at a line end (before the document's session for what the reader refuses, before the query
  for the rest) and its standard error up to the same point, then the refusal.
  (`refuse_xcoordmin` of PS16d is the `widened` row `x.h_qfmin.doms` since PS16-Q1, below.)

## Coordinates past NCBI's `int` (PS16-Q1, PD-LOSAT-NCBI-DEFECTS exception 11)

The Owner's Q1 = B (2026-11-14): where rpsbproc 0.5.1's signed `int` coordinate arithmetic overflows (undefined; its
two builds wrap it), LOSAT computes the value. Each case at the limit sits beside its nearby input, whose coordinates
are moved away from the limit by a known shift; the nearby row is a `compare` row (NCBI's bytes), and LOSAT's bytes of
the limit equal them moved back by that shift (the input file's name aside). Made 2026-11-14 by `gen_oracle.py` and
`gen_fixrows.py` of the task folder (`psiblast-20261025/PR/tmp/impl16q1/`, `oracle/INDEX.md` there); NCBI's binary and
bioconda's gave the same bytes in every case.

- `widened/<case>.{out,err,rc}`: LOSAT's bytes for the rows of `use` `widened`, where NCBI's builds print the wrap:
  - `ntplusmax.doms.std`, `ntplusmax.doms.rep.f` (`arch/ntplusmax.asn`: `nt67.asn` whose first alignment, plus strand,
    ends at base 2^31, `to 2147483647`; NCBI prints the end `-2147483648`, LOSAT 2147483648): nearby
    `ntplusmax3.*` (`to 2147483644`, end 2147483645; shift 3).
  - `ntminusmax.doms.std` (`nt67.asn` whose first minus-strand alignment ends at `to 2147483647`: `alignedLen` of
    715,827,882 codons, the end printed 2147483648, NCBI `-2147483648`): nearby `ntminusmax1.doms.std` (`to
    2147483646`, the same `alignedLen`, end 2147483647; shift 1).
  - `ntmotifmax.feats.std` (`nt68.asn` whose first sd00045 alignment is moved 2,147,483,047 bases, so that its first
    ANK motif starts at base 2^31; data `data_with/ntmot/`, `nt` without sd00045's `oligomer interface` site, so that
    the moved hit prints its motifs only; NCBI prints that motif's start `-2147483648`, LOSAT 2147483648): nearby
    `ntmotifmax3.feats.std` (moved 3 bases less: every number of the three motifs 3 lower).
  - `x.h_qfmin.doms` (`arch/x_h_p_coordmin.xml`, `Hsp_query-from` -2147483648; NCBI's `replace` throws at
    position 2147483647, LOSAT's at 18446744071562067967): nearby `x.h_qfmin1.doms` (`x_h_p_qfmin1.xml`,
    -2147483647: position 18446744071562067968; shift 1).
  - `x.h_qtmin.rep` (`x_h_p_qtmin.xml`, `Hsp_query-to` -2147483648: NCBI's wrapped range to 2^31 - 1 makes the hit a
    representative; LOSAT's ends below its start, so it prints no DOMAINS): nearby `x.h_qtmin1.rep` (-2147483647, the
    same bytes).
  - `x.h_hfmin.std` (`x_h_p_hfmin.xml`, `Hsp_hit-from` -2147483648: NCBI's wrapped domain start marks the hit `N`,
    LOSAT's -2^31 - 1 does not): nearby `x.h_hfmin1.std` (-2147483647, the same bytes).
- `compare` rows whose limit gives NCBI's bytes too (the overflow ends in an `unsigned int` or is not printed):
  `ntslen.both.std` / `ntslen1.both.std` (`nt67.asn` whose first subject interval runs from 2^32 - 5 to 2^31 - 2,
  `t - f + 1` = 2^31 + 4, and to 2^31 - 7, 2^31 - 1), `x.h_qtmin.full` / `x.h_qtmin1.full` (the hit of `x_h_p_qtmin.xml`
  in `-m full`: its end `-2147483648` in both), `x.h_hfmax.both` / `x.h_hfmax5.both` (`x_one.xml` with a query gap at
  column 11 and `Hsp_hit-from` 2147483638, so that `++hf` passes 2^31 - 1, and 2147483633), `x.h_span.doms` /
  `x.h_span1.doms` (a query of length 0, `Hsp_query-from` 0 and `Hsp_query-to` 2147483647: "see if NA"'s span 2^31,
  and 2147483646).
