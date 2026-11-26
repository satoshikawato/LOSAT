# blastdbcheck fixtures (row RG-D2b)

What they are: NCBI BLAST+ 2.17.0 `blastdbcheck` runs that `tests/blastdbcheck_fixtures.rs` repeats with
`losat blastdbcheck`, each row in its own fresh run directory.

## Origin
- Program: NCBI BLAST+ 2.17.0 `blastdbcheck` (the official `ncbi-blast-2.17.0+-x64-linux` build; its `makeblastdb` and
  `makeprofiledb` made the databases of `files/ord`, `files/ca`, `files/cdd` and `files/two`, in the D2b survey probe
  of 2026-11-24, with the wall clock frozen at 1700000000), run on the local development host (WSL2) by the RG-D2b
  oracle on 2026-11-25: `build_tree.py` and `run_ncbi.py` of `remaining-gaps-20261112/D/d2/d2b_oracle/` (outside the
  repository), under its `oracle.lock`. Each row ran in a fresh run directory holding only `in/`, `home/` and the
  row's case directory, as `flock ... env -i PATH=/usr/bin:/bin TZ=UTC LC_ALL=C BLAST_USAGE_REPORT=false
  HOME=<RUN>/home <row env> LD_PRELOAD=fixtime.so FIXTIME=<seed> prlimit --as=4294967296 --core=1 timeout 60
  blastdbcheck <args>` (timeout 5 s for `stride0` and `strideneg`), standard input `/dev/null`, working directory
  `<RUN>/cases/<dir>`. `fixtime.so` is an `LD_PRELOAD` shim whose `time`, `gettimeofday` and `clock_gettime` return
  `FIXTIME` seconds. `FIXTIME` is 1700000000, except where `cases.tsv` gives a `seed`. The rows without a seed were
  also run with `FIXTIME=1234567`: the same bytes, so their output does not depend on the time.
- Inputs: the base databases and lists are `tests/fixtures/blastdb_aliastool/in/` (its README names them); this
  directory does not copy them. `files/` holds only what `in/` lacks:
  `tax/<variant>/taxdb.{bti,btd}` (`good`, `no_9606`, `empty_fields`, `tabs0`, `high_bytes_utf8`, `bad_magic`: the
  taxonomy files of `tests/seqdb/taxdb*`); `edited/` (`p4.pin`/`p5.pin` with an empty date, `p4.pin` with total
  length 0, `p4.pin`/`p5.pin` with no sequences: the index header edited once and stored); `cdd/db.*` (a
  profile database of three profiles, `makeprofiledb -in l_p3.lst -out db`, the stand-in for an RPS database) and
  `two/{a1,b2}.{pin,pdb}` (the two files of the same database that differ when it is named `a1` or `b2`); `ord/o4.*`,
  `ord/o5.*` (a protein database without `-parse_seqids`, version 4 and 5); `ca/ca4.*`, `ca/ca5.*` (a few sequences, the first
  with two deflines joined by a control-A, `-parse_seqids`, taxid map only for the second defline: 9606); `alias/gi.gil`
  (the 20-byte binary GI list of `in/gi_ok.bgl`); `junk/` (an 8-byte and a 12-byte `.pin`); `negtax/neg.*` (one
  sequence of two definition lines joined by a control-A, `gi|5|ref|NP_000001.1|` without a taxid and
  `gi|99|ref|NP_000009.1|` with taxid -5, `makeblastdb -dbtype prot -parse_seqids -taxid_map -blastdb_version 4
  -title neg`, clock frozen at 1700000000, made on 2026-11-25 for review F3 of row RG-D2b).

## Rows
`cases.tsv`: `id`, `dir`, `env`, `seed`, `args`. `dir` is the case directory under `cases/` that is the working
directory; `env` is extra environment (`BLASTDB=../../in`) or `-`; `seed` is the `FIXTIME` of the rows whose output
depends on the time-seeded random sample (`-random N` with 0 < N < 200: `combo`, `fx_rnd3a`, `fx_rnd3b`, `fx_dir2`),
else `-`; the test skips the rows with a seed (a crate-internal test covers them); `args` are bash words with
single quotes only, `<RUN>` standing for the run directory.

Families (215 rows, 84 case directories): command line and frame (`help`, `xmlhelp`,
`noargs`, `arg_*`, `m_bad_lf`); verbosity, sampling and selection on `p5` (`t_p5*`, `t_*_full_v4`, `rnd*`, `ends*`,
`stride*`, `combo`, `fx_*`, `noisam`, `taxids_*`, `cdd`, `t_abs_p5`, `t_dir_v4`); MetaData (`nt_p5*`, `t_*` taxdb
variants, `lf0`, `lf1`, `*legacy*`, `d_*_date`, `d_p4_totlen0`, `d_*_zero_oids*`); damaged volumes (`d_p4_*`,
`d_p5_*`, `d_n4_*`, `d_n5_*`, `d_pmv_*`, `pp_*`, `d_p4_pnd_swap`); names, lists and directories (`m_*`, `r_*`,
`nest`, `m_blastdb_*`); alias checks (`a_*`, `a_dir`, `a_dir_v1`, `a_ov1`, `a_ov2`, `a_ov3`, `a_signs`, `a_giabs`,
`a_msk1`); review rows of 2026-11-25 (`mix1`, `mix2`: an alias file with the other end-of-line byte inside a line;
`dirslash`, `dirslash2`: `-dir` with trailing separators; `mix45_*`: alias files over version 4 and version 5
volumes, in both orders; `d_p4_pin_nucl_prot`, `d_p4_pin_nucl_dir`: an index of the other molecule type; `negtax`,
`negtax_v3`: `-must_have_taxids` with a negative taxid on a later definition line);
CDD (`c_*`, including `c_nucl`, `c_noloo`, `c_norps`, `c_noobsr`, `c_*magic`, `c_rps2`, `c_freq2`, `c_rps0`,
`c_two_diff`); widened overflow and shared-reader rows of 2026-11-26 (`ends_intmin`, `ends_intmin_near`,
`stride_intmax`, `d_p4_nonasc_seq`, `d_p4_nonasc_hdr`, `d_p4_badchar`, `d_p4_badchar_lf`, `d_p4_utf8`, `a_msk1`);
other versions (`ord4`, `ord5`, `ca4`, `ca5`, `t_pmv_full_v4`).

`legacy`, `legacy_atm`, `legacy_ab` and `nolegacy` each have one extra file next to `p5` (empty `p5.ptm`, non-empty
`p5.atm`, empty `p5.ab` (not matched: two letters), empty `p5.zzz`), so no output depends on the order in which the
file system lists a directory. Other rows with several files per volume are not affected: NCBI sorts what `-dir`
finds.

Rows that crash or never end in NCBI (`d_p4_psq_zero`, `d_p5_psq_zero`, `d_p4_psd_cut`, `r_garbage`, `r_garbage2`: exit
status 139; `stride0`, `strideneg`: the 5 s limit, status 124) are kept as NCBI gave them.

## Setup
`setup.tsv`: one operation per line, `dir`, `op`, arguments, in order; the harness builds `<RUN>/in` (a copy of
`blastdb_aliastool/in`), an empty `<RUN>/home` and the row's `<RUN>/cases/<dir>`:
- `db NAME [SUBDIR [PREFIX]]` copies every file `in/NAME.*` into the directory (or `SUBDIR`, created, `.` = the
  directory), `PREFIX` in front of each name; mode 0644.
- `tax VARIANT` copies `files/tax/VARIANT/taxdb.bti` and `taxdb.btd`.
- `put DEST SOURCE` copies `files/SOURCE` to `DEST` (parents created); `cp SRC DEST` copies a file of the directory.
- `rm NAME`; `trunc NAME SIZE` (a byte count, `-N` = drop the last N bytes, `half` = half the length rounded down);
  `poke NAME OFFSET HEX` (overwrite at the decimal offset); `text NAME TEXT` (`\n`, `\r`, `\t`, `\\` and `<RUN>`
  replaced; no TEXT = empty file); `mkdir SUBDIR`; `link NAME TARGET` (symbolic link, TARGET relative to the link's
  directory).
- `d_p4_pnd_swap` is `poke p4.pnd 15 02` (the OID of the sequence with GI 102 set to 2). The case directory without
  operations (`nt`) is an empty directory: no taxdb.

## How frozen
`ncbi/<id>.out`, `.err` and `.rc`: NCBI's standard output, standard error (raw, NCBI's words; only the run
directory's absolute path replaced by `<RUN>`) and exit status. The test maps the words to LOSAT's with
`losat::wording` (`losat_output_bytes`, `losat_bytes`). `ncbi/<id>.files`: each file below the run directory that the
run created or changed, `relative/path<TAB>size<TAB>FNV-1a 64 digest (16 hex digits)` of its bytes with the run
directory replaced by `<RUN>` (the rows with `-logfile`; empty for the others), and `ncbi/<id>.made/<path>` the bytes
themselves: their exception texts are NCBI's words, which the test maps to LOSAT's before it takes the size and the
digest. The freeze was run twice: the same bytes (the `.made` files were added on 2026-11-25 by a third run of the four
`-logfile` rows, whose other files came out unchanged). The 12 review rows of 2026-11-25 (above) were frozen once more
on that day, in a run of those rows only.

`-version` and `-version-full` are not rows here: `tests/ncbi_toolkit_version.rs` checks them with LOSAT's own version
line (the Owner, 2026-11-21).

## Rows that wait for LOSAT or follow approved exceptions
`LOSAT_REFUSES` in `tests/blastdbcheck_fixtures.rs` names the 8 rows where the test expects LOSAT's refusal
(exit status 1, `Error: ` ... `not supported by LOSAT`) instead of NCBI's bytes because NCBI crashes or never ends
(`d_p4_psq_zero`, `d_p5_psq_zero`, `c_rps0`, `d_p4_psd_cut`, `r_garbage`, `r_garbage2`, `stride0`, `strideneg`).
`APPROVED_WIDENED_OVERFLOW` (`a_ov1`, `a_ov2`, `a_dir`, `a_dir_v1`, `ends_intmin`) names rows covered by approved
exception 32 of `PD-LOSAT-NCBI-DEFECTS` (version 1.113, `Q-D2b-OV=A`), where LOSAT widens signed arithmetic and avoids
NCBI's signed overflow. NCBI's frozen result stays beside every row for provenance.
