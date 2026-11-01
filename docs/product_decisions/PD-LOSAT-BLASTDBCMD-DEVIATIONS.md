# Product Decision: blastdbcmd outputs that NCBI gets wrong

- Decision ID: `PD-LOSAT-BLASTDBCMD-DEVIATIONS`
- Version: 1.2
- Date: 2026-10-30 (this record; items 1 and 2 served by `losat blastdbcmd` from the pull
  request of line BDC's BDC2a). 1.1, 2026-11-01 (BDC2c): the evidence of items 1 and 2 on
  version 4 databases with parsed ids, whose ids `losat blastdbcmd` now looks up; the
  decision is unchanged. 1.2, 2026-11-01 (BDC3a: the taxonomy, the lists, the metadata and
  the IPGs): the explicit rejections of the scope; the four items are unchanged. The deviations were approved by the Owner on 2026-10-29
  ("これはあえて直してお出しするのは？TBLASTX みたいに", drop-in survey `DECISIONS.md` 12), as
  the TBLASTX local-subject genetic code was
  ([`PD-TBLASTX-LOCAL-SUBJECT-GENCODE`](PD-TBLASTX-LOCAL-SUBJECT-GENCODE.md)).
- Status: Accepted (approved deviations). Items 3 and 4 are approved and will be served with
  the masks (line BDC's BDC3); until then those inputs are rejected.

## Scope

`losat blastdbcmd` only. Four outputs of NCBI BLAST+ 2.17.0 `blastdbcmd` whose content is
wrong (the letters or the sequences do not belong to the header or the id asked for),
confirmed with NCBI's binary (line BDC's frozen probe, `dropin-survey-20261029/BDC/probe`).
Everything else `losat blastdbcmd` writes equals NCBI's bytes, under the root
[`AGENTS.md`](../../AGENTS.md) bit-perfect rule; what it does not serve is an explicit
rejection ("... is not supported by LOSAT's BLASTDBCMD"). NCBI's deterministic errors on
the same inputs (an exception text, an exit status) are reproduced, not excepted
([`PD-LOSAT-NCBI-DEFECTS`](PD-LOSAT-NCBI-DEFECTS.md) rule 2); its crash on `-range ""`
(exit status 139) and the `-taxids abc` exception that prints NCBI's build path (exit
status 255) are explicit rejections. So are (BDC3a, probes J and K of the same directory)
the other exceptions that print NCBI's build path (a `-taxids` / `-taxidlist` id that is
empty, blank, too large or not a number: J101, J102, J107, J108, J131; a database date that
NCBI's `CTime` does not read, for `-metadata`: J189, J190), NCBI's crash on an empty index
file of `-list` with `%t` (J168), the reversed protein letters NCBI writes for `-strand
minus` on a protein database read with a taxid list and `-dbtype guess` (J114: NCBI then
takes the database as nucleotide for the strand), `-ipg` / `-ipg_batch` on a database
without PIG files, where NCBI writes an OID it never set (J204, J205, J212, J214), and
`-ipg_batch` on a database opened with `-target_only` and a taxid list (K38). `-taxids
1` follows [`PD-LOSAT-NCBI-DEFECTS`](PD-LOSAT-NCBI-DEFECTS.md) exception 5.

## Decision

In each item the expected output is NCBI's own output for the neighbouring input that does
not reach the defect, so that the test oracle stays NCBI's binary.

1. **`-long_seqids` with `-range` on the minus strand (served).** NCBI's `CFastaOstream`
   path reverses the range location and sets `fReverseStrand` as well: the two reversals
   cancel (`seq_map_ci.hpp` `GetRefMinusStrand`) and the forward letters are written under
   the header `:c<to>-<from>` (probe S18 on a database with parsed ids; D16, D82, D87, D97,
   D101 on others). The short-header path of the same command without `-long_seqids`
   reverses once (E44, D83, D88). LOSAT writes the header of the `-long_seqids` path with
   the reverse complement of the range, the letters of the short-header run of the same
   range. Expected bytes: NCBI's `-long_seqids` header with the letters of NCBI's run
   without `-long_seqids` (`dropin-survey-20261029/BDC/probe/make_expected.py`,
   `LOSAT/tests/blastdbcmd/expected/D82 D87 D97 D101`).
2. **`gnl|BL_ORD_ID|n` in a version 4 database of several volumes (served).** NCBI
   resolves the id to the ordinal `n` without reading a file and then adds `n` to the first
   ordinal of every volume, keeping each sum that falls inside the database: `-entry
   'gnl|BL_ORD_ID|3'` on a database of four volumes writes four sequences (probe S01, and
   D57 with `-entry_batch -get_dups`; an alias of two databases, a partial alias and a
   database with parsed ids alike: I795, I816, I830, I849), while `%a` names the sequence
   `BL_ORD_ID:3` once.
   LOSAT writes the one sequence `n` of the whole database, as a database of one volume
   gives it. Expected bytes: the first record of NCBI's output (S01), or NCBI's run of the
   same batch without `-get_dups`, whose bulk lookup keeps the first ordinal alone (S02 for
   D57).
3. **`-mask_sequence_with` with `-strand minus` in FASTA (BDC3).** NCBI maps the mask to
   forward coordinates and lower-cases the output positions of the reversed letters (probe
   M04, S31); its own `-outfmt %s` path masks first and reverse-complements after (S32).
   LOSAT will mask and then reverse-complement, as `%s` does.
4. **Two or more `%m` in `-outfmt` (BDC3).** NCBI reads the algorithm id of the first `%m`
   only and leaves the digits of the second in the text (probe M03). LOSAT will write each
   `%m<id>` with its own algorithm, as the runs with one `%m` each give.

## Evidence

- Line BDC's frozen probe (NCBI 2.17.0 on the workstation, 2026-10-29, 2026-10-30 and
  2026-11-01):
  `dropin-survey-20261029/BDC/probe/out/{S01,S02,S18,S31,S32,M03,M04,D16,D57,D82,D83,D87,D88,D97,D101,I771,I795,I816,I830,I849}.*`,
  the code reading in `dropin-survey-20261029/BDC/notes/B_fasta_outfmt.md` (ORACLE 1 and 3)
  and `notes/A_id_resolution.md` (§2).
- The tests of items 1 and 2: `LOSAT/tests/blastdbcmd_retrieve.rs`
  (`long_ids_on_the_minus_strand_write_the_reverse_complement`,
  `ordinal_ids_of_a_database_in_volumes_are_database_wide`).
