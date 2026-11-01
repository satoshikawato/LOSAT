# Product Decision: blastdbcmd outputs that NCBI gets wrong

- Decision ID: `PD-LOSAT-BLASTDBCMD-DEVIATIONS`
- Version: 1.4
- Date: 2026-10-30 (this record; items 1 and 2 served by `losat blastdbcmd` from the pull
  request of line BDC's BDC2a). 1.1, 2026-11-01 (BDC2c): the evidence of items 1 and 2 on
  version 4 databases with parsed ids, whose ids `losat blastdbcmd` now looks up; the
  decision is unchanged. 1.2, 2026-11-01 (BDC3a: the taxonomy, the lists, the metadata and
  the IPGs): the explicit rejections of the scope; the four items are unchanged. 1.3,
  2026-11-01 (BDC3b, the masks): items 3 and 4 served, items 5 and 6 added on the same rule. The
  deviations of items 1 to 4 were approved by the Owner on 2026-10-29
  ("これはあえて直してお出しするのは？TBLASTX みたいに", drop-in survey `DECISIONS.md` 12), as
  the TBLASTX local-subject genetic code was
  ([`PD-TBLASTX-LOCAL-SUBJECT-GENCODE`](PD-TBLASTX-LOCAL-SUBJECT-GENCODE.md)), and those of
  items 5 and 6 on 2026-11-01 ("差 5・6 の承認: 承認します", drop-in survey `DECISIONS.md` 33).
  1.4, 2026-11-02 (BDC4, several names in `-db`, drop-in survey `DECISIONS.md` 34): the items hold
  on a list of names; item 6 covers one accession on two sequences. **That widening of item 6
  is proposed and waits for the Owner's approval.**
- Status: Accepted (approved deviations); the widening of item 6 in 1.4 proposed.

## Scope

`losat blastdbcmd` only. Six outputs of NCBI BLAST+ 2.17.0 `blastdbcmd` whose content is
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
`-ipg_batch` on a database opened with `-target_only` and a taxid list (K38). So are (BDC3b,
probe M) the exceptions of a damaged column file that print NCBI's build path (M701, M712,
M854, M864), `%m` digits that overflow NCBI's `int` (M143, M144, M145), `-mask_sequence_with`
below -1, which NCBI looks up with an uninitialised id (M225, M227), a sequence whose volume of
an alias has masks but not the algorithm, where NCBI reads the masks of the volume's id before
or an uninitialised one (M180-M184, M255-M270; on a list of names, probe N: N1269, N1306, N1313),
a list of names in which an alias and another name give the same volume, whose sequences NCBI's
`-info` and `-metadata` count twice, 39 for 21 (BDC4: N389, N393, N1627; the search programs
refuse it too), a volume with the big-endian data file of a
column alone, whose offsets NCBI reads as little-endian (M784), another column before the masks,
which NCBI reads as the masks (M824), and `-outfmt %s` without a range on a stored mask outside
the sequence, where NCBI crashes (M786, M886, M896). `-taxids
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
   database with parsed ids alike: I795, I816, I830, I849; a list of names: N846, N890), while
   `%a` names the sequence
   `BL_ORD_ID:3` once.
   LOSAT writes the one sequence `n` of the whole database, as a database of one volume
   gives it. Expected bytes: the first record of NCBI's output (S01), or NCBI's run of the
   same batch without `-get_dups`, whose bulk lookup keeps the first ordinal alone (S02 for
   D57).
3. **`-mask_sequence_with` with `-strand minus` in FASTA (served).** NCBI maps the mask to
   forward coordinates and lower-cases the output positions of the reversed letters (probe
   M04, S31, M348-M392); its own `-outfmt %s` path masks first and reverse-complements after
   (S32, M349-M393). LOSAT masks and then reverse-complements, as `%s` does, with
   `-long_seqids` too (M395-M406, with item 1). Expected bytes: NCBI's headers with the
   letters of NCBI's `%s` run of the same words.
4. **Two or more `%m` in `-outfmt` (served).** NCBI reads the algorithm id of the first `%m`
   only and leaves the digits of the second in the text (probe M03, M130). LOSAT writes each
   `%m<id>` with its own algorithm, as the runs with one `%m` each give; a `%m` without digits
   has the first one's, as in NCBI (M131), and an unknown id of a later `%m` is the error of
   its own run (M133, as M128). The first `%m` of the text may be the literal of `%%m<digits>`,
   where NCBI reads the id for the `%m` after it (M137: LOSAT writes `%%m11 %m12` with the
   masks of 12). Expected bytes: the runs with one `%m` each (M125, M126).
5. **`-mask_sequence_with` on a protein in FASTA without a range (served).** NCBI's mapper
   for a whole sequence (`sequence.cpp` `x_GetMaskingStates`) keeps the masks in the first third
   of a protein: `XP_100002.1` (190 residues), masked at 50-70, is lower-case at 50-64, and
   `P10001.1` (225 residues), masked at 80-105 and 165-186, not at all (probe M05, M285). The
   same sequences with a range are masked in full (M290), and so is `-outfmt %s` (M286). LOSAT
   masks the whole length. Expected bytes: NCBI's headers with the letters of NCBI's `%s` run.
6. **A batch line without masks after one with masks of the same sequence (served).** NCBI's
   FASTA formatter sets its mask only when the sequence has masks of the algorithm, so that
   `NM_100001.1 12`, which has none, is written with the masks of the line `NM_100001.1 11`
   before it (probe M415, M418). LOSAT writes it without masks, as the run of that line alone
   (M414) and `-outfmt %s` (M417) do. Expected bytes: NCBI's headers with the letters of
   NCBI's `%s` run. The formatter's mask stays for the same Seq-id on another sequence as well:
   `-db "db/nm5 db/m/nms" -entry NM_100001.1 -mask_sequence_with 12` writes the second sequence,
   which has no masks of 12, with the masks of the first (probe N, N1813; N1702 in a batch; 1.4).
   LOSAT writes it without masks, as `%s` does (N1814).

## Evidence

- Line BDC's frozen probe (NCBI 2.17.0 on the workstation, 2026-10-29, 2026-10-30 and
  2026-11-01):
  `dropin-survey-20261029/BDC/probe/out/{S01,S02,S18,S31,S32,M03,M04,M05,M125,M126,M128,M130,M133,M137,M285,M286,M290,M414,M415,M417,D16,D57,D82,D83,D87,D88,D97,D101,I771,I795,I816,I830,I849,N389,N846,N890,N1702,N1703,N1813,N1814}.*`,
  the code reading in `dropin-survey-20261029/BDC/notes/B_fasta_outfmt.md` (ORACLE 1 and 3)
  and `notes/A_id_resolution.md` (§2).
- The tests of items 1 and 2: `LOSAT/tests/blastdbcmd_retrieve.rs`
  (`long_ids_on_the_minus_strand_write_the_reverse_complement`,
  `ordinal_ids_of_a_database_in_volumes_are_database_wide`); of items 3 to 6:
  `LOSAT/tests/blastdbcmd_masks.rs` (`masks_are_ncbi_s_bytes`); on a list of names:
  `LOSAT/tests/blastdbcmd_several_names.rs` (`several_names_are_ncbi_s_bytes`).
