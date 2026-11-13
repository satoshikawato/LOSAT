# Mask Seq-ids matched as NCBI's `CSeq_id_Handle` (`-parse_seqids` databases)

Fixtures of the test `writedb::mask_data::tests_ids` (`src/writedb/mask_data/tests_ids.rs`): which
database sequence a mask file reaches, in NCBI BLAST+ 2.17.0 `makeblastdb -mask_data`, and the rules
of `seqid::mask_handle` that reproduce it.

## Files

- `cases.tsv`: one row per case, tab separated: name, FASTA file (in `fasta/`), kind (`n` nucleotide,
  `p` protein), `-parse_seqids` (1/0), mask files (comma separated, in `case/<name>/`, given to
  `-mask_data` in that order), the algorithm ids to read (`%m<id>`), the fate, NCBI's answer.
- `fasta/`: the sequences. `case/<name>/m*.asnt`: the mask files (`Blast-db-mask-info`, ASN.1 text), one Seq-id
  (or a few) in the Seq-loc of each.
- Fates: `match` = LOSAT gives NCBI's answer; `refuse` = NCBI's answer depends on memory addresses and
  LOSAT stops with an error that names what it cannot decide; `env:...` = NCBI's answer needs
  `OBJECTS_PACK_TEXTID=0` or `OBJECTS_PACK_GENERAL=0`, which LOSAT refuses (not run).
- The answer is `NOMATCH` (makeblastdb exits 0 and prints "No sequences matched any of the masks
  provided.") or the sequences that have a mask: `oid|ranges of each algorithm`, joined by a space; a
  range is `from-to;` as `blastdbcmd -outfmt %m` prints it (`to` is exclusive; `none` = no range).

## How they were made (2026-11-13, NCBI BLAST+ 2.17.0, `BLAST_USAGE_REPORT=false`)

- Cases `n*` and `p*` (98): hand-written mask files, each holding one Seq-id of one alternative of the
  `Seq-id` CHOICE (accession, name, version, release, case, packing, `general` tags, PDB chain and
  release date, patent numbers, `giim`, several ids of one sequence, two files), probing which
  database Seq-id of `fasta/nuc_main.fa`, `prot_main.fa` and the small files around them it reaches.
  Each was run as

  ```
  makeblastdb -in <fasta> -dbtype <nucl|prot> -parse_seqids -mask_data <files> -out <db>
  blastdbcmd -db <db> -entry all -outfmt '%o|%i|%m<algorithm id>'
  ```

  (`-blastdb_version` unset: version 5). The rows of `blastdbcmd` and the exit status and standard
  error of `makeblastdb` are the answers. `n53`-`n55` ran with `OBJECTS_PACK_TEXTID=0` /
  `OBJECTS_PACK_GENERAL=0` in the environment.
- Cases `v*` (28): one variant of a Seq-id per file (`lcl`, `gnl`, `ref`, `gb`, `emb`, `pdb`, `sp`, `gi`
  forms) against `fasta/idprobe_ids.fa`, same commands.
- `A_dust_parsed`, `B_nuc_dust_parsed`, `B_prot_seg_parsed`, `N_dust_plain`: the masks that
  `dustmasker -parse_seqids -outfmt maskinfo_asn1_text` and `segmasker -parse_seqids -outfmt
  maskinfo_asn1_text` write for `fasta/idprobe_ids.fa` and `fasta/idprobe2_*.fa`, and
  `dustmasker` without `-parse_seqids` (`N_dust_plain`, a database without `-parse_seqids`), with
  `makeblastdb -blastdb_version 4` and `5`.

Both format versions give the same answers; the test builds each case in version 5 and 4.
