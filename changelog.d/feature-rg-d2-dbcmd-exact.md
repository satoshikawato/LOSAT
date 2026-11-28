- `blastdbcmd -outfmt` with `%l`: LOSAT now reads a sequence's length only for a line it writes, at `%l`, as NCBI
  BLAST+ does. A damaged sequence whose length cannot be read no longer makes LOSAT refuse when none of its lines is
  written (every definition line removed by the alias files, or a header of no bytes); the other sequences' lines are
  written as in NCBI BLAST+. An unknown specifier before `%l` now gives NCBI BLAST+'s error after the part of the line
  before it.
- `blastdbcmd -range 0-` is an empty range, as in NCBI BLAST+: `-outfmt %s` of a sequence at offset 0 of an empty `.nsq`
  now gives NCBI BLAST+'s error `could not get sequence or range.` where LOSAT refused.
- `blastdbcheck`: a nucleotide at offset 0 of an empty `.nsq` is now reported as NCBI BLAST+ reports it (`sequence
  failed; oid=0: sequence data is empty`, exit status 1), where LOSAT refused the database.
