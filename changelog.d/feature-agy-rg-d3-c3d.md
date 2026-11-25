- `SERIAL_VERIFY_DATA_READ` and `SERIAL_VERIFY_DATA_GET` set to `NO`, `NEVER`, `DEFVALUE` or
  `DEFVALUE_ALWAYS` (any case), which turn off NCBI BLAST+'s checks of a missing mandatory member
  and of an unset member, are now taken into account by every program that reads such data, not
  only `rpsbproc`: an `-in_pssm` file of `psiblast` or `tblastn`, the profiles and binary
  definition lines (`ASN1_BlastDefLine`) of `makeprofiledb`, a `-mask_data` file of `makeblastdb`,
  and the definition lines of a database read by the search programs, `blastdbcmd`,
  `blastdbcheck`, `makeblastdb -input_type blastdb`, `dustmasker`, `segmasker`, `windowmasker` and
  `makembindex`. Well-formed inputs give the same output as before with every value. An input that
  fails a check these values turn off, where NCBI BLAST+ reads on, is refused where NCBI's default
  check stops (for `blastdbcheck`, where it fails the check of that OID), after the same output.
- Database searches: a definition line that lacks its Seq-id member altogether (not one with an
  empty list of Seq-ids) stops the search with NCBI BLAST+'s read error, as NCBI's does, instead of
  its hits being dropped without a message.
