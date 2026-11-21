- `makeblastdb -input_type blastdb` copies existing BLAST databases into a new one, as NCBI BLAST+
  2.17.0 does: one or several databases (alias names too), into format version 4 or 5, with their
  definition lines, PIGs, masks and other columns, the input's title by default, and `-taxid`,
  `-taxid_map`, `-mask_data`, `-hash_index`, `-oid_masks` and `-max_file_sz` as for FASTA input.
  LOSAT refused the input type.
