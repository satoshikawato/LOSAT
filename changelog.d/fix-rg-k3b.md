- FASTA records of 4,294,967,296 letters or more (gaps included), whose length NCBI BLAST+ 2.17.0
  counts in 32 bits: `losat makeblastdb` skips them, stops at their length or reports NCBI's error
  for a protein delta as NCBI does, and refuses only one whose letters NCBI would build (it no
  longer refuses a record that NCBI skips for its gi); `losat convert2blastmask` opens no interval
  at position 4,294,967,295, as NCBI; the refusals name the record.
