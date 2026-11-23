- `losat blastdbcmd -metadata`: when the JSON cannot be written (a full disk, a file size limit), to
  standard output or to `-out`, LOSAT now gives NCBI's messages, leaves the bytes NCBI leaves in the
  file and exits with NCBI's status, instead of an error.
