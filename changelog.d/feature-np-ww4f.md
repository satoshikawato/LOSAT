- LOSAT Web's check of a tool's command line now ends as the run ends for the files that the run
  opens before it reads one: an `-out` that cannot be opened (`nosub/m.asnb` in a folder that the run's
  folder does not have, `input.fa/../m.asnb`, `m.asnb/`, a folder) shows NCBI's
  `File is not accessible` error in the form for dustmasker, segmasker, convert2blastmask and
  windowmasker's masking, a missing or unreadable `-in` shows the run's error in the form for the
  same four, and so does dustmasker's `-hard_masking` with an output format other than FASTA.
  Before, the check accepted these command lines and the run then failed. makeblastdb's check shows the error
  of a `-logfile` that cannot be opened, and of the database's name and title, before the run.
- `losat makeblastdb -in in.fa/` (a file's name with a trailing `/`) runs as NCBI BLAST+ 2.17.0 runs
  it: the warning that the input does not appear to be FASTA, no sequences, "No volumes were
  created." and exit status 0, where LOSAT said that the file does not exist (exit status 1).
