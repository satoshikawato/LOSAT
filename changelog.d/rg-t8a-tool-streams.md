- `losat dustmasker`, `losat segmasker`, `losat windowmasker`, `losat convert2blastmask` and `losat blastdbcmd`
  write their output as NCBI's programs do: with `-out` naming the standard output and the messages in the same
  pipe, the messages come between the output's bytes where NCBI's do, also for records longer than 8 KiB; when the
  reader closes the pipe early, the messages NCBI writes before it stops are written too. `losat windowmasker
  -mk_counts -infmt blastdb` writes the error of a damaged taxonomy database once, as NCBI does, not once per pass
  over the input.
