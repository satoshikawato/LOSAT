- LOSAT Web reads the name of `-export_search_strategy` against the run's files, as it reads `-logfile`: a name in a
  folder that the run's files do not have (`out/s.asn` without a folder `out`) or in one of their files, and the name
  of one of their folders, is NCBI's `File is not accessible` error, as on the command line run in that folder; before,
  the Web GUI wrote the file. A file of the run that no option names (such as `.ncbirc`) and a file that the search
  reads (such as the alias file of a `-db` database, even where it is missing) are refused in the Web GUI, as for
  `-logfile`; before, the Web GUI wrote a file of that name, where the command line's file takes the place of the run's.
- LOSAT Web names the search strategy file of `-export_search_strategy` as the command line's file, as it names the
  file of `-logfile` and the PSSM files: `./s.asn`, `in/../s.asn` and `/proc/self/cwd/s.asn` give a file `s.asn`,
  shown and saved under that name; before, the Web GUI kept the name as it was typed.
- LOSAT Web gives the psiblast PSSM file of `-out_pssm` or `-out_ascii_pssm`, and the file of `-logfile`, whose name
  passes a folder of the run's files outside the run's directory: `/data/../tmp/p.asn` beside a query `/data/q.fa`
  gives the file `/tmp/p.asn`, and `-logfile /data/../tmp/run.log` the file `/tmp/run.log`; before, the Web GUI failed
  that run after the search with an error of its own (a file that the argv does not name).
- LOSAT Web names a file outside the run's directory whose name has `..` after a link as the file that the command
  line makes: `-export_search_strategy /var/run/../tmp/s.asn` gives `/tmp/s.asn` (`/var/run` is `/run`), as do
  `-out_pssm`, `-out_ascii_pssm` and the `-out` of `windowmasker` and `blast_formatter`; before, the Web GUI named it
  `/var/tmp/s.asn`, another file. `-logfile` keeps NCBI's own reading of its name (`/var/tmp/run.log` for
  `/var/run/../tmp/run.log`).
- LOSAT Web labels each file that a run makes besides its outputs by the option whose file it is, as the run that
  makes it says, also where the file's name is the value of `-logfile` as typed or the name of the `-logfile` file
  that the run does not make: with `-logfile -` (no file) and `-export_search_strategy ./-`, the file `-` is the
  search strategy file, and with `-logfile sub/../s.asn -export_search_strategy s.asn` and no folder `sub`, where the
  command line makes only the search strategy file, the file `s.asn` is that file; before, the Web GUI showed each as
  the `-logfile` file.
