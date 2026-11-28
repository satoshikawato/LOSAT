- LOSAT Web's tool tabs read every spelling of an `-out` or `-logfile` name as the file the run makes,
  as the command line opens it: `-out sub/../m.asnb` (through a folder of the run), `.//m.asnb`,
  `sub//m.asnb` and `/proc/self/cwd/m.asnb` now get the `-out` label of their file and its hand-overs
  ("Use as -mask_data in makeblastdb", "Use in BLASTN", "Mask with these counts"), where only a
  leading `./` was understood; `-logfile` is read as NCBI reads it (`sub/../run.log` is `run.log`).
  An `-out` file outside the run's folder (`../m.asnb`, `/tmp/m.asnb`) goes to makeblastdb's
  `-mask_data` under its name alone, which makeblastdb's form refused as a path before.
