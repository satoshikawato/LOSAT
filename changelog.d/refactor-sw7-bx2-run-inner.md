- Engine: the set-up of a BLASTX run's options (with its search strategy file and the
  requests of its archive), the read of its `-subject` file and the archive of a query batch
  (outfmt 11) move out of `run_inner` into their own module. No output changes.
