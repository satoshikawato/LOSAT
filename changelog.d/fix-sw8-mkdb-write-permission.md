- Fixed: `losat makeblastdb` with `-out` in an existing directory it cannot write to now stops
  before the log with makeblastdb's message "You do not have write permissions on '<dir>'" and
  exit status 255, as makeblastdb does, on Linux and macOS. It used to write the log header and
  fail at the first volume file with the system's message. On Windows it still fails at the first
  file (makeblastdb reads the directory's access rights there).
