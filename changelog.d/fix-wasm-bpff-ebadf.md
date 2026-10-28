- Wasm (threaded command on Node): `blastp -task blastp-fast` with `-out` and 2 or more threads failed
  with "Bad file descriptor (os error 8)": the query blocks searched two at a time are written by
  whichever thread finishes the block that completes the order, and each thread of the Node host
  (`LOSAT/tests/run_losat_wasi_threads.js`) had its own table of files opened for writing. In a
  command that requests 2 or more threads, the host now also opens a file that is opened for writing
  (a regular file, a new one, `/dev/null`, `/dev/stdout` when standard output is a file, a terminal
  or a pipe, or a named pipe) itself and keeps it in the table that every thread shares, as it does
  for files opened for reading. Results do not change; one-thread runs, native builds and LOSAT Web
  (which writes no files) are not affected.
