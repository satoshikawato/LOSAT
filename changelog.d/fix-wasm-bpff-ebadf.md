- Wasm (threaded command on Node): `blastp -task blastp-fast` with `-out` and 2 or more threads failed
  with "Bad file descriptor (os error 8)": the query blocks searched two at a time are written by
  whichever thread finishes the block that completes the order, and each thread of the Node host
  (`LOSAT/tests/run_losat_wasi_threads.js`) had its own table of files opened for writing. The host
  now opens the regular files, devices (such as `/dev/null`) and named pipes that the command opens
  for reading or writing itself, and keeps them in one table that every thread shares, as a native
  process has one table of descriptors. Outputs are those of the native command; native builds and
  LOSAT Web (which opens no files) are not affected.
- Wasm (serial and threaded commands on Node): `-query /dev/stdin` and `-out /dev/stdout` failed with
  "File is not accessible" when standard input or output was a pipe, at every thread count, since
  Node's WASI could not follow `/dev/stdin` to a pipe. Both runners now open such a file by its name
  as the operating system does, so they work as in the native command.
