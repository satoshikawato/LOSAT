- The NCBI C++ Toolkit options `-version-full`, `-version-full-xml` and `-version-full-json` run in every program, as
  `-version` does: NCBI BLAST+ 2.17.0's text, XML and JSON forms with the first line `blastn: LOSAT 0.2.0` (NCBI writes `blastn: 2.17.0+`) and LOSAT's
  package, without the lines that NCBI's build writes about itself. They were refused.
- `-version`, `-h` and `-help` no longer stop with a panic (exit status 134) when the standard output cannot be written
  (`-version > /dev/full`): the text is dropped and the program exits 0, as NCBI's does.
- The library API's search programs and the Web engine refuse `-version` and the full versions by name; `-version=1` is
  an unknown option, as in NCBI.
