- The NCBI C++ Toolkit option `-xmlhelp` runs on the command line of every program: it writes NCBI BLAST+ 2.17.0's XML
  description of the program's arguments without the options that LOSAT refuses by name (and the dependencies that name
  them), so a tool that builds its form or wrapper from the XML offers only what LOSAT runs, with NCBI's types, defaults
  and constraints. blastdbcmd, dustmasker, segmasker and convert2blastmask write NCBI's whole XML. As with `-help`, the
  arguments after it are not read. It was refused.
- The library API (and pylosat), the Web engine and the Web app refuse `-xmlhelp` by name.
