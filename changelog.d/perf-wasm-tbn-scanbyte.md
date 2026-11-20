- Web (WebAssembly): TBLASTN's subject scan does less work for each word and for each short list of hits: a
  1-thread search of one human protein against the human genome takes about 10 % less time in Node, and about
  2.5 % fewer cycles in the native x86-64 build (BUSCO's eukaryota proteins against the yeast genome: about 1 %
  fewer); the reports do not change.
