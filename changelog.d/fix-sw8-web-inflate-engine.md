- Compressed inputs: the gzip, bzip2 and zstd decoders of `blastinput::compressed` take the file's
  bytes in pieces of any size (`Inflater`, in the shape of NCBI's `CCompressionProcessor::Process`
  and `Finish`); the command line's reader runs them, and `PushInflater` lets LOSAT Web's adapter run
  the same decoders on the pieces it reads, so a `.gz` or `.zst` source of LOSAT Web will give the
  command line's bytes and refusals (the Web side follows in the same pull request). The sentence of a
  file that is not what its name says has one producer (`Compression::not_data`). The command line's
  output is unchanged.
