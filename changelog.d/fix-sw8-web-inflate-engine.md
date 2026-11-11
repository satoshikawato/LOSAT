- Compressed inputs: the gzip, bzip2 and zstd decoders of `blastinput::compressed` take the file's
  bytes in pieces of any size (`Inflater`, in the shape of NCBI's `CCompressionProcessor::Process`
  and `Finish`); the command line's reader runs them, and `PushInflater` lets LOSAT Web's adapter run
  the same decoders on the pieces it reads (LOSAT Web, below). The sentence of a file that is not what
  its name says has one producer (`Compression::not_data`). The command line's output is unchanged.
- LOSAT Web: a `.gz` or `.zst` source is refused exactly as the command line refuses it, with the
  command line's sentence (`'x.fa.zst' is named .zst but is not zstd data (a zstd frame's checksum
  does not match its data); LOSAT reads a file named .zst through zstd`), and otherwise inflates to
  the command line's bytes; the web app ships two fewer packages (fflate and fzstd). Before, a zstd
  file with a wrong checksum or content size was read without an error, and 20 of 87 gzip edge cases
  read differently from makeblastdb. The app inflates through the engine's decoders over three new
  exports of ABI v2 (`losat_web2_inflate_begin`, `losat_web2_inflate_chunk`,
  `losat_web2_inflate_end`; the inflated bytes on stream 5; the version stays 2), and the engine's
  name rule gives the compression and the name of the inflated file. A `.bz2` source is still
  refused by the web app.
