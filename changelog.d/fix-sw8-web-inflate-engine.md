- Compressed inputs: the gzip, bzip2 and zstd decoders of `blastinput::compressed` take the file's
  bytes in pieces of any size (`Inflater`, in the shape of NCBI's `CCompressionProcessor::Process`
  and `Finish`); the command line's reader runs them, and `PushInflater` lets LOSAT Web's adapter run
  the same decoders on the pieces it reads (LOSAT Web, below). The sentence of a file that is not what
  its name says has one producer (`Compression::not_data`). The command line's output is unchanged,
  except that a read error in the middle of a `.zst` is no longer wrapped in ruzstd's own text (as for
  `.gz` and `.bz2` already). A zstd decoder reserves its window at the first frame, as at every later
  frame: `makeblastdb -in` of a `--long=27` `.zst` peaks at 102 MB instead of 140 MB.
- LOSAT Web: a `.gz` or `.zst` source is refused exactly as the command line refuses it, with the
  command line's sentence (`'x.fa.zst' is named .zst but is not zstd data (a zstd frame's checksum
  does not match its data); LOSAT reads a file named .zst through zstd`), and otherwise inflates to
  the command line's bytes; the web app ships two fewer packages (fflate and fzstd). Before, a zstd
  file with a wrong checksum or content size was read without an error, and 20 of 87 gzip edge cases
  read differently from makeblastdb. The app inflates through the engine's decoders over three new
  exports of ABI v2 (`losat_web2_inflate_begin`, `losat_web2_inflate_chunk`,
  `losat_web2_inflate_end`; the inflated bytes on stream 5; the version stays 2), and the engine's
  name rule gives the compression and the name of the inflated file. A downloaded Example set
  takes the same path (the engine's decoders in the Data worker instead of the browser's
  `DecompressionStream`): a set served as several gzip members, which every browser refused, now
  loads, and a damaged one fails with the engine's sentence.
  Each compressed file inflates on a reactor instance of its own, let go when it ends. On a 96 MB
  FASTA, inflation is 1.3 to 4.8 times faster than before in Chromium, Firefox 157 and WebKit (zstd
  with a long window up to 24 times) and holds less memory, except Chromium with zstd (peak +22 %,
  +30 % with `--long=27`).
