# Compressed input fixtures

- `mask_nucl.fna.zst`: `../seqdb/mask_nucl.fna` compressed by the Zstandard CLI v1.5.7
  (`zstd -19 -q -c ../seqdb/mask_nucl.fna > mask_nucl.fna.zst`): one frame of one compressed
  block (Huffman-coded literals and FSE-coded sequences) with its content size and an XXH64
  checksum. LOSAT only decodes zstd, and the tests' own frames
  (`blastinput::compressed::test_data::zstd_frames`) hold raw blocks without a checksum, so this
  is the only input whose intact form reaches the decoder's entropy stages and the checksum
  check of `ZstdFrames`.
- `prot.faa.zst`: `../seqdb/prot.faa` compressed by the same CLI from standard input
  (`zstd -19 -q -c < ../seqdb/prot.faa > prot.faa.zst`), the form that `cat f | zstd` writes: a
  frame without a content size (a Window_Descriptor instead, 8 MiB) and with an XXH64 checksum.

Both are read by the tests in `LOSAT/src/blastinput/damaged_input_tests.rs`.
