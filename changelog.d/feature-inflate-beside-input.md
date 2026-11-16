- A search program inflates a `-query` or `-subject` named `.gz`, `.bz2` or `.zst` into an unnamed
  temporary file beside the input, in its directory, instead of under `TMPDIR` (memory on a tmpfs
  `/tmp`); when that directory refuses the file or has no room, the copy goes under `TMPDIR` as
  before. A copy that cannot be written is reported as such, not as bad compressed data.
