- `losat makeblastdb -mask_data` reads a mask string with a byte outside the visible ASCII range as makeblastdb does: the
  byte is stored as `#` and NCBI's error for it is written (`Error: (803.7) [makeblastdb]` with the member's path, the line
  or byte and the part of the string read), and a UTF8String tag for such a string gets NCBI's warning once per run. These
  files were refused before.
