- `losat rpsbproc` reads an archive whose text strings hold bytes outside the visible ASCII range as NCBI's rpsbproc
  0.5.1 does, instead of refusing it. In a VisibleString, each such byte is read as `#`, and NCBI's error is posted for
  each one, before the session, with NCBI's path to the value, the line and the piece of the string. This covers every
  part of the archive, including the parts that rpsbproc does not use. A UTF8String, such as a DefLine's title, keeps
  its bytes.
