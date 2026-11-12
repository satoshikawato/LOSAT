- Library: the remaining engine modules (`algorithm`, `blastinput`, `cli`, `config`, `report`, `utils`, `seqdb` and
  `api::local_blast`) are private to the crate (`pub(crate)`) instead of public and hidden from the docs, so another
  crate can call only `LOSAT::api` (LOSAT Web's adapter keeps its hidden `web_adapter` module; line LIB, P-V5). The
  tests that named those modules now go through `api`, or live in the crate next to what they test. The code that no
  build used once the modules were private is gone (about 7,100 lines), with 85 tests that tested only that code. The
  command line's output and LOSAT Web behave as before.
