- WASI command builds (`wasm32-wasip1`, `wasm32-wasip1-threads`): blastx with `-outfmt` 5 (BLAST XML) and
  13 to 16 (XML2 and JSON) now runs as natively. Before, every such run stopped with
  `Command line argument error: Argument "query". File is not accessible` (exit status 1): in these formats
  LOSAT reads the query file once before the search, and it read a second handle on the file, which WASI
  cannot make. It now reads the open file itself and puts it back where it was; native builds give the same
  bytes as before. The WASI command checks (`LOSAT/tests/check_wasm_threading.py`) run every blastx row of
  these formats of the frozen `-subject` and `-db` fixtures through each WASI command (`--outfmts` of
  `check_losat.py` and `db_search_fixtures.py check`).
