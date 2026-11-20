- WASI command builds (`wasm32-wasip1`, `wasm32-wasip1-threads`): `makeblastdb -logfile FILE` now writes
  the log and the messages to the file as natively. Before, every such run stopped with
  `Command line argument error: Argument "logfile". File is not accessible` (exit status 1): makeblastdb
  wrote its log through a second handle of the file that its messages go to, which WASI cannot make. It now
  writes both through the one open file, in the order they are written; native builds give the same bytes
  as before. The WASI command checks (`LOSAT/tests/check_wasm_threading.py`) run every makeblastdb
  `-logfile` row of the drop-in corpus (`dropin_corpus.py check`) and two runs appending a log and a
  warning to one file through each WASI command.
