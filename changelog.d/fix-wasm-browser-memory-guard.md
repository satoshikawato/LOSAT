- The threaded Wasm artifacts are built with the shared-memory guard that the Node runner applied before it
  compiled them (`build_wasi_artifacts.py`, the threaded release artifact; LOSAT Web already guarded its
  reactor), so a browser host gets guarded bytes too. `node LOSAT/tests/wasi_shared_memory.js guard|verify`
  is the build step; it fails on an instruction encoding it cannot read. A guarded module carries a
  `losat.shared_memory_guard` custom section, and the runner uses such a module as it is and frees its bytes
  once it is compiled. No output change.
