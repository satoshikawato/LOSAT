- Engine: the parameters of each query's composition-based redo in a BLASTP call and the
  redo's three parallel paths (the windows, the matches of one query and the queries on the
  thread pool) move out of `run_call` into their own module. No output changes.
