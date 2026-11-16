- Engine: the plan of the composition-based redo in a BLASTP call (which of its paths runs)
  and its one-thread loop move out of `run_call` into the redo's module, where the plan and
  all four paths read one set of the call's values. No output changes.
