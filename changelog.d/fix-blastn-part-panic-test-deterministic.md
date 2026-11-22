- Test-only: the blastn tests that a panicking participant ends the partitioned seed stage
  (`a_panicking_participant_ends_the_stage` and its stage in pieces) make the participant whose fold
  panics take part instead of counting on a 1 ms sleep in every fold, so they no longer fail on a
  loaded machine (a helper could leave the stalled stage before its first fold and the owner then
  folded the whole stage: 28 of 1,000 runs on two loaded CPUs). No program output changes.
