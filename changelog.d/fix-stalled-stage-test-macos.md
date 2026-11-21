- BLASTN with two or more threads: an idle helper thread of the partitioned seed stage now reads the
  stage's progress on every idle turn after its first 64 spinning ones and leaves a stage that has not
  moved for 10 ms. It read it on every 64th turn, and on a loaded machine a turn that yields can give
  up the CPU for a scheduler quantum (up to 10 ms on macOS), so a stall of many times 10 ms could go
  unseen and a stage that waits for such a helper to leave waited that long. Reports do not change.
  Test-only: the test that a stalled stage lets the helpers leave holds the stage from the time a helper
  is in it until the helpers have left, so it no longer depends on when the OS runs them (it failed on
  the 3-CPU macOS runner of the platform checks).
