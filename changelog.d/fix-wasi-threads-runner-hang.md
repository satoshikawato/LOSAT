- Tests: the hand-written threaded Wasm commands of `LOSAT/tests/test_wasi_runners.js` now wait for
  their threads with the thread count they compared, not with a second load of it. A thread that
  added the last count between the two loads and notified before the wait began left the command
  waiting until the test killed it at 20 s, in about 1 of 500 to 3,000 runs of a 2-thread command
  (the guard cache index test failed so locally and in nightly run 44322855271). A new test fails
  every time with the old form. The threaded WASI runner was not at fault and is unchanged.
