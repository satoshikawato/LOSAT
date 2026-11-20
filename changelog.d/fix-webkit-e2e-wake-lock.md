- Web tests: in WebKit, the pages' Screen Wake Lock is a stand-in, and a lint rule keeps the e2e tests on the
  `test` that installs it. WebKit's GLib ports connect to the session bus for each lock, and GLib's setuid check
  on the way aborted the web process ("getauxval () failed: Interrupted system call") when WebKit's thread
  suspend signal handler, which leaves errno = EINTR, ran inside it; JavaScriptCore suspends the main thread
  thousands of times while it compiles WebAssembly. The page then closed in the middle of a test, more often on
  a busy machine (3 of 8 local runs of the WebKit tests with the engine on an evening when other jobs shared the
  machine, 0 of 7 later at a lower load). Chromium and Firefox keep the real wake lock, and the application is
  unchanged.
