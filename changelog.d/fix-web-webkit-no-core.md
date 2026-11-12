- Tests: the Web E2E tests start WebKit on Linux with a core file limit of 1 byte
  (`web/app/tests/e2e/support/webkit-launcher.sh`), so a crashing WebKit process is reported at once
  ("Target crashed") instead of after the time the kernel spends piping its core to systemd-coredump on
  CI's runners (20 to 60 s in the crashes caught) or to wsl-capture-crash on WSL. The WebKit stalls caught
  on CI (a source stuck "indexing", a page that stopped answering) were such crashes; the other stall forms
  of browser-engine runs (a run stuck "Preparing", a describe stuck "loading") are probably the same. No app change.
