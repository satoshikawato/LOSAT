#!/bin/sh
# Starts the WebKit of the E2E tests on Linux (playwright.config.ts) with a core file limit of 1 byte, so that a
# WebKit process that crashes writes no core. On CI's runners and on WSL the kernel sends cores to a program through
# a pipe (systemd-coredump, wsl-capture-crash); a pipe cannot seek, so the kernel writes every page of the mappings
# the process used, unused ones as zeros (several GB per dump on WSL), and the page answers nothing until the dump
# ends: 20 to 60 s in the crashes caught on CI, where a test saw a source stuck "indexing" or a page that stopped
# answering (line WFL item PH). With a limit of 1 byte the kernel gives up a piped core at once, and Playwright
# reports the crash ("Target crashed") when it happens.
real=${LOSAT_WEB_WEBKIT_REAL:?LOSAT_WEB_WEBKIT_REAL names the WebKit launcher that this one starts}
if [ ! -x "$real" ]; then
  echo "webkit-launcher: $real does not exist or cannot run; install WebKit with: npx playwright install webkit" >&2
  exit 127
fi
# Where the hard limit is 0, no core is written anyway, and a soft limit of 1 cannot be set.
[ "$(ulimit -H -c)" = 0 ] && exec "$real" "$@"
exec prlimit --core=1: -- "$real" "$@"
