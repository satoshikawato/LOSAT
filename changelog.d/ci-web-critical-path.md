- CI: web.yml's longest path is shorter. `build_reactors.py` builds the two reactors at once (4.3 min
  instead of 7.6 on 4 CPUs, same bytes); the native CLI that V-ABI and V-BR compare with is built
  without fat LTO (2.7 min instead of 6.7); Chromium's engine tests run in two shards and WebKit's in
  three; and the FakeEngine application's E2E tests run in one job per browser.
