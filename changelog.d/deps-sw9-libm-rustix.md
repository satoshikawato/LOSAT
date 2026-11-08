- Dependencies: `libm`, declared but never called (the FP contract check forbids it), is no longer a
  dependency of LOSAT; rustix's exact `=1.1.6` requirement, which no record explained, is now the
  caret `1.1.6` like every other dependency (Cargo.lock still fixes 1.1.6). No engine change.
