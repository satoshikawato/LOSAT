- The variables `LOSAT_X_DPWINDOW` and `LOSAT_X_DPTRACEPOOL` are gone, and the X-drop DP keeps its memory by one rule
  without a band-start constant: where its rows would grow, it moves them down with the band when half of them or
  more lie below the band start. The shared pool of long traceback rows now follows `LOSAT_X_DPNIBBLE`, so
  `LOSAT_X_DPNIBBLE=0` gives one byte per cell and every worker its own trace, and the pool keeps one idle trace
  instead of two. Outputs do not change: <<NUMBERS>>.
