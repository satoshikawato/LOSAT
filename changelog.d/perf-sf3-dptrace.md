- The variables `LOSAT_X_DPWINDOW` and `LOSAT_X_DPTRACEPOOL` are gone, and the X-drop DP keeps its memory by one rule
  without a band-start constant: where its rows would grow, it moves them down with the band when half of them or
  more lie below the band start. The shared pool of long traceback rows now follows `LOSAT_X_DPNIBBLE`, so
  `LOSAT_X_DPNIBBLE=0` gives one byte per cell and every worker its own trace. Neither time nor memory changes:
  Sakai x EDL933 `-subject` at 16 threads, blastn and dc-megablast x0.993 and x0.995 the time and x1.015 and x1.005 the max RSS (rounds 0.93-1.12, as dev's own spread), at 4 threads x0.998 / x0.997 and x1.007 / x1.008; one-thread instructions x0.994-0.999; outputs equal to dev's and NCBI's.
