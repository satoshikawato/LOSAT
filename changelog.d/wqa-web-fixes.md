- LOSAT Web: the five mode tabs fit a 320 px window, so the search page no longer scrolls sideways and every mode name
  stays inside its tab (#1329).
- LOSAT Web: error messages are easier to read on the grey page and in a selected row; the error colour is a little
  darker and now meets 4.5:1 on every ground it is drawn on (#1330).
- LOSAT Web: an input that holds no sequence (an empty file, whitespace only, or a defline only) no longer gets the
  check mark and "Accepted"; its card says that the input adds nothing to the search (#1337).
- LOSAT Web: Back on the results returns to the Search tab with the form as it was, and Forward shows the results
  again, instead of leaving LOSAT and losing every run and input (#1339).
- LOSAT Web: when the browser does not let the page store result files, the note under the form says in LOSAT's words
  that results are kept in this tab's memory (at most 512 MB) and why, instead of printing the browser's error text;
  in Firefox with the site's data blocked, the page now starts (#1340).
- LOSAT Web: a number and its unit, such as "512 MB", "420 ms" or "40 nt", stay on one line (#1341).
- LOSAT Web: when a part of LOSAT cannot be loaded or stops, the page says so once, with what to do. It asks again by
  itself once the part can be loaded, keeping the pasted inputs, and no run or button hangs; Cancel also works while a
  run's results are being organized (#1343).
- LOSAT Web: a search runs to its results in Firefox with the site's data blocked, instead of failing with "The
  operation is insecure" (#1346).
