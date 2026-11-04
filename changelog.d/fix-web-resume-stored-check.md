- LOSAT Web: a tab shown again after it was hidden with searches in progress no longer says "The
  stored results are available." when the browser has removed them (Chromium did in an
  off-the-record context after about half an hour): the check reads back a byte of each run's
  stored results and names the runs that lost theirs ("The browser removed the stored results of
  runs 1–3. Run them again."). Opening such a run says "The browser removed the stored results of
  this run. Run it again.", and a run that cannot make the tab's storage again says "The browser
  removed this tab's storage. Reload the page.", without naming cleared site data as the cause.
