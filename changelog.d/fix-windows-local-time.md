- Windows: makeblastdb's log time (`Building a new DB, current time:`) and the database's date
  (blastdbcmd `-info`'s `Date:`, a report's `Posted date:`) are the local time, as NCBI BLAST+
  writes them on Windows, where they were UTC: LOSAT reads `TZ` in the form of Microsoft's C
  runtime (`tzn[+|-]hh[:mm[:ss]][dzn]`, with the United States daylight-time rules when a
  daylight name follows) or, when `TZ` is unset or empty, the Windows time zone with its
  daylight-time rules, as the C runtime's `localtime` does. Other platforms give the same bytes
  as before.
