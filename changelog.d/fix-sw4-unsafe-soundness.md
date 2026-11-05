- Internal: `BlastAaLookupTable::get_hits` (tblastx lookup) no longer reads with `get_unchecked`. It is a safe public
  function of a table with `pub` fields, so an index past the table could read outside it. It now returns the cell's
  hits through the bounds-checked `XCell::hits`; the hits it returns, the reports and the instruction count of the
  default path are unchanged (SWEBOK 2026-11-04 audit, KA4). The megablast table fill's arena write states where its
  bound comes from (`Fill::read`'s contract and its five callers).
