- LOSAT Web has a **Format Archive** tab (`blast_formatter`): it writes the reports of a search again
  from its `-outfmt 11` archive (a file, pasted text, or a search of the session through
  **Format archive…** beside **Export** on its Outputs) in `-outfmt` 8, 9 and 12, the same bytes as the
  command line, also for the archive of a `-db` search with the database's files under
  **Run directory**. The `losat` npm package runs `blast_formatter` with `engine.runTool`
  (ABI v2 `losat_web2_tool` and `losat_web2_validate`).
- Tests: `blast_formatter`'s outputs are compared byte for byte with the command line's in the V-ABI
  cases, the npm package, LOSAT Web's e2e (every `-outfmt` the tab offers) and `pylosat.blast_formatter`.
