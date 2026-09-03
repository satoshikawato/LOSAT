# Accepted parity exceptions

## Local-subject TBLASTX genetic code

For TBLASTX local `-s/--subject` searches, LOSAT honors `--db-gencode` during subject translation, search, and reporting for every non-default genetic code.

NCBI BLAST+ local `-subject` behavior may differ for this subject-genetic-code case. Do not count a difference caused only by that behavior as a LOSAT parity defect.

The exception does not permit differences in timing, ordering, scoring, filtering, statistics, pruning, formatting, query genetic code, or default-code searches. Record the exact option and affected field whenever this exception is applied.

No other exception is accepted unless `AGENTS.md` is updated with its scope and evidence.
