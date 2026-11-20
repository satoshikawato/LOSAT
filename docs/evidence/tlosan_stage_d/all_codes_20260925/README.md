# Inputs of the all-genetic-codes check (2026-09-25)

`fixtures/code<N>.faa` and `fixtures/code<N>.fna`: written by `generate_fixtures.py` from the Stage A query and subject
(`docs/evidence/tlosan_stage_a/fixtures/query.faa` and `subject_code1.fna`), with five spaced W sites changed to the amino
acid and codon where NCBI's `gc.prt` (pinned commit in `generation.manifest.txt`) differs from code 1 for code N;
`fixtures.tsv` lists the changes. Added in commit d05503f6b (2026-09-25).
