# Independent capacity evidence design review

Read-only `ncbi_parity_auditor` inspected accepted C-N4 and C-P1 owner lifetimes. Linear memory and RSS are not direct scratch-capacity evidence.

C-P1 resets trace_rows_used and clears rows while retaining capacity. Record every retained row's capacity, not only used rows. A diagnostic-only GapAlignScratch Drop observation captures each owner's retained capacity maximum because its rows/outer vectors do not shrink. Record trace-row payload sum/max, outer rows, offsets, DP and edit-op capacity separately. Per-owner maxima/sums must not be called simultaneous process peak. The new reserve call should record actual capacity before/after, band length, growth event count and byte increments; these are not allocator call counts. An identical repeated input cycle should demonstrate stable retained capacity.

BLASTP scratch owners live within subject/query/serial-search scopes and do not persist across searches. C-N4 arena/handle/output vectors are call-local; actual capacity and element width should be observed before materialization. collect() capacity is not assumed equal to input length. Deleted owned scripts are dropped immediately. Sort-internal allocation and whole-process peak are separate measures. These observations motivate the isolated capacity diagnostic; they are not measured capacity results.
