# Ownership boundary for the reuse investigation

These are source-level ownership observations, not measured allocator live-byte
counts or proof of a plateau over unlimited calls. The frozen integrated source
manifest binds the files; no production change was made for this investigation.

BLASTN `purge_hsps_for_subject_ex` in `purge_endpoints.rs:532` receives owned
input HSPs. The arena and handle vector are function-local owned vectors. The
arena does not grow after construction. Deleted payloads are taken and dropped
at the existing purge points; final live payloads are taken into the returned
vector at lines779–792. The handle vector is consumed and the remaining arena
is dropped on return. The returned payloads remain owned by the caller. There
is no new static, leaked owner or cross-search cache. Relevant NCBI semantics:
`blast_hits.c:2471–2478,2490–2499,2516–2525,2530–2535`.

BLASTP `GapAlignScratch` in `gapalign.rs:345` owns its DP vector, row vectors,
offset vector and reverse-operation vector. Row reuse resets its logical row
count (`gap_reset_traceback_state`, line547), clears a selected row before
writing (`gap_alloc_trace_row`, line577), and retains its capacity within that
scratch owner's lifetime. The candidate's lines1822–1823 borrow reserved spare
capacity; each written value is a byte, and the single `set_len` at line1910
exposes only the initialized prefix. It introduces no additional owner or
destructor-bearing payload. NCBI references are `blast_gapalign.h:69–80` and
`blast_gapalign.c:70–115,513–540,563–639`.

The engine creates Kappa scratch as local owned state for its parallel work
closures (`blast_engine.rs:5549,5651`) or the sequential search path (line5780).
Synchronous postprocessing borrows that state; the candidate does not change
those lifetimes. Search pools retain the existing scoped construction and
joining behavior (`utils/threading.rs:92–124`, NCBI `prelim_stage.cpp:145–188`).

These facts rule out a newly introduced escaping owner in the touched N/P
storage operations. They do not identify the cause of every observed Wasm
high-water increase. Actual retained capacities, process RSS and per-instance
linear-memory observations are reported separately. Additional diagnostics
must not turn the original six-call plateau failures into passing records.
