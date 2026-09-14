# Independent observed N/P retained-memory review

The read-only ncbi_parity_auditor supports the limited claim that the observed
post-call N/P retained allocation maxima and final values did not regress.

The auditor recalculated all 96 jobs / 16 series from actual allocator process
JSON. All eight matched pairs have equal maximum and final usable bytes and
block counts, reaching the observed maximum by call2. The boundary is after
runPair input deallocation and worker exit, while retaining the result. This is
not peak live memory during a search, allocated linear memory, or process RSS.

The actual fixed-window process JSON confirms 256 jobs and six of sixteen
last-eight plateau failures. Seven of eight candidate maxima are smaller than
baseline; P132 session0 is +655,360 bytes. Those original failures and the positive
delta remain visible. Actual original reuse records and usage.txt confirm all
four RSS pairs satisfy max(10% of baseline maximum, 16 MiB) and the 4 GiB budget.

All four hashes in runtime-memory-np-nonregression-evidence.json match. The
follow-up document correctly keeps finite observation limits and unproven cause
attribution explicit. Capacity remains pending; no unmeasured X retained-memory
or unlimited-stability claim is supported. Applying the existing Wasm evidence
to C-NA64 additionally requires the planned actual artifact byte identity check.
