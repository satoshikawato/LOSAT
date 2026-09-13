# P1 local candidate — REJECTED

Extracted only the NCBI small-gap predecessor scan (`link_hsps.c:703–746`) into a Wasm non-inlined helper, preserving predecessor order, strict score ties, boundaries, and f64 values. Native inlining was retained. The ordered search, chains and HSP set were untouched.

Fresh `p1-local-exploration-v2` diagnostics and all 3 exploration repetitions were raw-byte identical to the current NCBI oracle at n1/n8, with valid requested worker counts. n1 median: baseline 6.2401 s, candidate 6.2330 s. n8 median: baseline 6.5639 s, candidate 6.4468 s (1.78% improvement). This falls short of the predeclared 10% adoption threshold. No candidate code is integrated into the working tree. Broader adoption tests are not warranted for a rejected candidate.

The first experiment used an incomplete isolated source copy and is explicitly invalid (`../p1-local-exploration/INVALID_BUILD.md`). Its artifact and failure records remain available; none of its timings support this decision.

The next stage uses the original baseline artifacts and normal Node settings. TurboFan-only settings remain a separately measured P1 condition.
