# Independent native loop-alignment diagnostic audit

Reviewer: /root/audit_megablast (ncbi_parity_auditor), read-only.
All6 actual outputs (14,590,527 bytes) match their saved official oracles.
Inputs/argv/fixed order/non-overlap/exit0/timed=false and both589-file source
bindings pass. Bin/lib fingerprints retain the same Rust/features/profile,
adding only -C llvm-args=-align-loops=64 and
-C llvm-args=-max-bytes-for-alignment=63. Artifact SHA256:
a1ab28454b880fe913f88329056e852405c3ade8ba85e576be4a5a256b7624db.

Single diagnostic times, not adoption medians:
Mje B1/integrated/aligned16.730531/18.352228/16.679154s;
APself30.284110/32.364619/29.505726s.

Actual binary disassembly matches the saved assembly. The link function changes
from5639 instructions (5572 non-NOP+67NOP),33313 bytes, to5751 instructions
(5572 non-NOP+179NOP),34993 bytes. The observed loop head maps from0x1c7dd0
(mod64=16) to0x1e1c40 (mod64=0).
Independent normalization confirms5572 identical non-NOP instructions and
control flow. It preserves immediates/registers/ordinary memory operands and
external symbols+offsets, maps internal targets to instruction ordinals after
padding, and normalizes RIP relocations. It does not blindly delete numbers.

This supports the same non-padding operations/control flow in this function
with different placement, accompanied by shorter times in these single
searches. It does not prove a unique-loop cause, whole-binary instruction
equivalence or a general speedup rate. The flags affect multiple loops/functions.
Reproducing the flags/bytes with isolated ordinary Cargo configuration and then
qualifying the fresh native matrix is appropriate. Old native FAIL remains;
Wasm evidence cannot be reused before actual new-bundle identity is checked.
