# Independent source audit: C-X4 dead helper field

Read-only reviewer: /root/audit_megablast (ncbi_parity_auditor).
The reviewer supports removing maxsum1 and its dedicated maintenance as dead
state. NCBI core/unit_tests references are link_hsps.c105 (definition),
576/671/688/874 (writes), and850 (the only selection test), which is disabled:
if(0) if(H2_helper->maxsum1<=H_hsp_sum)break;

Rust uses are the field, sentinel writes, running_max accumulation dedicated
to that field, and a self-dependent max update. There is no read into candidate
selection, comparisons, statistics, output or diagnostics. The private ordinary
Rust type has no external ABI, serialization, raw-memory comparison or
transmute dependency. Removing only the field and these producers while keeping
sum/next_larger/usize domains/zero sentinels/visitation does not affect semantics.
Live DualMaximumTree/best_sum maxima are separate and must remain unchanged.

Actual target-specific layout sizes and new artifact performance/raw parity
still require measurement. The original integrated native FAIL is unchanged.

The reviewer subsequently checked the actual C-X4 patch/source and native
machine code. All589 bound files match; only linking.rs differs, and no C-X3
scanner is present. The saved patch matches the actual difference. Native
helper initialization uses index*5*8 and stride0x28 before, versus index<<5
and stride0x20 after. The n+2 initialization and two sentinel addresses follow
the same40-to32-byte stride change. The actual objdump output was checked
against both saved ASM files. This is artifact-specific native machine-code
evidence; it is not direct size_of output or evidence of the Wasm layout.
