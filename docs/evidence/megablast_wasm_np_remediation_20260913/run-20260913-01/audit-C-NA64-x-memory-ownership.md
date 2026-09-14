# Independent TBLASTX ownership boundary review

The read-only ncbi_parity_auditor found no new long-lived owner introduced by
C-X2. C-NA64 linking.rs is byte-identical to the previously audited integrated
Rust source (`9bc8d672...`). The limited claim is normal-return owner destruction,
not actual allocator/epoch reclamation time or release of linear-memory pages.

At linking.rs:807 the reference Vec borrows local results. After stable sorts at
808–809 it moves into replay_ncbi_output_list at811. Replay (839–893) returns owned
HSP clones in the existing order; the reference Vec, id map and index vector drop
at replay return, and original results drop at caller return. Empty/singleton
results return directly at804–805. The new Vec never escapes into TLS, global
storage or LAST_RESULT. NCBI correspondence is link_hsps.c:990–994 pointer sorts,
1080–1085 HSP restoration, and1087–1088 helper-array deallocation.

Existing BufferPools own lh_helpers and hsp_links. Parallel TLS POOLS reuse their
capacity between frame groups (683–710); sequential branches have local pools
(733–757 /764–788). DualMaximumTree remains group-local at1243. Search-scoped
worker state and build_scoped pools remain unchanged; Rayon scoped workers
terminate at search completion. The added reference Vec enters none of those
retained owners.

During replay, original results, references and new output coexist. Removing4N
payload copies therefore does not establish a lower simultaneous peak. Planned
capacity diagnostics measure the added reference Vec capacity and element width;
they omit sort workspace, the total existing pool and simultaneous output peak.

Actual X RSS nonregression, process/linear budgets, raw/worker reuse checks and
reference capacity measurements are still required. Together with this source
boundary they can support that added owners are not retained after a search and
observed resource use stays within the agreed limits. X allocator live bytes
remain unmeasured/N/A; N/P live-memory numbers must not be applied to X. No new
runtime change or extra allocator measurement was requested by this review.
