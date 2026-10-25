//! Defaults of the performance switches (`LOSAT_X_*`; list and descriptions in
//! docs/agents/diagnostics.md, "Performance switches").
//!
//! Until main 62102f7a every switch of the strict set was off unless its variable was set to any
//! value. The set is now on by default. Each switch keeps its variable, and `0` turns it off:
//!
//! | `LOSAT_X_<NAME>` | the switch is |
//! | --- | --- |
//! | `0` | off |
//! | not set, empty, or any other value | on (for a value switch, the number given or the built-in value) |
//!
//! White space around a value is ignored (`set LOSAT_X_BPSEED=0 && ...` in cmd.exe stores `0 `),
//! so a value of spaces only counts as empty. Only `0` means off: `00`, `false` or `off` are "any
//! other value".
//!
//! `LOSAT_X_DEFAULT`, which turned every switch of the set off at once (the previous paths), was
//! removed on 2026-10-23 (switch-tidy item 3); a run ignores it. The shadow and verify switches,
//! the counters and the tuning variables are not part of the set and keep "set to any value".
//!
//! A switch that graduated has no variable any more, so its previous path is not reachable: an
//! always-on switch lost that code, and a conditional or platform switch (`GRADUATED`) runs its
//! fallback only where its rule sends the work.
//!
//! Every switch of the set gives the reference result by construction, so no variable of the set
//! changes the output.
//!
//! No NCBI counterpart: these are LOSAT variables that choose between two LOSAT code paths with
//! the same result. NCBI BLAST+ has no equivalent, and this module computes no value NCBI
//! computes.

use std::ffi::OsStr;
use std::ops::RangeInclusive;

/// No NCBI counterpart: the switches of the strict set, the variables that are on by default,
/// sorted by name. A switch joins the set with its name in order here and its `export` line in
/// order in the block of `docs/agents/diagnostics.md` ("Performance switches"; a test below
/// compares the two), and its entry there says since when it is in the set and on what evidence.
/// The size of the set is written nowhere, and a joining switch adds nothing at the end of a
/// shared list, so two pull requests that each add a switch change different lines.
///
/// The switches that joined after the set became the default, until this rule (2026-10-17):
/// `LOSAT_X_PAIRLUT` and
/// `LOSAT_X_HSPINPLACE` were added after the change as switches that are off unless set and
/// joined the set once they were verified with it on. `LOSAT_X_SEEDBUCKET` was left out at the
/// change: it numbered the lookup hits of a chunk scan in a `u32`, exact only up to 2^32 hits
/// per scan. It now numbers them per flush (`tblastx/x_seed_bucket.rs`) and joined the set once
/// it was timed on the searches it acts on with the set on (the TBLASTX searches that
/// `LOSAT_X_PAIRPAR` leaves to the scan loop). `LOSAT_X_TBNOIDCHUNK` (one entry into the TBLASTN
/// search pool takes four chunks of subjects per thread, cut by residues as NCBI's database
/// reader cuts them for its search threads, `tblastn/stage_d_pipeline.rs`),
/// `LOSAT_X_SEEDRELEASE` (when TBLASTX frees the tables of its seed stage,
/// `tblastx/x_seed_release.rs`) were added after the change as switches that are off unless set
/// and joined the set once they were verified with it on. `LOSAT_X_GREEDYFUSE` (the fused loop
/// of megablast's score-only greedy alignment, `blastn/alignment/x_greedy_fused.rs`) and
/// `LOSAT_X_LINKSORT` (the index sorts of the TBLASTX linking call, `tblastx/x_link_sort.rs`) were
/// added after the change as switches that are off unless set and joined the set once they were
/// verified with it on. `LOSAT_X_NEWTONLANES8` (eight Newton problems per round in the lanes of
/// `LOSAT_X_NEWTONLANES`, `composition_adjustment/x_newton_lanes.rs`; it acts on aarch64 only)
/// was added after the change as a switch that is off unless set and joined the set once it was
/// verified on the arm64 runners with it on. `LOSAT_X_TBNSTATE` and `LOSAT_X_TBNXLATE` (TBLASTN's
/// kept scan state and table translation) joined the set together the same way. Both graduated
/// on 2026-10-22: they left `STRICT_SET` and their old code was deleted.
/// `LOSAT_X_TBNOIDNEXT` (the threads take the next chunk of subjects when they are free,
/// `tblastn/x_oid_flow.rs`) and `LOSAT_X_TBNOIDUNITS` (the units of a subject longer than a
/// chunk's budget on the pool, `tblastn/stage_d_pipeline.rs`) joined the same way, together,
/// after them. `LOSAT_X_TBNOIDCHUNK`,
/// `LOSAT_X_TBNOIDNEXT`, `LOSAT_X_TBNOIDUNITS`, `LOSAT_X_TBNSSIDE`, `LOSAT_X_TBNBATCH`,
/// `LOSAT_X_SEGSHARE` and `LOSAT_X_DEKKER` graduated on 2026-10-22 (see `GRADUATED`).
/// `LOSAT_X_ERFABSORB` (Spouge's ErfC returns 2.0 below -6 without evaluating ErfImpl,
/// `stats/spouge.rs`) and `LOSAT_X_SEGTRIMMEMO` (the SEG trim of a segment kept per thread for the
/// maskers that repeat segments, `utils/x_seg_trim_memo.rs`) graduated on 2026-10-26: they left
/// `STRICT_SET` and their previous code was deleted (the maskers that do not ask for the memo keep
/// trimming every segment, a rule of the masker and not a switch).
/// `LOSAT_X_LINKSLICES` (the groups of the TBLASTX linking call linked as slices of its one list,
/// `tblastx/x_link_slices.rs`) was added after the change as a switch that is off unless set and
/// joined the set once it was verified with it on.
/// `LOSAT_X_BPSUBJRUN` (BLASTP's parallel `-subject` loop searches the subjects in runs of
/// consecutive subjects, each run with one diagonal table that is cleared once,
/// `blastp/blast_engine/x_subject_runs.rs`) was added after the change as a switch that is off
/// unless set and joined the set once it was verified and timed with it on.
/// `LOSAT_X_TBNSTORE` (TBLASTN `-subject` records held packed, `tblastn/x_subject_store.rs`) was
/// added after the change as a switch that is off unless set and joined the set once it was
/// verified and timed with it on.
/// `LOSAT_X_DPNEON` (the NEON build of the vector X-drop kernels on aarch64, `utils/xdrop_simd.rs`)
/// was added after the change as a switch that is off unless set and joined the set once it was
/// verified with it on the Linux and macOS arm64 runners.
/// `LOSAT_X_DPNEON2` (the X-drop kernels of `LOSAT_X_DPAVX2` on aarch64 NEON in register pairs,
/// `utils/x_xdrop_neon2.rs`) was added on 2026-10-15 as a switch that is off unless set and joined
/// the set on 2026-10-24 (Owner 2026-10-23), after the arm64 platform checks hashed it in every
/// run since it was added.
/// `LOSAT_X_TBNBATCHRULE` (which query batches `LOSAT_X_TBNBATCH` searches side by side, only the
/// full waves of a small subject set, `tblastn/stage_d_pipeline.rs`) was added after the change as
/// a switch that is off unless set and joined the set once it was verified and timed with it on.
/// The eight blastp switches of the OrthoFinder work (`LOSAT_X_BPSEED`,
/// `LOSAT_X_BPUNGAP`, `LOSAT_X_BPGAP`, `LOSAT_X_BPLUT`, `LOSAT_X_BPPFREE`, `LOSAT_X_BPPCOMPACT`,
/// `LOSAT_X_BPREDOFIT` and `LOSAT_X_BPQBLOCK`) joined after a 49-job full-scale run with all of
/// them on: the outputs equal NCBI's and the max RSS is below NCBI's on every job.
/// `LOSAT_X_BNSTREAM` (BLASTN hands each subject's HSPs to the hit lists in OID order from a
/// bounded window while the threads search, `blastn/blast_engine/x_bnstream.rs`) was added after
/// the change as a switch that is off unless set and joined the set once it was verified and timed
/// with it on.
/// `LOSAT_X_TBXNEXT` (the threads of TBLASTX's parallel subject loop take the next subject when
/// they are free, `tblastx/x_subject_next.rs`) and `LOSAT_X_TBXFRAMES` (a TBLASTX subject's frames
/// made when they are read and freed when they are no longer read, `tblastx/subjects.rs`) were
/// added after the change as switches that are off unless set and joined the set together once
/// they were verified and timed with it on.
/// `LOSAT_X_TBXFRAMES` and `LOSAT_X_TBXUNGAP` (the byte table of TBLASTX's two-hit ungapped
/// extension, `tblastx/extension/x_two_hit_fast.rs`) graduated on 2026-10-26: they left
/// `STRICT_SET` and their previous code was deleted (the calls the byte table does not cover
/// still run `extend_hit_two_hit_reference`, a rule of the table walks and not a switch).
/// `LOSAT_X_TBNHITBF` (the branch-free two-hit loop of TBLASTN's unit path,
/// `tblastn/x_tbn_hitbf.rs`), `LOSAT_X_TBNLINKWIN` (the linking candidates of a head chosen by
/// subject position, `tblastn/stage_d_linking/x_link_window.rs`) and `LOSAT_X_TBNSCANBF` (the
/// lookup scan without a branch on the presence vector, `tblastn/x_scan_bf.rs`) graduated on
/// 2026-10-26: they left `STRICT_SET` and their previous code was deleted (NCBI's whole-window
/// loop of the linking call stays as the reference of a unit test, and `XTwoHit::pairs` stays as
/// the body that the branch-free loop hands the other hits to).
/// `LOSAT_X_BPDIAG` (the two-hit diagonal test of a lookup cell eight hits at a time,
/// `blastp/x_seed_fused.rs`) joined after a second run of the same jobs with it on as well.
/// `LOSAT_X_SEEDCOMPACT` (TBLASTX's unit path keeps the saved init HSPs of a unit as 28-byte
/// records copied out at their exact length, `tblastx/x_pair_par.rs`) was added after the change as
/// a switch that is off unless set and joined the set once it was verified and timed with it on.
/// `LOSAT_X_BPCLUT` (the compressed lookup table of `blastp-fast` built flat, in three passes, by
/// the search's pool at more than one thread, `tblastx/lookup/x_clut_direct.rs`) and
/// `LOSAT_X_BPLKPFLAT` (how that build is scheduled and where its vectors come from, which acts
/// only where `LOSAT_X_BPCLUT` builds the table on the pool, `tblastx/lookup/x_clut_par.rs`) were
/// added after the change as switches that are off unless set and joined the set together once
/// they were verified and timed with them on.
/// `LOSAT_X_BPDEFER` (the candidates of a lookup cell's AVX2 batches called after its last batch
/// of eight, `blastp/x_seed_fused.rs`) joined after a third run of the same jobs with it on as well.
/// `LOSAT_X_TBNOWNPOOLS` (the query batches `LOSAT_X_TBNBATCHRULE` leaves one after another
/// searched at the same time, each on a pool of its own, `tblastn/stage_d_pipeline.rs`) was added
/// after the change as a switch that is off unless set and joined the set once it was verified
/// and timed with it on.
/// `LOSAT_X_TBNFLOWHELP` (a thread that waits in TBLASTN's subject flow runs the pool's jobs
/// meanwhile, `tblastn/x_oid_flow.rs`) joined the set in the pull request that added it, after it
/// was verified and timed with it on and off.
/// `LOSAT_X_TBNUNITXLATE` (a long packed TBLASTN subject whose plan is made for its search alone
/// gets a plan without residues, and each unit translates the residues it reads,
/// `tblastn/x_subject_translate.rs`) joined the set in the pull request that added it, after it
/// was verified and timed with it on and off.
/// `LOSAT_X_BXSUBJFLOW` (BLASTX's subject loop of a query chunk as a flow of tasks of 64 subjects
/// searched ahead by the pool's other threads, `blastx/x_subject_flow.rs`) joined the set in the
/// pull request that added it (Owner 2026-10-17: speed-ups land default ON), after it was
/// verified and timed with it on and off.
/// `LOSAT_X_BNTABPAR` (BLASTN's outfmt 6 / 7 rows of a report formatted per query on the search's
/// pool while no LOSAT_X_QBATCHPAR slot exists, `blastn/x_tab_par.rs`) joined the set in the pull
/// request that added it, after it was verified and timed with it on and off.
/// `LOSAT_X_DCTWOSLOTS` (BLASTN's dc-megablast takes two LOSAT_X_QBATCHPAR slots of half the
/// threads where the slot model takes three or more and the two-slot rule holds,
/// `blastn/blast_engine/x_qbatch_par.rs`) joined the set in the pull request that added it, after
/// it was verified and timed with it on and off.
/// `LOSAT_X_BNEIGHTSLOTS` (BLASTN's `-task blastn` runs the slot model's best layout of at most
/// eight LOSAT_X_QBATCHPAR slots where the model takes more at 12 to 64 threads,
/// `blastn/blast_engine/x_qbatch_par.rs`) joined the set in the pull request that added it, after
/// it was verified and timed with it on and off.
/// `LOSAT_X_BNCSR` (the chains of BLASTN's contiguous lookup tables of word sizes 9 to 11 laid out
/// as CSR and scanned through a prefetch ring, `blastn/x_bn_csr.rs`) joined the set in the pull
/// request that added it, after it was verified and timed with it on and off.
/// `LOSAT_X_QGROUPPAR` (BLASTN's preliminary gapped loop and traceback of a subject searched on one
/// thread run per query on the pool, `blastn/blast_engine/x_qgroup_par.rs`) joined the set in the
/// pull request that added it, after it was verified and timed with it on and off.
/// `LOSAT_X_BXCHUNKFLOW` (BLASTX's chunks of a split query searched one after another, each chunk's
/// subjects as the subject flow on the whole pool, where each chunk has enough tasks for the pool,
/// `blastx/x_subject_flow.rs`) joined the set in the pull request that added it, after it was
/// verified and timed with it on and off.
/// `LOSAT_X_BXDIAGKEEP`, `LOSAT_X_BXLUTPOOL` and `LOSAT_X_BXUNITCAP` (BLASTX's memory at many
/// threads: the workers' diagonal tables kept one per thread, the lookup tables' neighbour lists
/// shared and their vectors kept between builds, at most K query batches or chunks searched at a
/// time, `blastx/x_bx_mem.rs`, `tblastx/lookup/x_lut_direct.rs`) joined the set in the pull
/// request that added them, after they were verified and timed with them on and off.
/// `LOSAT_X_PAIRPARGATE` (TBLASTX's unit path only for the query batches where it is expected to
/// be faster than the reference frame loop, `tblastx/x_pair_par.rs`) and `LOSAT_X_TBXLUT` (the
/// reference loop's full lookup table filled in place from per-cell counts,
/// `tblastx/lookup/x_lut_direct.rs`) were added after the change as switches that are off unless
/// set and joined the set together once they were verified and timed with them on and off.
/// `LOSAT_X_BXSCANBF` (BLASTX's seed scans find the subject words with lookup hits by tblastn's
/// branch-free candidate loop, `blastx/x_seed_scan.rs`) graduated on 2026-10-26: it left
/// `STRICT_SET` and the scan on the presence test was deleted (the compressed lookup table keeps
/// its own scan, a rule of the table and not a switch).
pub(crate) const STRICT_SET: &[&str] = &[
    "LOSAT_X_BNCSR",
    "LOSAT_X_BNDBCUT",
    "LOSAT_X_BNDBSTREAM",
    "LOSAT_X_BNDCONE",
    "LOSAT_X_BNEIGHTSLOTS",
    "LOSAT_X_BNFEWSUBJ",
    "LOSAT_X_BNLUTDROP",
    "LOSAT_X_BNMAPPED",
    "LOSAT_X_BNQDROP",
    "LOSAT_X_BNQPACK",
    "LOSAT_X_BNSTORE",
    "LOSAT_X_BNSTREAM",
    "LOSAT_X_BNTABPAR",
    "LOSAT_X_BNTBRANGE",
    "LOSAT_X_BNUNGAPCOUNT",
    "LOSAT_X_BPCLUT",
    "LOSAT_X_BPCLUTSCANBF",
    "LOSAT_X_BPCOMPOREUSE",
    "LOSAT_X_BPDBCUT",
    "LOSAT_X_BPDEFER",
    "LOSAT_X_BPDIAG",
    "LOSAT_X_BPFFILL",
    "LOSAT_X_BPFSEED",
    "LOSAT_X_BPFSTAGE",
    "LOSAT_X_BPGAP",
    "LOSAT_X_BPGAP16",
    "LOSAT_X_BPHITBF",
    "LOSAT_X_BPLISTEDLEAN",
    "LOSAT_X_BPLKDROP",
    "LOSAT_X_BPLKPAHEAD",
    "LOSAT_X_BPLKPFLAT",
    "LOSAT_X_BPLUT",
    "LOSAT_X_BPMASKED",
    "LOSAT_X_BPPCOMPACT",
    "LOSAT_X_BPPFREE",
    "LOSAT_X_BPQBLOCK",
    "LOSAT_X_BPQBLOCKFMT",
    "LOSAT_X_BPQBLOCKSA",
    "LOSAT_X_BPQBLOCKWEB",
    "LOSAT_X_BPQBLOCKXML",
    "LOSAT_X_BPQBLOCKXML2",
    "LOSAT_X_BPREDODBSUBJ",
    "LOSAT_X_BPREDOFIT",
    "LOSAT_X_BPREDOLPT",
    "LOSAT_X_BPREDOPAR",
    "LOSAT_X_BPREDOWINEXPAND",
    "LOSAT_X_BPSCANBF",
    "LOSAT_X_BPSEED",
    "LOSAT_X_BPSEGWIN",
    "LOSAT_X_BPSUBJMETA",
    "LOSAT_X_BPSUBJRUN",
    "LOSAT_X_BPUNGAP",
    "LOSAT_X_BPUNGAPSENT",
    "LOSAT_X_BXCHUNKFLOW",
    "LOSAT_X_BXCROSS",
    "LOSAT_X_BXCROSSDEFER",
    "LOSAT_X_BXDBSTREAM",
    "LOSAT_X_BXDEFER",
    "LOSAT_X_BXDIAG",
    "LOSAT_X_BXDIAGKEEP",
    "LOSAT_X_BXEVALEMPTY",
    "LOSAT_X_BXFLOW",
    "LOSAT_X_BXHITBF",
    "LOSAT_X_BXKAONCE",
    "LOSAT_X_BXLUTPOOL",
    "LOSAT_X_BXSCANBYTE",
    "LOSAT_X_BXSUBJFLOW",
    "LOSAT_X_BXUNGAPSENT",
    "LOSAT_X_BXUNITCAP",
    "LOSAT_X_CTXO1",
    "LOSAT_X_DBIDXREAD",
    "LOSAT_X_DCCSR",
    "LOSAT_X_DCTWOSLOTS",
    "LOSAT_X_DIAGLV",
    "LOSAT_X_DIAGLV2",
    "LOSAT_X_DISCLUTPAR",
    "LOSAT_X_DISPTABLE",
    "LOSAT_X_DPEMASK",
    "LOSAT_X_DPNEON2",
    "LOSAT_X_DPNIBBLE",
    "LOSAT_X_DPTRACEPOOL",
    "LOSAT_X_DPWINDOW",
    "LOSAT_X_FOLDLEAN",
    "LOSAT_X_GREEDYROWFUSE",
    "LOSAT_X_MBFILL",
    "LOSAT_X_MBLUTKEEP",
    "LOSAT_X_MBLUTONCE",
    "LOSAT_X_MBLUTPREFAULT",
    "LOSAT_X_PAIRPARDISC",
    "LOSAT_X_PAIRPARDISCARR",
    "LOSAT_X_PAIRPARFOLD",
    "LOSAT_X_PAIRPARGATE",
    "LOSAT_X_PAIRPARHOME",
    "LOSAT_X_PAIRPARMERGE",
    "LOSAT_X_PAIRPARN1",
    "LOSAT_X_PAIRPARNA",
    "LOSAT_X_PAIRWISEPAR",
    "LOSAT_X_QBATCHPAR",
    "LOSAT_X_QDUSTPAR",
    "LOSAT_X_QGROUPPAR",
    "LOSAT_X_QPREPPAR",
    "LOSAT_X_QSTREAM",
    "LOSAT_X_REDOCONTAINIDX",
    "LOSAT_X_RPTSTREAM",
    "LOSAT_X_SEEDCOMPACT",
    "LOSAT_X_SEEDUNMAP",
    "LOSAT_X_SUMPMEMO",
    "LOSAT_X_TBNBATCHRULE",
    "LOSAT_X_TBNDBSTREAM",
    "LOSAT_X_TBNFLOWHELP",
    "LOSAT_X_TBNFRAMEREUSE",
    "LOSAT_X_TBNKEEPFEW",
    "LOSAT_X_TBNLINKSETS",
    "LOSAT_X_TBNOWNPOOLS",
    "LOSAT_X_TBNREDOLPT",
    "LOSAT_X_TBNSCANBYTE",
    "LOSAT_X_TBNSTAGEDROP",
    "LOSAT_X_TBNSTORE",
    "LOSAT_X_TBNUNGAP",
    "LOSAT_X_TBNUNITXLATE",
    "LOSAT_X_TBNXLATESIMD",
    "LOSAT_X_TBXBAND",
    "LOSAT_X_TBXEXTBATCH",
    "LOSAT_X_TBXHITBF",
    "LOSAT_X_TBXKMERGE",
    "LOSAT_X_TBXLUT",
    "LOSAT_X_TBXMASKSEEK",
    "LOSAT_X_TBXNEXT",
    "LOSAT_X_TBXPIPE",
    "LOSAT_X_TBXPREFETCH",
    "LOSAT_X_TBXQFRAMES",
    "LOSAT_X_TBXQSLICE",
    "LOSAT_X_TBXRESOLVE",
    "LOSAT_X_TBXRPTPAR",
    "LOSAT_X_TBXSCANBF",
    "LOSAT_X_TBXSTORE",
    "LOSAT_X_TBXUNGAPSENT",
    "LOSAT_X_TBXUNITMASK",
];

/// No NCBI counterpart: the conditional and platform switches that left the strict set when they
/// graduated (on dev at least 7 days and in two full-set platform-checks runs with equal outputs
/// and no failed test, 2026-10-22), sorted by name. Their variables are gone; their rules (input
/// sizes, threads, options, memory, CPU) still choose between their own path and the fallback, so
/// the names stay as the keys of `x_reach!` (`utils/xreach.rs`) and of the rows of
/// `LOSAT/tests/switch_reach_expected.tsv`. An always-on switch that graduates leaves both lists:
/// nothing chooses any more.
pub(crate) const GRADUATED: &[&str] = &[
    "LOSAT_X_AHEAD",
    "LOSAT_X_BXBATCH",
    "LOSAT_X_BXCHUNK",
    "LOSAT_X_BXPAR",
    "LOSAT_X_BXPOOL",
    "LOSAT_X_DEKKER",
    "LOSAT_X_DPAVX2",
    "LOSAT_X_DPNEON",
    "LOSAT_X_DPW128",
    "LOSAT_X_GREEDYSPEC",
    "LOSAT_X_MBBATCH",
    "LOSAT_X_MBDELAY",
    "LOSAT_X_NEWTONLANES",
    "LOSAT_X_NEWTONLANES8",
    "LOSAT_X_PAIRLUT",
    "LOSAT_X_PAIRPAR",
    "LOSAT_X_SEEDBUCKET",
    "LOSAT_X_SEGSHARE",
    "LOSAT_X_TBNBATCH",
    "LOSAT_X_TBNOIDCHUNK",
    "LOSAT_X_TBNOIDNEXT",
    "LOSAT_X_TBNOIDUNITS",
    "LOSAT_X_TBNSSIDE",
];

/// No NCBI counterpart: the table in the module header, for one variable's value. A value that
/// is not UTF-8 is "any other value".
fn state(value: Option<&OsStr>) -> bool {
    !matches!(value.map(|v| v.to_str().map(str::trim)), Some(Some("0")))
}

/// No NCBI counterpart: whether the strict-set switch `name` is on. Reads the environment on
/// every call, as the direct reads it replaces did; callers keep the result where it matters.
pub(crate) fn on(name: &str) -> bool {
    debug_assert!(
        STRICT_SET.contains(&name),
        "{name} is not a strict-set switch"
    );
    state(std::env::var_os(name).as_deref())
}

/// No NCBI counterpart: a strict-set switch that carries a number. `None` when the switch is off;
/// otherwise the number given if it lies in `range`, and `builtin` if the variable is not set or
/// holds anything else.
pub(crate) fn value(name: &str, range: RangeInclusive<usize>, builtin: usize) -> Option<usize> {
    debug_assert!(
        STRICT_SET.contains(&name),
        "{name} is not a strict-set switch"
    );
    let raw = std::env::var_os(name);
    if !state(raw.as_deref()) {
        return None;
    }
    Some(given(raw.as_deref(), &range).unwrap_or(builtin))
}

/// No NCBI counterpart: the number in a value switch's variable, if it is one within `range`.
fn given(value: Option<&OsStr>, range: &RangeInclusive<usize>) -> Option<usize> {
    value?
        .to_str()?
        .trim()
        .parse::<usize>()
        .ok()
        .filter(|n| range.contains(n))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn os(text: &str) -> Option<&OsStr> {
        Some(OsStr::new(text))
    }

    #[test]
    fn a_switch_follows_its_variable() {
        for off in ["0", "0 ", " 0\t"] {
            assert!(!state(os(off)), "{off:?}");
        }
        for on in ["1", "16", "on", "00", "false", "", "  "] {
            assert!(state(os(on)), "{on:?}");
        }
        assert!(state(None));
    }

    #[test]
    fn a_value_switch_takes_a_number_in_range_only() {
        let range = 2..=65_536;
        assert_eq!(given(os("32"), &range), Some(32));
        assert_eq!(given(os("1"), &range), None);
        assert_eq!(given(os("65537"), &range), None);
        assert_eq!(given(os(" 32 "), &range), Some(32));
        assert_eq!(given(os("3 2"), &range), None);
        assert_eq!(given(os("x"), &range), None);
        assert_eq!(given(os(""), &range), None);
        assert_eq!(given(None, &range), None);
    }

    #[test]
    fn the_strict_set_names_are_distinct_and_sorted() {
        assert!(STRICT_SET.windows(2).all(|pair| pair[0] < pair[1]));
        assert!(STRICT_SET.iter().all(|name| name.starts_with("LOSAT_X_")));
    }

    #[test]
    fn the_graduated_names_are_distinct_sorted_and_outside_the_set() {
        assert!(GRADUATED.windows(2).all(|pair| pair[0] < pair[1]));
        assert!(GRADUATED.iter().all(|name| name.starts_with("LOSAT_X_")));
        assert!(GRADUATED.iter().all(|name| !STRICT_SET.contains(name)));
    }

    /// The block of `docs/agents/diagnostics.md` that writes the set as `export` lines names the
    /// switches of `STRICT_SET`, one per line, in the same order.
    #[test]
    fn the_documented_block_lists_the_strict_set() {
        let doc = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/agents/diagnostics.md"
        ));
        let anchor = doc
            .find("The strict set (the variables that are on;")
            .expect("the block's heading in diagnostics.md");
        let block = &doc[anchor..];
        let start = block.find("```bash\n").expect("the block") + "```bash\n".len();
        let end = start + block[start..].find("```").expect("the end of the block");
        let documented: Vec<&str> = block[start..end]
            .lines()
            .map(|line| {
                let name = line
                    .strip_prefix("export ")
                    .and_then(|rest| rest.split_once('='))
                    .map(|(name, _)| name);
                name.unwrap_or_else(|| panic!("not one `export NAME=value` line: {line:?}"))
            })
            .collect();
        assert_eq!(documented, STRICT_SET);
    }
}
