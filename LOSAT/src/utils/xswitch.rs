//! Defaults of the performance switches (`LOSAT_X_*`; list and descriptions in
//! docs/agents/diagnostics.md, "Performance switches").
//!
//! Until main 62102f7a every switch of the strict set was off unless its variable was set to any
//! value. The set is now on by default. Each switch keeps its variable, and one variable,
//! `LOSAT_X_DEFAULT`, returns all of them to the previous default:
//!
//! | `LOSAT_X_<NAME>` | `LOSAT_X_DEFAULT` | the switch is |
//! | --- | --- | --- |
//! | `0` | any | off |
//! | any other non-empty value | any | on (for a value switch, the number given or the built-in value) |
//! | not set, or empty | `0` | off, the default before the change |
//! | not set, or empty | not set, empty, or any other value | on |
//!
//! White space around a value is ignored (`set LOSAT_X_DEFAULT=0 && ...` in cmd.exe stores `0 `),
//! so a value of spaces only counts as empty. Only `0` means off: `00`, `false` or `off` are "any
//! other value".
//!
//! `LOSAT_X_DEFAULT=0` alone runs the previous paths. `LOSAT_X_DEFAULT=0` with some variables set
//! gives the run that those variables gave before the change. The shadow and verify switches, the
//! counters and the tuning variables are not part of the set and keep "set to any value".
//!
//! Every switch of the set gives the reference result by construction, so none of these variables
//! changes the output.
//!
//! No NCBI counterpart: these are LOSAT variables that choose between two LOSAT code paths with
//! the same result. NCBI BLAST+ has no equivalent, and this module computes no value NCBI
//! computes.

use std::ffi::OsStr;
use std::ops::RangeInclusive;
use std::sync::OnceLock;

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
/// verified on the arm64 runners with it on. `LOSAT_X_TBNSTATE` (the scan state of TBLASTN's unit
/// path kept by a thread from one subject to the next, `tblastn/search_gapped.rs`) and
/// `LOSAT_X_TBNXLATE` (the translated frames of a TBLASTN subject read from tables,
/// `tblastn/x_subject_translate.rs`) were added after the change as switches that are off unless
/// set and joined the set together once they were verified and timed with it on.
/// `LOSAT_X_TBNOIDNEXT` (the threads take the next chunk of subjects when they are free,
/// `tblastn/x_oid_flow.rs`) and `LOSAT_X_TBNOIDUNITS` (the units of a subject longer than a
/// chunk's budget on the pool, `tblastn/stage_d_pipeline.rs`) joined the same way, together,
/// after them.
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
/// `LOSAT_X_BXSUBJFLOW` (BLASTX's subject loop of a query chunk as a flow of tasks of 64 subjects
/// searched ahead by the pool's other threads, `blastx/x_subject_flow.rs`) joined the set in the
/// pull request that added it (Owner 2026-10-17: speed-ups land default ON), after it was
/// verified and timed with it on and off.
/// `LOSAT_X_QGROUPPAR` (BLASTN's preliminary gapped loop and traceback of a subject searched on one
/// thread run per query on the pool, `blastn/blast_engine/x_qgroup_par.rs`) joined the set in the
/// pull request that added it, after it was verified and timed with it on and off.
/// `LOSAT_X_PAIRPARGATE` (TBLASTX's unit path only for the query batches where it is expected to
/// be faster than the reference frame loop, `tblastx/x_pair_par.rs`) and `LOSAT_X_TBXLUT` (the
/// reference loop's full lookup table filled in place from per-cell counts,
/// `tblastx/lookup/x_lut_direct.rs`) were added after the change as switches that are off unless
/// set and joined the set together once they were verified and timed with them on and off.
pub(crate) const STRICT_SET: &[&str] = &[
    "LOSAT_X_AHEAD",
    "LOSAT_X_BNDBCUT",
    "LOSAT_X_BNDBSTREAM",
    "LOSAT_X_BNFEWSUBJ",
    "LOSAT_X_BNSTORE",
    "LOSAT_X_BNSTREAM",
    "LOSAT_X_BPCLUT",
    "LOSAT_X_BPCOMPOREUSE",
    "LOSAT_X_BPDBCUT",
    "LOSAT_X_BPDEFER",
    "LOSAT_X_BPDIAG",
    "LOSAT_X_BPGAP",
    "LOSAT_X_BPHITBF",
    "LOSAT_X_BPLISTEDLEAN",
    "LOSAT_X_BPLKPAHEAD",
    "LOSAT_X_BPLKPFLAT",
    "LOSAT_X_BPLUT",
    "LOSAT_X_BPPCOMPACT",
    "LOSAT_X_BPPFREE",
    "LOSAT_X_BPQBLOCK",
    "LOSAT_X_BPQBLOCKFMT",
    "LOSAT_X_BPQBLOCKSA",
    "LOSAT_X_BPQBLOCKXML",
    "LOSAT_X_BPQBLOCKXML2",
    "LOSAT_X_BPREDODBSUBJ",
    "LOSAT_X_BPREDOFIT",
    "LOSAT_X_BPREDOLPT",
    "LOSAT_X_BPREDOPAR",
    "LOSAT_X_BPREDOWINEXPAND",
    "LOSAT_X_BPSCANBF",
    "LOSAT_X_BPSEED",
    "LOSAT_X_BPSUBJMETA",
    "LOSAT_X_BPSUBJRUN",
    "LOSAT_X_BPUNGAP",
    "LOSAT_X_BXBATCH",
    "LOSAT_X_BXCHUNK",
    "LOSAT_X_BXDIAG",
    "LOSAT_X_BXFLOW",
    "LOSAT_X_BXHITBF",
    "LOSAT_X_BXLAZYCTX",
    "LOSAT_X_BXLEAN",
    "LOSAT_X_BXLUT",
    "LOSAT_X_BXPAR",
    "LOSAT_X_BXPOOL",
    "LOSAT_X_BXSCAN",
    "LOSAT_X_BXSCANBF",
    "LOSAT_X_BXSEGMEMO",
    "LOSAT_X_BXSUBJFLOW",
    "LOSAT_X_CODONFAST",
    "LOSAT_X_COMPFAST",
    "LOSAT_X_CTXFAST",
    "LOSAT_X_CTXO1",
    "LOSAT_X_DBIDXREAD",
    "LOSAT_X_DCCSR",
    "LOSAT_X_DEKKER",
    "LOSAT_X_DIAGLV",
    "LOSAT_X_DISCLUTPAR",
    "LOSAT_X_DPAVX2",
    "LOSAT_X_DPEMASK",
    "LOSAT_X_DPFAST",
    "LOSAT_X_DPNEON",
    "LOSAT_X_DPW128",
    "LOSAT_X_DUSTFAST",
    "LOSAT_X_DUSTRING",
    "LOSAT_X_ENVCACHE",
    "LOSAT_X_ERFABSORB",
    "LOSAT_X_ERFMEMO",
    "LOSAT_X_GREEDYFAST",
    "LOSAT_X_GREEDYFUSE",
    "LOSAT_X_GREEDYPRELIM",
    "LOSAT_X_GREEDYSPEC",
    "LOSAT_X_HSPINPLACE",
    "LOSAT_X_IDEALMEMO",
    "LOSAT_X_ITREEFAST",
    "LOSAT_X_KARLINFAST",
    "LOSAT_X_LINKSLICES",
    "LOSAT_X_LINKSORT",
    "LOSAT_X_LUTARENA",
    "LOSAT_X_LUTSPLIT",
    "LOSAT_X_MBBATCH",
    "LOSAT_X_MBDELAY",
    "LOSAT_X_MBLUTKEEP",
    "LOSAT_X_MBLUTONCE",
    "LOSAT_X_MBLUTPREFAULT",
    "LOSAT_X_NEWTONEXACT",
    "LOSAT_X_NEWTONLANES",
    "LOSAT_X_NEWTONLANES8",
    "LOSAT_X_PAIRLUT",
    "LOSAT_X_PAIRPAR",
    "LOSAT_X_PAIRPARDISC",
    "LOSAT_X_PAIRPARDISCARR",
    "LOSAT_X_PAIRPARFOLD",
    "LOSAT_X_PAIRPARGATE",
    "LOSAT_X_PAIRPARNA",
    "LOSAT_X_PAIRWISEPAR",
    "LOSAT_X_PURGEFAST",
    "LOSAT_X_QBATCHPAR",
    "LOSAT_X_QDUSTPAR",
    "LOSAT_X_QGROUPPAR",
    "LOSAT_X_QPREPPAR",
    "LOSAT_X_QSTREAM",
    "LOSAT_X_RPTSTREAM",
    "LOSAT_X_SEEDBUCKET",
    "LOSAT_X_SEEDCOMPACT",
    "LOSAT_X_SEEDRELEASE",
    "LOSAT_X_SEEDUNMAP",
    "LOSAT_X_SEGFAST",
    "LOSAT_X_SEGMEMO",
    "LOSAT_X_SEGSHARE",
    "LOSAT_X_SEGTRIMMEMO",
    "LOSAT_X_TBNBATCH",
    "LOSAT_X_TBNBATCHRULE",
    "LOSAT_X_TBNEVENTS",
    "LOSAT_X_TBNFRAMEREUSE",
    "LOSAT_X_TBNHITBF",
    "LOSAT_X_TBNKEEPFEW",
    "LOSAT_X_TBNLINKWIN",
    "LOSAT_X_TBNOIDCHUNK",
    "LOSAT_X_TBNOIDNEXT",
    "LOSAT_X_TBNOIDUNITS",
    "LOSAT_X_TBNOWNPOOLS",
    "LOSAT_X_TBNPAR",
    "LOSAT_X_TBNQSIDE",
    "LOSAT_X_TBNSCANBF",
    "LOSAT_X_TBNSCANBYTE",
    "LOSAT_X_TBNSSIDE",
    "LOSAT_X_TBNSTATE",
    "LOSAT_X_TBNSTORE",
    "LOSAT_X_TBNXLATE",
    "LOSAT_X_TBNXLATESIMD",
    "LOSAT_X_TBXBAND",
    "LOSAT_X_TBXEXTBATCH",
    "LOSAT_X_TBXFRAMES",
    "LOSAT_X_TBXHITBF",
    "LOSAT_X_TBXKMERGE",
    "LOSAT_X_TBXLUT",
    "LOSAT_X_TBXNEXT",
    "LOSAT_X_TBXPIPE",
    "LOSAT_X_TBXPREFETCH",
    "LOSAT_X_TBXQFRAMES",
    "LOSAT_X_TBXQSLICE",
    "LOSAT_X_TBXRPTPAR",
    "LOSAT_X_TBXSCANBF",
    "LOSAT_X_TBXUNGAP",
    "LOSAT_X_TBXUNITMASK",
];

/// No NCBI counterpart: the table in the module header, for one variable's value. A value that
/// is not UTF-8 is "any other value".
fn state(value: Option<&OsStr>, default_on: bool) -> bool {
    match value.map(|v| v.to_str().map(str::trim)) {
        Some(Some("0")) => false,
        Some(Some("")) | None => default_on,
        Some(_) => true,
    }
}

/// No NCBI counterpart: reads `LOSAT_X_DEFAULT` once; `0` makes the switches that are not set
/// off.
fn default_on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| state(std::env::var_os("LOSAT_X_DEFAULT").as_deref(), true))
}

/// No NCBI counterpart: whether the strict-set switch `name` is on. Reads the environment on
/// every call, as the direct reads it replaces did; callers keep the result where it matters.
pub(crate) fn on(name: &str) -> bool {
    debug_assert!(
        STRICT_SET.contains(&name),
        "{name} is not a strict-set switch"
    );
    state(std::env::var_os(name).as_deref(), default_on())
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
    if !state(raw.as_deref(), default_on()) {
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
    fn a_switch_follows_its_variable_then_the_default() {
        for default_on in [true, false] {
            assert!(!state(os("0"), default_on));
            assert!(!state(os("0 "), default_on));
            assert!(!state(os(" 0\t"), default_on));
            assert!(state(os("1"), default_on));
            assert!(state(os("16"), default_on));
            assert!(state(os("on"), default_on));
            assert!(state(os("00"), default_on));
            assert!(state(os("false"), default_on));
            assert_eq!(state(os(""), default_on), default_on);
            assert_eq!(state(os("  "), default_on), default_on);
            assert_eq!(state(None, default_on), default_on);
        }
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
