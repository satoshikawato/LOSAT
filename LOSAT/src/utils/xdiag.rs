//! The diagnostic build: the trace, debug, dump and profile variables that no test, script or CI
//! step sets are read only when the `diagnostics` feature is on.
//!
//! A default build ignores these variables: `var` and `var_os` return "not set" without reading the
//! environment, so the code behind them (the tracing statics and their printing paths) is dead and
//! the compiler drops it, and `main` writes one notice line to stderr for each of them that is set
//! (`ignored_notices`). Build with `cargo build --release --features diagnostics` to use them (the
//! WASI command adds `--target wasm32-wasip1 --no-default-features`). Both builds compile the same
//! code; there are no `#[cfg(feature = "diagnostics")]` blocks, only the constant `ON`.
//!
//! The rule for a new diagnostic variable: one that no test, script or CI step sets is read through
//! this module and listed in `NAMES` and in the list of `docs/agents/diagnostics.md` ("Diagnostic
//! build"); one that they set is read directly with `std::env`, because a test of a default build
//! would not see it otherwise. The tests below check the list, the docs block, and that no other
//! file reads a listed name directly.
//!
//! No NCBI counterpart: these are LOSAT debugging variables and NCBI BLAST+ has no equivalent;
//! this module only decides whether LOSAT reads them, and it does not change any value NCBI
//! computes.

use std::env::VarError;
use std::ffi::OsString;

/// Whether this build reads the diagnostic variables.
pub const ON: bool = cfg!(feature = "diagnostics");

/// The diagnostic variables, sorted. `docs/agents/diagnostics.md` lists the same names.
pub const NAMES: &[&str] = &[
    "LOSAT_DEBUG_BLASTP_TRACE_HSP",
    "LOSAT_DEBUG_BLASTP_WSSV",
    "LOSAT_DEBUG_CHAINING",
    "LOSAT_DEBUG_CUTOFFS_ALL",
    "LOSAT_DEBUG_EXTENSION",
    "LOSAT_DEBUG_HSP_SAVING",
    "LOSAT_DEBUG_SCAN_SOFF",
    "LOSAT_DEBUG_SCAN_WINDOW",
    "LOSAT_DUMP_TBLASTX_STAGE",
    "LOSAT_LINK_WIDE_KEYS",
    "LOSAT_STARTUP_TRACE",
    "LOSAT_TRACE_CHAIN_HSP",
    "LOSAT_TRACE_HSP",
    "LOSAT_TRACE_HSP_MASKS",
    "LOSAT_TRACE_LINK_SELECTIONS",
    "LOSAT_X_ADJMEMO_STATS",
    "LOSAT_X_DPCAPTURE",
    "LOSAT_X_DPCAPTURE_EVERY",
    "LOSAT_X_NEWTONEXACT_PROF",
];

/// `std::env::var_os(name)` in a diagnostic build, `None` otherwise.
///
/// No NCBI counterpart: reads a LOSAT debugging variable; it does not change any value NCBI
/// computes.
#[inline]
pub fn var_os(name: &str) -> Option<OsString> {
    debug_assert!(NAMES.contains(&name), "{name} is not in xdiag::NAMES");
    if ON {
        std::env::var_os(name)
    } else {
        None
    }
}

/// `std::env::var(name)` in a diagnostic build, `Err(VarError::NotPresent)` otherwise.
///
/// No NCBI counterpart: reads a LOSAT debugging variable; it does not change any value NCBI
/// computes.
#[inline]
pub fn var(name: &str) -> Result<String, VarError> {
    debug_assert!(NAMES.contains(&name), "{name} is not in xdiag::NAMES");
    if ON {
        std::env::var(name)
    } else {
        Err(VarError::NotPresent)
    }
}

/// The stderr line for a diagnostic variable that a build without the feature ignores.
///
/// No NCBI counterpart: LOSAT message about a LOSAT variable; it does not change any value NCBI
/// computes.
pub fn notice(name: &str) -> String {
    format!(
        "{name} is ignored: this losat was built without the diagnostics feature \
         (docs/agents/diagnostics.md, \"Diagnostic build\")"
    )
}

/// One notice per diagnostic variable that is set, in `NAMES` order; empty in a diagnostic build.
/// The loop asks `getenv` for each name, as a diagnostic build would, so the match rules of the
/// platform (case-insensitive names on Windows) are the same.
///
/// No NCBI counterpart: LOSAT message about LOSAT variables; it does not change any value NCBI
/// computes.
pub fn ignored_notices() -> Vec<String> {
    if ON {
        return Vec::new();
    }
    NAMES
        .iter()
        .filter(|name| std::env::var_os(name).is_some())
        .map(|name| notice(name))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::xswitch::{GRADUATED, STRICT_SET};
    use std::path::{Path, PathBuf};

    /// Every file under `dir` (recursively) with one of `extensions`, sorted.
    fn files_under(dir: &Path, extensions: &[&str]) -> Vec<PathBuf> {
        let mut found = Vec::new();
        let mut pending = vec![dir.to_path_buf()];
        while let Some(current) = pending.pop() {
            let Ok(entries) = std::fs::read_dir(&current) else {
                continue;
            };
            for entry in entries.flatten() {
                let path = entry.path();
                if path.is_dir() {
                    pending.push(path);
                } else if path
                    .extension()
                    .and_then(|ext| ext.to_str())
                    .is_some_and(|ext| extensions.contains(&ext))
                {
                    found.push(path);
                }
            }
        }
        found.sort();
        found
    }

    /// The text of `path` with all whitespace removed (rustfmt splits long calls over lines), or
    /// `None` for a file over 4 MiB (a fixture) or one that is not UTF-8. Reads the file once.
    fn squeezed(path: &Path) -> Option<String> {
        if std::fs::metadata(path).ok()?.len() > 4 << 20 {
            return None;
        }
        Some(squeeze(&std::fs::read_to_string(path).ok()?))
    }

    fn squeeze(text: &str) -> String {
        text.chars().filter(|c| !c.is_whitespace()).collect()
    }

    fn is_xdiag_source(path: &Path) -> bool {
        path.to_string_lossy()
            .replace('\\', "/")
            .ends_with("utils/xdiag.rs")
    }

    /// Whether `squeezed` (whitespace removed) sets or passes the variable `name` to a program:
    /// `"NAME"` or `NAME=` is left in it once the `env_remove("NAME")` calls are removed (also
    /// with the trailing comma of a call that rustfmt split). The match is of the whole name (the
    /// closing quote or `=` follows it), so `LOSAT_TRACE_HSP` is not found inside
    /// `LOSAT_TRACE_HSP_MASKS`.
    fn sets_name(squeezed: &str, name: &str) -> bool {
        let kept = squeezed
            .replace(&format!("env_remove(\"{name}\")"), "")
            .replace(&format!("env_remove(\"{name}\",)"), "");
        kept.contains(&format!("\"{name}\"")) || kept.contains(&format!("{name}="))
    }

    #[test]
    fn the_names_are_distinct_sorted_and_outside_the_switches() {
        assert!(NAMES.windows(2).all(|pair| pair[0] < pair[1]));
        for name in NAMES {
            assert!(name.starts_with("LOSAT_"), "{name}");
            assert!(!STRICT_SET.contains(name), "{name} is in STRICT_SET");
            assert!(!GRADUATED.contains(name), "{name} is in GRADUATED");
            assert!(!name.ends_with("SHADOW"), "{name} is a shadow switch");
            assert!(!name.ends_with("VERIFY"), "{name} is a verify switch");
        }
    }

    /// The block of `docs/agents/diagnostics.md` under "Diagnostic build" names the variables of
    /// `NAMES`, one per line, in the same order.
    #[test]
    fn the_documented_block_lists_the_names() {
        let doc = include_str!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/agents/diagnostics.md"
        ));
        let anchor = doc
            .find("The variables of the diagnostic build (")
            .expect("the block's heading in diagnostics.md");
        let block = &doc[anchor..];
        let start = block.find("```text\n").expect("the block") + "```text\n".len();
        let end = start + block[start..].find("```").expect("the end of the block");
        let documented: Vec<&str> = block[start..end]
            .lines()
            .map(|line| line.trim_end_matches('\r'))
            .collect();
        assert_eq!(documented, NAMES);
    }

    /// The calls that may take a listed name as a string literal: this module's two readers and
    /// the tblastx helper that passes its argument to `var`.
    const READERS: [&str; 3] = [
        "xdiag::var(",
        "xdiag::var_os(",
        "parse_trace_hsp_target_env(",
    ];

    /// No file under `src` other than this one names a listed variable except as the argument of
    /// one of `READERS` (so none reads it with `std::env`, which would bypass `ON`, or through
    /// a const or another helper), and every listed name is still used (a removed variable leaves
    /// `NAMES`).
    #[test]
    fn only_this_module_reads_the_names_from_the_environment() {
        let src = Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/src"));
        if !src.is_dir() {
            eprintln!("skipped: {} does not exist", src.display());
            return;
        }
        let mut failures = Vec::new();
        let mut mentioned = vec![false; NAMES.len()];
        for path in files_under(src, &["rs"]) {
            if is_xdiag_source(&path) {
                continue;
            }
            let Some(text) = squeezed(&path) else {
                continue;
            };
            for (index, name) in NAMES.iter().enumerate() {
                let quoted = format!("\"{name}\"");
                for (at, _) in text.match_indices(&quoted) {
                    mentioned[index] = true;
                    if !READERS.iter().any(|reader| text[..at].ends_with(reader)) {
                        failures.push(format!("{}: {quoted} is not read by xdiag", path.display()));
                    }
                }
            }
        }
        for (index, name) in NAMES.iter().enumerate() {
            if !mentioned[index] {
                failures.push(format!("{name} is used in no file under src"));
            }
        }
        assert!(failures.is_empty(), "{failures:#?}");
    }

    /// No test, script or CI step sets a listed name: a default build ignores it, so a test that
    /// did would pass in CI (`--all-features`) and fail in a default local run. Reads every `.rs`,
    /// `.py`, `.sh`, `.tsv`, `.yml` and `.yaml` file under `LOSAT/tests`, `scripts`, `tests` and
    /// `.github/workflows` (a directory that does not exist is skipped); `env_remove` is fine. The
    /// two checks of the notice itself set `LOSAT_STARTUP_TRACE` on purpose: `cli_smoke.rs` passes
    /// in both builds, `ci_fast_regressions.py` checks the default release binary.
    #[test]
    fn no_test_script_or_workflow_sets_a_name() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let mut failures = Vec::new();
        for dir in ["LOSAT/tests", "scripts", "tests", ".github/workflows"] {
            let dir = root.join(dir);
            if !dir.is_dir() {
                eprintln!("skipped: {} does not exist", dir.display());
                continue;
            }
            for path in files_under(&dir, &["rs", "py", "sh", "tsv", "yml", "yaml"]) {
                let slashed = path.to_string_lossy().replace('\\', "/");
                let notice_test = slashed.ends_with("tests/cli_smoke.rs")
                    || slashed.ends_with("tests/ci_fast_regressions.py");
                let Some(text) = squeezed(&path) else {
                    continue;
                };
                for name in NAMES {
                    if notice_test && *name == "LOSAT_STARTUP_TRACE" {
                        continue;
                    }
                    if sets_name(&text, name) {
                        failures.push(format!("{}: {name}", path.display()));
                    }
                }
            }
        }
        assert!(failures.is_empty(), "{failures:#?}");
    }

    #[test]
    fn the_matching_of_the_reverse_guard() {
        let sets = |text: &str, name: &str| sets_name(&squeeze(text), name);
        assert!(sets("(\"LOSAT_TRACE_HSP\", \"1\"),", "LOSAT_TRACE_HSP"));
        assert!(sets("env=LOSAT_TRACE_HSP=1 losat", "LOSAT_TRACE_HSP"));
        assert!(!sets(
            "cmd.env_remove(\"LOSAT_TRACE_HSP\");",
            "LOSAT_TRACE_HSP"
        ));
        assert!(!sets(
            "cmd.env_remove(\n    \"LOSAT_TRACE_HSP\",\n);",
            "LOSAT_TRACE_HSP"
        ));
        assert!(!sets(
            "(\"LOSAT_TRACE_HSP_MASKS\", \"1\"), LOSAT_TRACE_HSP_MASKS=1",
            "LOSAT_TRACE_HSP"
        ));
        assert!(sets(
            "(\"LOSAT_TRACE_HSP_MASKS\", \"1\")",
            "LOSAT_TRACE_HSP_MASKS"
        ));
    }

    #[test]
    fn the_notice_names_the_variable_and_the_docs() {
        assert_eq!(
            notice("LOSAT_TRACE_HSP"),
            "LOSAT_TRACE_HSP is ignored: this losat was built without the diagnostics feature \
             (docs/agents/diagnostics.md, \"Diagnostic build\")"
        );
    }
}
