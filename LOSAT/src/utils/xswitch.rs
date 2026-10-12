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

/// No NCBI counterpart: the switches of the strict set, the variables that are on by default.
/// `LOSAT_X_SEEDBUCKET` is not among them. When the set became the default it numbered the lookup
/// hits of a chunk scan in a `u32`, exact only up to 2^32 hits per scan; it now numbers them per
/// flush (`tblastx/x_seed_bucket.rs`) and stays opt-in until it is measured on the searches it
/// would act on with the set on (the fallbacks of `LOSAT_X_PAIRPAR`).
pub(crate) const STRICT_SET: [&str; 43] = [
    "LOSAT_X_DPFAST",
    "LOSAT_X_COMPFAST",
    "LOSAT_X_PURGEFAST",
    "LOSAT_X_ITREEFAST",
    "LOSAT_X_CTXFAST",
    "LOSAT_X_GREEDYFAST",
    "LOSAT_X_GREEDYSPEC",
    "LOSAT_X_AHEAD",
    "LOSAT_X_ENVCACHE",
    "LOSAT_X_DUSTFAST",
    "LOSAT_X_DUSTRING",
    "LOSAT_X_MBDELAY",
    "LOSAT_X_MBBATCH",
    "LOSAT_X_SEGFAST",
    "LOSAT_X_SEGMEMO",
    "LOSAT_X_SEGSHARE",
    "LOSAT_X_KARLINFAST",
    "LOSAT_X_IDEALMEMO",
    "LOSAT_X_ERFMEMO",
    "LOSAT_X_DEKKER",
    "LOSAT_X_LUTARENA",
    "LOSAT_X_LUTSPLIT",
    "LOSAT_X_TBNPAR",
    "LOSAT_X_TBNEVENTS",
    "LOSAT_X_TBNQSIDE",
    "LOSAT_X_TBNSSIDE",
    "LOSAT_X_TBNBATCH",
    "LOSAT_X_BXPAR",
    "LOSAT_X_BXLAZYCTX",
    "LOSAT_X_BXLEAN",
    "LOSAT_X_BXPOOL",
    "LOSAT_X_BXCHUNK",
    "LOSAT_X_BXBATCH",
    "LOSAT_X_CODONFAST",
    "LOSAT_X_NEWTONEXACT",
    "LOSAT_X_DPAVX2",
    "LOSAT_X_DPW128",
    "LOSAT_X_BXLUT",
    "LOSAT_X_BXSCAN",
    "LOSAT_X_BXSEGMEMO",
    "LOSAT_X_PAIRPAR",
    "LOSAT_X_NEWTONLANES",
    "LOSAT_X_GREEDYPRELIM",
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
    fn the_strict_set_names_are_distinct() {
        let mut names = STRICT_SET.to_vec();
        names.sort_unstable();
        names.dedup();
        assert_eq!(names.len(), STRICT_SET.len());
        assert!(names.iter().all(|name| name.starts_with("LOSAT_X_")));
    }
}
