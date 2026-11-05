//! Diagnostic (`LOSAT_X_REACH`): which side of each conditional switch of the strict set, or of
//! one that graduated from it, a run took.
//!
//! No NCBI counterpart: it records which LOSAT code path ran; it changes no output and computes
//! no value NCBI computes.
//!
//! A conditional switch (class `conditional` in `LOSAT/tests/switch_reach_expected.tsv`) takes
//! its own path only where a rule holds (input sizes, thread count, options, memory); elsewhere
//! the reference path runs although the switch is on. A graduated switch (`GRADUATED` in
//! `utils/xswitch.rs`) has no variable and is always on; its rule and its marks stay. With
//! `LOSAT_X_REACH=<absolute directory>` each of its decision points records the side it took
//! (`x_reach!`), and the run writes one file
//! `reach-<pid>-<nanoseconds>.tsv` into the directory at its end (`write_report`, called by
//! `main`). `LOSAT/tests/switch_reach_check.py` folds the files of the CI fixtures and fails
//! where a side the expected file requires was not reached by any fixture.
//!
//! A decision point records nothing while its switch is off, and without the variable it costs
//! one load of a cached flag. It sits where the choice is made (per search, batch, subject, unit
//! or call), never inside a per-hit loop.

use std::ffi::OsString;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU8, Ordering};
use std::sync::OnceLock;

use super::xswitch::{GRADUATED, STRICT_SET};

/// The switch's own path ran.
const NEW: u8 = 1;
/// The reference path ran with the switch on (its rule did not hold).
const REF: u8 = 2;

/// The switches that can record a side: those of `STRICT_SET`, then those of `GRADUATED`.
const SWITCHES: usize = STRICT_SET.len() + GRADUATED.len();

static SIDES: [AtomicU8; SWITCHES] = [const { AtomicU8::new(0) }; SWITCHES];

/// No NCBI counterpart: the position of `name` in `STRICT_SET`, else `STRICT_SET.len()` plus its
/// position in `GRADUATED`. Evaluated in a constant by `x_reach!`, so a name outside both lists
/// does not compile.
#[doc(hidden)]
pub const fn index(name: &str) -> usize {
    let mut i = 0;
    while i < STRICT_SET.len() {
        if same(STRICT_SET[i].as_bytes(), name.as_bytes()) {
            return i;
        }
        i += 1;
    }
    let mut g = 0;
    while g < GRADUATED.len() {
        if same(GRADUATED[g].as_bytes(), name.as_bytes()) {
            return STRICT_SET.len() + g;
        }
        g += 1;
    }
    panic!("x_reach!: not in STRICT_SET or GRADUATED (utils/xswitch.rs)");
}

const fn same(a: &[u8], b: &[u8]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let mut i = 0;
    while i < a.len() {
        if a[i] != b[i] {
            return false;
        }
        i += 1;
    }
    true
}

/// No NCBI counterpart: `LOSAT_X_REACH` read once (`reach_directory`).
fn setting() -> &'static Result<Option<PathBuf>, OsString> {
    static SETTING: OnceLock<Result<Option<PathBuf>, OsString>> = OnceLock::new();
    SETTING.get_or_init(|| reach_directory(std::env::var_os("LOSAT_X_REACH")))
}

/// No NCBI counterpart: the directory of `LOSAT_X_REACH`; `None` when the variable is not set or
/// empty, or is not an absolute path (`write_report` then says so on stderr, once, at the end of
/// the run: a decision point may run on a pool thread, which must not take the stderr lock).
fn directory() -> Option<&'static PathBuf> {
    setting().as_ref().ok().and_then(Option::as_ref)
}

/// No NCBI counterpart: the directory a value of `LOSAT_X_REACH` names. `Ok(None)` when it is
/// not set or empty; `Err(value)` when it is not an absolute path. A relative value would be a
/// directory under the run's current directory: on 2026-10-31 runs with the values `0` and `1`
/// left `LOSAT/0/` and `LOSAT/1/` in a checkout (12 files, removed by #715).
fn reach_directory(value: Option<OsString>) -> Result<Option<PathBuf>, OsString> {
    let Some(value) = value.filter(|value| !value.is_empty()) else {
        return Ok(None);
    };
    let directory = PathBuf::from(&value);
    if directory.is_absolute() {
        Ok(Some(directory))
    } else {
        Err(value)
    }
}

/// No NCBI counterpart: records that the switch at `index` took its own path (`own_path`) or
/// the reference path. Called through `x_reach!`. Never inlined: a decision point keeps only the
/// call, so the code around it compiles as it did without the mark.
#[doc(hidden)]
#[inline(never)]
pub fn mark(index: usize, own_path: bool) {
    if directory().is_some() {
        SIDES[index].fetch_or(if own_path { NEW } else { REF }, Ordering::Relaxed);
    }
}

/// Records the side a conditional switch of the strict set (or a graduated one) took, while the
/// switch is on:
/// `x_reach!("LOSAT_X_<NAME>", own_path)`, with `own_path` true where the switch's path runs
/// and false where its rule sends the work to the reference path.
macro_rules! x_reach {
    ($name:literal, $own_path:expr) => {{
        const INDEX: usize = $crate::utils::xreach::index($name);
        $crate::utils::xreach::mark(INDEX, $own_path)
    }};
}
pub(crate) use x_reach;

/// No NCBI counterpart: writes the sides this run took into the `LOSAT_X_REACH` directory, one
/// row per switch that recorded a side (`switch`, `new`, `ref`, `program`: 1 = reached). A
/// failure to write is reported on stderr and does not change the exit status; so is a value
/// that is not an absolute path, which writes nothing.
pub fn write_report() {
    let directory = match setting() {
        Ok(Some(directory)) => directory,
        Ok(None) => return,
        Err(value) => {
            eprintln!(
                "[X_REACH] ignored: LOSAT_X_REACH={} is not an absolute path",
                value.to_string_lossy()
            );
            return;
        }
    };
    let program = crate::cli::entry_argv(std::env::args_os())
        .get(1)
        .map_or_else(String::new, |word| word.to_string_lossy().into_owned());
    let mut text = String::from("switch\tnew\tref\tprogram\n");
    for (name, sides) in STRICT_SET.iter().chain(GRADUATED).zip(SIDES.iter()) {
        let sides = sides.load(Ordering::Relaxed);
        if sides != 0 {
            text.push_str(&format!(
                "{name}\t{}\t{}\t{program}\n",
                u8::from(sides & NEW != 0),
                u8::from(sides & REF != 0)
            ));
        }
    }
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |elapsed| elapsed.as_nanos());
    #[cfg(not(target_os = "wasi"))]
    let pid = std::process::id();
    #[cfg(target_os = "wasi")]
    let pid = 0;
    let path = directory.join(format!("reach-{pid}-{nanos}.tsv"));
    let written = std::fs::create_dir_all(directory).and_then(|()| std::fs::write(&path, text));
    if let Err(error) = written {
        eprintln!("[X_REACH] cannot write {}: {error}", path.display());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_finds_every_name_of_the_set_then_the_graduated_ones() {
        for (i, name) in STRICT_SET.iter().chain(GRADUATED).enumerate() {
            assert_eq!(index(name), i);
        }
    }
}
