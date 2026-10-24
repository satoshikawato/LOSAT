//! The sides of the conditional switches that a run took (`LOSAT_X_REACH`,
//! `LOSAT/src/utils/xreach.rs`): a test passes `Reach::var()` with the variables of a run and
//! asks `Reach::sides` whether the switch's own path ran on its input. It replaced the shadow
//! counters as the proof that a path engaged (switch-tidy, 2026-10-24).
#![allow(dead_code)]

use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

/// A directory that the runs given its `var()` write their `reach-*.tsv` files into (one
/// directory per call of `new`); removed on drop.
pub struct Reach {
    dir: PathBuf,
}

impl Reach {
    /// A new empty directory per call: `temp_dir()/losat-reach-<label>-<pid>-<n>`. The tests of
    /// one binary run at the same time, so no two calls share a directory.
    pub fn new(label: &str) -> Reach {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let n = NEXT.fetch_add(1, Ordering::Relaxed);
        let dir =
            std::env::temp_dir().join(format!("losat-reach-{label}-{}-{n}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("reach directory");
        Reach { dir }
    }

    /// `("LOSAT_X_REACH", <directory>)`, to add to the variables of a run.
    pub fn var(&self) -> (&'static str, &str) {
        ("LOSAT_X_REACH", self.dir.to_str().expect("UTF-8 path"))
    }

    /// `(new, ref)` of `switch` over every `reach-*.tsv` of the directory (the OR of its rows);
    /// `(false, false)` when the runs recorded no row for it; panics when no run wrote a reach
    /// file (a run that got no `LOSAT_X_REACH`, or could not write).
    pub fn sides(&self, switch: &str) -> (bool, bool) {
        let (mut new, mut reference, mut files) = (false, false, 0usize);
        for entry in std::fs::read_dir(&self.dir).expect("reach directory") {
            let path = entry.expect("reach entry").path();
            let name = path
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("");
            if !(name.starts_with("reach-") && name.ends_with(".tsv")) {
                continue;
            }
            files += 1;
            let text = std::fs::read_to_string(&path).expect("reach file");
            for row in text.lines().skip(1) {
                let mut fields = row.split('\t');
                if fields.next() == Some(switch) {
                    new |= fields.next() == Some("1");
                    reference |= fields.next() == Some("1");
                }
            }
        }
        assert!(
            files > 0,
            "no reach-*.tsv in {}: the run did not get LOSAT_X_REACH ({switch}) or could not write",
            self.dir.display()
        );
        (new, reference)
    }
}

impl Drop for Reach {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}
