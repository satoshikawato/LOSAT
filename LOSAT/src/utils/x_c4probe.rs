//! PROBE (v4-C s4, not for merge): call counts and repeated inputs of the blastp redo's
//! "other" functions. Active only when `LOSAT_C4PROBE` is set; the counts go to stderr at exit.
use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::sync::Mutex;
use std::sync::OnceLock;

pub fn on() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("LOSAT_C4PROBE").is_some())
}

#[derive(Default)]
struct Site {
    calls: u64,
    units: u64,
    units_first: u64,
    distinct: HashSet<u64>,
}

static SITES: Mutex<Option<HashMap<(&'static str, &'static str), Site>>> = Mutex::new(None);

pub fn key_of<T: Hash + ?Sized>(value: &T) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

pub fn key_f64s(values: &[f64]) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    for v in values {
        v.to_bits().hash(&mut hasher);
    }
    hasher.finish()
}

/// Count one call at (group, site) with `units` of work and an input key.
pub fn note(group: &'static str, site: &'static str, units: u64, key: u64) {
    if !on() {
        return;
    }
    let mut guard = SITES.lock().expect("c4probe lock");
    let map = guard.get_or_insert_with(HashMap::new);
    for s in [site, "*"] {
        let entry = map.entry((group, s)).or_default();
        entry.calls += 1;
        entry.units += units;
        if entry.distinct.insert(key) {
            entry.units_first += units;
        }
    }
}

pub fn dump() {
    if !on() {
        return;
    }
    let mut guard = SITES.lock().expect("c4probe lock");
    let Some(map) = guard.take() else {
        return;
    };
    let mut rows: Vec<_> = map.iter().collect();
    rows.sort_by(|a, b| a.0.cmp(b.0));
    for ((group, site), s) in rows {
        eprintln!(
            "[C4PROBE] {group}\t{site}\tcalls={}\tunits={}\tunits_first={}\tdistinct={}",
            s.calls,
            s.units,
            s.units_first,
            s.distinct.len()
        );
    }
}
