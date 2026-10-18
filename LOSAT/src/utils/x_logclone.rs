//! `log_slice`: the natural logarithm of every element, through libm.
//!
//! PROBE BRANCH exp/c2-log-probe (not for merge): with `C2PROBE_DUMP=<dir>` every distinct
//! input bit pattern of `log_slice` is collected and written to `<dir>/<pid>.bin` (little-endian
//! u64) by `print_shadow_stats` at exit, together with the number of calls and values.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/optimize_target_freq.c:636,679
//! ```c
//! temp = log(x[k] / q[k]);
//! ...
//!     scores[k] = log(target_freqs[k] / (row_freqs[i] * col_freqs[j]));
//! ```

use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock};

static VALUES: AtomicU64 = AtomicU64::new(0);
static CALLS: AtomicU64 = AtomicU64::new(0);

fn dump_set() -> Option<&'static Mutex<rustc_hash::FxHashSet<u64>>> {
    static SET: OnceLock<Option<Mutex<rustc_hash::FxHashSet<u64>>>> = OnceLock::new();
    SET.get_or_init(|| {
        std::env::var_os("C2PROBE_DUMP").map(|_| Mutex::new(rustc_hash::FxHashSet::default()))
    })
    .as_ref()
}

/// Writes the collected inputs (probe).
pub fn print_shadow_stats() {
    let Some(set) = dump_set() else { return };
    let dir = std::env::var_os("C2PROBE_DUMP").unwrap();
    let set = set.lock().unwrap();
    let mut path = std::path::PathBuf::from(dir);
    path.push(format!("{}.bin", std::process::id()));
    let mut f = std::io::BufWriter::new(std::fs::File::create(&path).unwrap());
    for &b in set.iter() {
        f.write_all(&b.to_le_bytes()).unwrap();
    }
    f.flush().unwrap();
    let mut meta = path.clone();
    meta.set_extension("txt");
    std::fs::write(
        meta,
        format!(
            "calls={} values={} distinct={}\n",
            CALLS.load(Ordering::Relaxed),
            VALUES.load(Ordering::Relaxed),
            set.len()
        ),
    )
    .unwrap();
}

/// `out[k] = ln(src[k])`.
#[inline]
pub(crate) fn log_slice(src: &[f64], out: &mut [f64]) {
    debug_assert_eq!(src.len(), out.len());
    if let Some(set) = dump_set() {
        CALLS.fetch_add(1, Ordering::Relaxed);
        VALUES.fetch_add(src.len() as u64, Ordering::Relaxed);
        let mut set = set.lock().unwrap();
        for &x in src {
            set.insert(x.to_bits());
        }
    }
    for (o, &x) in out.iter_mut().zip(src) {
        *o = x.ln();
    }
}
