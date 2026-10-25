//! EXPERIMENT (line C s7, branch exp/c8-bx-newton-probe, not for merging). No NCBI counterpart: a
//! log of the Newton problems that run exactly (`LOSAT_C8NEWTONLOG=<path>`: one file per thread,
//! `<path>.<pid>.<n>`), each with the BLASTX match (one `redo_list_observed` call) and frame
//! (context) it belongs to and the time of its exact solve; and the ignored test `c8_newton_replay`
//! that times the same problems solved one at a time (`optimize_target_frequencies_exact`) and in
//! lane batches grouped as a prefetch could group them (per frame call, per match, per k matches).
//! It does not change any value NCBI computes.

use std::cell::{Cell, RefCell};
use std::io::Write;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::OnceLock;

const N: usize = 20;
const NN: usize = N * N;
/// Bytes of one record: match, frame, status, pad (u32 each), ns (u64), relative entropy, q, row, col.
const REC: usize = 4 * 4 + 8 + 8 * (1 + NN + 2 * N);

static PATH: OnceLock<Option<String>> = OnceLock::new();
static THREADS: AtomicU32 = AtomicU32::new(0);

thread_local! {
    static MATCH: Cell<u32> = const { Cell::new(0) };
    static FRAME: Cell<u32> = const { Cell::new(0) };
    static FILE: RefCell<Option<std::fs::File>> = const { RefCell::new(None) };
}

pub(crate) fn on() -> bool {
    PATH.get_or_init(|| std::env::var("LOSAT_C8NEWTONLOG").ok().filter(|v| !v.is_empty()))
        .is_some()
}

pub(crate) fn next_match() {
    MATCH.with(|m| m.set(m.get() + 1));
}

pub(crate) fn set_frame(frame: usize) {
    FRAME.with(|f| f.set(frame as u32));
}

pub(crate) fn record(q: &[f64], row: &[f64; N], col: &[f64; N], re: f64, ns: u64, status: i32) {
    let Some(Some(path)) = PATH.get() else {
        return;
    };
    let mut buf = Vec::with_capacity(REC);
    buf.extend_from_slice(&MATCH.with(|m| m.get()).to_le_bytes());
    buf.extend_from_slice(&FRAME.with(|f| f.get()).to_le_bytes());
    buf.extend_from_slice(&status.to_le_bytes());
    buf.extend_from_slice(&0u32.to_le_bytes());
    buf.extend_from_slice(&ns.to_le_bytes());
    buf.extend_from_slice(&re.to_le_bytes());
    for v in q.iter().chain(row.iter()).chain(col.iter()) {
        buf.extend_from_slice(&v.to_le_bytes());
    }
    assert_eq!(buf.len(), REC);
    FILE.with(|f| {
        let mut f = f.borrow_mut();
        if f.is_none() {
            let n = THREADS.fetch_add(1, Ordering::Relaxed);
            let name = format!("{path}.{}.{n}", std::process::id());
            *f = Some(std::fs::File::create(name).expect("create the Newton log"));
        }
        f.as_mut()
            .unwrap()
            .write_all(&buf)
            .expect("write the Newton log");
    });
}

#[cfg(test)]
mod tests {
    use super::super::x_newton_exact::optimize_target_frequencies_exact;
    use super::super::x_newton_lanes::{solve_batch, NewtonInput, NewtonOutput};
    use super::*;
    use std::time::Instant;

    struct Rec {
        m: u32,
        frame: u32,
        status: i32,
        ns: u64,
        input: NewtonInput,
    }

    fn f64s(b: &[u8]) -> impl Iterator<Item = f64> + '_ {
        b.chunks_exact(8)
            .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
    }

    fn read(path: &str) -> Vec<Rec> {
        let bytes = std::fs::read(path).expect("read C8_NEWTON_LOG");
        assert_eq!(bytes.len() % REC, 0);
        bytes
            .chunks_exact(REC)
            .map(|r| {
                let u = |i: usize| u32::from_le_bytes(r[4 * i..4 * i + 4].try_into().unwrap());
                let ns = u64::from_le_bytes(r[16..24].try_into().unwrap());
                let mut it = f64s(&r[24..]);
                let re = it.next().unwrap();
                let mut q = [0.0; NN];
                let mut row = [0.0; N];
                let mut col = [0.0; N];
                for v in q.iter_mut().chain(row.iter_mut()).chain(col.iter_mut()) {
                    *v = it.next().unwrap();
                }
                Rec {
                    m: u(0),
                    frame: u(1),
                    status: u(2) as i32,
                    ns,
                    input: NewtonInput {
                        q,
                        row,
                        col,
                        relative_entropy: re,
                    },
                }
            })
            .collect()
    }

    fn hist(sizes: &[usize]) -> String {
        let mut h = [0usize; 10];
        for &s in sizes {
            h[s.min(9)] += 1;
        }
        let n = sizes.len().max(1);
        let mean = sizes.iter().sum::<usize>() as f64 / n as f64;
        format!(
            "groups={} mean={mean:.2} 1:{} 2:{} 3:{} 4:{} 5:{} 6:{} 7:{} 8:{} 9+:{}",
            sizes.len(),
            h[1],
            h[2],
            h[3],
            h[4],
            h[5],
            h[6],
            h[7],
            h[8],
            h[9]
        )
    }

    /// Groups of consecutive records with the same key.
    fn groups(recs: &[Rec], key: impl Fn(&Rec) -> u64) -> Vec<std::ops::Range<usize>> {
        let mut out = Vec::new();
        let mut start = 0;
        for i in 1..=recs.len() {
            if i == recs.len() || key(&recs[i]) != key(&recs[start]) {
                out.push(start..i);
                start = i;
            }
        }
        out
    }

    #[test]
    #[ignore = "needs C8_NEWTON_LOG=<log file>"]
    fn c8_newton_replay() {
        let Ok(path) = std::env::var("C8_NEWTON_LOG") else {
            return;
        };
        let reps: usize = std::env::var("C8_REPS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(1);
        let recs = read(&path);
        let n = recs.len();
        let logged_ns: u64 = recs.iter().map(|r| r.ns).sum();
        // duplicates (same input bits)
        let mut seen = std::collections::HashSet::new();
        let mut dup = 0usize;
        for r in &recs {
            let mut k = Vec::with_capacity(441);
            k.push(r.input.relative_entropy.to_bits());
            k.extend(r.input.q.iter().map(|v| v.to_bits()));
            k.extend(r.input.row.iter().map(|v| v.to_bits()));
            k.extend(r.input.col.iter().map(|v| v.to_bits()));
            if !seen.insert(k) {
                dup += 1;
            }
        }
        let by_frame = groups(&recs, |r| ((r.m as u64) << 8) | r.frame as u64);
        let by_match = groups(&recs, |r| r.m as u64);
        println!(
            "C8N problems={n} matches_with_newton={} duplicates={dup} logged_exact_us_per_problem={:.2}",
            by_match.len(),
            logged_ns as f64 / n.max(1) as f64 / 1e3
        );
        println!(
            "C8N per frame call: {}",
            hist(&by_frame.iter().map(|g| g.len()).collect::<Vec<_>>())
        );
        println!(
            "C8N per match: {}",
            hist(&by_match.iter().map(|g| g.len()).collect::<Vec<_>>())
        );
        // exact, one at a time
        let mut want = vec![NewtonOutput::new(); n];
        let t = Instant::now();
        for _ in 0..reps {
            for (r, o) in recs.iter().zip(want.iter_mut()) {
                o.status = optimize_target_frequencies_exact(
                    &mut o.x,
                    &r.input.q,
                    &r.input.row,
                    &r.input.col,
                    r.input.relative_entropy,
                );
            }
        }
        let exact_s = t.elapsed().as_secs_f64() / reps as f64;
        for (r, o) in recs.iter().zip(&want) {
            assert_eq!(r.status, o.status, "replayed status");
        }
        println!(
            "C8N exact: {:.2} us/problem",
            exact_s * 1e6 / n.max(1) as f64
        );
        let inputs: Vec<NewtonInput> = recs
            .iter()
            .map(|r| NewtonInput {
                q: r.input.q,
                row: r.input.row,
                col: r.input.col,
                relative_entropy: r.input.relative_entropy,
            })
            .collect();
        let mut run = |label: &str, gs: &[std::ops::Range<usize>]| {
            let mut out = vec![NewtonOutput::new(); n];
            let t = Instant::now();
            for _ in 0..reps {
                for g in gs {
                    solve_batch(&inputs[g.clone()], &mut out[g.clone()]);
                }
            }
            let s = t.elapsed().as_secs_f64() / reps as f64;
            for (i, (a, b)) in out.iter().zip(&want).enumerate() {
                assert_eq!(a.status, b.status, "{label}: status of problem {i}");
                assert!(
                    a.x.iter().zip(b.x.iter()).all(|(u, v)| u.to_bits() == v.to_bits()),
                    "{label}: x of problem {i}"
                );
            }
            println!(
                "C8N lanes {label}: {:.2} us/problem, x{:.3} of exact",
                s * 1e6 / n.max(1) as f64,
                s / exact_s
            );
        };
        run("per frame call", &by_frame);
        run("per match", &by_match);
        for k in [4usize, 16, 64] {
            let gs: Vec<_> = by_match
                .chunks(k)
                .map(|c| c[0].start..c[c.len() - 1].end)
                .collect();
            run(&format!("per {k} matches"), &gs);
        }
        let gs: Vec<_> = (0..n).step_by(32).map(|s| s..(s + 32).min(n)).collect();
        run("32 at a time", &gs);
    }
}
