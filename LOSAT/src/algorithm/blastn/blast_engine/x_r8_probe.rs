//! EXPERIMENT (R8_PROBE; line R's R8 probe of plan v4, `v4-R/DESIGN_R8.md`; not for merge): the
//! cost of `QueryContextIndex::context_for_offset` per lookup hit of the direct-table scan
//! (blastn-short), measured by removing it and by adding it once more, and exact counts of the
//! scan's hits. Every arm returns the context `context_for_offset` returns for the offset of a
//! lookup hit (that offset lies inside a context's sequence), so no value changes.
//!
//! `R8_PROBE` is a comma-separated list:
//! - `ctx=div`: `bsearch_context_info`, NCBI's BSearchContextInfo with its b/e narrowing
//!   (blast_query_info.c:219-243): one division and no loop when every context has one length.
//! - `ctx=tab<k>` (k = 2..12): the first run after every 2^k offsets plus a forward walk
//!   (R4's table at k = 8).
//! - `ctx=dbl`: `context_for_offset`, then once more on an opaque copy of the offset.
//! - `cnt`: count lookup hits, hits rejected by the diagonal, ungapped calls and hits at or above
//!   the cutoff; one `[R8] cnt` line per scan on stderr.
//!
//! Unset or empty = the dev code path.
use super::QueryContextIndex;
use std::sync::OnceLock;

#[derive(Clone, Copy, Debug, Default)]
pub(super) struct Probe {
    /// 0 dev, 1 `div`, 2 `tab<k>`, 3 `dbl`.
    pub(super) ctx: u8,
    pub(super) tab_k: u32,
    pub(super) cnt: bool,
}

pub(super) fn probe() -> Probe {
    static P: OnceLock<Probe> = OnceLock::new();
    *P.get_or_init(|| {
        let mut p = Probe::default();
        if let Ok(v) = std::env::var("R8_PROBE") {
            for t in v.split(',').map(str::trim) {
                match t {
                    "" => {}
                    "ctx=div" => p.ctx = 1,
                    "ctx=dbl" => p.ctx = 3,
                    "cnt" => p.cnt = true,
                    _ if t.starts_with("ctx=tab") => {
                        let k: u32 = t[7..].parse().expect("R8_PROBE ctx=tab<k>");
                        assert!((2..=12).contains(&k), "R8_PROBE ctx=tab<k>: k in 2..=12");
                        p.ctx = 2;
                        p.tab_k = k;
                    }
                    _ => panic!("R8_PROBE: unknown item {t}"),
                }
            }
            eprintln!("[R8_PROBE] {p:?}");
        }
        p
    })
}

/// The run index after every `2^k` offsets: `after[b]` = the number of runs that start at or
/// before `b << k`.
pub(super) struct CtxTab {
    after: Vec<u32>,
    k: u32,
}

impl CtxTab {
    pub(super) fn build(q: &QueryContextIndex, k: u32) -> Self {
        let blocks = (q.x_total_len >> k) + 2;
        let mut after = Vec::with_capacity(blocks);
        let mut a = 0usize;
        for b in 0..blocks {
            let start = b << k;
            while a < q.x_ranges.len() && q.x_ranges[a].0 <= start {
                a += 1;
            }
            after.push(a as u32);
        }
        Self { after, k }
    }

    /// `q.context_for_offset(n)`.
    #[inline(always)]
    fn context(&self, q: &QueryContextIndex, n: usize) -> usize {
        if n >= q.x_total_len || !q.direct_map.is_empty() {
            return q.context_for_offset(n);
        }
        let r = &q.x_ranges;
        let mut after = self.after[n >> self.k] as usize;
        while after < r.len() && r[after].0 <= n {
            after += 1;
        }
        if after == 0 {
            return 0;
        }
        let (_, end, context) = r[after - 1];
        if n < end {
            context as usize
        } else {
            0
        }
    }
}

/// The context of the lookup hit at concatenated query offset `n` by the probe's arm `ctx`
/// (`ctx` != 0).
#[inline(always)]
pub(super) fn context(ctx: u8, tab: Option<&CtxTab>, q: &QueryContextIndex, n: usize) -> usize {
    match ctx {
        1 => super::bsearch_context_info(n, q),
        2 => tab.expect("R8_PROBE ctx=tab table").context(q, n),
        _ => {
            let c = q.context_for_offset(n);
            let c2 = q.context_for_offset(std::hint::black_box(n));
            std::hint::black_box(c2);
            c
        }
    }
}

/// Counts of one scan (`cnt`).
#[derive(Default)]
pub(super) struct Counts {
    pub(super) hits: u64,
    pub(super) rejected: u64,
    pub(super) ungapped: u64,
    pub(super) passed: u64,
}

impl Counts {
    pub(super) fn print(&self, contexts: usize, s_len: usize) {
        eprintln!(
            "[R8] cnt hits={} rejected={} ungapped={} passed={} contexts={} s_len={}",
            self.hits, self.rejected, self.ungapped, self.passed, contexts, s_len
        );
    }
}
