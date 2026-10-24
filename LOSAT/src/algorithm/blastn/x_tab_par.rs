//! EXPERIMENT (LOSAT_X_BNTABPAR): BLASTN's outfmt 6 / 7 rows of a report formatted on the search's
//! pool.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:828-832
//! ```c
//! ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) {
//!         const CSeq_align& s = **itr;
//!         tabinfo.SetFields(s, *m_Scope, &m_ScoringMatrix);
//!         tabinfo.Print();
//! }
//! ```
//! NCBI prints the rows of a query's alignments one after another, and the queries of a batch one
//! after another, on the thread of the application's loop. `write_output_blastn_hitlists_to_writer`
//! (`hsp.rs`) does the same on the thread of `run_in_pool`'s batch loop, while the other threads of
//! the search's pool wait: 31,168 cDNAs against 3 Mb with `-task blastn` (1.02 M rows, 90 query
//! batches written one at a time by LOSAT_X_RPTSTREAM, LOSAT_X_QBATCHPAR=0) spent 0.88 s in the
//! report calls at 8, 16 and 32 threads alike (ws, 2026-10-24), 87 % of it in the rows
//! (`write_hit_fields_bytes`, local perf).
//!
//! What the switch changes (read once per process, `mode`), when the caller gives a pool of two
//! threads or more (`pool`: the search's pool while no LOSAT_X_QBATCHPAR slot exists), there is no
//! format observer (`FormatProbe`, which flushes around every row; the Web adapter) and the
//! BLASTN trace is off (`tracing::enabled`; it logs rows in the serial order): the queries are
//! taken in windows of at least `ROWS_PER_THREAD` rows per thread of the pool (at least one query),
//! and the pool formats each query of a window into a buffer of its own with
//! `write_query_rows`, the serial function of one query; the calling thread then writes the
//! window's buffers in query order, each query's warnings (`QueryWarnings`) before its buffer,
//! before the next window is formatted. The last line (`write_tabular_epilog`) stays on the
//! calling thread.
//!
//! Why the bytes are those of the serial writer: each buffer holds what `write_query_rows` writes
//! for that query, into a buffer instead of the writer, and the serial writer is that call for each
//! query in order with the query's warnings before it (`QueryWarnings::before_query`, which flushes
//! the report written so far and then writes to its own sink, as it does here). A query's rows read
//! only its own hit list, the ids, the titles and the fixed configuration; no value is carried from
//! one query to the next (the running HSP index is read only by the format observer, which keeps
//! the serial writer). Writing to a buffer cannot fail, so a failed write is one of the writer's,
//! which stops the report as before.
//!
//! Off (the serial writer): the switch off (`LOSAT_X_BNTABPAR=0`); no pool or one of one thread
//! (`-num_threads 1`, builds without the `parallel` feature, Wasm without threads); while
//! LOSAT_X_QBATCHPAR slots exist (their pools would share the CPUs; their reports already overlap
//! the other slots' searches); a format observer; the BLASTN trace.
//!
//! Memory: the formatted rows of one window (`ROWS_PER_THREAD` x threads rows, about 6 MiB at 32
//! threads with the ~100-byte rows of the cDNA input) plus the buffer of the window's last query.
//!
//! With LOSAT_X_REACH the decision is recorded once per report call (`x_reach!`).
//!
//! No NCBI counterpart for the windows and buffers: they choose where the same bytes are made and
//! change no value NCBI computes.

use std::io::{self, Write};
use std::sync::{Arc, OnceLock};

use crate::api::local_blast::FormatProbe;
use crate::report::query_warnings::QueryWarnings;
use crate::report::OutputConfig;
use crate::utils::threading::SearchPool;

use super::hsp::{
    write_output_blastn_hitlists_to_writer, write_query_rows, write_tabular_epilog, BlastnHitList,
    BlastnOutputFormat,
};

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;

/// Rows per thread of the pool in one window of queries.
const ROWS_PER_THREAD: usize = 2048;

/// Bytes reserved per row of a query's buffer (a tabular row of the cDNA input is ~100 bytes).
const ROW_BYTES: usize = 96;

// No NCBI counterpart: reads the switch once (LOSAT_X_BNTABPAR is a strict-set switch: on unless
// it is `0`). It chooses between two LOSAT code paths with the same result.
fn mode() -> bool {
    static MODE: OnceLock<bool> = OnceLock::new();
    *MODE.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNTABPAR"))
}

/// The pool the tabular rows are formatted on (`Some` with the switch on and a pool of two
/// threads or more), or `None` for the serial writer.
pub(crate) fn pool<'p, 'q>(pool: &'p SearchPool<'q>) -> Option<&'p SearchPool<'q>> {
    const PARALLEL: bool = cfg!(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ));
    (PARALLEL && mode() && pool.enabled() && pool.threads() > 1).then_some(pool)
}

/// `write_output_blastn_hitlists_to_writer`'s report: its rows formatted on `pool` (from `pool`)
/// where the switch's rule allows, else by the serial writer.
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_output<W: Write>(
    pool: Option<&SearchPool<'_>>,
    hit_lists: &[Option<BlastnHitList>],
    writer: &mut W,
    query_ids: &[Arc<[u8]>],
    subject_ids: &[Arc<[u8]>],
    output_format: BlastnOutputFormat,
    query_titles: &[Arc<[u8]>],
    subject_title: &str,
    unsearched: &[bool],
    epilog: Option<usize>,
    probe: Option<&mut FormatProbe<'_>>,
    warnings: Option<&mut QueryWarnings<'_>>,
) -> io::Result<()> {
    match pool {
        Some(pool) if probe.is_none() && !super::tracing::enabled() => {
            crate::utils::xreach::x_reach!("LOSAT_X_BNTABPAR", true);
            write_on_pool(
                pool,
                hit_lists,
                writer,
                query_ids,
                subject_ids,
                output_format,
                query_titles,
                subject_title,
                unsearched,
                epilog,
                warnings,
            )
        }
        _ => {
            if mode() {
                crate::utils::xreach::x_reach!("LOSAT_X_BNTABPAR", false);
            }
            write_output_blastn_hitlists_to_writer(
                hit_lists,
                writer,
                query_ids,
                subject_ids,
                output_format,
                query_titles,
                subject_title,
                unsearched,
                epilog,
                probe,
                warnings,
            )
        }
    }
}

/// The rows of `hit_list` (its HSPs; 0 without a list).
fn rows_of(hit_list: Option<&BlastnHitList>) -> usize {
    hit_list.map_or(0, |list| {
        list.hsplist_array.iter().map(|hsps| hsps.hsps.len()).sum()
    })
}

/// The serial writer's bytes, with each window's queries formatted on `pool` (module docs).
#[allow(clippy::too_many_arguments)]
fn write_on_pool<W: Write>(
    pool: &SearchPool<'_>,
    hit_lists: &[Option<BlastnHitList>],
    writer: &mut W,
    query_ids: &[Arc<[u8]>],
    subject_ids: &[Arc<[u8]>],
    output_format: BlastnOutputFormat,
    query_titles: &[Arc<[u8]>],
    subject_title: &str,
    unsearched: &[bool],
    epilog: Option<usize>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
) -> io::Result<()> {
    let config = OutputConfig::ncbi_compat();
    let window_rows = ROWS_PER_THREAD.saturating_mul(pool.threads().max(1));
    let format = |q_idx: usize| -> io::Result<Vec<u8>> {
        let hit_list = hit_lists[q_idx].as_ref();
        let mut buffer = Vec::with_capacity(rows_of(hit_list) * ROW_BYTES);
        // The running HSP index is the format observer's, which this path never has.
        let mut hsp_index = 0;
        write_query_rows(
            &mut buffer,
            q_idx,
            hit_list,
            query_ids,
            subject_ids,
            output_format,
            query_titles,
            subject_title,
            unsearched,
            &config,
            &mut hsp_index,
            None,
        )?;
        Ok(buffer)
    };
    let mut start = 0;
    while start < hit_lists.len() {
        let mut end = start;
        let mut rows = 0;
        while end < hit_lists.len() && (end == start || rows < window_rows) {
            rows += rows_of(hit_lists[end].as_ref());
            end += 1;
        }
        let buffers = map(pool, start..end, &format);
        for (q_idx, buffer) in (start..end).zip(buffers) {
            if let Some(warnings) = warnings.as_deref_mut() {
                warnings.before_query(q_idx, &mut *writer)?;
            }
            writer.write_all(&buffer?)?;
        }
        start = end;
    }
    write_tabular_epilog(writer, output_format, epilog)
}

/// `f` of every query of `queries`, in order, on `pool`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn map<R: Send>(
    pool: &SearchPool<'_>,
    queries: std::ops::Range<usize>,
    f: &(impl Fn(usize) -> R + Sync),
) -> Vec<R> {
    pool.install(|| queries.into_par_iter().map(f).collect())
}

/// `f` of every query of `queries`, in order (builds without threads, where `pool` gives `None`
/// and this is never called).
#[cfg(not(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
)))]
fn map<R: Send>(
    _pool: &SearchPool<'_>,
    queries: std::ops::Range<usize>,
    f: &(impl Fn(usize) -> R + Sync),
) -> Vec<R> {
    queries.map(f).collect()
}

#[cfg(all(
    test,
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod tests {
    use super::super::hsp::{BlastnHsp, BlastnHspList};
    use super::*;

    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:125-148
    // ```c
    // typedef struct BlastHSP {
    //    Int4 score;
    //    double evalue;
    //    BlastSeg query;
    //    BlastSeg subject;
    //    Int4 context;
    // } BlastHSP;
    // ```
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n.max(1)
        }
    }

    fn hsp(rng: &mut Rng, q_idx: u32, s_idx: u32) -> BlastnHsp {
        let length = 1 + rng.below(3000) as usize;
        let num_ident = rng.below(length as u64 + 1) as usize;
        let e_value = match rng.below(5) {
            0 => 0.0,
            1 => 1e-180 * (1 + rng.below(9)) as f64,
            2 => (rng.below(1000) as f64) * 1e-5,
            3 => 1.0 + rng.below(9) as f64 / 7.0,
            _ => 10f64.powi(-(rng.below(300) as i32)),
        };
        let q_start = 1 + rng.below(5000) as usize;
        let s_start = 1 + rng.below(3_000_000) as usize;
        let minus = rng.below(2) == 1;
        BlastnHsp {
            identity: 100.0 * num_ident as f64 / length as f64,
            length,
            mismatch: length - num_ident,
            gapopen: rng.below(9) as usize,
            q_start,
            q_end: q_start + length - 1,
            s_start,
            s_end: if minus {
                s_start.saturating_sub(length - 1).max(1)
            } else {
                s_start + length - 1
            },
            e_value,
            bit_score: rng.below(6000) as f64 / 3.0 + 0.5,
            num_ident,
            query_frame: if minus { -1 } else { 1 },
            query_length: q_start + length,
            q_idx,
            s_idx,
            raw_score: rng.below(3000) as i32,
            internal_q_offset_0: q_start - 1,
            internal_q_end_0: q_start + length - 1,
            internal_s_offset_0: s_start - 1,
            internal_s_end_0: s_start + length - 1,
            internal_query_context_offset: 0,
            gap_info: None,
            num_positives: num_ident,
        }
    }

    /// Hit lists of `queries` queries: some without a list, some with empty HSP lists, the others
    /// with up to `max_rows` rows over up to three subjects.
    fn hit_lists(rng: &mut Rng, queries: usize, max_rows: u64) -> Vec<Option<BlastnHitList>> {
        (0..queries as u32)
            .map(|q_idx| {
                if rng.below(4) == 0 {
                    return None;
                }
                let subjects = rng.below(4) as u32;
                let mut list = BlastnHitList::new(subjects.max(1) as usize);
                for s_idx in 0..subjects {
                    let rows = rng.below(max_rows + 1);
                    let hsps = (0..rows).map(|_| hsp(rng, q_idx, s_idx)).collect();
                    list.hsplist_array.push(BlastnHspList {
                        oid: s_idx,
                        query_index: q_idx,
                        hsps,
                        best_evalue: 0.0,
                    });
                }
                list.hsplist_count = list.hsplist_array.len();
                Some(list)
            })
            .collect()
    }

    fn ids(prefix: &str, n: usize) -> Vec<Arc<[u8]>> {
        (0..n)
            .map(|i| Arc::<[u8]>::from(format!("{prefix}{i}").into_bytes()))
            .collect()
    }

    struct Case {
        hit_lists: Vec<Option<BlastnHitList>>,
        query_ids: Vec<Arc<[u8]>>,
        subject_ids: Vec<Arc<[u8]>>,
        query_titles: Vec<Arc<[u8]>>,
        unsearched: Vec<bool>,
        warnings: Vec<Vec<u8>>,
    }

    fn case(rng: &mut Rng, queries: usize, max_rows: u64) -> Case {
        Case {
            hit_lists: hit_lists(rng, queries, max_rows),
            query_ids: ids("query_", queries),
            subject_ids: ids("lcl|subject_", 3),
            query_titles: ids("query title ", queries),
            unsearched: (0..queries).map(|_| rng.below(9) == 0).collect(),
            warnings: (0..queries)
                .map(|q| {
                    if rng.below(5) == 0 {
                        format!("Warning: [blastn] Query_{q} warning\n").into_bytes()
                    } else {
                        Vec::new()
                    }
                })
                .collect(),
        }
    }

    /// The report and the warnings' sink of the serial writer (`pool` `None`) or the pool's.
    fn run(
        case: &Case,
        format: BlastnOutputFormat,
        epilog: Option<usize>,
        pool: Option<&SearchPool<'_>>,
    ) -> (Vec<u8>, Vec<u8>) {
        let mut report = Vec::new();
        let mut sink = Vec::new();
        let mut warnings = QueryWarnings {
            before: &case.warnings,
            sink: &mut sink,
        };
        let args = (
            &case.hit_lists[..],
            &case.query_ids[..],
            &case.subject_ids[..],
            &case.query_titles[..],
            &case.unsearched[..],
        );
        match pool {
            Some(pool) => write_on_pool(
                pool,
                args.0,
                &mut report,
                args.1,
                args.2,
                format,
                args.3,
                "User specified sequence set (Input: subject.fa)",
                args.4,
                epilog,
                Some(&mut warnings),
            ),
            None => write_output_blastn_hitlists_to_writer(
                args.0,
                &mut report,
                args.1,
                args.2,
                format,
                args.3,
                "User specified sequence set (Input: subject.fa)",
                args.4,
                epilog,
                None,
                Some(&mut warnings),
            ),
        }
        .unwrap();
        (report, sink)
    }

    #[test]
    fn rows_on_the_pool_are_the_serial_writers_bytes() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        let cases: Vec<Case> = [
            (0usize, 0u64),
            (1, 0),
            (1, 5),
            (7, 3),
            (60, 40),
            (400, 9),
            (25, 3000),
        ]
        .into_iter()
        .map(|(queries, max_rows)| case(&mut rng, queries, max_rows))
        .collect();
        for threads in [2usize, 3, 8] {
            crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
                for case in &cases {
                    let queries = case.hit_lists.len();
                    for format in [
                        BlastnOutputFormat::Tabular,
                        BlastnOutputFormat::TabularWithComments,
                    ] {
                        for epilog in [None, Some(queries + 3)] {
                            let serial = run(case, format, epilog, None);
                            let pooled = run(case, format, epilog, Some(pool));
                            assert_eq!(
                                serial, pooled,
                                "threads {threads} queries {queries} {format:?} {epilog:?}"
                            );
                        }
                    }
                }
                Ok(())
            })
            .unwrap();
        }
    }

    #[test]
    fn a_query_larger_than_a_window_and_many_small_ones() {
        // The first query has more rows than a window of two threads, the others a few each.
        let mut rng = Rng(7);
        let mut case = case(&mut rng, 50, 4);
        let mut big = BlastnHitList::new(1);
        big.hsplist_array.push(BlastnHspList {
            oid: 1,
            query_index: 0,
            hsps: (0..(ROWS_PER_THREAD * 2 + 17))
                .map(|_| hsp(&mut rng, 0, 1))
                .collect(),
            best_evalue: 0.0,
        });
        big.hsplist_count = 1;
        case.hit_lists[0] = Some(big);
        let serial = run(&case, BlastnOutputFormat::Tabular, None, None);
        crate::utils::threading::with_search_pool(2, "blastn", |pool| {
            assert_eq!(
                serial,
                run(&case, BlastnOutputFormat::Tabular, None, Some(pool))
            );
            Ok(())
        })
        .unwrap();
        assert!(serial.0.len() > ROWS_PER_THREAD * 2 * 40);
    }
}
