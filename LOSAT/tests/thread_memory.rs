//! The shape of BLASTN's memory over threads: what a `-subject` search of four threads holds
//! beyond the same search of one thread (the extra) must not grow with the input, and must stay
//! below a constant per added thread.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3224-3236
//! ```c++
//!         // This is temporarily ignored (per SB-635)
//!         if (args.Exist(kArgSubject) && args[kArgSubject].HasValue() &&
//!             m_NumThreads != CThreadable::kMinNumThreads) {
//!             m_NumThreads = CThreadable::kMinNumThreads;
//! ```
//! No NCBI counterpart for the property: NCBI searches a `-subject` on one thread whatever
//! `-num_threads` says, and LOSAT honors `-num_threads` there (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES`),
//! so the extra is LOSAT's own. A workspace made per piece of work instead of per thread, kept at the
//! size of the largest piece it has seen, or holding a whole stage's results gives an extra that
//! grows with the input; workspaces of a constant size per thread do not.
//!
//! What is measured. The global allocator of this binary counts the bytes of the live blocks and
//! their most at once (`HELD`); a search runs through the library API (`LOSAT::api::run`, outfmt 6)
//! on a thread of its own, and its peak is the most bytes held during the call beyond those held
//! before it. On 64-bit Linux with glibc LOSAT asks for 32 MiB and one page, rounded up to a whole
//! element, for a buffer that it wants glibc to map (`utils/x_mapped.rs` `X_MAPPED_BYTES`,
//! `algorithm/blastn/x_mb_lookup_par.rs` `lines_capacity`) and never touches the pages beyond what
//! it uses. A block of that size (`RESERVE`) is counted apart (`RESERVED`), not in the peak: the
//! allocator cannot see how much of it is used. On dev at these inputs that leaves out the lookup
//! fill's arena of a search of two threads or more and, at 64 units, the ungapped hits of the
//! one-thread search (1.4 MiB used, measured with `LOSAT_X_BNMAPPED=0`), so the 64-unit blastn extra
//! below is about 1.4 MiB above what the two searches use.
//!
//! The inputs (`inputs`), one query and one subject (the pair path of one chunk), are units of
//! 25,000 random letters each (fixed seed): in each unit a 6,000-letter piece of the query is planted
//! in the subject with 6 % substitutions and 0.6 % indels (a long gapped alignment), eight pieces of
//! 150-599 letters with 4-15 % substitutions (every second one reverse-complemented), and one
//! 300-letter element is planted six times in each sequence with 3 % substitutions (dense seeds and
//! 36 HSPs). A larger input has more units of the same kind, never longer alignments, so a
//! workspace of a constant size per thread gives the same extra at every size. Sizes 4, 16 and 64
//! units (100 kb, 400 kb, 1.6 Mb): blastn grows mainly from 16 to 64 units, dc-megablast from 4 to
//! 16, so one step of four would miss one of them.
//!
//! The checks (`check`), on 64-bit x86 Linux with glibc, where the table below was measured
//! (`MEASURED_HERE`; elsewhere the checks print the table only): the extra of the next size is at
//! most `GROWTH` times that of the size before plus `SLACK` (it does not grow with the input), and
//! at every size at most `PER_THREAD` per added thread. Everywhere: the four-thread reports equal the
//! one-thread report byte for byte and have at least `ROWS_PER_UNIT` rows per unit, so the inputs
//! reach the parallel stages and the searches find what was planted. The one-thread peak was the
//! same to the byte in every repetition measured (one search allocates the same blocks in the same
//! order); the four-thread peak depends on which thread takes which piece of work, so it is the
//! median of `REPEATS` searches. A search of each task at the first size runs first, unmeasured,
//! so that the process's one-time tables are not in a measured peak.
//!
//! `TASKS` holds the shape of each task on dev (`Grows`: the extra grows with the input). A change
//! that gives a task's workspaces a constant size per thread turns its row to `Flat`, and from then
//! on the row guards it; the test fails when a row no longer holds, either way.
//!
//! Measured on dev 8a5656f7f (test profile, default features; peaks in KiB; extra = T4 - T1;
//! growth = this extra / the previous size's):
//!
//! | task | units | T1 | T4 (median of 3) | extra | extra / 3 | growth |
//! |---|---|---|---|---|---|---|
//! | blastn | 4 | 22,237 | 22,763 | 525 | 175 | - |
//! | blastn | 16 | 23,348 | 24,324 | 976 | 325 | 1.86 |
//! | blastn | 64 | 35,455 | 40,863 | 5,407 | 1,802 | 5.54 |
//! | dc-megablast | 4 | 17,884 | 19,663 | 1,778 | 592 | - |
//! | dc-megablast | 16 | 21,293 | 25,736 | 4,443 | 1,481 | 2.50 |
//! | dc-megablast | 64 | 34,746 | 39,971 | 5,225 | 1,741 | 1.18 |
//! | megablast | 4 | 22,298 | 22,554 | 255 | 85 | - |
//! | megablast | 16 | 53,584 | 53,919 | 335 | 111 | 1.31 |
//! | megablast | 64 | 65,740 | 66,300 | 560 | 186 | 1.67 |
//!
//! blastn grows (5,407 > 1.25 x 976 + 1,024) and dc-megablast grows (4,443 > 1.25 x 1,778 +
//! 1,024); megablast is flat. With `--all-features` (as CI) and with `--release` the verdicts are the
//! same (blastn at 64 units: extra 6,253 and 6,049 KiB; the four-thread median moves by up to 0.9 MiB
//! between runs). No size is above 2.1 MiB per added thread. The binary takes 4 s on an idle
//! workstation of 32 CPUs and up to 11 s under load.
//!
//! The file is built with the `parallel` feature only: without it `-num_threads 4` is rejected
//! (`threading.rs:121-130`).

#![cfg(feature = "parallel")]

use std::alloc::{GlobalAlloc, Layout, System};
use std::io::Write;
use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use LOSAT::api::{read_fasta, run, validate, FastaRecord, Inputs, Program, Role};

// ---------------------------------------------------------------------------------------------
// The counting allocator
// ---------------------------------------------------------------------------------------------

/// No NCBI counterpart: bytes of live blocks and the most of them at once since `start`.
struct Counter {
    live: AtomicUsize,
    peak: AtomicUsize,
}

impl Counter {
    const fn new() -> Self {
        Self {
            live: AtomicUsize::new(0),
            peak: AtomicUsize::new(0),
        }
    }

    fn add(&self, bytes: usize) {
        let live = self.live.fetch_add(bytes, Ordering::Relaxed) + bytes;
        self.peak.fetch_max(live, Ordering::Relaxed);
    }

    fn sub(&self, bytes: usize) {
        self.live.fetch_sub(bytes, Ordering::Relaxed);
    }

    /// Starts a window: the peak becomes the bytes live now, which it returns.
    fn start(&self) -> usize {
        let live = self.live.load(Ordering::Relaxed);
        self.peak.store(live, Ordering::Relaxed);
        live
    }

    /// The most bytes live at once since `start` returned `live`, beyond `live`.
    fn peak_since(&self, live: usize) -> usize {
        self.peak.load(Ordering::Relaxed).saturating_sub(live)
    }
}

/// No NCBI counterpart: every block outside `RESERVE`.
static HELD: Counter = Counter::new();
/// No NCBI counterpart: the blocks of `RESERVE`.
static RESERVED: Counter = Counter::new();

/// No NCBI counterpart: the sizes of the blocks that LOSAT asks glibc to map on 64-bit Linux:
/// `(32 << 20) + 4096` bytes (`X_MAPPED_BYTES`, `lines_capacity`) rounded up to a whole element of
/// at most 64 bytes (BLASTN's elements: letters, ungapped hits of 44 bytes, lines of 64). A copy of
/// private constants of the library, which this file cannot read: `X_MAPPED_BYTES`
/// (`src/utils/x_mapped.rs:37`) and `lines_capacity` (`src/algorithm/blastn/x_mb_lookup_par.rs:403-405`).
const RESERVE: Range<usize> = (32 << 20) + 4096..(32 << 20) + 4096 + 64;

fn counter(bytes: usize) -> &'static Counter {
    if RESERVE.contains(&bytes) {
        &RESERVED
    } else {
        &HELD
    }
}

/// No NCBI counterpart: the system allocator, counting (`HELD`, `RESERVED`).
struct Counting;

// SAFETY: every call goes to `System` with the caller's arguments unchanged; the counters are
// atomics and allocate nothing.
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        // SAFETY: the caller's contract for `alloc` is `System.alloc`'s.
        let block = unsafe { System.alloc(layout) };
        if !block.is_null() {
            counter(layout.size()).add(layout.size());
        }
        block
    }

    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        // SAFETY: as `alloc`.
        let block = unsafe { System.alloc_zeroed(layout) };
        if !block.is_null() {
            counter(layout.size()).add(layout.size());
        }
        block
    }

    unsafe fn dealloc(&self, block: *mut u8, layout: Layout) {
        // SAFETY: `block` came from this allocator, so from `System`, with `layout`.
        unsafe { System.dealloc(block, layout) };
        counter(layout.size()).sub(layout.size());
    }

    unsafe fn realloc(&self, block: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        // SAFETY: as `dealloc`, and `new_size` is the caller's.
        let moved = unsafe { System.realloc(block, layout, new_size) };
        if !moved.is_null() {
            counter(layout.size()).sub(layout.size());
            counter(new_size).add(new_size);
        }
        moved
    }
}

#[global_allocator]
static COUNTING: Counting = Counting;

/// No NCBI counterpart: the tests of this binary measure one at a time (`cargo test` runs them on
/// threads of one process; nextest runs each in a process of its own).
static ONE_AT_A_TIME: Mutex<()> = Mutex::new(());

// ---------------------------------------------------------------------------------------------
// The inputs
// ---------------------------------------------------------------------------------------------

/// No NCBI counterpart: splitmix64, so that the inputs are the same on every machine.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn base(&mut self) -> u8 {
        b"ACGT"[(self.next() >> 62) as usize]
    }

    fn letters(&mut self, len: usize) -> Vec<u8> {
        (0..len).map(|_| self.base()).collect()
    }
}

/// A copy of `source` with `substitutions` and `indels` per thousand letters (half deletions,
/// half insertions).
fn mutated(rng: &mut Rng, source: &[u8], substitutions: usize, indels: usize) -> Vec<u8> {
    let mut copy = Vec::with_capacity(source.len() + source.len() / 50);
    for &letter in source {
        let roll = rng.below(1000);
        if roll < indels / 2 {
            continue;
        }
        if roll < indels {
            copy.push(rng.base());
            copy.push(letter);
        } else if roll < indels + substitutions {
            let at = b"ACGT".iter().position(|&b| b == letter).unwrap_or(0);
            copy.push(b"ACGT"[(at + 1 + rng.below(3)) % 4]);
        } else {
            copy.push(letter);
        }
    }
    copy
}

fn reverse_complement(source: &[u8]) -> Vec<u8> {
    source
        .iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

/// Writes `piece` over `target` from `at`, as far as `target` goes.
fn plant(target: &mut [u8], at: usize, piece: &[u8]) {
    let end = (at + piece.len()).min(target.len());
    target[at..end].copy_from_slice(&piece[..end - at]);
}

/// The letters of a unit of the query and of the subject.
const UNIT: usize = 25_000;
/// The query piece planted in the subject whole (the long alignment).
const LONG: usize = 6_000;

/// One unit (module header), appended to `query` and `subject`.
fn unit(rng: &mut Rng, query: &mut Vec<u8>, subject: &mut Vec<u8>) {
    let mut q = rng.letters(UNIT);
    let mut s = rng.letters(UNIT);
    let from = rng.below(UNIT - LONG);
    let long = mutated(rng, &q[from..from + LONG], 60, 6);
    plant(&mut s, rng.below(UNIT - LONG), &long);
    for k in 0..8 {
        let len = 150 + rng.below(450);
        let from = rng.below(UNIT - len);
        let divergence = 40 + rng.below(110);
        let mut piece = mutated(rng, &q[from..from + len], divergence, 4);
        if k % 2 == 1 {
            piece = reverse_complement(&piece);
        }
        plant(&mut s, rng.below(UNIT - len), &piece);
    }
    let element = rng.letters(300);
    for _ in 0..6 {
        let copy = mutated(rng, &element, 30, 2);
        plant(&mut q, rng.below(UNIT - 300), &copy);
        let copy = mutated(rng, &element, 30, 2);
        plant(&mut s, rng.below(UNIT - 300), &copy);
    }
    query.extend_from_slice(&q);
    subject.extend_from_slice(&s);
}

fn fasta(name: &str, letters: &[u8]) -> Vec<u8> {
    let mut text = format!(">{name}\n").into_bytes();
    for line in letters.chunks(80) {
        text.extend_from_slice(line);
        text.push(b'\n');
    }
    text
}

/// The query and subject records of `units` units; a larger input starts with a smaller one.
fn inputs(units: usize) -> (Vec<FastaRecord>, Vec<FastaRecord>) {
    let mut rng = Rng(0x005E_ED0F_6E0E);
    let (mut query, mut subject) = (Vec::new(), Vec::new());
    for _ in 0..units {
        unit(&mut rng, &mut query, &mut subject);
    }
    (
        read_fasta(Program::Blastn, Role::Query, &fasta("q", &query)).expect("query records"),
        read_fasta(Program::Blastn, Role::Subject, &fasta("s", &subject)).expect("subject records"),
    )
}

// ---------------------------------------------------------------------------------------------
// The searches and the checks
// ---------------------------------------------------------------------------------------------

/// The sizes, in units.
const SIZES: [usize; 3] = [4, 16, 64];
/// Four-thread searches per size (the median is taken).
const REPEATS: usize = 3;
/// The extra of the next size may be this many times that of the size before, plus `SLACK`.
const GROWTH: f64 = 1.25;
const SLACK: i64 = 1 << 20;
/// The extra per added thread, at every size.
const PER_THREAD: i64 = 4 << 20;
/// The rows of a report, at least, per unit (each unit plants about 45 HSPs).
const ROWS_PER_UNIT: usize = 20;
/// Where `TASKS` was measured: the checks of the shape run only there.
const MEASURED_HERE: bool = cfg!(all(
    target_os = "linux",
    target_env = "gnu",
    target_arch = "x86_64"
));

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Shape {
    /// The extra stays within `GROWTH` and `SLACK` from each size to the next.
    Flat,
    /// The extra grows with the input.
    Grows,
}

/// The shape of each task on dev (module header). A `Grows` row guards only the cap per added
/// thread (`PER_THREAD`): the test passes while the extra grows. dc-megablast's `Grows` rests on
/// the step from 4 to 16 units (4,443 KiB against a threshold of about 3,240 KiB), so a later
/// change can turn it either way. A change of memory per thread should read the printed table, not
/// only pass or fail.
const TASKS: [(&str, Shape); 3] = [
    ("blastn", Shape::Grows),
    ("dc-megablast", Shape::Grows),
    ("megablast", Shape::Flat),
];

/// One search: its peak (`HELD`), the reserve blocks' peak (`RESERVED`) and its report.
struct Searched {
    peak: usize,
    reserve: usize,
    report: Vec<u8>,
}

fn search(task: &str, threads: usize, query: &[FastaRecord], subject: &[FastaRecord]) -> Searched {
    let threads = threads.to_string();
    let search =
        validate(Program::Blastn, &["-task", task, "-num_threads", &threads]).expect("the options");
    // On a new thread, so that every search starts without the thread-local caches of the last.
    std::thread::scope(|scope| {
        scope
            .spawn(|| {
                let mut report: Vec<u8> = Vec::new();
                let mut diagnostics: Vec<u8> = Vec::new();
                let (held, reserved) = (HELD.start(), RESERVED.start());
                {
                    let mut reports: [(&str, &mut (dyn Write + Send)); 1] = [("6", &mut report)];
                    run(
                        &search,
                        Inputs::new(query, subject),
                        &mut reports,
                        &mut diagnostics,
                    )
                    .expect("the search");
                }
                let (peak, reserve) = (HELD.peak_since(held), RESERVED.peak_since(reserved));
                assert!(
                    diagnostics.is_empty(),
                    "{}",
                    String::from_utf8_lossy(&diagnostics)
                );
                Searched {
                    peak,
                    reserve,
                    report,
                }
            })
            .join()
            .expect("the search thread")
    })
}

/// One size: the one-thread peak and the median four-thread peak, and the reserves'.
struct Row {
    units: usize,
    one: usize,
    four: usize,
    one_reserve: usize,
    four_reserve: usize,
}

impl Row {
    fn extra(&self) -> i64 {
        self.four as i64 - self.one as i64
    }
}

fn kib(bytes: i64) -> i64 {
    bytes / 1024
}

fn table(task: &str, rows: &[Row]) -> String {
    let mut text = format!(
        "{task}: peaks in KiB (reserve blocks apart), extra = T4 - T1\n\
         units  T1  T4(median of {REPEATS})  extra  extra/3  growth  reserve(T1, T4)\n"
    );
    for (i, row) in rows.iter().enumerate() {
        let growth = match i {
            0 => "-".to_string(),
            _ => format!("{:.2}", row.extra() as f64 / rows[i - 1].extra() as f64),
        };
        text.push_str(&format!(
            "{}  {}  {}  {}  {}  {growth}  {}, {}\n",
            row.units,
            kib(row.one as i64),
            kib(row.four as i64),
            kib(row.extra()),
            kib(row.extra() / 3),
            kib(row.one_reserve as i64),
            kib(row.four_reserve as i64),
        ));
    }
    text
}

fn check(task: &str, expected: Shape) {
    let _one = ONE_AT_A_TIME
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    let (query, subject) = inputs(SIZES[0]);
    search(task, 1, &query, &subject);
    let mut rows = Vec::new();
    for units in SIZES {
        let (query, subject) = inputs(units);
        let one = search(task, 1, &query, &subject);
        let rows_found = one.report.iter().filter(|&&b| b == b'\n').count();
        assert!(
            rows_found >= ROWS_PER_UNIT * units,
            "{task}, {units} units: {rows_found} rows, fewer than {ROWS_PER_UNIT} per unit"
        );
        let mut fours: Vec<Searched> = (0..REPEATS)
            .map(|_| search(task, 4, &query, &subject))
            .collect();
        for four in &fours {
            assert!(
                four.report == one.report,
                "{task}, {units} units: the four-thread report differs from the one-thread report"
            );
        }
        fours.sort_by_key(|four| four.peak);
        let median = &fours[REPEATS / 2];
        rows.push(Row {
            units,
            one: one.peak,
            four: median.peak,
            one_reserve: one.reserve,
            four_reserve: median.reserve,
        });
    }
    let table = table(task, &rows);
    println!("{table}");
    if !MEASURED_HERE {
        return;
    }
    for row in &rows {
        assert!(
            row.extra() <= 3 * PER_THREAD,
            "{task}, {} units: {} KiB per added thread, more than {} KiB\n\
             If a padded block's size or count changed (`X_MAPPED_BYTES`, `lines_capacity`), check \
             `RESERVE` first: a block it no longer matches is counted in the peak.\n{table}",
            row.units,
            kib(row.extra() / 3),
            kib(PER_THREAD)
        );
    }
    let grows = rows
        .windows(2)
        .any(|pair| pair[1].extra() as f64 > GROWTH * pair[0].extra() as f64 + SLACK as f64);
    let measured = if grows { Shape::Grows } else { Shape::Flat };
    match (expected, measured) {
        (Shape::Flat, Shape::Grows) => panic!(
            "{task}: a regression of memory per thread: the extra, flat on dev, now grows with the \
             input. Do not edit the row of TASKS; find the change that made the extra grow.\n{table}"
        ),
        (Shape::Grows, Shape::Flat) => panic!(
            "{task}: the extra no longer grows with the input. Set the task's row in TASKS to \
             `Flat`, and update the table and the prose of the module header with the new \
             numbers.\n{table}"
        ),
        _ => {}
    }
}

fn expected(task: &str) -> Shape {
    TASKS
        .iter()
        .find(|(name, _)| *name == task)
        .map(|(_, shape)| *shape)
        .expect("a task of TASKS")
}

#[test]
fn blastn_extra_per_thread() {
    check("blastn", expected("blastn"));
}

#[test]
fn dc_megablast_extra_per_thread() {
    check("dc-megablast", expected("dc-megablast"));
}

#[test]
fn megablast_extra_per_thread() {
    check("megablast", expected("megablast"));
}
