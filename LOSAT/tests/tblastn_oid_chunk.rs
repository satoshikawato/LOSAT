//! How TBLASTN hands its subjects to the search pool (`LOSAT_X_TBNOIDCHUNK`,
//! `src/algorithm/tblastn/stage_d_pipeline.rs`), on whole searches of many subjects.
//!
//! NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:697-707
//! ```c++
//!         Int8 tot_length = m_Atlas.GetSliceSize() / (4*m_NumThreads) + 1;
//!
//!         res.length = vol->GetSequence(vol_oid++, &seq);
//!         if (res.length < 0) return;
//!         // must return at least one sequence
//!         do {
//!             tot_length -= res.length;
//!             res.address = seq;
//!             buffer->results.push_back(res);
//!             res.length = vol->GetSequence(vol_oid++, &seq);
//!         } while (res.length >= 0 && tot_length >= res.length && vol_oid < m_RestrictEnd);
//! ```
//! With two or more threads each NCBI search thread takes a chunk of consecutive subjects whose
//! bases stay within `min(2^30, F) / (4 * threads) + 1`, at least one subject; F is the largest
//! file of the database, and for `-subject` LOSAT takes ceil(bases / 4), the `.nsq` that
//! `makeblastdb` writes for the same subjects. Without the switch one entry into LOSAT's pool
//! takes one subject per thread; with it, four such chunks per thread. The switch is in the strict
//! set (`src/utils/xswitch.rs`): it is on by default, off with `LOSAT_X_TBNOIDCHUNK=0`, and off
//! with `LOSAT_X_DEFAULT=0` unless its variable is set. It must change the entries and nothing in
//! the report, so the tests compare the report with the one-thread search, which has no pool and
//! searches the subjects in NCBI's loop order, and read the `pool_entries` field of the
//! `[losat-thread-activity]` line (`LOSAT_WASI_THREADS_DEBUG=1`) to check which entries ran.
//!
//! The last test checks `LOSAT_X_TBNSTATE` (`src/algorithm/tblastn/search_gapped.rs`, in the
//! strict set), the scan state of the unit path kept by a thread from one subject to the next, on
//! the same searches.
//!
//! The fixture suites have no search of many subjects, so the tests cut
//! `tests/fasta/MelaMJNV.fasta` into subjects of 1,000 bases (240 of them, or 200 with one subject
//! of 48,000 bases first, in the middle or last) and search them with the first ten proteins of
//! `tests/fasta/PemoMJNVA.faa`.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

const SUBJECTS: usize = 240;
const SUBJECT_LENGTH: usize = 1_000;
const QUERIES: usize = 10;
const THREADS: usize = 4;
/// The sets with one long subject: `SHORT_SUBJECTS` pieces of `SUBJECT_LENGTH` bases and one of
/// `LONG_SUBJECT` bases, which is longer than a chunk's budget at every thread count tested.
const SHORT_SUBJECTS: usize = 200;
const LONG_SUBJECT: usize = 48_000;

/// The records of a FASTA file as (title line, letters).
fn records(path: &Path) -> Vec<(String, String)> {
    let text = std::fs::read_to_string(path).expect("fixture");
    let mut out: Vec<(String, String)> = Vec::new();
    for line in text.lines() {
        match line.strip_prefix('>') {
            Some(title) => out.push((title.to_string(), String::new())),
            None => out
                .last_mut()
                .expect("title line first")
                .1
                .push_str(line.trim()),
        }
    }
    out
}

/// The query and subject files of the tests, written once: the first `QUERIES` proteins, and
/// the genome cut into `SUBJECTS` pieces of `SUBJECT_LENGTH` bases.
fn inputs() -> &'static (PathBuf, PathBuf) {
    static FILES: OnceLock<(PathBuf, PathBuf)> = OnceLock::new();
    FILES.get_or_init(|| {
        let fasta = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
        let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("tblastn_oid_chunk");
        std::fs::create_dir_all(&dir).expect("work directory");
        let mut query = String::new();
        for (title, letters) in records(&fasta.join("PemoMJNVA.faa")).iter().take(QUERIES) {
            query.push_str(&format!(">{title}\n{letters}\n"));
        }
        let genome = &records(&fasta.join("MelaMJNV.fasta"))[0].1;
        assert!(
            genome.len() >= SUBJECTS * SUBJECT_LENGTH,
            "genome too short"
        );
        let mut subjects = String::new();
        for index in 0..SUBJECTS {
            let piece = &genome[index * SUBJECT_LENGTH..(index + 1) * SUBJECT_LENGTH];
            subjects.push_str(&format!(">piece{index}\n{piece}\n"));
        }
        let files = (dir.join("query.faa"), dir.join("subjects.fna"));
        std::fs::write(&files.0, query).expect("query file");
        std::fs::write(&files.1, subjects).expect("subject file");
        files
    })
}

/// Where the long subject stands in its set.
#[derive(Clone, Copy, Debug)]
enum Long {
    First,
    Middle,
    Last,
}

/// The subject file with one long subject (`LONG_SUBJECT` bases, the start of the genome) and
/// `SHORT_SUBJECTS` pieces of `SUBJECT_LENGTH` bases after it in the genome; the long subject is
/// first, after half of the pieces, or last. The three files are written once for the whole test
/// binary: the tests run side by side, and a test that wrote a file again while another one's
/// search read it gave that search an empty subject file.
fn long_subject_set(position: Long) -> PathBuf {
    static FILES: OnceLock<Vec<PathBuf>> = OnceLock::new();
    FILES.get_or_init(|| {
        [Long::First, Long::Middle, Long::Last]
            .into_iter()
            .map(write_long_subject_set)
            .collect()
    })[position as usize]
        .clone()
}

/// Writes the subject file of `long_subject_set`.
fn write_long_subject_set(position: Long) -> PathBuf {
    let fasta = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fasta");
    let genome = &records(&fasta.join("MelaMJNV.fasta"))[0].1;
    assert!(
        genome.len() >= LONG_SUBJECT + SHORT_SUBJECTS * SUBJECT_LENGTH,
        "genome too short"
    );
    let long = format!(">long\n{}\n", &genome[..LONG_SUBJECT]);
    let at = match position {
        Long::First => 0,
        Long::Middle => SHORT_SUBJECTS / 2,
        Long::Last => SHORT_SUBJECTS,
    };
    let mut subjects = String::new();
    for index in 0..SHORT_SUBJECTS {
        if index == at {
            subjects.push_str(&long);
        }
        let start = LONG_SUBJECT + index * SUBJECT_LENGTH;
        let piece = &genome[start..start + SUBJECT_LENGTH];
        subjects.push_str(&format!(">piece{index}\n{piece}\n"));
    }
    if at == SHORT_SUBJECTS {
        subjects.push_str(&long);
    }
    let path = inputs()
        .1
        .with_file_name(format!("long_{position:?}.fna").to_lowercase());
    std::fs::write(&path, subjects).expect("subject file");
    path
}

/// One TBLASTN run of the queries against `subjects` (`-outfmt 6` on standard output) with the
/// given variables as its only `LOSAT_*` variables, except `LOSAT_BLASTX_WORKER_LOG`: a build
/// with the worker probe (`--all-features`) writes to that file as soon as a pool starts, so the
/// caller's value is kept, or a file of the tests is named. Returns the report and the standard
/// error.
fn tblastn_on(subjects: &Path, threads: usize, variables: &[(&str, &str)]) -> (Vec<u8>, String) {
    const WORKER_LOG: &str = "LOSAT_BLASTX_WORKER_LOG";
    let query = &inputs().0;
    let mut command = Command::new(env!("CARGO_BIN_EXE_LOSAT"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if name_text.starts_with("LOSAT_") && name_text != WORKER_LOG {
            command.env_remove(name);
        }
    }
    if std::env::var_os(WORKER_LOG).is_none() {
        command.env(WORKER_LOG, query.with_file_name("worker.log"));
    }
    let output = command
        .arg("tblastn")
        .arg("-query")
        .arg(query)
        .arg("-subject")
        .arg(subjects)
        .args(["-outfmt", "6", "-num_threads", &threads.to_string()])
        .envs(variables.iter().copied())
        .output()
        .expect("run LOSAT");
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(output.status.success(), "{variables:?}: {stderr}");
    (output.stdout, stderr)
}

/// `tblastn_on` with the `SUBJECTS` pieces.
fn tblastn(threads: usize, variables: &[(&str, &str)]) -> (Vec<u8>, String) {
    tblastn_on(&inputs().1, threads, variables)
}

/// The subjects whose (frame, chunk) units ran on the pool: the `[losat-thread-stage]` lines of
/// the unit stage with `parallel_selected=true`, one per subject searched by units.
fn units_on_the_pool(stderr: &str) -> usize {
    stderr
        .lines()
        .filter(|line| {
            line.starts_with("[losat-thread-stage] program=tblastn stage=prelim_units ")
                && line.contains(" parallel_selected=true ")
        })
        .count()
}

/// The `pool_entries` field of the one `[losat-thread-activity]` line of the subject stage.
fn pool_entries(stderr: &str) -> usize {
    let lines: Vec<_> = stderr
        .lines()
        .filter(|line| {
            line.starts_with("[losat-thread-activity] program=tblastn stage=subject_core ")
        })
        .collect();
    assert_eq!(lines.len(), 1, "one query batch, one line:\n{stderr}");
    lines[0]
        .split_whitespace()
        .find_map(|field| field.strip_prefix("pool_entries="))
        .unwrap_or_else(|| panic!("no pool_entries in {}", lines[0]))
        .parse()
        .expect("number of entries")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1413
// ```c
//     itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
//
//     /* iterate over all subject sequences */
//     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//            != BLAST_SEQSRC_EOF) {
// ```
// The one-thread search runs this loop over the subjects in order. The pool's entries of four
// chunks per thread (the default, and `LOSAT_X_DEFAULT=0` with `LOSAT_X_TBNOIDCHUNK=1`) and of
// one subject per thread (`LOSAT_X_TBNOIDCHUNK=0`, and `LOSAT_X_DEFAULT=0`) must give its
// report.
#[test]
fn the_report_is_the_one_thread_report_by_default_and_with_the_switch_off() {
    let (expected, _) = tblastn(1, &[]);
    let rows = expected
        .split(|&byte| byte == b'\n')
        .filter(|row| !row.is_empty());
    let hit_subjects: std::collections::BTreeSet<_> = rows
        .map(|row| {
            row.split(|&byte| byte == b'\t')
                .nth(1)
                .expect("subject column")
                .to_vec()
        })
        .collect();
    assert!(
        hit_subjects.len() > 1,
        "the search must hit several subjects"
    );
    for variables in [
        &[][..],
        &[("LOSAT_X_TBNOIDCHUNK", "0")],
        &[("LOSAT_X_DEFAULT", "0")],
        &[("LOSAT_X_DEFAULT", "0"), ("LOSAT_X_TBNOIDCHUNK", "1")],
    ] {
        let (report, _) = tblastn(THREADS, variables);
        assert!(
            report == expected,
            "{variables:?} at {THREADS} threads differs from one thread"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:697
// ```c++
//         Int8 tot_length = m_Atlas.GetSliceSize() / (4*m_NumThreads) + 1;
// ```
// 240 subjects of 1,000 bases: F = 60,000, and at `THREADS` = 4 a chunk may hold 3,751 bases,
// three subjects. With no variable an entry takes 4 * `THREADS` chunks (5 entries). With
// `LOSAT_X_TBNOIDCHUNK=0`, and with `LOSAT_X_DEFAULT=0` and no variable for the switch, it takes
// `THREADS` subjects (60 entries); `LOSAT_X_DEFAULT=0` with `LOSAT_X_TBNOIDCHUNK=1` gives the
// chunks again. The report is the same in all four.
#[test]
fn the_switch_is_on_by_default_and_follows_its_variable_and_the_global_one() {
    let debug = ("LOSAT_WASI_THREADS_DEBUG", "1");
    let chunks = SUBJECTS / (4 * THREADS * 3);
    let subjects = SUBJECTS / THREADS;
    let cases: [(&[(&str, &str)], usize); 4] = [
        (&[debug], chunks),
        (&[debug, ("LOSAT_X_TBNOIDCHUNK", "0")], subjects),
        (&[debug, ("LOSAT_X_DEFAULT", "0")], subjects),
        (
            &[
                debug,
                ("LOSAT_X_DEFAULT", "0"),
                ("LOSAT_X_TBNOIDCHUNK", "1"),
            ],
            chunks,
        ),
    ];
    let mut reports = Vec::new();
    for (variables, entries) in cases {
        let (report, stderr) = tblastn(THREADS, variables);
        assert_eq!(pool_entries(&stderr), entries, "{variables:?}");
        reports.push(report);
    }
    assert!(
        reports.iter().all(|report| *report == reports[0]),
        "the reports differ"
    );
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:697-707
// ```c++
//         Int8 tot_length = m_Atlas.GetSliceSize() / (4*m_NumThreads) + 1;
// ...
//         } while (res.length >= 0 && tot_length >= res.length && vol_oid < m_RestrictEnd);
// ```
// A subject longer than a chunk's budget is a chunk of its own, wherever it stands; the chunks
// before and after it are cut short or start anew. With the long subject first, in the middle
// and last the report at 2, 3 and 8 threads must be the one-thread report. The 248,000 bases
// give F = 62,000; at 2 threads a chunk may hold 7,751 bases (seven short subjects) and an entry
// takes eight chunks, so by hand: first, {long}, 28 chunks of 7 and {4} (30 chunks, 4 entries);
// in the middle, 14 of 7, {2}, {long}, 14 of 7 and {2} (31 chunks, 4 entries); last, 28 of 7,
// {4} and {long} (30 chunks, 4 entries).
#[test]
fn a_long_subject_first_in_the_middle_or_last_gives_the_one_thread_report() {
    let debug = ("LOSAT_WASI_THREADS_DEBUG", "1");
    for (position, entries_at_two) in [(Long::First, 4), (Long::Middle, 4), (Long::Last, 4)] {
        let subjects = long_subject_set(position);
        let (expected, _) = tblastn_on(&subjects, 1, &[]);
        assert!(
            expected
                .split(|&byte| byte == b'\n')
                .any(|row| row.split(|&byte| byte == b'\t').nth(1) == Some(&b"long"[..])),
            "{position:?}: the search must hit the long subject"
        );
        for threads in [2, 3, 8] {
            let (report, stderr) = tblastn_on(&subjects, threads, &[debug]);
            assert!(
                report == expected,
                "{position:?}: {threads} threads differ from one thread"
            );
            if threads == 2 {
                assert_eq!(pool_entries(&stderr), entries_at_two, "{position:?}");
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:333-344
// ```c++
//     begin_chunk = * state_obj;
//
//     // fill the cache for all sequence in mmaped slice
//     if (m_NumThreads) {
//         SSeqResBuffer * buffer = m_CachedSeqs[cacheID];
//         x_FillSeqBuffer(buffer, begin_chunk);
// ```
// With `LOSAT_X_TBNOIDNEXT` (off unless set) the threads take the next chunk when they are free,
// as NCBI's do, and the lists come back in OID order (`src/algorithm/tblastn/x_oid_flow.rs`): the
// search is one entry into the pool and the line names the chunks. The report must be the
// one-thread report on the 240 pieces and on the sets with one long subject, at 2, 3 and 8
// threads; with `LOSAT_X_TBNOIDCHUNK=0` there are no chunks and the switch does nothing.
#[test]
fn the_chunks_taken_when_a_thread_is_free_give_the_one_thread_report() {
    let debug = ("LOSAT_WASI_THREADS_DEBUG", "1");
    let next = ("LOSAT_X_TBNOIDNEXT", "1");
    let mut sets = vec![(None, inputs().1.clone())];
    for position in [Long::First, Long::Middle, Long::Last] {
        sets.push((Some(position), long_subject_set(position)));
    }
    for (position, subjects) in sets {
        let (expected, _) = tblastn_on(&subjects, 1, &[]);
        for threads in [2, 3, 8] {
            let (report, stderr) = tblastn_on(&subjects, threads, &[debug, next]);
            assert!(
                report == expected,
                "{position:?}: {threads} threads with the flow differ from one thread"
            );
            assert_eq!(
                pool_entries(&stderr),
                1,
                "{position:?} at {threads} threads"
            );
            assert!(
                stderr.contains(" flow_chunks="),
                "{position:?} at {threads} threads: no flow fields"
            );
        }
    }
    let (report, stderr) = tblastn(THREADS, &[debug, next, ("LOSAT_X_TBNOIDCHUNK", "0")]);
    assert!(
        report == tblastn(1, &[]).0,
        "LOSAT_X_TBNOIDCHUNK=0 with the flow"
    );
    assert_eq!(pool_entries(&stderr), SUBJECTS / THREADS);
    assert!(!stderr.contains(" flow_chunks="));
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-482
// ```c
//     while (TRUE) {
//         status = s_GetNextSubjectChunk(subject, &backup, kNucleotide,
//                                        dbseq_chunk_overlap);
//
//         if (status == SUBJECT_SPLIT_DONE) break;
// ```
// NCBI's thread searches the chunks of its subject one after another. With
// `LOSAT_X_TBNOIDUNITS` (off unless set) the (frame, chunk) units of a subject longer than a
// chunk's budget go to the pool, with the entries of `LOSAT_X_TBNOIDCHUNK` and with the flow of
// `LOSAT_X_TBNOIDNEXT`. The report at 2, 3 and 8 threads must be the one-thread report with the
// long subject first, in the middle and last. Only the long subject (48,000 bases; a chunk may
// hold 7,751, 5,167 and 1,938 bases at 2, 3 and 8 threads) runs its units on the pool, and no
// subject does without the switch.
#[test]
fn the_units_of_a_long_subject_on_the_pool_give_the_one_thread_report() {
    let debug = ("LOSAT_WASI_THREADS_DEBUG", "1");
    let units = ("LOSAT_X_TBNOIDUNITS", "1");
    let next = ("LOSAT_X_TBNOIDNEXT", "1");
    for position in [Long::First, Long::Middle, Long::Last] {
        let subjects = long_subject_set(position);
        let (expected, _) = tblastn_on(&subjects, 1, &[]);
        for threads in [2, 3, 8] {
            for (variables, on_the_pool) in [
                (&[debug, units][..], 1),
                (&[debug, units, next][..], 1),
                (&[debug][..], 0),
            ] {
                let (report, stderr) = tblastn_on(&subjects, threads, variables);
                assert!(
                    report == expected,
                    "{position:?}: {threads} threads with {variables:?} differ from one thread"
                );
                assert_eq!(
                    units_on_the_pool(&stderr),
                    on_the_pool,
                    "{position:?} at {threads} threads with {variables:?}"
                );
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1002-1003,1412-1413
// ```c
//     if ((status = BlastExtendWordNew(query->length, word_params,
//                                     &aux_struct->ewp)) != 0)
// ...
//     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//            != BLAST_SEQSRC_EOF) {
// ```
// NCBI's thread keeps its diagonal table from one subject to the next and only moves its offset
// on. With `LOSAT_X_TBNSTATE` (in the strict set: on by default) LOSAT's threads keep the scan
// state of the unit path the same way. The report at 1, 2, 3 and 8 threads must be the one-thread
// report without it (`LOSAT_X_TBNSTATE=0`), on the 240 pieces and on the sets with one long
// subject, also with the long subject's units on the pool (`LOSAT_X_TBNOIDUNITS`). With
// `LOSAT_X_TBNPARSHADOW` the run compares the unit path with the serial reference on every subject
// and stops if they differ.
#[test]
fn the_scan_state_kept_by_a_thread_gives_the_one_thread_report() {
    let without = ("LOSAT_X_TBNSTATE", "0");
    let units = ("LOSAT_X_TBNOIDUNITS", "1");
    let mut sets = vec![(None, inputs().1.clone())];
    for position in [Long::First, Long::Middle, Long::Last] {
        sets.push((Some(position), long_subject_set(position)));
    }
    for (position, subjects) in sets {
        let (expected, _) = tblastn_on(&subjects, 1, &[without]);
        for threads in [1, 2, 3, 8] {
            for variables in [&[][..], &[units][..]] {
                let (report, _) = tblastn_on(&subjects, threads, variables);
                assert!(
                    report == expected,
                    "{position:?}: {threads} threads with {variables:?} differ from one thread"
                );
            }
        }
    }
    let (expected, _) = tblastn(1, &[without]);
    for threads in [1, THREADS] {
        let (report, _) = tblastn(threads, &[("LOSAT_X_TBNPARSHADOW", "1")]);
        assert!(
            report == expected,
            "the shadow run at {threads} threads differs from one thread"
        );
    }
}
