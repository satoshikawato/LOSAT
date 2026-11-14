//! TBLASTX's merge of the HSP lists of two overlapping 5,000,000-letter parts of a subject frame
//! (`merge_tblastx_subject_chunk_hits`, `src/algorithm/tblastx/blast_engine/run_impl.rs`) with the
//! one-hit word finder (`-window_size 0`), whose HSPs may have a length of 0 or less, on whole
//! searches against NCBI BLAST+ 2.17.0.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:578-584
//! ```c
//!         Blast_HSPListAdjustOffsets(hsp_list, backup.offset);
//!         overlap = (backup.offset == backup.hard_ranges[backup.hm_index].left) ?
//!                   0 : dbseq_chunk_overlap;
//!         status = Blast_HSPListsMerge(&hsp_list, &combined_hsp_list,
//!                      kHspNumMax, &(backup.offset), INT4_MIN,
//!                      overlap, score_options->gapped_calculation,
//!                      Blast_ProgramIsMapping(program_number));
//! ```
//! The subject is written at test time: 15,000,300 letters of A (its frames read K and F only,
//! which give no word hit with the options below) with a segment copied at nt 420 and at nt
//! 14,999,760, so each copy lies in the overlap of the two parts of the minus frames and of the
//! plus frames (`MAX_DBSEQ_LEN` 5,000,000, `DBSEQ_CHUNK_OVERLAP` 100 letters). Both parts find the
//! HSPs there, and the merge joins the copies.
//! - `c3` (`tests/fasta/outfmt0/x2b_c3_*`, PAM30 `-seg no -threshold 13 -evalue 1e30`): HSPs of
//!   length -1 in both overlaps. LOSAT prints NCBI's bytes at 1 and 4 threads
//!   (`tests/fixtures/tblastx_chunk_merge/c3_big.ncbi_2.17.0.fmt6.out`: `tblastx -query
//!   tests/fasta/outfmt0/x2b_c3_q.fna -subject <the subject> ... -outfmt 6 -num_threads 1`, run
//!   from LOSAT/ with the subject's name `c3_big`).
//! - `r7_0148` (`tests/fasta/outfmt0/x2b_r7_0148_*`, PAM30 `-threshold 8 -evalue 1e30`): an HSP
//!   of length 0 that both parts find, whose merge NCBI computes by dividing by a sum of lengths of
//!   0 (its result is undefined): LOSAT refuses the search, exit 1 with its message
//!   (PD-LOSAT-NCBI-DEFECTS 1.19, rule 3).

use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

/// The subject's length and the two places of the segment (nt, 0-based).
const SUBJECT_LENGTH: usize = 15_000_300;
const SEGMENT_STARTS: [usize; 2] = [14_999_760, 420];

fn manifest_path(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(relative)
}

/// A path in the temporary directory that no other test of this process uses.
fn temp_path(tag: &str) -> PathBuf {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    std::env::temp_dir().join(format!(
        "losat_tblastx_chunk_merge_{}_{}_{tag}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ))
}

/// The subject of 15,000,300 letters named `name`, with the sequence of the FASTA file `segment`
/// at `SEGMENT_STARTS`, in a directory removed when dropped.
struct Subject {
    dir: PathBuf,
    path: PathBuf,
}

impl Subject {
    fn new(segment: &str, name: &str) -> Self {
        let text = std::fs::read_to_string(manifest_path(segment)).expect("read the segment");
        let segment: Vec<u8> = text
            .lines()
            .filter(|line| !line.starts_with('>'))
            .flat_map(|line| line.trim().bytes())
            .collect();
        let mut sequence = vec![b'A'; SUBJECT_LENGTH];
        for start in SEGMENT_STARTS {
            sequence[start..start + segment.len()].copy_from_slice(&segment);
        }
        let mut fasta = Vec::with_capacity(SUBJECT_LENGTH + SUBJECT_LENGTH / 80 + 64);
        fasta.extend_from_slice(format!(">{name}\n").as_bytes());
        for line in sequence.chunks(80) {
            fasta.extend_from_slice(line);
            fasta.push(b'\n');
        }
        let dir = temp_path(name);
        std::fs::create_dir_all(&dir).expect("create the directory");
        let path = dir.join(format!("{name}_s.fna"));
        std::fs::write(&path, fasta).expect("write the subject");
        Subject { dir, path }
    }
}

impl Drop for Subject {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// `losat tblastx` from LOSAT/ with no `LOSAT_*` variable but `LOSAT_BLASTX_WORKER_LOG` (a binary
/// built with the `blastx-worker-probe` feature needs that file name once a search has more than
/// one thread; it selects nothing).
fn tblastx(query: &str, subject: &Subject, args: &[&str], threads: u32) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let text = name.to_string_lossy();
        if text.starts_with("LOSAT_")
            || matches!(
                text.as_ref(),
                "BATCH_SIZE" | "CHUNK_SIZE" | "OVERLAP_CHUNK_SIZE" | "RAYON_NUM_THREADS"
            )
        {
            command.env_remove(name);
        }
    }
    command
        .current_dir(env!("CARGO_MANIFEST_DIR"))
        .env("LOSAT_BLASTX_WORKER_LOG", temp_path("worker.log"))
        .arg("tblastx")
        .args(["-query", query, "-subject"])
        .arg(&subject.path)
        .args(args)
        .args(["-outfmt", "6", "-num_threads", &threads.to_string()])
        .output()
        .expect("run LOSAT")
}

/// `c3`: the merges of HSPs of length -1 with their copies found by the other part give NCBI's
/// report at 1 and 4 threads.
#[test]
fn one_hit_hsps_of_negative_length_merge_across_parts_as_ncbi() {
    let subject = Subject::new("tests/fasta/outfmt0/x2b_c3_s.fna", "c3_big");
    let expected = std::fs::read(manifest_path(
        "tests/fixtures/tblastx_chunk_merge/c3_big.ncbi_2.17.0.fmt6.out",
    ))
    .expect("read NCBI's report");
    for threads in [1, 4] {
        let output = tblastx(
            "tests/fasta/outfmt0/x2b_c3_q.fna",
            &subject,
            &[
                "-matrix",
                "PAM30",
                "-seg",
                "no",
                "-threshold",
                "13",
                "-window_size",
                "0",
                "-evalue",
                "1e30",
            ],
            threads,
        );
        assert!(
            output.status.success(),
            "threads {threads}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(output.stderr.is_empty(), "threads {threads}: stderr");
        assert_eq!(
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&expected),
            "threads {threads}"
        );
    }
}

/// `r7_0148`: the HSP of length 0 that both parts find would be merged with a sum of lengths of
/// 0; the search is refused at 1 and 4 threads.
#[test]
fn a_merge_with_lengths_adding_up_to_zero_is_refused() {
    let subject = Subject::new("tests/fasta/outfmt0/x2b_r7_0148_s.fna", "r7_0148_big");
    for threads in [1, 4] {
        let output = tblastx(
            "tests/fasta/outfmt0/x2b_r7_0148_q.fna",
            &subject,
            &[
                "-matrix",
                "PAM30",
                "-threshold",
                "8",
                "-window_size",
                "0",
                "-evalue",
                "1e30",
            ],
            threads,
        );
        assert_eq!(output.status.code(), Some(1), "threads {threads}");
        assert!(output.stdout.is_empty(), "threads {threads}: no report");
        assert_eq!(
            String::from_utf8_lossy(&output.stderr),
            "Error: TBLASTX -window_size 0: a one-hit ungapped alignment whose length is 0 or less \
             meets another in the overlap of two 5,000,000-letter parts of a subject frame, and \
             their lengths add up to 0 (NCBI BLAST+ divides by that sum when it merges them; its \
             result is undefined), which LOSAT does not support\n",
            "threads {threads}"
        );
    }
}
