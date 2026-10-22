#![allow(warnings, clippy::all)]

//! outfmt 8 (ASN.1 text), 9 (ASN.1 binary) and 12 (JSON) of the LOSAT CLI against NCBI BLAST+ 2.17.0's output
//! for the same `-subject` searches, byte for byte, at one and four threads: BLASTN's three
//! tasks (megablast, blastn, dc-megablast), BLASTP, BLASTX (gapped and `-ungapped`), TBLASTN
//! (with and without `-subject_loc`) and TBLASTX on the outfmt sample inputs
//! (`tests/fixtures/seqalign/in/`; outfmt-20261019 `ncbi-samples/fixture/`).
//!
//! The translated searches' files: `ncbi-samples/out/<name>.out` of `blastx -query
//! fixture/qnucl.fa -subject fixture/sprot.fa`, `tblastn -query fixture/qprot.fa -subject
//! fixture/snucl.fa` and `tblastx -query fixture/qnucl.fa -subject fixture/snucl.fa`
//! (`-num_threads 1 -outfmt 8|9|12`), and the OF3 design probes (`S/tmp/probes/<name>.out`) with
//! `-ungapped -comp_based_stats 0` (BLASTX) and `-subject_loc 31-400` (TBLASTN). TBLASTX's
//! JSON (149 KB) is left to outfmt-20261019 `S/scripts/compare_samples_of3b2.py`.
//!
//! The expected files are NCBI's output, unchanged (`ncbi-samples/out/<name>.out` of
//! `blastn -task <task> -query fixture/qnucl.fa -subject fixture/snucl.fa -num_threads 1
//! -outfmt 8|9|12` and `blastp -query fixture/qprot.fa -subject fixture/sprot.fa ...`); NCBI
//! writes nothing to standard error and exits with 0. A `-subject` search's documents name no
//! path, so the searches run from the fixture folder with its own file names. (The `-db`
//! searches are compared by outfmt-20261019 `S/scripts/compare_samples_of3b1.py`; their frozen
//! fixtures belong to `db_search_manifest.tsv`.)
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-698
//! ```c++
//!     if (m_FormatType == CFormattingArgs::eAsnText || m_FormatType == CFormattingArgs::eJsonSeqalign) {
//!         if (results.HasAlignments()) {
//!     	    CRef<CSeq_align_set> aln_set (new CSeq_align_set);
//!             CBlastFormatUtil::PruneSeqalign(*(results.GetSeqAlign()), *aln_set, m_HitlistSize);
//!             if(m_FormatType == CFormattingArgs::eAsnText)
//!             	m_Outfile << MSerial_AsnText << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
//!             else
//!             	m_Outfile << MSerial_Json << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
//!         }
//! ```

use std::path::PathBuf;
use std::process::Command;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/seqalign")
}

/// The CLI without the caller's LOSAT_* variables (as `cli_smoke.rs`; the BLASTX worker log
/// of an `--all-features` build stays).
fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        let name_text = name.to_string_lossy();
        if name_text.starts_with("LOSAT_")
            && name != "LOSAT_X_DEFAULT"
            && name != "LOSAT_BLASTX_WORKER_LOG"
        {
            command.env_remove(name);
        }
    }
    command.current_dir(fixtures());
    command
}

#[test]
fn subject_searches_write_ncbis_seq_annots() {
    let searches: [(&str, &[&str]); 4] = [
        (
            "blastn-megablast_subject",
            &[
                "blastn",
                "-task",
                "megablast",
                "-query",
                "in/qnucl.fa",
                "-subject",
                "in/snucl.fa",
            ],
        ),
        (
            "blastn-blastn_subject",
            &[
                "blastn",
                "-task",
                "blastn",
                "-query",
                "in/qnucl.fa",
                "-subject",
                "in/snucl.fa",
            ],
        ),
        (
            "blastn-dcmegablast_subject",
            &[
                "blastn",
                "-task",
                "dc-megablast",
                "-query",
                "in/qnucl.fa",
                "-subject",
                "in/snucl.fa",
            ],
        ),
        (
            "blastp_subject",
            &["blastp", "-query", "in/qprot.fa", "-subject", "in/sprot.fa"],
        ),
    ];
    for (name, argv) in searches {
        for outfmt in ["8", "9", "12"] {
            let expected = std::fs::read(fixtures().join(format!("{name}_fmt{outfmt}.out")))
                .expect("NCBI's output");
            for threads in ["1", "4"] {
                let output = losat()
                    .args(argv)
                    .args(["-outfmt", outfmt, "-num_threads", threads])
                    .output()
                    .expect("run LOSAT");
                let case = format!("{name} -outfmt {outfmt} -num_threads {threads}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stderr),
                    "",
                    "{case}: standard error"
                );
                assert!(output.status.success(), "{case}: {:?}", output.status);
                assert!(
                    output.stdout == expected,
                    "{case}: the documents differ from NCBI's"
                );
            }
        }
    }
}

/// NCBI writes no document for a query without hits, and nothing at all when no query has
/// hits (blast_format.cpp:690-698 above: `results.HasAlignments()`); NCBI BLAST+ 2.17.0 writes
/// 0 bytes to both streams and exits with 0 for this search (outfmt-20261019 S, 2026-10-19).
#[test]
fn a_search_without_hits_writes_nothing() {
    for outfmt in ["8", "9", "12"] {
        let output = losat()
            .args([
                "blastp",
                "-query",
                "in/qprot.fa",
                "-subject",
                "in/sprot.fa",
                "-evalue",
                "1e-300",
            ])
            .args(["-outfmt", outfmt])
            .output()
            .expect("run LOSAT");
        assert!(output.status.success(), "{:?}", output.status);
        assert!(
            output.stdout.is_empty() && output.stderr.is_empty(),
            "outfmt {outfmt}"
        );
    }
}

/// The Std-seg forms of the translated programs: one Std-seg per edit-script operation with
/// the translated row in nucleotides (BLASTX's query, TBLASTN's subject), and the `diags`
/// Seq-aligns of an ungapped search (TBLASTX, BLASTX `-ungapped`) with each HSP's scores in
/// its Std-seg (blast_seqalign.cpp:306-368, 1323-1447).
#[test]
fn translated_subject_searches_write_ncbis_seq_annots() {
    let searches: [(&str, &[&str], &[&str]); 5] = [
        (
            "blastx_subject",
            &["blastx", "-query", "in/qnucl.fa", "-subject", "in/sprot.fa"],
            &["8", "9", "12"],
        ),
        (
            "blastx_subject_ungapped",
            &[
                "blastx",
                "-query",
                "in/qnucl.fa",
                "-subject",
                "in/sprot.fa",
                "-ungapped",
                "-comp_based_stats",
                "0",
            ],
            &["8"],
        ),
        (
            "tblastn_subject",
            &[
                "tblastn",
                "-query",
                "in/qprot.fa",
                "-subject",
                "in/snucl.fa",
            ],
            &["8", "9", "12"],
        ),
        (
            "tblastn_subject_sloc",
            &[
                "tblastn",
                "-query",
                "in/qprot.fa",
                "-subject",
                "in/snucl.fa",
                "-subject_loc",
                "31-400",
            ],
            &["8"],
        ),
        (
            "tblastx_subject",
            &[
                "tblastx",
                "-query",
                "in/qnucl.fa",
                "-subject",
                "in/snucl.fa",
            ],
            &["8", "9"],
        ),
    ];
    for (name, argv, formats) in searches {
        for &outfmt in formats {
            let expected = std::fs::read(fixtures().join(format!("{name}_fmt{outfmt}.out")))
                .expect("NCBI's output");
            for threads in ["1", "4"] {
                let output = losat()
                    .args(argv)
                    .args(["-outfmt", outfmt, "-num_threads", threads])
                    .output()
                    .expect("run LOSAT");
                let case = format!("{name} -outfmt {outfmt} -num_threads {threads}");
                assert_eq!(
                    String::from_utf8_lossy(&output.stderr),
                    "",
                    "{case}: standard error"
                );
                assert!(output.status.success(), "{case}: {:?}", output.status);
                assert!(
                    output.stdout == expected,
                    "{case}: the documents differ from NCBI's"
                );
            }
        }
    }
}
