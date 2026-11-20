//! TBLASTX `-subject` records held packed (`LOSAT_X_TBXSTORE`, `src/algorithm/tblastx/blast_engine/
//! run_impl.rs` `read_subject_set`, `src/blastinput/x_packed_subjects.rs`) on whole searches.
//!
//! The switch is in the strict set (`src/utils/xswitch.rs`): on unless it is `0`. On the command
//! line each record's letters are packed as it is read and the search reads the stores as it reads
//! a database's sequences (`Subjects::Stored`), the reports by record. A set with a byte without an
//! ncbi4na code would keep every record's letters, but the reader gives none (it reads `U` as `T`);
//! `run_local`'s records are never packed.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:515-520
//! ```c
//!     virtual SBlastSequence GetCompressedPlusStrand() {
//!         SBlastSequence retval(size());
//!         string ncbi4na = kEmptyStr;
//!         m_SeqVector.GetSeqData(m_SeqVector.begin(), m_SeqVector.end(), ncbi4na);
//!         s_Ncbi4naToNcbi2na(ncbi4na, size(), retval.data.get());
//! ```
//! A store gives the frames of its record (the unit test of `subjects.rs` compares them in one
//! process), so the property under test is that the reports and messages do not depend on the
//! switch; that they equal NCBI's is the business of the fixture comparisons. `LOSAT_X_REACH`
//! shows the side each run took. The input (the first 6 kb of MjeNMV against the MelaMJNV genome
//! with IUPAC codes, lower-case letters and N blocks, one across a translation window, a record
//! without letters and the PemoMJNVA genome) gives HSPs in both genomes.

use std::path::{Path, PathBuf};
use std::process::Command;

mod reach_support;

const SWITCH: &str = "LOSAT_X_TBXSTORE";
const OFF: (&str, &str) = ("LOSAT_X_TBXSTORE", "0");

fn fasta(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fasta")
        .join(name)
}

/// The binary with the given variables as its only `LOSAT_*` variables that choose a code path.
/// `LOSAT_BLASTX_WORKER_LOG` stays: a binary built with the `blastx-worker-probe` feature
/// (`cargo test --all-features`, as in CI) needs that file name as soon as a search has more than
/// one thread, and it selects nothing.
fn losat(switches: &[(&str, &str)]) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.envs(switches.iter().copied());
    command
}

fn letters(name: &str) -> Vec<u8> {
    std::fs::read_to_string(fasta(name))
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(|line| line.trim().bytes())
        .collect()
}

fn record(defline: &str, letters: &[u8]) -> String {
    let mut text = format!(">{defline}\n");
    for line in letters.chunks(70) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    text
}

/// The query and subject files, in a directory of the test's own: `subjects` (MelaMJNV with
/// ambiguity codes, a record without letters, PemoMJNVA), `genomes` (the two genomes only, for
/// `-subject_loc`, whose interval would start past the end of the empty record and stop the
/// reading there) and `rna` (PemoMJNVA's `T` written `U` and `u` in part).
struct Inputs {
    dir: PathBuf,
    query: PathBuf,
    subjects: PathBuf,
    genomes: PathBuf,
    rna: PathBuf,
}

impl Inputs {
    fn new(test: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("losat-tblastx-store-{test}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let query = dir.join("q.fna");
        std::fs::write(&query, record("mje6k", &letters("MjeNMV.fasta")[..6000])).unwrap();
        let mut mela = letters("MelaMJNV.fasta");
        let iupac = b"NRYKMSWBDHVnrykmswbdhv-";
        for (n, at) in (1_000..mela.len()).step_by(997).enumerate() {
            mela[at] = iupac[n % iupac.len()];
        }
        // Lower-case letters where most HSPs are.
        mela[2_000..8_000].make_ascii_lowercase();
        // N blocks, one across the border of the first translation window (196,608 bases).
        for (start, len) in [(20_000, 300), (196_500, 250)] {
            mela[start..start + len].fill(b'N');
        }
        let mela = record("mela MelaMJNV with ambiguity codes", &mela);
        let pemo_letters = letters("PemoMJNVA.fasta");
        let pemo = record("pemo PemoMJNVA", &pemo_letters);
        let subjects = dir.join("s.fna");
        std::fs::write(&subjects, mela.clone() + ">empty no letters\n" + &pemo).unwrap();
        let genomes = dir.join("g.fna");
        std::fs::write(&genomes, mela.clone() + &pemo).unwrap();
        let mut rna_letters = pemo_letters;
        for (at, base) in rna_letters[..50_000].iter_mut().enumerate() {
            if *base == b'T' {
                *base = if at < 40_000 { b'U' } else { b'u' };
            }
        }
        let rna = dir.join("r.fna");
        std::fs::write(&rna, mela + &record("pemo PemoMJNVA as RNA", &rna_letters)).unwrap();
        Self {
            dir,
            query,
            subjects,
            genomes,
            rna,
        }
    }

    /// One tblastx search of `subjects` with `options`; returns its standard output and error and
    /// the sides of the switch the run took.
    fn tblastx(
        &self,
        subjects: &Path,
        options: &[&str],
        switches: &[(&str, &str)],
    ) -> (Vec<u8>, Vec<u8>, (bool, bool)) {
        let reach = reach_support::Reach::new("tblastx_store");
        let mut variables = switches.to_vec();
        variables.push(reach.var());
        let output = losat(&variables)
            .arg("tblastx")
            .arg("-query")
            .arg(&self.query)
            .arg("-subject")
            .arg(subjects)
            .args(options)
            .output()
            .expect("run LOSAT");
        assert!(
            output.status.success(),
            "{switches:?} {options:?}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        (output.stdout, output.stderr, reach.sides(SWITCH))
    }
}

impl Drop for Inputs {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.dir);
    }
}

/// The packed records give `=0`'s reports and messages: outfmt 0 (the subject set's numbers, the
/// empty record's warning, the alignments), outfmt 6 at 2 threads (the records packed on a second
/// thread) and outfmt 5 with `-subject_loc` (the records cut before they are packed, the subject
/// lengths); the default takes the new side only, `=0` records no side.
#[test]
fn packed_records_give_the_same_reports() {
    let inputs = Inputs::new("same");
    for (subjects, options) in [
        (&inputs.subjects, &["-outfmt", "0", "-num_threads", "1"][..]),
        (&inputs.subjects, &["-outfmt", "6", "-num_threads", "2"][..]),
        (
            &inputs.genomes,
            &["-outfmt", "5", "-subject_loc", "1500-150000"][..],
        ),
    ] {
        let (off, off_messages, off_side) = inputs.tblastx(subjects, options, &[OFF]);
        let text = String::from_utf8_lossy(&off);
        for subject in ["mela", "pemo"] {
            assert!(text.contains(subject), "{options:?}: no hit in {subject}");
        }
        assert_eq!(off_side, (false, false), "{options:?}: =0's side");
        let (report, messages, reached) = inputs.tblastx(subjects, options, &[]);
        assert_eq!(reached, (true, false), "{options:?}: default's side");
        assert_eq!(report, off, "{options:?}: report against =0");
        assert_eq!(messages, off_messages, "{options:?}: messages against =0");
    }
    let (_, messages, _) = inputs.tblastx(&inputs.subjects, &["-outfmt", "6"], &[]);
    assert!(
        String::from_utf8_lossy(&messages).contains("Subject sequence contains no data"),
        "the empty record's warning"
    );
}

/// A record with `U` and `u` is packed too (the reader gives `T` and `t`, `fasta.cpp`'s
/// `CSeqportUtil::Pack`) and gives `=0`'s report.
#[test]
fn a_record_with_u_is_packed_as_t() {
    let inputs = Inputs::new("rna");
    let options = ["-outfmt", "6"];
    let (off, _, _) = inputs.tblastx(&inputs.rna, &options, &[OFF]);
    let (report, _, reached) = inputs.tblastx(&inputs.rna, &options, &[]);
    assert_eq!(reached, (true, false), "default's side");
    assert_eq!(report, off, "against =0");
}
