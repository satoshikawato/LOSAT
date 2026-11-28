//! A prefetch stop after query batches that the program still holds (review of line RG series D's
//! D-SO-5, L3; 2026-11-28): two sequences hold the stopping subject's accession under gis of
//! their own (`gi|8001|ref|NP_000002.1|` and `gi|8002|ref|NP_000002.1|`, made here), the first
//! batch reports the first, the second batch the second, and the third batch reports the
//! accession's own sequence (`ref|NP_000002.1|` of `multidb/ldp`), whose lookup fails in a list
//! with a database that has no accession table (`multidb/lda`). NCBI BLAST+ 2.17.0 does not look
//! that accession up: its data loader holds it from the second batch's load, the last one, and
//! NCBI writes the second holder's letters for the subject (exit status 0; the review's probe,
//! D/logs/tbx_stream_off/l3probe). LOSAT refuses the batch, naming that holder (OID 1), in BLASTP
//! and TBLASTN, which hold their batches, and in BLASTN both as it streams its reports and with
//! `LOSAT_X_RPTSTREAM=0`, which holds them (`common::database::TaxdbRun::loads_at_stop`: each held
//! batch a note of its own, as the stream records them).

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn losat() -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_losat"));
    for (name, _) in std::env::vars_os() {
        if name.to_string_lossy().starts_with("LOSAT_") && name != "LOSAT_BLASTX_WORKER_LOG" {
            command.env_remove(name);
        }
    }
    command.env_remove("BLASTDB");
    command
}

fn temp_dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "losat-db-stop-held-batches-{tag}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// `count` random residues after a methionine (no stop codon in their back-translation).
fn protein(seed: u64, count: usize) -> String {
    let mut state = seed;
    std::iter::once('M')
        .chain((1..count).map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            char::from(b"ACDEFGHIKLMNPQRSTVWY"[((state >> 33) % 20) as usize])
        }))
        .collect()
}

/// The letters of `protein` coded by one codon each.
fn coding(protein: &str) -> String {
    protein
        .chars()
        .map(|residue| {
            let index = "ACDEFGHIKLMNPQRSTVWY".find(residue).unwrap();
            [
                "GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA", "CTG", "ATG", "AAT",
                "CCG", "CAG", "CGT", "TCT", "ACC", "GTT", "TGG", "TAT",
            ][index]
        })
        .collect()
}

fn record(defline: &str, letters: &str) -> String {
    let mut text = format!(">{defline}\n");
    for line in letters.as_bytes().chunks(60) {
        text.push_str(std::str::from_utf8(line).unwrap());
        text.push('\n');
    }
    text
}

/// The record of `fasta` whose definition line starts with `id`.
fn record_of(fasta: &Path, id: &str) -> String {
    let text = format!("\n{}", std::fs::read_to_string(fasta).unwrap());
    let start = text.find(&format!("\n>{id}")).unwrap() + 1;
    let end = text[start..]
        .find("\n>")
        .map_or(text.len(), |end| start + end + 1);
    text[start..end].to_owned()
}

#[test]
fn a_stop_after_held_batches_names_the_holder_of_the_later_batch() {
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    let multidb = manifest.join("tests/db_search/multidb");
    let dir = temp_dir("names");
    let database = |name: &str, defline: &str, letters: &str, dbtype: &str| {
        let fasta = dir.join(format!("{name}.fa"));
        std::fs::write(&fasta, record(defline, letters)).unwrap();
        let db = dir.join(name).to_str().unwrap().to_owned();
        let output = losat()
            .arg("makeblastdb")
            .arg("-in")
            .arg(&fasta)
            .args(["-dbtype", dbtype, "-parse_seqids", "-blastdb_version", "5"])
            .args(["-out", &db])
            .output()
            .unwrap();
        assert!(output.status.success(), "{output:?}");
        db
    };
    let (one, two) = (protein(81, 200), protein(82, 200));
    // The names `h1*` and `h2*` sort before `lda*` and `ldp*`, so that the sequences' OIDs follow
    // the list's order wherever the run's directory is (NCBI BLAST+ and LOSAT gave `multidb`'s
    // volumes the first OIDs for names `one*` and `two*` in a directory under the home's `.cache`).
    let holders = |accession: &str, one: &str, two: &str, dbtype: &str| {
        [("8001", "one", one), ("8002", "two", two)].map(|(gi, which, letters)| {
            database(
                &format!("h{}{dbtype}", &gi[3..]),
                &format!("gi|{gi}|ref|{accession}| holder {which} of the stopping accession"),
                letters,
                dbtype,
            )
        })
    };
    let [pa, pb] = holders("NP_000002.1", &one, &two, "prot");
    let [na, nb] = holders("NM_000002.1", &coding(&one), &coding(&two), "nucl");
    let queries = |name: &str, one: &str, two: &str, third: String| {
        let path = dir.join(name);
        std::fs::write(
            &path,
            record("q1 of holder one", one) + &record("q2 of holder two", two) + &third,
        )
        .unwrap();
        path
    };
    let proteins = queries(
        "q.faa",
        &one,
        &two,
        record_of(&multidb.join("mix2_q.faa"), "P69905.2"),
    );
    let nucleotides = queries(
        "q.fna",
        &coding(&one),
        &coding(&two),
        record_of(&multidb.join("mix2_q.fna"), "ENST00000320868.9"),
    );
    let list = |a: &str, b: &str, suffix: &str| {
        format!(
            "{a} {b} {} {}",
            multidb.join(format!("lda{suffix}")).display(),
            multidb.join(format!("ldp{suffix}")).display()
        )
    };
    let run = |program: &str, query: &Path, dbs: &str, stream: bool| -> Output {
        let mut command = losat();
        command
            .arg(program)
            .arg("-query")
            .arg(query)
            .args(["-db", dbs, "-evalue", "1e-10", "-outfmt", "6"])
            .env("BATCH_SIZE", "1");
        if !stream {
            command.env("LOSAT_X_RPTSTREAM", "0");
        }
        command.output().unwrap()
    };
    for (program, query, dbs, stream, accession) in [
        ("blastp", &proteins, list(&pa, &pb, ""), true, "NP_000002.1"),
        (
            "tblastn",
            &proteins,
            list(&na, &nb, "n"),
            true,
            "NM_000002.1",
        ),
        (
            "blastn",
            &nucleotides,
            list(&na, &nb, "n"),
            true,
            "NM_000002.1",
        ),
        (
            "blastn",
            &nucleotides,
            list(&na, &nb, "n"),
            false,
            "NM_000002.1",
        ),
    ] {
        let output = run(program, query, &dbs, stream);
        let case = format!("{program} stream={stream}");
        assert_eq!(output.status.code(), Some(1), "{case}: {output:?}");
        let stderr = String::from_utf8_lossy(&output.stderr);
        let last = stderr.lines().last().unwrap_or_default();
        assert!(
            last.starts_with(&format!(
                "Error: database sequence 3 has the Seq-align id ref|{accession}|, whose lookup in the database fails"
            )) && last.contains(
                "but which NCBI's prefetch may not look up: another sequence of the database (OID 1) holds an id"
            ),
            "{case}: {stderr}"
        );
        // BLASTN's stream writes the two batches before; the programs that hold them refuse
        // before the first byte.
        let rows: Vec<String> = String::from_utf8_lossy(&output.stdout)
            .lines()
            .map(|line| line.split('\t').take(2).collect::<Vec<_>>().join(" "))
            .collect();
        let want: &[&str] = match (program, stream) {
            ("blastn", true) => &["q1 NM_000002.1", "q2 NM_000002.1"],
            _ => &[],
        };
        assert_eq!(rows, want, "{case}");
    }
    let _ = std::fs::remove_dir_all(&dir);
}
