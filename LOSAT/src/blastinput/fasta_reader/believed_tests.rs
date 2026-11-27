use super::super::{read_all, FastaInputSource, ReaderConfig};
use super::*;
use crate::seqid::id_scope::IdScope;

fn read(fasta: &[u8], protein: bool) -> Result<Vec<super::super::FastaRecord>, String> {
    let config = ReaderConfig::query("BLASTN", protein, false).with_parse_deflines(true);
    let mut source = FastaInputSource::from_bytes(fasta, config);
    read_all(&mut source, &mut |_| Ok(())).map_err(|error| format!("{error:#}"))
}

fn ids_of(record: &super::super::FastaRecord) -> &BelievedIds {
    record.believed().expect("-parse_deflines keeps the ids")
}

// NCBI BLAST+ 2.17.0 -parse_deflines (probe `DI7/ws/pd_probe.sh`, `pdo/*`): `Query=`, `qseqid`,
// `qacc`, `qaccver` and `qgi` of each defline.
#[test]
fn the_ids_of_a_defline_follow_the_ncbi_oracle() {
    for (defline, line, qseqid, qacc, qaccver, gi) in [
        (
            &b"gi|42|ref|NP_0042.1| t"[..],
            &b"NP_0042.1 t"[..],
            &b"gi|42|ref|NP_0042.1|"[..],
            &b"NP_0042"[..],
            &b"NP_0042.1"[..],
            42,
        ),
        (
            b"ref|NP_0043.1| t",
            b"NP_0043.1 t",
            b"ref|NP_0043.1|",
            b"NP_0043",
            b"NP_0043.1",
            0,
        ),
        (b"lcl|q3 t", b"q3 t", b"q3", b"q3", b"q3", 0),
        (
            b"qplain4 t",
            b"qplain4 t",
            b"qplain4",
            b"qplain4",
            b"qplain4",
            0,
        ),
        (
            b"gnl|qdb|q5 t",
            b"qdb:q5 t",
            b"gnl|qdb|q5",
            b"qdb:q5",
            b"qdb:q5",
            0,
        ),
        (b"gi|0 t", b"gi|0 t", b"gi|0", b"0", b"0", 0),
        (
            b"prf||123456A t",
            b"prf||123456A t",
            b"prf||123456A",
            b"123456A",
            b"123456A",
            0,
        ),
        (
            b"pdb|1ABC|A t",
            b"1ABC_A t",
            b"pdb|1ABC|A",
            b"1ABC_A",
            b"1ABC_A",
            0,
        ),
        (b"q1", b"q1", b"q1", b"q1", b"q1", 0),
        (b"a,b t", b"a_b t", b"a_b", b"a_b", b"a_b", 0),
    ] {
        let fasta = [b">", defline, b"\nACGTACGT\n"].concat();
        let records = read(&fasta, false).unwrap();
        let ids = ids_of(&records[0]);
        let text = String::from_utf8_lossy(defline).into_owned();
        assert_eq!(&*records[0].query_line(), line, "{text}");
        assert_eq!(ids.tabular[0], qseqid, "{text}");
        assert_eq!(ids.tabular[1], qacc, "{text}");
        assert_eq!(ids.tabular[2], qaccver, "{text}");
        assert_eq!(ids.gi, gi, "{text}");
    }
}

// NCBI BLAST+ 2.17.0 (probe `edge3`, `blastn.gen.*`): the generated ids count the records
// without ids, not the records.
#[test]
fn generated_ids_count_the_records_without_ids() {
    let records = read(b">has1 a\nACGT\n>\nACGT\n>has3 b\nACGT\n>\nACGT\n", false).unwrap();
    let lines: Vec<_> = records
        .iter()
        .map(|record| record.query_line().into_owned())
        .collect();
    assert_eq!(lines, [&b"has1 a"[..], b"Query_1", b"has3 b", b"Query_2"]);
}

// NCBI BLAST+ 2.17.0 (probe `edge2`, `range.*`): a range replaces the ids by a generated one.
#[test]
fn a_range_gives_a_generated_id_and_keeps_the_title() {
    let records = read(b">rq:10-40 t\nACGT\n", false).unwrap();
    assert_eq!(&*records[0].query_line(), b"Query_1 t");
    assert_eq!(records[0].title, b"t");
}

/// The scope of the query records `fasta` of `program` (`QueryScopes`), read as one batch.
fn one_scope(fasta: &[u8], program: &'static str) -> Result<(), String> {
    let config = ReaderConfig::query(program, false, false).with_parse_deflines(true);
    let mut source = FastaInputSource::from_bytes(fasta, config);
    let records = read_all(&mut source, &mut |_| Ok(())).map_err(|e| e.to_string())?;
    crate::blastinput::query_batch::QueryScopes::of_records(
        crate::blastinput::query_batch::QueryScopeSpan::Batch,
        program,
        &records,
    )
    .map(drop)
    .map_err(|error| error.to_string())
}

// Two queries of one scope (series Q of line RG, Q1c): the reader reads them, the batch that
// holds both refuses them, two batches do not (NCBI empties its scope after a batch) but with
// `-subject` and outfmt 5 (`QueryScopeSpan::Run`).
#[test]
fn two_records_with_one_id_are_refused_in_one_scope() {
    for fasta in [
        &b">dup1 a\nACGT\n>dup1 b\nACGT\n"[..],
        b">abc a\nACGT\n>ABC b\nACGT\n",
        b">gb|AB000001.1| a\nACGT\n>emb|AB000001.1| b\nACGT\n",
        b">gi|5|ref|NM_000001.1| a\nACGT\n>gi|5 b\nACGT\n",
        b">gb|AAA00007.1|NAMEA a\nACGT\n>gb||NAMEA b\nACGT\n",
        b">gb||NAMEA a\nACGT\n>gb|AAA00007.1|NAMEA b\nACGT\n",
    ] {
        let records = read(fasta, false).unwrap();
        let error = one_scope(fasta, "BLASTP").unwrap_err();
        assert!(error.contains("two query records"), "{error}");
        use crate::blastinput::query_batch::{QueryScopeSpan, QueryScopes};
        let mut batches = QueryScopes::new(QueryScopeSpan::Batch, "BLASTP");
        batches.add_batch(&records[..1], &[false], 0).unwrap();
        batches.add_batch(&records[1..], &[false], 1).unwrap();
        let mut run = QueryScopes::new(QueryScopeSpan::Run, "BLASTP");
        run.add_batch(&records[..1], &[false], 0).unwrap();
        let error = run
            .add_batch(&records[1..], &[false], 1)
            .unwrap_err()
            .to_string();
        assert!(error.contains("two query records"), "{error}");
        let config = ReaderConfig::subject("BLASTN", false, false).with_parse_deflines(true);
        let mut source = FastaInputSource::from_bytes(fasta, config);
        let error = read_all(&mut source, &mut |_| Ok(()))
            .unwrap_err()
            .to_string();
        assert!(error.contains("two subject records"), "{error}");
    }
}

// NCBI BLAST+ 2.17.0 (probes of the review): NCBI's scope finds these as two sequences.
#[test]
fn a_version_or_a_chain_tells_two_ids_apart() {
    for fasta in [
        &b">NM_000001.1 a\nACGT\n>NM_000001.2 b\nACGT\n"[..],
        b">pdb|1ABC|A a\nACGT\n>pdb|1ABC|B b\nACGT\n",
        b">gb|AAA00007.1|NAMEA a\nACGT\n>gb|AAA00008.1|NAMEA b\nACGT\n",
        b">sp|P12345.2|NAME_X a\nACGT\n>sp|Q99999.1|NAME_X b\nACGT\n",
    ] {
        one_scope(fasta, "BLASTN").unwrap();
    }
}

// NCBI BLAST+ 2.17.0 (`db.pd.*.ver`): the same letters and reports, but BLASTN warns.
#[test]
fn an_accession_with_and_without_a_version_is_refused_by_blastn_only() {
    let fasta = b">NM_000005 a\nACGT\n>NM_000005.1 b\nACGT\n";
    let error = one_scope(fasta, "BLASTN").unwrap_err();
    assert!(error.contains("with and without a version"), "{error}");
    one_scope(fasta, "TBLASTX").unwrap();
}

// fasta.cpp:1005, 2117-2125: the message names the molecule that the accession named.
#[test]
fn the_invalid_letters_message_names_the_molecule_of_the_accession() {
    let message = |fasta: &[u8], protein: bool| {
        let config = ReaderConfig::query("BLASTN", protein, false).with_parse_deflines(true);
        let mut source = FastaInputSource::from_bytes(fasta, config);
        let mut text = Vec::new();
        read_all(&mut source, &mut |line| {
            text.extend_from_slice(line);
            Ok(())
        })
        .unwrap();
        String::from_utf8(text).unwrap()
    };
    assert!(message(b">NM_000009.1 x\nACGTACGTAC!GTACGTACGT\n", false)
        .starts_with("FASTA-Reader: Ignoring invalid nucleotide residues"));
    assert!(message(b">NP_000009.1 x\nMKLVMKLVMK!LVMKLVMKLV\n", true)
        .starts_with("FASTA-Reader: Ignoring invalid protein residues"));
    assert!(message(b">plain x\nACGTACGTAC!GTACGTACGT\n", false)
        .starts_with("FASTA-Reader: Ignoring invalid residues"));
}

#[test]
fn a_geninfo_import_id_is_refused() {
    let error = read(b">gim|123 a\nACGT\n", false).unwrap_err();
    assert!(error.contains("Geninfo import id"), "{error}");
}

#[test]
fn a_database_ordinal_is_refused() {
    let error = read(b">gnl|BL_ORD_ID|5 a\nACGT\n", false).unwrap_err();
    assert!(error.contains("database ordinal"), "{error}");
}

// `xSetSeqMol` + `AssignMolType` + `CSeqportUtil::Pack`: NCBI reads each data line with the
// molecule the accession names (warning on letters invalid for that molecule), then encodes the
// record for the program's molecule (line RG-XT item X7-5).
#[test]
fn other_molecule_accession_follows_ncbis_reading() {
    let read_with_warnings = |fasta: &[u8], protein: bool| {
        let config = ReaderConfig::query("BLASTN", protein, false).with_parse_deflines(true);
        let mut source = FastaInputSource::from_bytes(fasta, config);
        let mut warnings = Vec::new();
        let records = read_all(&mut source, &mut |line| {
            warnings.extend_from_slice(line);
            Ok(())
        })
        .unwrap();
        (records, String::from_utf8(warnings).unwrap())
    };
    // Protein program reading a nucleotide-named accession: non-IUPACna letters are dropped with
    // NCBI's nucleotide warning; kept nucleotide letters (including U/u) stay as protein letters.
    let (records, warnings) = read_with_warnings(b">NM_000001.1 a\nMKLVEFu\n", true);
    assert_eq!(records[0].seq(), b"MKVu");
    assert_eq!(
        warnings,
        "FASTA-Reader: Ignoring invalid nucleotide residues at position(s): On line 2: 3, 5-6\n"
    );
    // Nucleotide program reading a protein-named accession: all protein letters are accepted
    // without warning; `CSeqportUtil::Pack` maps U->T and non-IUPACna letters to N, keeping masks
    // across `>?` gap lines and on ranged accessions (`xSetSeqMol` runs before `PostProcessIDs`).
    let nuc_cases: [(&[u8], &[u8]); 5] = [
        (b">NP_000001.1 a\nACGUefILPQXZ*\n", b"ACGTnnNNNNNNN"),
        (b">NP_000001.1 a\nACGTNNRY\n", b"ACGTNNRY"),
        (b">NP_000001.1 a\nACGU\n", b"ACGT"),
        (b">NP_000001.1 a\nae\n>?2\nfgA\n", b"annnngA"),
        (b">NP_000001.1:10-40 a\nACGEF\n", b"ACGNN"),
    ];
    for (fasta, expected) in nuc_cases {
        let (records, warnings) = read_with_warnings(fasta, false);
        assert_eq!(records[0].seq(), expected);
        assert_eq!(records[0].believed().unwrap().named_molecule, Some(true));
        assert!(warnings.is_empty(), "{warnings}");
    }
    // Protein program with a nucleotide-named accession: lowercase mask stays open across a `>?`
    // gap line (`X`), and ranged accessions still set `named_molecule` to nucleotide.
    let (records, warnings) = read_with_warnings(b">NM_000001.1 a\naeu\n>?2\nfgcA\n", true);
    assert_eq!(records[0].seq(), b"auxxgcA");
    assert!(
        warnings.contains("On line 2: 2") && warnings.contains("On line 4: 1"),
        "{warnings}"
    );
    let (records, warnings) = read_with_warnings(b">NM_000001.1:10-40 a\nMKLVEF\n", true);
    assert_eq!(records[0].seq(), b"MKV");
    assert_eq!(records[0].believed().unwrap().named_molecule, Some(false));
    assert!(
        warnings.contains("invalid nucleotide residues"),
        "{warnings}"
    );
    assert!(read(b">NM_000001.1 a\nACGT\n", false).is_ok());
    assert!(read(b">NP_000001.1 a\nMKLVEF\n", true).is_ok());
}

// fasta.cpp:386-397: without fDLOptional a defline is required.
#[test]
fn a_record_without_a_defline_is_nothing_read() {
    let error = read(b"ACGT\n", false).unwrap_err();
    assert!(
        error.contains("Input doesn't start with a defline"),
        "{error}"
    );
}

#[test]
fn the_archive_form_follows_the_id() {
    let local = read(b">lcl|q a\nACGT\n", false).unwrap();
    assert!(local[0].archive_seq_loc().is_none());
    let general = read(b">gnl|db|q a\nACGT\n", false).unwrap();
    assert!(general[0].archive_seq_loc().is_none());
    let mut text = read(b">NM_000001.1 a\nACGTACGT\n", false).unwrap();
    let loc = text[0].archive_seq_loc().unwrap();
    assert_eq!((loc.from, loc.to), (0, 7));
    text[0].cut_in_place(2, 6);
    let loc = text[0].archive_seq_loc().unwrap();
    assert_eq!((loc.from, loc.to), (2, 5));
}

/// The records of `fasta` read by the reader of `role` of `program` in the search `scope`.
fn read_in(
    scope: &std::sync::Arc<IdScope>,
    role: &str,
    program: &'static str,
    protein: bool,
    fasta: &[u8],
) -> Result<Vec<super::super::FastaRecord>, String> {
    let config = if role == "query" {
        ReaderConfig::query(program, protein, false)
    } else {
        ReaderConfig::subject(program, protein, false)
    };
    let mut source = FastaInputSource::from_bytes(fasta, config.with_parse_deflines(true));
    source.set_id_scope(Some(scope.clone()));
    read_all(&mut source, &mut |_| Ok(())).map_err(|error| format!("{error:#}"))
}

// A query and a subject share no id, nor one accession with and without a version, unless they
// are one record: defline and letters (DECISIONS 12, 15, 19).
#[test]
fn a_query_and_a_subject_share_an_id_only_as_one_record() {
    const QUERY: &[u8] =
        b">q1 t\nACGT\n>NM_000001 t\nACGT\n>ref|NM_000002.1| t\nACGT\n>NM_000003.1 t\nACGT\n\
              >gb|AAA00007.1|NAMEA t\nACGT\n>gb||NAMEB t\nACGT\n";
    for (subject, refused) in [
        (&b">s1 t\nACGT\n"[..], false),
        (b">NM_000002.2 t\nACGT\n", false),
        (b">q1 t\nACGT\n", false),
        (b">ref|NM_000002.1| t\nACGT\n", false),
        (b">lcl|q1 t\nACGT\n", true),
        (b">q1  t\nACGT\n", true),
        (b">q1 t\nAC\nGT\n", false),
        (b">q1 t\nacgu\n", false),
        (b">q1 other\nACGT\n", true),
        (b">q1\nACGT\n", true),
        (b">q1 t\nACGA\n", true),
        (b">NM_000001.1 t\nACGT\n", true),
        (b">NM_000002 t\nACGT\n", true),
        (b">NM_000002.1 t\nACGT\n", true),
        (b">NM_000003 t\nACGT\n", true),
        (b">gb|AAA00007.1|NAMEA t\nACGT\n", false),
        (b">gb||NAMEB t\nACGT\n", false),
        (b">gb|AAA00008.1|NAMEA t\nACGT\n", false),
        (b">gb||NAMEA t\nACGT\n", true),
        (b">gb|AAA00009.1|NAMEB t\nACGT\n", true),
    ] {
        let scope = std::sync::Arc::new(IdScope::default());
        read_in(&scope, "query", "BLASTN", false, QUERY).unwrap();
        let read = read_in(&scope, "subject", "BLASTN", false, subject);
        let text = String::from_utf8_lossy(subject).into_owned();
        assert_eq!(read.is_err(), refused, "{text}");
        if refused {
            assert!(
                read.unwrap_err().contains("a query and a subject"),
                "{text}"
            );
        }
    }
}

// The order of the readers does not change it, and the second reader's record knows.
#[test]
fn one_record_is_marked_with_the_lower_case_of_either_copy() {
    for (first, second, lowercase) in [
        (&b">q1 t\nACGT\n"[..], &b">q1 t\nACGT\n"[..], false),
        (b">q1 t\nacGT\n", b">q1 t\nACGT\n", true),
        (b">q1 t\nACGT\n", b">q1 t\nacGT\n", true),
    ] {
        let scope = std::sync::Arc::new(IdScope::default());
        let subjects = read_in(&scope, "subject", "BLASTN", false, first).unwrap();
        assert_eq!(ids_of(&subjects[0]).shared, None);
        let queries = read_in(&scope, "query", "BLASTN", false, second).unwrap();
        assert_eq!(ids_of(&queries[0]).shared, Some(SharedRecord { lowercase }));
    }
}

// BLASTX and TBLASTN: a nucleotide and a protein reader never read one record (DECISIONS 19).
#[test]
fn a_nucleotide_and_a_protein_reader_do_not_read_one_record() {
    let scope = std::sync::Arc::new(IdScope::default());
    read_in(&scope, "query", "TBLASTN", true, b">q1 t\nACGT\n").unwrap();
    let error = read_in(&scope, "subject", "TBLASTN", false, b">q1 t\nACGT\n").unwrap_err();
    assert!(error.contains("a query and a subject"), "{error}");
}

// DECISIONS 16: the ids belong to one search, not to the process.
#[test]
fn the_ids_of_a_search_do_not_reach_the_next_search() {
    let first = std::sync::Arc::new(IdScope::default());
    read_in(&first, "query", "BLASTN", false, b">kept1 t\nACGT\n").unwrap();
    read_in(&first, "subject", "BLASTN", false, b">kept1 u\nACGT\n").unwrap_err();
    let second = std::sync::Arc::new(IdScope::default());
    read_in(&second, "subject", "BLASTN", false, b">kept1 u\nACGT\n").unwrap();
}

// Masks of a shared query: the masks of the query, or with `-lcase_masking` the lower case of
// either copy (DECISIONS 19).
#[test]
fn outfmt_0_refuses_the_masks_of_a_shared_record_only() {
    let scope = std::sync::Arc::new(IdScope::default());
    read_in(&scope, "subject", "BLASTN", false, b">same1 t\nACgt\n").unwrap();
    let queries = read_in(
        &scope,
        "query",
        "BLASTN",
        false,
        b">other1 t\nACGT\n>same1 t\nACGT\n",
    )
    .unwrap();
    let refused = |masked: usize, lcase: bool| {
        refuse_shared_masks(&queries, |index| index == masked, lcase, "BLASTN").is_err()
    };
    assert!(!refused(0, false));
    assert!(refused(1, false));
    assert!(refused(0, true));
}

// NCBI BLAST+ 2.17.0 (probe `DI7/ws/pd_sub.py` `subjname`, `noid`): a `-subject` record with a
// local id `Subject_N`, written or generated, is shown as one read without the option.
#[test]
fn subject_local_ids_are_shown_as_records() {
    for (defline, local) in [
        (&b">lcl|Subject_1 x"[..], true),
        (b">Subject_x x", true),
        (b">", true),
        (b">Subjectx x", false),
        (b">plain1 x", false),
        (b">NM_000003.1 x", false),
    ] {
        let mut fasta = defline.to_vec();
        fasta.extend_from_slice(b"\nACGT\n");
        let config = ReaderConfig::subject("BLASTN", false, false).with_parse_deflines(true);
        let mut source = FastaInputSource::from_bytes(&fasta, config);
        let records = read_all(&mut source, &mut |_| Ok(())).unwrap();
        let ids = records[0].believed().unwrap();
        assert_eq!(
            ids.subject_local(),
            local,
            "{}",
            String::from_utf8_lossy(defline)
        );
    }
}
