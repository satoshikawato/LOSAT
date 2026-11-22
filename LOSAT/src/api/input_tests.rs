use super::*;

// Test protein: N-terminal part of E. coli AsnA (RefSeq NP_418200.1), Cys51 changed to Ala.
const PROTEIN: &[u8] =
    b">p1 a protein\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ\n";

#[test]
fn records_carry_the_readers_local_ids_and_titles() {
    let query = read_fasta(Program::Blastp, Role::Query, PROTEIN).unwrap();
    let subject = read_fasta(Program::Blastp, Role::Subject, PROTEIN).unwrap();
    assert_eq!(query.len(), 1);
    assert_eq!(query[0].local_id, "Query_1");
    assert_eq!(subject[0].local_id, "Subject_1");
    assert_eq!(query[0].title(), b"p1 a protein");
    assert_eq!(query[0].seq().len(), 66);
    // An input of white space only has no record; an empty input neither.
    assert!(read_fasta(Program::Blastn, Role::Query, b" \n\t\n")
        .unwrap()
        .is_empty());
    assert!(read_fasta(Program::Blastn, Role::Query, b"")
        .unwrap()
        .is_empty());
}

#[test]
fn the_readers_flags_follow_the_programs_molecules() {
    assert!(!is_protein(Program::Blastn, Role::Query));
    assert!(!is_protein(Program::Blastn, Role::Subject));
    assert!(is_protein(Program::Blastp, Role::Query));
    assert!(is_protein(Program::Blastp, Role::Subject));
    assert!(!is_protein(Program::Blastx, Role::Query));
    assert!(is_protein(Program::Blastx, Role::Subject));
    assert!(is_protein(Program::Tblastn, Role::Query));
    assert!(!is_protein(Program::Tblastn, Role::Subject));
    assert!(!is_protein(Program::Tblastx, Role::Query));
    assert!(!is_protein(Program::Tblastx, Role::Subject));
}

#[test]
fn a_reader_error_is_the_command_lines_message() {
    // A first line that is neither a defline nor plausible data (`!`, `#` and `;` open
    // comment lines): NCBI's CFastaReader error, with the CLI's prefix.
    let error = read_fasta(Program::Blastn, Role::Query, b"xyz{}[]|~`^&*()=+\n");
    match error {
        Err(error @ Error::Input(_)) => {
            let message = error.to_string();
            assert!(
                message.starts_with("LOSAT query error: CFastaReader: Near line 1"),
                "{message}"
            );
            assert!(!message.ends_with('\n'));
        }
        other => panic!("expected a reader error, got {other:?}"),
    }
}

#[test]
fn a_file_reads_as_the_command_line_opens_it() {
    use std::io::Write;
    let directory = crate::utils::test_dir::temp_dir().join(format!(
        "losat-read-input-file-{}",
        crate::utils::test_dir::process_id()
    ));
    std::fs::create_dir_all(&directory).unwrap();
    let plain = directory.join("q.fa");
    std::fs::write(&plain, PROTEIN).unwrap();
    assert_eq!(read_input_file(Role::Query, &plain).unwrap(), PROTEIN);
    let packed = directory.join("q.fa.gz");
    let mut encoder = flate2::write::GzEncoder::new(
        std::fs::File::create(&packed).unwrap(),
        flate2::Compression::default(),
    );
    encoder.write_all(PROTEIN).unwrap();
    encoder.finish().unwrap();
    assert_eq!(read_input_file(Role::Subject, &packed).unwrap(), PROTEIN);
    let missing = directory.join("none.fa");
    let error = read_input_file(Role::Subject, &missing).unwrap_err();
    assert_eq!(
        error,
        Error::Input(
            format!(
                "Command line argument error: Argument \"subject\". File is not accessible:  `{}'",
                missing.display()
            )
            .into()
        )
    );
    // A name that is not UTF-8 keeps its bytes, as on the command line (line LIB-B).
    #[cfg(unix)]
    {
        use std::os::unix::ffi::OsStrExt;
        let missing = directory.join(std::ffi::OsStr::from_bytes(b"none\xE9.fa"));
        let error = read_input_file(Role::Query, &missing).unwrap_err();
        let expected = [
            b"Command line argument error: Argument \"query\". File is not accessible:  `"
                .as_slice(),
            missing.as_os_str().as_bytes(),
            b"'",
        ]
        .concat();
        assert_eq!(error, Error::Input(expected));
    }
    std::fs::remove_dir_all(&directory).unwrap();
}

#[test]
fn inputs_default_names_are_query_and_subject() {
    let records = read_fasta(Program::Blastp, Role::Query, PROTEIN).unwrap();
    let inputs = Inputs::new(&records, &records);
    assert_eq!(
        (inputs.query_name, inputs.subject_name),
        ("query", "subject")
    );
    let named = inputs.named("q.faa", "s.faa");
    assert_eq!((named.query_name, named.subject_name), ("q.faa", "s.faa"));
}
