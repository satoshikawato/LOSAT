//! `XmlRefusal::of_taxonomy`: only a refusal of the taxonomy lookups drops a held document;
//! a read error passes through, and the run writes what it holds (`HeldState::drop`).
use super::*;
use crate::cli::NativeError;
use crate::seqdb::SeqDbError;
use std::path::PathBuf;

fn through(error: SeqDbError) -> anyhow::Error {
    XmlRefusal::of_taxonomy(crate::algorithm::common::database::taxdb_error_for_tests(
        error,
    ))
}

#[test]
fn a_refused_record_or_name_is_a_refusal_with_its_message() {
    let error = through(SeqDbError::Unsupported {
        path: PathBuf::from("taxdb.btd"),
        message: "the taxonomy record of taxid 9606 reaches past the mapped end".into(),
    });
    assert!(error.is::<XmlRefusal>());
    assert!(error.to_string().contains("taxid 9606"), "{error}");
}

#[test]
fn a_read_error_is_not_a_refusal() {
    let io = SeqDbError::Io {
        path: PathBuf::from("taxdb.btd"),
        source: std::io::Error::other("Input/output error"),
    };
    let error = through(io);
    assert!(!error.is::<XmlRefusal>(), "{error}");
    let format = SeqDbError::Format {
        path: PathBuf::from("db.phr"),
        message: "byte range 0..9 passes the end of the file (4 bytes)".into(),
    };
    let error = through(format);
    assert!(!error.is::<XmlRefusal>(), "{error}");
    // An error of NCBI's own (exit status 2) is passed on whole.
    let native = anyhow::Error::new(NativeError {
        exit: 2,
        message: "BLAST Database error: x\n".into(),
    });
    assert!(XmlRefusal::of_taxonomy(native).is::<NativeError>());
}
