use super::*;
use crate::utils::test_dir::TempDir;

// The file makeblastdb 2.17.0 wrote for the three proteins of survey A's
// f1_tiny.fa (`-title t -out db`, version 5).
#[test]
fn the_json_layout_is_ncbi_s() {
    let metadata = DbMetadata {
        dbname: b"db".to_vec(),
        dbtype: "Protein",
        db_version: 5,
        description: b"t".to_vec(),
        number_of_letters: 120,
        number_of_sequences: 3,
        last_updated: "2026-10-29T20:32:00".to_string(),
        number_of_volumes: 1,
        number_of_taxids: None,
        bytes_total: 37723,
        bytes_to_cache: 228,
        files: [
            "db.pdb", "db.phr", "db.pin", "db.pot", "db.psq", "db.ptf", "db.pto",
        ]
        .map(|name| name.as_bytes().to_vec())
        .to_vec(),
    };
    let expected =
        "{\n  \"version\": \"1.2\",\n  \"dbname\": \"db\",\n  \"dbtype\": \"Protein\",\n  \
             \"db-version\": 5,\n  \"description\": \"t\",\n  \"number-of-letters\": 120,\n  \
             \"number-of-sequences\": 3,\n  \"last-updated\": \"2026-10-29T20:32:00\",\n  \
             \"number-of-volumes\": 1,\n  \"bytes-total\": 37723,\n  \"bytes-to-cache\": 228,\n  \
             \"files\": [\n    \"db.pdb\",\n    \"db.phr\",\n    \"db.pin\",\n    \"db.pot\",\n    \
             \"db.psq\",\n    \"db.ptf\",\n    \"db.pto\"\n  ]\n}\n";
    assert_eq!(String::from_utf8(metadata.to_json()).unwrap(), expected);
}

// makeblastdb 2.17.0 with `-title` = `q"b\s<TAB>t<0x01>c é end ` (é in UTF-8).
#[test]
fn strings_are_escaped_as_ncbi_s_json_writer_does() {
    let escaped = |value: &[u8]| {
        let mut out = JsonText::default();
        let Ok(()) = write_json_string(&mut out, value);
        out.bytes
    };
    assert_eq!(
        escaped("q\"b\\s\tt\u{1}c é end ".as_bytes()),
        b"\"q\\\"b\\\\s\\u0009t\\u0001c \xc3\xa9 end \"".to_vec()
    );
    assert_eq!(
        escaped(b"a\x1bb\x7f\xff\x00c"),
        b"\"a\\u001Bb\x7f\xff\"".to_vec()
    );
}

/// The path of the database `db` in `dir`, as the text the metadata writer takes.
fn db_path(dir: &Path) -> String {
    dir.join("db").to_string_lossy().into_owned()
}

/// The metadata of a protein database of `volumes` volumes of four files each, with a
/// description of `pad` bytes.
fn metadata(volumes: usize, pad: usize) -> DbMetadata {
    let mut files: Vec<Vec<u8>> = (0..volumes)
        .flat_map(|v| ["phr", "pin", "pog", "psq"].map(|e| format!("db.{v:03}.{e}").into_bytes()))
        .collect();
    files.extend(
        ["db.pal", "db.pdb", "db.pos", "db.pot", "db.ptf", "db.pto"]
            .map(|name| name.as_bytes().to_vec()),
    );
    DbMetadata {
        dbname: b"db".to_vec(),
        dbtype: "Protein",
        db_version: 5,
        description: vec![b'x'; pad],
        number_of_letters: 7800,
        number_of_sequences: 130,
        last_updated: "2026-11-12T12:56:00".to_string(),
        number_of_volumes: volumes as i32,
        number_of_taxids: None,
        bytes_total: 412_004,
        bytes_to_cache: 61_380,
        files,
    }
}

/// The metadata file's text (when it is a regular file: `/dev/full` reads without end), its
/// messages and the exception that ends makeblastdb.
fn written(path: &str, metadata: &DbMetadata) -> (Option<Vec<u8>>, String, Option<&'static str>) {
    let mut err = Vec::new();
    let escaped = write_metadata_file(path.as_bytes(), true, metadata, &mut err);
    let pjs = format!("{path}.pjs");
    let text = std::fs::metadata(&pjs)
        .is_ok_and(|meta| meta.is_file())
        .then(|| std::fs::read(&pjs).ok())
        .flatten();
    (text, String::from_utf8(err).unwrap(), escaped)
}

// The object is written through the serializer's buffer: what reaches the file is the text
// of `to_json`, whatever its length.
#[test]
fn the_metadata_file_holds_the_object_and_a_line_feed() {
    for volumes in [1, 130] {
        let tmp = TempDir::unique(&format!("metadata-ok{volumes}"));
        let metadata = metadata(volumes, 3);
        let (text, err, escaped) = written(&db_path(&tmp), &metadata);
        assert_eq!(text, Some(metadata.to_json()), "{volumes}");
        assert_eq!((err.as_str(), escaped), ("", None));
    }
}

// A file that cannot be opened (a directory of its name): the serializer's `Close` finds the
// write fault at the end of what it holds. The text up to 4096 bytes is written at
// `Flush`; a longer one fails at its first write, when the 4097th byte does not fit.
// makeblastdb 2.17.0: `byte 479` (b_p5_pjs) and `byte 4096` (b_m_pjs), exit status 0.
#[test]
fn a_metadata_file_that_cannot_be_opened_reports_a_write_fault() {
    for (volumes, at) in [(1, None), (130, Some(4096))] {
        let tmp = TempDir::unique(&format!("metadata-blocked{volumes}"));
        std::fs::create_dir(format!("{}.pjs", db_path(&tmp))).unwrap();
        let metadata = metadata(volumes, 3);
        let at = at.unwrap_or(metadata.to_json().len() - 1);
        let (_, err, escaped) = written(&db_path(&tmp), &metadata);
        assert_eq!(
            err,
            format!(
                "Error: (802.5) [makeblastdb] CObjectOStream: error at byte {at}: stack is \
                     empty: write fault\n"
            )
        );
        assert_eq!(escaped, None);
    }
}

// The buffer is written when a call does not fit: a line feed at byte 4092 fits, its
// indentation of 4 spaces does not, so the first write holds 4093 bytes.
#[test]
fn the_buffer_is_written_where_a_call_does_not_fit() {
    let plain = metadata(130, 0).to_json();
    let eol = plain[..4097]
        .windows(5)
        .rposition(|w| w == b"\n    ")
        .unwrap();
    let metadata = metadata(130, 4092 - eol);
    assert_eq!(&metadata.to_json()[4092..4097], b"\n    ");
    let tmp = TempDir::unique("metadata-straddle");
    std::fs::create_dir(format!("{}.pjs", db_path(&tmp))).unwrap();
    let (_, err, _) = written(&db_path(&tmp), &metadata);
    assert_eq!(
        err,
        "Error: (802.5) [makeblastdb] CObjectOStream: error at byte 4093: stack is empty: \
             write fault\n"
    );
}

// A device that takes no byte: a text of up to 4096 bytes sits in the stream's buffer of
// 8191 bytes until `out.flush()` in `COStreamBuffer::Flush` fails, `eFlush` leaves
// makeblastdb's block, and the serializer's `Close` fails again and reports (f_p5_pjs,
// exit status 6). A longer text gets its first 4096 bytes into the stream's buffer; the
// next 4096 do not fit the 4095 left and their write fails (f_m_pjs: `byte 8192`, exit
// status 0).
#[cfg(target_os = "linux")]
#[test]
fn a_metadata_file_on_a_full_device_cannot_be_written() {
    let tmp = TempDir::unique("metadata-full");
    std::os::unix::fs::symlink("/dev/full", format!("{}.pjs", db_path(&tmp))).unwrap();
    let small = metadata(1, 3);
    let at = small.to_json().len() - 1;
    let (_, err, escaped) = written(&db_path(&tmp), &small);
    // NCBI's text in LOSAT's words (`crate::wording::losat_text`).
    assert_eq!(
        err,
        crate::wording::losat_text(&format!(
            "Error: (802.5) [makeblastdb] CObjectOStream: error at byte {at}: stack is empty: \
                 cannot close output stream\n\
                 Error: (CIOException::eFlush) [makeblastdb] COStreamBuffer::Flush: failed\n\
                 Error: (CSerialException::eIoError) [makeblastdb] byte {at}: cannot close output \
                 stream\n\
                 Error: (802.4) [makeblastdb] Cannot close serializing output stream \
                 (CSerialException::eIoError) byte {at}: cannot close output stream\n"
        ))
    );
    assert_eq!(escaped, Some("COStreamBuffer::Flush: failed"));
    let (_, err, escaped) = written(&db_path(&tmp), &metadata(130, 3));
    assert_eq!(
        err,
        "Error: (802.5) [makeblastdb] CObjectOStream: error at byte 8192: stack is empty: \
             write fault\n"
    );
    assert_eq!(escaped, None);
}
