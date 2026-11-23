//! The list of blast_formatter's refusals: the formats it writes are not refused, every other
//! format of NCBI's blast_formatter is, and the texts are LOSAT's (no NCBI footer).

use super::*;

#[test]
fn the_written_formats_are_not_refused() {
    for number in WRITTEN_FORMATS {
        assert_eq!(output_format(number), None, "{number}");
    }
    for number in NCBI_FORMATS {
        assert_eq!(
            output_format(number).is_none(),
            WRITTEN_FORMATS.contains(&number),
            "{number}"
        );
    }
}

#[test]
fn the_documented_refusals_are_one_of_each_kind() {
    let documented = documented();
    assert_eq!(documented[0], Refusal::Rid);
    let formats: Vec<i32> = documented
        .iter()
        .filter_map(|refusal| match refusal {
            Refusal::OutputFormat(number) => Some(*number),
            _ => None,
        })
        .collect();
    assert_eq!(
        formats,
        [0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 13, 14, 15, 16, 17, 18, 20]
    );
    assert!(documented.contains(&Refusal::BinaryArchive));
    assert!(documented.contains(&Refusal::XmlArchive));
    assert!(documented.contains(&Refusal::FatalMessage));
    for refusal in &documented {
        assert!(!refusal.manual_words().is_empty(), "{refusal:?}");
    }
}

#[test]
fn the_texts_are_losat_s() {
    assert_eq!(
        Refusal::OutputFormat(6).reason(),
        "output format 6 is not supported by LOSAT's BLAST_FORMATTER (it writes 8, 9 and 12)"
    );
    assert_eq!(
        Refusal::Rid.reason(),
        "-rid (the results of a search on NCBI's servers) is not supported by LOSAT's \
         BLAST_FORMATTER"
    );
    assert_eq!(
        Refusal::Database("nt".to_string()).reason(),
        "the database 'nt' of the BLAST archive cannot be opened here; NCBI BLAST+ then asks \
         NCBI's servers for its sequences, which is not supported by LOSAT's BLAST_FORMATTER"
    );
    // A refusal travels as an error to `catch_all`, which writes its reason.
    let error: anyhow::Error = Refusal::FatalMessage.into();
    assert_eq!(format!("{error:#}"), Refusal::FatalMessage.reason());
    for refusal in documented() {
        let reason = refusal.reason();
        assert!(reason.contains("LOSAT's BLAST_FORMATTER"), "{reason}");
        assert!(!reason.contains("user manual"), "{reason}");
    }
}
