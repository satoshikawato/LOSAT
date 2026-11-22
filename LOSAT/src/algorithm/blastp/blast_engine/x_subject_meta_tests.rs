use super::super::tabular_subject_ids;
use super::*;

/// Subject records whose IDs take every branch of `tabular_subject_id`: a first word, the local
/// ID of an empty title, `Subject_` first words that `GenerateDefline` keeps, rewrites or
/// cannot decode (`\x81`, no encoding NCBI guesses: `None`), and an empty record.
fn records() -> Vec<FastaRecord> {
    let titles: [&[u8]; 9] = [
        b"sp|P1|A first protein",
        b"",
        b"Subject_1 desc",
        b"Subject_\x81 desc",
        b"   ",
        b"Subject_\xe9 desc",
        b"lcl|Subject_7 foo",
        b"Subject_\x81",
        b"x",
    ];
    titles
        .iter()
        .enumerate()
        .map(|(index, title)| {
            FastaRecord::new(format!("Subject_{}", index + 1), title, b"MKVLAAGIVA")
        })
        .collect()
}

// The stored values are those that every block makes: each subject's shown ID (none where it
// fails) and the titles.
#[test]
fn stored_values_are_the_per_block_values() {
    let records = records();
    let meta = XSubjectMeta::new(&records);
    assert_eq!(meta.shown().len(), records.len());
    for (index, record) in records.iter().enumerate() {
        let direct = crate::report::defline::tabular_subject_id(
            record.local_id.as_bytes(),
            &record.title,
            true,
        )
        .ok();
        assert_eq!(meta.shown()[index].as_deref(), direct.as_deref(), "{index}");
        assert_eq!(&*meta.titles()[index], record.title(), "{index}");
    }
    // Two of the records cannot be written (`Subject_\x81`), the others can.
    let failing: Vec<usize> = (0..records.len())
        .filter(|&index| meta.shown()[index].is_none())
        .collect();
    assert_eq!(failing, vec![3, 7]);
    // No reported subject: `tabular_subject_ids` shows every ID, none where it fails.
    assert_eq!(meta.shown(), tabular_subject_ids(&records));
    assert_eq!(meta.shown()[0].as_deref(), Some(&b"sp|P1|A"[..]));
}
