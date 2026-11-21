use super::*;

fn id(text: &str) -> SetKey {
    use crate::seqid::accguide::parse_flags::{ANY_LOCAL, PARTIAL_OK, RAW_TEXT};
    let ids = crate::seqid::parse::parse_ids(
        text.as_bytes(),
        PARTIAL_OK | ANY_LOCAL | RAW_TEXT,
        &mut |_| {},
    )
    .unwrap();
    SetKey::Id(ids.into_iter().next().unwrap().unwrap())
}

// Review 1 of line DBID pssmrows, M1: only a shown subject that a previous key can be one
// handle with is refused (blast_format.cpp:541-553, cited at `refuse_unsure_shown`).
#[test]
fn only_a_shown_unsure_subject_is_refused() {
    let found = FoundBefore {
        found: vec![true, false, false],
        unsure: [(2, (id("gb|ABC12.1|"), id("emb|ABC12.1|")))]
            .into_iter()
            .collect(),
    };
    found.refuse_unsure_shown([0, 1]).unwrap();
    let error = found.refuse_unsure_shown([0, 2]).unwrap_err().to_string();
    assert!(error.contains("gb|ABC12.1| and emb|ABC12.1|"), "{error}");
}
