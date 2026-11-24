use super::*;
use crate::dbindex::stream::Bases;

/// A chunk whose end passes 2^32 (here with chunks of 2^32 - 10 letters, 5 apart) is
/// refused: NCBI's `TSeqPos` end wraps and its index reads past the sequence (line RG
/// DECISIONS KC-11d; makembindex 2.17.0 on ws, remaining-gaps-20261112/K/kc/p11d/p3 l93n,
/// l94n: SIGSEGV).
#[test]
fn a_chunk_whose_end_wraps_is_refused() {
    let options = Options {
        chunk_size: (1 << 32) - 10,
        chunk_overlap: (1 << 32) - 15,
        ..Options::default_options()
    };
    let mut map = SubjectMap::new(&options);
    map.new_sequence_init(
        Some(Rc::new(SeqData {
            bases: Bases::of_letters(&b"ACGT".repeat(5)),
            masks: Vec::new(),
            id: None,
            texts: Vec::new(),
        })),
        0,
    );
    assert!(matches!(map.add_sequence_chunk(), Ok((true, false))));
    assert!(matches!(map.add_sequence_chunk(), Ok((true, false))));
    let Err(StreamError::Refused(what)) = map.add_sequence_chunk() else {
        panic!("not refused")
    };
    assert_eq!(
        what,
        "a sequence of 11 letters or more, whose last chunk of 4,294,967,286 letters NCBI's \
             makembindex reads past the sequence's end (its end wraps at 32 bits) is not \
             supported by LOSAT's MAKEMBINDEX"
    );
}
