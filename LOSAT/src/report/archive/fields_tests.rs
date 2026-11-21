use super::*;

#[test]
fn value_types_are_numbered_as_ncbi_e_choice() {
    assert_eq!(ValueType::BigInteger.which(), 1);
    assert_eq!(ValueType::Cutoff.which(), 4);
    assert_eq!(ValueType::Integer.which(), 5);
    assert_eq!(ValueType::String.which(), 12);
    assert_eq!(ValueType::IntegerList.which(), 17);
    assert_eq!(ValueType::QueryMask.which(), 27);
    for (tag, ty) in ValueType::ALL.iter().enumerate() {
        assert_eq!(usize::from(ty.tag()), tag);
        assert_eq!(ValueType::of_tag(tag as u16), Some(*ty));
    }
    assert_eq!(ValueType::of_tag(27), None);
}

#[test]
fn match_needs_name_and_type() {
    assert!(CULLING_LIMIT.matches(b"Culling", ValueType::Integer));
    assert!(!CULLING_LIMIT.matches(b"Culling", ValueType::Boolean));
    assert!(CULLING.matches(b"Culling", ValueType::Boolean));
    assert!(!EVALUE_THRESHOLD.matches(b"EvalueThreshold", ValueType::Real));
    assert!(!WORD_SIZE.matches(b"wordsize", ValueType::Integer));
}
