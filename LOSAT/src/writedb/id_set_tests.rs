use super::*;

#[test]
fn a_label_is_in_the_set_only_when_its_bytes_are() {
    let mut set = LabelSet::default();
    assert!(set.insert(b"REF|NP_000549.1"));
    assert!(!set.insert(b"REF|NP_000549.1"));
    assert!(set.insert(b"REF|NP_000549"));
    assert!(set.insert(b"REF|NP_000549.2"));
    assert!(set.insert(b""));
    assert!(!set.insert(b""));
    assert_eq!(set.len(), 4);
    set.clear();
    assert_eq!(set.len(), 0);
    assert!(set.insert(b"REF|NP_000549.1"));
}

#[test]
fn many_labels_grow_the_table_and_chunks_without_losing_one() {
    let mut set = LabelSet::default();
    for i in 0..200_000u32 {
        assert!(set.insert(format!("LCL|UNIREF50_A{i:010}").as_bytes()));
    }
    // A label longer than a chunk takes a chunk of its own.
    let long = vec![b'X'; CHUNK + 5];
    assert!(set.insert(&long));
    assert!(!set.insert(&long));
    for i in (0..200_000u32).step_by(7) {
        assert!(!set.insert(format!("LCL|UNIREF50_A{i:010}").as_bytes()));
    }
    assert_eq!(set.len(), 200_001);
    // 23-byte labels: the bytes, one length byte, and the table.
    let per_label = (set.heap_bytes() - CHUNK - 5) as f64 / 200_000.0;
    assert!(per_label < 48.0, "{per_label} bytes per label");
}

#[test]
fn a_key_map_keeps_the_last_value_of_a_key() {
    let mut map = KeyMap::default();
    assert!(map.is_empty());
    assert_eq!(map.get(b"A"), None);
    map.set(b"A", 1);
    map.set(b"B", -2);
    map.set(b"A", 3);
    map.set(b"", 7);
    assert_eq!(map.get(b"A"), Some(3));
    assert_eq!(map.get(b"B"), Some(-2));
    assert_eq!(map.get(b""), Some(7));
    assert_eq!(map.get(b"a"), None);
    assert_eq!(map.len(), 3);
    // Through several doublings of the table and more than one chunk.
    for i in 0..100_000i32 {
        map.set(format!("K{i}").as_bytes(), i);
    }
    let long = vec![b'x'; CHUNK + 3];
    map.set(&long, 42);
    for i in (0..100_000i32).step_by(13) {
        assert_eq!(map.get(format!("K{i}").as_bytes()), Some(i));
        map.set(format!("K{i}").as_bytes(), -i);
    }
    for i in (0..100_000i32).step_by(13) {
        assert_eq!(map.get(format!("K{i}").as_bytes()), Some(-i));
    }
    assert_eq!(map.get(&long), Some(42));
    assert_eq!(map.len(), 100_004);
}

#[test]
fn labels_with_one_hash_are_told_apart_by_their_bytes() {
    // Every label lands in a table of 1024 slots; equal low bits must not merge them.
    let mut set = LabelSet::default();
    let labels: Vec<Vec<u8>> = (0..700u32).map(|i| i.to_le_bytes().to_vec()).collect();
    for label in &labels {
        assert!(set.insert(label));
    }
    for label in &labels {
        assert!(!set.insert(label));
    }
    assert_eq!(set.len(), 700);
}
