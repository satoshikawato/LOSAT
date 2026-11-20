//! The set of Seq-id labels of one volume, for the duplicate check of
//! `CWriteDB_Volume::WriteSequence` (`m_IdSet`,
//! `c++/src/objtools/blast/seqdb_writer/writedb_volume.cpp:228-253`).
//!
//! NCBI keeps a `set<string>`: a tree node and a `std::string` for every label, a heap
//! buffer more for a label of more than 15 bytes (about 112 bytes for a label such as
//! `LCL|UNIREF50_A0A5A9P0L4` with glibc and libstdc++). Here the labels are packed into an
//! arena of 1 MiB chunks (a LEB128 length and the bytes) and found through an open-addressed
//! table of 64-bit slots (the arena position and 16 bits of the label's hash). A label is in
//! the set only when its bytes equal a stored label's bytes, so the set answers as NCBI's
//! does, without false positives or negatives; the hash only chooses where to look.
//! Per label: its bytes plus one or two length bytes, and 8 bytes per slot of a table at
//! most three quarters full (10.7 to 16 bytes per label after the table doubles).
//!
//! The same table, with a 4-byte value after each key, holds the `-taxid_map` map
//! (`KeyMap`, `CTaxIdSet::m_TaxIdMap`, a `map<string, TTaxId>` in NCBI).

/// The bytes of one arena chunk.
const CHUNK_BITS: u32 = 20;
const CHUNK: usize = 1 << CHUNK_BITS;
/// The slots of a table before its first key.
const MIN_SLOTS: usize = 1 << 10;
/// The arena position (+ 1) in a slot.
const POS_MASK: u64 = (1 << 48) - 1;

/// Byte strings packed into an arena of 1 MiB chunks (a LEB128 length, the bytes, then
/// `VALUE` bytes of a value) and found through an open-addressed table of 64-bit slots
/// (the arena position and 16 bits of the key's hash). Equality is decided by the bytes.
#[derive(Debug, Default, Clone)]
struct ArenaTable<const VALUE: usize> {
    /// The keys, never split across chunks.
    chunks: Vec<Vec<u8>>,
    /// The chunk being filled.
    current: usize,
    /// 0 for an empty slot; else the key's position + 1 (low 48 bits) and the top 16
    /// bits of its hash.
    slots: Vec<u64>,
    len: usize,
}

/// A 64-bit hash of a label (64-bit words through a multiply-rotate step, then
/// MurmurHash3's finalizer): it spreads the labels over the slots; equality is decided by
/// the bytes.
fn label_hash(bytes: &[u8]) -> u64 {
    let mut h: u64 = (bytes.len() as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    let mut words = bytes.chunks_exact(8);
    for word in &mut words {
        let w = u64::from_le_bytes(word.try_into().expect("8 bytes"));
        h = (h ^ w).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(29);
    }
    let rest = words.remainder();
    if !rest.is_empty() {
        let mut last = [0u8; 8];
        last[..rest.len()].copy_from_slice(rest);
        h = (h ^ u64::from_le_bytes(last))
            .wrapping_mul(0x9E37_79B9_7F4A_7C15)
            .rotate_left(29);
    }
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
    h ^ (h >> 33)
}

impl<const VALUE: usize> ArenaTable<VALUE> {
    /// Empties the table; its memory is kept.
    fn clear(&mut self) {
        for chunk in &mut self.chunks {
            chunk.clear();
        }
        self.current = 0;
        self.slots.iter_mut().for_each(|slot| *slot = 0);
        self.len = 0;
    }

    /// The chunk, the start and the length of the key stored at `pos`.
    fn locate(&self, pos: u64) -> (usize, usize, usize) {
        let index = (pos >> CHUNK_BITS) as usize;
        let chunk = &self.chunks[index];
        let mut at = (pos as usize) & (CHUNK - 1);
        let mut length = 0usize;
        let mut shift = 0;
        loop {
            let byte = chunk[at];
            at += 1;
            length |= usize::from(byte & 0x7f) << shift;
            if byte & 0x80 == 0 {
                break;
            }
            shift += 7;
        }
        (index, at, length)
    }

    /// The bytes of the key stored at `pos`.
    fn key_at(&self, pos: u64) -> &[u8] {
        let (index, at, length) = self.locate(pos);
        &self.chunks[index][at..at + length]
    }

    /// The value bytes of the key stored at `pos`.
    fn value_at(&self, pos: u64) -> &[u8] {
        let (index, at, length) = self.locate(pos);
        &self.chunks[index][at + length..at + length + VALUE]
    }

    fn value_at_mut(&mut self, pos: u64) -> &mut [u8] {
        let (index, at, length) = self.locate(pos);
        &mut self.chunks[index][at + length..at + length + VALUE]
    }

    /// Stores `key` and `value` in the arena and returns the position.
    fn push(&mut self, key: &[u8], value: &[u8; VALUE]) -> u64 {
        let mut head = [0u8; 10];
        let mut head_len = 0;
        let mut length = key.len();
        loop {
            let byte = (length & 0x7f) as u8;
            length >>= 7;
            if length == 0 {
                head[head_len] = byte;
                head_len += 1;
                break;
            }
            head[head_len] = byte | 0x80;
            head_len += 1;
        }
        let size = head_len + key.len() + VALUE;
        // A chunk takes a key that ends within its 1 MiB, or any key while it is empty:
        // every position is then below 2^20 within its chunk.
        if self.chunks.is_empty() {
            self.chunks.push(Vec::with_capacity(CHUNK.max(size)));
        }
        if !self.chunks[self.current].is_empty() && self.chunks[self.current].len() + size > CHUNK {
            self.current += 1;
            if self.current == self.chunks.len() {
                self.chunks.push(Vec::with_capacity(CHUNK.max(size)));
            }
        }
        let chunk = &mut self.chunks[self.current];
        let pos = ((self.current as u64) << CHUNK_BITS) | chunk.len() as u64;
        chunk.extend_from_slice(&head[..head_len]);
        chunk.extend_from_slice(key);
        chunk.extend_from_slice(value);
        pos
    }

    /// Doubles the table (or makes the first one) and puts every key in its new slot.
    fn grow(&mut self) {
        let size = (self.slots.len() * 2).max(MIN_SLOTS);
        let old = std::mem::replace(&mut self.slots, vec![0; size]);
        let mask = size - 1;
        for slot in old.into_iter().filter(|&slot| slot != 0) {
            let pos = (slot & POS_MASK) - 1;
            let mut at = (label_hash(self.key_at(pos)) as usize) & mask;
            while self.slots[at] != 0 {
                at = (at + 1) & mask;
            }
            self.slots[at] = slot;
        }
    }

    /// The position of `key`, or the empty slot where it goes and the tag of its hash.
    fn find(&self, key: &[u8]) -> Result<u64, (usize, u64)> {
        let hash = label_hash(key);
        let tag = hash >> 48;
        if self.slots.is_empty() {
            return Err((0, tag));
        }
        let mask = self.slots.len() - 1;
        let mut at = (hash as usize) & mask;
        loop {
            let slot = self.slots[at];
            if slot == 0 {
                return Err((at, tag));
            }
            if slot >> 48 == tag && self.key_at((slot & POS_MASK) - 1) == key {
                return Ok((slot & POS_MASK) - 1);
            }
            at = (at + 1) & mask;
        }
    }

    /// Adds `key` with `value` when it is not in the table: `Ok(position)` of the new key,
    /// or `Err(position)` of the key already there (its value unchanged).
    fn insert(&mut self, key: &[u8], value: &[u8; VALUE]) -> Result<u64, u64> {
        if (self.len + 1) * 4 > self.slots.len() * 3 {
            self.grow();
        }
        let (at, tag) = match self.find(key) {
            Ok(pos) => return Err(pos),
            Err(empty) => empty,
        };
        let pos = self.push(key, value);
        self.slots[at] = (tag << 48) | (pos + 1);
        self.len += 1;
        Ok(pos)
    }

    /// The bytes this table holds (arena capacity and slots), for the memory figures.
    #[cfg(test)]
    fn heap_bytes(&self) -> usize {
        self.chunks.iter().map(Vec::capacity).sum::<usize>() + self.slots.capacity() * 8
    }
}

/// `set<string> m_IdSet` of one volume.
#[derive(Debug, Default)]
pub(crate) struct LabelSet {
    table: ArenaTable<0>,
}

impl LabelSet {
    /// The number of labels (`m_IdSet.size()`).
    pub(crate) fn len(&self) -> usize {
        self.table.len
    }

    /// `m_IdSet.clear()`: the set is emptied; its memory is kept for the next volume.
    pub(crate) fn clear(&mut self) {
        self.table.clear();
    }

    /// `m_IdSet.insert(label).second`: `true` when the label was not in the set and is
    /// now.
    pub(crate) fn insert(&mut self, label: &[u8]) -> bool {
        self.table.insert(label, &[]).is_ok()
    }

    /// The bytes this set holds (arena capacity and table), for the memory figures.
    #[cfg(test)]
    pub(crate) fn heap_bytes(&self) -> usize {
        self.table.heap_bytes()
    }
}

/// A `map<string, int>` in the same form (`CTaxIdSet::m_TaxIdMap`): the key's bytes, one or
/// two length bytes and the 4-byte value in the arena, 8 bytes per slot of a table at most
/// three quarters full. While the table doubles, the old and the new table are held (16 to
/// 32 bytes of slots per key at that moment).
#[derive(Debug, Default, Clone)]
pub(crate) struct KeyMap {
    table: ArenaTable<4>,
}

impl KeyMap {
    /// `m_TaxIdMap.empty()`.
    pub(crate) fn is_empty(&self) -> bool {
        self.table.len == 0
    }

    /// The number of keys.
    #[cfg(test)]
    pub(crate) fn len(&self) -> usize {
        self.table.len
    }

    /// `m_TaxIdMap[key] = value`: a key already there takes the new value.
    pub(crate) fn set(&mut self, key: &[u8], value: i32) {
        if let Err(pos) = self.table.insert(key, &value.to_le_bytes()) {
            self.table
                .value_at_mut(pos)
                .copy_from_slice(&value.to_le_bytes());
        }
    }

    /// `m_TaxIdMap.find(key)`.
    pub(crate) fn get(&self, key: &[u8]) -> Option<i32> {
        let pos = self.table.find(key).ok()?;
        Some(i32::from_le_bytes(
            self.table.value_at(pos).try_into().expect("4 bytes"),
        ))
    }

    /// The bytes this map holds (arena capacity and table), for the memory figures.
    #[cfg(test)]
    pub(crate) fn heap_bytes(&self) -> usize {
        self.table.heap_bytes()
    }
}

#[cfg(test)]
mod tests {
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
}
