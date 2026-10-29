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

/// The bytes of one arena chunk.
const CHUNK_BITS: u32 = 20;
const CHUNK: usize = 1 << CHUNK_BITS;
/// The slots of a table before its first label.
const MIN_SLOTS: usize = 1 << 10;

/// `set<string> m_IdSet` of one volume.
#[derive(Debug, Default)]
pub(crate) struct LabelSet {
    /// The labels, each as LEB128 length + bytes, never split across chunks.
    chunks: Vec<Vec<u8>>,
    /// The chunk being filled.
    current: usize,
    /// 0 for an empty slot; else the label's position + 1 (low 48 bits) and the top 16
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

impl LabelSet {
    /// The number of labels (`m_IdSet.size()`).
    pub(crate) fn len(&self) -> usize {
        self.len
    }

    /// `m_IdSet.clear()`: the set is emptied; its memory is kept for the next volume.
    pub(crate) fn clear(&mut self) {
        for chunk in &mut self.chunks {
            chunk.clear();
        }
        self.current = 0;
        self.slots.iter_mut().for_each(|slot| *slot = 0);
        self.len = 0;
    }

    /// The bytes of the label stored at `pos`.
    fn label_at(&self, pos: u64) -> &[u8] {
        let chunk = &self.chunks[(pos >> CHUNK_BITS) as usize];
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
        &chunk[at..at + length]
    }

    /// Stores `label` in the arena and returns its position.
    fn push_label(&mut self, label: &[u8]) -> u64 {
        let mut head = [0u8; 10];
        let mut head_len = 0;
        let mut length = label.len();
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
        let size = head_len + label.len();
        // A chunk takes a label that ends within its 1 MiB, or any label while it is
        // empty: every position is then below 2^20 within its chunk.
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
        chunk.extend_from_slice(label);
        pos
    }

    /// Doubles the table (or makes the first one) and puts every label in its new slot.
    fn grow(&mut self) {
        let size = (self.slots.len() * 2).max(MIN_SLOTS);
        let old = std::mem::replace(&mut self.slots, vec![0; size]);
        let mask = size - 1;
        for slot in old.into_iter().filter(|&slot| slot != 0) {
            let pos = (slot & ((1 << 48) - 1)) - 1;
            let mut at = (label_hash(self.label_at(pos)) as usize) & mask;
            while self.slots[at] != 0 {
                at = (at + 1) & mask;
            }
            self.slots[at] = slot;
        }
    }

    /// `m_IdSet.insert(label).second`: `true` when the label was not in the set and is
    /// now.
    pub(crate) fn insert(&mut self, label: &[u8]) -> bool {
        if (self.len + 1) * 4 > self.slots.len() * 3 {
            self.grow();
        }
        let hash = label_hash(label);
        let tag = hash >> 48;
        let mask = self.slots.len() - 1;
        let mut at = (hash as usize) & mask;
        loop {
            let slot = self.slots[at];
            if slot == 0 {
                break;
            }
            if slot >> 48 == tag && self.label_at((slot & ((1 << 48) - 1)) - 1) == label {
                return false;
            }
            at = (at + 1) & mask;
        }
        let pos = self.push_label(label);
        self.slots[at] = (tag << 48) | (pos + 1);
        self.len += 1;
        true
    }

    /// The bytes this set holds (arena capacity and table), for the memory figures.
    #[cfg(test)]
    pub(crate) fn heap_bytes(&self) -> usize {
        self.chunks.iter().map(Vec::capacity).sum::<usize>() + self.slots.capacity() * 8
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
