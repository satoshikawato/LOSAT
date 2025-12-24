use rustc_hash::FxHashMap;

/// Type alias for nucleotide k-mer lookup table
/// Maps encoded k-mer to list of (query_index, position) pairs
pub type NuclKmerLookup = FxHashMap<u64, Vec<(u32, u32)>>;

/// Encode a nucleotide k-mer to a 2-bit encoded integer
///
/// Encoding: T/U=0, C=1, A=2, G=3
/// Returns None if the k-mer contains ambiguous bases
#[inline]
pub fn encode_nucl_kmer(seq: &[u8], start: usize, k: usize) -> Option<u64> {
    if start + k > seq.len() {
        return None;
    }

    let mut code: u64 = 0;
    for i in 0..k {
        let base = seq[start + i];
        let bits = match base {
            b'T' | b't' | b'U' | b'u' => 0,
            b'C' | b'c' => 1,
            b'A' | b'a' => 2,
            b'G' | b'g' => 3,
            _ => return None, // Ambiguous base
        };
        code = (code << 2) | bits;
    }
    Some(code)
}

/// Decode a 2-bit encoded k-mer back to nucleotide sequence
pub fn decode_nucl_kmer(code: u64, k: usize) -> Vec<u8> {
    let mut seq = vec![0u8; k];
    let mut c = code;
    for i in (0..k).rev() {
        let bits = (c & 3) as u8;
        seq[i] = match bits {
            0 => b'T',
            1 => b'C',
            2 => b'A',
            3 => b'G',
            _ => unreachable!(),
        };
        c >>= 2;
    }
    seq
}

/// Build a k-mer lookup table from query sequences
///
/// Returns a hash map mapping encoded k-mers to their positions
pub fn build_nucl_lookup<S: AsRef<[u8]>>(
    queries: &[(String, S)],
    word_size: usize,
) -> NuclKmerLookup {
    let mut lookup: NuclKmerLookup = FxHashMap::default();

    for (query_idx, (_, seq)) in queries.iter().enumerate() {
        let seq = seq.as_ref();
        if seq.len() < word_size {
            continue;
        }

        for pos in 0..=(seq.len() - word_size) {
            if let Some(code) = encode_nucl_kmer(seq, pos, word_size) {
                lookup
                    .entry(code)
                    .or_default()
                    .push((query_idx as u32, pos as u32));
            }
        }
    }

    lookup
}

/// Build a k-mer lookup table from bio::io::fasta::Record sequences
pub fn build_nucl_lookup_from_records(
    records: &[bio::io::fasta::Record],
    word_size: usize,
) -> NuclKmerLookup {
    let mut lookup: NuclKmerLookup = FxHashMap::default();

    for (query_idx, record) in records.iter().enumerate() {
        let seq = record.seq();
        if seq.len() < word_size {
            continue;
        }

        for pos in 0..=(seq.len() - word_size) {
            if let Some(code) = encode_nucl_kmer(seq, pos, word_size) {
                lookup
                    .entry(code)
                    .or_default()
                    .push((query_idx as u32, pos as u32));
            }
        }
    }

    lookup
}

/// Configuration for nucleotide word finding
#[derive(Debug, Clone, Copy)]
pub struct NuclWordFinderConfig {
    /// Word size (k-mer length)
    pub word_size: usize,
    /// Maximum hits per k-mer (for filtering repetitive k-mers)
    pub max_hits_per_kmer: Option<usize>,
}

impl Default for NuclWordFinderConfig {
    fn default() -> Self {
        Self {
            word_size: 11,
            max_hits_per_kmer: None,
        }
    }
}

impl NuclWordFinderConfig {
    pub fn megablast() -> Self {
        Self {
            word_size: 28,
            max_hits_per_kmer: None,
        }
    }

    pub fn blastn() -> Self {
        Self {
            word_size: 11,
            max_hits_per_kmer: None,
        }
    }
}

/// Nucleotide word finder for seed detection
pub struct NuclWordFinder {
    lookup: NuclKmerLookup,
    config: NuclWordFinderConfig,
}

impl NuclWordFinder {
    /// Create a new word finder from query sequences
    pub fn new<S: AsRef<[u8]>>(queries: &[(String, S)], config: NuclWordFinderConfig) -> Self {
        let mut lookup = build_nucl_lookup(queries, config.word_size);

        // Filter out over-represented k-mers if configured
        if let Some(max_hits) = config.max_hits_per_kmer {
            lookup.retain(|_, hits| hits.len() <= max_hits);
        }

        Self { lookup, config }
    }

    /// Create a new word finder from FASTA records
    pub fn from_records(records: &[bio::io::fasta::Record], config: NuclWordFinderConfig) -> Self {
        let mut lookup = build_nucl_lookup_from_records(records, config.word_size);

        if let Some(max_hits) = config.max_hits_per_kmer {
            lookup.retain(|_, hits| hits.len() <= max_hits);
        }

        Self { lookup, config }
    }

    /// Find all seed hits for a subject sequence
    pub fn find_seeds(&self, subject: &[u8]) -> Vec<SeedHit> {
        let mut hits = Vec::new();
        let word_size = self.config.word_size;

        if subject.len() < word_size {
            return hits;
        }

        for s_pos in 0..=(subject.len() - word_size) {
            if let Some(code) = encode_nucl_kmer(subject, s_pos, word_size) {
                if let Some(query_hits) = self.lookup.get(&code) {
                    for &(query_idx, q_pos) in query_hits {
                        hits.push(SeedHit {
                            query_idx: query_idx as usize,
                            query_pos: q_pos as usize,
                            subject_pos: s_pos,
                            word_size,
                        });
                    }
                }
            }
        }

        hits
    }

    /// Get the word size
    pub fn word_size(&self) -> usize {
        self.config.word_size
    }

    /// Get the number of unique k-mers in the index
    pub fn num_kmers(&self) -> usize {
        self.lookup.len()
    }
}

/// A seed hit representing a k-mer match
#[derive(Debug, Clone, Copy)]
pub struct SeedHit {
    /// Index of the query sequence
    pub query_idx: usize,
    /// Position in the query sequence (0-based)
    pub query_pos: usize,
    /// Position in the subject sequence (0-based)
    pub subject_pos: usize,
    /// Length of the seed match
    pub word_size: usize,
}

impl SeedHit {
    /// Get the diagonal of this hit (query_pos - subject_pos)
    pub fn diagonal(&self) -> isize {
        self.query_pos as isize - self.subject_pos as isize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_kmer() {
        let seq = b"ACGTACGT";
        let code = encode_nucl_kmer(seq, 0, 4).unwrap();
        let decoded = decode_nucl_kmer(code, 4);
        assert_eq!(&decoded, b"ACGT");
    }

    #[test]
    fn test_encode_ambiguous() {
        let seq = b"ACNGT";
        assert!(encode_nucl_kmer(seq, 0, 4).is_none());
    }

    #[test]
    fn test_build_lookup() {
        let queries = vec![
            ("q1".to_string(), b"ACGTACGT".to_vec()),
            ("q2".to_string(), b"TACGTACG".to_vec()),
        ];

        let lookup = build_nucl_lookup(&queries, 4);

        // ACGT should be found in both sequences
        let acgt_code = encode_nucl_kmer(b"ACGT", 0, 4).unwrap();
        assert!(lookup.contains_key(&acgt_code));
    }

    #[test]
    fn test_word_finder() {
        let queries = vec![("q1".to_string(), b"ACGTACGTACGT".to_vec())];

        let config = NuclWordFinderConfig {
            word_size: 4,
            max_hits_per_kmer: None,
        };

        let finder = NuclWordFinder::new(&queries, config);
        let subject = b"NNNNACGTNNNN";
        let hits = finder.find_seeds(subject);

        // Should find ACGT at position 4 in subject
        assert!(!hits.is_empty());
        assert!(hits.iter().any(|h| h.subject_pos == 4));
    }
}
