//! Genetic Code Singleton
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/gencode_singleton.c
//!
//! This module provides genetic code translation tables for converting
//! nucleotide codons to amino acids.

/// Genetic code translation table
pub struct GeneticCode {
    pub table: [u8; 64],
}

impl GeneticCode {
    /// Translate a codon to an amino acid
    pub fn get(&self, codon: &[u8]) -> u8 {
        if codon.len() != 3 {
            return b'X';
        }
        let mut idx = 0;
        for &b in codon {
            idx <<= 2;
            match b.to_ascii_uppercase() {
                b'T' | b'U' => idx |= 0,
                b'C' => idx |= 1,
                b'A' => idx |= 2,
                b'G' => idx |= 3,
                _ => return b'X',
            }
        }
        self.table[idx]
    }

    /// Create a genetic code from an NCBI genetic code ID
    ///
    /// Supported IDs: 1-6, 9-16, 21-31, 33
    /// Reference: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    pub fn from_id(id: u8) -> Self {
        let table_str = match id {
            1 => b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            2 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSS**VVVVAAAADDEEGGGG",
            3 => b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            4 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            5 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
            6 => b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            9 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
            10 => b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            11 => b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            12 => b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            13 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
            14 => b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
            15 => b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            16 => b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            21 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
            22 => b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            23 => b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            24 => b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
            25 => b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            26 => b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            27 => b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            28 => b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            29 => b"FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            30 => b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            31 => b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            33 => b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
            _ => {
                eprintln!(
                    "Warning: Genetic code {} not implemented or invalid, using Standard (1).",
                    id
                );
                b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
            }
        };

        let mut table = [0u8; 64];
        table.copy_from_slice(table_str);
        GeneticCode { table }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_standard_code() {
        let code = GeneticCode::from_id(1);

        // ATG -> M (Methionine, start codon)
        assert_eq!(code.get(b"ATG"), b'M');

        // TAA, TAG, TGA -> * (Stop codons)
        assert_eq!(code.get(b"TAA"), b'*');
        assert_eq!(code.get(b"TAG"), b'*');
        assert_eq!(code.get(b"TGA"), b'*');

        // TTT, TTC -> F (Phenylalanine)
        assert_eq!(code.get(b"TTT"), b'F');
        assert_eq!(code.get(b"TTC"), b'F');
    }

    #[test]
    fn test_ambiguous_codon() {
        let code = GeneticCode::from_id(1);
        assert_eq!(code.get(b"NNN"), b'X');
        assert_eq!(code.get(b"ATN"), b'X');
    }

    #[test]
    fn test_invalid_length() {
        let code = GeneticCode::from_id(1);
        assert_eq!(code.get(b"AT"), b'X');
        assert_eq!(code.get(b"ATGA"), b'X');
    }
}
