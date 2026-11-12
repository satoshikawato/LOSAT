/// Scoring specification for nucleotide alignments
#[derive(Debug, Clone, Copy)]
pub struct NuclScoringSpec {
    pub reward: i32,
    pub penalty: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for NuclScoringSpec {
    fn default() -> Self {
        Self {
            reward: 1,
            penalty: -2,
            gap_open: 0,
            gap_extend: 0,
        }
    }
}

/// Scoring specification for protein alignments
#[derive(Debug, Clone, Copy)]
pub struct ProteinScoringSpec {
    pub matrix: ScoringMatrix,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for ProteinScoringSpec {
    fn default() -> Self {
        Self {
            matrix: ScoringMatrix::Blosum62,
            gap_open: 11,
            gap_extend: 1,
        }
    }
}

/// The nine built-in protein matrices of NCBI BLAST+ (`util/tables/sm_*.c`, each with a
/// Karlin-Altschul table in `blast_stat.c`).
///
/// NCBI reference: c++/src/algo/blast/core/blast_stat.c:2958-2992
/// ```c
///    matrix_info = MatrixInfoNew("BLOSUM80", blosum80_values, blosum80_prefs, BLOSUM80_VALUES_MAX);
/// ...
///    if (!standard_only) {
///        matrix_info = MatrixInfoNew("IDENTITY", prot_idenity_values,
///                                    prot_identity_prefs,
///                                    PROT_IDENTITY_VALUES_MAX);
///        ListNodeAddPointer(&retval, 0, matrix_info);
///    }
/// ```
/// (BLOSUM62_20 is compiled only under `BLOSUM62_20_ENABLE`, which the 2.17.0 build does
/// not define.)
/// `Identity` is the last variant so that the others keep their discriminants.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ScoringMatrix {
    Blosum45,
    Blosum50,
    #[default]
    Blosum62,
    Blosum80,
    Blosum90,
    Pam30,
    Pam70,
    Pam250,
    Identity,
}

impl ScoringMatrix {
    /// Every matrix, in the order of the variants.
    pub const ALL: [ScoringMatrix; 9] = [
        ScoringMatrix::Blosum45,
        ScoringMatrix::Blosum50,
        ScoringMatrix::Blosum62,
        ScoringMatrix::Blosum80,
        ScoringMatrix::Blosum90,
        ScoringMatrix::Pam30,
        ScoringMatrix::Pam70,
        ScoringMatrix::Pam250,
        ScoringMatrix::Identity,
    ];

    /// The name of the matrix's table (`MatrixInfo.name`).
    pub const fn ncbi_name(self) -> &'static str {
        match self {
            ScoringMatrix::Blosum45 => "BLOSUM45",
            ScoringMatrix::Blosum50 => "BLOSUM50",
            ScoringMatrix::Blosum62 => "BLOSUM62",
            ScoringMatrix::Blosum80 => "BLOSUM80",
            ScoringMatrix::Blosum90 => "BLOSUM90",
            ScoringMatrix::Pam30 => "PAM30",
            ScoringMatrix::Pam70 => "PAM70",
            ScoringMatrix::Pam250 => "PAM250",
            ScoringMatrix::Identity => "IDENTITY",
        }
    }
}

impl std::str::FromStr for ScoringMatrix {
    type Err = String;

    /// The matrix whose table name equals `s` ignoring the case, as NCBI finds a matrix's
    /// table.
    ///
    /// NCBI reference: c++/src/algo/blast/core/blast_stat.c:3593-3605
    /// ```c
    ///    vnp = head = BlastLoadMatrixValues(standard_only);
    ///    while (vnp)
    ///    {
    ///       matrix_info = vnp->ptr;
    ///       if (strcasecmp(matrix_info->name, matrix_name) == 0)
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        ScoringMatrix::ALL
            .into_iter()
            .find(|matrix| matrix.ncbi_name().eq_ignore_ascii_case(s))
            .ok_or_else(|| format!("Unknown scoring matrix: {}", s))
    }
}
