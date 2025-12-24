/// Compatibility mode for BLAST output
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CompatMode {
    /// BLEMIR default behavior (current implementation)
    #[default]
    Blemir,
    /// NCBI BLAST compatible behavior
    Ncbi,
}

impl std::str::FromStr for CompatMode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "blemir" | "default" => Ok(CompatMode::Blemir),
            "ncbi" | "blast" => Ok(CompatMode::Ncbi),
            _ => Err(format!(
                "Unknown compatibility mode: {}. Use 'blemir' or 'ncbi'",
                s
            )),
        }
    }
}

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

/// Supported scoring matrices
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
}

impl std::str::FromStr for ScoringMatrix {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "BLOSUM45" => Ok(ScoringMatrix::Blosum45),
            "BLOSUM50" => Ok(ScoringMatrix::Blosum50),
            "BLOSUM62" => Ok(ScoringMatrix::Blosum62),
            "BLOSUM80" => Ok(ScoringMatrix::Blosum80),
            "BLOSUM90" => Ok(ScoringMatrix::Blosum90),
            "PAM30" => Ok(ScoringMatrix::Pam30),
            "PAM70" => Ok(ScoringMatrix::Pam70),
            "PAM250" => Ok(ScoringMatrix::Pam250),
            _ => Err(format!("Unknown scoring matrix: {}", s)),
        }
    }
}

/// HSP post-processing configuration
#[derive(Debug, Clone, Copy)]
pub struct PostProcessConfig {
    /// Whether to chain nearby HSPs into longer alignments
    pub enable_chaining: bool,
    /// Overlap threshold for filtering redundant HSPs (0.0-1.0)
    pub overlap_threshold: f64,
    /// Maximum gap allowed between HSPs for chaining
    pub max_gap: usize,
    /// Maximum diagonal drift allowed for chaining
    pub max_diag_drift: isize,
}

impl Default for PostProcessConfig {
    fn default() -> Self {
        Self {
            enable_chaining: true,
            overlap_threshold: 0.5,
            max_gap: 1000,
            max_diag_drift: 100,
        }
    }
}

impl PostProcessConfig {
    /// Configuration for NCBI BLAST compatible behavior (no chaining)
    pub fn ncbi_compat() -> Self {
        Self {
            enable_chaining: false,
            overlap_threshold: 0.0,
            max_gap: 0,
            max_diag_drift: 0,
        }
    }
}

/// Seeding configuration
#[derive(Debug, Clone, Copy)]
pub struct SeedConfig {
    /// Word size for seeding
    pub word_size: usize,
    /// Threshold for neighborhood word finding (protein only)
    pub threshold: i32,
    /// Whether to use neighborhood word finding (NCBI style)
    pub use_neighborhood: bool,
}

impl Default for SeedConfig {
    fn default() -> Self {
        Self {
            word_size: 3,
            threshold: 13,
            use_neighborhood: false,
        }
    }
}

impl SeedConfig {
    /// Configuration for NCBI BLAST compatible seeding
    pub fn ncbi_compat(word_size: usize, threshold: i32) -> Self {
        Self {
            word_size,
            threshold,
            use_neighborhood: true,
        }
    }
}

/// Extension configuration
#[derive(Debug, Clone, Copy)]
pub struct ExtensionConfig {
    /// X-drop for ungapped extension
    pub x_drop_ungapped: i32,
    /// X-drop for gapped extension
    pub x_drop_gapped: i32,
    /// Band width for banded alignment
    pub band_width: usize,
    /// Maximum extension length
    pub max_extension: usize,
    /// Whether to compute full traceback for exact statistics
    pub full_traceback: bool,
}

impl Default for ExtensionConfig {
    fn default() -> Self {
        Self {
            x_drop_ungapped: 20,
            x_drop_gapped: 30,
            band_width: 32,
            max_extension: 10000,
            full_traceback: false,
        }
    }
}

impl ExtensionConfig {
    /// Configuration for NCBI BLAST compatible extension
    pub fn ncbi_compat() -> Self {
        Self {
            x_drop_ungapped: 20,
            x_drop_gapped: 30,
            band_width: 32,
            max_extension: 10000,
            full_traceback: true,
        }
    }
}

/// Statistics configuration
#[derive(Debug, Clone, Copy)]
pub struct StatsConfig {
    /// Whether to use length adjustment for effective search space
    pub use_length_adjustment: bool,
}

impl Default for StatsConfig {
    fn default() -> Self {
        Self {
            use_length_adjustment: false,
        }
    }
}

impl StatsConfig {
    /// Configuration for NCBI BLAST compatible statistics
    pub fn ncbi_compat() -> Self {
        Self {
            use_length_adjustment: true,
        }
    }
}

/// Complete BLASTN configuration
#[derive(Debug, Clone)]
pub struct BlastnConfig {
    pub compat_mode: CompatMode,
    pub scoring: NuclScoringSpec,
    pub seed: SeedConfig,
    pub extension: ExtensionConfig,
    pub post_process: PostProcessConfig,
    pub stats: StatsConfig,
}

impl Default for BlastnConfig {
    fn default() -> Self {
        Self {
            compat_mode: CompatMode::Blemir,
            scoring: NuclScoringSpec::default(),
            seed: SeedConfig::default(),
            extension: ExtensionConfig::default(),
            post_process: PostProcessConfig::default(),
            stats: StatsConfig::default(),
        }
    }
}

impl BlastnConfig {
    pub fn ncbi_compat(scoring: NuclScoringSpec, word_size: usize) -> Self {
        Self {
            compat_mode: CompatMode::Ncbi,
            scoring,
            seed: SeedConfig {
                word_size,
                threshold: 0,
                use_neighborhood: false,
            },
            extension: ExtensionConfig::ncbi_compat(),
            post_process: PostProcessConfig::ncbi_compat(),
            stats: StatsConfig::ncbi_compat(),
        }
    }
}

/// Complete TBLASTX configuration
#[derive(Debug, Clone)]
pub struct TblastxConfig {
    pub compat_mode: CompatMode,
    pub scoring: ProteinScoringSpec,
    pub seed: SeedConfig,
    pub extension: ExtensionConfig,
    pub post_process: PostProcessConfig,
    pub stats: StatsConfig,
}

impl Default for TblastxConfig {
    fn default() -> Self {
        Self {
            compat_mode: CompatMode::Blemir,
            scoring: ProteinScoringSpec::default(),
            seed: SeedConfig::default(),
            extension: ExtensionConfig {
                x_drop_ungapped: 20,
                x_drop_gapped: 25,
                band_width: 32,
                max_extension: 500,
                full_traceback: false,
            },
            post_process: PostProcessConfig::default(),
            stats: StatsConfig::default(),
        }
    }
}

impl TblastxConfig {
    pub fn ncbi_compat(scoring: ProteinScoringSpec, word_size: usize, threshold: i32) -> Self {
        Self {
            compat_mode: CompatMode::Ncbi,
            scoring,
            seed: SeedConfig::ncbi_compat(word_size, threshold),
            extension: ExtensionConfig {
                x_drop_ungapped: 20,
                x_drop_gapped: 25,
                band_width: 32,
                max_extension: 500,
                full_traceback: true,
            },
            post_process: PostProcessConfig::ncbi_compat(),
            stats: StatsConfig::ncbi_compat(),
        }
    }
}
