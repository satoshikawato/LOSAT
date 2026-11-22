//! The Blast4 parameters: the one registry of their names and value types (`CBlast4Field`), which
//! the request writer (`options`, `writer`: what NCBI's `Set*` calls record) and the strategy
//! reader (`read::options`: what NCBI's `CBlastOptionsBuilder` applies) share.
//!
//! NCBI reference (598d8ae6): c++/src/objects/blast/names.cpp:367-373
//! A field is a name with the `Blast4-value` alternative that NCBI's table gives it
//! (names.cpp:53-337, `CBlast4Field::Get`); a parameter is the field's only when both agree
//! ([`Field::matches`]). One name has two fields: `Culling` is the culling limit as an integer
//! ([`CULLING_LIMIT`]) and the builder's culling switch as a boolean ([`CULLING`]).
//!
//! A [`Field`] is made here only, so every parameter that LOSAT writes or reads names a field of
//! this table.

/// An alternative of `Blast4-value` (blast.asn:512-551), in the module's order.
///
/// NCBI reference (598d8ae6): c++/src/objects/blast/blast.asn:512-526
/// The discriminants are NCBI's `CBlast4_value::E_Choice` (`e_not_set` is 0, the first
/// alternative 1), the number that `Which()` returns and NCBI's messages print.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ValueType {
    BigInteger = 1,
    Bioseq,
    Boolean,
    Cutoff,
    Integer,
    Matrix,
    Real,
    SeqAlign,
    SeqId,
    SeqLoc,
    StrandType,
    String,
    BigIntegerList,
    BioseqList,
    BooleanList,
    CutoffList,
    IntegerList,
    MatrixList,
    RealList,
    SeqAlignList,
    SeqIdList,
    SeqLocList,
    StrandTypeList,
    StringList,
    BioseqSet,
    SeqAlignSet,
    QueryMask,
}

impl ValueType {
    /// Every alternative, in the module's order.
    pub const ALL: [ValueType; 27] = [
        ValueType::BigInteger,
        ValueType::Bioseq,
        ValueType::Boolean,
        ValueType::Cutoff,
        ValueType::Integer,
        ValueType::Matrix,
        ValueType::Real,
        ValueType::SeqAlign,
        ValueType::SeqId,
        ValueType::SeqLoc,
        ValueType::StrandType,
        ValueType::String,
        ValueType::BigIntegerList,
        ValueType::BioseqList,
        ValueType::BooleanList,
        ValueType::CutoffList,
        ValueType::IntegerList,
        ValueType::MatrixList,
        ValueType::RealList,
        ValueType::SeqAlignList,
        ValueType::SeqIdList,
        ValueType::SeqLocList,
        ValueType::StrandTypeList,
        ValueType::StringList,
        ValueType::BioseqSet,
        ValueType::SeqAlignSet,
        ValueType::QueryMask,
    ];

    /// `CBlast4_value::Which()`: the alternative's `E_Choice` number.
    pub fn which(self) -> i32 {
        self as i32
    }

    /// The alternative's tag (its position from 0), which the encoders take and the generic
    /// reader gives (`Value::Choice`).
    pub fn tag(self) -> u8 {
        self as u8 - 1
    }

    /// The alternative of tag `tag`.
    pub fn of_tag(tag: u16) -> Option<ValueType> {
        Self::ALL.get(usize::from(tag)).copied()
    }

    /// The alternative's name in the module.
    pub fn variant(self) -> &'static str {
        match self {
            ValueType::BigInteger => "big-integer",
            ValueType::Bioseq => "bioseq",
            ValueType::Boolean => "boolean",
            ValueType::Cutoff => "cutoff",
            ValueType::Integer => "integer",
            ValueType::Matrix => "matrix",
            ValueType::Real => "real",
            ValueType::SeqAlign => "seq-align",
            ValueType::SeqId => "seq-id",
            ValueType::SeqLoc => "seq-loc",
            ValueType::StrandType => "strand-type",
            ValueType::String => "string",
            ValueType::BigIntegerList => "big-integer-list",
            ValueType::BioseqList => "bioseq-list",
            ValueType::BooleanList => "boolean-list",
            ValueType::CutoffList => "cutoff-list",
            ValueType::IntegerList => "integer-list",
            ValueType::MatrixList => "matrix-list",
            ValueType::RealList => "real-list",
            ValueType::SeqAlignList => "seq-align-list",
            ValueType::SeqIdList => "seq-id-list",
            ValueType::SeqLocList => "seq-loc-list",
            ValueType::StrandTypeList => "strand-type-list",
            ValueType::StringList => "string-list",
            ValueType::BioseqSet => "bioseq-set",
            ValueType::SeqAlignSet => "seq-align-set",
            ValueType::QueryMask => "query-mask",
        }
    }
}

/// A Blast4 parameter's field (`CBlast4Field`): its name and the value type of NCBI's table.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Field {
    name: &'static str,
    value_type: ValueType,
}

impl Field {
    const fn new(name: &'static str, value_type: ValueType) -> Field {
        Field { name, value_type }
    }

    /// `CBlast4Field::GetName`.
    pub fn name(self) -> &'static str {
        self.name
    }

    /// `CBlast4Field::GetType`.
    pub fn value_type(self) -> ValueType {
        self.value_type
    }

    /// `CBlast4Field::Match` (names.cpp:367-373): a parameter named `name` whose value is the
    /// alternative `value_type` is this field's.
    pub fn matches(self, name: &[u8], value_type: ValueType) -> bool {
        value_type == self.value_type && name == self.name.as_bytes()
    }
}

// NCBI reference (598d8ae6): c++/src/objects/blast/names.cpp:53-60
// The constants below are the cases of `CBlast4Field::Get` in names.cpp's order, each named after
// its `EBlastOptIdx` (`eBlastOpt_WordThreshold` is `WORD_THRESHOLD`); the commented-out cases
// (`Program`, `LookupTableType`, `XDropoff`, ...) have no field. The cases that neither LOSAT's
// writer nor NCBI's `CBlastOptionsBuilder` names are left out (`ComplexityAdjustMode`, `MaskLevel`,
// `UseRealDbSize` and the `Web_` fields but `Web_RunPsiBlast` and `Web_StepNumber`): the builder
// takes a parameter of those names as one it does not know, which needs no field.
use ValueType as T;

pub const WORD_THRESHOLD: Field = Field::new("WordThreshold", T::Integer);
pub const WORD_SIZE: Field = Field::new("WordSize", T::Integer);
pub const MB_TEMPLATE_LENGTH: Field = Field::new("MBTemplateLength", T::Integer);
pub const MB_TEMPLATE_TYPE: Field = Field::new("MBTemplateType", T::Integer);
pub const FILTER_STRING: Field = Field::new("FilterString", T::String);
pub const MASK_AT_HASH: Field = Field::new("MaskAtHash", T::Boolean);
pub const DUST_FILTERING: Field = Field::new("DustFiltering", T::Boolean);
pub const DUST_FILTERING_LEVEL: Field = Field::new("DustFilteringLevel", T::Integer);
pub const DUST_FILTERING_WINDOW: Field = Field::new("DustFilteringWindow", T::Integer);
pub const DUST_FILTERING_LINKER: Field = Field::new("DustFilteringLinker", T::Integer);
pub const SEG_FILTERING: Field = Field::new("SegFiltering", T::Boolean);
pub const SEG_FILTERING_WINDOW: Field = Field::new("SegFilteringWindow", T::Integer);
pub const SEG_FILTERING_LOCUT: Field = Field::new("SegFilteringLocut", T::Real);
pub const SEG_FILTERING_HICUT: Field = Field::new("SegFilteringHicut", T::Real);
pub const REPEAT_FILTERING: Field = Field::new("RepeatFiltering", T::Boolean);
pub const REPEAT_FILTERING_DB: Field = Field::new("RepeatFilteringDB", T::String);
pub const STRAND_OPTION: Field = Field::new("StrandOption", T::StrandType);
pub const QUERY_GENETIC_CODE: Field = Field::new("QueryGeneticCode", T::Integer);
pub const WINDOW_SIZE: Field = Field::new("WindowSize", T::Integer);
pub const GAP_X_DROPOFF: Field = Field::new("GapXDropoff", T::Real);
pub const GAP_X_DROPOFF_FINAL: Field = Field::new("GapXDropoffFinal", T::Real);
pub const GAP_TRIGGER: Field = Field::new("GapTrigger", T::Real);
pub const GAP_EXTN_ALGORITHM: Field = Field::new("GapExtnAlgorithm", T::Integer);
pub const HITLIST_SIZE: Field = Field::new("HitlistSize", T::Integer);
/// `eBlastOpt_CullingLimit`: named `Culling`, as [`CULLING`] is.
pub const CULLING_LIMIT: Field = Field::new("Culling", T::Integer);
pub const EVALUE_THRESHOLD: Field = Field::new("EvalueThreshold", T::Cutoff);
pub const CUTOFF_SCORE: Field = Field::new("CutoffScore", T::Cutoff);
pub const PERCENT_IDENTITY: Field = Field::new("PercentIdentity", T::Real);
pub const QUERY_COV_HSP_PERC: Field = Field::new("QueryCovHspPerc", T::Real);
pub const MAX_HSPS_PER_SUBJECT: Field = Field::new("MaxHspsPerSubject", T::Integer);
pub const SUM_STATISTICS_MODE: Field = Field::new("SumStatistics", T::Boolean);
pub const LONGEST_INTRON_LENGTH: Field = Field::new("LongestIntronLength", T::Integer);
/// `eBlastOpt_GappedMode`: "same as !ungapped mode" (names.cpp:153).
pub const GAPPED_MODE: Field = Field::new("UngappedMode", T::Boolean);
pub const MATRIX_NAME: Field = Field::new("MatrixName", T::String);
pub const MATCH_REWARD: Field = Field::new("MatchReward", T::Integer);
pub const MISMATCH_PENALTY: Field = Field::new("MismatchPenalty", T::Integer);
pub const GAP_OPENING_COST: Field = Field::new("GapOpeningCost", T::Integer);
pub const GAP_EXTENSION_COST: Field = Field::new("GapExtensionCost", T::Integer);
pub const OUT_OF_FRAME_MODE: Field = Field::new("OutOfFrameMode", T::Boolean);
pub const DB_LENGTH: Field = Field::new("DbLength", T::BigInteger);
pub const EFFECTIVE_SEARCH_SPACE: Field = Field::new("EffectiveSearchSpace", T::BigInteger);
pub const DB_GENETIC_CODE: Field = Field::new("DbGeneticCode", T::Integer);
pub const PHI_PATTERN: Field = Field::new("PHIPattern", T::String);
pub const INCLUSION_THRESHOLD: Field = Field::new("InclusionThreshold", T::Real);
/// `eBlastOpt_PseudoCount`.
pub const PSEUDO_COUNT: Field = Field::new("PseudoCountWeight", T::Integer);
pub const GAP_TRACEBACK_ALGORITHM: Field = Field::new("GapTracebackAlgorithm", T::Integer);
pub const COMPOSITION_BASED_STATS: Field = Field::new("CompositionBasedStats", T::Integer);
pub const SMITH_WATERMAN_MODE: Field = Field::new("SmithWatermanMode", T::Boolean);
pub const UNIFIED_P: Field = Field::new("UnifiedP", T::Integer);
pub const WINDOW_MASKER_DATABASE: Field = Field::new("WindowMaskerDatabase", T::String);
pub const WINDOW_MASKER_TAX_ID: Field = Field::new("WindowMaskerTaxId", T::Integer);
pub const FORCE_MB_INDEX: Field = Field::new("ForceMbIndex", T::Boolean);
pub const MB_INDEX_NAME: Field = Field::new("MbIndexName", T::String);
pub const BEST_HIT_SCORE_EDGE: Field = Field::new("BestHitScoreEdge", T::Real);
pub const BEST_HIT_OVERHANG: Field = Field::new("BestHitOverhang", T::Real);
/// `eBlastOpt_SubjectBestHit`.
pub const SUBJECT_BEST_HIT: Field = Field::new("UseSubjectBestHit", T::Boolean);
pub const IGNORE_MSA_MASTER: Field = Field::new("IgnoreMsaMaster", T::Boolean);
pub const DOMAIN_INCLUSION_THRESHOLD: Field = Field::new("DomainInclusionThreshold", T::Real);
/// `eBlastOpt_Culling`: named `Culling`, as [`CULLING_LIMIT`] is.
pub const CULLING: Field = Field::new("Culling", T::Boolean);
pub const ENTREZ_QUERY: Field = Field::new("EntrezQuery", T::String);
pub const FINAL_DB_SEQ: Field = Field::new("FinalDbSeq", T::Integer);
pub const FIRST_DB_SEQ: Field = Field::new("FirstDbSeq", T::Integer);
pub const GI_LIST: Field = Field::new("GiList", T::IntegerList);
pub const DB_FILTERING_ALGORITHM_ID: Field = Field::new("DbFilteringAlgorithmId", T::Integer);
pub const DB_FILTERING_ALGORITHM_KEY: Field = Field::new("DbFilteringAlgorithmKey", T::String);
pub const SUBJECT_MASKING_TYPE: Field = Field::new("SubjectMaskingType", T::Integer);
pub const HSP_RANGE_MAX: Field = Field::new("HspRangeMax", T::Integer);
pub const L_CASE_MASK: Field = Field::new("LCaseMask", T::QueryMask);
pub const MATRIX_TABLE: Field = Field::new("MatrixTable", T::Matrix);
pub const NEGATIVE_GI_LIST: Field = Field::new("NegativeGiList", T::IntegerList);
pub const REQUIRED_END: Field = Field::new("RequiredEnd", T::Integer);
pub const REQUIRED_START: Field = Field::new("RequiredStart", T::Integer);
pub const TAXID_LIST: Field = Field::new("TaxidList", T::IntegerList);
pub const NEGATIVE_TAXID_LIST: Field = Field::new("NegativeTaxidList", T::IntegerList);
pub const WEB_RUN_PSI_BLAST: Field = Field::new("Web_RunPsiBlast", T::Boolean);
pub const WEB_STEP_NUMBER: Field = Field::new("Web_StepNumber", T::Integer);

#[cfg(test)]
#[path = "fields_tests.rs"]
mod tests;
