//! The report records of the serialized formats (outfmt 8, 9, 12): per query, borrowed
//! views of a program's final hit lists in NCBI's terms, built only when such a format is
//! requested (`ReportNeeds::seqalign`), into buffers reused from query to query, and dropped
//! after the query's document is written. Nothing here flows through the search's pools,
//! pipes or heaps; each program's adapter (`algorithm/<program>/seqalign.rs`) builds them
//! from its own HSPs.
//!
//! The fields are those of NCBI's `BlastHSP` that the Seq-align conversion reads.
//!
//! NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:126-148
//! ```c
//! typedef struct BlastHSP {
//!    Int4 score;           /**< This HSP's raw score */
//!    Int4 num_ident;       /**< Number of identical base pairs in this HSP */
//!    double bit_score;     /**< Bit score, calculated from score */
//!    double evalue;        /**< This HSP's e-value */
//!    BlastSeg query;       /**< Query sequence info. */
//!    BlastSeg subject;     /**< Subject sequence info. */
//!    Int4     context;     /**< Context number of query */
//!    GapEditScript* gap_info;/**< ALL gapped alignment is here */
//!    Int4 num;             /**< How many HSP's are linked together for sum
//!                               statistics evaluation? If unset (0), this HSP is
//!                               not part of a linked set, i.e. value 0 is treated
//!                               the same way as 1. */
//!    Int2		comp_adjustment_method;  /**< which mode of composition
//!                                               adjustment was used; relevant
//!                                               only for blastp and tblastn */
//!    SPHIHspInfo* pat_info; /**< In PHI BLAST, information about this pattern
//!                                  match. */
//!    Int4 num_positives;
//!
//!    BlastHSPMappingInfo* map_info;
//!
//! } BlastHSP;
//! ```
//! NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:96-101
//! ```c
//! typedef struct BlastSeg {
//!    Int2 frame;  /**< Translation frame */
//!    Int4 offset; /**< Start of hsp */
//!    Int4 end;    /**< End of hsp */
//!    Int4 gapped_start;/**< Where the gapped extension started. */
//! } BlastSeg;
//! ```

use crate::blastinput::app::{Program, ReportFormat};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::common::GapEditOp;

use super::seqid::SeqId;

/// What the requested formats of a search need from its reports beyond outfmt 0, 6 and 7,
/// decided once per search where the program reads its formats; the search itself never
/// reads it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ReportNeeds {
    /// A Seq-annot of Seq-aligns is written (outfmt 8, 12): the reports build
    /// `HspRecord`s, and a database's name and title are checked
    /// (`report::serial::seqalign::check_database`).
    pub seqalign: bool,
}

impl ReportNeeds {
    /// The needs of the requested formats.
    pub fn of(formats: impl IntoIterator<Item = ReportFormat>) -> Self {
        let mut needs = Self::default();
        for format in formats {
            needs.seqalign |= crate::report::serial::SerialFormat::of(format).is_some();
        }
        needs
    }
}

/// One HSP of a query's final hit list, as NCBI's `BlastHSP` holds it when the Seq-aligns
/// are made (`BlastHitList2SeqAlign_OMF`).
#[derive(Clone, Copy, Debug)]
pub struct HspRecord<'a> {
    /// The index of the HSP's subject in the query's `SubjectRecord`s. The HSPs of one
    /// subject are consecutive, in NCBI's Seq-align order.
    pub subject: u32,
    /// `query.offset` and `query.end` (exclusive): positions in the HSP's query context
    /// (the reverse complement for a minus-strand BLASTN HSP; the frame's protein for a
    /// translated query), in the searched letters.
    pub q_offset: i32,
    pub q_end: i32,
    /// `subject.offset` and `subject.end` (exclusive), likewise for the subject.
    pub s_offset: i32,
    pub s_end: i32,
    /// `query.frame` and `subject.frame`: 0 for a protein, +1/-1 for a nucleotide strand,
    /// +-1..3 for a translation frame.
    pub q_frame: i8,
    pub s_frame: i8,
    /// `score` (the raw score).
    pub score: i32,
    /// `num`: the size of the HSP's linked set (0 or 1 when it is not linked).
    pub num: i32,
    /// `evalue` (the sum e-value of a linked set) and `bit_score`.
    pub evalue: f64,
    pub bit_score: f64,
    /// `num_ident` and `num_positives`: the engine's counts (not a display recount).
    pub num_ident: i32,
    pub num_positives: i32,
    /// `comp_adjustment_method`: 0, 1 or 2.
    pub comp_adjust: u8,
    /// `gap_info`: the edit script, in the context's direction. An empty script stands
    /// for one substitution over the query's extent.
    pub ops: &'a [GapEditOp],
}

/// A query of the search, as the Seq-aligns name and measure it.
#[derive(Clone, Copy, Debug)]
pub struct QueryRecord<'a> {
    /// `local str "Query_N"`, N = the record's position in the input.
    pub id: SeqId<'a>,
    /// The length of the searched letters (with `-query_loc`, of the interval;
    /// `query_data.GetSeqLength`).
    pub searched_length: i32,
    /// The start of the `-query_loc` interval in the record (0 without one), which
    /// `RemapToQueryLoc` adds to the query row.
    pub offset: i32,
}

/// A subject with HSPs, as the Seq-aligns name and measure it.
#[derive(Clone, Copy, Debug)]
pub struct SubjectRecord<'a> {
    /// `local str "Subject_N"` (`-subject`) or `general { db "BL_ORD_ID", tag id <oid> }`
    /// (`-db`).
    pub id: SeqId<'a>,
    /// The length of the searched letters (with `-subject_loc`, of the interval).
    pub searched_length: i32,
    /// The start of the `-subject_loc` interval in the record (0 without one, and always 0
    /// for `-db`), which `s_RemapToSubjectLoc` adds to the subject row.
    pub offset: i32,
}

/// What a whole search contributes to every document.
#[derive(Clone, Copy, Debug)]
pub struct RunRecord<'a> {
    /// The program, which names the "Blast Type" and selects the segment form.
    pub program: Program,
    /// The database title (`-db`), or `None` for a `-subject` search.
    pub db_title: Option<&'a [u8]>,
    /// Whether the search is gapped (`BLASTHspListToSeqAlign`) or ungapped
    /// (`BLASTUngappedHspListToSeqAlign`).
    pub gapped: bool,
}

impl RunRecord<'_> {
    /// Whether the Seq-aligns are remapped to their subjects' locations: a `-subject` search,
    /// whose sequence source can return partial sequences (a database's cannot).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1630-1645
    /// ```c++
    /// 	if (seqinfo_src->CanReturnPartialSequence() == true)
    /// 	{
    ///         	CConstRef<CSeq_loc> subj_loc = seqinfo_src->GetSeqLoc(kOid);
    ///         	NON_CONST_ITERATE(vector<CRef<CSeq_align > >, iter, hit_align) {
    ///              	   RemapToQueryLoc(*iter, query_loc);
    ///              	   if ( !is_ooframe )
    ///                    	s_RemapToSubjectLoc(*iter, *subj_loc);
    ///              	   seq_aligns->Set().push_back(*iter);
    /// 		}
    ///         }
    /// 	else
    /// 	{
    /// 		NON_CONST_ITERATE(vector<CRef<CSeq_align > >, iter, hit_align) {
    ///            	    RemapToQueryLoc(*iter, query_loc);
    ///            	    seq_aligns->Set().push_back(*iter);
    /// 		}
    /// 	}
    /// ```
    pub fn remaps_subjects(&self) -> bool {
        self.db_title.is_none()
    }
}

/// How a program's reports name its subjects: the records of a `-subject` set
/// (`Subject_N`), or the ordinal IDs of a database's reported sequences by the rank that
/// the final hit lists use, with the database's title (NCBI's `m_DbInfo.front().definition`).
#[derive(Clone, Copy)]
pub enum SubjectNames<'a> {
    Records(&'a [FastaRecord]),
    Database {
        ordinals: &'a [u32],
        title: &'a [u8],
    },
}

impl<'a> SubjectNames<'a> {
    /// The database title of `RunRecord::db_title`.
    pub fn db_title(&self) -> Option<&'a [u8]> {
        match *self {
            SubjectNames::Records(_) => None,
            SubjectNames::Database { title, .. } => Some(title),
        }
    }

    /// The record of subject `index` (by rank for a database): its Seq-id, the searched
    /// length and the `-subject_loc` start `offset`. A database's sequences are searched
    /// whole and have no subject remap (`CanReturnPartialSequence` is false for CSeqDB,
    /// blast_seqalign.cpp:1630-1645); their length is `database_length` (the caller's, where
    /// a row reads it).
    pub fn record(&self, index: usize, offset: usize, database_length: i32) -> SubjectRecord<'a> {
        match *self {
            SubjectNames::Records(records) => SubjectRecord {
                id: SeqId::Local(&records[index].local_id),
                searched_length: records[index].seq().len() as i32,
                offset: offset as i32,
            },
            SubjectNames::Database { ordinals, .. } => SubjectRecord {
                id: SeqId::BlOrdId(ordinals[index]),
                searched_length: database_length,
                offset: 0,
            },
        }
    }
}

/// The queries of a program's reports: the searched records (with their local IDs
/// `Query_N`), their searched lengths, and where each lies in its record (`-query_loc`).
pub struct QueryNames<'a> {
    pub records: &'a [FastaRecord],
    pub searched_lengths: &'a [usize],
    pub placements: &'a Placements,
}

impl<'a> QueryNames<'a> {
    /// The record of query `index`.
    pub fn record(&self, index: usize) -> QueryRecord<'a> {
        QueryRecord {
            id: SeqId::Local(&self.records[index].local_id),
            searched_length: self.searched_lengths[index] as i32,
            offset: self.placements.offset(index) as i32,
        }
    }
}
