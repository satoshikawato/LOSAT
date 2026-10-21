//! The report records of outfmt 5: per query, borrowed views of a program's final hit lists in
//! the terms of NCBI's `CBlastOutput` model (objects/blastxml/blastxml.asn), built at write time
//! into buffers reused from query to query. Nothing here flows through the search's pools,
//! pipes or heaps; each program's adapter (`algorithm/<program>/blast_xml.rs`) builds them.

use crate::common::GapEditOp;

use super::matrix::FormattingMatrix;

/// What a whole search contributes to every document (`BlastOutput` and `Parameters`).
pub struct XmlRun<'a> {
    /// `BlastOutput_program`: `Blast_ProgramNameFromType` (`blastn`, `blastp`, ...).
    pub program: &'static str,
    /// `BlastOutput_version`: `BLASTP 2.17.0+`.
    pub version: &'a str,
    /// `BlastOutput_reference`, before the XML escapes (`Sch&auml;ffer` stands as written).
    pub reference: &'static str,
    /// `BlastOutput_db`: the `-db` value; empty for `-subject`.
    pub db_name: &'a [u8],
    /// `Parameters_matrix`; none for BLASTN.
    pub matrix: Option<&'a str>,
    /// `Parameters_expect`.
    pub expect: f64,
    /// `Parameters_gap-open` and `Parameters_gap-extend`.
    pub gap_open: i32,
    pub gap_extend: i32,
    /// `Parameters_sc-match` and `Parameters_sc-mismatch` (BLASTN's reward and penalty).
    pub reward_penalty: Option<(i32, i32)>,
    /// `Parameters_filter` (`filter::filter_string`).
    pub filter: &'a [u8],
    /// `Statistics_db-num` and `Statistics_db-len`: the database's sequences and letters; 0 and
    /// 0 for `-subject`.
    pub db_num: i32,
    pub db_len: i64,
    /// The positives table of the recount; `None` for BLASTN (a nucleotide search shows `|`).
    pub positives: Option<&'a FormattingMatrix>,
    /// The letter of a masked query position: `N` for BLASTN, `X` otherwise
    /// (blastxml_format.cpp:280-283).
    pub mask_char: u8,
}

/// The `Statistics` of a query.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct XmlStats {
    pub hsp_len: i32,
    pub eff_space: i64,
    pub kappa: f64,
    pub lambda: f64,
    pub entropy: f64,
}

impl XmlStats {
    /// A query without a valid context (an all-X protein, an all-N nucleotide, an empty
    /// record) or in a batch that NCBI did not search: `CBlastAncillaryData` gives -1 blocks
    /// (api/local_blast.cpp:177-217), the data class 0 for the lengths.
    pub const NONE: XmlStats = XmlStats {
        hsp_len: 0,
        eff_space: 0,
        kappa: -1.0,
        lambda: -1.0,
        entropy: -1.0,
    };
}

/// A query of the search (an `Iteration`).
pub struct XmlQuery<'a> {
    /// `Iteration_query-ID`: `Query_N` (`CBlastFormatUtil::GetSeqIdString`).
    pub id: &'a [u8],
    /// `GenerateDefline` of the query, as generated: empty for an empty nucleotide title.
    /// `Iteration_query-def` is `No definition line` then; `Search_query-title` (XML2) is
    /// absent (blastxml2_format.cpp:590-592).
    pub title: &'a [u8],
    /// `Iteration_query-len`: the length of the searched interval.
    pub len: i32,
    pub stats: XmlStats,
    /// `Iteration_message` (`names::xml_message`); none when empty.
    pub message: Option<&'a [u8]>,
    /// Whether NCBI searched the query's batch (it has a valid context in it). The XML2
    /// `Search` of a query that it did not search has no `hits` member and adds
    /// ` \nNo hits found` to its message (blastxml2_format.cpp:614-622); outfmt 5 does not
    /// read it.
    pub searched: bool,
    /// The query's masked regions in NCBI's order and frames (`TMaskedQueryRegions`):
    /// outfmt 5 masks an HSP's row with the regions of the HSP's own frame, XML2 lists every
    /// region of a plus frame in `Search_query-masking` and masks with the regions of frame 0
    /// or the HSP's translated frame.
    pub masks: &'a [QueryMask],
}

/// A subject with HSPs (a `Hit`).
pub struct XmlHit<'a> {
    pub id: &'a [u8],
    /// The subject's title as generated or stored: empty for an empty nucleotide `-subject`
    /// title or an empty stored title. `Hit_def` is `No definition line` then; XML2's
    /// `HitDescr_title` is empty.
    pub title: &'a [u8],
    pub accession: &'a [u8],
    /// `Hit_len`: the length of the subject's record.
    pub len: i32,
    /// The Blast-def-line's taxid of a `-db` sequence, 0 for none (XML2's `HitDescr_taxid`).
    pub taxid: i32,
}

/// A masked region of the query in NCBI's form (`CSeqLocInfo`): inclusive 0-based positions in
/// the query record and the frame of the context it belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct QueryMask {
    pub from: i32,
    pub to: i32,
    pub frame: i32,
}

/// An HSP (an `Hsp`) with its alignment rows.
pub struct XmlHsp<'a> {
    /// `Hsp_score` (the Seq-align's raw score).
    pub score: i32,
    /// The Seq-align's bit score and e-value (the e-value is 0 below 1e-180; a linked set's
    /// sum e-value); the writer rounds both to six digits.
    pub bit_score: f64,
    pub evalue: f64,
    pub query_from: i32,
    pub query_to: i32,
    pub hit_from: i32,
    pub hit_to: i32,
    pub query_frame: i32,
    pub hit_frame: i32,
    /// The rows as the alignment shows them: the unmasked, upper-case query letters and the
    /// subject letters, `-` for a gap, equally long.
    pub qrow: &'a [u8],
    pub srow: &'a [u8],
    /// The edit script of the rows (`None`: one segment of substitutions over the whole row), which
    /// gives the Dense-seg segments that the mask loop walks.
    pub ops: Option<&'a [GapEditOp]>,
    /// The 0-based record position of the query row's first residue and how the residues
    /// follow (`hsp::mask_query_row`): 1 for a protein or plus-strand nucleotide row, 3 for a
    /// translated plus frame, -3 for a minus frame (`query_first` is then the highest
    /// nucleotide of the first residue).
    pub query_first: i32,
    pub query_step: i32,
}
