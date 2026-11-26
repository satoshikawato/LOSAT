//! The report records of outfmt 5: per query, borrowed views of a program's final hit lists in
//! the terms of NCBI's `CBlastOutput` model (objects/blastxml/blastxml.asn), built at write time
//! into buffers reused from query to query. Nothing here flows through the search's pools,
//! pipes or heaps; each program's adapter (`algorithm/<program>/blast_xml.rs`) builds them.

use crate::common::GapEditOp;

use super::matrix::FormattingMatrix;
use super::subjects::SubjectNames;

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
    /// `Statistics_db-num` and `Statistics_db-len`: the database's sequences and letters (the
    /// sums over the names of the `-db` list, `PairwiseDatabase`); 0 and 0 for `-subject`. The
    /// number is the `Int8` sum that XML2 and JSON write (data4xml2format.hpp:215-217); BLAST
    /// XML sums NCBI's `int` (data4xmlformat.hpp:276), equal where that sum does not overflow,
    /// which LOSAT refuses (`XmlRun::check_sums`).
    pub db_num: i64,
    pub db_len: i64,
    /// Whether the sum of `db_len` passes NCBI's `Int8` (`PairwiseDatabase::total_letters_overflow`),
    /// which LOSAT refuses when it writes a document (`XmlRun::check_sums`).
    pub db_len_overflow: bool,
    /// The positives table of the recount; `None` for BLASTN (a nucleotide search shows `|`).
    pub positives: Option<&'a FormattingMatrix>,
    /// The letter of a masked query position: `N` for BLASTN, `X` otherwise
    /// (blastxml_format.cpp:280-283).
    pub mask_char: u8,
}

impl XmlRun<'_> {
    /// LOSAT's refusal of the database counts that a report of `outfmt` (5, 13 to 16) of this
    /// run would print where NCBI's sums overflow (`common::database::check_sums_of`), made by
    /// the documents' `StructuredXml::start` before their first byte.
    pub(crate) fn check_sums(&self, outfmt: i32) -> anyhow::Result<()> {
        crate::algorithm::common::database::check_sums_of(
            self.db_num,
            self.db_len_overflow,
            outfmt,
            &self.program.to_ascii_uppercase(),
        )
    }
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
    /// Whether NCBI's result of the query has an alignment set: its batch was searched (it has a
    /// valid context in it), or the merge of several RPS volumes made one for every query
    /// (rpsblast_local.cpp:139-185). The XML2 `Search` of a query without one has no `hits`
    /// member and adds ` \nNo hits found` to its message (blastxml2_format.cpp:614-622);
    /// outfmt 5 does not read it.
    pub searched: bool,
    /// The query's masked regions in NCBI's order and frames (`TMaskedQueryRegions`):
    /// outfmt 5 masks an HSP's row with the regions of the HSP's own frame, XML2 lists every
    /// region of a plus frame in `Search_query-masking` and masks with the regions of frame 0
    /// or the HSP's translated frame.
    pub masks: &'a [QueryMask],
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2363-2389
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xml2format.cpp:95-110
    /// Whether the result joins the held report of the queries before it where the XML2 and JSON
    /// documents hold their reports (`blast_xml2::Xml2Document::holds_back`; the one owner's
    /// answer, `query_batch::XmlShow::joins_held`): NCBI adds the result's `Search` to that
    /// report, shown with the id, title and length of the report's first `Search` and masked by
    /// its masks. Read by no other document.
    pub joins_held: bool,
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
    /// What the taxonomy database said about the taxid (XML2's `HitDescr_sciname`); outfmt 5
    /// does not read it.
    pub taxonomy: &'a HitTaxonomy,
    /// The `HitDescr` list of a database sequence with parsed Seq-ids (`SubjectNames::
    /// descriptions`; empty for the others).
    pub descriptions: &'a [crate::report::db_subject::HitDescription],
    /// `GenerateDefline` threw for the `-subject` record (`SubjectNames::defline_thrown`):
    /// outfmt 5's `Hit_len` is 0.
    pub defline_thrown: bool,
    /// The exception that the lookup of the database sequence's Seq-align id throws where the
    /// `Hit` loads it (`SubjectNames::lookup_throws`; the names are those of NCBI's `catch`,
    /// `SubjectNames::of_thrown`): outfmt 5 writes the `catch`'s `Hit` and HSPs
    /// (`writer::Xml5Hits`), the XML2 and JSON documents raise it (`CaughtException`).
    pub lookup_throws: Option<&'a crate::ncbi_exception::NcbiException>,
}

impl<'a> XmlHit<'a> {
    /// The `Hit` of a subject with the names `names` and the record length `len`.
    pub(crate) fn of(names: &'a SubjectNames, len: i32) -> Self {
        XmlHit {
            id: &names.id,
            title: &names.title,
            accession: &names.accession,
            len,
            taxid: names.taxid,
            taxonomy: &names.taxonomy,
            descriptions: &names.descriptions,
            defline_thrown: names.defline_thrown,
            lookup_throws: names.lookup_throws.as_ref(),
        }
    }
}

/// What NCBI's XML2 and JSON formats do with the taxid of a `Hit` once the taxonomy database has
/// been asked (`s_SeqAlignToXMLHit`, blastxml2_format.cpp:380-386; DECISIONS 28).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:380-386
/// (`GetTaxInfo` throws when the taxid is not found; the `catch` of the hit replaces its
/// description by the sequence id alone.)
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum HitTaxonomy {
    /// No scientific name is written: no taxid, `-subject`, or no `taxdb.bti` on the search
    /// path.
    #[default]
    Unnamed,
    /// The scientific name, written after the taxid (an empty one too).
    Named(Vec<u8>),
    /// The taxonomy database does not know the taxid (or is damaged): the description is the id
    /// alone. The texts are the `ERR_POST` messages of this hit's lookup (NCBI writes them again
    /// at every hit).
    NotFound(Vec<String>),
}

/// The `taxonomy` of a hit without a taxonomy lookup.
#[cfg(test)]
pub const UNNAMED: &HitTaxonomy = &HitTaxonomy::Unnamed;

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
    /// The edit script of the rows (`None`: an ungapped HSP, one segment of substitutions over the
    /// whole row and no gaps), which gives the Dense-seg segments that the mask loop walks and the
    /// gaps (`writer::count_hsp`).
    pub ops: Option<&'a [GapEditOp]>,
    /// The 0-based record position of the query row's first residue and how the residues
    /// follow (`hsp::mask_query_row`): 1 for a protein or plus-strand nucleotide row, 3 for a
    /// translated plus frame, -3 for a minus frame (`query_first` is then the highest
    /// nucleotide of the first residue).
    pub query_first: i32,
    pub query_step: i32,
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:810-849
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1199-1203, 1491-1495
    /// The Seq-align's `num_ident` as `GetAlnScores` reads it: the engine's identities, -1 where
    /// the Seq-align has no `num_ident` score (psiblast's HSP of 0 identities, whose count NCBI
    /// sets to -1 before it writes the Seq-align). Outfmt 5 writes it in the `catch` of an HSP
    /// whose subject's lookup throws; every other HSP counts its identities again from the rows.
    pub alignment_identities: i32,
}

impl XmlHsp<'_> {
    /// An HSP of a `Hit` whose writer does not read its alignment (`HspSink::reads_alignment`:
    /// the HSPs of a `Hit` whose subject's lookup throws in outfmt 5): its Seq-align scores alone
    /// (the raw score, the bit score, the e-value as the Seq-align has it, and
    /// `alignment_identities`), without coordinates, frames or rows.
    pub fn unaligned(score: i32, bit_score: f64, evalue: f64, alignment_identities: i32) -> Self {
        XmlHsp {
            score,
            bit_score,
            evalue,
            query_from: 0,
            query_to: 0,
            hit_from: 0,
            hit_to: 0,
            query_frame: 0,
            hit_frame: 0,
            qrow: &[],
            srow: &[],
            ops: None,
            query_first: 0,
            query_step: 1,
            alignment_identities,
        }
    }
}
