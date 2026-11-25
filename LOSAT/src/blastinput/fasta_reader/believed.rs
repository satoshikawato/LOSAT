//! `-parse_deflines` (`CBlastInputSourceConfig::SetBelieveDeflines`, `fParseRawID`): a record's
//! first word is its Seq-ids, and the rest of its defline its title.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:319-322
//! The ids are parsed by the reader of `makeblastdb -parse_seqids` (`process_ids`, with the same
//! `fQuickIDCheck` and `fAssumeNuc` or `fAssumeProt`).
//!
//! What NCBI's object manager does with these ids is not reproduced where it shows another
//! sequence than the record's (`dropin-survey-20261029/DI7/DECISIONS.md` 8): LOSAT refuses
//! such input. NCBI keeps the queries of a batch in one scope (every batch's with `-subject` and in
//! outfmt 5) and finds a sequence again by its id, so of two queries with one id there it searches
//! the first's letters twice, and it shows a database sequence with the id of such a query as that
//! query (probe of 2026-10-30; per batch: series Q of line RG, Q1c, `query_batch::QueryScopes`).
//! When an id's accession names the other molecule (`named_molecule`), LOSAT reproduces NCBI's
//! deterministic reading (`fasta_reader::reader::scan_data_line`, line RG-XT item X7-5).

use std::collections::HashSet;

use crate::report::archive::bioseq::{QuerySeqLoc, SeqLocStrand};
use crate::seqid::accguide::{identify_seq_id, parse_flags, GENERAL, LOCAL};
use crate::seqid::defline_set::{scope_keys, tabular_id, TabularId};
use crate::seqid::id_scope::{name_keys, version_keys, SharedRecord};
use crate::seqid::label;
use crate::seqid::model::{Dbtag, ObjectId, ParsedSeqId};
use crate::seqid::rank::{best_rank_score, find_best, text_score};
use crate::writedb::fasta_source::defline_ids::{parse_range, post_seq_id_warning_of, process_ids};
use crate::writedb::fasta_source::reader::ReaderWarning;

/// The Seq-ids of a record that `-parse_deflines` read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BelievedIds {
    /// `CBioseq::GetId()`: the Bioseq of the archive.
    pub ids: Vec<ParsedSeqId>,
    /// `FindBestChoice` of `ids` with `BestRank`: the id of the record's Seq-aligns.
    pub best: ParsedSeqId,
    /// `qseqid`, `qacc` and `qaccver` (`s_GetSeqIdListString`).
    pub tabular: [Vec<u8>; 3],
    /// `qgi`.
    pub gi: i64,
    /// `id->GetSeqIdString()`: the id of the record in NCBI's warnings.
    pub warning_id: Vec<u8>,
    /// `IsLocalId` of `best` (blast_aux.cpp:1094-1108): the archive then holds the Bioseqs.
    pub local: bool,
    /// The strand of the record's Seq-loc.
    pub strand: SeqLocStrand,
    /// Whether the record is one that a reader of the other role read before it
    /// (`crate::seqid::id_scope::IdScope::add_record`).
    pub shared: Option<SharedRecord>,
    /// Whether the reports show these ids as `-parse_deflines` shows a record's ids: false for the
    /// ids of a Bioseq that the reports show as a record without them, which only the record's
    /// Seq-id takes (`FastaRecord::keeping_bioseq_ids`, `RecordNames::seq_id`).
    pub shown: bool,
    /// `CFastaReader::xSetSeqMol`: whether the defline's ids named protein (`Some(true)`),
    /// nucleotide (`Some(false)`), or neither (`None`) when the record was read.
    pub named_molecule: Option<bool>,
}

impl BelievedIds {
    /// The ids `ids` of a record.
    pub(crate) fn new(ids: Vec<ParsedSeqId>, strand: SeqLocStrand) -> Self {
        let best = find_best(&ids, best_rank_score)
            .cloned()
            .expect("a defline parse gives at least one id");
        Self::with_alignment_id(ids, best, strand)
    }

    /// Records the molecule `CFastaReader::xSetSeqMol` inferred from the defline's parsed ids
    /// (before `PostProcessIDs` replaces ranged ids with generated ones).
    pub(crate) fn with_named_molecule(mut self, named_molecule: Option<bool>) -> Self {
        self.named_molecule = named_molecule;
        self
    }

    /// The ids `ids` (at least one) of a Bioseq that is the query as it is (the PSSM's of
    /// `-in_pssm`): its Seq-aligns and warnings take the first, of which NCBI makes the query's
    /// Seq-loc (psiblast_app.cpp:540, where `SetWhole().Assign` takes the first id of the PSSM's query
    /// Seq), where a FASTA query's Seq-loc takes `FindBestChoice` of the ids it read
    /// (`new`).
    pub(crate) fn of_bioseq(ids: Vec<ParsedSeqId>, strand: SeqLocStrand) -> Self {
        let first = ids.first().cloned().expect("a Bioseq with a Seq-id");
        Self::with_alignment_id(ids, first, strand)
    }

    fn with_alignment_id(ids: Vec<ParsedSeqId>, best: ParsedSeqId, strand: SeqLocStrand) -> Self {
        let field = |kind| {
            let mut out = Vec::new();
            tabular_id(&ids, kind, &mut out);
            out
        };
        let tabular = [
            field(TabularId::Full),
            field(TabularId::Accession),
            field(TabularId::AccessionVersion),
        ];
        let gi = label::first_gi(&ids);
        let mut warning_id = Vec::new();
        label::content(&best, false, &mut warning_id);
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1102-1106
        let info = identify_seq_id(
            &best,
            parse_flags::ANY_RAW | parse_flags::ANY_LOCAL,
            &mut |_| {},
        );
        let local = matches!(best, ParsedSeqId::Local(_)) || info == LOCAL || info == GENERAL;
        Self {
            ids,
            best,
            tabular,
            gi,
            warning_id,
            local,
            strand,
            shared: None,
            shown: true,
            named_molecule: None,
        }
    }

    /// `CAlignFormatUtil::GetSeqIdString` of `ids` with `true`: the id of `Query=` and of the XML query.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4944-4956
    /// The record as the Seq-loc of the archive's `seq-loc-list`, `None` for a local id.
    pub(crate) fn seq_loc(&self, from: i32, to: i32) -> Option<QuerySeqLoc<'_>> {
        (!self.local).then_some(QuerySeqLoc {
            from,
            to,
            strand: self.strand,
            id: &self.best,
        })
    }

    /// Whether these are the ids of a `-subject` record that the reports show as a record without
    /// `-parse_deflines` shows (a local id `Subject_N`, written or generated): the branch of
    /// `CShowBlastDefline::GetSeqIdList` on `lcl|Subject_` puts the first word of the title in its
    /// place. A `gnl|BL_ORD_ID|n` subject takes that branch too, and `duplicate_id` refuses it.
    pub(crate) fn subject_local(&self) -> bool {
        matches!(
            self.ids.as_slice(),
            [ParsedSeqId::Local(ObjectId::Str(name))] if name.starts_with(b"Subject_")
        )
    }

    pub(crate) fn bare_id(&self) -> String {
        let mut out = Vec::new();
        if let Some(id) = find_best(&self.ids, text_score) {
            label::bare(id, &mut out);
        }
        String::from_utf8_lossy(&out).into_owned()
    }
}

/// A record's defline as `-parse_deflines` reads it: where its title starts, and its ids.
#[derive(Debug)]
pub(crate) struct BelievedDefline {
    /// The position of the first byte after the first word (the title's search starts here).
    pub title_search: usize,
    /// The ids of the first word; empty for a blank defline (`GenerateID`).
    pub ids: Vec<ParsedSeqId>,
    pub range: bool,
}

/// `CFastaDeflineReader::ParseDefline` with `fParseRawID` and `fQuickIDCheck`. `Err` is the
/// message of the `CObjReaderParseException`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta_reader_utils.cpp:157-207
/// `defline[pos]` is a signed `char` on the platforms compared: a byte from 0x80 up ends the id.
pub(crate) fn parse_believed_defline(
    defline: &[u8],
    protein: bool,
    program: &str,
    line_number: u64,
    warn: &mut dyn FnMut(&[u8]),
) -> Result<BelievedDefline, String> {
    use crate::ncbi_str::{is_blank, is_space};
    let len = defline.len();
    if len <= 1 || is_blank(&defline[1..]) {
        return Ok(BelievedDefline {
            title_search: len,
            ids: Vec::new(),
            range: false,
        });
    }
    let start = 1 + defline[1..]
        .iter()
        .position(|&byte| !is_space(byte))
        .unwrap_or(len - 1);
    let mut pos = start;
    while pos < len && (defline[pos] as i8) > (b' ' as i8) {
        pos += 1;
    }
    let (range_len, _, _) = parse_range(&defline[start..pos]);
    let id_string = &defline[start..pos - range_len];
    if is_blank(id_string) {
        return Err("Unable to locate sequence id in definition line".to_string());
    }
    // A warning of the Seq-id parser names the application (`post_seq_id_warning_of`); the
    // reader's own warnings are lines without a prefix.
    let name = program.to_ascii_lowercase();
    let mut warn_named = |warning: ReaderWarning<'_>| match warning {
        ReaderWarning::Line(line) => warn(line),
        ReaderWarning::SeqId(warning) => post_seq_id_warning_of(&name, warning, &mut *warn),
    };
    let ids = process_ids(id_string, protein, true, true, line_number, &mut warn_named)
        .map_err(|error| error.message)?;
    Ok(BelievedDefline {
        title_search: pos,
        ids,
        range: range_len > 0,
    })
}

/// LOSAT's refusal of an id of a record, or `None`: a `gnl|BL_ORD_ID|n` id and a Geninfo import id.
/// `origin` names where the id comes from (`-parse_deflines`, the PSSM of `-in_pssm`, or the
/// master sequence of psiblast's `-in_msa` alignment).
pub(crate) fn refused_id(
    ids: &[ParsedSeqId],
    role: &str,
    origin: &str,
    program: &str,
) -> Option<anyhow::Error> {
    let fasta = |id: &ParsedSeqId| {
        let mut out = Vec::new();
        label::fasta(id, &mut out);
        String::from_utf8_lossy(&out).into_owned()
    };
    for id in ids {
        if matches!(id, ParsedSeqId::General(Dbtag { db, .. }) if db.as_slice() == b"BL_ORD_ID") {
            return Some(anyhow::anyhow!(
                "the {role} Seq-id {} of {origin} names a database ordinal; this is not \
                 supported by LOSAT's {program}",
                fasta(id)
            ));
        }
    }
    for id in ids {
        // NCBI stops with `Cannot resolve gim|n` (probe of the review of DECISIONS 15).
        if matches!(id, ParsedSeqId::Giim(_)) {
            return Some(anyhow::anyhow!(
                "the {role} Seq-id {} of {origin} is a Geninfo import id, which NCBI BLAST+ \
                 cannot resolve; this is not supported by LOSAT's {program}",
                fasta(id)
            ));
        }
    }
    None
}

/// LOSAT's refusal of the ids of a record, or `None`: an id that a record of the role read before
/// has: the subjects of a reader, the queries of a scope (`query_batch::QueryScopes`). The readers
/// refuse what `refused_id` refuses first. An accession with a version and the same accession without one
/// count as one id for two subjects, and in BLASTN two queries (`dropin-survey-20261029/DI7/
/// DECISIONS.md` 11, 12). Adds the keys of the ids. A record of the other role is
/// `crate::seqid::id_scope::IdScope::add_record`'s.
pub(crate) fn duplicate_id(
    ids: &[ParsedSeqId],
    seen: &mut HashSet<u64>,
    role: &str,
    program: &str,
) -> Option<anyhow::Error> {
    let fasta = |id: &ParsedSeqId| {
        let mut out = Vec::new();
        label::fasta(id, &mut out);
        String::from_utf8_lossy(&out).into_owned()
    };
    // Where two queries are in one scope (`query_batch::QueryScopeSpan`): psiblast's batches
    // hold one query, after each of which NCBI's formatter empties the scope against a database
    // (`CBlastFormat::ResetScopeHistory`), so only `-subject` keeps two, where NCBI BLAST+ 2.17.0
    // searches the first's letters for both or crashes (review of RG Q1, S1; series Q, Q1c; review
    // 5 of line WQA D10, L3: g2a with `-subject` exits 1 with the reader's error after two searched
    // batches, g2b crashes, exit 139).
    let scope = match (role, program) {
        ("query", "PSIBLAST") => " with -subject",
        ("query", _) => " in one query batch, or in any with -subject or outfmt 5",
        _ => "",
    };
    let mut keys = Vec::new();
    for id in ids {
        keys.clear();
        scope_keys(id, &mut keys);
        if keys.iter().any(|key| seen.contains(key)) {
            let why = match (role, program) {
                ("query", "PSIBLAST") => {
                    "NCBI BLAST+'s psiblast searches the letters of the first for both, or crashes"
                }
                ("query", _) => {
                    "NCBI BLAST+'s scope holds both and searches the letters of the first for both"
                }
                _ => "NCBI BLAST+ shows the first for both",
            };
            return Some(anyhow::anyhow!(
                "two {role} records have the Seq-id {} under -parse_deflines{scope} ({why}); this \
                 is not supported by LOSAT's {program}",
                fasta(id)
            ));
        }
        seen.extend(keys.iter().copied());
        if let Some((own, other)) = name_keys(id) {
            if seen.contains(&other) {
                return Some(anyhow::anyhow!(
                    "two {role} records have the name of {} with and without an accession under \
                     -parse_deflines{scope} (NCBI BLAST+ shows them as one sequence); this is not \
                     supported by LOSAT's {program}",
                    fasta(id)
                ));
            }
            seen.insert(own);
        }
        if program == "BLASTN" || role == "subject" {
            if let Some((own, other)) = version_keys(id) {
                if seen.contains(&other) {
                    let why = if role == "query" {
                        "NCBI BLAST+'s BLASTN writes a warning of its scope"
                    } else {
                        "NCBI BLAST+ shows them as one sequence"
                    };
                    return Some(anyhow::anyhow!(
                        "two {role} records have the accession of {} with and without a version \
                         under -parse_deflines{scope} ({why}); this is not supported by LOSAT's \
                         {program}",
                        fasta(id)
                    ));
                }
                seen.insert(own);
            }
        }
    }
    None
}

/// With outfmt 0, LOSAT's refusal of a query that is also a subject (`BelievedIds::shared`) and has
/// masks (`masked`, by query; with `-lcase_masking` the lower-case letters of either copy too):
/// `CDisplaySeqalign::x_FillLocList` finds a mask's row by its id, and the two rows have one id, so
/// it draws the mask at the other row's places (`dropin-survey-20261029/DI7/DECISIONS.md` 19).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4292-4306
pub(crate) fn refuse_shared_masks(
    queries: &[super::FastaRecord],
    masked: impl Fn(usize) -> bool,
    lcase_masking: bool,
    program: &str,
) -> anyhow::Result<()> {
    for (index, record) in queries.iter().enumerate() {
        let Some(ids) = record.believed() else {
            continue;
        };
        let Some(shared) = ids.shared else {
            continue;
        };
        if masked(index) || (lcase_masking && shared.lowercase) {
            let mut fasta = Vec::new();
            label::fasta(&ids.best, &mut fasta);
            anyhow::bail!(
                "the record {} is both a query and a subject under -parse_deflines and has masks \
                 (NCBI BLAST+ draws them at the other record's places in outfmt 0); this is not \
                 supported by LOSAT's {program}",
                String::from_utf8_lossy(&fasta)
            );
        }
    }
    Ok(())
}

/// `refuse_shared_masks` as the refusal of outfmt 0 alone (`FormatRefusal::check`), in each
/// outfmt 0 output of a call (`refusals`), before any format of the call is written: the command
/// line's run of outfmt 0 stops with it, and a call of several formats in which outfmt 0 fails
/// alone (LOSAT Web's outfmt 0 when the argv asks for another) writes the others, as the command
/// line's run of each of them does (line WQA's sweep, D7, 2026-11-19). BLASTN, BLASTP and TBLASTX
/// make it so.
pub(crate) fn refuse_shared_masks_in<'r>(
    refusals: impl IntoIterator<Item = &'r mut crate::api::local_blast::FormatRefusal>,
    queries: &[super::FastaRecord],
    masked: impl Fn(usize) -> bool,
    lcase_masking: bool,
    program: &str,
) -> anyhow::Result<()> {
    for refusal in refusals {
        refusal.check(refuse_shared_masks(
            queries,
            &masked,
            lcase_masking,
            program,
        ))?;
    }
    Ok(())
}

/// `CFastaReader::xSetSeqMol`: whether the first id whose accession names a molecule names a
/// protein, `None` when none does.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:547-562
pub(crate) fn named_molecule(ids: &[ParsedSeqId]) -> Option<bool> {
    use crate::seqid::accguide::{F_NUC, F_PROT};
    for id in ids {
        let info = identify_seq_id(
            id,
            parse_flags::ANY_RAW | parse_flags::ANY_LOCAL,
            &mut |_| {},
        );
        if info & F_NUC != 0 {
            return Some(false);
        }
        if info & F_PROT != 0 {
            return Some(true);
        }
    }
    None
}

/// `CFastaReader::x_NucOrProt`, in the reader's message about invalid letters: the molecule the
/// accession named, nothing without one.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:2117-2125
pub(crate) fn nuc_or_prot(named: Option<bool>) -> &'static str {
    match named {
        Some(true) => "protein ",
        Some(false) => "nucleotide ",
        None => "",
    }
}

/// LOSAT's refusal of a first word with a range (`x:10-40`) in a run that writes the BLAST archive:
/// NCBI's Bioseq holds a Seq-hist alignment to the old id (`CFastaReader::xCreateAlignment`), which
/// LOSAT's archive does not write.
pub(crate) fn refuse_range_in_archive(role: &str, program: &str) -> anyhow::Error {
    anyhow::anyhow!(
        "a {role} Seq-id with a range under -parse_deflines in the BLAST archive is not supported \
         by LOSAT's {program}"
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:633-640
/// Whether the defline `defline` read under `-parse_deflines` has a first word with a range:
/// what `refuse_range_in_archive` refuses, for a caller whose reader did not (`archive_known`).
pub(crate) fn defline_has_range(defline: &[u8], protein: bool, program: &'static str) -> bool {
    parse_believed_defline(defline, protein, program, 0, &mut |_| {})
        .is_ok_and(|parsed| parsed.range)
}

/// The local id `id_prefix` + `n` that `GenerateID` makes for a record without ids.
pub(crate) fn generated_ids(id: &str) -> Vec<ParsedSeqId> {
    vec![ParsedSeqId::Local(ObjectId::Str(id.as_bytes().to_vec()))]
}

#[cfg(test)]
#[path = "believed_tests.rs"]
mod tests;
