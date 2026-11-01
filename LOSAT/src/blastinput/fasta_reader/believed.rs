//! `-parse_deflines` (`CBlastInputSourceConfig::SetBelieveDeflines`, `fParseRawID`): a record's
//! first word is its Seq-ids, and the rest of its defline its title.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:319-322
//! ```c++
//!     CFastaReader::TFlags flags = m_Config.GetBelieveDeflines() ?
//!                                     CFastaReader::fParseRawID:
//!                                     (CFastaReader::fNoParseID |
//!                                      CFastaReader::fDLOptional);
//! ```
//! The ids are parsed by the reader of `makeblastdb -parse_seqids` (`process_ids`, with the same
//! `fQuickIDCheck` and `fAssumeNuc` or `fAssumeProt`).
//!
//! What NCBI's object manager does with these ids is not reproduced where it shows another
//! sequence than the record's (`dropin-survey-20261029/DI7/DECISIONS.md` 7 and 8): LOSAT refuses
//! such input. NCBI keeps every query in one scope and finds a sequence again by its id, so of
//! two queries with one id it searches the first's letters twice, and it shows a database
//! sequence with the id of a query as that query (probe of 2026-10-30).

use std::collections::HashSet;

use crate::report::archive::bioseq::{QuerySeqLoc, SeqLocStrand};
use crate::seqid::accguide::{identify_seq_id, parse_flags, GENERAL, LOCAL};
use crate::seqid::defline_set::{scope_keys, tabular_id, TabularId};
use crate::seqid::id_scope::{name_keys, version_keys, SharedRecord};
use crate::seqid::label;
use crate::seqid::model::{Dbtag, ObjectId, ParsedSeqId};
use crate::seqid::rank::{best_rank_score, find_best, text_score};
use crate::writedb::fasta_source::defline_ids::{parse_range, process_ids};

/// The Seq-ids of a record that `-parse_deflines` read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BelievedIds {
    /// `CBioseq::GetId()`: the Bioseq of the archive.
    pub ids: Vec<ParsedSeqId>,
    /// `FindBestChoice(ids, BestRank)`: the id of the record's Seq-aligns.
    pub best: ParsedSeqId,
    /// `qseqid`, `qacc` and `qaccver` (`s_GetSeqIdListString`).
    pub tabular: [Vec<u8>; 3],
    /// `qgi`.
    pub gi: i64,
    /// `id->GetSeqIdString()`: the id of the record in NCBI's warnings.
    pub warning_id: Vec<u8>,
    /// `IsLocalId(best)` (blast_aux.cpp:1094-1108): the archive then holds the Bioseqs.
    pub local: bool,
    /// The strand of the record's Seq-loc.
    pub strand: SeqLocStrand,
    /// Whether the record is one that a reader of the other role read before it
    /// (`crate::seqid::id_scope::IdScope::add_record`).
    pub shared: Option<SharedRecord>,
}

impl BelievedIds {
    /// The ids `ids` of a record.
    pub(crate) fn new(ids: Vec<ParsedSeqId>, strand: SeqLocStrand) -> Self {
        let best = find_best(&ids, best_rank_score)
            .cloned()
            .expect("a defline parse gives at least one id");
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
        // ```c++
        //     CSeq_id::EAccessionInfo id_info = seqid->IdentifyAccession();
        //     if (seqid->IsLocal() || id_info == CSeq_id::eAcc_local  ||
        //     	id_info == CSeq_id::eAcc_general) {
        //         retval = true;
        //     }
        // ```
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
        }
    }

    /// `CAlignFormatUtil::GetSeqIdString(ids, true)`: the id of `Query=` and of the XML query.
    ///
    /// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4944-4956
    /// ```c++
    /// string CAlignFormatUtil::GetBareId(const CSeq_id& id)
    /// {
    ///     string retval;
    ///
    ///     if (id.IsGi() || id.IsPrf() || id.IsPir()) {
    ///         retval = id.AsFastaString();
    ///     }
    ///     else {
    ///         retval = id.GetSeqIdString(true);
    ///     }
    /// ```
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
/// ```c++
///     const size_t len = defline.length();
///     if (len <= 1 ||
///         NStr::IsBlank(defline.substr(1))) {
///         return;
///     }
///     ...
///     // ignore spaces between '>' and the sequence ID
///     size_t start;
///     for(start = 1 ; start < len; ++start ) {
///         if( ! isspace(defline[start]) ) {
///             break;
///         }
///     }
///     ...
///         pos = start;
///         while (pos < len && defline[pos] > ' ') {
///             pos++;
///         }
///
///         if ( ! (fFastaFlags & CFastaReader::fDisableParseRange) ) {
///             range_len = ParseRange(defline.substr(start, pos - start),
///                     data.range_start, data.range_end, pMessageListener);
///         }
///
///         auto id_string = defline.substr(start, pos - start - range_len);
///         if (NStr::IsBlank(id_string)) {
///             NCBI_THROW2(CObjReaderParseException, eFormat,
///                 "Unable to locate sequence id in definition line", 0);
///         }
/// ```
/// `defline[pos]` is a signed `char` on the platforms compared: a byte from 0x80 up ends the id.
pub(crate) fn parse_believed_defline(
    defline: &[u8],
    protein: bool,
    program: &str,
    line_number: u64,
    warn: &mut dyn FnMut(&[u8]),
) -> Result<BelievedDefline, String> {
    let space = crate::seqid::parse::is_space;
    let len = defline.len();
    if len <= 1 || defline[1..].iter().all(|&byte| space(byte)) {
        return Ok(BelievedDefline {
            title_search: len,
            ids: Vec::new(),
            range: false,
        });
    }
    let start = 1 + defline[1..]
        .iter()
        .position(|&byte| !space(byte))
        .unwrap_or(len - 1);
    let mut pos = start;
    while pos < len && (defline[pos] as i8) > (b' ' as i8) {
        pos += 1;
    }
    let (range_len, _, _) = parse_range(&defline[start..pos]);
    let id_string = &defline[start..pos - range_len];
    if id_string.iter().all(|&byte| space(byte)) {
        return Err("Unable to locate sequence id in definition line".to_string());
    }
    // NCBI's diagnostics name the application; `process_ids` is makeblastdb's.
    let name = program.to_ascii_lowercase();
    let mut warn_named = |message: &[u8]| {
        let text = String::from_utf8_lossy(message)
            .replace("[makeblastdb]", &format!("[{name}]"))
            .into_bytes();
        warn(&text);
    };
    let ids = process_ids(id_string, protein, line_number, &mut warn_named)
        .map_err(|error| error.message)?;
    Ok(BelievedDefline {
        title_search: pos,
        ids,
        range: range_len > 0,
    })
}

/// LOSAT's refusal of the ids of a record, or `None`: a `gnl|BL_ORD_ID|n` id and an id that a record
/// of the role read before has. An accession with a version and the same accession without one
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
    for id in ids {
        if matches!(id, ParsedSeqId::General(Dbtag { db, .. }) if db.as_slice() == b"BL_ORD_ID") {
            return Some(anyhow::anyhow!(
                "the {role} Seq-id {} of -parse_deflines names a database ordinal; this is not \
                 supported by LOSAT's {program}",
                fasta(id)
            ));
        }
    }
    for id in ids {
        // NCBI stops with `Cannot resolve gim|n` (probe of the review of DECISIONS 15).
        if matches!(id, ParsedSeqId::Giim(_)) {
            return Some(anyhow::anyhow!(
                "the {role} Seq-id {} of -parse_deflines is a Geninfo import id, which NCBI BLAST+ \
                 cannot resolve; this is not supported by LOSAT's {program}",
                fasta(id)
            ));
        }
    }
    let mut keys = Vec::new();
    for id in ids {
        keys.clear();
        scope_keys(id, &mut keys);
        if keys.iter().any(|key| seen.contains(key)) {
            return Some(if role == "query" {
                anyhow::anyhow!(
                    "two query records have the Seq-id {} under -parse_deflines (in one query \
                     batch NCBI BLAST+ searches the letters of the first for both); this is not \
                     supported by LOSAT's {program}",
                    fasta(id)
                )
            } else {
                anyhow::anyhow!(
                    "two subject records have the Seq-id {} under -parse_deflines (NCBI BLAST+ \
                     shows the first for both); this is not supported by LOSAT's {program}",
                    fasta(id)
                )
            });
        }
        seen.extend(keys.iter().copied());
        if let Some((own, other)) = name_keys(id) {
            if seen.contains(&other) {
                return Some(anyhow::anyhow!(
                    "two {role} records have the name of {} with and without an accession under \
                     -parse_deflines (NCBI BLAST+ shows them as one sequence); this is not \
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
                         under -parse_deflines ({why}); this is not supported by LOSAT's {program}",
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
/// ```c++
/// void CDisplaySeqalign::x_FillLocList(TSAlnSeqlocInfoList& loc_list,
///                                      const list< CRef<CSeqLocInfo> >*  masks) const
/// ```
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

/// `CFastaReader::xSetSeqMol`: whether the first id whose accession names a molecule names a
/// protein, `None` when none does.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:547-562
/// ```c++
/// bool CFastaReader::xSetSeqMol(const list<CRef<CSeq_id>>& ids, CSeq_inst_Base::EMol& mol) {
///
///     for (auto pId : ids) {
///         const CSeq_id::EAccessionInfo acc_info = pId->IdentifyAccession();
///         if (acc_info & CSeq_id::fAcc_nuc) {
///             mol = CSeq_inst::eMol_na;
///             return true;
///         }
///
///         if (acc_info & CSeq_id::fAcc_prot) {
///             mol = CSeq_inst::eMol_aa;
///             return true;
///         }
///     }
///     return false;
/// }
/// ```
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
/// ```c++
/// std::string CFastaReader::x_NucOrProt(void) const
/// {
///     if( m_CurrentSeq && m_CurrentSeq->IsSetInst() &&
///         m_CurrentSeq->GetInst().IsSetMol() )
///     {
///         return ( m_CurrentSeq->GetInst().IsAa() ? "protein " : "nucleotide " );
///     } else {
///         return kEmptyStr;
///     }
/// }
/// ```
pub(crate) fn nuc_or_prot(named: Option<bool>) -> &'static str {
    match named {
        Some(true) => "protein ",
        Some(false) => "nucleotide ",
        None => "",
    }
}

/// Whether both readings of a data line give the same letters: those of a molecule named by the
/// accession and those of the program's. The nucleotide reading drops what it does not know and
/// reads `U` as `T`; every letter here reads the same in both, and blanks and digits are skipped by
/// both. A line with another byte is LOSAT's refusal (`dropin-survey-20261029/DI7/DECISIONS.md` 7).
pub(crate) fn reads_the_same_in_both(line: &[u8]) -> bool {
    line.iter().all(|&byte| {
        byte.is_ascii_whitespace()
            || byte.is_ascii_digit()
            || matches!(
                byte.to_ascii_uppercase(),
                b'A' | b'C'
                    | b'G'
                    | b'T'
                    | b'R'
                    | b'Y'
                    | b'M'
                    | b'K'
                    | b'W'
                    | b'S'
                    | b'B'
                    | b'D'
                    | b'H'
                    | b'V'
                    | b'N'
            )
    })
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

/// The local id `id_prefix` + `n` that `GenerateID` makes for a record without ids.
pub(crate) fn generated_ids(id: &str) -> Vec<ParsedSeqId> {
    vec![ParsedSeqId::Local(ObjectId::Str(id.as_bytes().to_vec()))]
}

#[cfg(test)]
mod tests {
    use super::super::{read_all, FastaInputSource, ReaderConfig};
    use super::*;
    use crate::seqid::id_scope::IdScope;

    fn read(fasta: &[u8], protein: bool) -> Result<Vec<super::super::FastaRecord>, String> {
        let config = ReaderConfig::query("BLASTN", protein, false).with_parse_deflines(true);
        let mut source = FastaInputSource::from_bytes(fasta, config);
        read_all(&mut source, &mut |_| Ok(())).map_err(|error| format!("{error:#}"))
    }

    fn ids_of(record: &super::super::FastaRecord) -> &BelievedIds {
        record.believed().expect("-parse_deflines keeps the ids")
    }

    // NCBI BLAST+ 2.17.0 -parse_deflines (probe `DI7/ws/pd_probe.sh`, `pdo/*`): `Query=`, `qseqid`,
    // `qacc`, `qaccver` and `qgi` of each defline.
    #[test]
    fn the_ids_of_a_defline_follow_the_ncbi_oracle() {
        for (defline, line, qseqid, qacc, qaccver, gi) in [
            (
                &b"gi|42|ref|NP_0042.1| t"[..],
                &b"NP_0042.1 t"[..],
                &b"gi|42|ref|NP_0042.1|"[..],
                &b"NP_0042"[..],
                &b"NP_0042.1"[..],
                42,
            ),
            (
                b"ref|NP_0043.1| t",
                b"NP_0043.1 t",
                b"ref|NP_0043.1|",
                b"NP_0043",
                b"NP_0043.1",
                0,
            ),
            (b"lcl|q3 t", b"q3 t", b"q3", b"q3", b"q3", 0),
            (
                b"qplain4 t",
                b"qplain4 t",
                b"qplain4",
                b"qplain4",
                b"qplain4",
                0,
            ),
            (
                b"gnl|qdb|q5 t",
                b"qdb:q5 t",
                b"gnl|qdb|q5",
                b"qdb:q5",
                b"qdb:q5",
                0,
            ),
            (b"gi|0 t", b"gi|0 t", b"gi|0", b"0", b"0", 0),
            (
                b"prf||123456A t",
                b"prf||123456A t",
                b"prf||123456A",
                b"123456A",
                b"123456A",
                0,
            ),
            (
                b"pdb|1ABC|A t",
                b"1ABC_A t",
                b"pdb|1ABC|A",
                b"1ABC_A",
                b"1ABC_A",
                0,
            ),
            (b"q1", b"q1", b"q1", b"q1", b"q1", 0),
            (b"a,b t", b"a_b t", b"a_b", b"a_b", b"a_b", 0),
        ] {
            let fasta = [b">", defline, b"\nACGTACGT\n"].concat();
            let records = read(&fasta, false).unwrap();
            let ids = ids_of(&records[0]);
            let text = String::from_utf8_lossy(defline).into_owned();
            assert_eq!(&*records[0].query_line(), line, "{text}");
            assert_eq!(ids.tabular[0], qseqid, "{text}");
            assert_eq!(ids.tabular[1], qacc, "{text}");
            assert_eq!(ids.tabular[2], qaccver, "{text}");
            assert_eq!(ids.gi, gi, "{text}");
        }
    }

    // NCBI BLAST+ 2.17.0 (probe `edge3`, `blastn.gen.*`): the generated ids count the records
    // without ids, not the records.
    #[test]
    fn generated_ids_count_the_records_without_ids() {
        let records = read(b">has1 a\nACGT\n>\nACGT\n>has3 b\nACGT\n>\nACGT\n", false).unwrap();
        let lines: Vec<_> = records
            .iter()
            .map(|record| record.query_line().into_owned())
            .collect();
        assert_eq!(lines, [&b"has1 a"[..], b"Query_1", b"has3 b", b"Query_2"]);
    }

    // NCBI BLAST+ 2.17.0 (probe `edge2`, `range.*`): a range replaces the ids by a generated one.
    #[test]
    fn a_range_gives_a_generated_id_and_keeps_the_title() {
        let records = read(b">rq:10-40 t\nACGT\n", false).unwrap();
        assert_eq!(&*records[0].query_line(), b"Query_1 t");
        assert_eq!(records[0].title, b"t");
    }

    #[test]
    fn two_records_with_one_id_are_refused() {
        for fasta in [
            &b">dup1 a\nACGT\n>dup1 b\nACGT\n"[..],
            b">abc a\nACGT\n>ABC b\nACGT\n",
            b">gb|AB000001.1| a\nACGT\n>emb|AB000001.1| b\nACGT\n",
            b">gi|5|ref|NM_000001.1| a\nACGT\n>gi|5 b\nACGT\n",
            b">gb|AAA00007.1|NAMEA a\nACGT\n>gb||NAMEA b\nACGT\n",
            b">gb||NAMEA a\nACGT\n>gb|AAA00007.1|NAMEA b\nACGT\n",
        ] {
            let error = read(fasta, false).unwrap_err();
            assert!(error.contains("two query records"), "{error}");
        }
    }

    // NCBI BLAST+ 2.17.0 (probes of the review): NCBI's scope finds these as two sequences.
    #[test]
    fn a_version_or_a_chain_tells_two_ids_apart() {
        for fasta in [
            &b">NM_000001.1 a\nACGT\n>NM_000001.2 b\nACGT\n"[..],
            b">pdb|1ABC|A a\nACGT\n>pdb|1ABC|B b\nACGT\n",
            b">gb|AAA00007.1|NAMEA a\nACGT\n>gb|AAA00008.1|NAMEA b\nACGT\n",
            b">sp|P12345.2|NAME_X a\nACGT\n>sp|Q99999.1|NAME_X b\nACGT\n",
        ] {
            read(fasta, false).unwrap();
        }
    }

    // NCBI BLAST+ 2.17.0 (`db.pd.*.ver`): the same letters and reports, but BLASTN warns.
    #[test]
    fn an_accession_with_and_without_a_version_is_refused_by_blastn_only() {
        let fasta = b">NM_000005 a\nACGT\n>NM_000005.1 b\nACGT\n";
        let error = read(fasta, false).unwrap_err();
        assert!(error.contains("with and without a version"), "{error}");
        let config = ReaderConfig::query("TBLASTX", false, false).with_parse_deflines(true);
        let mut source = FastaInputSource::from_bytes(fasta, config);
        read_all(&mut source, &mut |_| Ok(())).unwrap();
    }

    // fasta.cpp:1005, 2117-2125: the message names the molecule that the accession named.
    #[test]
    fn the_invalid_letters_message_names_the_molecule_of_the_accession() {
        let message = |fasta: &[u8], protein: bool| {
            let config = ReaderConfig::query("BLASTN", protein, false).with_parse_deflines(true);
            let mut source = FastaInputSource::from_bytes(fasta, config);
            let mut text = Vec::new();
            read_all(&mut source, &mut |line| {
                text.extend_from_slice(line);
                Ok(())
            })
            .unwrap();
            String::from_utf8(text).unwrap()
        };
        assert!(message(b">NM_000009.1 x\nACGTACGTAC!GTACGTACGT\n", false)
            .starts_with("FASTA-Reader: Ignoring invalid nucleotide residues"));
        assert!(message(b">NP_000009.1 x\nMKLVMKLVMK!LVMKLVMKLV\n", true)
            .starts_with("FASTA-Reader: Ignoring invalid protein residues"));
        assert!(message(b">plain x\nACGTACGTAC!GTACGTACGT\n", false)
            .starts_with("FASTA-Reader: Ignoring invalid residues"));
    }

    #[test]
    fn a_geninfo_import_id_is_refused() {
        let error = read(b">gim|123 a\nACGT\n", false).unwrap_err();
        assert!(error.contains("Geninfo import id"), "{error}");
    }

    #[test]
    fn a_database_ordinal_is_refused() {
        let error = read(b">gnl|BL_ORD_ID|5 a\nACGT\n", false).unwrap_err();
        assert!(error.contains("database ordinal"), "{error}");
    }

    // `xSetSeqMol`: NCBI reads a record by the molecule its accession names. The letters of a
    // protein are dropped by the nucleotide reading (180 aa -> 121, probe `pdo`); a nucleotide
    // record under a protein accession reads the same.
    #[test]
    fn only_a_molecule_that_changes_the_letters_is_refused() {
        let error = read(b">NM_000001.1 a\nMKLVEF\n", true).unwrap_err();
        assert!(error.contains("names the other molecule"), "{error}");
        assert!(read(b">NP_000001.1 a\nACGTNNRY\n", false).is_ok());
        let error = read(b">NP_000001.1 a\nACGU\n", false).unwrap_err();
        assert!(error.contains("names the other molecule"), "{error}");
        assert!(read(b">NM_000001.1 a\nACGT\n", false).is_ok());
        assert!(read(b">NP_000001.1 a\nMKLVEF\n", true).is_ok());
    }

    // fasta.cpp:386-397: without fDLOptional a defline is required.
    #[test]
    fn a_record_without_a_defline_is_nothing_read() {
        let error = read(b"ACGT\n", false).unwrap_err();
        assert!(
            error.contains("Input doesn't start with a defline"),
            "{error}"
        );
    }

    #[test]
    fn the_archive_form_follows_the_id() {
        let local = read(b">lcl|q a\nACGT\n", false).unwrap();
        assert!(local[0].archive_seq_loc().is_none());
        let general = read(b">gnl|db|q a\nACGT\n", false).unwrap();
        assert!(general[0].archive_seq_loc().is_none());
        let mut text = read(b">NM_000001.1 a\nACGTACGT\n", false).unwrap();
        let loc = text[0].archive_seq_loc().unwrap();
        assert_eq!((loc.from, loc.to), (0, 7));
        text[0].cut_in_place(2, 6);
        let loc = text[0].archive_seq_loc().unwrap();
        assert_eq!((loc.from, loc.to), (2, 5));
    }

    /// The records of `fasta` read by the reader of `role` of `program` in the search `scope`.
    fn read_in(
        scope: &std::sync::Arc<IdScope>,
        role: &str,
        program: &'static str,
        protein: bool,
        fasta: &[u8],
    ) -> Result<Vec<super::super::FastaRecord>, String> {
        let config = if role == "query" {
            ReaderConfig::query(program, protein, false)
        } else {
            ReaderConfig::subject(program, protein, false)
        };
        let mut source = FastaInputSource::from_bytes(fasta, config.with_parse_deflines(true));
        source.set_id_scope(Some(scope.clone()));
        read_all(&mut source, &mut |_| Ok(())).map_err(|error| format!("{error:#}"))
    }

    // A query and a subject share no id, nor one accession with and without a version, unless they
    // are one record: defline and letters (DECISIONS 12, 15, 19).
    #[test]
    fn a_query_and_a_subject_share_an_id_only_as_one_record() {
        const QUERY: &[u8] =
            b">q1 t\nACGT\n>NM_000001 t\nACGT\n>ref|NM_000002.1| t\nACGT\n>NM_000003.1 t\nACGT\n\
              >gb|AAA00007.1|NAMEA t\nACGT\n>gb||NAMEB t\nACGT\n";
        for (subject, refused) in [
            (&b">s1 t\nACGT\n"[..], false),
            (b">NM_000002.2 t\nACGT\n", false),
            (b">q1 t\nACGT\n", false),
            (b">ref|NM_000002.1| t\nACGT\n", false),
            (b">lcl|q1 t\nACGT\n", true),
            (b">q1  t\nACGT\n", true),
            (b">q1 t\nAC\nGT\n", false),
            (b">q1 t\nacgu\n", false),
            (b">q1 other\nACGT\n", true),
            (b">q1\nACGT\n", true),
            (b">q1 t\nACGA\n", true),
            (b">NM_000001.1 t\nACGT\n", true),
            (b">NM_000002 t\nACGT\n", true),
            (b">NM_000002.1 t\nACGT\n", true),
            (b">NM_000003 t\nACGT\n", true),
            (b">gb|AAA00007.1|NAMEA t\nACGT\n", false),
            (b">gb||NAMEB t\nACGT\n", false),
            (b">gb|AAA00008.1|NAMEA t\nACGT\n", false),
            (b">gb||NAMEA t\nACGT\n", true),
            (b">gb|AAA00009.1|NAMEB t\nACGT\n", true),
        ] {
            let scope = std::sync::Arc::new(IdScope::default());
            read_in(&scope, "query", "BLASTN", false, QUERY).unwrap();
            let read = read_in(&scope, "subject", "BLASTN", false, subject);
            let text = String::from_utf8_lossy(subject).into_owned();
            assert_eq!(read.is_err(), refused, "{text}");
            if refused {
                assert!(
                    read.unwrap_err().contains("a query and a subject"),
                    "{text}"
                );
            }
        }
    }

    // The order of the readers does not change it, and the second reader's record knows.
    #[test]
    fn one_record_is_marked_with_the_lower_case_of_either_copy() {
        for (first, second, lowercase) in [
            (&b">q1 t\nACGT\n"[..], &b">q1 t\nACGT\n"[..], false),
            (b">q1 t\nacGT\n", b">q1 t\nACGT\n", true),
            (b">q1 t\nACGT\n", b">q1 t\nacGT\n", true),
        ] {
            let scope = std::sync::Arc::new(IdScope::default());
            let subjects = read_in(&scope, "subject", "BLASTN", false, first).unwrap();
            assert_eq!(ids_of(&subjects[0]).shared, None);
            let queries = read_in(&scope, "query", "BLASTN", false, second).unwrap();
            assert_eq!(ids_of(&queries[0]).shared, Some(SharedRecord { lowercase }));
        }
    }

    // BLASTX and TBLASTN: a nucleotide and a protein reader never read one record (DECISIONS 19).
    #[test]
    fn a_nucleotide_and_a_protein_reader_do_not_read_one_record() {
        let scope = std::sync::Arc::new(IdScope::default());
        read_in(&scope, "query", "TBLASTN", true, b">q1 t\nACGT\n").unwrap();
        let error = read_in(&scope, "subject", "TBLASTN", false, b">q1 t\nACGT\n").unwrap_err();
        assert!(error.contains("a query and a subject"), "{error}");
    }

    // DECISIONS 16: the ids belong to one search, not to the process.
    #[test]
    fn the_ids_of_a_search_do_not_reach_the_next_search() {
        let first = std::sync::Arc::new(IdScope::default());
        read_in(&first, "query", "BLASTN", false, b">kept1 t\nACGT\n").unwrap();
        read_in(&first, "subject", "BLASTN", false, b">kept1 u\nACGT\n").unwrap_err();
        let second = std::sync::Arc::new(IdScope::default());
        read_in(&second, "subject", "BLASTN", false, b">kept1 u\nACGT\n").unwrap();
        let mut keys = Vec::new();
        for id in read_in(&first, "subject", "BLASTN", false, b">s t\nACGT\n").unwrap()[0]
            .believed()
            .unwrap()
            .ids
            .iter()
        {
            scope_keys(id, &mut keys);
        }
        assert!(!first.query_has_key(&keys));
        keys.clear();
        scope_keys(
            &ParsedSeqId::Local(ObjectId::Str(b"kept1".to_vec())),
            &mut keys,
        );
        assert!(first.query_has_key(&keys));
        assert!(!second.query_has_key(&keys));
    }

    // Masks of a shared query: the masks of the query, or with `-lcase_masking` the lower case of
    // either copy (DECISIONS 19).
    #[test]
    fn outfmt_0_refuses_the_masks_of_a_shared_record_only() {
        let scope = std::sync::Arc::new(IdScope::default());
        read_in(&scope, "subject", "BLASTN", false, b">same1 t\nACgt\n").unwrap();
        let queries = read_in(
            &scope,
            "query",
            "BLASTN",
            false,
            b">other1 t\nACGT\n>same1 t\nACGT\n",
        )
        .unwrap();
        let refused = |masked: usize, lcase: bool| {
            refuse_shared_masks(&queries, |index| index == masked, lcase, "BLASTN").is_err()
        };
        assert!(!refused(0, false));
        assert!(refused(1, false));
        assert!(refused(0, true));
    }

    // NCBI BLAST+ 2.17.0 (probe `DI7/ws/pd_sub.py` `subjname`, `noid`): a `-subject` record with a
    // local id `Subject_N`, written or generated, is shown as one read without the option.
    #[test]
    fn subject_local_ids_are_shown_as_records() {
        for (defline, local) in [
            (&b">lcl|Subject_1 x"[..], true),
            (b">Subject_x x", true),
            (b">", true),
            (b">Subjectx x", false),
            (b">plain1 x", false),
            (b">NM_000003.1 x", false),
        ] {
            let mut fasta = defline.to_vec();
            fasta.extend_from_slice(b"\nACGT\n");
            let config = ReaderConfig::subject("BLASTN", false, false).with_parse_deflines(true);
            let mut source = FastaInputSource::from_bytes(&fasta, config);
            let records = read_all(&mut source, &mut |_| Ok(())).unwrap();
            let ids = records[0].believed().unwrap();
            assert_eq!(
                ids.subject_local(),
                local,
                "{}",
                String::from_utf8_lossy(defline)
            );
        }
    }
}
