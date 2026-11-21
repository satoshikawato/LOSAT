//! A database sequence with parsed Seq-ids in the reports: what each format prints of its
//! definition lines (`seqid::DeflineSet`), chosen and printed by the functions of `seqid`.
//!
//! NCBI's formatters find a subject again by its Seq-align id: the scope gives the Bioseq of
//! the definition line that holds the id (`DeflineSet::bioseq_defline`, the "target") with
//! the whole set as a user object. The ids and title of the target are the Bioseq's (tabular
//! `sseqid`, `sacc`, `saccver`, the description titles of outfmt 0, the XML `Hit_accession`);
//! the formats that show every definition line read the set (the outfmt 0 alignment heading,
//! the XML `Hit_def`, the XML2 and JSON `HitDescr` list, `stitle`).
//!
//! A database without parsed Seq-ids keeps its own path (the title and `gnl|BL_ORD_ID|n`,
//! `algorithm::common::database::stored_title`); `SeqDb::read_header` tells the two apart.
//!
//! LOSAT writes the names of the sequence it searched. NCBI's formatters look the Seq-align id
//! up again in the database, and where the lookup finds another sequence first (an accession
//! without a version after the same accession with one, `gb|X.1|` and `emb|X.1|`, `lcl|X` after
//! `gb|X.1|`: PLAN §9.6 of line DBID) they write that sequence's names and letters, or stop.
//! LOSAT looks the id up as NCBI does (`seqdb::lookup::SeqidLookup`, `found_instead`) and refuses
//! such a subject where a format reads its Bioseq (`fetched`). The two-letter PDB chains of one
//! molecule (`A`, `AA`, `B`) are not such a case (probe of the ID2a review).

use crate::seqdb::DbSequence;
use crate::seqid::defline_set::{
    alignment_line_id, compare, description_id, seq_id_list, tabular_id, Sic, TabularId,
};
use crate::seqid::model::PatentNumber;
use crate::seqid::rank::{find_best, text_score};
use crate::seqid::{label, BlastDefline, DeflineSet, ParsedSeqId};

use super::blast_xml::HitTaxonomy;
use super::defline::{
    bioseq_title_from_source, generate_bioseq_defline, BioseqFacts, Defline, SourceTaxname,
    UnknownEncoding,
};

/// A reported database sequence with parsed Seq-ids: its definition lines, its Seq-align id
/// and the definition line of its Bioseq.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParsedSubject {
    set: DeflineSet,
    alignment_id: ParsedSeqId,
    /// `DeflineSet::bioseq_defline` of `alignment_id`.
    target: usize,
    /// Whether NCBI's formatters find no Bioseq by the Seq-align id whatever its lookup gives
    /// (`new`, `fetched`).
    unfetchable: bool,
    /// What NCBI's lookup of the Seq-align id in the database finds where that is not this
    /// sequence (`found_instead`, `fetched`).
    found_instead: Option<FoundInstead>,
    /// Another sequence that NCBI's scope can show for this one by a weak match where its
    /// formatters read the batch's sequences one by one (`weak_match`, `fetched`).
    weak_match: Option<WeakMatch>,
    /// Another sequence (`DbSequence`) that holds an id NCBI's scope takes for the Seq-align id,
    /// besides the sequence that the lookup gives first or where it gives none
    /// (`FoundInstead::Sequence`, `Nothing`): NCBI's scope can read it for this one where it has
    /// read it before (`also_held`, `pssm_read`).
    also_held: Option<DbSequence>,
    /// `-show_gis`: the outfmt 0 description line and heading show it (`description`, `heading`).
    show_gis: bool,
    /// The name of the organism of the subject Bioseq's first source, from which NCBI makes its
    /// title (`with_source`).
    source: SourceTaxname,
    /// Whether this is a `-subject` record (`of_record`): a Bioseq without definition lines.
    record: bool,
    /// Whether it is a record that the reports show as a record without `-parse_deflines`
    /// (`BelievedIds::subject_local`): only its tabular ids are its own (`shown_as_record`).
    local_record: bool,
    /// The strings of the Seq-aligns' `use_this_seqid` (`under_positive_list`); empty without a
    /// positive list.
    use_this_seqid: Vec<Vec<u8>>,
    /// `tabular_saccver`, kept for the reports of later query batches.
    saccver: std::sync::OnceLock<Option<std::sync::Arc<[u8]>>>,
    /// `query_scope_keys`, made at the first report that asks.
    query_scope_keys: std::sync::OnceLock<QueryScopeKeys>,
}

/// The keys by which NCBI's scope finds a query for the Seq-align id of a reported database
/// sequence (`ParsedSubject::query_scope_keys`; `common::database::refuse_id_of_a_query`), and
/// for its display id, by which the description table looks its title up again
/// (`common::database::refuse_display_id_of_a_query`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QueryScopeKeys {
    /// Its own keys (`defline_set::scope_keys`).
    pub own: Vec<u64>,
    /// The keys of the ids that it matches weakly (`id_scope::weak_keys`).
    pub weak: Vec<u64>,
    /// The keys of the display id (`ParsedSubject::display_id`) and of the ids that it matches
    /// weakly, as `own` and `weak` of the Seq-align id: the same lookup of NCBI's scope.
    pub display: Vec<u64>,
}

/// The `use_this_seqid` strings of a database subject, `subject` when it has parsed Seq-ids: a
/// sequence without them has a `general` id of `BL_ORD_ID`, which is neither a text Seq-id nor a
/// GI (`ParsedSubject::under_positive_list`).
pub fn use_this_seqid(subject: Option<&ParsedSubject>) -> &[Vec<u8>] {
    subject.map_or(&[], ParsedSubject::use_this_seqid)
}

/// `AdjustSeqIdType` with `eGetId_ForceAcc` (`ParsedSubject::sam_id`) of a Bioseq with the ids
/// `ids` and the Seq-align id `alignment_id`.
pub fn forced_accession<'a>(
    alignment_id: &'a ParsedSeqId,
    ids: &'a [ParsedSeqId],
) -> &'a ParsedSeqId {
    let acc_ver =
        |id: &crate::seqid::model::TextseqId| id.accession.is_some() && id.version.is_some();
    if alignment_id.textseq().is_some_and(acc_ver) {
        return alignment_id;
    }
    let mut found = None;
    for id in ids {
        if let Some(text) = id.textseq() {
            found = Some(id);
            if acc_ver(text) {
                break;
            }
        }
    }
    found.unwrap_or(alignment_id)
}

/// What NCBI's lookup of a subject's Seq-align id in its database finds where that is not the
/// subject (`ParsedSubject::found_instead`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FoundInstead {
    /// Another sequence of the database (`DbSequence`) first: NCBI's formatters show it where
    /// its ids compare equal to the Seq-align id, else no sequence.
    Sequence(DbSequence),
    /// The sequence itself first (or the same record), and later another sequence of the
    /// database (`DbSequence`) that holds an id NCBI's scope takes for the Seq-align id: where
    /// both are reported, NCBI's formatters can show that one for it.
    AlsoHeld(DbSequence),
    /// Another sequence of the database (`DbSequence`) with the same definition lines but for
    /// their taxids and the same letters (the same record in two databases made with other
    /// taxids), that the lookup gives first or that also holds the Seq-align id: NCBI's
    /// formatters can show it for the subject, which only the reports that write the subject's
    /// taxids show apart (`fetched`).
    OtherTaxids(DbSequence),
    /// Such a sequence, also with the same taxids in another volume, that every report that
    /// reads the Bioseq can show apart, for the reason given.
    CopyApart(DbSequence, CopyApart),
    /// No sequence.
    Nothing,
    /// The lookup fails on a file of the database (`LookupFailure`): NCBI's formatters, and its
    /// prefetch, stop with the `CSeqDBException`.
    Failed(LookupFailure),
    /// LOSAT's own list of the other sequences that hold the Seq-align id fails on a file that
    /// NCBI's lookup does not read (`why`: what LOSAT cannot list, with the error), so whether
    /// another sequence holds an id that NCBI's scope takes for it is unknown: the list of every
    /// OID of the lookup (`SeqidLookup::seqid_to_oids` with `multi`) on an ISAM of a volume after
    /// the first one that holds the id (a version 4 string ISAM that does not open), where NCBI's
    /// lookup has returned (review 3 of line DBID loadord, H3), or the searches of the keys that
    /// differ from the Seq-align id's in the case of its letters (`SeqidLookup::case_folded_oids`
    /// and `case_variant_oids`: version 5). `instead`: the sequence that the lookup gives where it
    /// is another one (as `Sequence`), else the lookup gives the subject or the same record.
    HoldersUnknown {
        instead: Option<DbSequence>,
        why: String,
    },
    /// The lookup throws on the id itself before it reads a file (`SeqidLookup::throws_on_id`: a
    /// `gnl|BL_ORD_ID|0` that makeblastdb stores with a text tag), with this message: NCBI's
    /// formatters stop, its prefetch catches the exception with its warning
    /// (`ParsedSubject::prefetch_fails`).
    ThrowsOnId(String),
}

/// What NCBI's scope can show for a subject by a weak match where its formatters read the
/// reported sequences of a batch one by one (`PRE_FETCH_SEQS_LIMIT` turns the batch's prefetch
/// off; the prefetch loads each by its own id): `ParsedSubject::weak_match`, refused by
/// `fetched` with `one_by_one`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum WeakMatch {
    /// Another sequence of the database (`DbSequence`) that holds an id the Seq-align id
    /// matches weakly (`seqid::id_scope::weak_keys`): where that one is read first, NCBI's scope
    /// shows it for the subject.
    Held(DbSequence),
    /// A Seq-align id with weak keys that a version 5 database cannot look up: one with a name,
    /// or a `pir` or `prf` one (`common::database::HeaderReader::weakly_held`).
    Unlisted,
    /// The lookup of the sequences that could match weakly failed, with this message (LOSAT's
    /// lookup of the name alone; NCBI's scope makes none).
    Failed(String),
}

/// What NCBI's prefetch of a query batch reads for a subject that no report shows, by its
/// Seq-align id (`ParsedSubject::prefetched`): psiblast's traceback hit lists beyond
/// `-max_target_seqs`, and the Seq-annots with `PRE_FETCH_SEQS_LIMIT=2147483647`
/// (`common::database::prefetch_reads`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Prefetched {
    /// The sequence itself, or the same record whose lookups of its taxids write the same texts
    /// (`ParsedSubject::fetched` without the taxonomy formats).
    Itself,
    /// No sequence: the prefetch reads nothing for it.
    Nothing,
    /// Another sequence (`DbSequence`) where LOSAT knows it, whose taxids the prefetch can look up
    /// instead: `FoundInstead::Sequence` and `HoldersUnknown` with another sequence that the
    /// lookup gives (that sequence where one of its ids compares equal to the Seq-align id, else
    /// none), `AlsoHeld` (`or_itself`: the subject or that sequence, as the ids of the blobs read
    /// before in the run decide, `CBlastDbDataLoader::m_Ids`), and the same record whose lookups
    /// write other texts (`None`); `or_none`: or nothing, where the Bioseq of another sequence
    /// that holds the id, which the scope read before, can stand for it without a load (one that
    /// `also_held` names, or one that LOSAT cannot list, `HoldersUnknown`); `reason` as `fetched`
    /// gives it.
    Another {
        other: Option<DbSequence>,
        or_itself: bool,
        or_none: bool,
        reason: String,
    },
    /// The sequence itself, or nothing where the Bioseq of another sequence that holds the id,
    /// which the scope read before and which LOSAT cannot list, stands for it
    /// (`FoundInstead::HoldersUnknown` where the lookup gives the subject or the same record);
    /// `reason` as `fetched` gives it.
    ItselfOrNothing { reason: String },
    /// A lookup that fails (the reason): on the id itself NCBI's prefetch catches the exception
    /// and writes the warning `Error pre-fetching sequence data`
    /// (`ParsedSubject::prefetch_fails`, `common::database::batch_prefetch`); on a file it
    /// rethrows the `CSeqDBException` and stops.
    Failed(String),
}

/// A lookup of a Seq-align id that fails on a file of the database (`FoundInstead::Failed`): the
/// message, and whether it is NCBI's own (`SeqDbError::ncbi_message`: the text of the
/// `CSeqDBException` that NCBI writes after `BLAST Database error: `, as "DB contains no accession
/// info." of a list that puts version 5 databases with and without the accession table side by
/// side), else LOSAT's text of an error NCBI does not have.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LookupFailure {
    message: String,
    ncbi: NcbiFailure,
}

/// What NCBI's lookup does where LOSAT's fails on a file ([`LookupFailure`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum NcbiFailure {
    /// It throws a `CSeqDBException` with the message (NCBI's own text).
    Message,
    /// It throws, with a message that LOSAT does not reproduce.
    Other,
    /// It crashes on the file (`SeqDbError::ncbi_crashes`).
    Crashes,
}

impl LookupFailure {
    /// What NCBI does where it reads a subject whose lookup fails on a file: its prefetch rethrows
    /// the `CSeqDBException` (`BLAST Database error`, exit status 2); where `PRE_FETCH_SEQS_LIMIT`
    /// turns the prefetch off, its formatters catch it and write `Unknown`, then stop or not as
    /// the format decides (DECISIONS 53 (1) of line DBID).
    pub const NCBI_READS: &'static str = "NCBI BLAST+ 2.17.0's prefetch stops with a BLAST Database error, and where the prefetch is off its formatters fail to read it";

    /// What NCBI does where it reads a subject whose lookup throws NCBI's exception at a place
    /// LOSAT knows ([`LookupFailure::exception`]), and why LOSAT refuses it in this run: LOSAT
    /// writes each format's failure as NCBI does (`FormatThrow::of`) but for the formats whose
    /// writers it does not make so, and a call of the library that mixes formats that stop at
    /// other places after it (`FormatThrow::and`): NCBI's reference for such a call is one run of
    /// each format, which stop at different places with different exit statuses and standard
    /// errors, where the call is one search with one result (review 2 of line DBID pfoff-3, M2).
    pub const NCBI_THROWS: &'static str = "NCBI BLAST+ 2.17.0's prefetch stops with a BLAST Database error, and where the prefetch is off its formatters fail to read it and stop or go on as each format does; LOSAT writes that for outfmt 0 to 4, 6, 7, 10, 13 to 16, 18 and 20, but not for outfmt 5 or 17, nor in a call of LOSAT's library that mixes outfmt 0, outfmt 1 to 4 and the other formats, which NCBI writes in runs of their own that stop at different places, where the call is one search with one result";

    /// What NCBI does where its lookup reads a file that it crashes on (an ISAM whose header does
    /// not fit its files, as an empty `.psi` or `.psd`: SIGSEGV; review 4 of line DBID loadord):
    /// a run whose lookups never reach it (the sequence was loaded before by another id) goes on.
    pub const NCBI_CRASHES: &'static str =
        "NCBI BLAST+ 2.17.0 crashes where its lookup reads that file";

    /// NCBI's own message.
    pub fn ncbi(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            ncbi: NcbiFailure::Message,
        }
    }

    /// LOSAT's text of an error that NCBI does not have.
    pub fn losat(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
            ncbi: NcbiFailure::Other,
        }
    }

    /// The failure of a lookup with the error `error`: NCBI's message where it is one of NCBI's
    /// (`SeqDbError::ncbi_message`), NCBI's crash where LOSAT refuses a file that NCBI crashes on
    /// (`SeqDbError::ncbi_crashes`), else LOSAT's text.
    pub fn of_error(error: &crate::seqdb::SeqDbError) -> Self {
        match error.ncbi_message() {
            Some(message) => Self::ncbi(message),
            None if error.ncbi_crashes() => Self {
                message: error.to_string(),
                ncbi: NcbiFailure::Crashes,
            },
            None => Self::losat(error.to_string()),
        }
    }

    /// The text NCBI writes after `BLAST Database error: `, where the message is NCBI's own.
    pub fn ncbi_message(&self) -> Option<&str> {
        (self.ncbi == NcbiFailure::Message).then_some(self.message.as_str())
    }

    /// The `CSeqDBException` that NCBI's lookup throws, where the message is NCBI's own
    /// (`ncbi_message`) and is one that NCBI's lookup of a Seq-align id was seen to throw, with
    /// its place (`seqdb::lookup_throw_site`: "DB contains no accession info." of
    /// `CBlastLMDBManager::CBlastEnv::GetDbi`, at a list that puts version 5 databases with and
    /// without the accession table side by side); `None` for every other text (also the other
    /// texts of NCBI's database reader whose places LOSAT knows, `seqdb::throw_site`, which no
    /// lookup was seen to throw), where only the prefetch's report of it (`BLAST Database error:
    /// <message>`) is known.
    pub fn exception(&self) -> Option<crate::ncbi_exception::NcbiException> {
        let message = self.ncbi_message()?;
        let site = crate::seqdb::lookup_throw_site(message)?;
        Some(crate::ncbi_exception::NcbiException::new(site, message))
    }

    /// What NCBI does where it reads the subject ([`LookupFailure::NCBI_READS`],
    /// [`LookupFailure::NCBI_THROWS`] where its exception is known,
    /// [`LookupFailure::NCBI_CRASHES`]), for the reason of LOSAT's refusal.
    pub fn what_ncbi_does(&self) -> &'static str {
        match self.ncbi {
            NcbiFailure::Message if self.exception().is_some() => Self::NCBI_THROWS,
            NcbiFailure::Message | NcbiFailure::Other => Self::NCBI_READS,
            NcbiFailure::Crashes => Self::NCBI_CRASHES,
        }
    }
}

impl std::fmt::Display for LookupFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

/// Why the reports can show apart a copy of a subject with the same definition lines but for
/// their taxids (`FoundInstead::CopyApart`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CopyApart {
    /// A title can be made from the taxids (`ParsedSubject::title_reads_source`), which only
    /// the reports that show a generated title show apart
    /// (`ParsedSubject::generated_title_fetched`).
    Title,
    /// The taxonomy database writes a text at a lookup that the volume of the copy makes for
    /// its taxids.
    TaxonomyTexts,
}

/// One `HitDescr` of the XML2 and JSON formats (`s_SeqAlignToXMLHit`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HitDescription {
    pub id: Vec<u8>,
    /// Not written when empty.
    pub accession: Vec<u8>,
    /// Not written when the definition line has no title member.
    pub title: Option<Vec<u8>>,
    /// Not written when 0.
    pub taxid: i32,
    /// What the taxonomy database said about `taxid` (`Unnamed` where nothing was asked); a
    /// `Named` scientific name is written after the taxid.
    pub taxonomy: HitTaxonomy,
}

impl ParsedSubject {
    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:390-395
    // ```c++
    //     int oid = -1;
    //
    //     if (! m_BlastDb->SeqidToOid(*seqid, oid)) {
    //         _TRACE("FAILED to find '" << seqid->AsFastaString() << "'");
    //         return -1;
    //     }
    // ```
    // NCBI's formatters find the subject's sequence again from its Seq-align id, by the lookup that
    // `found_instead` ports, whatever kind of id it is: in version 5 `bbs|123`, `nat|X.1|` and
    // `lcl|777` by their LMDB keys, in version 4 the ISAM files that makeblastdb writes for neither
    // `bbs|` nor `nat|`, a `gnl|PIG|` id in the PIG files that it fills from no FASTA id, and a
    // `gnl|BL_ORD_ID|n` as the stored n plus the first OID of its volume
    // (`SeqDb::defline_set`), looked up from the start of each volume, so as its own sequence where
    // the stored n is its own number in its volume, also beside another id (probe of line DBID N2b,
    // `blastdb-ids-20261021/R/unfetch/p1`: 317 cases of NCBI 2.17.0, databases of NCBI's and
    // LOSAT's makeblastdb). The lookup does not decide one case: where the definition line of the
    // Bioseq is found by the content of a local id alone (`lcl|777` before the definition line of
    // the Seq-align id `gi|777`, `DeflineSet::bioseq_defline`), NCBI 2.17.0 finds no Bioseq, in
    // version 4 and 5 (outfmt 6 stops with a CCoreException, XML2 writes the Seq-align id alone;
    // probe of the ID2a review). LOSAT refuses that subject instead of repeating a defect of
    // NCBI's.
    /// The subject of the definition lines `set` of a database. `Err` names what LOSAT does not
    /// write as NCBI: a stored title or Seq-id with a byte outside ' ' ..= '~' (NCBI writes an
    /// error message each time it reads the definition line) and a definition line without a
    /// Seq-id (NCBI's formatters take the best of none; makeblastdb writes none such).
    pub fn new(set: DeflineSet) -> Result<Self, String> {
        if set
            .deflines
            .iter()
            .any(|defline| defline.replaced_title_bytes > 0 || defline.replaced_id_bytes > 0)
        {
            return Err("has a stored title or Seq-id with a byte outside ' ' to '~' (NCBI BLAST+ shows it as '#' and writes an error message each time it reads the definition line)".into());
        }
        if set.deflines.iter().any(|defline| defline.ids.is_empty()) {
            return Err("has a definition line without a Seq-id".into());
        }
        let Some(alignment_id) = set.alignment_id().cloned() else {
            return Err("has no Seq-id".into());
        };
        let Some(target) = set.bioseq_defline(&alignment_id) else {
            return Err("has no definition line that holds its Seq-align id".into());
        };
        // A giim Seq-align id is never found again: NCBI's handles never take it for one with an
        // equal id (`algorithm::psiblast::rounds::HandleIndex`), so its scope finds no Bioseq by it
        // (NCBI BLAST+ 2.17.0 writes `Unknown`; review 2 of line DBID pssmrows, H1).
        let unfetchable = matches!(alignment_id, ParsedSeqId::Giim(_))
            || set
                .display_defline(&alignment_id)
                .is_some_and(|shown| shown.index != target);
        Ok(Self {
            set,
            alignment_id,
            target,
            unfetchable,
            found_instead: None,
            weak_match: None,
            also_held: None,
            show_gis: false,
            source: SourceTaxname::Unnamed,
            record: false,
            local_record: false,
            use_this_seqid: Vec::new(),
            saccver: std::sync::OnceLock::new(),
            query_scope_keys: std::sync::OnceLock::new(),
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:458-466
    // ```c++
    //     // set ID
    //     retval->SetInt().SetId().Assign(*FindBestChoice(itr->GetId(), CSeq_id::BestRank));
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:272-283
    // ```c++
    //     CRef<CBlast_def_line_set> bdlRef = CSeqDB::ExtractBlastDefline(handle);
    //     const list< CRef< CBlast_def_line > > &bdl = (bdlRef.Empty()) ? list< CRef< CBlast_def_line > >() : bdlRef->Get();
    //     ...
    //     if(bdl.empty()){ //no blast defline struct, should be no such case now
    //         //actually not so fast...as we now fetch from entrez even when it's not in blast db
    //         //there is no blast defline in such case.
    //         seqid = GetSeqIdListString(ids, show_gi);
    //         defline = sequence::CDeflineGenerator().GenerateDefline(handle);
    // ```
    /// The `-subject` record with the ids `ids` and the FASTA title `title` under
    /// `-parse_deflines`, as one definition line: its Bioseq has no Blast-def-line user object, so
    /// its titles are generated and `stitle` and `salltitles` are `N/A`, and NCBI's scope holds
    /// it, so it is always found.
    pub fn of_record(
        ids: &crate::blastinput::fasta_reader::believed::BelievedIds,
        title: &[u8],
    ) -> Self {
        Self {
            set: DeflineSet {
                deflines: vec![BlastDefline {
                    title: (!title.is_empty()).then(|| title.to_vec()),
                    ids: ids.ids.clone(),
                    ..BlastDefline::default()
                }],
            },
            alignment_id: ids.best.clone(),
            target: 0,
            unfetchable: false,
            found_instead: None,
            weak_match: None,
            also_held: None,
            show_gis: false,
            source: SourceTaxname::Unnamed,
            record: true,
            local_record: ids.subject_local(),
            use_this_seqid: Vec::new(),
            saccver: std::sync::OnceLock::new(),
            query_scope_keys: std::sync::OnceLock::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:484-498
    // ```c++
    //     if (sid_in->IsLocal()) {
    //         string id_token;
    //         vector<string> title_tokens;
    //         title_tokens = NStr::Split(CAlignFormatUtil::GetTitle(bh), " ", title_tokens);
    //         ...
    //         if (id_token == NcbiEmptyString || parse_local) {
    //             const CObject_id& obj_id = sid_in->GetLocal();
    // ```
    /// Whether the outfmt 0, XML and archive reports show this `-subject` record as one read
    /// without `-parse_deflines` (a local id `Subject_N`, `BelievedIds::subject_local`); its
    /// `sall*` are still its own ids (`s_ReplaceLocalId` keeps a local id with `parse_local`).
    pub fn shown_as_record(&self) -> bool {
        self.local_record
    }

    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:393-416
    // ```c++
    //     if (! m_BlastDb->SeqidToOid(*seqid, oid)) {
    //         _TRACE("FAILED to find '" << seqid->AsFastaString() << "'");
    //         return -1;
    //     }
    //     _TRACE("Found '" << seqid->AsFastaString() << "' at OID " << oid);
    //
    //     // Test for deflines.  If the filtering eliminates the Seq-id we
    //     // are interested in, we just pretend we don't know anything about
    //     // this Seq-id.  If there are other data loaders installed, they
    //     // will have an opportunity to resolve the Seq-id.
    //
    //     bool found = false;
    //
    //     IBlastDbAdapter::TSeqIdList filtered = m_BlastDb->GetSeqIDs(oid);
    //
    //     ITERATE(IBlastDbAdapter::TSeqIdList, id, filtered) {
    //         if (seqid->Compare(**id) == CSeq_id::e_YES) {
    //             found = true;
    //             break;
    //         }
    //     }
    // ```
    // NCBI's data loader serves the Bioseq of the first OID that `SeqidToOid` gives for the
    // Seq-align id (`seqdb::lookup::SeqidLookup::seqid_to_oid`), and only where one of that
    // sequence's ids compares equal to it: a version 5 database finds `gb|X|` under the key `X`,
    // which makeblastdb also writes for `gb|X.1|`, and `gb|X.1|` compares equal to `gb|X|` and
    // to `emb|X.1|`. NCBI 2.17.0 then writes the other sequence's names and counts the
    // alignment's identities against its letters (outfmt 6 `AAA00001.1 100.000 ... 0` for the
    // subject `gb|AAA00001|` of 99.561 % identity), or, for `lcl|X` after `gb|X.1|`, finds no
    // Bioseq and stops with a CCoreException in outfmt 6 (probe of line DBID N1,
    // `blastdb-ids-20261021/R/accoid/p1`; version 4 alike).
    /// The subject as NCBI's formatters find it again: `found` is what NCBI's lookup of the
    /// Seq-align id in the database finds where that is not this sequence (another one or
    /// none), which the formats that read the subject's Bioseq refuse (`fetched`).
    pub fn found_instead(mut self, found: Option<FoundInstead>) -> Self {
        self.found_instead = found;
        self
    }

    /// The subject with `weak`, what NCBI's scope can show for it by a weak match where its
    /// formatters read the batch's sequences one by one (`fetched`'s `one_by_one`).
    pub fn weak_match(mut self, weak: Option<WeakMatch>) -> Self {
        self.weak_match = weak;
        self
    }

    /// The subject with `holder`, another sequence that holds an id NCBI's scope takes for the
    /// Seq-align id where the lookup gives a third sequence first or none (`pssm_read`).
    pub fn also_held(mut self, holder: Option<DbSequence>) -> Self {
        self.also_held = holder;
        self
    }

    /// Whether NCBI's formatters find no Bioseq by the Seq-align id whatever its lookup gives
    /// (see `new`), so that no lookup is needed. Read by tests only: the engine asks
    /// `common::database::scope_bioseq`, which looks it up first.
    #[cfg(test)]
    pub fn unfetchable(&self) -> bool {
        self.unfetchable
    }

    /// What NCBI's lookup of the Seq-align id finds where that is not this sequence
    /// (`found_instead`).
    pub fn found(&self) -> Option<&FoundInstead> {
        self.found_instead.as_ref()
    }

    /// The subject as the outfmt 0 report of a run with or without `-show_gis` shows it.
    pub fn showing_gis(mut self, show_gis: bool) -> Self {
        self.show_gis = show_gis;
        self
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1146-1168
    // ```c++
    //         bool found_taxid_in_taxonomy_blastdb = true;
    //         ...
    //         try {
    //             found_taxid_in_taxonomy_blastdb = tax_info.GetTaxNames(*taxid, namerecord);
    //         }
    //         catch (CSeqDBException &) {
    //             found_taxid_in_taxonomy_blastdb = false;
    //         }
    //         ...
    //         if (found_taxid_in_taxonomy_blastdb) {
    //             org->SetTaxname().swap(namerecord.scientific_name);
    // ```
    // A title made from the source (`x_SetTitleFromBioSrc`, of the first source of the Bioseq,
    // whose volume puts the target's taxids first) names the organism only where the taxonomy
    // database gives a name for the taxid; else NCBI writes the title of a Bioseq without one
    // (review 1 of line DBID N1b L1, `blastdb-ids-20261021/R/titlesrc/review1/r2.tsv`).
    // A RefSeq chromosome, mRNA or ncRNA record gets its title from the source too: NCBI's own
    // functions for them make nothing of a database Bioseq (`report::defline::source_title`).
    // The end of a protein title ending in `]` is rewritten from the same name
    // (`x_AdjustProteinTitleSuffix`, `adjusted_title_shown`; probe of line RG T5c).
    /// The subject with the name of its Bioseq's first source (`names` of its taxid,
    /// `first_source_taxid`, asked only where `title_from_taxonomy` or, of a `protein` database,
    /// `adjusted_title_shown`), from which its generated title is made (`source`,
    /// `bioseq_title`). Other subjects ask nothing.
    pub fn with_source(mut self, protein: bool, names: impl FnOnce(i32) -> SourceTaxname) -> Self {
        if self.title_from_taxonomy() || self.adjusted_title_shown(protein) {
            if let Some(taxid) = self.first_source_taxid() {
                self.source = names(taxid);
            }
        }
        self
    }

    /// The taxid of the subject Bioseq's first source where NCBI fetches it by its Seq-align id
    /// (`DeflineSet::first_source_taxid`, the first of `volume_lookups`).
    pub fn first_source_taxid(&self) -> Option<i32> {
        self.set.first_source_taxid(Some(&self.alignment_id))
    }

    /// The name of the organism of the subject Bioseq's first source (`with_source`).
    pub fn source(&self) -> &SourceTaxname {
        &self.source
    }

    /// The subject as a database with a positive list (`positive_list`,
    /// `SeqDb::has_positive_list`) reports it: its Seq-aligns name the ids of the definition
    /// lines the list keeps. `use_gi` is whether the Seq-align id is a GI.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqinfosrc_aux.cpp:113-150
    /// ```c++
    /// void GetFilteredRedundantSeqids(const IBlastSeqInfoSrc & seqinfo_src,
    ///                              int                      oid,
    ///                              vector<string>            & seqids,
    /// 			     bool			use_gi)
    /// {
    ///     seqids.resize(0);
    ///
    ///     if (! seqinfo_src.HasGiList()) {
    ///         return;
    ///     }
    ///
    ///     list< CRef<CSeq_id> > seqid_list = seqinfo_src.GetId(oid);
    ///
    ///     ITERATE(list< CRef<CSeq_id> >, id, seqid_list) {
    /// 	if (use_gi)
    /// 	{
    ///         	if ((**id).IsGi())
    /// 		{
    /// 			string sid_string = "gi:" + (*id)->GetSeqIdString(true);
    /// 			seqids.push_back(sid_string);
    /// 		}
    /// 	}
    /// 	else
    /// 	{
    /// 		const CTextseq_id* tsid = (*id)->GetTextseq_Id();
    /// 		if (tsid && tsid->IsSetAccession())
    /// 		{
    /// 			string sid_string = "seqid:" + (*id)->GetSeqIdString(true);
    /// 			seqids.push_back(sid_string);
    /// 		}
    /// 		else if ((*id)->IsPdb())
    /// 		{
    /// 			string sid_string = "seqid:" + (*id)->GetSeqIdString(true);
    /// 			seqids.push_back(sid_string);
    /// 		}
    /// 	}
    ///     }
    /// }
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1604-1605
    /// ```c++
    ///         // Get SeqIds for entrez query restriction.
    /// 	vector<string> seqid_list;
    ///         GetFilteredRedundantSeqids(*seqinfo_src, hsp_list->oid, seqid_list, subject_id->IsGi());
    /// ```
    pub fn under_positive_list(mut self, positive_list: bool) -> Self {
        if !positive_list {
            return self;
        }
        let use_gi = self.alignment_id.gi().is_some();
        for id in self.set.ids() {
            let prefix: &[u8] = if use_gi {
                if id.gi().is_none() {
                    continue;
                }
                b"gi:"
            } else {
                if !(id.textseq().is_some_and(|text| text.accession.is_some())
                    || matches!(id, ParsedSeqId::Pdb(_)))
                {
                    continue;
                }
                b"seqid:"
            };
            let mut string = prefix.to_vec();
            crate::seqid::label::content(id, true, &mut string);
            self.use_this_seqid.push(string);
        }
        self
    }

    /// The strings of the Seq-aligns' `use_this_seqid` (`under_positive_list`).
    pub fn use_this_seqid(&self) -> &[Vec<u8>] {
        &self.use_this_seqid
    }

    /// `Err` for a subject whose sequence NCBI's formatters cannot find by its Seq-align id
    /// (see `new` and `found_instead`): the formats that read the subject's Bioseq (outfmt 0, 5,
    /// 6, 7, 10, 13 to 16 and SAM, 17) refuse it; the Seq-annots and the archive (8, 9, 11, 12)
    /// write the Seq-align id alone. `taxonomy`: whether a report writes the subject's taxids
    /// or their names (outfmt 13 to 16, the tabular `staxid`, `staxids` and the names), which
    /// alone refuse a sequence found as the same record with other taxids
    /// (`FoundInstead::OtherTaxids`; also where a generated title can be made from them, which
    /// the reports that show it refuse, `generated_title_fetched`). `one_by_one`: whether NCBI's
    /// formatters may read the subject one by one (a batch that reports it has its prefetch
    /// off, `common::database::ParsedShown::one_by_one`), which alone refuses a sequence that
    /// another can stand for by a weak match (`WeakMatch`), after what the lookup finds.
    /// `load`: what the reports do at the subject's `Hit` (`HitLoads::load`, the one answer): a
    /// subject whose `Hit` throws (`HitLoad::Throws`: its lookup fails on a file with NCBI's
    /// exception in a batch whose prefetch is off, and every format of the run writes it as NCBI
    /// does, `FormatThrow::of_run`: the XML2 and JSON documents of outfmt 13 to 16, the tables of
    /// outfmt 6, 7, 10 and 20, the organism report of outfmt 18) is not refused for what the lookup
    /// finds (it finds nothing: it throws before any definition line is read), and one that no
    /// report loads (`HitLoad::Unread`: shown only after such a `Hit`) is not refused at all;
    /// every other subject (`HitLoad::Loaded`) is refused as above, a lookup that fails on a file
    /// among them (outfmt 5 here, before any byte is written).
    pub fn fetched(&self, taxonomy: bool, one_by_one: bool, load: HitLoad) -> Result<(), String> {
        match load {
            HitLoad::Unread => return Ok(()),
            HitLoad::Throws => {}
            HitLoad::Loaded => self.found_fetched(taxonomy)?,
        }
        self.weakly_fetched(one_by_one)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:357-360, 412-417
    // ```c++
    //     const CSeq_id& kSeqId = kAlignSet.Get().front()->GetSeq_id(1);
    //
    //     try {
    //         const CBioseq_Handle& subj_handle = scope->GetBioseqHandle(kSeqId);
    //     ...
    //     } catch (const CException&) {
    //     	CRef<blastxml2::CHitDescr> hit_exp(new blastxml2::CHitDescr);
    //     	hit_exp->SetId(kSeqId.AsFastaString());
    //         hit->SetDescription().push_back(hit_exp);
    //         hit->SetLen(sequence::GetLength(kSeqId, scope));
    //     };
    // ```
    // Where the batch's prefetch is off, the `Hit` loads its subject by the Seq-align id: the
    // lookup throws on the file there, and again in the `catch` (`GetLength` loads by the same id),
    // outside any `catch` of the hit; the report function of outfmt 15 and 16 catches it, 13 and 14
    // end the program with it (probe of line DBID mixdb, `PU.off`; DECISIONS 62 (2)).
    /// The exception that NCBI's lookup of the Seq-align id throws where a formatter loads the
    /// subject one by one (a batch whose prefetch is off): a lookup that fails on a file with
    /// NCBI's `CSeqDBException` whose place LOSAT knows (`FoundInstead::Failed`,
    /// `LookupFailure::exception`). Where the batch prefetches, its prefetch throws first and
    /// stops the run (`prefetch_stops`).
    pub fn lookup_throws(&self) -> Option<crate::ncbi_exception::NcbiException> {
        match &self.found_instead {
            Some(FoundInstead::Failed(failure)) => failure.exception(),
            _ => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:440-450 (quoted at
    // `algorithm::psiblast::pssm::input_msa`)
    /// LOSAT's refusal of what psiblast's PSSM reads by the Seq-align id of a subject that it
    /// includes (the `CSeqVector` of each Seq-align's range in NCBI's scope), where LOSAT does not
    /// read it as NCBI: what the scope gets is `common::database::scope_bioseq`'s, and this
    /// refuses another sequence that also holds an id the scope takes for the Seq-align id, which
    /// the scope reads for it where it has read it before (`FoundInstead::AlsoHeld` where the
    /// lookup gives the subject first, `also_held` where it gives another or none), a lookup that
    /// fails on a file with LOSAT's error (`FoundInstead::Failed`; with NCBI's exception, which the
    /// read catches, `prefetch_stops`: the caller reads no letters), the same record whose read writes other
    /// taxonomy texts (`CopyApart::TaxonomyTexts`), and with `one_by_one` a weak match
    /// (`fetched`). The reasons are what the lookup finds, also for a subject whose definition
    /// line the formatters do not find (`unfetchable`), which the PSSM reads by the scope as any
    /// other (review 2 of line DBID pssmread, M2).
    pub fn pssm_read(&self, one_by_one: bool) -> Result<(), String> {
        let id = || String::from_utf8_lossy(&self.alignment_fasta()).into_owned();
        match &self.found_instead {
            Some(
                found @ (FoundInstead::AlsoHeld(_)
                | FoundInstead::HoldersUnknown { .. }
                | FoundInstead::Failed(_)
                | FoundInstead::CopyApart(_, CopyApart::TaxonomyTexts)),
            ) => return Err(self.found_reason(found, false).unwrap_or_default()),
            Some(FoundInstead::Sequence(_) | FoundInstead::Nothing) => {
                if let Some(holder) = self.also_held {
                    return Err(format!(
                        "has the Seq-align id {}, which database sequence {holder} also holds (an id that NCBI BLAST+'s scope takes for it) besides what the lookup gives: NCBI BLAST+ 2.17.0 reads that sequence for it where it has read it before, as its reads and the blobs' layout decide",
                        id()
                    ));
                }
            }
            _ => {}
        }
        self.weakly_fetched(one_by_one)
    }

    /// Another sequence (`DbSequence`) that holds an id NCBI's scope takes for the Seq-align id
    /// besides what the lookup gives (`also_held`, `FoundInstead::AlsoHeld`), whose Bioseq the
    /// scope gets for this one where it holds it from an earlier read
    /// (`common::database::scope_bioseq`).
    pub fn other_holder(&self) -> Option<DbSequence> {
        match &self.found_instead {
            Some(FoundInstead::AlsoHeld(other)) => Some(*other),
            Some(FoundInstead::Sequence(_) | FoundInstead::Nothing) => self.also_held,
            _ => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1437-1439
    // ```c++
    //     seqids = defline->GetSeqid();
    //     ...
    //     bioseq->SetId().swap(seqids);
    // ```
    /// Whether the Bioseq that `GetBioseqNoData` builds of this sequence for its Seq-align id
    /// holds that id: the ids of the definition line its filter keeps (`target`,
    /// `DeflineSet::bioseq_defline`) hold one that compares equal to it (the scope's test on the
    /// Bioseq's ids, `common::database::scope_bioseq`). Not where that line is found by the
    /// content of a local id alone (`lcl|777` before `gi|777`), nor for a giim id.
    pub fn bioseq_holds_alignment_id(&self) -> bool {
        !matches!(self.alignment_id, ParsedSeqId::Giim(_))
            && self
                .target()
                .ids
                .iter()
                .any(|id| compare(id, &self.alignment_id) == Sic::Yes)
    }

    /// `fetched`'s refusal of a weak match where NCBI's formatters may read the subject one by
    /// one (`one_by_one`).
    fn weakly_fetched(&self, one_by_one: bool) -> Result<(), String> {
        let id = || String::from_utf8_lossy(&self.alignment_fasta()).into_owned();
        match self.weak_match.as_ref().filter(|_| one_by_one) {
            None => Ok(()),
            Some(WeakMatch::Held(other)) => Err(format!(
                "has the Seq-align id {}, which can match an id of database sequence {other} weakly (an accession with a version for one without, a name with or without an accession, a pir or prf accession as a name): in a batch whose prefetch PRE_FETCH_SEQS_LIMIT turns off, where NCBI BLAST+ 2.17.0 reads the reported sequences one by one and that sequence first, it can write that sequence's names and letters for it",
                id()
            )),
            Some(WeakMatch::Unlisted) => Err(format!(
                "has the Seq-align id {}, which ids of other sequences can match weakly in any letter case (a name with or without an accession, a pir or prf accession as a name) and which a version 5 database does not look up: in a batch whose prefetch PRE_FETCH_SEQS_LIMIT turns off, where NCBI BLAST+ 2.17.0 reads the reported sequences one by one, it can write the names and letters of such a sequence read before",
                id()
            )),
            Some(WeakMatch::Failed(message)) => Err(format!(
                "has the Seq-align id {}, for which LOSAT cannot look up the sequences that can match it weakly ({message}): in a batch whose prefetch PRE_FETCH_SEQS_LIMIT turns off, where NCBI BLAST+ 2.17.0 reads the reported sequences one by one, it can write the names and letters of such a sequence read before",
                id()
            )),
        }
    }

    /// `fetched` but for the weak matches: what the lookup of the Seq-align id finds
    /// (`found_instead`), whatever the batch's prefetch.
    fn found_fetched(&self, taxonomy: bool) -> Result<(), String> {
        if self.unfetchable {
            return Err(format!(
                "has the Seq-align id {}, by which NCBI BLAST+'s formatters do not find the sequence (NCBI BLAST+ 2.17.0 writes Unknown or another sequence's names, or stops)",
                String::from_utf8_lossy(&self.alignment_fasta())
            ));
        }
        match &self.found_instead {
            None => Ok(()),
            Some(found) => self.found_reason(found, taxonomy).map_or(Ok(()), Err),
        }
    }

    /// The reason of a refusal of what the lookup of the Seq-align id finds, `found`, for a format
    /// that reads the Bioseq (`fetched`; `taxonomy` as there), or for psiblast's PSSM
    /// (`pssm_read`); `None` for a copy of the record that such a report shows as the subject.
    fn found_reason(&self, found: &FoundInstead, taxonomy: bool) -> Option<String> {
        let id = String::from_utf8_lossy(&self.alignment_fasta()).into_owned();
        Some(match found {
            FoundInstead::OtherTaxids(_) | FoundInstead::CopyApart(_, CopyApart::Title)
                if !taxonomy =>
            {
                return None
            }
            FoundInstead::OtherTaxids(other) | FoundInstead::CopyApart(other, CopyApart::Title) => format!(
                "has the Seq-align id {id}, which database sequence {other} also holds with the same definition lines and letters but other taxids (NCBI BLAST+ 2.17.0 can write that sequence's taxids and their names for it)"
            ),
            FoundInstead::CopyApart(other, CopyApart::TaxonomyTexts) => format!(
                "has the Seq-align id {id}, which database sequence {other} also holds with the same definition lines and letters: where NCBI BLAST+ 2.17.0 shows that sequence for it, the volume of that sequence looks its taxids up, and the taxonomy database writes a text at such a lookup"
            ),
            FoundInstead::Sequence(other) => format!(
                "has the Seq-align id {id}, by which NCBI BLAST+'s lookup gives database sequence {other} first (NCBI BLAST+ 2.17.0 writes that sequence's names and letters, or stops)"
            ),
            FoundInstead::AlsoHeld(other) => format!(
                "has the Seq-align id {id}, which database sequence {other} also holds (where both are reported, NCBI BLAST+ 2.17.0 can write that sequence's names and letters for it)"
            ),
            FoundInstead::HoldersUnknown { instead, why } => format!(
                "has the Seq-align id {id}, {}for which LOSAT cannot list {why}: where such a sequence is also reported, NCBI BLAST+ 2.17.0 can write its names and letters for it",
                instead.map_or_else(String::new, |other| format!(
                    "by which NCBI BLAST+'s lookup gives database sequence {other} first (NCBI BLAST+ 2.17.0 writes that sequence's names and letters, or stops), and "
                ))
            ),
            FoundInstead::Failed(failure) => format!(
                "has the Seq-align id {id}, whose lookup in the database fails ({failure}; {})",
                failure.what_ncbi_does()
            ),
            FoundInstead::ThrowsOnId(message) => format!(
                "has the Seq-align id {id}, whose lookup in the database fails ({message}; NCBI BLAST+ 2.17.0 stops)"
            ),
            FoundInstead::Nothing => format!(
                "has the Seq-align id {id}, by which NCBI BLAST+'s formatters do not find the sequence (NCBI BLAST+ 2.17.0 writes Unknown or another sequence's names, or stops)"
            ),
        })
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:770-786
    // ```c++
    //        s_ExtractSeqidsAndRanges(results, ids, ranges);
    //        _TRACE("Prefetching " << ids.size() << " sequence lengths");
    //        LoadSequencesToScope(ids, ranges, scope);
    // 	}
    //     catch(const CSeqDBException & e) {
    //         _TRACE(e.GetMsg());
    //     	throw;
    //     }
    //     catch (const CException& e) {
    // 	   if(s_IsUsingRemoteBlastDbDataLoader()) {
    //           ERR_POST(Warning << "Error fetching sequence data from BLAST databases at NCBI, "
    //                      "please try again later");
    // 	   }
    // 	   else {
    //           ERR_POST(Warning << "Error pre-fetching sequence data");
    // 	   }
    //    }
    // ```
    // The prefetch loads each subject's Bioseq by its Seq-align id, as the formatters find it
    // (`fetched`), by the exact lookup only (the weak matches come where it is off): an id that
    // the lookup finds nowhere loads nothing, without a message (probe of DBID P1 review 1:
    // `gnl|PIG|124`, `gnl|BL_ORD_ID|40` of no sequence); a lookup that fails (a `gnl|BL_ORD_ID|0`
    // stored with a text tag) ends the prefetch with that warning; `lcl|555` that an earlier
    // sequence holds loads that sequence, whose taxids it looks up.
    // The loader reads by the Seq-align id what its lookup gives, also for a subject whose
    // definition line the formatters do not find (`unfetchable`): its own blob where the lookup
    // gives it, whose taxids it looks up (review 2 of DBID P1: `lcl|777` before `gi|777`,
    // `bad_magic`), another sequence's where it gives that one first (review 2 of line DBID
    // pssmread, H3: `gi|777|gb|YYY00001.1|` of taxid 99990 before `lcl|777 ^A gi|777`).
    /// What NCBI's prefetch reads for this subject where no report shows it (`Prefetched`).
    pub fn prefetched(&self) -> Prefetched {
        let reason = || {
            self.found_instead
                .as_ref()
                .and_then(|found| self.found_reason(found, false))
                .unwrap_or_default()
        };
        match &self.found_instead {
            None
            | Some(FoundInstead::OtherTaxids(_) | FoundInstead::CopyApart(_, CopyApart::Title)) => {
                Prefetched::Itself
            }
            Some(FoundInstead::Nothing) => Prefetched::Nothing,
            // The prefetch catches the exception of a failed lookup and writes its warning, or
            // rethrows a database error (`CSeqDBException`) and stops; the formats that read the
            // Bioseq stop (`fetched`).
            Some(FoundInstead::ThrowsOnId(message)) => Prefetched::Failed(format!(
                "has the Seq-align id {}, whose lookup in the database fails ({message}; NCBI BLAST+ 2.17.0's prefetch writes Error pre-fetching sequence data, and a report that reads it stops)",
                String::from_utf8_lossy(&self.alignment_fasta())
            )),
            Some(FoundInstead::Failed(failure)) => Prefetched::Failed(format!(
                "has the Seq-align id {}, whose lookup in the database fails ({failure}; {})",
                String::from_utf8_lossy(&self.alignment_fasta()),
                failure.what_ncbi_does()
            )),
            Some(FoundInstead::Sequence(other)) => Prefetched::Another {
                other: Some(*other),
                or_itself: false,
                or_none: self.also_held.is_some(),
                reason: reason(),
            },
            Some(FoundInstead::HoldersUnknown {
                instead: Some(other),
                ..
            }) => Prefetched::Another {
                other: Some(*other),
                or_itself: false,
                or_none: true,
                reason: reason(),
            },
            Some(FoundInstead::HoldersUnknown { instead: None, .. }) => {
                Prefetched::ItselfOrNothing { reason: reason() }
            }
            Some(FoundInstead::AlsoHeld(other)) => Prefetched::Another {
                other: Some(*other),
                or_itself: true,
                or_none: false,
                reason: reason(),
            },
            Some(FoundInstead::CopyApart(_, CopyApart::TaxonomyTexts)) => Prefetched::Another {
                other: None,
                or_itself: false,
                or_none: false,
                reason: reason(),
            },
        }
    }

    /// Whether the lookup of this subject's Seq-align id throws on the id where NCBI's prefetch
    /// loads it (`FoundInstead::ThrowsOnId`): the batch's prefetch catches it and writes `Error
    /// pre-fetching sequence data` (`common::database::batch_prefetch`). A lookup that fails on a
    /// file stops NCBI (`FoundInstead::Failed`).
    pub fn prefetch_fails(&self) -> bool {
        matches!(self.found_instead, Some(FoundInstead::ThrowsOnId(_)))
    }

    /// NCBI's message where the lookup of this subject's Seq-align id fails on a file of the
    /// database with NCBI's own `CSeqDBException` (`FoundInstead::Failed`,
    /// `LookupFailure::ncbi_message`; "DB contains no accession info." of a list that puts
    /// version 5 databases with and without the accession table side by side): NCBI's prefetch of
    /// a batch that reports the subject rethrows it, which stops the run with `BLAST Database
    /// error: <message>` (exit status 2) before the batch's reports
    /// (`common::database::batch_prefetch`).
    pub fn prefetch_stops(&self) -> Option<&str> {
        match &self.found_instead {
            Some(FoundInstead::Failed(failure)) => failure.ncbi_message(),
            _ => None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2250-2255
    // ```c++
    // 		if(bdl->IsSetTitle()){
    // 			alnDispParams->title = bdl->GetTitle();
    // 		}
    //         if(alnDispParams->title.empty()) {
    //             alnDispParams->title = CDeflineGenerator().GenerateDefline(bsp_handle);
    //         }
    // ```
    // NCBI's formatters generate the title of the Bioseq they fetch (`CDeflineGenerator`, from
    // the sources of that sequence's volume) for the outfmt 0 description (showdefline.cpp:498)
    // and heading, for the first word that replaces an id in the tabular formats
    // (showdefline.cpp:227) and for the XML2 and JSON `title` (blastxml2_format.cpp:405); the
    // tabular `stitle` and `salltitles` (tabular.cpp:415-434) and the XML `Hit_def`
    // (`GetBioseqHandleDeflineAndId`, showdefline.cpp:283-299) write the titles of the
    // definition lines.
    /// `Err` for a subject found as the same record with other taxids from which NCBI's
    /// formatters can make the title they generate for it (`FoundInstead::CopyApart` with
    /// `CopyApart::Title`), for the reports that show a generated title: outfmt 0
    /// (`check_shown_parsed_subject`) and the tabular word of a replaced id
    /// (`check_printed_word`); outfmt 13 to 16 write the taxids (`fetched`). The callers ask
    /// it where the title they show can depend on the source (`pairwise_title_reads_source`, the
    /// word of a nucleotide or of a protein title ending in `]`).
    pub fn generated_title_fetched(&self) -> Result<(), String> {
        match &self.found_instead {
            Some(FoundInstead::CopyApart(other, CopyApart::Title)) => Err(format!(
                "has the Seq-align id {}, which database sequence {other} also holds with the same definition lines and letters but other taxids, from which NCBI BLAST+ 2.17.0 can make that sequence's title for it (a title made from the taxids with the taxonomy database, or the end of a protein title rewritten)",
                String::from_utf8_lossy(&self.alignment_fasta())
            )),
            _ => Ok(()),
        }
    }

    /// The Seq-align id (outfmt 8, 9, 11, 12; `DeflineSet::alignment_id`).
    pub fn alignment_id(&self) -> &ParsedSeqId {
        &self.alignment_id
    }

    // NCBI reference (598d8ae6): c++/src/objtools/format/sam_formatter.cpp:131-141
    // ```c++
    // void CSAM_CIGAR_Formatter::AdjustSeqIdType(CConstRef<CSeq_id>& id)
    // {
    //     ...
    //     CSeq_id_Handle forced_id = sequence::GetId(*id, *scope, force_type);
    //     if (forced_id) {
    //         id.Reset(forced_id.GetSeqId());
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:2743-2757 (`GetAccVer`, which
    // `GetId(eGetId_ForceAcc)` calls: the id itself where it is an accession with a version,
    // else `x_GetAccVer` of the Bioseq's ids)
    // ```c++
    //     if ( !(flags & CScope::fForceLoad) ) {
    //         if ( idh.IsAccVer() ) {
    //             return idh;
    //         }
    //     }
    //     ...
    //             CSeq_id_Handle ret = CScope::x_GetAccVer(info->GetIds());
    // ```
    // NCBI reference (598d8ae6): c++/src/objmgr/scope.cpp:425-445
    // ```c++
    //     ITERATE ( TIds, iter, ids ) {
    //         if ( iter->IsGi() ) {
    //             // optimization to avoid calling GetSeqId() (it may be expensive)
    //             continue;
    //         }
    //         if ( auto seq_id = iter->GetSeqId() ) { // need to hold CConstRef<>
    //             if ( auto text_id = seq_id->GetTextseq_Id() ) {
    //                 // found some text id, it can be returned
    //                 ret = *iter;
    //                 if ( text_id->IsSetAccession() && text_id->IsSetVersion() ) {
    //                     // found acc.ver, stop searching
    //                     break;
    //                 }
    // ```
    /// The Seq-id that NCBI's SAM formatter (outfmt 17) names the subject by
    /// (`AdjustSeqIdType` with `eGetId_ForceAcc`): the Seq-align id where it is an accession
    /// with a version, else the first such id of the subject Bioseq (its gis passed over), else
    /// its last Textseq-id, else the Seq-align id.
    pub fn sam_id(&self) -> &ParsedSeqId {
        forced_accession(&self.alignment_id, &self.target().ids)
    }

    /// The label of the subject's rows in a query-anchored report (outfmt 1 to 4): the Bioseq's
    /// Seq-ids as `x_FillSeqid` reads them, with `-show_gis` (`report::query_anchored::fill_seqid`).
    pub(crate) fn anchored_label(&self, out: &mut Vec<u8>) {
        crate::report::query_anchored::fill_seqid(&self.target().ids, self.show_gis, out);
    }

    /// The definition lines.
    pub fn deflines(&self) -> &[BlastDefline] {
        &self.set.deflines
    }

    /// The definition line of the subject Bioseq.
    fn target(&self) -> &BlastDefline {
        &self.set.deflines[self.target]
    }

    /// What `GenerateDefline` reads of the subject Bioseq: its Seq-ids, whether it has a source
    /// and the name of its first source (none where LOSAT refuses it, `source`).
    fn bioseq_facts(&self) -> BioseqFacts {
        let mut facts = bioseq_facts(&self.target().ids);
        facts.has_source = self.has_source();
        facts.taxname = self.source.name().ok().flatten();
        facts
    }

    /// `GenerateDefline` of the subject Bioseq: for the alignment heading
    /// (`leave_prefix_suffix` false) or the description table (true).
    pub fn bioseq_title(
        &self,
        protein: bool,
        leave_prefix_suffix: bool,
    ) -> Result<Defline, UnknownEncoding> {
        let title = self.target().title.as_deref().unwrap_or_default();
        generate_bioseq_defline(title, protein, leave_prefix_suffix, &self.bioseq_facts())
    }

    /// Whether the subject Bioseq has a source: some definition line has a positive taxid
    /// (`CSeqDBVol::x_GetTaxonomy` makes a source of each, seqdbvol.cpp:1136).
    fn has_source(&self) -> bool {
        self.first_source_taxid().is_some()
    }

    /// The lookups of the volume's taxonomy cache when NCBI fetches the subject Bioseq by its
    /// Seq-align id: the positive taxids of the definition lines in the order of
    /// `CSeqDBVol::x_GetTaxonomy` (`x_GetTaxDefline` puts the line that holds the id first,
    /// `DeflineSet::taxonomy_lookups`), a taxid of several lines once for each.
    pub fn volume_lookups(&self) -> Vec<i32> {
        self.set.taxonomy_lookups(Some(&self.alignment_id))
    }

    /// The positive taxids that the volume looks up when it fetches the subject Bioseq, each once,
    /// in the order of their first lookup ([`ParsedSubject::volume_lookups`]).
    pub fn volume_taxids(&self) -> Vec<i32> {
        let mut taxids = Vec::new();
        for taxid in self.volume_lookups() {
            if !taxids.contains(&taxid) {
                taxids.push(taxid);
            }
        }
        taxids
    }

    /// Whether NCBI makes the subject Bioseq's title from its source (`bioseq_title_from_source`)
    /// and that source has a taxid: the title then depends on the taxonomy database.
    pub fn title_from_taxonomy(&self) -> bool {
        let title = self.target().title.as_deref().unwrap_or_default();
        bioseq_title_from_source(title, || self.bioseq_facts()) && self.has_source()
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:3964-3972
    // ```c++
    //     // adjust protein partial/organelle/taxname suffix, if necessary
    //     if ( m_IsAA && (! m_MainTitle.empty()) && (! (flags & fLeavePrefixSuffix)) ) {
    //         if (m_Idx) {
    //             x_AdjustProteinTitleSuffixIdx (bsh);
    //         } else {
    //             x_AdjustProteinTitleSuffix (bsh);
    //         }
    //     }
    // ```
    // `x_AdjustProteinTitleSuffix` (`report::defline::adjust_protein_title_suffix`) rewrites the
    // end of a protein title ending in `]` from the Bioseq's first source, also without a name
    // (a ` [ ]` end goes): `alpha, partial [Bacillus sp.]` is written `alpha [Bacillus sp.]`
    // (probe of line RG T5c, `T/t5c/out/ncbi/ppar.blastp.f0.tx`).
    /// Whether the subject Bioseq's title can lose its end to `x_AdjustProteinTitleSuffix`
    /// (`report::defline::protein_title_suffix_adjustable` of the target's title).
    fn title_ends_in_bracket(&self) -> bool {
        let title = self.target().title.as_deref().unwrap_or_default();
        crate::report::defline::protein_title_suffix_adjustable(title)
    }

    /// Whether some definition line has an id that `seq_id_list` replaces by the first word of
    /// the generated title (`replaces_an_id`; `all_tabular_ids`, `lookup_failure`).
    fn has_replaced_id(&self) -> bool {
        self.set
            .deflines
            .iter()
            .any(|defline| replaces_an_id(&defline.ids))
    }

    /// Whether a report shows the subject Bioseq's generated title (flags 0) where
    /// `x_AdjustProteinTitleSuffix` can rewrite its end from the source: a protein with a source
    /// whose title ends in `]`, with a line without a title (the outfmt 0 heading, outfmt 18) or
    /// an id that the title's first word replaces (the tabular ids). `with_source` names the
    /// source then.
    fn adjusted_title_shown(&self, protein: bool) -> bool {
        protein
            && self.title_ends_in_bracket()
            && self.has_source()
            && (self.heading_generates_title() || self.has_replaced_id())
    }

    /// Whether a title that the reports write of the subject can depend on its Bioseq's
    /// source, the taxids of its definition lines, whatever they are (`title_from_taxonomy`,
    /// `adjusted_title_shown`): in outfmt 0 (`pairwise_title_reads_source`), and the tabular word
    /// of an id that `seq_id_list` replaces, of a protein title ending in `]` (a title made from
    /// the source with the taxonomy database is also shown in outfmt 0; of a protein the word is
    /// `unnamed` either way).
    pub fn title_reads_source(&self, protein: bool, taxdb_found: bool) -> bool {
        self.pairwise_title_reads_source(protein, taxdb_found)
            || (protein && self.title_ends_in_bracket() && self.has_replaced_id())
    }

    /// Whether a title that outfmt 0 shows of the subject can depend on its Bioseq's source
    /// (`title_reads_source`): a title made from the source where `taxdb.bti` is found
    /// (`taxdb_found`; without it NCBI 2.17.0 writes the title of a Bioseq without a source,
    /// probe of line DBID N1b, `blastdb-ids-20261021/R/copytaxid/p1`) in the description and the
    /// heading, and a protein title ending in `]` in the heading of a set with a definition line
    /// without a title, which shows the generated title without `fLeavePrefixSuffix` (the
    /// description keeps the stored end: probe of line DBID N1b L1,
    /// `blastdb-ids-20261021/R/titlesrc/p1`).
    pub fn pairwise_title_reads_source(&self, protein: bool, taxdb_found: bool) -> bool {
        let title = self.target().title.as_deref().unwrap_or_default();
        (taxdb_found && bioseq_title_from_source(title, || self.bioseq_facts()))
            || (protein && self.title_ends_in_bracket() && self.heading_generates_title())
    }

    /// Whether the outfmt 0 alignment heading shows the subject Bioseq's generated title: some
    /// definition line has no title (`heading`).
    pub fn heading_generates_title(&self) -> bool {
        self.set
            .deflines
            .iter()
            .any(|defline| defline.title.as_deref().is_none_or(<[u8]>::is_empty))
    }

    /// The first word of the subject Bioseq's generated title, which `seq_id_list` puts in the
    /// place of a `gnl|BL_ORD_ID` or `lcl|Subject_` id of the Bioseq; empty when no such id is
    /// there (the title is then not made). `Err` where NCBI's `x_CleanAndCompress` reads past the
    /// end of the title (NCBI crashes) and the other refusals of `generated_title_word`.
    /// `Ok(None)` where `GenerateDefline` throws for that title (a title in no encoding that
    /// `CUtf8::GuessEncoding` recognises): `GetSeqIdList` makes no ids (`lookup_failure`).
    fn title_word(&self, protein: bool) -> Result<Option<Vec<u8>>, &'static str> {
        if !replaces_an_id(&self.target().ids) {
            return Ok(Some(Vec::new()));
        }
        self.generated_title_word(protein)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:222-230
    // ```c++
    //         if (((*itr)->IsGeneral() &&
    //             (*itr)->AsFastaString().find("gnl|BL_ORD_ID")
    //             != string::npos) ||
    // 		(*itr)->AsFastaString().find("lcl|Subject_") != string::npos) {
    //             vector<string> title_tokens;
    //             string defline = sequence::CDeflineGenerator().GenerateDefline(bh);
    //             if (defline != NcbiEmptyString) {
    //                 id_token =
    //                     NStr::Split(defline, " ", title_tokens)[0];
    // ```
    // The first word of a title that NCBI makes from the source with the taxonomy database is
    // the organism's for a nucleotide (`Homo` for an empty title with taxid 9606) and `unnamed`
    // for a protein, as without a source; of a protein title ending in `]`, the word before the
    // end that `x_AdjustProteinTitleSuffix` rewrites where the Bioseq has a source (probe of line
    // DBID N1b L1, `blastdb-ids-20261021/R/titlesrc/p1`). NCBI prints of the ids that
    // `GetSeqIdList` makes the GI and the best one (`GetSeqIdListString`, showdefline.cpp:155-187),
    // so a replaced id beside a better one is not printed.
    /// `Err` where a tabular report prints (`printed`) an id that the first word of a title
    /// NCBI can make from the source replaces: a title made from the taxids of the same record
    /// found for the subject (`generated_title_fetched`) for a nucleotide or a protein title
    /// ending in `]`. The subject's own title made from its source is made with its name
    /// (`source`); the other refusals of the word are `tabular_ids`' and `all_tabular_ids`'.
    pub fn check_printed_word(
        &self,
        protein: bool,
        printed: super::tabular::PrintedIds,
    ) -> Result<(), &'static str> {
        let problem = if self.generated_title_fetched().is_err()
            && (!protein || self.title_ends_in_bracket())
        {
            "a Seq-align id that another database sequence also holds with the same definition lines and letters but other taxids, from which NCBI BLAST+ 2.17.0 can make the title whose first word replaces a printed id (a title made from the taxids with the taxonomy database, or the end of a protein title rewritten)"
        } else {
            return Ok(());
        };
        // Whether the printed ids of `ids` change with the word.
        let shows_word = |ids: &[ParsedSeqId], kinds: &[TabularId]| {
            let (one, other) = (seq_id_list(ids, b"x"), seq_id_list(ids, b"y"));
            kinds.iter().any(|&kind| {
                let (mut a, mut b) = (Vec::new(), Vec::new());
                tabular_id(&one, kind, &mut a);
                tabular_id(&other, kind, &mut b);
                a != b
            })
        };
        let subject = printed.subject
            && shows_word(
                &self.target().ids,
                &[
                    TabularId::Full,
                    TabularId::Accession,
                    TabularId::AccessionVersion,
                ],
            );
        let all =
            printed.all
                && !self.record
                && self.set.deflines.iter().any(|defline| {
                    shows_word(&defline.ids, &[TabularId::Full, TabularId::Accession])
                });
        if subject || all {
            return Err(problem);
        }
        Ok(())
    }

    /// The first word of the subject Bioseq's generated title (`title_word` without the check
    /// of the target's ids): `Err` also where LOSAT refuses the source of a title made from it
    /// (`SourceTaxname::Refused`); `Ok(None)` where `GenerateDefline` throws.
    fn generated_title_word(&self, protein: bool) -> Result<Option<Vec<u8>>, &'static str> {
        self.source.name()?;
        let Ok(defline) = self.bioseq_title(protein, false) else {
            return Ok(None);
        };
        // A nucleotide title stops at the end of the string (approved exception 2 of
        // PD-LOSAT-NCBI-DEFECTS for BLASTN, TBLASTN and TBLASTX).
        if protein && defline.reads_past_end {
            return Err(concat!(
                "a title ",
                crate::report::defline::punctuation_title!()
            ));
        }
        let end = defline
            .text
            .iter()
            .position(|&byte| byte == b' ')
            .unwrap_or(defline.text.len());
        Ok(Some(defline.text[..end].to_vec()))
    }

    /// The tabular `sseqid`, `sacc` and `saccver` (`tabular_id` of the subject Bioseq's ids after
    /// `seq_id_list`); `Ok(None)` where `SetSubjectId`'s `GetSeqIdList` throws (`lookup_failure`).
    pub fn tabular_ids(&self, protein: bool) -> Result<Option<[Vec<u8>; 3]>, &'static str> {
        let Some(word) = self.title_word(protein)? else {
            return Ok(None);
        };
        let ids = seq_id_list(&self.target().ids, &word);
        let field = |kind| {
            let mut out = Vec::new();
            tabular_id(&ids, kind, &mut out);
            out
        };
        Ok(Some([
            field(TabularId::Full),
            field(TabularId::Accession),
            field(TabularId::AccessionVersion),
        ]))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:862-869
    // ```c++
    //         } catch (const CException&) {
    //             list<CRef<CSeq_id> > subject_ids;
    //             CRef<CSeq_id> id(new CSeq_id());
    //             id->Assign(align.GetSeq_id(1));
    //             subject_ids.push_back(id);
    //             SetSubjectId(subject_ids);
    // ```
    /// Which of NCBI's lookups of the subject's tabular ids throw (`report::tabular::FailedLookup`):
    /// `SetSubjectId` where the subject Bioseq's ids hold one that `GetSeqIdList` replaces
    /// (`lcl|Subject_`, `gnl|BL_ORD_ID`) and `GenerateDefline` throws, `x_SetSubjectIds` where a
    /// definition line holds one (a `-subject` record's ids are not looked up there); `None`
    /// where neither throws. `Err` as `tabular_ids`'.
    pub fn lookup_failure(
        &self,
        protein: bool,
    ) -> Result<Option<super::tabular::LookupFailure>, &'static str> {
        let subject_id = replaces_an_id(&self.target().ids);
        let all_ids = !self.record && self.has_replaced_id();
        if !subject_id && !all_ids {
            return Ok(None);
        }
        Ok(self
            .generated_title_word(protein)?
            .is_none()
            .then_some(super::tabular::LookupFailure {
                subject_id,
                all_ids,
                bioseq: false,
            }))
    }

    /// The Seq-align id as `sallseqid` and `sallacc` print it (`tabular_id` of the id alone,
    /// `eFullId` and `eAccession`): what the catch of `SetFields` pushes where a lookup throws
    /// (`report::tabular::FailedLookup::full`, `accession`).
    pub fn alignment_tabular_ids(&self) -> [Vec<u8>; 2] {
        let ids = std::slice::from_ref(&self.alignment_id);
        let field = |kind| {
            let mut out = Vec::new();
            tabular_id(ids, kind, &mut out);
            out
        };
        [field(TabularId::Full), field(TabularId::Accession)]
    }

    /// The Seq-align id as `sallgi` prints it (`FindGi` of the id alone): its GI, 0 for an id that
    /// is none.
    pub fn alignment_gi(&self) -> i64 {
        self.alignment_id.gi().unwrap_or(0)
    }

    /// The keys of the ids of queries that NCBI's scope would show for this sequence
    /// (`QueryScopeKeys`), made once.
    pub fn query_scope_keys(&self) -> &QueryScopeKeys {
        self.query_scope_keys.get_or_init(|| {
            let mut own = Vec::new();
            scope_keys(self, &mut own);
            let mut weak = Vec::new();
            crate::seqid::id_scope::weak_keys(&self.alignment_id, &mut weak);
            let mut display = Vec::new();
            if let Some(id) = self.display_id() {
                crate::seqid::defline_set::scope_keys(id, &mut display);
                crate::seqid::id_scope::weak_keys(id, &mut display);
            }
            QueryScopeKeys { own, weak, display }
        })
    }

    /// The `saccver` of `tabular_ids`, made once for a subject (`protein` is its database's);
    /// `None` where the lookup throws, as `tabular_ids` (the programs' id lists,
    /// `common::tabular::SubjectId`).
    pub fn tabular_saccver(
        &self,
        protein: bool,
    ) -> Result<Option<std::sync::Arc<[u8]>>, &'static str> {
        if let Some(saccver) = self.saccver.get() {
            return Ok(saccver.clone());
        }
        let saccver = self
            .tabular_ids(protein)?
            .map(|[_, _, saccver]| std::sync::Arc::from(saccver));
        Ok(self.saccver.get_or_init(|| saccver).clone())
    }

    /// `Err` as `tabular_ids`', which it does not make.
    pub fn check_tabular_ids(&self, protein: bool) -> Result<(), &'static str> {
        self.title_word(protein).map(drop)
    }

    /// The tabular `sgi`: the first GI of the subject Bioseq's ids (`FindGi`; `seq_id_list` only
    /// replaces ids that are no GI), 0 for none.
    pub fn tabular_gi(&self) -> i64 {
        self.target().gi()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:544-556
    // ```c++
    //     if (bdlRef.NotEmpty() && bdlRef->CanGet() && bdlRef->IsSet() && !bdlRef->Get().empty()){
    //         vector< CConstRef<CSeq_id> > original_seqids;
    //
    //         ITERATE(CBlast_def_line_set::Tdata, itr, bdlRef->Get()) {
    //             original_seqids.clear();
    //             ITERATE(CBlast_def_line::TSeqid, id, (*itr)->GetSeqid()) {
    //                 original_seqids.push_back(*id);
    //             }
    //             list<CRef<objects::CSeq_id> > next_seqid_list;
    //             // Next call replaces BL_ORD_ID if found.
    //             CShowBlastDefline::GetSeqIdList(bh,original_seqids,next_seqid_list);
    //             m_SubjectIds.push_back(next_seqid_list);
    //         }
    // ```
    /// The tabular ids of every definition line (`x_SetSubjectIds`: `sallseqid`, `sallgi`,
    /// `sallacc`), in the stored order: each line's ids after `seq_id_list` with the first word
    /// of the subject Bioseq's title, as `tabular_ids` makes them. Where `GenerateDefline` throws,
    /// the ids of the lines before the first line with an id to replace, which `x_SetSubjectIds`
    /// pushed before it threw (`AllTabularIds::failed`). `Err` as `tabular_ids`'.
    pub fn all_tabular_ids(&self, protein: bool) -> Result<AllTabularIds, &'static str> {
        // `GetSeqIdList` replaces the ids of each definition line, so the word is made when any
        // line has an id to replace, not only the target's.
        let replaced = !self.record && self.has_replaced_id();
        let word = if replaced {
            self.generated_title_word(protein)?
        } else {
            Some(Vec::new())
        };
        let Some(word) = word else {
            let before = self
                .set
                .deflines
                .iter()
                .take_while(|defline| !replaces_an_id(&defline.ids))
                .map(|defline| defline_tabular_ids(defline, b""))
                .collect();
            return Ok(AllTabularIds {
                ids: before,
                failed: true,
            });
        };
        Ok(AllTabularIds {
            ids: self
                .set
                .deflines
                .iter()
                .map(|defline| defline_tabular_ids(defline, &word))
                .collect(),
            failed: false,
        })
    }

    /// The title of the first definition line, empty for none (`stitle`; the writers show an
    /// empty one as `N/A`).
    pub fn first_title(&self) -> &[u8] {
        self.set
            .deflines
            .first()
            .and_then(|defline| defline.title.as_deref())
            .unwrap_or_default()
    }

    /// The start of the outfmt 0 description line (`description_id` of `GetDisplayIds`' id) and
    /// its title: `GenerateDefline` of the subject Bioseq with `fLeavePrefixSuffix` (`Unknown`
    /// where NCBI's `HtmlDecode` throws), after two spaces where `gi|N|` was written.
    pub fn description(&self, protein: bool) -> (Vec<u8>, Defline) {
        let (start, wrote_gi) = self.description_id();
        let mut title = self.bioseq_title(protein, true).unwrap_or(Defline {
            text: b"Unknown".to_vec(),
            reads_past_end: false,
        });
        if wrote_gi {
            title.text.splice(0..0, *b"  ");
        }
        (start, title)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1576-1604
    // ```c++
    //     sdl->id = id;
    //     sdl->defline = "Unknown";
    //     ...
    //     try{
    //         const CBioseq_Handle& handle = m_ScopeRef->GetBioseqHandle(*id);
    //         x_FillDeflineAndId(handle, *id, use_this_seqid, sdl, blast_rank);
    //     } catch (const CException&){
    //         sdl->defline = "Unknown";
    //         ...
    //         if((*id).Which() == CSeq_id::e_Gi){
    //             sdl->gi = (*id).GetGi();
    //         } else {
    //             sdl->gi = ZERO_GI;
    //         }
    //         sdl->id = id;
    // ```
    /// The outfmt 0 description line (`description`) of a subject whose lookup throws there
    /// (`HitLoad::Throws`): the `catch` keeps the Seq-align id and the title `Unknown`, and the
    /// GI of a GI id; the line's start is made of them as of any other (`description_id`: the bare
    /// id and a space, or with `-show_gis` the FASTA form of an id that is not a GI, which no space
    /// follows: `ref|NP_100001.1|Unknown`). A GI Seq-align id goes the same way (`description_id`,
    /// which the probes check with the loaded lines and the thrown accession, gnl, lcl and giim
    /// ones), though no database makes it throw: NCBI finds a GI by its ISAM index, never by the
    /// accession table whose absence throws (review 2 of pfoff-3, L3).
    pub fn thrown_description(&self) -> (Vec<u8>, Defline) {
        let mut start = Vec::new();
        let gi = self.alignment_id.gi().unwrap_or(0);
        let wrote_gi = description_id(&self.alignment_id, gi, self.show_gis, &mut start);
        let mut title = b"Unknown".to_vec();
        if wrote_gi {
            title.splice(0..0, *b"  ");
        }
        (
            start,
            Defline {
                text: title,
                reads_past_end: false,
            },
        )
    }

    /// Whether the Seq-align id is a local id, which `CAlignFormatUtil::GetDbType` does not look
    /// up (`AfterThrow::alignment_end`).
    pub fn alignment_id_is_local(&self) -> bool {
        matches!(self.alignment_id, ParsedSeqId::Local(_))
    }

    /// The start of the outfmt 0 description line (`description`).
    pub fn description_start(&self) -> Vec<u8> {
        self.description_id().0
    }

    /// `description`'s start and whether `gi|N|` was written.
    fn description_id(&self) -> (Vec<u8>, bool) {
        let mut start = Vec::new();
        let mut wrote_gi = false;
        if let Some((defline, found)) = self.display_defline() {
            if let Some(id) = defline.display_id() {
                let gi = if found { defline.gi() } else { 0 };
                wrote_gi = description_id(id, gi, self.show_gis, &mut start);
            }
        }
        (start, wrote_gi)
    }

    /// The definition line whose id the description line shows (`DeflineSet::display_defline`
    /// of the Seq-align id) and whether it holds the Seq-align id.
    fn display_defline(&self) -> Option<(&BlastDefline, bool)> {
        let shown = self.set.display_defline(&self.alignment_id)?;
        Some((&self.set.deflines[shown.index], shown.found))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:432-433, 497-498
    // ```c++
    //     //get id (sdl->id, sdl-gi)
    //     sdl->id = CAlignFormatUtil::GetDisplayIds(handle,aln_id,use_this_seqid,&sdl->gi,&sdl->taxid,&sdl->textSeqID);
    //     ...
    //     //get defline
    //     sdl->defline = CDeflineGenerator().GenerateDefline(m_ScopeRef->GetBioseqHandle(*(sdl->id)), sequence::CDeflineGenerator::fLeavePrefixSuffix);
    // ```
    /// The display id: the id of the description line (`GetDisplayIds`' `wid`, the `WorstRank`
    /// best id of `display_defline`), by which the description table (outfmt 0 to 4) looks the
    /// title up again in NCBI's scope, which finds a query of the batch with that id before the
    /// database sequence (`common::database::refuse_display_id_of_a_query`).
    pub fn display_id(&self) -> Option<&ParsedSeqId> {
        self.display_defline()?.0.display_id()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2380-2384, 2451-2461
    // ```c++
    //                 if(alnDispParams) {
    //                     numBdl++;
    //                     if(isFirst){
    //                         out << ">";
    //                     } else{
    //                         out << " ";
    //     ...
    //                     if (out.tellp() > 1L) {
    //                         out << " ";
    //                     }
    // 					if(!alnDispParams->title.empty()) {
    //                         s_WrapOutputLine(out, (m_AlignOption&eHtml) ?
    //                                          CHTMLHelper::
    //                                          HTMLEncode(alnDispParams->title) :
    //                                          alnDispParams->title);
    //                     }
    //                     out<<"\n";
    //                     isFirst = false;
    // ```
    // `out` holds the heading of this subject only; the title of a definition line is its own,
    // or the subject Bioseq's generated title when it is empty (`x_FillDispParams`).
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:3166-3180
    // ```c++
    //     ITERATE(CSeq_align_set::Tdata, iter, alnset.Get()){
    //         const CSeq_align::TSegs& seg = (*iter)->GetSegs();
    //         if(seg.Which() == CSeq_align::C_Segs::e_Std){
    //             if(seg.GetStd().size() > 1){
    //                 //has more than one stdseg. Need to seperate as each
    //                 //is a distinct HSP
    //                 ITERATE (CSeq_align::C_Segs::TStd, iterStdseg, seg.GetStd()){
    //                     CRef<CSeq_align> aln(new CSeq_align);
    //                     if((*iterStdseg)->IsSetScores()){
    //                         aln->SetScore() = (*iterStdseg)->GetScores();
    //                     }
    //                     aln->SetSegs().SetStd().push_back(*iterStdseg);
    //                     alnSetRef->Set().push_back(aln);
    //                 }
    //
    //             } else {
    //                 alnSetRef->Set().push_back(*iter);
    //             }
    // ```
    /// The definition-line lines of the outfmt 0 alignment of the subject (each ending with a
    /// newline; `Length=` follows): the definition lines (`heading_shows`), the first after `>`
    /// and the others after a space, with its id (`alignment_line_id`) and its title wrapped at
    /// 60 bytes. `every_line`: the Seq-align has no `ext` when the report is made, so nothing
    /// limits the lines: a Seq-align of an ungapped search with two or more HSPs, which the
    /// formatter splits into new ones (`PrepareBlastUngappedSeqalign` above). The generated title
    /// is made once, for the first line without a title; where it throws, `Err` with the Seq-id
    /// text of NCBI's line in place of each alignment (`heading_skipped_id`).
    pub fn heading(&self, protein: bool, every_line: bool) -> Result<Vec<u8>, Vec<u8>> {
        let mut out = Vec::new();
        let mut generated = None;
        let mut first = true;
        for defline in &self.set.deflines {
            if !every_line && !self.heading_shows(defline) {
                continue;
            }
            let start = out.len();
            out.push(if first { b'>' } else { b' ' });
            alignment_line_id(defline, self.show_gis, &mut out);
            if !first || out.len() - start > 1 {
                out.push(b' ');
            }
            first = false;
            let title = match defline.title.as_deref() {
                Some(title) if !self.heading_generates(title) => title,
                _ => match &generated {
                    Some(text) => text,
                    None => match self.bioseq_title(protein, false) {
                        Ok(defline) => generated.insert(defline.text),
                        Err(UnknownEncoding) => return Err(self.skipped_id()),
                    },
                },
            };
            wrap_title(title, &mut out);
            out.push(b'\n');
        }
        Ok(out)
    }

    /// Whether the heading shows the generated title for a definition line with the title
    /// `title` (`x_FillDispParams`: an empty one; every line of a `-subject` record, whose
    /// title is the record's own, as `heading_skipped_id` has it).
    fn heading_generates(&self, title: &[u8]) -> bool {
        title.is_empty() || self.record
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1976-1981
    // ```c++
    //                 } catch (const CException&){
    //                     out << "Sequence with id "
    //                         << (avRef->GetSeqId(1)).GetSeqIdString().c_str()
    //                         <<" no longer exists in database...alignment skipped\n";
    //                     continue;
    //                 }
    // ```
    // `x_DisplayAlnvecInfo` fills each shown definition line (`x_FillDispParams`) before it
    // writes the heading, and a line without a title gets the Bioseq's generated title, whose
    // `NStr::HtmlDecode` throws where it cannot guess the encoding (a name of the taxonomy
    // database, line RG item T5; a `-subject` record's title, which its heading always generates:
    // probe `T/t5/sub`, `AB999999` for `gb|AB999999.1|`): no heading and no alignment are written,
    // the message instead, for each Seq-align of the subject (`previousId` is not set, so the
    // next one throws too).
    /// The Seq-id text (`CSeq_id::GetSeqIdString` of the Seq-align id) that outfmt 0 writes in
    /// place of the alignments of a subject whose heading (`heading`, with `every_line`) shows the
    /// generated title and NCBI cannot make it (`UnknownEncoding`); `None` otherwise.
    pub fn heading_skipped_id(&self, protein: bool, every_line: bool) -> Option<Vec<u8>> {
        let generated = self.set.deflines.iter().any(|defline| {
            (every_line || self.heading_shows(defline))
                && self.heading_generates(defline.title.as_deref().unwrap_or_default())
        });
        if !generated || self.bioseq_title(protein, false).is_ok() {
            return None;
        }
        Some(self.skipped_id())
    }

    /// `CSeq_id::GetSeqIdString` of the Seq-align id: the id of NCBI's line in place of each
    /// alignment of a subject whose heading it cannot make (`heading_skipped_id`).
    pub fn skipped_id(&self) -> Vec<u8> {
        let mut id = Vec::new();
        label::content(&self.alignment_id, false, &mut id);
        id
    }

    /// Whether the `Score`-best id of the subject Bioseq is a local id (XML2 then writes no
    /// accession of a `-subject` record, `SubjectNames::of_believed`).
    pub fn has_local_best_id(&self) -> bool {
        matches!(
            find_best(&self.target().ids, text_score),
            Some(ParsedSeqId::Local(_))
        )
    }

    /// The XML `Hit_id` and `Hit_accession` of a `-subject` record (`of_record`; `Hit_def` is the
    /// generated title, `bioseq_title`).
    pub fn record_xml_ids(&self) -> (Vec<u8>, Vec<u8>) {
        let mut id = Vec::new();
        label::seqid_list_string(&self.target().ids, true, &mut id);
        let mut accession = Vec::new();
        if let Some(best) = find_best(&self.target().ids, text_score) {
            label::format_label(best, false, &mut accession);
        }
        (id, accession)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2196-2202
    // ```c++
    //     CRef<CSeq_id> wid = FindBestChoice(ids, CSeq_id::WorstRank);
    //     TGi gi_in_use_this_gi = ZERO_GI;
    //     bool isGiList = false;
    //     bool match = CAlignFormatUtil::MatchSeqInSeqList(gi, wid, use_this_seqid,&isGiList);
    //     if(match && isGiList) gi_in_use_this_gi = gi;
    //
    //     if(use_this_seqid.empty() || match) {
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4648-4664
    // ```c++
    // bool CAlignFormatUtil::MatchSeqInSeqList(TGi cur_gi, CRef<CSeq_id> &seqID, list<string> &use_this_seq,bool *isGiList)
    // {
    //     bool found = false;
    //     bool isGi = false;
    //
    //     string curSeqID = CAlignFormatUtil::GetLabel(seqID,true); //uses GetSeqIdString(true)
    //     ITERATE(list<string>, iter_seq, use_this_seq){
    //         isGi = false;
    //         string useThisSeq = s_UseThisSeqToTextSeqID(*iter_seq, isGi);
    //         if((isGi && cur_gi == NStr::StringToNumeric<TGi>((useThisSeq))) || (!isGi && curSeqID == useThisSeq)){
    //             found = true;
    //             break;
    //          }
    //     }
    // ```
    /// Whether the outfmt 0 alignment heading shows `defline`: under a positive list only a line
    /// that `use_this_seqid` names (a line of a local id alone is not shown, nor one without a
    /// GI when the strings are GIs; probe G05 of line DI-8).
    fn heading_shows(&self, defline: &BlastDefline) -> bool {
        if self.use_this_seqid.is_empty() {
            return true;
        }
        // `CAlignFormatUtil::GetLabel`: the tag of a `general` id, else `GetSeqIdString(true)`.
        let mut label_of_id = Vec::new();
        match defline.display_id() {
            Some(id @ ParsedSeqId::General(_)) => {
                label::content_general_db_is_content(id, &mut label_of_id);
            }
            Some(id) => label::content(id, true, &mut label_of_id),
            None => {}
        }
        // `NStr::StringToNumeric<TGi>(useThisSeq)` (`NStr::StringToInt8`, flags 0) of a GI
        // that `under_positive_list` wrote in decimal.
        let gi = defline.gi();
        self.use_this_seqid
            .iter()
            .any(|string| match string.strip_prefix(b"gi:") {
                Some(value) => crate::ncbi_str::string_to_int8(value, 0).is_ok_and(|n| n == gi),
                None => string.strip_prefix(b"seqid:") == Some(label_of_id.as_slice()),
            })
    }

    /// The outfmt 5 `Hit_id` and `Hit_def` (`DeflineSet::hit_id_and_def` with the GIs shown;
    /// an empty `Hit_def` is written `No definition line`) and `Hit_accession` (the label of
    /// the `Score`-best id of the subject Bioseq, without its version; a `general` id's tag).
    pub fn xml_hit(&self) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
        let (id, def) = self.set.hit_id_and_def(true);
        let mut accession = Vec::new();
        if let Some(best) = find_best(&self.target().ids, text_score) {
            label::format_label(best, false, &mut accession);
        }
        (id, def, accession)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:365-388
    // ```c++
    //         	ITERATE(list<CRef<CBlast_def_line> >, itr, bdlRef->Get()) {
    //         		const CBlast_def_line & defline = **itr;
    //         		CRef<blastxml2::CHitDescr> hit_exp(new blastxml2::CHitDescr);
    //         		hit_exp->SetId(CShowBlastDefline::GetSeqIdListString(defline.GetSeqid(), true));
    //
    //         		 CRef<CSeq_id> best_id = FindBestChoice(defline.GetSeqid(), CSeq_id::Score);
    //         		 CSeq_id_Handle id_handle = CSeq_id_Handle::GetHandle(*best_id);
    //         		 string accession = CAlignFormatUtil::GetLabel(id_handle.GetSeqId());
    //         		 if(accession != kEmptyStr)
    //         			 hit_exp->SetAccession(accession);
    //
    //         		if(defline.IsSetTitle())
    //         			hit_exp->SetTitle(defline.GetTitle());
    //
    //         		if(defline.IsSetTaxid() && defline.GetTaxid() != ZERO_TAX_ID) {
    //         			TTaxId tax_id = defline.GetTaxid();
    //        				hit_exp->SetTaxid(tax_id);
    // ```
    /// The `HitDescr` of each definition line (XML2 and JSON), without the scientific names.
    pub fn hit_descriptions(&self) -> Vec<HitDescription> {
        self.set
            .deflines
            .iter()
            .map(|defline| {
                let mut id = Vec::new();
                label::seqid_list_string(&defline.ids, true, &mut id);
                let mut accession = Vec::new();
                if let Some(best) = defline.display_id() {
                    label::format_label(best, false, &mut accession);
                }
                HitDescription {
                    id,
                    accession,
                    title: defline.title.clone(),
                    taxid: defline.taxid.unwrap_or(0),
                    taxonomy: HitTaxonomy::Unnamed,
                }
            })
            .collect()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/taxFormat.cpp:935-976
    // ```c++
    // CTaxFormat::SSeqInfo *CTaxFormat::x_FillTaxDispParams(const CRef< CBlast_def_line > &bdl,
    //                                      const CBioseq_Handle& bsp_handle,
    //                                      double bits,
    //                                      double evalue,
    // 		                             list<string>& use_this_seqid)
    // {
    //     SSeqInfo *seqInfo = NULL;
    //
    // 	const list<CRef<CSeq_id> > ids = bdl->GetSeqid();
    // 	TGi gi =  CAlignFormatUtil::GetGiForSeqIdList(ids);
    //     CRef<CSeq_id> wid = FindBestChoice(ids, CSeq_id::WorstRank);
    //     bool match = CAlignFormatUtil::MatchSeqInSeqList(gi, wid, use_this_seqid);
    //
    //     if(use_this_seqid.empty() || match) {
    // 		...
    // 		seqInfo->seqID = FindBestChoice(ids, CSeq_id::WorstRank);
    // 		seqInfo->label =  CAlignFormatUtil::GetLabel(seqInfo->seqID);//Just accession without db part like ref| or pdbd|
    // ```
    /// The definition lines of the organism report (outfmt 18, `report::tax_report`): of a
    /// database sequence, each line that a positive list keeps (`heading_shows`) with the label of
    /// its `WorstRank`-best id without version, its title (the Bioseq's generated title where it
    /// has none) and its taxid; of a `-subject` record (no `Blast-def-line`), the Bioseq's
    /// `WorstRank`-best id, its generated title and taxid 0 (`TaxDefline::of_record`). A generated
    /// title that LOSAT does not write as NCBI is a `TaxTitle::Refused` (those of the outfmt 0
    /// heading, which shows the same title, `common::database::check_shown_parsed_subject`; the
    /// one made from a source that LOSAT refuses, `SourceTaxname::Refused`, only where the report
    /// shows it). A title made from the source names the organism of the Bioseq's first source
    /// (`with_source`). NCBI BLAST+ 2.17.0 keeps the positive list in an ungapped search too
    /// (probe `T/t2/fix1` `f1.d3.*.sl*`: tblastx and `blastn -ungapped` with `-seqidlist` show the
    /// listed line alone for a subject of 2 to 5 HSPs).
    pub fn tax_deflines(&self, protein: bool) -> Vec<crate::report::tax_report::TaxDefline> {
        use crate::report::tax_report::{generated_title, TaxDefline, TaxTitle};
        let label_of = |id: Option<&ParsedSeqId>| {
            let mut label = Vec::new();
            if let Some(id) = id {
                label::format_label(id, false, &mut label);
            }
            label
        };
        let refused = |reason: &str, shown_only: bool| TaxTitle::Refused {
            reason: reason.to_owned(),
            shown_only,
        };
        let generated = || -> TaxTitle {
            if let Err(reason) = self.source.name() {
                return refused(&format!("has {reason}"), true);
            }
            if let Err(reason) = self.generated_title_fetched() {
                return refused(&reason, false);
            }
            let title = self.target().title.as_deref().unwrap_or_default();
            generated_title(title, protein, &self.bioseq_facts())
        };
        if self.record {
            return vec![TaxDefline {
                label: label_of(find_best(&self.target().ids, text_score)),
                title: generated(),
                taxid: 0,
            }];
        }
        self.set
            .deflines
            .iter()
            .filter(|defline| self.heading_shows(defline))
            .map(|defline| TaxDefline {
                label: label_of(defline.display_id()),
                title: match defline.title.as_deref() {
                    Some(title) if !title.is_empty() => TaxTitle::Text(title.to_vec()),
                    _ => generated(),
                },
                taxid: defline.taxid.unwrap_or(0),
            })
            .collect()
    }

    /// The FASTA form of the Seq-align id (`CSeq_id::AsFastaString`): the id of the `HitDescr`
    /// that NCBI writes in place of the rest when the taxonomy database does not know a taxid
    /// (blastxml2_format.cpp:422-428), and the id that LOSAT's refusals name.
    pub fn alignment_fasta(&self) -> Vec<u8> {
        let mut out = Vec::new();
        label::fasta(&self.alignment_id, &mut out);
        out
    }
}

/// The ids of every definition line of a subject (`ParsedSubject::all_tabular_ids`); where
/// `failed`, those of the lines before the one whose `GetSeqIdList` throws.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AllTabularIds {
    pub ids: Vec<super::tabular::DeflineIds>,
    pub failed: bool,
}

/// Whether `GetSeqIdList` replaces one of `ids` with the first word of the generated title
/// (`seq_id_list`: a `gnl|BL_ORD_ID` or `lcl|Subject_` id).
fn replaces_an_id(ids: &[ParsedSeqId]) -> bool {
    matches!(seq_id_list(ids, b"x"), std::borrow::Cow::Owned(_))
}

/// The tabular ids of one definition line after `seq_id_list` with `word`.
fn defline_tabular_ids(defline: &BlastDefline, word: &[u8]) -> super::tabular::DeflineIds {
    let ids = seq_id_list(&defline.ids, word);
    let field = |kind| {
        let mut out = Vec::new();
        tabular_id(&ids, kind, &mut out);
        out
    };
    super::tabular::DeflineIds {
        id: field(TabularId::Full),
        gi: defline.gi(),
        accession: field(TabularId::Accession),
    }
}

/// The subjects of a `-subject` run under `-parse_deflines`, by subject index: `records` gives the
/// ids (`None`: shown as a record without the option) and title of each. LOSAT's refusals of them
/// in outfmt 0 are `check_pairwise_records`.
pub fn parsed_records<'a>(
    records: impl Iterator<
        Item = (
            Option<&'a crate::blastinput::fasta_reader::believed::BelievedIds>,
            &'a [u8],
        ),
    >,
    show_gis: bool,
) -> Vec<Option<std::sync::Arc<ParsedSubject>>> {
    records
        .map(|(ids, title)| {
            ids.map(|ids| {
                std::sync::Arc::new(ParsedSubject::of_record(ids, title).showing_gis(show_gis))
            })
        })
        .collect()
}

/// LOSAT's refusals of outfmt 0 of the subjects of `parsed_records`, in their order: those that
/// outfmt 0 shows by their parsed ids (`ParsedSubject::shown_as_record`) go through
/// `check_pairwise_record`. A call passes the result to the refusal of each of its formats that
/// shows the titles (`FormatRefusal::check`), made right after `parsed_records`; `program` names
/// the program.
pub fn check_pairwise_records(
    subjects: &[Option<std::sync::Arc<ParsedSubject>>],
    protein: bool,
    program: &str,
) -> anyhow::Result<()> {
    subjects
        .iter()
        .flatten()
        .filter(|subject| !subject.shown_as_record())
        .try_for_each(|subject| check_pairwise_record(subject, protein, program))
}

/// LOSAT's refusals of a `-subject` record in outfmt 0, as of a database sequence
/// (`algorithm::common::database::check_shown_parsed_subject`): a protein title that NCBI's
/// `x_CleanAndCompress` reads past the end of, and a description id longer than 65 bytes.
fn check_pairwise_record(
    subject: &ParsedSubject,
    protein: bool,
    program: &str,
) -> anyhow::Result<()> {
    let (start, title) = subject.description(protein);
    let name = || String::from_utf8_lossy(&start).trim_end().to_string();
    if protein && title.reads_past_end {
        anyhow::bail!(
            "the subject record {} has a title of punctuation that NCBI BLAST+'s \
             x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes in outfmt 0), \
             which is not supported by LOSAT's {program}",
            name()
        );
    }
    if start.len() > 65 {
        anyhow::bail!(
            "the subject record {} has an outfmt 0 description id longer than 65 bytes (NCBI \
             BLAST+ 2.17.0 pads the description line without end), which is not supported by \
             LOSAT's {program}",
            name()
        );
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/accguide2.txt (the rows of the RefSeq
// kinds that `x_SetFlags` tests; `3+6` is a prefix of three characters and six digits)
// ```
// 3+6  AC_      refseq_chromosome
// 3+6  NC_      refseq_chromosome
// 3+6  NM_      refseq_mrna
// 3+9  NM_      refseq_mrna
// 3+6  NR_      refseq_ncrna
// 3+6  XM_      refseq_mrna_predicted
// 3+9  XM_      refseq_mrna_predicted
// ```
/// The RefSeq kind that `CSeq_id::IdentifyAccession` gives `accession` where `x_SetFlags` tests
/// it: a chromosome or ncRNA (`BioseqFacts::refseq_chromosome_or_ncrna`) or an mRNA
/// (`BioseqFacts::refseq_mrna`).
fn refseq_accession_kind(accession: &[u8]) -> Option<RefseqKind> {
    if accession.len() < 3 || accession[2] != b'_' {
        return None;
    }
    let digits = &accession[3..];
    if !digits.iter().all(u8::is_ascii_digit) {
        return None;
    }
    let prefix = [
        accession[0].to_ascii_uppercase(),
        accession[1].to_ascii_uppercase(),
    ];
    match &prefix {
        b"AC" | b"NC" | b"NR" if digits.len() == 6 => Some(RefseqKind::ChromosomeOrNcrna),
        b"NM" | b"XM" if digits.len() == 6 || digits.len() == 9 => Some(RefseqKind::Mrna),
        _ => None,
    }
}

/// A RefSeq kind whose empty title `GenerateDefline` makes with its own function
/// (`refseq_accession_kind`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum RefseqKind {
    /// `m_IsNC` or `m_IsNR`.
    ChromosomeOrNcrna,
    /// `m_IsNM`.
    Mrna,
}

impl BioseqFacts {
    /// Sets the flag of the RefSeq kind of `accession` (`refseq_accession_kind`), if any.
    fn add_refseq_accession(&mut self, accession: Option<&[u8]>) {
        match accession.and_then(refseq_accession_kind) {
            Some(RefseqKind::ChromosomeOrNcrna) => self.refseq_chromosome_or_ncrna = true,
            Some(RefseqKind::Mrna) => self.refseq_mrna = true,
            None => {}
        }
    }
}

/// `BioseqFacts` of the Seq-ids `ids` of a Bioseq.
pub(crate) fn bioseq_facts(ids: &[ParsedSeqId]) -> BioseqFacts {
    let mut flags = BioseqFacts::default();
    for id in ids {
        match id {
            ParsedSeqId::Tpg(text) | ParsedSeqId::Tpe(text) | ParsedSeqId::Tpd(text) => {
                flags.third_party = true;
                flags.add_refseq_accession(text.accession.as_deref());
            }
            ParsedSeqId::Other(text)
            | ParsedSeqId::Genbank(text)
            | ParsedSeqId::Embl(text)
            | ParsedSeqId::Ddbj(text) => {
                flags.add_refseq_accession(text.accession.as_deref());
            }
            ParsedSeqId::Pdb(pdb) => {
                let (chain_id, chain) = flags.pdb.get_or_insert_with(|| (Vec::new(), 0));
                if let Some(id) = &pdb.chain_id {
                    *chain_id = id.clone();
                } else if let Some(value) = pdb.chain {
                    *chain = value;
                }
            }
            ParsedSeqId::Patent(patent) => {
                let number = match &patent.cit.id {
                    PatentNumber::Number(number) | PatentNumber::AppNumber(number) => {
                        number.clone()
                    }
                };
                flags.patent = Some((patent.seqid, patent.cit.country.clone(), number));
            }
            _ => {}
        }
    }
    flags
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:336-359
// ```c++
// static void s_WrapOutputLine(CNcbiOstream& out, const string& str)
// {
//     const int line_len = 60;
//     bool do_wrap = false;
//     int length = (int) str.size();
//     if (length > line_len) {
//         for (int i = 0; i < length; i ++){
//             if(i > 0 && i % line_len == 0){
//                 do_wrap = true;
//             }
//             out << str[i];
//             if(do_wrap && isspace((unsigned char) str[i])){
//                 out << "\n";
//                 do_wrap = false;
//             }
//         }
//     } else {
//         out << str;
//     }
// }
// ```
/// `s_WrapOutputLine`: `title`, with a newline after the first white space at or after each
/// multiple of 60 bytes.
pub fn wrap_title(title: &[u8], out: &mut Vec<u8>) {
    let mut do_wrap = false;
    for (i, &c) in title.iter().enumerate() {
        if i > 0 && i % 60 == 0 {
            do_wrap = true;
        }
        out.push(c);
        if do_wrap && (c.is_ascii_whitespace() || c == 11) {
            out.push(b'\n');
            do_wrap = false;
        }
    }
}

/// Whether NCBI's XML and outfmt 0 reports join the hits of `a` and `b` when they follow each
/// other in a hit list: their Seq-align ids match (`CSeq_id::Match`, which is
/// `Compare(..) == e_YES`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:478-487
/// ```c++
///             CSeq_align_set one_subject_alnset;
///             CConstRef<CSeq_id> current_id(s_GetSubjectId(*(*iter)));
///             for ( ; iter != alnset.Get().end(); ++iter) {
///                 CConstRef<CSeq_id> next_id(s_GetSubjectId(*(*iter)));
///                 if (!current_id->Match(*next_id)) {
///                     break;
///                 }
///                 one_subject_alnset.Set().push_back(*iter);
///             }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1750-1757
/// ```c++
///     ITERATE(CSeq_align_set::Tdata, iter, source.Get()){
///         const CSeq_id& cur_id = (*iter)->GetSeq_id(1);
///         if(previous_id.Empty()) {
///             temp =  new CSeq_align_set;
///             temp->Set().push_back(*iter);
///             target.push_back(temp);
///         } else if (cur_id.Match(*previous_id)){
///             temp->Set().push_back(*iter);
/// ```
pub fn joined_in_reports(a: &ParsedSubject, b: &ParsedSubject) -> bool {
    compare(&a.alignment_id, &b.alignment_id) == Sic::Yes
}

// No NCBI counterpart: LOSAT joins the hits that NCBI joins only where they are copies of one
// record; it refuses the others (`common::database::refuse_joined_hits`).
/// Whether `later`, which follows `earlier` in a hit list, is shown in `earlier`'s hit by the
/// reports that join following hits (`joined_in_reports`) as LOSAT writes them: the two are
/// copies of one record (the same definition lines but for their taxids). NCBI shows that hit
/// with the names of the copy its lookup gives first, and each alignment with the letters of the
/// copy its scope holds for the Seq-align id; the reports that read the Bioseq refuse either
/// subject where another sequence holding its Seq-align id is not shown as it (`fetched`: other
/// letters, other taxids where they are written), so the copies' names and letters are the same.
pub fn joins_as_copy(earlier: &ParsedSubject, later: &ParsedSubject) -> bool {
    let without_taxid = |defline: &BlastDefline| BlastDefline {
        taxid: None,
        ..defline.clone()
    };
    joined_in_reports(earlier, later)
        && earlier.set.deflines.len() == later.set.deflines.len()
        && earlier
            .set
            .deflines
            .iter()
            .zip(&later.set.deflines)
            .all(|(a, b)| without_taxid(a) == without_taxid(b))
}

/// `joins_as_copy` of the subjects of ranks `earlier` and `later` of `parsed` (the subjects with
/// parsed Seq-ids by rank): false where either has none.
pub fn rank_joins_as_copy(
    parsed: &[Option<std::sync::Arc<ParsedSubject>>],
    earlier: usize,
    later: usize,
) -> bool {
    match (
        parsed.get(earlier).and_then(Option::as_deref),
        parsed.get(later).and_then(Option::as_deref),
    ) {
        (Some(earlier), Some(later)) => joins_as_copy(earlier, later),
        _ => false,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:476-506
// ```c++
//     while (iter != alnset->Get().end()) {
//         CRef<blastxml2::CHit> new_hit(new blastxml2::CHit);
//         ...
//         if ((*iter)->GetSegs().IsDisc()) {
//             s_SeqAlignToXMLHit(new_hit, *(*iter), scope, matrix, mask_info,
//                                ungapped, master_gentice_code, slave_genetic_code, hasTaxDB);
//             ++iter;
//         } else {
//             ...
//             for ( ; iter != alnset->Get().end(); ++iter) {
//                 CConstRef<CSeq_id> next_id(s_GetSubjectId(*(*iter)));
//                 if (!current_id->Match(*next_id)) {
//                     break;
//                 }
//                 one_subject_alnset.Set().push_back(*iter);
//             }
//             ...
//             s_SeqAlignToXMLHit(new_hit, disc_align_wrap, scope, matrix,
//                                mask_info, ungapped, master_gentice_code, slave_genetic_code, hasTaxDB);
//         }
//
//         hits.push_back(new_hit);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:357-360, 412-417
// (quoted at `ParsedSubject::lookup_throws`: the `Hit` loads its subject by the Seq-align id of
// its first Seq-align, and a lookup that throws there throws again outside the hit's `catch`)
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:750-760
// ```c++
// BlastXML2_FormatReport(const IBlastXML2ReportData* data, CNcbiOstream *out_stream)
// {
// 	blastxml2::CBlastOutput2 bxmlout;
// 	try {
// 		s_FillBlastOutput(bxmlout, data);
// 		s_WriteXML2ObjectNoHeader(bxmlout, out_stream);
// 	}
// 	catch(CException &e){
// 	    ERR_POST(Error << e.GetMsg() << e.what() );
// 	    return;
// 	}
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:833-870 (quoted at
// `report::tabular::TabularSpec::row_throw`: a table's row whose subject block catches the throw)
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:891-898 (quoted at
// `report::tabular::TabularSpec::reads_subject`: `slen` read outside any `catch`)
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:1586-1592 (quoted at
// `ParsedSubject::thrown_description`: the description line's `catch`)
// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4188-4197
// ```c++
// CAlignFormatUtil::DbType CAlignFormatUtil::GetDbType(const CSeq_align_set& actual_aln_list, CScope & scope)
// {
//     //determine if the database has gi by looking at the 1st hit.
//     //Could be wrong but simple for now
//     DbType type = eDbTypeNotSet;
//     CRef<CSeq_align> first_aln = actual_aln_list.Get().front();
//     const CSeq_id& subject_id = first_aln->GetSeq_id(1);
//
//     if (subject_id.Which() != CSeq_id::e_Local){
//         const CBioseq_Handle& handleTemp  = scope.GetBioseqHandle(subject_id);
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1881 (`x_InitAlignParams`,
// called by `DisplaySeqalign` before the alignments of outfmt 0 to 4, without a `catch`)
// ```c++
//     m_CanRetrieveSeq = CAlignFormatUtil::GetDbType(actual_aln_list,m_Scope) == CAlignFormatUtil::eDbTypeNotSet ? false : true;
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1950-1981 (quoted at
// `report::pairwise::write_subject_alignments_part`: outfmt 0's `catch` of each alignment, which
// writes the skipped line), 1998-2003 (quoted at `report::x_pairwise_par::write_serial`: the merge
// of outfmt 1 to 4 loads each shown alignment's subject without a `catch`)
/// What a report does after a subject whose lookup throws where its batch's prefetch is off
/// (`ParsedSubject::lookup_throws`; DECISIONS 62 of line DBID): the one answer of each writer,
/// which [`HitLoads::of`] takes for the run ([`FormatThrow::of_run`]). The first three are ordered
/// by what the reports load after the throw (each loads at least what the one before it loads);
/// the pairwise reports' two (`is_pairwise`) load the subjects of their description lines and of
/// their alignments, an order of their own (sorted by `-sorthits`), which is in no such order with
/// the others ([`FormatThrow::and`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum AfterThrow {
    /// NCBI's application ends the run with the exception (`CATCH_ALL`'s `BLAST Database error:
    /// <message>`, exit status 2; `blast_xml::writer::run_ending_exception`): outfmt 13 and 14,
    /// and a table whose rows read the subject only after its subject block (`slen`, the letters
    /// and frames), at the subject's first row that reads it.
    EndsRun,
    /// The report function catches it and leaves the query's report out (outfmt 15 and 16,
    /// `BlastXML2_FormatReport`, `BlastJSON_FormatReport`); the next query's report goes on.
    EndsReport,
    /// The `catch` of the table's subject block ends the row's fields, the row is written with
    /// what the query's formatter kept (`report::tabular::TabularWriter::row`), and the report
    /// goes on: a table with a field of that block (outfmt 6, 7, 10 and 20); the organism report
    /// (outfmt 18, `CTaxFormat::x_InitTaxInfoMap`'s `catch`) leaves the HSP out and goes on.
    WritesRow,
    /// The pairwise report, outfmt 0: the `catch` of the description line writes it with the
    /// Seq-align id and `Unknown` (`ParsedSubject::thrown_description`) and the report goes on;
    /// then `GetDbType` loads the subject of the first alignment shown, without a `catch` for an
    /// id that is not local: where it throws, the run ends after the query's description lines
    /// (`alignment_end`); otherwise the `catch` of each alignment writes NCBI's skipped line in
    /// place of each alignment of a subject that throws (`report::pairwise::DatabaseHeading::
    /// Thrown`), and the report goes on.
    SkipsAlignment,
    /// The query-anchored reports, outfmt 1 to 4: the description lines as outfmt 0's
    /// (`SkipsAlignment`); then `GetDbType` and the merge load the subject of each alignment shown,
    /// in order, without a `catch`: the run ends after the query's description lines at the first
    /// that throws (`alignment_end`).
    EndsAtAlignment,
}

impl AfterThrow {
    /// Whether the writer is the pairwise report's (outfmt 0 to 4: `SkipsAlignment`,
    /// `EndsAtAlignment`), whose reports read what their writers show ([`PairwiseReads`]).
    pub fn is_pairwise(self) -> bool {
        matches!(
            self,
            AfterThrow::SkipsAlignment | AfterThrow::EndsAtAlignment
        )
    }

    /// Where the alignments of a pairwise report (outfmt 0 to 4) end the run: the position, among
    /// the subjects of the alignments it shows, in order (`aligned`), of the one whose load ends it
    /// (`throws`: its load throws; `local`: its Seq-align id is a local id, which `GetDbType` does
    /// not look up). `SkipsAlignment`: the first, where it throws and is not local; `EndsAtAlignment`:
    /// the first that throws; `None` for the other writers and where none ends it. The one decision
    /// of both the model of the loads ([`HitLoads::of`], with the lookups) and the writer
    /// (`report::x_pairwise_par`, with the loads).
    pub fn alignment_end<K: Copy>(
        self,
        aligned: impl IntoIterator<Item = K>,
        throws: impl Fn(K) -> bool,
        local: impl Fn(K) -> bool,
    ) -> Option<usize> {
        let mut aligned = aligned.into_iter();
        match self {
            AfterThrow::SkipsAlignment => aligned
                .next()
                .filter(|&first| throws(first) && !local(first))
                .map(|_| 0),
            AfterThrow::EndsAtAlignment => aligned.position(throws),
            AfterThrow::EndsRun | AfterThrow::EndsReport | AfterThrow::WritesRow => None,
        }
    }
}

/// The counts of the description lines and the alignments that a pairwise report (outfmt 0 to 4)
/// shows of each query (`-num_descriptions`, `-num_alignments`;
/// `report::x_pairwise_par::described_subjects`, `QueryGroup::aligned`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Shown {
    pub descriptions: usize,
    pub alignments: usize,
}

/// What one query's pairwise report (outfmt 0 to 4) reads where its batch's prefetch is off, as
/// its writer reads it (`report::x_pairwise_par::QueryGroup::read_order`, for [`HitLoads::of`] and
/// the texts of its reads): its hits in its order, each the subjects that it shows (its own
/// first, then those of the copies of its record shown in it, each once), the hits whose
/// description lines it shows (indices of `hits`, in the order of the lines:
/// `x_pairwise_par::described_subjects`), and the number of its first hits whose alignments it
/// shows (`QueryGroup::aligned`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct PairwiseReads {
    pub hits: Vec<Vec<usize>>,
    pub described: Vec<usize>,
    pub aligned: usize,
}

impl PairwiseReads {
    /// The subjects in the order that the report reads them, each once: those of its description
    /// lines, then those of its alignments (the texts of the reads one by one,
    /// `common::database::subject_read_texts`). A subject of a hit that neither shows is never read.
    pub fn subjects(&self) -> Vec<usize> {
        let mut read = vec![false; self.hits.len()];
        let mut subjects = Vec::new();
        for hit in self.described.iter().copied().chain(0..self.aligned) {
            if !std::mem::replace(&mut read[hit], true) {
                subjects.extend(&self.hits[hit]);
            }
        }
        subjects
    }
}

/// The subjects of one query's report for [`HitLoads::of`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ReportKeys {
    /// The subject key of each of the report's entries (HSP lists or HSPs) in its order (`None`
    /// for an entry that shows nothing), of which the entries make `Hit`s as the XML writers make
    /// them (`blast_xml::hit_ranges`); of a pairwise report (outfmt 0 to 4), those of a batch that
    /// prefetches, which loads every subject, in any order.
    Entries(Vec<Option<usize>>),
    /// The reads of a pairwise report whose batch reads one by one ([`PairwiseReads`]).
    Pairwise(PairwiseReads),
}

impl From<Vec<Option<usize>>> for ReportKeys {
    fn from(entries: Vec<Option<usize>>) -> Self {
        ReportKeys::Entries(entries)
    }
}

/// What the reports of a run do after a subject whose lookup throws where its batch's prefetch is
/// off ([`FormatThrow::of_run`]): `after` ([`AfterThrow`]), at the rows that load the subject,
/// `rows`: every row, or for a table that reads the subject for `nident` alone the rows of HSPs
/// with identities (`report::tabular::SubjectFetch::WithIdentities`, `TabularSpec::throw_rows`);
/// a row of no identities there reads nothing of the subject, so its load neither throws nor ends
/// anything (it prints the 0 of `x_ResetFields`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RunThrow {
    pub after: AfterThrow,
    pub rows: super::tabular::SubjectFetch,
}

impl RunThrow {
    /// What a format whose every row (or `Hit`) loads its subject does: `after`.
    pub fn every_row(after: AfterThrow) -> Self {
        RunThrow {
            after,
            rows: super::tabular::SubjectFetch::Every,
        }
    }
}

/// What one format of a run does with a subject whose lookup throws where its batch's prefetch is
/// off ([`RunThrow`]), as [`FormatThrow::of_run`] takes it for the run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum FormatThrow {
    /// LOSAT refuses such a subject (`ParsedSubject::fetched`): a format whose writer LOSAT does
    /// not make write it as NCBI does (outfmt 5, 17).
    Refused,
    /// The format reads nothing of the subject (a table of fields that read nothing of it; the
    /// Seq-annots and the archive), so NCBI never looks it up there.
    #[default]
    Unread,
    /// The format writes it as NCBI does, then does as the [`RunThrow`] says.
    Writes(RunThrow),
}

impl FormatThrow {
    /// What the format `format` does, with its field list `spec` (`None` for a format without
    /// one, and for the default outfmt 6 and 7, whose `saccver` is a field of the subject block),
    /// in `program`'s reports.
    pub fn of(
        format: crate::blastinput::app::ReportFormat,
        spec: Option<&super::tabular::TabularSpec>,
        program: crate::blastinput::app::Program,
    ) -> Self {
        use crate::blastinput::app::ReportFormat;
        if let Some(kind) = super::blast_xml2::Xml2Kind::of(format) {
            return FormatThrow::Writes(RunThrow::every_row(kind.after_throw()));
        }
        if let Some(view) = super::pairwise::PairwiseView::of(format) {
            return FormatThrow::Writes(RunThrow::every_row(view.after_throw()));
        }
        match format {
            ReportFormat::Tabular
            | ReportFormat::TabularWithComments
            | ReportFormat::CommaSeparatedValues
            | ReportFormat::CommaSeparatedValuesWithHeader
            | ReportFormat::TaxFormat => spec.map_or(
                FormatThrow::Writes(RunThrow::every_row(AfterThrow::WritesRow)),
                |spec| spec.format_throw(program),
            ),
            ReportFormat::AsnText
            | ReportFormat::AsnBinary
            | ReportFormat::ArchiveFormat
            | ReportFormat::JsonSeqalign => FormatThrow::Unread,
            _ => FormatThrow::Refused,
        }
    }

    /// What the `-outfmt` value `outfmt` does in `program`'s reports, as NCBI reads it
    /// (`app::parse_formatting_string`; a value that NCBI refuses is refused).
    pub fn of_outfmt(outfmt: &str, program: crate::blastinput::app::Program) -> Self {
        use crate::blastinput::app::{parse_formatting_string, ReportFormat};
        let Some((format, choice)) = parse_formatting_string(outfmt).ok().and_then(|choice| {
            ReportFormat::from_number(choice.number).map(|format| (format, choice))
        }) else {
            return FormatThrow::Refused;
        };
        let spec = super::tabular::TabularSpec::of(format, &choice);
        Self::of(format, spec.as_ref(), program)
    }

    /// The answer of two formats of one run: a refusal by either refuses, a format that reads
    /// nothing leaves the other's answer, and of two that write, the one that loads more after the
    /// throw (a subject is loaded where some report loads it), at the rows that either loads at.
    /// A pairwise report (outfmt 0 to 4, `AfterThrow::is_pairwise`) with a format of another
    /// writer refuses: it loads in an order of its own (the library's calls of several formats).
    pub fn and(self, other: Self) -> Self {
        use super::tabular::SubjectFetch;
        match (self, other) {
            (FormatThrow::Refused, _) | (_, FormatThrow::Refused) => FormatThrow::Refused,
            (FormatThrow::Unread, answer) | (answer, FormatThrow::Unread) => answer,
            (FormatThrow::Writes(one), FormatThrow::Writes(other))
                if one.after != other.after
                    && (one.after.is_pairwise() || other.after.is_pairwise()) =>
            {
                FormatThrow::Refused
            }
            (FormatThrow::Writes(one), FormatThrow::Writes(other)) => {
                FormatThrow::Writes(RunThrow {
                    after: one.after.max(other.after),
                    rows: if one.rows == other.rows {
                        one.rows
                    } else {
                        SubjectFetch::Every
                    },
                })
            }
        }
    }

    /// What the reports of a run do after a subject whose lookup throws ([`FormatThrow::and`] of
    /// its formats): `None` where a format refuses it or no format reads it, so that every subject
    /// is loaded ([`HitLoads::of`]). The one answer of every program.
    pub fn of_run(formats: impl IntoIterator<Item = Self>) -> Option<RunThrow> {
        match formats
            .into_iter()
            .fold(FormatThrow::Unread, FormatThrow::and)
        {
            FormatThrow::Writes(run) => Some(run),
            FormatThrow::Refused | FormatThrow::Unread => None,
        }
    }

    /// [`FormatThrow::of_run`] of the `-outfmt` values `outfmts` of `program`'s run.
    pub fn of_outfmts<'a>(
        outfmts: impl IntoIterator<Item = &'a str>,
        program: crate::blastinput::app::Program,
    ) -> Option<RunThrow> {
        Self::of_run(
            outfmts
                .into_iter()
                .map(|outfmt| Self::of_outfmt(outfmt, program)),
        )
    }
}

/// What NCBI's reports do at the `Hit` (outfmt 13 to 16), the rows (outfmt 6, 7, 10 and 20) or
/// the organism report's HSPs (outfmt 18) of a subject: the one answer to "this subject's lookup
/// throws there" and "this subject is never loaded", which the refusals
/// (`ParsedSubject::fetched`, the tables' rows at which NCBI has no result,
/// `report::tabular::check_failed_lookups`), the names of the `Hit`s
/// (`blast_xml::subjects::HitNames`), the tables' rows (`report::tabular::SubjectLookup::Throws`)
/// and their names and taxonomy (`common::tabular::tabular_taxonomy`, the taxonomy database's
/// lookups), the organism report (`tax_report::SubjectTable`), the taxonomy reads
/// (`common::database::read_reported_taxonomy`, the texts of the subjects that a report reads one
/// by one) and the run's record of the loads (`TaxdbRun::note_loaded`, the stop's `before`) take
/// ([`HitLoads::of`]).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HitLoad {
    /// A report loads the subject (or the batch's prefetch does), as every subject of a run that
    /// writes another format.
    #[default]
    Loaded,
    /// The subject's lookup throws where a report loads it (`lookup_throws`): NCBI loads nothing
    /// for it, its rows are written as NCBI's formatter leaves them, and the report goes on or
    /// ends at its first row that loads it as the run's [`RunThrow`] says.
    Throws,
    /// No report loads the subject: it is shown only after a subject that throws where the
    /// reports end there, or only in rows that read nothing of it (a table of `nident` alone, HSPs
    /// of no identities) where its lookup would throw.
    Unread,
}

/// [`HitLoad`] of each subject of a call's reports, by the key that the caller gives its subjects
/// (a rank, or an ordinal id): those whose lookup throws where a report reaches them, with the
/// exception, and those that no report loads; every other subject is loaded. Empty (every subject
/// loaded) where a format of the run refuses a lookup that throws, or no report reads one by one.
#[derive(Debug, Default)]
pub struct HitLoads {
    throws: std::collections::HashMap<usize, crate::ncbi_exception::NcbiException>,
    unread: std::collections::HashSet<usize>,
}

impl HitLoads {
    /// The loads of the query reports `reports`, in the order of the reports: whether the
    /// query's batch has its prefetch off, so that its report loads each subject where it reaches
    /// it (`one_by_one`, `common::database::BatchPrefetch::reads_one_by_one`), and the subject key
    /// of each of its HSP lists (or HSPs) in the report's order (`None` for an entry that shows
    /// nothing). The entries make `Hit`s as the XML writers make them (`blast_xml::hit_ranges`: a
    /// copy of one record that follows another joins its `Hit`, `joins_as_copy`; the entries of
    /// one subject that follow each other are one). In a report that reads one by one, a `Hit`
    /// whose subject's lookup throws (`lookup_throws` of `parsed(key)`) throws (its subjects are
    /// not loaded) at its first row that loads the subject: every row for the run's `rows`
    /// `Every`, the rows of HSPs with identities for `WithIdentities` (`identities(report, entry)`:
    /// the largest of the identities of the HSPs of entry `entry` of report `report`, the
    /// Seq-align's `num_ident`, asked only there); a `Hit` none of whose rows loads it is not loaded
    /// and ends nothing (`HitLoad::Unread`, unless another report loads it or throws there). What
    /// comes after a `Hit` that throws is the run's `after` ([`FormatThrow::of_run`]):
    /// - `WritesRow` (a table of a field of the subject block, the organism report): each row looks
    ///   its own Seq-align id up, so every subject whose lookup throws throws, every other one
    ///   loads;
    /// - `EndsReport` (outfmt 15 and 16): the `Hit`s after it in the report load nothing, the next
    ///   report goes on;
    /// - `EndsRun` (outfmt 13 and 14, a table that reads the subject after its subject block):
    ///   nothing after it in the run loads (`ended`, carried from the reports before: see
    ///   [`HitLoads::after`]);
    /// - `SkipsAlignment` and `EndsAtAlignment` (outfmt 0 to 4): a report that reads one by one
    ///   gives its reads as its writer makes them (`ReportKeys::Pairwise`,
    ///   `x_pairwise_par::QueryGroup::read_order`): each description line loads its subject or
    ///   throws, then the alignments load theirs up to the one that ends the run, if any
    ///   (`pairwise_report`); one whose batch prefetches gives its entries in any order.
    ///
    /// A report whose batch prefetches loads every subject (a lookup that throws stops the
    /// prefetch first, `common::database::batch_prefetch_stop`). A subject is loaded where some
    /// report loads it; one that throws in a report throws in every report that reaches it.
    /// `run` `None`: a format of the run refuses the throw or none reads the subjects, and every
    /// subject is loaded (`ParsedSubject::fetched` refuses such a subject where a format reads it).
    pub fn of<'p, K: Into<ReportKeys>>(
        run: Option<RunThrow>,
        reports: impl IntoIterator<Item = (bool, K)>,
        parsed: impl Fn(usize) -> Option<&'p ParsedSubject>,
        identities: impl Fn(usize, usize) -> usize,
    ) -> Self {
        Self::after(
            &std::cell::Cell::new(false),
            run,
            reports,
            parsed,
            identities,
        )
    }

    /// [`HitLoads::of`] of the reports `reports` that follow those whose loads ended the run where
    /// `ended` holds true (an `EndsRun` subject that threw), which it holds after these reports
    /// too ([`RunLoads::loads`]).
    fn after<'p, K: Into<ReportKeys>>(
        ended: &std::cell::Cell<bool>,
        run: Option<RunThrow>,
        reports: impl IntoIterator<Item = (bool, K)>,
        parsed: impl Fn(usize) -> Option<&'p ParsedSubject>,
        identities: impl Fn(usize, usize) -> usize,
    ) -> Self {
        let mut loads = HitLoads::default();
        let Some(run) = run else {
            return loads;
        };
        let mut loaded = std::collections::HashSet::new();
        let mut passed = std::collections::HashSet::new();
        for (report, (one_by_one, keys)) in reports.into_iter().enumerate() {
            let subjects = match keys.into() {
                ReportKeys::Pairwise(reads) => {
                    debug_assert!(one_by_one && run.after.is_pairwise());
                    loads.pairwise_report(
                        run.after,
                        &reads,
                        &parsed,
                        ended,
                        (&mut loaded, &mut passed),
                    );
                    continue;
                }
                ReportKeys::Entries(subjects) => subjects,
            };
            if run.after.is_pairwise() {
                // A pairwise report whose batch prefetches loads every subject that it shows,
                // unless the run ended before.
                debug_assert!(!one_by_one);
                let keys = subjects.into_iter().flatten();
                if ended.get() {
                    passed.extend(keys);
                } else {
                    loaded.extend(keys);
                }
                continue;
            }
            // A table and the organism report look each row's Seq-align id up.
            let hits = if run.after == AfterThrow::WritesRow {
                (0..subjects.len()).map(|entry| entry..entry + 1).collect()
            } else {
                super::blast_xml::hit_ranges(
                    &subjects,
                    |subject| *subject,
                    |earlier, later| match (parsed(earlier), parsed(later)) {
                        (Some(earlier), Some(later)) => joins_as_copy(earlier, later),
                        _ => false,
                    },
                )
            };
            // The report ended at a `Hit` that threw (`EndsReport`), or the run did (`EndsRun`).
            let mut thrown = ended.get();
            for range in hits {
                let Some(first) = subjects[range.start] else {
                    continue;
                };
                let keys = subjects[range.clone()].iter().flatten().copied();
                let exception = parsed(first)
                    .filter(|_| one_by_one && !thrown)
                    .and_then(ParsedSubject::lookup_throws);
                let Some(exception) = exception else {
                    if thrown {
                        passed.extend(keys);
                    } else {
                        loaded.extend(keys);
                    }
                    continue;
                };
                // The `Hit`'s first row that loads its subject (`rows`).
                let loads_it = run.rows == super::tabular::SubjectFetch::Every
                    || range.clone().any(|entry| {
                        subjects[entry].is_some() && run.rows.reads(identities(report, entry))
                    });
                if !loads_it {
                    passed.extend(keys);
                    continue;
                }
                loads.throws.insert(first, exception);
                thrown = run.after != AfterThrow::WritesRow;
                if run.after == AfterThrow::EndsRun {
                    ended.set(true);
                }
                passed.extend(keys.filter(|&key| key != first));
            }
        }
        loads.unread = passed
            .into_iter()
            .filter(|key| !loaded.contains(key) && !loads.throws.contains_key(key))
            .collect();
        loads
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:686-700
    // The description lines load the subjects that they show, in their order: those of
    // `x_pairwise_par::described_subjects` (psiblast's later rounds, `two_part_lines`).
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:1014-1040
    // (`PruneSeqalign` keeps the alignments of the first `m_NumAlignments` hits, quoted at
    // `report::x_pairwise_par::QueryGroup::aligned`)
    /// The loads of one pairwise report (outfmt 0 to 4, `AfterThrow::is_pairwise`, whose run does
    /// `after` after a throw) in a batch whose prefetch is off, whose reads are `reads` (its
    /// writer's own answers, `x_pairwise_par::QueryGroup::read_order`): unless the run ended
    /// before (`ended`), its description lines load the subjects of the hits that they show, then
    /// its alignments those of its first `reads.aligned` hits up to the one whose load ends the
    /// run (`AfterThrow::alignment_end`, which `ended` then records). A subject whose lookup throws
    /// (`lookup_throws`) loads nothing there, nor do the copies in its hit; every subject that
    /// neither loads is passed.
    fn pairwise_report<'p>(
        &mut self,
        after: AfterThrow,
        reads: &PairwiseReads,
        parsed: &impl Fn(usize) -> Option<&'p ParsedSubject>,
        ended: &std::cell::Cell<bool>,
        (loaded, passed): (
            &mut std::collections::HashSet<usize>,
            &mut std::collections::HashSet<usize>,
        ),
    ) {
        let first = |hit: usize| reads.hits[hit].first().copied().and_then(parsed);
        let throws = |hit: usize| first(hit).and_then(ParsedSubject::lookup_throws);
        let mut reached = vec![false; reads.hits.len()];
        if !ended.get() {
            let end = after.alignment_end(
                0..reads.aligned,
                |hit| throws(hit).is_some(),
                |hit| first(hit).is_some_and(ParsedSubject::alignment_id_is_local),
            );
            let aligned = 0..end.map_or(reads.aligned, |end| end + 1);
            for hit in reads.described.iter().copied().chain(aligned) {
                reached[hit] = true;
            }
            ended.set(end.is_some());
        }
        for (hit, keys) in reads.hits.iter().enumerate() {
            match (reached[hit], throws(hit), keys.split_first()) {
                (true, Some(exception), Some((&first, copies))) => {
                    self.throws.insert(first, exception);
                    passed.extend(copies);
                }
                (true, None, _) => loaded.extend(keys),
                _ => passed.extend(keys),
            }
        }
    }

    /// These loads with each subject's key renamed by `key` (an ordinal id to its rank; a key
    /// that `key` does not name is dropped, which reads as loaded).
    pub fn renamed(self, key: impl Fn(usize) -> Option<usize>) -> Self {
        HitLoads {
            throws: self
                .throws
                .into_iter()
                .filter_map(|(old, exception)| Some((key(old)?, exception)))
                .collect(),
            unread: self.unread.into_iter().filter_map(&key).collect(),
        }
    }

    /// What the reports do at the `Hit` of the subject `key`.
    pub fn load(&self, key: usize) -> HitLoad {
        if self.throws.contains_key(&key) {
            HitLoad::Throws
        } else if self.unread.contains(&key) {
            HitLoad::Unread
        } else {
            HitLoad::Loaded
        }
    }

    /// Whether a report (or a batch's prefetch) loads the subject `key`.
    pub fn loaded(&self, key: usize) -> bool {
        self.load(key) == HitLoad::Loaded
    }

    /// The exception that the `Hit` of the subject `key` throws (`HitLoad::Throws`).
    pub fn thrown(&self, key: usize) -> Option<&crate::ncbi_exception::NcbiException> {
        self.throws.get(&key)
    }

    /// Whether the load of some subject throws (`HitLoad::Throws`).
    pub fn any_throws(&self) -> bool {
        !self.throws.is_empty()
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203 (the `CATCH_ALL` of
// every application, quoted at `cli::NativeError::exception`): an exception that leaves a report
// ends the application; no later report, batch or load runs.
/// The run's record of what its reports do after a subject whose lookup throws where a batch's
/// prefetch is off (`run`, [`FormatThrow::of_run`]) and of whether a report's load ended the run
/// (an `EndsRun` subject that threw): the one record of the run's end where a program makes each
/// query's report's loads apart from the others (the texts of BLASTN's and TBLASTX's reports of a
/// batch whose prefetch is off, their `report_loads`), kept across the run's batches, so that no
/// report after the end loads anything and no batch after it is prefetched (review 2 of line DBID
/// pfoff-2, L1).
#[derive(Debug, Default)]
pub struct RunLoads {
    run: Option<RunThrow>,
    ended: std::cell::Cell<bool>,
}

impl RunLoads {
    /// The record of a run whose reports do `run` after a throw, before any report.
    pub fn new(run: Option<RunThrow>) -> Self {
        RunLoads {
            run,
            ended: std::cell::Cell::new(false),
        }
    }

    /// What the run's reports do after a subject whose lookup throws (`None`: a report refuses
    /// it, or none reads one by one).
    pub fn run(&self) -> Option<RunThrow> {
        self.run
    }

    /// Whether a report's load ended the run (NCBI's application ends there).
    pub fn ended(&self) -> bool {
        self.ended.get()
    }

    /// [`HitLoads::of`] of the reports `reports` after those of the run before (none after its
    /// end), which records where these reports end the run.
    pub fn loads<'p, K: Into<ReportKeys>>(
        &self,
        reports: impl IntoIterator<Item = (bool, K)>,
        parsed: impl Fn(usize) -> Option<&'p ParsedSubject>,
        identities: impl Fn(usize, usize) -> usize,
    ) -> HitLoads {
        HitLoads::after(&self.ended, self.run, reports, parsed, identities)
    }
}

// No NCBI counterpart: LOSAT's check that no two sequences of a database are joined
// (`common::database::reports_may_refuse`).
/// Adds to `keys` the keys of the Seq-align id of `subject` (`defline_set::match_keys`): two
/// subjects that the reports join (`joined_in_reports`) share one.
pub fn join_keys(subject: &ParsedSubject, keys: &mut Vec<u64>) {
    crate::seqid::defline_set::match_keys(&subject.alignment_id, keys);
}

/// The keys of the Seq-align id of `subject` as NCBI's scope finds it
/// (`defline_set::scope_keys`).
pub fn scope_keys(subject: &ParsedSubject, keys: &mut Vec<u64>) {
    crate::seqid::defline_set::scope_keys(&subject.alignment_id, keys);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seqdb::{SeqDb, SeqType};
    use crate::seqid::model::ObjectId;
    use crate::seqid::model::{Dbtag, PdbSeqId};

    /// The subjects of `tests/db_search/ids/prot_v4` (`make_ids.py`), by the FASTA form of the
    /// first id of their first definition line.
    fn subject(first: &str) -> ParsedSubject {
        let base =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/ids/prot_v4");
        let db = SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap();
        for ordinal in 0..db.num_seqs() {
            let set = db.defline_set(ordinal).unwrap();
            let mut fasta = Vec::new();
            label::fasta(&set.deflines[0].ids[0], &mut fasta);
            if fasta == first.as_bytes() {
                return ParsedSubject::new(set).unwrap();
            }
        }
        panic!("no sequence {first}");
    }

    fn text(bytes: &[u8]) -> &str {
        std::str::from_utf8(bytes).unwrap()
    }

    // Approved exception 2 of PD-LOSAT-NCBI-DEFECTS: BLASTN, TBLASTN and TBLASTX stop an outfmt 0
    // title made only of punctuation at the end of the string (NCBI's `x_CleanAndCompress` reads
    // past it and crashes); a nucleotide subject with parsed Seq-ids gets the title of the same
    // stored title without them, and so does the tabular title word. A protein subject is refused.
    #[test]
    fn a_nucleotide_title_of_punctuation_stops_at_its_end() {
        use crate::algorithm::common::database::check_shown_parsed_subject;
        let subject = |id: ParsedSeqId| {
            let defline = BlastDefline {
                title: Some(b", ,".to_vec()),
                ids: vec![id],
                ..BlastDefline::default()
            };
            ParsedSubject::new(DeflineSet {
                deflines: vec![defline],
            })
            .unwrap()
        };
        let stored = crate::report::defline::generate_defline(b", ,", false, true).unwrap();
        assert!(stored.reads_past_end);
        let accession = subject(ParsedSeqId::Genbank(crate::seqid::model::TextseqId {
            name: None,
            accession: Some(b"AF000050".to_vec()),
            release: None,
            version: Some(1),
        }));
        let (start, title) = accession.description(false);
        assert_eq!(text(&start), "AF000050.1 ");
        assert_eq!(title.text, stored.text);
        check_shown_parsed_subject(
            &accession,
            false,
            false,
            || crate::seqdb::Oid::of_unrestricted(0),
            "BLASTN",
        )
        .expect("a nucleotide title of punctuation is shown");
        assert!(check_shown_parsed_subject(
            &accession,
            true,
            false,
            || crate::seqdb::Oid::of_unrestricted(0),
            "BLASTP"
        )
        .is_err());
        let local = subject(ParsedSeqId::Local(ObjectId::Str(b"Subject_51".to_vec())));
        assert!(local.tabular_ids(true).is_err());
        assert!(local.tabular_saccver(true).is_err());
        assert!(local.check_tabular_ids(true).is_err());
        assert!(local.check_tabular_ids(false).is_ok());
        let [id, _, _] = local.tabular_ids(false).unwrap().unwrap();
        let word = stored.text.split(|&byte| byte == b' ').next().unwrap();
        assert_eq!(id, word);
    }

    // NCBI BLAST+ 2.17.0 `blastn -db tests/db_search/ids/nucl_v5 -outfmt 17` names the subjects
    // so (`LOSAT/tests/fixtures/db_search/db.id.blastn.17.v5.out`, ordinal order): the Seq-align
    // id `gi|654321` and `gi|2020` (a GI in the second definition line) give the accession of
    // the Bioseq's definition line, `AB000009` has no version, and `lcl|Subject_31` no Textseq-id.
    #[test]
    fn sam_names_the_subject_by_an_accession_of_its_bioseq() {
        let base =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/ids/nucl_v5");
        let db = SeqDb::open(base.to_str().unwrap(), SeqType::Nucleotide).unwrap();
        let names: Vec<String> = (0..db.num_seqs())
            .map(|ordinal| {
                let subject = ParsedSubject::new(db.defline_set(ordinal).unwrap()).unwrap();
                let mut name = Vec::new();
                crate::report::seqid::SeqId::Db(subject.sam_id()).sam_name(&mut name);
                String::from_utf8(name).unwrap()
            })
            .collect();
        assert_eq!(
            names,
            [
                "NM_000001.1",
                "AF000002.3",
                "AJ000003.1",
                "AB000004.1",
                "NC_000005.11",
                "mynuc6",
                "mydb:nuc7",
                "NM_000008.2",
                "AB000009",
                "mynuc_10",
                "BK000012.1",
                "AF000013.1",
                "US1234567_14",
                "AF000015.1",
                "Subject_30",
                "Subject_31",
                "BK000018.1",
                "NC_000019.1",
                "AF000020.1",
                "AF000021.1",
                "AF000022.1",
                "NR_000023.1",
                "ti:123",
                "BN000025.1",
            ]
        );
    }

    fn description(first: &str) -> String {
        let (mut line, title) = subject(first).description(true);
        line.extend_from_slice(&title.text);
        String::from_utf8(line).unwrap()
    }

    // The strings are NCBI BLAST+ 2.17.0's `blastp -db` (outfmt 0 and 6) for these records
    // (`LOSAT/tests/fixtures/db_search/db.id.blastp.0.v4.out`, `db.id.blastp.6c.v4.out`).
    #[test]
    fn titles_that_ncbi_makes_from_the_seq_ids() {
        assert_eq!(
            description("tpg|DAA10001.1|"),
            "DAA10001.1 tpa titled protein"
        );
        assert_eq!(
            description("tpg|DAA10002.1|"),
            "DAA10002.1 unnamed protein product"
        );
        assert_eq!(
            text(&subject("tpg|DAA10002.1|").heading(true, false).unwrap()),
            ">DAA10002.1 TPA: unnamed protein product\n"
        );
        assert_eq!(description("pdb|2XYZ|B"), "2XYZ_B Chain B");
        assert_eq!(description("pdb|2XYY| "), "2XYY unnamed protein product");
        assert_eq!(
            description("pat|US|7654321|3"),
            "US7654321_3 Sequence 3 from Patent US 7654321"
        );
        assert_eq!(description("ref|NC_100013.1|"), "NC_100013.1 ");
        assert_eq!(
            text(&subject("ref|NC_100013.1|").heading(true, false).unwrap()),
            ">NC_100013.1 \n"
        );
        assert_eq!(
            description("ref|WP_100014.1|"),
            "WP_100014.1 unnamed protein product"
        );
        assert_eq!(description("lcl|Subject_27"), "subject seven title");
        assert_eq!(
            text(&subject("lcl|Subject_30").heading(true, false).unwrap()),
            ">Subject_30 unnamed protein product\n"
        );
        assert_eq!(description("pdb|5ABC|D"), "5ABC_D Chain D");
        assert_eq!(
            text(&subject("pdb|5ABC|D").heading(true, false).unwrap()),
            ">5ABC_D Chain D\n XP_100015.1 second\n"
        );
        assert_eq!(
            text(&subject("lcl|Subject_36").heading(true, false).unwrap()),
            ">Subject_36 unnamed protein product\n DAA10016.1 tpa second\n"
        );
        // A title from the source: the taxonomy database can name its organism.
        assert!(!subject("ref|NP_100006.1|").title_from_taxonomy());
        assert!(!subject("pdb|2XYZ|B").title_from_taxonomy());
        assert!(subject("gi|1818").title_from_taxonomy());
    }

    #[test]
    fn tabular_ids_replace_subject_ids_by_the_title_word() {
        let ids = |first: &str| {
            subject(first)
                .tabular_ids(true)
                .unwrap()
                .unwrap()
                .map(|id| String::from_utf8(id).unwrap())
        };
        assert_eq!(ids("lcl|Subject_36"), ["unnamed", "unnamed", "unnamed"]);
        assert_eq!(ids("lcl|Subject_27"), ["subject", "subject", "subject"]);
        assert_eq!(
            ids("lcl|Subject_28"),
            ["tpg|DAA10008.1|", "DAA10008", "DAA10008.1"]
        );
        assert_eq!(ids("lcl|Subject_29"), ["pdb|3ABC|C", "3ABC_C", "3ABC_C"]);
        assert!(subject("lcl|Subject_29").check_tabular_ids(true).is_ok());
        // The `saccver`, made once.
        let once = subject("lcl|Subject_28");
        let first = once.tabular_saccver(true).unwrap().unwrap();
        assert_eq!(&*first, b"DAA10008.1");
        assert!(std::sync::Arc::ptr_eq(
            &first,
            &once.tabular_saccver(true).unwrap().unwrap()
        ));
    }

    fn one_defline(ids: Vec<ParsedSeqId>) -> DeflineSet {
        DeflineSet {
            deflines: vec![BlastDefline {
                title: Some(b"t".to_vec()),
                ids,
                ..BlastDefline::default()
            }],
        }
    }

    fn gb(accession: &str, version: i32) -> ParsedSeqId {
        ParsedSeqId::Genbank(crate::seqid::model::TextseqId {
            accession: Some(accession.as_bytes().to_vec()),
            version: Some(version),
            ..Default::default()
        })
    }

    fn lines(lines: Vec<Vec<ParsedSeqId>>) -> ParsedSubject {
        let deflines = lines
            .into_iter()
            .map(|ids| BlastDefline {
                title: Some(b"t".to_vec()),
                ids,
                ..BlastDefline::default()
            })
            .collect();
        ParsedSubject::new(DeflineSet { deflines }).unwrap()
    }

    fn strings(subject: &ParsedSubject) -> Vec<String> {
        subject
            .use_this_seqid()
            .iter()
            .map(|string| String::from_utf8(string.clone()).unwrap())
            .collect()
    }

    // Review 3 of line RG T3a, D1 (NCBI BLAST+ 2.17.0 probes `T/t3/review3_t3a`, `T/t3/fix3_t3a`):
    // the description line shows the display id, the `WorstRank` best id of the line that holds
    // the Seq-align id (`gi|777|lcl|Query_1`: the Seq-align id `gi|777`, the display id
    // `lcl|Query_1`; `ref|NP_000777.1|lcl|Query_1`: `ref|NP_000777.1|` for both), and looks the
    // title up again by it: its keys and those of the ids it matches weakly are `display`
    // (`ref|NP_000781|` finds a query `ref|NP_000781.1|`, `ref|NP_000782.1|` not one
    // `ref|NP_000782|`), apart from the Seq-align id's `own` and `weak`.
    #[test]
    fn the_display_id_has_its_own_scope_keys() {
        let local = |id: &str| ParsedSeqId::Local(ObjectId::Str(id.as_bytes().to_vec()));
        let refseq = |accession: &str, version: Option<i32>| {
            ParsedSeqId::Other(crate::seqid::model::TextseqId {
                accession: Some(accession.as_bytes().to_vec()),
                version,
                ..Default::default()
            })
        };
        let fasta = |id: Option<&ParsedSeqId>| {
            let mut out = Vec::new();
            label::fasta(id.unwrap(), &mut out);
            String::from_utf8(out).unwrap()
        };
        let keys = |id: &ParsedSeqId| {
            let mut keys = Vec::new();
            crate::seqid::defline_set::scope_keys(id, &mut keys);
            keys
        };
        let query = local("Query_1");
        let gi = lines(vec![vec![ParsedSeqId::Gi(777), local("Query_1")]]);
        assert_eq!(text(&gi.alignment_fasta()), "gi|777");
        assert_eq!(fasta(gi.display_id()), "lcl|Query_1");
        let scope = gi.query_scope_keys();
        assert!(keys(&query).iter().all(|key| scope.display.contains(key)));
        assert!(!keys(&query).iter().any(|key| scope.own.contains(key)));
        let refq = lines(vec![vec![refseq("NP_000777", Some(1)), local("Query_1")]]);
        assert_eq!(text(&refq.alignment_fasta()), "ref|NP_000777.1|");
        assert_eq!(fasta(refq.display_id()), "ref|NP_000777.1|");
        assert!(!keys(&query)
            .iter()
            .any(|key| refq.query_scope_keys().display.contains(key)));
        let versioned_query = |accession: &str, version: Option<i32>| {
            crate::seqid::id_scope::version_keys(&refseq(accession, version))
                .map(|(own, _)| own)
                .unwrap()
        };
        let weak = lines(vec![vec![ParsedSeqId::Gi(781), refseq("NP_000781", None)]]);
        assert_eq!(fasta(weak.display_id()), "ref|NP_000781|");
        assert!(weak
            .query_scope_keys()
            .display
            .contains(&versioned_query("NP_000781", Some(1))));
        let versioned = lines(vec![vec![
            ParsedSeqId::Gi(782),
            refseq("NP_000782", Some(1)),
        ]]);
        assert!(!versioned
            .query_scope_keys()
            .display
            .contains(&versioned_query("NP_000782", None)));
    }

    // Line DBID N1b: the titles whose bytes the taxids of the record can change, whatever the
    // record's taxids are (`title_reads_source`).
    #[test]
    fn a_title_reads_the_source_where_it_is_made_from_it_or_its_end_is_rewritten() {
        let titled = |titles: &[&str]| {
            let deflines = titles
                .iter()
                .enumerate()
                .map(|(line, title)| BlastDefline {
                    title: Some(title.as_bytes().to_vec()),
                    ids: vec![gb(&format!("AF00000{line}"), 1)],
                    ..BlastDefline::default()
                })
                .collect();
            ParsedSubject::new(DeflineSet { deflines }).unwrap()
        };
        // A title made from the source names its organism with the taxonomy database only.
        let empty = titled(&[""]);
        assert!(empty.title_reads_source(true, true));
        assert!(empty.title_reads_source(false, true));
        assert!(!empty.title_reads_source(true, false));
        // A protein title ending in ']' where the heading generates the title (another line
        // without a title), with or without the taxonomy database; not for a nucleotide one.
        let bracket = titled(&["alpha, partial [Homo sapiens]"]);
        assert!(!bracket.title_reads_source(true, true));
        let generated = titled(&["alpha, partial [Homo sapiens]", ""]);
        assert!(generated.title_reads_source(true, false));
        assert!(!generated.title_reads_source(false, false));
        assert!(!titled(&["alpha", ""]).title_reads_source(true, true));
    }

    // Probes G01, G05 of line DI-8 (`GetFilteredRedundantSeqids`, `x_FillAlnDispParams`).
    #[test]
    fn a_positive_list_names_the_text_ids_or_the_gis() {
        let local = || ParsedSeqId::Local(ObjectId::Str(b"loc".to_vec()));
        // Text ids of every line, whichever comes first; a local id is not named.
        let subject = lines(vec![vec![local()], vec![gb("AF000001", 1)]]);
        assert!(subject.use_this_seqid().is_empty());
        let subject = subject.under_positive_list(true);
        assert_eq!(strings(&subject), ["seqid:AF000001.1"]);
        // The heading leaves out the line the strings do not name, and `>` goes to the next.
        assert_eq!(subject.heading(true, false).unwrap(), b">AF000001.1 t\n");
        // An ungapped Seq-align of two HSPs reaches the report without the strings.
        assert_eq!(
            subject.heading(true, true).unwrap(),
            b">loc t\n AF000001.1 t\n"
        );
        // A local id alone: no string, every line in the heading.
        let subject = lines(vec![vec![local()]]).under_positive_list(true);
        assert!(subject.use_this_seqid().is_empty());
        assert_eq!(subject.heading(true, false).unwrap(), b">loc t\n");
        // A Seq-align id that is a GI: only GIs are named, and a line without one is left out.
        let subject = lines(vec![
            vec![ParsedSeqId::Gi(888), gb("GG000001", 1)],
            vec![gb("GG000002", 1)],
            vec![ParsedSeqId::Gi(890), gb("GG000003", 1)],
        ])
        .under_positive_list(true);
        assert_eq!(strings(&subject), ["gi:888", "gi:890"]);
        assert_eq!(
            subject.heading(true, false).unwrap(),
            b">GG000001.1 t\n GG000003.1 t\n"
        );
        // A name without an accession is not named; a PDB id is.
        let name_only = ParsedSeqId::Prf(crate::seqid::model::TextseqId {
            name: Some(b"NAME".to_vec()),
            ..Default::default()
        });
        let pdb = ParsedSeqId::Pdb(crate::seqid::model::PdbSeqId {
            mol: b"3DEF".to_vec(),
            chain: None,
            rel: None,
            chain_id: Some(b"B".to_vec()),
        });
        let subject = lines(vec![vec![name_only], vec![pdb]]).under_positive_list(true);
        assert_eq!(strings(&subject), ["seqid:3DEF_B"]);
        // Without a positive list nothing changes.
        let subject =
            lines(vec![vec![local()], vec![gb("AF000001", 1)]]).under_positive_list(false);
        assert!(subject.use_this_seqid().is_empty());
        assert_eq!(
            subject.heading(true, false).unwrap(),
            b">loc t\n AF000001.1 t\n"
        );
    }

    #[test]
    fn the_kind_of_a_seq_align_id_refuses_nothing_before_its_lookup() {
        // NCBI's lookup decides whether its formatters find these (`found_instead`; probe of line
        // DBID N2b), not the kind of id.
        for id in [
            ParsedSeqId::Local(ObjectId::Id(123)),
            ParsedSeqId::Gibbsq(123),
            ParsedSeqId::NamedAnnotTrack(crate::seqid::model::TextseqId {
                accession: Some(b"NA000001".to_vec()),
                version: Some(1),
                ..Default::default()
            }),
            ParsedSeqId::General(Dbtag {
                db: b"PIG".to_vec(),
                tag: ObjectId::Id(123),
            }),
            ParsedSeqId::General(Dbtag {
                db: b"BL_ORD_ID".to_vec(),
                tag: ObjectId::Id(5),
            }),
        ] {
            let subject = ParsedSubject::new(one_defline(vec![id])).unwrap();
            assert!(!subject.unfetchable(), "{subject:?}");
            assert!(subject.fetched(false, true, HitLoad::Loaded).is_ok());
        }
        // A local id that alone decides the definition line of the Bioseq: no lookup finds it.
        let shown = lines(vec![
            vec![ParsedSeqId::Local(ObjectId::Id(777))],
            vec![ParsedSeqId::Gi(777)],
        ]);
        assert!(shown.unfetchable());
        let error = shown.fetched(false, true, HitLoad::Loaded).unwrap_err();
        assert!(
            error.contains("Seq-align id gi|777, by which NCBI BLAST+'s formatters do not find"),
            "{error}"
        );
        // A GI in the set makes the GI the Seq-align id.
        let with_gi = one_defline(vec![ParsedSeqId::Gibbsq(123), ParsedSeqId::Gi(7)]);
        assert!(ParsedSubject::new(with_gi)
            .unwrap()
            .fetched(false, true, HitLoad::Loaded)
            .is_ok());
        let replaced = DeflineSet {
            deflines: vec![BlastDefline {
                ids: vec![ParsedSeqId::Gi(7)],
                replaced_title_bytes: 1,
                ..BlastDefline::default()
            }],
        };
        assert!(ParsedSubject::new(replaced).is_err());
        assert!(ParsedSubject::new(DeflineSet::default()).is_err());
    }

    #[test]
    fn a_seq_align_id_by_which_ncbi_finds_another_sequence_or_none() {
        let subject = ParsedSubject::new(one_defline(vec![gb("AAA00001", 0)])).unwrap();
        assert!(subject
            .clone()
            .found_instead(None)
            .fetched(false, true, HitLoad::Loaded)
            .is_ok());
        let error = subject
            .clone()
            .found_instead(Some(FoundInstead::Sequence(DbSequence::of_unrestricted(4))))
            .fetched(false, true, HitLoad::Loaded)
            .unwrap_err();
        assert!(
            error.contains("lookup gives database sequence 4 first"),
            "{error}"
        );
        let error = subject
            .clone()
            .found_instead(Some(FoundInstead::Nothing))
            .fetched(false, true, HitLoad::Loaded)
            .unwrap_err();
        assert!(error.contains("do not find the sequence"), "{error}");
        let error = subject
            .found_instead(Some(FoundInstead::Failed(LookupFailure::ncbi(
                "DB contains no accession info.",
            ))))
            .fetched(false, true, HitLoad::Loaded)
            .unwrap_err();
        assert!(
            error.contains("whose lookup in the database fails (DB contains no accession info.;"),
            "{error}"
        );
    }

    // Line DBID pfoff-1 (DECISIONS 62): where a batch's prefetch is off, the XML2 and JSON `Hit`
    // loads the subject by its Seq-align id and the lookup throws NCBI's `CSeqDBException` (blastp
    // 2.17.0, outfmt 15 of `pP pU`, `R/mixdb2/review1/out/ncbi/bp.lim1.U2.15.err`); only the text
    // seen there gives that exception (review 1, M2).
    #[test]
    fn only_the_lookup_text_ncbi_was_seen_to_throw_gives_its_exception() {
        let failed = LookupFailure::ncbi(crate::seqdb::tax_files::NO_ACCESSION_INFO);
        let exception = failed.exception().unwrap();
        assert_eq!(exception.get_msg(), "DB contains no accession info.");
        assert_eq!(
            exception.what(),
            format!(
                "NCBI C++ Exception:\n    T0 \"{}/src/objtools/blast/seqdb_reader/seqdb_lmdb.cpp\", line 142: Error: (CSeqDBException::eInvalid) BLASTDB::ncbi::CBlastLMDBManager::CBlastEnv::GetDbi() - DB contains no accession info.\n",
                crate::ncbi_exception::NCBI_BUILD_ROOT
            )
        );
        // LOSAT's own error, NCBI's text of a place that LOSAT does not know (the text of
        // `CATCH_ALL` for `eTooManyOpenFiles`), and the texts of NCBI's reader whose places LOSAT
        // knows but that no lookup was seen to throw (`seqdb::throw_site`): no exception.
        for failure in [
            LookupFailure::losat("byte range 0..9 passes the end of the file"),
            LookupFailure::ncbi("Too many open files, please raise the open file limit"),
            LookupFailure::ncbi("Error: File (x.pos) not found."),
            LookupFailure::ncbi("Database name is required."),
            LookupFailure::ncbi(crate::seqdb::tax_files::NO_VOLUME_INFO),
        ] {
            assert!(failure.exception().is_none(), "{failure:?}");
        }
    }

    // Line DBID pfoff-1 (DECISIONS 62; review 1, H1 and M1): the one answer of what the XML2 and
    // JSON reports load (`HitLoads`): in a report read one by one, the `Hit`s before the first
    // whose lookup throws load their subjects, that one throws, and those after it load nothing
    // (NCBI's report function ends the report there); a subject is loaded where some report
    // loads it; a report whose batch prefetches loads every subject; a run that writes another
    // format loads every subject (which refuses the throw). The refusals take that answer.
    #[test]
    fn the_hits_after_one_whose_lookup_throws_load_nothing() {
        let throws = || {
            Some(FoundInstead::Failed(LookupFailure::ncbi(
                crate::seqdb::tax_files::NO_ACCESSION_INFO,
            )))
        };
        let plain =
            |accession: &str| ParsedSubject::new(one_defline(vec![gb(accession, 1)])).unwrap();
        // 0 and 4 throw; 2 and 3 are copies of one record (one `Hit`); 1, 5 and 6 are found.
        let copy = plain("AAA00003");
        let parsed: Vec<Option<ParsedSubject>> = vec![
            Some(plain("AAA00001").found_instead(throws())),
            Some(plain("AAA00002")),
            Some(copy.clone()),
            Some(copy),
            Some(plain("AAA00005").found_instead(throws())),
            Some(plain("AAA00006")),
            None,
        ];
        assert!(joins_as_copy(
            parsed[2].as_ref().unwrap(),
            parsed[3].as_ref().unwrap()
        ));
        let subject = |key: usize| parsed[key].as_ref();
        let none = |_: usize, _: usize| 0;
        let loads = HitLoads::of(
            Some(RunThrow::every_row(AfterThrow::EndsReport)),
            [
                // Read one by one: 1 loads, 0 throws, the copies 2 and 3, 4 and 6 load nothing.
                (
                    true,
                    vec![Some(1), None, Some(0), Some(2), Some(3), Some(4), Some(6)],
                ),
                // 6 loads in a later report; 5 only in a batch that prefetches.
                (true, vec![Some(6)]),
                (false, vec![Some(5)]),
            ],
            subject,
            none,
        );
        assert_eq!(loads.load(1), HitLoad::Loaded);
        assert_eq!(loads.load(0), HitLoad::Throws);
        assert_eq!(
            loads
                .thrown(0)
                .map(crate::ncbi_exception::NcbiException::get_msg),
            Some("DB contains no accession info.")
        );
        for unread in [2, 3, 4] {
            assert_eq!(loads.load(unread), HitLoad::Unread, "{unread}");
            assert!(loads.thrown(unread).is_none());
        }
        assert!(loads.loaded(5) && loads.loaded(6));
        // A key the reports never show reads as loaded; renamed keys keep their answers.
        assert!(loads.loaded(7));
        let renamed = loads.renamed(|key| (key != 1).then_some(key + 10));
        assert_eq!(renamed.load(10), HitLoad::Throws);
        assert_eq!(renamed.load(14), HitLoad::Unread);
        assert_eq!(renamed.load(0), HitLoad::Loaded);
        // A report whose batch prefetches, and a run that writes another format, load every
        // subject.
        for (after, one_by_one) in [
            (Some(RunThrow::every_row(AfterThrow::EndsReport)), false),
            (None, true),
        ] {
            let every = HitLoads::of(after, [(one_by_one, vec![Some(0), Some(1)])], subject, none);
            assert!(every.loaded(0) && every.loaded(1));
        }
        // The refusals take the answer: a `Hit` that throws is not refused for what the lookup
        // finds, a subject that no report loads not at all, every other one as before (a lookup
        // that fails on a file too).
        let thrown = subject(0).unwrap();
        assert!(thrown.fetched(true, true, HitLoad::Throws).is_ok());
        assert!(thrown.fetched(true, true, HitLoad::Unread).is_ok());
        let error = thrown.fetched(false, true, HitLoad::Loaded).unwrap_err();
        assert!(
            error.contains("whose lookup in the database fails (DB contains no accession info.;"),
            "{error}"
        );
        let nothing = plain("AAA00007").found_instead(Some(FoundInstead::Nothing));
        assert!(nothing.fetched(false, true, HitLoad::Loaded).is_err());
        assert!(nothing.fetched(false, true, HitLoad::Unread).is_ok());
        // A weak match is refused where the `Hit` throws too (NCBI's scope may find it first).
        let weak = thrown
            .clone()
            .weak_match(Some(WeakMatch::Held(DbSequence::of_unrestricted(5))));
        assert!(weak.fetched(false, true, HitLoad::Throws).is_err());
    }

    // Line DBID pfoff-2 (DECISIONS 62 (12), probe `R/pfoff/p2`): what the tables load after a
    // subject whose lookup throws (`AfterThrow`, `FormatThrow`): a table of a field of the subject
    // block writes each row and goes on, every throwing subject throws and every other loads; a
    // table that reads the subject after that block (`slen`, the letters) and outfmt 13 and 14
    // end the run, so nothing after the first throw loads, in that report or a later one; a
    // table of other fields reads nothing; the run's answer is its formats' (`FormatThrow::and`).
    #[test]
    fn the_formats_say_what_loads_after_a_throw() {
        use crate::blastinput::app::Program;
        let throws = || {
            Some(FoundInstead::Failed(LookupFailure::ncbi(
                crate::seqdb::tax_files::NO_ACCESSION_INFO,
            )))
        };
        let plain =
            |accession: &str| ParsedSubject::new(one_defline(vec![gb(accession, 1)])).unwrap();
        // 1 and 3 throw; 0, 2 and 4 are found.
        let parsed: Vec<ParsedSubject> = vec![
            plain("AAA00001"),
            plain("AAA00002").found_instead(throws()),
            plain("AAA00003"),
            plain("AAA00004").found_instead(throws()),
            plain("AAA00005"),
        ];
        let subject = |key: usize| parsed.get(key);
        let reports = || {
            [
                (true, vec![Some(0), Some(1), Some(2)]),
                (true, vec![Some(3), Some(4)]),
            ]
        };
        let every = |after| Some(RunThrow::every_row(after));
        let none = |_: usize, _: usize| 0;
        let rows = HitLoads::of(every(AfterThrow::WritesRow), reports(), subject, none);
        assert!(rows.loaded(0) && rows.loaded(2) && rows.loaded(4));
        assert_eq!(rows.load(1), HitLoad::Throws);
        assert_eq!(rows.load(3), HitLoad::Throws);
        assert!(rows.any_throws());
        let run = HitLoads::of(every(AfterThrow::EndsRun), reports(), subject, none);
        assert!(run.loaded(0));
        assert_eq!(run.load(1), HitLoad::Throws);
        // The run ended at 1: the next report reads nothing (3 does not throw).
        for unread in [2, 3, 4] {
            assert_eq!(run.load(unread), HitLoad::Unread, "{unread}");
        }
        let report = HitLoads::of(every(AfterThrow::EndsReport), reports(), subject, none);
        assert_eq!(report.load(2), HitLoad::Unread);
        assert_eq!(report.load(3), HitLoad::Throws);
        assert_eq!(report.load(4), HitLoad::Unread);
        // The formats' answers.
        let of = |outfmt: &str| FormatThrow::of_outfmt(outfmt, Program::Blastp);
        let writes = |after| FormatThrow::Writes(RunThrow::every_row(after));
        for rows in [
            "6",
            "7 qseqid sseqid evalue",
            "10 qseqid stitle slen",
            "20",
            "20 qseqid sacc",
            "18",
        ] {
            assert_eq!(of(rows), writes(AfterThrow::WritesRow), "{rows}");
        }
        assert_eq!(of("6 qseqid slen"), writes(AfterThrow::EndsRun));
        assert_eq!(of("6 qseqid pident"), writes(AfterThrow::EndsRun));
        // `nident` alone reads the subject only at a row with identities.
        assert_eq!(
            of("6 qseqid nident"),
            FormatThrow::Writes(RunThrow {
                after: AfterThrow::EndsRun,
                rows: super::super::tabular::SubjectFetch::WithIdentities,
            })
        );
        assert_eq!(of("6 qseqid evalue bitscore"), FormatThrow::Unread);
        assert_eq!(of("13"), writes(AfterThrow::EndsRun));
        assert_eq!(of("15"), writes(AfterThrow::EndsReport));
        assert_eq!(of("8"), FormatThrow::Unread);
        for refused in ["5", "17"] {
            assert_eq!(of(refused), FormatThrow::Refused, "{refused}");
        }
        // Pairwise: the first shown alignment ends the run, a later one is skipped; query-anchored:
        // any shown alignment ends the run.
        assert_eq!(of("0"), writes(AfterThrow::SkipsAlignment));
        for anchored in ["1", "2", "3", "4"] {
            assert_eq!(
                of(anchored),
                writes(AfterThrow::EndsAtAlignment),
                "{anchored}"
            );
        }
        // BLASTN's `pident` and `nident` read nothing of the subject.
        for unread in ["6 qseqid pident", "6 qseqid nident"] {
            assert_eq!(
                FormatThrow::of_outfmt(unread, Program::Blastn),
                FormatThrow::Unread,
                "{unread}"
            );
        }
        // A run's answer: a refusal refuses, a format that reads nothing leaves the other's.
        let after =
            |outfmts: &[&str]| FormatThrow::of_outfmts(outfmts.iter().copied(), Program::Blastp);
        assert_eq!(
            after(&["6", "8"]),
            Some(RunThrow::every_row(AfterThrow::WritesRow))
        );
        // A pairwise format with another whose throws end elsewhere: refused.
        assert_eq!(after(&["6", "0"]), None);
        assert_eq!(after(&["0", "3"]), None);
        assert_eq!(
            after(&["0", "8"]),
            Some(RunThrow::every_row(AfterThrow::SkipsAlignment))
        );
        assert_eq!(after(&["6 qseqid evalue"]), None);
        assert_eq!(
            after(&["6 qseqid slen", "15"]),
            Some(RunThrow::every_row(AfterThrow::EndsReport))
        );
        // Rows of `nident` alone and rows of every HSP load at every row.
        assert_eq!(
            after(&["6 qseqid nident", "6 qseqid slen"]),
            Some(RunThrow::every_row(AfterThrow::EndsRun))
        );
    }

    // Review 2 of line DBID pfoff-2, M1 and L1 (probe `R/pfoff/p2/review2`, `zP`): where `nident`
    // alone reads the subject (`SubjectFetch::WithIdentities`), a throwing subject whose rows have
    // no identities reads nothing (its rows print the 0 of `x_ResetFields` and the run goes on);
    // the run ends at the first throwing subject with a row of identities, and nothing after it
    // loads, in that report, a later one, or one that a later call makes apart
    // (`RunLoads::loads`).
    #[test]
    fn a_throw_ends_the_run_only_at_a_row_that_reads_the_subject() {
        let throws = || {
            Some(FoundInstead::Failed(LookupFailure::ncbi(
                crate::seqdb::tax_files::NO_ACCESSION_INFO,
            )))
        };
        let plain =
            |accession: &str| ParsedSubject::new(one_defline(vec![gb(accession, 1)])).unwrap();
        // 1 and 3 throw; 0, 2 and 4 are found.
        let parsed: Vec<ParsedSubject> = vec![
            plain("AAA00001"),
            plain("AAA00002").found_instead(throws()),
            plain("AAA00003"),
            plain("AAA00004").found_instead(throws()),
            plain("AAA00005"),
        ];
        let subject = |key: usize| parsed.get(key);
        let nident = Some(RunThrow {
            after: AfterThrow::EndsRun,
            rows: super::super::tabular::SubjectFetch::WithIdentities,
        });
        // The identities of each report's entries: 1's rows have none, 3's have some.
        let identities = |report: usize, entry: usize| [[5, 0, 7], [2, 4, 0]][report][entry];
        let reports = || {
            [
                (true, vec![Some(0), Some(1), Some(2)]),
                (true, vec![Some(3), Some(4)]),
            ]
        };
        let loads = HitLoads::of(nident, reports(), subject, identities);
        assert_eq!(loads.load(1), HitLoad::Unread);
        assert!(loads.thrown(1).is_none());
        assert!(loads.loaded(0) && loads.loaded(2));
        assert_eq!(loads.load(3), HitLoad::Throws);
        assert_eq!(loads.load(4), HitLoad::Unread);
        // Every row reads the subject: 1 throws and ends the run.
        let every = HitLoads::of(
            Some(RunThrow::every_row(AfterThrow::EndsRun)),
            reports(),
            subject,
            identities,
        );
        assert_eq!(every.load(1), HitLoad::Throws);
        assert_eq!(every.load(3), HitLoad::Unread);
        // The reports made apart, one call each, carry the end of the run.
        let run = RunLoads::new(nident);
        let [first, second] = reports();
        let first = run.loads([first], subject, |_, entry| [5, 0, 7][entry]);
        assert!(!run.ended());
        assert_eq!(first.load(1), HitLoad::Unread);
        let second = run.loads([second], subject, |_, entry| [2, 4][entry]);
        assert!(run.ended());
        assert_eq!(second.load(3), HitLoad::Throws);
        let later = run.loads([(true, vec![Some(0), Some(2)])], subject, |_, _| 9);
        assert_eq!(later.load(0), HitLoad::Unread);
        assert_eq!(later.load(2), HitLoad::Unread);
        // A run without a throw's answer loads every subject and never ends.
        let none = RunLoads::new(None);
        assert!(none.loads(reports(), subject, identities).loaded(3));
        assert!(!none.ended());
    }

    // Review 1 of line DBID pfoff-2, H1 and M1 (probe `R/pfoff/p2/review1`, `*.H.*`, `*.H2.*`):
    // the tables' taxonomy, from which their names are looked up, is made only of the subjects
    // that a report loads (`common::tabular::tabular_taxonomy` with `HitLoads`): a throwing
    // subject's taxid (9606 here), which NCBI never reads, opens no taxonomy database and is not
    // looked up, so its record is never refused and the texts of opening wait for a loaded row.
    #[test]
    fn the_names_of_the_tables_take_only_the_loaded_subjects() {
        let throws = || {
            Some(FoundInstead::Failed(LookupFailure::ncbi(
                crate::seqdb::tax_files::NO_ACCESSION_INFO,
            )))
        };
        let with_taxid = |accession: &str, taxid: i32| {
            let mut set = one_defline(vec![gb(accession, 1)]);
            set.deflines[0].taxid = Some(taxid);
            ParsedSubject::new(set).unwrap()
        };
        // 0 throws (taxid 9606); 1 is found (taxid 99990).
        let parsed: Vec<Option<std::sync::Arc<ParsedSubject>>> = vec![
            Some(std::sync::Arc::new(
                with_taxid("NM_000001", 9606).found_instead(throws()),
            )),
            Some(std::sync::Arc::new(with_taxid("XM_000002", 99990))),
        ];
        let loads = HitLoads::of(
            Some(RunThrow::every_row(AfterThrow::WritesRow)),
            [(true, vec![Some(0), Some(1)])],
            |rank| parsed[rank].as_deref(),
            |_, _| 0,
        );
        let base =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/db_search/ids/prot_v4");
        let db = SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap();
        let taxonomy = crate::algorithm::common::tabular::tabular_taxonomy(
            &db,
            &[0, 1],
            &parsed,
            [0, 1],
            &loads,
            "BLASTP",
        )
        .unwrap();
        assert_eq!(
            taxonomy[0],
            super::super::tabular::TabularTaxonomy::default()
        );
        assert_eq!(taxonomy[1].taxid, 99990);
        // Every subject loaded (no report reads one by one): both taxids.
        let every = crate::algorithm::common::tabular::tabular_taxonomy(
            &db,
            &[0, 1],
            &parsed,
            [0, 1],
            &HitLoads::default(),
            "BLASTP",
        )
        .unwrap();
        assert_eq!(every[0].taxid, 9606);
        // The names of the rows: no lookup asks a taxid of a subject that is not loaded, so the
        // damaged taxonomy database (`bad_magic`, whose open writes a text) is not opened for it.
        let spec = super::super::tabular::TabularSpec::of(
            crate::blastinput::app::ReportFormat::Tabular,
            &crate::blastinput::app::FormatChoice {
                number: 6,
                spec: "qseqid sseqid ssciname".to_string(),
                delimiter: String::new(),
            },
        )
        .unwrap();
        let damaged =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxdb/bad_magic");
        let names = |taxonomy: &super::super::tabular::TabularTaxonomy| {
            let search =
                crate::seqdb::SearchPath::new(&damaged, Some(damaged.as_os_str()), None).unwrap();
            super::super::tabular::TabularNames::look_up(
                &spec,
                &crate::seqdb::taxdb::LazyTaxDb::new(search),
                [taxonomy],
                "BLASTP",
            )
            .unwrap()
        };
        assert!(!names(&taxonomy[0]).looks_up());
        assert!(names(&every[0]).looks_up());
    }

    // Line DBID pssmread (probe `R/pssmread`, reviews 1 and 2): LOSAT's refusals of what
    // psiblast's PSSM reads by the Seq-align id (what the scope gets is
    // `common::database::scope_bioseq`'s).
    #[test]
    fn the_pssm_refuses_reads_that_losat_does_not_make_as_ncbi() {
        let subject = ParsedSubject::new(one_defline(vec![gb("AAA00001", 0)])).unwrap();
        let read = |subject: &ParsedSubject| subject.pssm_read(false);
        assert_eq!(read(&subject), Ok(()));
        assert!(subject.bioseq_holds_alignment_id());
        assert_eq!(subject.other_holder(), None);
        let other = subject
            .clone()
            .found_instead(Some(FoundInstead::Sequence(DbSequence::of_unrestricted(4))));
        assert_eq!(read(&other), Ok(()));
        for found in [
            FoundInstead::Nothing,
            FoundInstead::ThrowsOnId("Error: ... is not a valid ordinal id".into()),
            FoundInstead::OtherTaxids(DbSequence::of_unrestricted(2)),
        ] {
            assert_eq!(read(&subject.clone().found_instead(Some(found))), Ok(()));
        }
        // Refused: another holder besides what the lookup gives, a lookup failing on a file, the
        // subject first with another holder.
        let held = other
            .clone()
            .also_held(Some(DbSequence::of_unrestricted(6)));
        assert_eq!(held.other_holder(), Some(DbSequence::of_unrestricted(6)));
        assert!(read(&held)
            .unwrap_err()
            .contains("which database sequence 6 also holds"));
        let none_held = subject
            .clone()
            .found_instead(Some(FoundInstead::Nothing))
            .also_held(Some(DbSequence::of_unrestricted(6)));
        assert!(read(&none_held).is_err());
        for found in [
            FoundInstead::AlsoHeld(DbSequence::of_unrestricted(6)),
            FoundInstead::Failed(LookupFailure::ncbi("DB contains no accession info.")),
        ] {
            assert!(read(&subject.clone().found_instead(Some(found))).is_err());
        }
        // A weak match refuses where the batch reads one by one, whatever the lookup gives.
        let weak = other
            .clone()
            .weak_match(Some(WeakMatch::Held(DbSequence::of_unrestricted(7))));
        assert_eq!(weak.pssm_read(false), Ok(()));
        assert!(weak.pssm_read(true).is_err());
        // A definition line that the formatters do not find: its Bioseq does not hold the id, and
        // the reasons are what the lookup finds (review 2, M2).
        let shown = lines(vec![
            vec![ParsedSeqId::Local(ObjectId::Id(777))],
            vec![ParsedSeqId::Gi(777)],
        ]);
        assert!(shown.unfetchable());
        assert!(!shown.bioseq_holds_alignment_id());
        assert_eq!(read(&shown), Ok(()));
        let also = shown
            .clone()
            .found_instead(Some(FoundInstead::AlsoHeld(DbSequence::of_unrestricted(2))));
        assert_eq!(also.other_holder(), Some(DbSequence::of_unrestricted(2)));
        let error = read(&also).unwrap_err();
        assert!(
            error.contains("which database sequence 2 also holds"),
            "{error}"
        );
        assert!(!error.contains("formatters do not find"), "{error}");
        let failed = shown.found_instead(Some(FoundInstead::Failed(LookupFailure::ncbi(
            "no accession info",
        ))));
        assert!(read(&failed)
            .unwrap_err()
            .contains("whose lookup in the database fails"));
    }

    // Line DBID weakbatch (probe w0): NCBI's scope matches a Seq-align id weakly only where the
    // formatters read the batch's subjects one by one, its prefetch off; the prefetch loads each
    // subject by its own id.
    #[test]
    fn a_weak_match_is_refused_only_where_the_batch_reads_one_by_one() {
        let subject = ParsedSubject::new(one_defline(vec![gb("AAA00001", 0)])).unwrap();
        for weak in [
            WeakMatch::Held(DbSequence::of_unrestricted(4)),
            WeakMatch::Unlisted,
            WeakMatch::Failed("no accession info".into()),
        ] {
            let subject = subject.clone().weak_match(Some(weak.clone()));
            assert!(
                subject.fetched(true, false, HitLoad::Loaded).is_ok(),
                "{weak:?}"
            );
            let error = subject.fetched(false, true, HitLoad::Loaded).unwrap_err();
            assert!(
                error.contains("in a batch whose prefetch PRE_FETCH_SEQS_LIMIT turns off"),
                "{error}"
            );
            assert_eq!(subject.prefetched(), Prefetched::Itself, "{weak:?}");
            // What the lookup finds is refused whatever the batch: a copy with other taxids in
            // the reports of taxids, before the weak match.
            let copy = subject.found_instead(Some(FoundInstead::OtherTaxids(
                DbSequence::of_unrestricted(2),
            )));
            assert!(copy.fetched(false, false, HitLoad::Loaded).is_ok());
            assert!(copy
                .fetched(true, false, HitLoad::Loaded)
                .unwrap_err()
                .contains("other taxids"));
            assert!(copy
                .fetched(false, true, HitLoad::Loaded)
                .unwrap_err()
                .contains("weakly"));
        }
        let other =
            subject.found_instead(Some(FoundInstead::AlsoHeld(DbSequence::of_unrestricted(4))));
        assert!(other.fetched(false, false, HitLoad::Loaded).is_err());
    }

    // Line DBID P1: what NCBI's prefetch reads for a subject that no report shows (psiblast's
    // traceback hit list beyond -max_target_seqs), from what its lookup finds.
    #[test]
    fn the_prefetch_of_an_unshown_subject_reads_what_the_lookup_finds() {
        let subject = ParsedSubject::new(one_defline(vec![gb("AAA00001", 0)])).unwrap();
        let prefetched =
            |found: Option<FoundInstead>| subject.clone().found_instead(found).prefetched();
        for itself in [
            None,
            Some(FoundInstead::OtherTaxids(DbSequence::of_unrestricted(2))),
            Some(FoundInstead::CopyApart(
                DbSequence::of_unrestricted(2),
                CopyApart::Title,
            )),
        ] {
            assert_eq!(prefetched(itself.clone()), Prefetched::Itself, "{itself:?}");
        }
        assert_eq!(prefetched(Some(FoundInstead::Nothing)), Prefetched::Nothing);
        // A file's error stops NCBI's prefetch; the exception of the id is its warning.
        assert!(matches!(
            prefetched(Some(FoundInstead::Failed(LookupFailure::ncbi("no accession info")))),
            Prefetched::Failed(reason) if reason.contains("prefetch stops with a BLAST Database error")
        ));
        let throws = Some(FoundInstead::ThrowsOnId("Invalid choice selection".into()));
        assert!(matches!(
            prefetched(throws.clone()),
            Prefetched::Failed(reason) if reason.contains("writes Error pre-fetching sequence data")
        ));
        assert!(subject.clone().found_instead(throws).prefetch_fails());
        assert!(!subject
            .clone()
            .found_instead(Some(FoundInstead::Failed(LookupFailure::ncbi(
                "no accession info"
            ))))
            .prefetch_fails());
        assert!(matches!(
            prefetched(Some(FoundInstead::Sequence(DbSequence::of_unrestricted(4)))),
            Prefetched::Another {
                other: Some(other),
                or_itself: false,
                ..
            } if other == DbSequence::of_unrestricted(4)
        ));
        assert!(matches!(
            prefetched(Some(FoundInstead::AlsoHeld(DbSequence::of_unrestricted(4)))),
            Prefetched::Another {
                other: Some(other),
                or_itself: true,
                ..
            } if other == DbSequence::of_unrestricted(4)
        ));
        assert!(matches!(
            prefetched(Some(FoundInstead::CopyApart(
                DbSequence::of_unrestricted(4),
                CopyApart::TaxonomyTexts
            ))),
            Prefetched::Another { other: None, .. }
        ));
        // Other holders that LOSAT cannot list (line DBID v4late): what the lookup gives, or
        // nothing where such a holder's Bioseq read before stands for it.
        assert!(matches!(
            prefetched(Some(FoundInstead::HoldersUnknown {
                instead: Some(DbSequence::of_unrestricted(4)),
                why: String::new(),
            })),
            Prefetched::Another {
                other: Some(other),
                or_itself: false,
                or_none: true,
                ..
            } if other == DbSequence::of_unrestricted(4)
        ));
        assert!(matches!(
            prefetched(Some(FoundInstead::HoldersUnknown {
                instead: None,
                why: String::new(),
            })),
            Prefetched::ItselfOrNothing { .. }
        ));
        // A subject whose definition line the formatters do not find: what its lookup gives, as
        // any other (review 2 of line DBID pssmread, H3).
        let shown = lines(vec![
            vec![ParsedSeqId::Local(ObjectId::Id(777))],
            vec![ParsedSeqId::Gi(777)],
        ]);
        assert_eq!(shown.prefetched(), Prefetched::Itself);
        assert!(matches!(
            shown.clone().found_instead(Some(FoundInstead::Sequence(DbSequence::of_unrestricted(3)))).prefetched(),
            Prefetched::Another { other: Some(other), or_itself: false, or_none: false, reason }
                if other == DbSequence::of_unrestricted(3) && reason.contains("gives database sequence 3 first")
        ));
        assert_eq!(
            shown
                .found_instead(Some(FoundInstead::Nothing))
                .prefetched(),
            Prefetched::Nothing
        );
    }

    // Line DBID N1b L1: a copy with other taxids from which NCBI can make the subject's title is
    // refused by the reports that show a generated title (outfmt 0, the tabular word of a
    // replaced id) and by those of taxids, not by the others (tabular ids, `stitle`, XML).
    #[test]
    fn a_title_made_from_a_copy_s_taxids_is_refused_where_a_generated_title_is_shown() {
        use crate::algorithm::common::database::check_shown_parsed_subject;
        let title_apart = Some(FoundInstead::CopyApart(
            DbSequence::of_unrestricted(3),
            CopyApart::Title,
        ));
        let subject = ParsedSubject::new(one_defline(vec![gb("AAA00001", 1)]))
            .unwrap()
            .found_instead(title_apart.clone());
        assert!(subject.fetched(false, true, HitLoad::Loaded).is_ok());
        let error = subject.fetched(true, true, HitLoad::Loaded).unwrap_err();
        assert!(error.contains("database sequence 3 also holds"), "{error}");
        assert!(error.contains("taxids and their names for it"), "{error}");
        assert!(subject.check_tabular_ids(true).is_ok());
        assert!(subject.all_tabular_ids(true).is_ok());
        // Outfmt 0 where it shows a title that can depend on the source: an empty title with
        // the taxonomy database, a heading generated of a title ending in `]`.
        assert!(check_shown_parsed_subject(
            &subject,
            true,
            true,
            || crate::seqdb::Oid::of_unrestricted(0),
            "BLASTP"
        )
        .is_ok());
        let titled = |titles: &[&str]| {
            let deflines = titles
                .iter()
                .enumerate()
                .map(|(line, title)| BlastDefline {
                    title: Some(title.as_bytes().to_vec()),
                    ids: vec![gb(&format!("AF00000{line}"), 1)],
                    ..BlastDefline::default()
                })
                .collect();
            ParsedSubject::new(DeflineSet { deflines })
                .unwrap()
                .found_instead(title_apart.clone())
        };
        let error = check_shown_parsed_subject(
            &titled(&[""]),
            false,
            true,
            || crate::seqdb::Oid::of_unrestricted(0),
            "BLASTN",
        )
        .unwrap_err()
        .to_string();
        assert!(error.contains("can make that sequence's title"), "{error}");
        assert!(check_shown_parsed_subject(
            &titled(&[""]),
            false,
            false,
            || crate::seqdb::Oid::of_unrestricted(0),
            "BLASTN"
        )
        .is_ok());
        let heading = titled(&["alpha, partial [Homo sapiens]", ""]);
        assert!(check_shown_parsed_subject(
            &heading,
            true,
            false,
            || crate::seqdb::Oid::of_unrestricted(0),
            "BLASTP"
        )
        .is_err());
        let bracket = titled(&["alpha, partial [Homo sapiens]"]);
        assert!(check_shown_parsed_subject(
            &bracket,
            true,
            true,
            || crate::seqdb::Oid::of_unrestricted(0),
            "BLASTP"
        )
        .is_ok());
        // A local id `Subject_N` that the title's first word replaces, of the target or of
        // another definition line: the word can depend on the source for a nucleotide and for a
        // protein title ending in `]`; refused where a printed id shows the word.
        use crate::report::tabular::PrintedIds;
        let (subject_ids, all_ids) = (
            PrintedIds {
                subject: true,
                all: false,
            },
            PrintedIds {
                subject: false,
                all: true,
            },
        );
        let local = || ParsedSeqId::Local(ObjectId::Str(b"Subject_1".to_vec()));
        let replaced = ParsedSubject::new(one_defline(vec![local()]))
            .unwrap()
            .found_instead(title_apart.clone());
        assert!(replaced.check_tabular_ids(false).is_ok());
        let error = replaced.check_printed_word(false, subject_ids).unwrap_err();
        assert!(
            error.contains("whose first word replaces a printed id"),
            "{error}"
        );
        assert!(replaced.check_printed_word(false, all_ids).is_err());
        assert!(replaced
            .check_printed_word(false, PrintedIds::default())
            .is_ok());
        assert!(replaced.check_printed_word(true, subject_ids).is_ok());
        let second = ParsedSubject::new(DeflineSet {
            deflines: vec![
                one_defline(vec![gb("AAA00001", 1)]).deflines.remove(0),
                one_defline(vec![local()]).deflines.remove(0),
            ],
        })
        .unwrap()
        .found_instead(title_apart.clone());
        assert!(second.check_printed_word(false, subject_ids).is_ok());
        assert!(second.check_printed_word(false, all_ids).is_err());
        assert!(second.check_printed_word(true, all_ids).is_ok());
        // A replaced id beside a better one is not printed (`GetSeqIdListString`).
        let beside = ParsedSubject::new(one_defline(vec![gb("AAA00001", 1), local()]))
            .unwrap()
            .found_instead(title_apart.clone());
        assert!(beside.check_printed_word(false, subject_ids).is_ok());
        assert!(beside.check_printed_word(false, all_ids).is_ok());
        let bracket = ParsedSubject::new(DeflineSet {
            deflines: vec![BlastDefline {
                title: Some(b"alpha, partial [Homo sapiens]".to_vec()),
                ids: vec![local()],
                ..BlastDefline::default()
            }],
        })
        .unwrap()
        .found_instead(title_apart);
        assert!(bracket.check_printed_word(true, subject_ids).is_err());
        // The other copies are refused as before: the taxonomy texts by every report that reads
        // the Bioseq, the taxids by those that write them.
        let texts = ParsedSubject::new(one_defline(vec![gb("AAA00001", 1)]))
            .unwrap()
            .found_instead(Some(FoundInstead::CopyApart(
                DbSequence::of_unrestricted(3),
                CopyApart::TaxonomyTexts,
            )));
        assert!(texts.fetched(false, true, HitLoad::Loaded).is_err());
        assert!(texts.generated_title_fetched().is_ok());
        let taxids = ParsedSubject::new(one_defline(vec![gb("AAA00001", 1)]))
            .unwrap()
            .found_instead(Some(FoundInstead::OtherTaxids(
                DbSequence::of_unrestricted(3),
            )));
        assert!(taxids.fetched(false, true, HitLoad::Loaded).is_ok());
        assert!(taxids.generated_title_fetched().is_ok());
    }

    // Line DBID N1b L1 probe: of a nucleotide `lcl|Subject_1` without a title and with taxid 9606,
    // NCBI's tabular id is the first word of the title it makes from the source where the
    // taxonomy database names the taxid (`Homo`); a protein's is `unnamed` either way (probe of
    // line RG T5, `T/t5/out/ncbi/s.*.ids`: `Homo`, `Uncultured`, `unnamed`). A source that
    // LOSAT refuses (a lookup that crashes NCBI, a named RefSeq record) refuses the word.
    #[test]
    fn a_tabular_word_of_a_title_made_from_the_source_names_the_organism() {
        use crate::report::tabular::PrintedIds;
        let subject = |first: &str, taxid: Option<i32>, source: SourceTaxname| {
            ParsedSubject::new(DeflineSet {
                deflines: vec![BlastDefline {
                    ids: vec![match first.strip_prefix("ref|") {
                        Some(accession) => ParsedSeqId::Other(crate::seqid::model::TextseqId {
                            accession: Some(accession.as_bytes().to_vec()),
                            version: Some(1),
                            ..Default::default()
                        }),
                        None => ParsedSeqId::Local(ObjectId::Str(first.as_bytes().to_vec())),
                    }],
                    taxid,
                    ..BlastDefline::default()
                }],
            })
            .unwrap()
            .with_source(true, |taxid| {
                assert_eq!(taxid, 9606);
                source
            })
        };
        let named = || SourceTaxname::Named(std::sync::Arc::from(&b"uncultured bacterium"[..]));
        let printed = PrintedIds::of(None);
        let local = subject("Subject_1", Some(9606), named());
        assert!(local.check_printed_word(false, printed).is_ok());
        let [id, _, _] = local.tabular_ids(false).unwrap().unwrap();
        assert_eq!(text(&id), "Uncultured");
        let [id, _, _] = local.tabular_ids(true).unwrap().unwrap();
        assert_eq!(text(&id), "unnamed");
        assert_eq!(
            text(&local.bioseq_title(true, true).unwrap().text),
            "unnamed protein product [uncultured bacterium]"
        );
        assert_eq!(
            text(&local.bioseq_title(false, false).unwrap().text),
            "Uncultured bacterium"
        );
        let unnamed = subject("Subject_1", Some(9606), SourceTaxname::Unnamed);
        assert!(unnamed.tabular_ids(false).is_ok());
        assert_eq!(text(&unnamed.bioseq_title(false, false).unwrap().text), "");
        let failed = subject(
            "Subject_1",
            Some(9606),
            SourceTaxname::Refused("a failed lookup"),
        );
        assert_eq!(failed.tabular_ids(false).unwrap_err(), "a failed lookup");
        assert_eq!(failed.tabular_ids(true).unwrap_err(), "a failed lookup");
        // Without a taxid the Bioseq has no source: nothing is looked up.
        let none = ParsedSubject::new(DeflineSet {
            deflines: vec![BlastDefline {
                ids: vec![ParsedSeqId::Local(ObjectId::Str(b"Subject_1".to_vec()))],
                ..BlastDefline::default()
            }],
        })
        .unwrap()
        .with_source(true, |_| panic!("no source to name"));
        assert_eq!(none.source(), &SourceTaxname::Unnamed);
        // A named RefSeq chromosome gets the name (`x_SetTitleFromNC` makes nothing of a database
        // Bioseq), a protein's not upper-cased and not `unnamed protein product`; an unnamed one
        // keeps its empty title.
        let refseq = subject("ref|NC_000001", Some(9606), named());
        assert!(matches!(refseq.source(), SourceTaxname::Named(_)));
        assert_eq!(
            text(&refseq.bioseq_title(false, false).unwrap().text),
            "Uncultured bacterium"
        );
        assert_eq!(
            text(&refseq.bioseq_title(true, false).unwrap().text),
            "uncultured bacterium"
        );
        let refseq = subject("ref|NC_000001", Some(9606), SourceTaxname::Unnamed);
        assert_eq!(refseq.source(), &SourceTaxname::Unnamed);
    }

    // NCBI writes `Sequence with id AB000010 no longer exists in database...alignment skipped` in
    // place of each alignment of a sequence whose heading's generated title it cannot decode (probe
    // of line RG T5, `pt.100009.*.f0`); a title it can decode, or a stored one, is shown.
    #[test]
    fn a_heading_that_does_not_decode_names_the_seq_align_id() {
        let subject = |title: Option<&[u8]>, name: &[u8]| {
            ParsedSubject::new(DeflineSet {
                deflines: vec![BlastDefline {
                    ids: vec![gb("AB000010", 1)],
                    title: title.map(<[u8]>::to_vec),
                    taxid: Some(100009),
                    ..BlastDefline::default()
                }],
            })
            .unwrap()
            .with_source(false, |_| SourceTaxname::Named(std::sync::Arc::from(name)))
        };
        let bad = subject(None, b"Bad \x81 name");
        assert_eq!(
            bad.heading_skipped_id(false, false).as_deref(),
            Some(&b"AB000010"[..])
        );
        // The heading's own making of the title gives the same id (one GenerateDefline).
        assert_eq!(
            bad.heading(false, false).unwrap_err(),
            bad.heading_skipped_id(false, false).unwrap()
        );
        assert_eq!(text(&bad.description(false).1.text), "Unknown");
        assert_eq!(
            subject(None, b"Homo sapiens").heading_skipped_id(false, false),
            None
        );
        assert_eq!(
            text(
                &subject(None, b"Homo sapiens")
                    .heading(false, false)
                    .unwrap()
            ),
            ">AB000010.1 Homo sapiens\n"
        );
        assert_eq!(
            subject(Some(b"stored"), b"Bad \x81 name").heading_skipped_id(true, false),
            None
        );
        assert!(subject(Some(b"stored"), b"Bad \x81 name")
            .heading(true, false)
            .is_ok());
    }

    #[test]
    fn pdb_ids_without_a_chain_match_any_chain_of_the_molecule() {
        let pdb = |chain_id: Option<&str>| {
            let id = ParsedSeqId::Pdb(PdbSeqId {
                mol: b"1ABC".to_vec(),
                chain: chain_id.map(|id| i32::from(id.as_bytes()[0])),
                rel: None,
                chain_id: chain_id.map(|id| id.as_bytes().to_vec()),
            });
            ParsedSubject::new(one_defline(vec![id])).unwrap()
        };
        assert!(joined_in_reports(&pdb(None), &pdb(Some("AA"))));
        assert!(!joined_in_reports(&pdb(Some("A")), &pdb(Some("B"))));
    }

    // Line DBID N3: the copies of one record in the databases of an alias (three of r1, two of
    // r3, the b ones with taxid 9606 in pjoin_abt) are joined; r6 (r1's letters under another id)
    // and a Seq-align id that matches another one's without being its copy (pfound_v5's
    // gb|DDD00004| and gb|DDD00004.1|, pdb|2XYZ| and pdb|2XYZ|AA) are not.
    #[test]
    fn copies_of_one_record_are_joined_and_other_matching_ids_are_not() {
        let open = |name: &str| {
            let base = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/db_search/ids")
                .join(name);
            let db = SeqDb::open(base.to_str().unwrap(), SeqType::Protein).unwrap();
            (0..db.num_seqs())
                .map(|ordinal| ParsedSubject::new(db.defline_set(ordinal).unwrap()).unwrap())
                .collect::<Vec<_>>()
        };
        // a: r1 r2 r3, b: r1 r3 r4 r6, c: r1.
        for name in ["pjoin_abc", "pjoin_abt"] {
            let s = open(name);
            assert_eq!(s.len(), 8);
            for (earlier, later) in [(0, 3), (3, 0), (3, 7), (7, 0), (2, 4)] {
                assert!(
                    joins_as_copy(&s[earlier], &s[later]),
                    "{name} {earlier} {later}"
                );
            }
            for (earlier, later) in [(7, 6), (6, 3), (0, 1), (2, 5)] {
                assert!(
                    !joined_in_reports(&s[earlier], &s[later]),
                    "{name} {earlier} {later}"
                );
                assert!(
                    !joins_as_copy(&s[earlier], &s[later]),
                    "{name} {earlier} {later}"
                );
            }
        }
        let s = open("pfound_v5");
        for prefix in ["gb|DDD00004", "pdb|2XYZ|"] {
            let pair: Vec<&ParsedSubject> = s
                .iter()
                .filter(|subject| subject.alignment_fasta().starts_with(prefix.as_bytes()))
                .collect();
            assert_eq!(pair.len(), 2, "{prefix}");
            assert!(joined_in_reports(pair[0], pair[1]), "{prefix}");
            assert!(!joins_as_copy(pair[0], pair[1]), "{prefix}");
        }
        let parsed: Vec<Option<std::sync::Arc<ParsedSubject>>> = open("pjoin_abc")
            .into_iter()
            .map(|subject| Some(std::sync::Arc::new(subject)))
            .chain([None])
            .collect();
        assert!(rank_joins_as_copy(&parsed, 3, 0));
        assert!(!rank_joins_as_copy(&parsed, 0, 8));
        assert!(!rank_joins_as_copy(&parsed, 8, 0));
    }

    #[test]
    fn refseq_accessions_whose_empty_titles_come_from_the_source() {
        use RefseqKind::{ChromosomeOrNcrna, Mrna};
        for (accession, kind) in [
            ("NC_000001", Some(ChromosomeOrNcrna)),
            ("AC_000001", Some(ChromosomeOrNcrna)),
            ("NM_000001", Some(Mrna)),
            ("NM_000000001", Some(Mrna)),
            ("XM_000000001", Some(Mrna)),
            ("NR_000001", Some(ChromosomeOrNcrna)),
            ("NR_000000001", None),
            ("NC_0000001", None),
            ("NP_000001", None),
            ("WP_000000001", None),
            ("XR_000001", None),
        ] {
            assert_eq!(
                refseq_accession_kind(accession.as_bytes()),
                kind,
                "{accession}"
            );
        }
    }
    /// The subject of a `-subject` record read with `-parse_deflines` (`ParsedSubject::of_record`).
    fn record_subject(defline: &str, protein: bool) -> ParsedSubject {
        record_subject_of(defline.as_bytes(), protein)
    }

    /// `record_subject` of a defline of any bytes.
    fn record_subject_of(defline: &[u8], protein: bool) -> ParsedSubject {
        use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};
        let mut fasta = defline.to_vec();
        fasta.extend_from_slice(if protein { b"\nMKV\n" } else { b"\nACGT\n" });
        let config = ReaderConfig::subject("BLASTN", protein, false).with_parse_deflines(true);
        let mut source = FastaInputSource::from_bytes(&fasta, config);
        let record = read_all(&mut source, &mut |_| Ok(())).unwrap().remove(0);
        ParsedSubject::of_record(record.believed().unwrap(), &record.title)
    }

    // NCBI BLAST+ 2.17.0 -parse_deflines -subject (probe `T/t7/review` of line RG, `C.lclS9`,
    // `J.pS9tpa`): a parsed `lcl|Subject_9` keeps its id (`s_ReplaceLocalId` with `parse_local`),
    // and `GetSeqIdList` puts the first word of the title generated from the rest of the defline
    // in its place: `s7` of `TPA: s7 third ~`; where that title is in no encoding NCBI guesses the
    // lookup throws (`sseqid` Unknown, `sall*` the Seq-align id `Subject_9`), a record's
    // definition line is never looked up.
    #[test]
    fn parsed_subject_ids_of_a_title_ncbi_cannot_decode_fail_the_lookup() {
        let decodable = record_subject(">lcl|Subject_9 TPA: s7 third ~", false);
        assert_eq!(
            decodable.tabular_ids(false).unwrap().unwrap()[2],
            b"s7".to_vec()
        );
        let bad = record_subject_of(b">lcl|Subject_9 bad x\x81y title", false);
        assert_eq!(bad.tabular_ids(false).unwrap(), None);
        assert_eq!(bad.tabular_saccver(false).unwrap(), None);
        assert!(bad.check_tabular_ids(false).is_ok());
        assert_eq!(
            bad.lookup_failure(false).unwrap(),
            Some(crate::report::tabular::LookupFailure {
                subject_id: true,
                all_ids: false,
                bioseq: false,
            })
        );
        assert_eq!(
            bad.alignment_tabular_ids(),
            [b"Subject_9".to_vec(), b"Subject_9".to_vec()]
        );
        assert_eq!(bad.skipped_id(), b"Subject_9");
        let all = bad.all_tabular_ids(false).unwrap();
        assert!(!all.failed);
        assert_eq!(text(&all.ids[0].id), "Subject_9");
        // A record of a decodable title has no failure.
        assert_eq!(decodable.lookup_failure(false).unwrap(), None);
    }

    // NCBI BLAST+ 2.17.0 -parse_deflines -subject (probe `DI7/ws/pd_sub.py`, 2026-11-01): `sseqid`,
    // `sacc`, `saccver`, `sallseqid`, the outfmt 0 heading and description of a record.
    #[test]
    fn records_read_with_parse_deflines_follow_the_ncbi_oracle() {
        for (defline, protein, ids, all, heading, local) in [
            (
                ">NM_000003.1 x y",
                false,
                ["ref|NM_000003.1|", "NM_000003", "NM_000003.1"],
                "ref|NM_000003.1|",
                ">NM_000003.1 x y\n",
                false,
            ),
            (
                ">plain1 a",
                false,
                ["plain1", "plain1", "plain1"],
                "plain1",
                ">plain1 a\n",
                false,
            ),
            (
                ">tpg|BK000001.1| y",
                false,
                ["tpg|BK000001.1|", "BK000001", "BK000001.1"],
                "tpg|BK000001.1|",
                ">BK000001.1 TPA: y\n",
                false,
            ),
            (">s1", false, ["s1", "s1", "s1"], "s1", ">s1 \n", false),
            (
                ">lcl|Subject_1 x y",
                false,
                ["x", "x", "x"],
                "Subject_1",
                "",
                true,
            ),
            (
                ">",
                false,
                ["Subject_1", "Subject_1", "Subject_1"],
                "Subject_1",
                "",
                true,
            ),
            (
                ">",
                true,
                ["unnamed", "unnamed", "unnamed"],
                "Subject_1",
                "",
                true,
            ),
        ] {
            let subject = record_subject(defline, protein);
            let shown = subject.tabular_ids(protein).unwrap().unwrap();
            assert_eq!(
                shown.map(|id| String::from_utf8(id).unwrap()),
                ids,
                "{defline}"
            );
            let all_ids = subject.all_tabular_ids(protein).unwrap().ids;
            assert_eq!(text(&all_ids[0].id), all, "{defline}");
            assert_eq!(subject.shown_as_record(), local, "{defline}");
            if !local {
                assert_eq!(
                    text(&subject.heading(protein, false).unwrap()),
                    heading,
                    "{defline}"
                );
                assert!(subject.first_title().is_empty() || defline.contains(' '));
            }
        }
    }

    // NCBI BLAST+ 2.17.0 `blastp -num_threads 1` on the version 5 database `psfx` of line RG
    // T5c (`tests/db_search/t5/make_t5.py --suffix-only`; rows `db.t5.sfx.blastp.*`): the heading
    // of an untitled line, outfmt 18 and the tabular word of `lcl|Subject_9nn` show the target's
    // title with its end rewritten from the first source (`x_AdjustProteinTitleSuffix`); the
    // description (`fLeavePrefixSuffix`) keeps the stored end. Only those subjects ask for the
    // source's name.
    #[test]
    fn a_protein_title_ending_in_a_bracket_gets_its_end_from_the_source() {
        use crate::algorithm::common::database::check_shown_parsed_subject;
        let set = |lines: &[(ParsedSeqId, &str, Option<i32>)]| DeflineSet {
            deflines: lines
                .iter()
                .map(|(id, title, taxid)| BlastDefline {
                    ids: vec![id.clone()],
                    title: (!title.is_empty()).then(|| title.as_bytes().to_vec()),
                    taxid: *taxid,
                    ..BlastDefline::default()
                })
                .collect(),
        };
        let bacillus = || SourceTaxname::Named(std::sync::Arc::from(&b"Bacillus sp."[..]));
        let named = |lines: &[(ParsedSeqId, &str, Option<i32>)], protein| {
            ParsedSubject::new(set(lines))
                .unwrap()
                .with_source(protein, |taxid| {
                    assert_eq!(taxid, 100002);
                    bacillus()
                })
        };
        // The first source comes from the second line.
        let heading = named(
            &[
                (gb("AC400001", 1), "alpha, partial [Bacillus sp.]", None),
                (gb("AD400001", 1), "", Some(100002)),
            ],
            true,
        );
        assert_eq!(heading.source(), &bacillus());
        assert_eq!(
            text(&heading.heading(true, false).unwrap()),
            ">AC400001.1 alpha, partial [Bacillus sp.]\n AD400001.1 alpha [Bacillus sp.]\n"
        );
        assert_eq!(
            text(&heading.description(true).1.text),
            "alpha, partial [Bacillus sp.]"
        );
        assert!(check_shown_parsed_subject(
            &heading,
            true,
            true,
            || crate::seqdb::Oid::of_unrestricted(0),
            "BLASTP"
        )
        .is_ok());
        let deflines = heading.tax_deflines(true);
        assert!(matches!(
            &deflines[1].title,
            crate::report::tax_report::TaxTitle::Text(title) if title == b"alpha [Bacillus sp.]"
        ));
        // A nucleotide is not adjusted, and asks nothing.
        let nucleotide = ParsedSubject::new(set(&[
            (gb("AC400001", 1), "alpha, partial [Bacillus sp.]", None),
            (gb("AD400001", 1), "", Some(100002)),
        ]))
        .unwrap()
        .with_source(false, |_| panic!("a nucleotide title is not adjusted"));
        assert_eq!(
            text(&nucleotide.heading(false, false).unwrap()),
            ">AC400001.1 alpha, partial [Bacillus sp.]\n AD400001.1 alpha, partial [Bacillus sp.]\n"
        );
        // The tabular word of `lcl|Subject_9nn`.
        let local = |tag: &str| ParsedSeqId::Local(ObjectId::Str(tag.as_bytes().to_vec()));
        let word = named(
            &[(
                local("Subject_901"),
                "alpha, partial [Bacillus sp.]",
                Some(100002),
            )],
            true,
        );
        let [id, accession, version] = word.tabular_ids(true).unwrap().unwrap();
        assert_eq!([text(&id), text(&accession), text(&version)], ["alpha"; 3]);
        // Without a name (no taxonomy database) a ` [ ]` end goes.
        let nameless = ParsedSubject::new(set(&[(local("Subject_903"), "a, [ ]", Some(77777))]))
            .unwrap()
            .with_source(true, |_| SourceTaxname::Unnamed);
        assert_eq!(text(&nameless.tabular_ids(true).unwrap().unwrap()[0]), "a");
        // A title whose generated form no report shows (no untitled line, no replaced id), a
        // title not ending in `]` and a Bioseq without a source ask nothing.
        for lines in [
            vec![(
                gb("AC400002", 1),
                "alpha, partial [Bacillus sp.]",
                Some(100002),
            )],
            vec![
                (gb("AC400003", 1), "alpha protein", Some(100002)),
                (gb("AD400003", 1), "", Some(100002)),
            ],
            vec![
                (gb("AC400031", 1), "gg, partial [Bacillus sp.]", None),
                (gb("AD400031", 1), "", None),
            ],
        ] {
            let subject = ParsedSubject::new(set(&lines))
                .unwrap()
                .with_source(true, |_| panic!("no name is asked for {lines:?}"));
            assert_eq!(subject.source(), &SourceTaxname::Unnamed);
        }
        // An adjusted title of punctuation is read past its end (NCBI BLAST+ 2.17.0 crashes:
        // probe of line RG T5c, `pedge_e3`): refused.
        let edge = ParsedSubject::new(set(&[
            (gb("AG400003", 1), ", , [ ]", Some(77777)),
            (gb("AH400003", 1), "", None),
        ]))
        .unwrap()
        .with_source(true, |_| SourceTaxname::Unnamed);
        let error = check_shown_parsed_subject(
            &edge,
            true,
            true,
            || crate::seqdb::Oid::of_unrestricted(0),
            "BLASTP",
        )
        .unwrap_err()
        .to_string();
        assert!(
            error.contains("x_CleanAndCompress reads past the end"),
            "{error}"
        );
    }

    // DBID pfoff-3 (DECISIONS 62 (13); review 2 of the item, L7): where a pairwise report ends the
    // run is decided on the alignments that it shows, after the sort and the cut: NCBI's probe
    // `R/pfoff/p3` `bp.U1.b1.dflt.0` ends (exit 2) where `bp.U1.b1.sh4.0`, whose `-sorthits 4` puts
    // a loaded subject first, writes skipped lines (exit 0), and `bp.U2.b1.dflt.3` (outfmt 3) ends
    // where `bp.U2.b1.dflt.0` does not; a local Seq-align id is not looked up there, and
    // `-num_alignments 0` shows none. The model of the loads (`HitLoads::of` of the reads that
    // `QueryGroup::read_order` makes before the rows) and the writer (`aligned_end` on its group
    // and the loads) end at the same alignment: the subjects aligned before it load, it throws,
    // those after it are never read (with `-num_descriptions 0`, no description line reads them).
    #[test]
    fn the_model_and_the_writer_end_the_run_at_the_same_alignment() {
        use crate::report::hit_sort::ReportSort;
        use crate::report::x_pairwise_par::{aligned_end, writer_reads, QueryGroup};
        let throws = || {
            Some(FoundInstead::Failed(LookupFailure::ncbi(
                crate::seqdb::tax_files::NO_ACCESSION_INFO,
            )))
        };
        let plain =
            |accession: &str| ParsedSubject::new(one_defline(vec![gb(accession, 1)])).unwrap();
        let local = ParsedSubject::new(one_defline(vec![ParsedSeqId::Local(ObjectId::Str(
            b"ENST1".to_vec(),
        ))]))
        .unwrap();
        // Subject 0 throws and comes first by E-value; 1 loads and comes first by bit score; 2
        // throws. With `local_first`, 0's Seq-align id is a local one.
        let parsed = |local_first: bool| {
            vec![
                if local_first {
                    local.clone().found_instead(throws())
                } else {
                    plain("AAA00001").found_instead(throws())
                },
                plain("AAA00002"),
                plain("AAA00003").found_instead(throws()),
            ]
        };
        let rows: Vec<crate::report::pairwise::PairwiseHit> =
            [(0, 30.0, 1e-20), (1, 90.0, 1e-10), (2, 20.0, 1e-5)]
                .into_iter()
                .map(|(s_idx, bit_score, e_value)| {
                    crate::report::pairwise::PairwiseHit::from(crate::common::Hit {
                        identity: 100.0,
                        length: 10,
                        mismatch: 0,
                        gapopen: 0,
                        q_start: 1,
                        q_end: 10,
                        s_start: 1,
                        s_end: 10,
                        e_value,
                        bit_score,
                        num_ident: 10,
                        query_frame: 1,
                        query_length: 0,
                        q_idx: 0,
                        s_idx,
                        raw_score: 25,
                        sort_query_offset: 0,
                        sort_query_end: 0,
                        sort_subject_offset: 0,
                        sort_subject_end: 0,
                        has_sort_offsets: false,
                        gap_info: None,
                        num_positives: 10,
                    })
                })
                .collect();
        let exception = LookupFailure::ncbi(crate::seqdb::tax_files::NO_ACCESSION_INFO)
            .exception()
            .unwrap();
        // Where the model and the writer end the run, as the position among the shown alignments.
        let end = |after: AfterThrow, hits: Option<u8>, alignments: usize, local_first: bool| {
            let parsed = parsed(local_first);
            let reading = crate::report::pairwise::PairwiseReading {
                sort: ReportSort { hits, hsps: None },
                shown: Shown {
                    descriptions: 0,
                    alignments,
                },
            };
            let reads = QueryGroup::read_order(
                rows.iter().map(|row| (row.hit.s_idx, row)),
                |row| crate::report::hit_sort::HspScores::of(row),
                |_, _| false,
                reading,
                false,
                None,
            );
            let run = RunLoads::new(Some(RunThrow::every_row(after)));
            let loads = run.loads(
                [(true, ReportKeys::Pairwise(reads))],
                |key| parsed.get(key),
                |_, _| 0,
            );
            let writer = writer_reads(&rows, |_, _| false, reading, false, None);
            let aligned: Vec<u32> = (0..writer.aligned)
                .map(|hit| writer.hits[hit][0] as u32)
                .collect();
            let writer_end = aligned_end(after, &aligned, |s_idx| {
                loads
                    .thrown(s_idx as usize)
                    .map(|_| (&exception, parsed[s_idx as usize].alignment_id_is_local()))
            })
            .map(|(at, _)| at);
            assert_eq!(run.ended(), writer_end.is_some());
            for (position, &s_idx) in aligned.iter().enumerate() {
                let expected = match writer_end {
                    Some(at) if position > at => HitLoad::Unread,
                    _ if parsed[s_idx as usize].lookup_throws().is_some() => HitLoad::Throws,
                    _ => HitLoad::Loaded,
                };
                assert_eq!(loads.load(s_idx as usize), expected, "{s_idx}");
            }
            writer_end
        };
        let every = usize::MAX;
        // Outfmt 0 ends at its first alignment where it throws, not after `-sorthits` moves a
        // loaded one first, nor with no alignment shown.
        assert_eq!(end(AfterThrow::SkipsAlignment, None, every, false), Some(0));
        assert_eq!(end(AfterThrow::SkipsAlignment, Some(1), every, false), None);
        assert_eq!(end(AfterThrow::SkipsAlignment, None, 0, false), None);
        assert_eq!(end(AfterThrow::SkipsAlignment, None, every, true), None);
        // Outfmt 1 to 4 end at the first shown alignment that throws, a local one too, unless
        // the cut leaves none.
        assert_eq!(
            end(AfterThrow::EndsAtAlignment, Some(1), every, false),
            Some(1)
        );
        assert_eq!(end(AfterThrow::EndsAtAlignment, None, every, true), Some(0));
        assert_eq!(end(AfterThrow::EndsAtAlignment, Some(1), 1, false), None);
        for other in [
            AfterThrow::EndsRun,
            AfterThrow::EndsReport,
            AfterThrow::WritesRow,
        ] {
            assert!(!other.is_pairwise());
            assert_eq!(other.alignment_end([0], |_| true, |_| false), None);
        }
    }

    // DBID pfoff-3: the description lines load each shown subject in order (one that throws
    // loads nothing), then the alignments load theirs up to the one that ends the run; after the
    // run's end nothing loads, in a later report either; a batch that prefetches loads all; a
    // subject that neither the description lines nor the alignments show is never read (review 1
    // of pfoff-3, H1: psiblast's later round shows the subjects of its first HSPs).
    #[test]
    fn a_pairwise_report_loads_its_descriptions_then_its_alignments() {
        let throws = || {
            Some(FoundInstead::Failed(LookupFailure::ncbi(
                crate::seqdb::tax_files::NO_ACCESSION_INFO,
            )))
        };
        let plain =
            |accession: &str| ParsedSubject::new(one_defline(vec![gb(accession, 1)])).unwrap();
        // 1 and 3 throw.
        let parsed: Vec<ParsedSubject> = vec![
            plain("AAA00001"),
            plain("AAA00002").found_instead(throws()),
            plain("AAA00003"),
            plain("AAA00004").found_instead(throws()),
            plain("AAA00005"),
            plain("AAA00006"),
        ];
        let subject = |key: usize| parsed.get(key);
        let none = |_: usize, _: usize| 0;
        let run = |after| Some(RunThrow::every_row(after));
        // A report of one subject per hit whose description lines show its first `descriptions`
        // hits and whose alignments its first `alignments`; a batch that prefetches gives its
        // subjects alone.
        let report = |one_by_one: bool, keys: &[usize], descriptions: usize, alignments: usize| {
            let keys = if one_by_one {
                ReportKeys::Pairwise(PairwiseReads {
                    hits: keys.iter().map(|&key| vec![key]).collect(),
                    described: (0..keys.len().min(descriptions)).collect(),
                    aligned: keys.len().min(alignments),
                })
            } else {
                ReportKeys::Entries(keys.iter().map(|&key| Some(key)).collect())
            };
            (one_by_one, keys)
        };
        let loads = |run, reports: Vec<(bool, ReportKeys)>| {
            let loads = HitLoads::of(run, reports, subject, none);
            (0..6).map(|key| loads.load(key)).collect::<Vec<_>>()
        };
        use HitLoad::{Loaded as L, Throws as T, Unread as U};
        let every = usize::MAX;
        let three = |descriptions, alignments, keys: [&[usize]; 3]| {
            keys.iter()
                .map(|keys| report(true, keys, descriptions, alignments))
                .collect::<Vec<_>>()
        };
        // Outfmt 0: the second report's first alignment throws; its description lines loaded 4
        // before; the third report loads nothing.
        let ends_second: [&[usize]; 3] = [&[0, 1, 2], &[3, 4], &[5]];
        assert_eq!(
            loads(
                run(AfterThrow::SkipsAlignment),
                three(every, every, ends_second)
            ),
            [L, T, L, T, L, U]
        );
        // A first alignment that loads: the throwing ones are skipped and the run goes on.
        let skips: [&[usize]; 3] = [&[0, 1, 2], &[4, 3], &[5]];
        assert_eq!(
            loads(run(AfterThrow::SkipsAlignment), three(every, every, skips)),
            [L, T, L, T, L, L]
        );
        // `-num_alignments 0`: no alignment section, nothing ends.
        assert_eq!(
            loads(
                run(AfterThrow::SkipsAlignment),
                three(every, 0, ends_second)
            ),
            [L, T, L, T, L, L]
        );
        // Outfmt 1 to 4: the first report's second alignment ends the run.
        assert_eq!(
            loads(run(AfterThrow::EndsAtAlignment), three(every, every, skips)),
            [L, T, L, U, U, U]
        );
        // `-num_alignments 1`: only the first alignment is shown, which loads; 3 throws in the
        // description lines alone.
        assert_eq!(
            loads(run(AfterThrow::EndsAtAlignment), three(every, 1, skips)),
            [L, T, L, T, L, L]
        );
        // `-num_descriptions 0`: only the alignments up to the end load, and 2 is never reached.
        assert_eq!(
            loads(run(AfterThrow::EndsAtAlignment), three(0, every, skips)),
            [L, T, U, U, U, U]
        );
        // The description lines that a report shows are its writer's (`PairwiseReads::described`,
        // in the order of the lines): psiblast's later round shows those of its first HSPs, here
        // 4 alone (the part of the subjects found before) then 2, so 0, 5, 1 and 3, which neither
        // the lines nor the alignments (`-num_alignments 0`) show, are never read (1 and 3, whose
        // lookups throw, throw nothing).
        let two_part = (
            true,
            ReportKeys::Pairwise(PairwiseReads {
                hits: vec![vec![0], vec![2], vec![4], vec![5], vec![1], vec![3]],
                described: vec![2, 1],
                aligned: 0,
            }),
        );
        assert_eq!(
            loads(run(AfterThrow::SkipsAlignment), vec![two_part]),
            [U, U, L, U, L, U]
        );
        // A batch that prefetches loads every subject it shows, where `-num_descriptions 0
        // -num_alignments 1` reads the first alignment alone of a batch that does not.
        assert_eq!(
            loads(
                run(AfterThrow::SkipsAlignment),
                vec![
                    report(false, &[0, 2, 4], 0, 1),
                    report(true, &[5, 1, 3], 0, 1)
                ]
            ),
            [L, U, L, U, L, L]
        );
        // After the run's end a batch that prefetches loads nothing either.
        assert_eq!(
            loads(
                run(AfterThrow::SkipsAlignment),
                vec![
                    report(true, &[1, 0], every, every),
                    report(false, &[2, 3, 4, 5], every, every)
                ]
            ),
            [L, T, U, U, U, U]
        );
        // A hit that shows copies of its record: a thrown first subject loads none of them (1 and
        // 2, in hits that are not shown, are not read).
        let copies = (
            true,
            ReportKeys::Pairwise(PairwiseReads {
                hits: vec![vec![4, 5], vec![3, 0], vec![1], vec![2]],
                described: vec![0, 1],
                aligned: 0,
            }),
        );
        assert_eq!(
            loads(run(AfterThrow::SkipsAlignment), vec![copies]),
            [U, U, U, T, L, L]
        );
        // A later call's reports read nothing after the end (`RunLoads`).
        let calls = RunLoads::new(run(AfterThrow::SkipsAlignment));
        let first = calls.loads([report(true, &[3, 4], every, every)], subject, none);
        assert_eq!(first.load(3), HitLoad::Throws);
        assert!(calls.ended());
        let later = calls.loads([report(true, &[5], every, every)], subject, none);
        assert_eq!(later.load(5), HitLoad::Unread);
    }

    // DBID pfoff-3, NCBI 2.17.0 (probe `R/pfoff/p3`, `bp.U2.b1.dflt.0`, `bp.U2.b1.gis.0`): the
    // description line of a subject whose load throws is the Seq-align id and `Unknown`
    // (`x_GetDeflineInfo`'s catch), the id's FASTA form with `-show_gis` (no space), and its
    // alignments' line names the accession without its version.
    #[test]
    fn a_thrown_description_is_the_alignment_id_and_unknown() {
        let refseq = ParsedSubject::new(one_defline(vec![ParsedSeqId::Other(
            crate::seqid::model::TextseqId {
                accession: Some(b"NP_100001".to_vec()),
                version: Some(1),
                ..Default::default()
            },
        )]))
        .unwrap();
        let line = |subject: &ParsedSubject| {
            let (mut start, title) = subject.thrown_description();
            start.extend_from_slice(&title.text);
            String::from_utf8(start).unwrap()
        };
        assert_eq!(line(&refseq), "NP_100001.1 Unknown");
        assert_eq!(
            line(&refseq.clone().showing_gis(true)),
            "ref|NP_100001.1|Unknown"
        );
        assert_eq!(refseq.skipped_id(), b"NP_100001");
        assert!(!refseq.alignment_id_is_local());
        // A GI Seq-align id is a reading of the model, not NCBI's bytes: no database makes a GI's
        // lookup throw (NCBI finds a GI by its ISAM index, not by the accession table whose
        // absence throws; review 1 of pfoff-3, L3), so no probe or fixture shows the line. What
        // can be tested is that it is made as the line of a loaded subject of that id
        // (`description_start`, which the probes check), then `Unknown`, after two spaces where
        // `-show_gis` wrote `gi|N|`.
        let gi = ParsedSubject::new(one_defline(vec![ParsedSeqId::Gi(1234)])).unwrap();
        for subject in [gi.clone(), gi.showing_gis(true)] {
            let (start, title) = subject.thrown_description();
            assert_eq!(start, subject.description_start());
            let gap: &[u8] = if subject.show_gis { b"  " } else { b"" };
            assert_eq!(title.text, [gap, b"Unknown"].concat());
        }
    }
}
