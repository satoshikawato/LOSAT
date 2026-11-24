//! The definition lines of one database sequence (`Blast-def-line-set`) and the choices that
//! every report makes from them: the Seq-align id (`alignment_id`), the definition line of the
//! subject Bioseq (`bioseq_defline`) and of the shown id (`display_defline`), and the joins of
//! the formats that show every definition line (the XML `Hit_def`, the lines of an outfmt 0
//! alignment block; the tabular `sall*`, `salltitles` and `staxids` are made by
//! `report::db_subject::ParsedSubject`, `report::tabular::TabularTaxonomy` and the tabular
//! writer from these definition lines).
//!
//! NCBI's formatters read the definition lines through the object manager: the Seq-align id
//! is resolved to the sequence again, whose Bioseq holds the ids and title of one definition
//! line (the first that holds the id) and the whole set as a user object
//! (`CSeqDB::ExtractBlastDefline`). LOSAT knows the sequence and keeps the set; these
//! functions give the parts NCBI takes from each.

use std::borrow::Cow;
use std::collections::BTreeSet;

use super::label::{bare, content, fasta, first_gi, format_label, number, seqid_list_string};
use super::model::{Dbtag, ObjectId, ParsedSeqId, PatentNumber, PatentSeqId, PdbSeqId, TextseqId};
use super::rank::{best_rank_score, find_best, text_score};

/// One `Blast-def-line`, every member decoded (`seqid::ber::decode_full_set`).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BlastDefline {
    /// `title`, after NCBI's repair of the stored bytes; `None` when the member is absent
    /// (makeblastdb writes an empty title instead).
    pub title: Option<Vec<u8>>,
    /// The number of stored title bytes outside ' ' ..= '~' that read as '#' (NCBI writes an
    /// error message for each time it reads the definition line).
    pub replaced_title_bytes: usize,
    /// `seqid`, in the stored order.
    pub ids: Vec<ParsedSeqId>,
    /// The number of bytes of the Seq-id strings that read as '#'.
    pub replaced_id_bytes: usize,
    /// `taxid` (NCBI's `TTaxId`, an `int`; `seqid::ber::Reader::int`).
    pub taxid: Option<i32>,
    pub memberships: Option<Vec<i32>>,
    /// `links`: further taxids of the definition line (`GetTaxIds`).
    pub links: Option<Vec<i32>>,
    pub other_info: Option<Vec<i32>>,
}

/// The `Blast-def-line-set` of one sequence.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DeflineSet {
    pub deflines: Vec<BlastDefline>,
    /// Whether the read of the stored record met a UTF8String where a VisibleString is
    /// expected (`seqid::ber::Reader::visible_string`): NCBI reads it, with its warning once per
    /// process at its first read of such a record (`seqid::ber::utf8_string_path`). A set that
    /// was not read from a record has none.
    pub utf8: bool,
}

impl From<Vec<BlastDefline>> for DeflineSet {
    fn from(deflines: Vec<BlastDefline>) -> Self {
        DeflineSet {
            deflines,
            utf8: false,
        }
    }
}

/// `CSeq_id::E_SIC`, the result of `CSeq_id::Compare`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sic {
    /// The ids are of types that cannot be compared (`e_DIFF`).
    Diff,
    /// Comparable ids that do not match (`e_NO`).
    No,
    /// The ids match (`e_YES`).
    Yes,
}

/// The definition line whose id an NCBI formatter shows (`display_defline`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DisplayDefline {
    /// The index of the definition line.
    pub index: usize,
    /// Whether it holds the Seq-align id. If none does, NCBI shows the id of the last
    /// definition line and no GI (`gi` is then 0).
    pub found: bool,
    /// The `taxid` of the last definition line up to `index` that has one, 0 for none
    /// (`GetDisplayIds`' `taxid`: it need not be the shown line's).
    pub taxid: i32,
}

/// The kind of id of a tabular field (`CBlastTabularInfo::ESeqIdType`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TabularId {
    /// `sseqid` (`eFullId`).
    Full,
    /// `sacc` (`eAccession`).
    Accession,
    /// `saccver` (`eAccVersion`).
    AccessionVersion,
}

impl BlastDefline {
    // NCBI reference (598d8ae6): c++/src/objects/blastdb/Blast_def_line.cpp:97-125
    /// `CBlast_def_line::GetTaxIds`: the taxid and the links, without 0 unless it is alone.
    pub fn tax_ids(&self) -> BTreeSet<i32> {
        let mut ids: BTreeSet<i32> = self.taxid.into_iter().collect();
        if let Some(links) = &self.links {
            ids.extend(links.iter().copied());
        }
        if ids.len() > 1 {
            ids.remove(&0);
        }
        ids
    }

    /// The first GI of the definition line, 0 for none (`FindGi`).
    pub fn gi(&self) -> i64 {
        first_gi(&self.ids)
    }

    /// The id of the definition line that NCBI's displays choose: the first `WorstRank`-best.
    pub fn display_id(&self) -> Option<&ParsedSeqId> {
        find_best(&self.ids, text_score)
    }
}

impl DeflineSet {
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:956-1019
    /// The title of the Bioseq: the first definition line's title, then ` >` with the FASTA
    /// strings of the ids of each later one (joined by `|`) and its title.
    /// Each line's part (NCBI's `oss`) is appended to the title as it is written.
    pub fn bioseq_title(&self) -> Vec<u8> {
        let mut title = Vec::new();
        for (line, defline) in self.deflines.iter().enumerate() {
            if !title.is_empty() {
                title.push(b' ');
            }
            let mut wrote_seqids = false;
            if line > 0 {
                for (at, id) in defline.ids.iter().enumerate() {
                    title.push(if at == 0 { b'>' } else { b'|' });
                    fasta(id, &mut title);
                    wrote_seqids = true;
                }
            }
            if let Some(text) = &defline.title {
                if wrote_seqids {
                    title.push(b' ');
                }
                title.extend_from_slice(text);
            }
        }
        title
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1844-1864
    // `x_GetFilteredHeader` drops definition lines only for a membership bit or a GI list,
    // which LOSAT's searches do not use.
    /// Every Seq-id of every definition line, in order (`CSeqDB::GetSeqIDs`).
    pub fn ids(&self) -> impl Iterator<Item = &ParsedSeqId> {
        self.deflines.iter().flat_map(|defline| defline.ids.iter())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqinfosrc_aux.cpp:51-69
    /// The subject Seq-id of the sequence's Seq-aligns (outfmt 8, 9, 11, 12 print it; NCBI's
    /// formatters find the subject by it): the first GI of any definition line, else the first
    /// id of the first definition line. `None` for a set without ids.
    pub fn alignment_id(&self) -> Option<&ParsedSeqId> {
        let best = find_best(self.ids(), best_rank_score)?;
        if best.gi().is_some() {
            Some(best)
        } else {
            self.ids().next()
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1328-1359
    // The data loader asks for the Bioseq of the Seq-align id (cached_sequence.cpp:58):
    // `target` is the `alignment_id` of the same set.
    /// The definition line of the subject Bioseq that NCBI's formatters read for `target` (the
    /// Seq-align id): the first that holds it (`s_SeqDB_SeqIdIn`, or a local id of its content
    /// label). Its ids and title are the Bioseq's (tabular `sseqid`, `sacc`, ...). `None` when
    /// no definition line holds it (NCBI throws).
    // NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2898-2916
    // `TruncateSpaces_Unsafe` (`eTrunc_Both`) is `ncbi_str::truncate_spaces`; `SetStrOrId` gives
    // an integer tag for a label that reads as one, which `compare` treats as the same tag.
    pub fn bioseq_defline(&self, target: &ParsedSeqId) -> Option<usize> {
        let mut label = Vec::new();
        content(target, true, &mut label);
        let acc = crate::ncbi_str::truncate_spaces(&label);
        let local = ParsedSeqId::Local(ObjectId::Str(acc.to_vec()));
        self.deflines
            .iter()
            .position(|defline| seqid_in(&defline.ids, target) || seqid_in(&defline.ids, &local))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1069-1092
    // The data loader asks for the Bioseq of the id it was asked for (cached_sequence.cpp:58), so
    // `target` is the Seq-align id by which the sequence is read; `x_GetTaxonomy` walks the lines in
    // this order and asks its volume's cache for each positive taxid (seqdbvol.cpp:1129-1153).
    /// The positive taxids of the definition lines in the order in which NCBI's volume looks them
    /// up where it loads this sequence's Bioseq for the id `target` (`CSeqDBVol::x_GetTaxonomy` of
    /// `x_GetTaxDefline`): the first line that holds `target` (`s_SeqDB_SeqIdIn`, without the
    /// local form that `bioseq_defline` also takes) first, then the others in their stored order;
    /// a taxid of several lines comes once for each (each is a lookup of the volume's cache).
    /// Without a target (`GetBioseq` of the oid: blastdbcmd without `-target_only`, the maskers, a
    /// sequence read by its ordinal id) the stored order.
    pub fn taxonomy_lookups(&self, target: Option<&ParsedSeqId>) -> Vec<i32> {
        self.taxonomy_order(target).collect()
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:1004-1011
    /// The taxid of the Bioseq's first source, whose name `CDeflineGenerator` makes an empty title
    /// of (`report::defline::BioseqFacts::taxname`): the first of `taxonomy_lookups(target)`, the
    /// one owner of that order for every reader (the search's headers, blastdbcmd, the maskers).
    /// With `-target_only` it is taken from every definition line the filter keeps, not from the
    /// target's line alone that the Bioseq's ids and title come from.
    pub fn first_source_taxid(&self, target: Option<&ParsedSeqId>) -> Option<i32> {
        self.taxonomy_order(target).next()
    }

    /// `taxonomy_lookups` as an iterator (the maskers' reader walks it without a list).
    pub(crate) fn taxonomy_order(&self, target: Option<&ParsedSeqId>) -> TaxonomyOrder<'_> {
        let front = target.and_then(|target| {
            self.deflines
                .iter()
                .position(|defline| seqid_in(&defline.ids, target))
        });
        TaxonomyOrder {
            deflines: &self.deflines,
            front,
            front_given: false,
            next: 0,
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4741-4785
    /// The definition line whose `WorstRank`-best id NCBI shows for the Seq-align id `aln_id`
    /// (`CAlignFormatUtil::GetDisplayIds`: outfmt 0 descriptions): the first that holds
    /// `aln_id` or a `general` id of the same database, else the last one. `None` for a set
    /// without definition lines.
    pub fn display_defline(&self, aln_id: &ParsedSeqId) -> Option<DisplayDefline> {
        let last = self.deflines.len().checked_sub(1)?;
        let found = self.deflines.iter().position(|defline| {
            defline.ids.iter().any(|id| {
                compare(id, aln_id) == Sic::Yes
                    || matches!((id, aln_id),
                        (ParsedSeqId::General(a), ParsedSeqId::General(b)) if a.db == b.db)
            })
        });
        let index = found.unwrap_or(last);
        let taxid = self.deflines[..=index]
            .iter()
            .rev()
            .find_map(|defline| defline.taxid)
            .unwrap_or(0);
        Some(DisplayDefline {
            index,
            found: found.is_some(),
            taxid,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:266-304
    // `PutTargetGiFirst` with -1 (the default `this_gi_first`) leaves the set as it is
    // (Blast_def_line_set.cpp:212-237).
    /// `CShowBlastDefline::GetBioseqHandleDeflineAndId` (outfmt 5 `Hit_id` and `Hit_def`):
    /// the id list (`seqid_list_string`) of the first definition line with a title (of the
    /// last one when none has a title), and that title followed by ` >[gi|N|]<fasta> <title>`
    /// for each later definition line with a title.
    pub fn hit_id_and_def(&self, show_gi: bool) -> (Vec<u8>, Vec<u8>) {
        let (mut id, mut title) = (Vec::new(), Vec::new());
        let mut is_first = true;
        for defline in &self.deflines {
            if is_first {
                id.clear();
                seqid_list_string(&defline.ids, show_gi, &mut id);
            }
            let Some(text) = &defline.title else {
                continue;
            };
            if !is_first {
                title.extend_from_slice(b" >");
                let gi = defline.gi();
                if show_gi && gi > 0 {
                    title.extend_from_slice(b"gi|");
                    number(&mut title, gi);
                    title.push(b'|');
                }
                if let Some(best) = defline.display_id() {
                    fasta(best, &mut title);
                }
                title.push(b' ');
            }
            title.extend_from_slice(text);
            is_first = false;
        }
        (id, title)
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:165-197
// NCBI dereferences the best id of an empty list (accession kinds); every definition line
// that makeblastdb writes has an id, and LOSAT writes `Unknown` for none. `eGi` prints a
// number: the tabular `sgi` and `sallgi` are `BlastDefline::gi` of the definition lines, as
// `seq_id_list` replaces no GI (`report::db_subject::ParsedSubject::tabular_gi`,
// `all_tabular_ids`).
/// The tabular id of `kind` of the id list `ids` (`s_GetSeqIdListString`): `Unknown` for an
/// empty result. The subject fields (`sseqid`, `sacc`, `saccver`) take the ids of the subject
/// Bioseq's definition line (`bioseq_defline`) after `seq_id_list`.
pub fn tabular_id(ids: &[ParsedSeqId], kind: TabularId, out: &mut Vec<u8>) {
    let start = out.len();
    match kind {
        TabularId::Full => seqid_list_string(ids, true, out),
        TabularId::Accession | TabularId::AccessionVersion => {
            if let Some(best) = find_best(ids, text_score) {
                content(best, kind == TabularId::AccessionVersion, out);
            }
        }
    }
    if out.len() == start {
        out.extend_from_slice(b"Unknown");
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:212-243
// The tabular formats call it for the subject Bioseq's ids (`SetSubjectId`, tabular.cpp:
// 522-532; its `s_ReplaceLocalId` keeps a local id's own text for a `-db` search, where
// `m_ParseSubjectDefline` is false) and for every definition line (`x_SetSubjectIds`).
/// `CShowBlastDefline::GetSeqIdList`: `ids` with each `general` id whose FASTA form holds
/// `gnl|BL_ORD_ID`, and each id whose FASTA form holds `lcl|Subject_`, replaced by `local str
/// "<title_word>"`; with an empty `title_word` the ids stay. `title_word` is the first word (the
/// bytes before the first space) of the subject Bioseq's generated title, `GenerateDefline` of a
/// `CDeflineGenerator` on the Bioseq handle, which the caller makes: the Bioseq's title is its
/// definition line's (seqdbvol.cpp:1445-1463; `unnamed protein product` for an empty protein title,
/// empty for an empty nucleotide one), and its Seq-ids can add a prefix (`TPA: ` for a `tpg`, `tpe`
/// or `tpd` id) or a generated title (a `patent` or `pdb` id with an empty title).
pub fn seq_id_list<'a>(ids: &'a [ParsedSeqId], title_word: &[u8]) -> Cow<'a, [ParsedSeqId]> {
    let replaced = |id: &ParsedSeqId| {
        let mut fasta_form = Vec::new();
        fasta(id, &mut fasta_form);
        (matches!(id, ParsedSeqId::General(_)) && contains(&fasta_form, b"gnl|BL_ORD_ID"))
            || contains(&fasta_form, b"lcl|Subject_")
    };
    if title_word.is_empty() || !ids.iter().any(replaced) {
        return Cow::Borrowed(ids);
    }
    Cow::Owned(
        ids.iter()
            .map(|id| {
                if replaced(id) {
                    ParsedSeqId::Local(ObjectId::Str(title_word.to_vec()))
                } else {
                    id.clone()
                }
            })
            .collect(),
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2417-2437
// `seqID` is the `WorstRank`-best id of the definition line, `gi` its first GI, `label`
// `CAlignFormatUtil::GetLabel` of `seqID` with `CSeq_id::eContent` (the enum is the `bool`
// parameter `with_version`: true) (showalign.cpp:2205-2209). `m_UseLongSeqIds` is the
// `[BLAST] LONG_SEQID` registry value, false by default.
/// The id that an outfmt 0 alignment block prints for `defline` (after `>` or ` `, before
/// ` <title>`; `CDisplaySeqalign::x_PrintDefLine`): `gi|N|` and the FASTA form with
/// `-show_gis` and a GI, else the bare id (a local id as its label).
pub fn alignment_line_id(defline: &BlastDefline, show_gis: bool, out: &mut Vec<u8>) {
    let Some(id) = defline.display_id() else {
        return;
    };
    let gi = defline.gi();
    if show_gis && gi > 0 && id.gi().is_none() {
        out.extend_from_slice(b"gi|");
        number(out, gi);
        out.push(b'|');
    }
    let mut fasta_form = Vec::new();
    fasta(id, &mut fasta_form);
    if !contains(&fasta_form, b"gnl|BL_ORD_ID") || contains(&fasta_form, b"lcl|Subject_") {
        if fasta_form.starts_with(b"lcl|") {
            format_label(id, true, out);
        } else if show_gis && gi > 0 {
            out.extend_from_slice(&fasta_form);
        } else {
            bare(id, out);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:886-913
// `sdl->id` and `sdl->gi` are the id and GI of `GetDisplayIds` (`display_defline`; the GI
// is 0 unless the definition line holds the Seq-align id).
/// The start of an outfmt 0 description line (`CShowBlastDefline`, text mode): `gi|N|` with
/// `-show_gis` and a GI, then the id (`GetBareId` and a space; the FASTA form with `-show_gis`
/// unless the id is the GI; a local id without `lcl|`). Returns whether `gi|N|` was written:
/// NCBI then puts two spaces before the title.
pub fn description_id(id: &ParsedSeqId, gi: i64, show_gis: bool, out: &mut Vec<u8>) -> bool {
    let wrote_gi = show_gis && gi > 0;
    if wrote_gi {
        out.extend_from_slice(b"gi|");
        number(out, gi);
        out.push(b'|');
    }
    let mut fasta_form = Vec::new();
    fasta(id, &mut fasta_form);
    if !(contains(&fasta_form, b"gnl|BL_ORD_ID") || contains(&fasta_form, b"lcl|Subject_")) {
        let mut id_str = Vec::new();
        if show_gis && id.gi().is_none() {
            id_str.extend_from_slice(&fasta_form);
        } else {
            bare(id, &mut id_str);
            id_str.push(b' ');
        }
        if id_str.starts_with(b"lcl|") {
            out.extend_from_slice(&fasta_form[4..]);
        } else {
            out.extend_from_slice(&id_str);
        }
    }
    wrote_gi
}

fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    haystack
        .windows(needle.len())
        .any(|window| window == needle)
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1031-1052
/// `s_SeqDB_SeqIdIn`: whether the first id of `ids` that is comparable with `target` matches.
fn seqid_in(ids: &[ParsedSeqId], target: &ParsedSeqId) -> bool {
    for id in ids {
        match compare(id, target) {
            Sic::Yes => return true,
            Sic::No => return false,
            Sic::Diff => {}
        }
    }
    false
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:392-482
/// `CSeq_id::Compare`.
pub fn compare(a: &ParsedSeqId, b: &ParsedSeqId) -> Sic {
    let yes_no = |matches: bool| if matches { Sic::Yes } else { Sic::No };
    if a.alternative() != b.alternative() {
        let comparable = |id: &ParsedSeqId| {
            matches!(
                id,
                ParsedSeqId::Genbank(_)
                    | ParsedSeqId::Embl(_)
                    | ParsedSeqId::Ddbj(_)
                    | ParsedSeqId::Tpg(_)
                    | ParsedSeqId::Tpe(_)
                    | ParsedSeqId::Tpd(_)
                    | ParsedSeqId::Gpipe(_)
                    | ParsedSeqId::NamedAnnotTrack(_)
            )
        };
        if !comparable(a) || !comparable(b) {
            return Sic::Diff;
        }
        return yes_no(textseq_match(a.textseq().unwrap(), b.textseq().unwrap()));
    }
    match (a, b) {
        (ParsedSeqId::Local(x), ParsedSeqId::Local(y)) => yes_no(object_id_match(x, y)),
        (ParsedSeqId::Gibbsq(x), ParsedSeqId::Gibbsq(y))
        | (ParsedSeqId::Gibbmt(x), ParsedSeqId::Gibbmt(y)) => yes_no(x == y),
        (ParsedSeqId::Giim(x), ParsedSeqId::Giim(y)) => yes_no(x.id == y.id),
        (ParsedSeqId::Patent(x), ParsedSeqId::Patent(y)) => yes_no(patent_match(x, y)),
        (ParsedSeqId::General(x), ParsedSeqId::General(y)) => {
            if dbtag_match(x, y) {
                Sic::Yes
            } else if !x.db.eq_ignore_ascii_case(&y.db) {
                Sic::Diff
            } else {
                Sic::No
            }
        }
        (ParsedSeqId::Gi(x), ParsedSeqId::Gi(y)) => yes_no(x == y),
        (ParsedSeqId::Pdb(x), ParsedSeqId::Pdb(y)) => yes_no(pdb_match(x, y)),
        _ => yes_no(textseq_match(a.textseq().unwrap(), b.textseq().unwrap())),
    }
}

// No NCBI counterpart: a necessary condition of `compare`'s `Yes`, for LOSAT's check that no two
// sequences of a database can be taken for one where a report joins following hits
// (`report::db_subject::join_keys`).
/// Adds to `keys` the keys of `id`: two ids that `compare` finds the same (`Sic::Yes`) share at
/// least one key, so ids that share none never do (ids that share one need not). A text Seq-id
/// has its accession and its name in upper case, whatever its alternative (a superset of the
/// alternatives `compare` compares across); a `pdb` id its molecule; the other alternatives their
/// alternative and their value, an object id as its number or its string in upper case
/// (`object_id_match`). The keys are 64-bit hashes of those values.
pub fn match_keys(id: &ParsedSeqId, keys: &mut Vec<u64>) {
    use std::hash::{Hash, Hasher};
    let mut key = |parts: &[&[u8]]| {
        let mut hasher = std::hash::DefaultHasher::new();
        parts.hash(&mut hasher);
        keys.push(hasher.finish());
    };
    let tag_key = |tag: &ObjectId| -> Vec<u8> {
        match object_id_number(tag) {
            Some(number) => format!("n{number}").into_bytes(),
            None => match tag {
                ObjectId::Str(text) => [b"s".as_slice(), &text.to_ascii_uppercase()].concat(),
                ObjectId::Id(_) => unreachable!("an integer tag is a number"),
            },
        }
    };
    match id {
        ParsedSeqId::Local(tag) => key(&[b"local", &tag_key(tag)]),
        ParsedSeqId::Gibbsq(number) => key(&[b"gibbsq", &number.to_le_bytes()]),
        ParsedSeqId::Gibbmt(number) => key(&[b"gibbmt", &number.to_le_bytes()]),
        ParsedSeqId::Giim(giim) => key(&[b"giim", &giim.id.to_le_bytes()]),
        ParsedSeqId::Patent(patent) => {
            let number = match &patent.cit.id {
                PatentNumber::Number(number) | PatentNumber::AppNumber(number) => number,
            };
            key(&[
                b"patent",
                &patent.seqid.to_le_bytes(),
                &patent.cit.country.to_ascii_uppercase(),
                &number.to_ascii_uppercase(),
            ])
        }
        ParsedSeqId::General(dbtag) => key(&[
            b"general",
            &dbtag.db.to_ascii_uppercase(),
            &tag_key(&dbtag.tag),
        ]),
        ParsedSeqId::Gi(gi) => key(&[b"gi", &gi.to_le_bytes()]),
        ParsedSeqId::Pdb(pdb) => key(&[b"pdb", &pdb.mol]),
        _ => {
            let text = id
                .textseq()
                .expect("the other alternatives are text Seq-ids");
            if let Some(accession) = &text.accession {
                key(&[b"accession", &accession.to_ascii_uppercase()]);
            }
            if let Some(name) = &text.name {
                key(&[b"name", &name.to_ascii_uppercase()]);
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objmgr/data_source.cpp:1502
// The handle of a text Seq-id is its accession and version in any alternative, in any case (probe
// of the DI-7 review: `NM_1.1` and `NM_1.2`, `NM_1` and `NM_1.1`, `pdb|1ABC|A` and `|B` are two
// sequences to NCBI's scope; `gb|X.1|` and `emb|X.1|`, `abc` and `ABC` are one).
/// Adds to `keys` the keys of `id` as NCBI's scope finds a sequence by it: two ids that the scope
/// finds as one share a key. Other than `match_keys`, a version and a PDB chain tell ids apart.
pub fn scope_keys(id: &ParsedSeqId, keys: &mut Vec<u64>) {
    use std::hash::{Hash, Hasher};
    let mut key = |parts: &[&[u8]]| {
        let mut hasher = std::hash::DefaultHasher::new();
        parts.hash(&mut hasher);
        keys.push(hasher.finish());
    };
    if let ParsedSeqId::Pdb(pdb) = id {
        key(&[
            b"pdb",
            &pdb.mol.to_ascii_uppercase(),
            &pdb.chain.unwrap_or(32).to_le_bytes(),
            pdb.chain_id.as_deref().unwrap_or_default(),
        ]);
    } else if let Some(text) = id.textseq() {
        let version = text
            .version
            .map_or(Vec::new(), |v| v.to_le_bytes().to_vec());
        match (&text.accession, &text.name) {
            (Some(accession), _) => {
                key(&[b"accession", &accession.to_ascii_uppercase(), &version]);
            }
            (None, Some(name)) => key(&[b"name", &name.to_ascii_uppercase()]),
            (None, None) => {}
        }
    } else {
        match_keys(id, keys);
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1064-1095
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1713-1719
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:1141-1148 (`FindInfo`; the
// packed key of an accession, `ParseAcc` at 825-880, holds its prefix without case, its digits
// and its version, not the alternative)
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.hpp:767-768 (`CSeq_id_Local_Tree`),
// 920-925 (`CSeq_id_General_Tree`) and 999-1008 (`CSeq_id_Patent_Tree`)
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_id_tree.cpp:3467-3481 (the PDB tree's key;
// `CSeq_id_PDB_Info::Normalize` upper-cases the molecule after its first letter and sets
// `chain_id` from `chain`, then the tree compares the normalized ids, `Equals`)
// No NCBI counterpart: LOSAT's sets of `CSeq_id_Handle` (psiblast's inclusion sets and the rows
// of its PSSM, `algorithm::psiblast::rounds::SetKey`) take two ids for one where they are equal
// (`ParsedSeqId`'s `Eq`, which the handle of each id follows, but for a `giim` id, which those
// sets refuse) and refuse two other ids that share this key, where the handles can be one or two.
/// Adds to `out` the key of `id` that every id which NCBI's `CSeq_id_Handle` takes for the same
/// sequence shares: a text Seq-id's accession and name in upper case, its version and release,
/// with GenBank, EMBL and DDBJ as one alternative (their tree; an accession that packs is one
/// handle in the three); a local id's tag and a `general` id's database and tag as `match_keys`
/// takes them (a number, or a string in upper case); a patent's country, number (or application
/// number, as one) in upper case and the sequence's number; a PDB id's molecule in upper case and
/// its chain as written (`chain_id`, else `chain`, read as its default); the others' values. Ids
/// that NCBI keeps apart can share it (`gb|ABC1.1|` and `emb|abc1.1|`).
pub fn handle_superset_key(id: &ParsedSeqId, out: &mut Vec<u8>) {
    let mut part = |bytes: &[u8]| {
        out.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
        out.extend_from_slice(bytes);
    };
    let tag_key = |tag: &ObjectId| -> Vec<u8> {
        match object_id_number(tag) {
            Some(number) => format!("n{number}").into_bytes(),
            None => match tag {
                ObjectId::Str(text) => [b"s".as_slice(), &text.to_ascii_uppercase()].concat(),
                ObjectId::Id(_) => unreachable!("an integer tag is a number"),
            },
        }
    };
    let option = |value: Option<&[u8]>| -> Vec<u8> {
        value.map_or(b"-".to_vec(), |value| [b"+".as_slice(), value].concat())
    };
    match id {
        ParsedSeqId::Local(tag) => {
            part(b"local");
            part(&tag_key(tag));
        }
        ParsedSeqId::Gibbsq(number) => {
            part(b"gibbsq");
            part(&number.to_le_bytes());
        }
        ParsedSeqId::Gibbmt(number) => {
            part(b"gibbmt");
            part(&number.to_le_bytes());
        }
        ParsedSeqId::Giim(giim) => {
            part(b"giim");
            part(&giim.id.to_le_bytes());
        }
        ParsedSeqId::Patent(patent) => {
            let number = match &patent.cit.id {
                PatentNumber::Number(number) | PatentNumber::AppNumber(number) => number,
            };
            part(b"patent");
            part(&patent.cit.country.to_ascii_uppercase());
            part(&number.to_ascii_uppercase());
            part(&patent.seqid.to_le_bytes());
        }
        ParsedSeqId::General(dbtag) => {
            part(b"general");
            part(&dbtag.db.to_ascii_uppercase());
            part(&tag_key(&dbtag.tag));
        }
        ParsedSeqId::Gi(gi) => {
            part(b"gi");
            part(&gi.to_le_bytes());
        }
        ParsedSeqId::Pdb(pdb) => {
            let chain = match &pdb.chain_id {
                Some(chain_id) => chain_id.clone(),
                None => vec![pdb.chain_or_default() as u8],
            };
            part(b"pdb");
            part(&pdb.mol.to_ascii_uppercase());
            part(&chain);
        }
        _ => {
            let text = id
                .textseq()
                .expect("the other alternatives are text Seq-ids");
            let alternative = match id {
                ParsedSeqId::Genbank(_) | ParsedSeqId::Embl(_) | ParsedSeqId::Ddbj(_) => "gb",
                _ => id.alternative_name(),
            };
            part(alternative.as_bytes());
            part(&option(
                text.accession
                    .as_ref()
                    .map(|accession| accession.to_ascii_uppercase())
                    .as_deref(),
            ));
            part(&option(
                text.name
                    .as_ref()
                    .map(|name| name.to_ascii_uppercase())
                    .as_deref(),
            ));
            part(&option(
                text.version
                    .map(|version| version.to_le_bytes())
                    .as_ref()
                    .map(|bytes| bytes.as_slice()),
            ));
            part(&option(text.release.as_deref()));
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Patent_seq_id.cpp:53-56
// NCBI reference (598d8ae6): c++/src/objects/biblio/Id_pat.cpp:55-79
// `doc-type` takes no part.
fn patent_match(a: &PatentSeqId, b: &PatentSeqId) -> bool {
    a.seqid == b.seqid
        && a.cit.country.eq_ignore_ascii_case(&b.cit.country)
        && match (&a.cit.id, &b.cit.id) {
            (PatentNumber::Number(x), PatentNumber::Number(y))
            | (PatentNumber::AppNumber(x), PatentNumber::AppNumber(y)) => x.eq_ignore_ascii_case(y),
            _ => false,
        }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Textseq_id.cpp:143-174
fn textseq_match(a: &TextseqId, b: &TextseqId) -> bool {
    let versions = || match (a.version, b.version) {
        (Some(x), Some(y)) => x == y,
        _ => true,
    };
    if let (Some(x), Some(y)) = (&a.accession, &b.accession) {
        return x.eq_ignore_ascii_case(y) && versions();
    }
    if let (Some(x), Some(y)) = (&a.name, &b.name) {
        return x.eq_ignore_ascii_case(y) && versions();
    }
    false
}

// NCBI reference (598d8ae6): c++/src/objects/general/Dbtag.cpp:158-163
fn dbtag_match(a: &Dbtag, b: &Dbtag) -> bool {
    a.db.eq_ignore_ascii_case(&b.db) && object_id_match(&a.tag, &b.tag)
}

// NCBI reference (598d8ae6): c++/src/objects/general/Object_id.cpp:93-161
fn object_id_match(a: &ObjectId, b: &ObjectId) -> bool {
    match (object_id_number(a), object_id_number(b)) {
        (Some(x), Some(y)) => x == y,
        (None, None) => match (a, b) {
            (ObjectId::Str(x), ObjectId::Str(y)) => x.eq_ignore_ascii_case(y),
            _ => false,
        },
        _ => false,
    }
}

/// The value of an Object-id that `GetIdType` reads as an integer: an `id`, or a `str` that
/// `NStr::StringToInt8` with `fConvErr_NoThrow` and `fConvErr_NoErrMessage` reads
/// (`ncbi_str::string_to_int8` with no flags: no spaces, and the string ends at a 0 byte;
/// `NoThrow` gives 0 with `errno` set for an error, which is `e_Str`), without a `+` or leading
/// zeros; a 0 only as the whole string.
fn object_id_number(tag: &ObjectId) -> Option<i64> {
    match tag {
        ObjectId::Id(id) => Some(i64::from(*id)),
        ObjectId::Str(text) => {
            // `StringToInt8` with no flags fails at a first byte other than a digit or a sign
            // (no spaces): such a text is `e_Str` without the conversion and its exception.
            if !text
                .first()
                .is_some_and(|&c| c.is_ascii_digit() || c == b'+' || c == b'-')
            {
                return None;
            }
            match crate::ncbi_str::string_to_int8(text, 0).ok()? {
                0 => (text.len() == 1).then_some(0),
                value if value > 0 => (!matches!(text[0], b'0' | b'+')).then_some(value),
                value => (text[1] != b'0').then_some(value),
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:486-522
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Textseq_id.cpp:178-224
// NCBI reference (598d8ae6): c++/src/objects/general/Object_id.cpp:145-161
// NCBI reference (598d8ae6): c++/src/objects/general/Dbtag.cpp:176-183
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Patent_seq_id.cpp:60-67
// NCBI reference (598d8ae6): c++/src/objects/biblio/Id_pat.cpp:55-79
// NCBI reference (598d8ae6): c++/src/objects/seqloc/PDB_seq_id.cpp:70-91
/// A key of `id` that two Seq-ids share exactly when `CSeq_id::CompareOrdered` finds them equal,
/// as a map keyed by Seq-ids does (`SSeqIdKey`): the alternative, then its members as the
/// comparison reads them. Letters compared with `PNocase` go in lower case. A Textseq-id is its
/// accession and version where it has both (its name is not read then), else its accession,
/// version and name as far as each is set; an Object-id is its integer where `GetIdType` reads
/// one (`object_id_number`), else its string; a `general` id its database and tag; a `patent` id
/// its sequence number, country and number (`CId_pat::Match`; the document type is not read); a
/// `pdb` id its molecule and its chain id, or the chain's letter (case kept); a `giim` id its id.
pub(crate) fn ordered_key(id: &ParsedSeqId) -> Vec<u8> {
    fn text(key: &mut Vec<u8>, bytes: &[u8], nocase: bool) {
        key.extend_from_slice(&(bytes.len() as u64).to_le_bytes());
        key.extend(
            bytes
                .iter()
                .map(|&b| if nocase { b.to_ascii_lowercase() } else { b }),
        );
    }
    fn set_text(key: &mut Vec<u8>, bytes: Option<&Vec<u8>>) {
        match bytes {
            Some(bytes) => {
                key.push(1);
                text(key, bytes, true);
            }
            None => key.push(0),
        }
    }
    fn tag(key: &mut Vec<u8>, tag: &ObjectId) {
        match (object_id_number(tag), tag) {
            (Some(value), _) => {
                key.push(0);
                key.extend_from_slice(&value.to_le_bytes());
            }
            (None, ObjectId::Str(string)) => {
                key.push(1);
                text(key, string, true);
            }
            // `object_id_number` reads every `id` as its integer.
            (None, ObjectId::Id(value)) => {
                key.push(0);
                key.extend_from_slice(&i64::from(*value).to_le_bytes());
            }
        }
    }
    let mut key = vec![id.alternative()];
    match id {
        ParsedSeqId::Local(local) => tag(&mut key, local),
        ParsedSeqId::Gibbsq(value) | ParsedSeqId::Gibbmt(value) => {
            key.extend_from_slice(&value.to_le_bytes())
        }
        ParsedSeqId::Giim(giim) => key.extend_from_slice(&giim.id.to_le_bytes()),
        ParsedSeqId::Gi(gi) => key.extend_from_slice(&gi.to_le_bytes()),
        ParsedSeqId::General(dbtag) => {
            text(&mut key, &dbtag.db, true);
            tag(&mut key, &dbtag.tag);
        }
        ParsedSeqId::Patent(patent) => {
            key.extend_from_slice(&patent.seqid.to_le_bytes());
            text(&mut key, &patent.cit.country, true);
            match &patent.cit.id {
                PatentNumber::Number(number) => {
                    key.push(0);
                    text(&mut key, number, true);
                }
                PatentNumber::AppNumber(number) => {
                    key.push(1);
                    text(&mut key, number, true);
                }
            }
        }
        ParsedSeqId::Pdb(pdb) => {
            text(&mut key, &pdb.mol, true);
            // `s_GetEffectiveChain_id`: the chain id, else the chain as one letter, else nothing.
            let chain = pdb.chain.map(|chain| [chain as u8]);
            let effective: &[u8] = match (&pdb.chain_id, &chain) {
                (Some(chain_id), _) => chain_id,
                (None, Some(chain)) => chain,
                (None, None) => &[],
            };
            text(&mut key, effective, false);
        }
        _ => {
            if let Some(textseq) = id.textseq() {
                set_text(&mut key, textseq.accession.as_ref());
                match textseq.version {
                    Some(version) => {
                        key.push(1);
                        key.extend_from_slice(&version.to_le_bytes());
                    }
                    None => key.push(0),
                }
                if textseq.accession.is_none() || textseq.version.is_none() {
                    set_text(&mut key, textseq.name.as_ref());
                }
            }
        }
    }
    key
}

// NCBI reference (598d8ae6): c++/src/objects/seqloc/PDB_seq_id.cpp:53-67
fn pdb_match(a: &PdbSeqId, b: &PdbSeqId) -> bool {
    if let (Some(x), Some(y)) = (a.chain, b.chain) {
        if x != y {
            return false;
        }
    }
    if let (Some(x), Some(y)) = (&a.chain_id, &b.chain_id) {
        if x != y {
            return false;
        }
    }
    a.mol == b.mol
}

/// The positive taxids of a sequence's definition lines in the order of NCBI's lookups
/// (`DeflineSet::taxonomy_order`): the line `front` first, then the others in their stored order.
/// A plain loop over the lines (no iterator adapters: the readers ask its first item of every
/// sequence they read).
#[derive(Clone)]
pub(crate) struct TaxonomyOrder<'a> {
    deflines: &'a [BlastDefline],
    front: Option<usize>,
    front_given: bool,
    next: usize,
}

impl Iterator for TaxonomyOrder<'_> {
    type Item = i32;

    fn next(&mut self) -> Option<i32> {
        if !self.front_given {
            self.front_given = true;
            if let Some(taxid) = self
                .front
                .and_then(|index| self.deflines[index].taxid)
                .filter(|&taxid| taxid > 0)
            {
                return Some(taxid);
            }
        }
        while let Some(defline) = self.deflines.get(self.next) {
            let index = self.next;
            self.next += 1;
            if Some(index) == self.front {
                continue;
            }
            if let Some(taxid) = defline.taxid.filter(|&taxid| taxid > 0) {
                return Some(taxid);
            }
        }
        None
    }
}

#[cfg(test)]
#[path = "defline_set_tests.rs"]
mod tests;
