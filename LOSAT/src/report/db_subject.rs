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
//! `algorithm::common::database::read_title`); `SeqDb::read_header` tells the two apart.

use crate::seqid::defline_set::{
    alignment_line_id, compare, description_id, seq_id_list, tabular_id, Sic, TabularId,
};
use crate::seqid::model::{ObjectId, PatentNumber};
use crate::seqid::rank::{find_best, text_score};
use crate::seqid::{label, BlastDefline, DeflineSet, ParsedSeqId};

use super::defline::{
    bioseq_title_from_source, generate_bioseq_defline, BioseqIds, Defline, UnknownEncoding,
};

/// A reported database sequence with parsed Seq-ids: its definition lines, its Seq-align id
/// and the definition line of its Bioseq.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ParsedSubject {
    set: DeflineSet,
    alignment_id: ParsedSeqId,
    /// `DeflineSet::bioseq_defline` of `alignment_id`.
    target: usize,
    /// Whether NCBI's formatters cannot find the sequence by its Seq-align id (`fetched`).
    unfetchable: bool,
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
}

impl ParsedSubject {
    // NCBI reference (598d8ae6): c++/src/objtools/data_loaders/blastdb/bdbloader.cpp:406-415
    // ```c++
    //     int oid = -1;
    //     ...
    //     if (m_BlastDb->SeqidToOid(*seq_id, oid)) {
    // ```
    // NCBI's formatters find the subject's sequence again from its Seq-align id. For a gibbsq
    // id, a `general` id of the database `PIG`, a named-annotation-track id and a `general`
    // id `BL_ORD_ID` that makeblastdb stored as written, the lookup finds no sequence or
    // another one: NCBI 2.17.0 stops with a CCoreException in outfmt 6 (exit status 255) and
    // writes the title `Unknown`, or another sequence's, in outfmt 0 and 5 (probe of line
    // DBID, DECISIONS 8 (3); `LOSAT/tests/seqid/fasta/variants.fa`). In a version 5 database
    // (LMDB) a local integer id (`lcl|123`) is not found either (the same failures; probe of
    // ID2a, `blastdb-ids-20261021/R/probe`), where version 4 (ISAM) finds it. LOSAT refuses
    // these subjects instead of repeating a defect of NCBI's.
    /// The subject of the definition lines `set` of a database of format version 5
    /// (`version5`) or 4. `Err` names what LOSAT does not write as NCBI: a stored title or
    /// Seq-id with a byte outside ' ' ..= '~' (NCBI writes an error message each time it reads
    /// the definition line) and a set without a Seq-id.
    pub fn new(set: DeflineSet, version5: bool) -> Result<Self, String> {
        if set
            .deflines
            .iter()
            .any(|defline| defline.replaced_title_bytes > 0 || defline.replaced_id_bytes > 0)
        {
            return Err("has a stored title or Seq-id with a byte outside ' ' to '~' (NCBI BLAST+ shows it as '#' and writes an error message each time it reads the definition line)".into());
        }
        let Some(alignment_id) = set.alignment_id().cloned() else {
            return Err("has no Seq-id".into());
        };
        let unfetchable = unfetchable(&alignment_id, version5);
        let Some(target) = set.bioseq_defline(&alignment_id) else {
            return Err("has no definition line that holds its Seq-align id".into());
        };
        Ok(Self {
            set,
            alignment_id,
            target,
            unfetchable,
        })
    }

    /// `Err` for a subject whose sequence NCBI's formatters cannot find by its Seq-align id
    /// (see `new`): the formats that read the subject's Bioseq (outfmt 0, 5, 6, 7 and 13 to 16)
    /// refuse it; the Seq-annots and the archive (8, 9, 11, 12) write the Seq-align id alone.
    pub fn fetched(&self) -> Result<(), String> {
        if !self.unfetchable {
            return Ok(());
        }
        Err(format!(
            "has the Seq-align id {}, by which NCBI BLAST+'s formatters do not find the sequence (NCBI BLAST+ 2.17.0 writes Unknown or another sequence's names, or stops)",
            String::from_utf8_lossy(&self.alignment_fasta())
        ))
    }

    /// The Seq-align id (outfmt 8, 9, 11, 12; `DeflineSet::alignment_id`).
    pub fn alignment_id(&self) -> &ParsedSeqId {
        &self.alignment_id
    }

    /// The definition lines.
    pub fn deflines(&self) -> &[BlastDefline] {
        &self.set.deflines
    }

    /// The definition line of the subject Bioseq.
    fn target(&self) -> &BlastDefline {
        &self.set.deflines[self.target]
    }

    /// What `GenerateDefline` reads of the subject Bioseq's Seq-ids.
    fn bioseq_ids(&self) -> BioseqIds {
        bioseq_ids(&self.target().ids)
    }

    /// `GenerateDefline` of the subject Bioseq: for the alignment heading
    /// (`leave_prefix_suffix` false) or the description table (true).
    pub fn bioseq_title(
        &self,
        protein: bool,
        leave_prefix_suffix: bool,
    ) -> Result<Defline, UnknownEncoding> {
        let title = self.target().title.as_deref().unwrap_or_default();
        generate_bioseq_defline(title, protein, leave_prefix_suffix, &self.bioseq_ids())
    }

    /// Whether NCBI makes the subject Bioseq's title from its source (`bioseq_title_from_source`)
    /// and that source has a taxid: the title then depends on the taxonomy database.
    pub fn title_from_taxonomy(&self) -> bool {
        let title = self.target().title.as_deref().unwrap_or_default();
        bioseq_title_from_source(title, &self.bioseq_ids())
            && self
                .set
                .deflines
                .iter()
                .any(|defline| defline.taxid.is_some_and(|taxid| taxid != 0))
    }

    /// The first word of the subject Bioseq's generated title, which `seq_id_list` puts in the
    /// place of a `gnl|BL_ORD_ID` or `lcl|Subject_` id of the Bioseq; empty when no such id is
    /// there (the title is then not made). `Err` where NCBI's `x_CleanAndCompress` reads past the
    /// end of the title (NCBI crashes).
    fn title_word(&self, protein: bool) -> Result<Vec<u8>, &'static str> {
        let ids = &self.target().ids;
        if matches!(seq_id_list(ids, b"x"), std::borrow::Cow::Borrowed(_)) {
            return Ok(Vec::new());
        }
        let defline = self
            .bioseq_title(protein, false)
            .map_err(|UnknownEncoding| "a title in no encoding that NCBI BLAST+ recognises")?;
        if defline.reads_past_end {
            return Err("a title of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes)");
        }
        let end = defline
            .text
            .iter()
            .position(|&byte| byte == b' ')
            .unwrap_or(defline.text.len());
        Ok(defline.text[..end].to_vec())
    }

    /// The tabular `sseqid`, `sacc` and `saccver` (`tabular_id` of the subject Bioseq's ids after
    /// `seq_id_list`).
    pub fn tabular_ids(&self, protein: bool) -> Result<[Vec<u8>; 3], &'static str> {
        let word = self.title_word(protein)?;
        let ids = seq_id_list(&self.target().ids, &word);
        let field = |kind| {
            let mut out = Vec::new();
            tabular_id(&ids, kind, &mut out);
            out
        };
        Ok([
            field(TabularId::Full),
            field(TabularId::Accession),
            field(TabularId::AccessionVersion),
        ])
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

    /// The start of the outfmt 0 description line (`description_id` of `GetDisplayIds`' id,
    /// without `-show_gis`, which LOSAT refuses) and its title: `GenerateDefline` of the subject
    /// Bioseq with `fLeavePrefixSuffix` (`Unknown` where NCBI's `HtmlDecode` throws).
    pub fn description(&self, protein: bool) -> (Vec<u8>, Defline) {
        let mut start = Vec::new();
        if let Some(shown) = self.set.display_defline(&self.alignment_id) {
            let defline = &self.set.deflines[shown.index];
            if let Some(id) = defline.display_id() {
                let gi = if shown.found { defline.gi() } else { 0 };
                description_id(id, gi, false, &mut start);
            }
        }
        let title = self.bioseq_title(protein, true).unwrap_or(Defline {
            text: b"Unknown".to_vec(),
            reads_past_end: false,
        });
        (start, title)
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
    /// The definition-line lines of the outfmt 0 alignment of the subject (each ending with a
    /// newline; `Length=` follows): every definition line, the first after `>` and the others
    /// after a space, with its id (`alignment_line_id`) and its title wrapped at 60 bytes.
    pub fn heading(&self, protein: bool) -> Vec<u8> {
        let mut out = Vec::new();
        let mut generated = None;
        for (i, defline) in self.set.deflines.iter().enumerate() {
            let start = out.len();
            out.push(if i == 0 { b'>' } else { b' ' });
            alignment_line_id(defline, false, &mut out);
            if i > 0 || out.len() - start > 1 {
                out.push(b' ');
            }
            let title = match defline.title.as_deref() {
                Some(title) if !title.is_empty() => title,
                _ => generated.get_or_insert_with(|| {
                    self.bioseq_title(protein, false)
                        .map(|defline| defline.text)
                        .unwrap_or_default()
                }),
            };
            wrap_title(title, &mut out);
            out.push(b'\n');
        }
        out
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
                }
            })
            .collect()
    }

    /// The FASTA form of the Seq-align id (the `HitDescr` id that NCBI writes when it cannot
    /// make the descriptions).
    pub fn alignment_fasta(&self) -> Vec<u8> {
        let mut out = Vec::new();
        label::fasta(&self.alignment_id, &mut out);
        out
    }
}

/// Whether NCBI's formatters cannot find a sequence by the Seq-align id `id` in a database of
/// version 5 (`version5`) or 4 (`ParsedSubject::new`).
fn unfetchable(id: &ParsedSeqId, version5: bool) -> bool {
    match id {
        ParsedSeqId::Gibbsq(_) | ParsedSeqId::NamedAnnotTrack(_) => true,
        ParsedSeqId::Local(ObjectId::Id(_)) => version5,
        ParsedSeqId::General(dbtag) => {
            dbtag.db.eq_ignore_ascii_case(b"PIG") || dbtag.db == b"BL_ORD_ID"
        }
        _ => false,
    }
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
/// Whether `CSeq_id::IdentifyAccession` gives `accession` a RefSeq chromosome, mRNA or ncRNA
/// type (`BioseqIds::refseq_source_title`).
fn refseq_source_accession(accession: &[u8]) -> bool {
    if accession.len() < 3 || accession[2] != b'_' {
        return false;
    }
    let digits = &accession[3..];
    if !digits.iter().all(u8::is_ascii_digit) {
        return false;
    }
    let prefix = [
        accession[0].to_ascii_uppercase(),
        accession[1].to_ascii_uppercase(),
    ];
    match &prefix {
        b"AC" | b"NC" | b"NR" => digits.len() == 6,
        b"NM" | b"XM" => digits.len() == 6 || digits.len() == 9,
        _ => false,
    }
}

/// `BioseqIds` of the Seq-ids `ids` of a Bioseq.
fn bioseq_ids(ids: &[ParsedSeqId]) -> BioseqIds {
    let mut flags = BioseqIds::default();
    for id in ids {
        match id {
            ParsedSeqId::Tpg(text) | ParsedSeqId::Tpe(text) | ParsedSeqId::Tpd(text) => {
                flags.third_party = true;
                flags.refseq_source_title |= text
                    .accession
                    .as_deref()
                    .is_some_and(refseq_source_accession);
            }
            ParsedSeqId::Other(text)
            | ParsedSeqId::Genbank(text)
            | ParsedSeqId::Embl(text)
            | ParsedSeqId::Ddbj(text) => {
                flags.refseq_source_title |= text
                    .accession
                    .as_deref()
                    .is_some_and(refseq_source_accession);
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seqdb::{SeqDb, SeqType};
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
                return ParsedSubject::new(set, false).unwrap();
            }
        }
        panic!("no sequence {first}");
    }

    fn text(bytes: &[u8]) -> &str {
        std::str::from_utf8(bytes).unwrap()
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
            text(&subject("tpg|DAA10002.1|").heading(true)),
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
            text(&subject("ref|NC_100013.1|").heading(true)),
            ">NC_100013.1 \n"
        );
        assert_eq!(
            description("ref|WP_100014.1|"),
            "WP_100014.1 unnamed protein product"
        );
        assert_eq!(description("lcl|Subject_27"), "subject seven title");
        assert_eq!(
            text(&subject("lcl|Subject_30").heading(true)),
            ">Subject_30 unnamed protein product\n"
        );
        assert_eq!(description("pdb|5ABC|D"), "5ABC_D Chain D");
        assert_eq!(
            text(&subject("pdb|5ABC|D").heading(true)),
            ">5ABC_D Chain D\n XP_100015.1 second\n"
        );
        assert_eq!(
            text(&subject("lcl|Subject_36").heading(true)),
            ">Subject_36 unnamed protein product\n DAA10016.1 tpa second\n"
        );
        // A title from the source: the taxonomy database would name its organism.
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
                .map(|id| String::from_utf8(id).unwrap())
        };
        assert_eq!(ids("lcl|Subject_36"), ["unnamed", "unnamed", "unnamed"]);
        assert_eq!(ids("lcl|Subject_27"), ["subject", "subject", "subject"]);
        assert_eq!(
            ids("lcl|Subject_28"),
            ["tpg|DAA10008.1|", "DAA10008", "DAA10008.1"]
        );
        assert_eq!(ids("lcl|Subject_29"), ["pdb|3ABC|C", "3ABC_C", "3ABC_C"]);
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

    #[test]
    fn seq_align_ids_that_ncbi_does_not_find() {
        let local_id = one_defline(vec![ParsedSeqId::Local(ObjectId::Id(123))]);
        assert!(ParsedSubject::new(local_id.clone(), false)
            .unwrap()
            .fetched()
            .is_ok());
        let error = ParsedSubject::new(local_id, true)
            .unwrap()
            .fetched()
            .unwrap_err();
        assert!(error.contains("Seq-align id lcl|123,"), "{error}");
        for id in [
            ParsedSeqId::Gibbsq(123),
            ParsedSeqId::General(Dbtag {
                db: b"PIG".to_vec(),
                tag: ObjectId::Id(123),
            }),
            ParsedSeqId::General(Dbtag {
                db: b"BL_ORD_ID".to_vec(),
                tag: ObjectId::Id(5),
            }),
        ] {
            let subject = ParsedSubject::new(one_defline(vec![id]), false).unwrap();
            assert!(subject.fetched().is_err());
        }
        // A GI in the set makes the GI the Seq-align id.
        let with_gi = one_defline(vec![ParsedSeqId::Gibbsq(123), ParsedSeqId::Gi(7)]);
        assert!(ParsedSubject::new(with_gi, true).unwrap().fetched().is_ok());
        let replaced = DeflineSet {
            deflines: vec![BlastDefline {
                ids: vec![ParsedSeqId::Gi(7)],
                replaced_title_bytes: 1,
                ..BlastDefline::default()
            }],
        };
        assert!(ParsedSubject::new(replaced, false).is_err());
        assert!(ParsedSubject::new(DeflineSet::default(), false).is_err());
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
            ParsedSubject::new(one_defline(vec![id]), false).unwrap()
        };
        assert!(joined_in_reports(&pdb(None), &pdb(Some("AA"))));
        assert!(!joined_in_reports(&pdb(Some("A")), &pdb(Some("B"))));
    }

    #[test]
    fn refseq_accessions_whose_empty_titles_come_from_the_source() {
        for (accession, from_source) in [
            ("NC_000001", true),
            ("AC_000001", true),
            ("NM_000001", true),
            ("NM_000000001", true),
            ("XM_000000001", true),
            ("NR_000001", true),
            ("NR_000000001", false),
            ("NC_0000001", false),
            ("NP_000001", false),
            ("WP_000000001", false),
            ("XR_000001", false),
        ] {
            assert_eq!(
                refseq_source_accession(accession.as_bytes()),
                from_source,
                "{accession}"
            );
        }
    }
}
