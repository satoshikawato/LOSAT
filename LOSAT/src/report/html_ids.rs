//! The ids of `-html` reports: the link of a sequence's id (NCBI's `CAlignFormatUtil::GetIDUrlGen`),
//! the label of its score link and the anchor of its alignment heading. Reached only with `-html`:
//! the plain report never classifies an id here.
// The description rows (D5) and headings (D6, D7) call these next (T4a2 piece A, in progress).
#![cfg_attr(not(test), allow(dead_code))]

use crate::seqid::label;
use crate::seqid::model::ParsedSeqId;
use crate::seqid::rank::{find_best, text_score};

/// Where a link is written: a description row (`x_FillDeflineAndId`, `isAlignLink` false, NCBI's
/// empty RID) or an alignment heading (`x_InitSeqUrl`, `isAlignLink` true, RID `0`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LinkPlace {
    Description,
    Heading,
}

/// What the link of one id list depends on besides its ids (NCBI's `SSeqURLInfo` as the
/// applications fill it: no templates, no CSS information, no new window).
#[derive(Clone, Copy, Debug)]
pub struct SeqUrlInfo<'a> {
    /// `accession`: `GetLabel(eContent)` of the description's display id, or
    /// `CAlignFormatUtil::GetLabel(id, true)` of the heading's (`description_accession`,
    /// `heading_label`).
    pub accession: &'a [u8],
    /// `blast_rank`: the subject's rank in the query's hit list, from 1.
    pub rank: usize,
    /// `isDbNa`: the subjects are nucleotide sequences.
    pub nucleotide: bool,
    pub place: LinkPlace,
    /// `taxid` (the `<@taxid@>` of a local id's template; 0 for none).
    pub taxid: i32,
}

/// The formatter configuration that the links read: `[LOCAL_ID] TOOL_URL` of the file that
/// NCBI's `InitConfig` opens, the `-db` text as typed (empty for `-subject`) and the protocol of
/// `MapProtocol` (`[BLASTFMTUTIL] PROTOCOL` of that file, else `https:`). The registry
/// (`report/formatter_registry.rs`) gives it; `IdUrlConfig::DEFAULT` is a run without that file.
#[derive(Clone, Copy, Debug)]
pub struct IdUrlConfig<'a> {
    pub local_id_url: Option<&'a [u8]>,
    pub db_name: &'a [u8],
    pub protocol: &'a [u8],
}

impl IdUrlConfig<'static> {
    /// No configuration file: no local-id link, protocol `https:`.
    pub const DEFAULT: Self = IdUrlConfig {
        local_id_url: None,
        db_name: b"",
        protocol: b"https:",
    };
}

// NCBI reference (2.17.0): c++/src/objtools/align_format/align_format_util.cpp:181,187
// ```c++
//     { "ENTREZ", "<a title=\"Show report for <@acc@>\" <@cssInf@>href=\"https://www.ncbi.nlm.nih.gov/<@db@>/<@acc@>?report=genbank&log$=<@log@>&blast_rank=<@blast_rank@>&RID=<@rid@>\" <@target@>>" },
//     { "TRACE",  "<a title=\"Show report for <@val@>\" <@cssInf@>href=\"https://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=retrieve&dopt=fasta&val=<@val@>&RID=<@rid@>\">" },
// ```
const ENTREZ: &[u8] = b"<a title=\"Show report for <@acc@>\" <@cssInf@>href=\"https://www.ncbi.nlm.nih.gov/<@db@>/<@acc@>?report=genbank&log$=<@log@>&blast_rank=<@blast_rank@>&RID=<@rid@>\" <@target@>>";
const TRACE: &[u8] = b"<a title=\"Show report for <@val@>\" <@cssInf@>href=\"https://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=retrieve&dopt=fasta&val=<@val@>&RID=<@rid@>\">";

// NCBI reference (2.17.0): c++/src/objtools/align_format/align_format_util.cpp:3593-3599
// ```c++
// string CAlignFormatUtil::MapTemplate(string inpString,string tmplParamName,string templParamVal)
// {
//     string outString;
//     string tmplParam = "<@" + tmplParamName + "@>";
//     NStr::Replace(inpString,tmplParam,templParamVal,outString);
//     return outString;
// }
// ```
/// `MapTemplate`: every `<@name@>` of `text`, left to right, replaced by `value` (`NStr::Replace`,
/// corelib/ncbistr.cpp:3310-3390: a replaced value is not searched again).
pub fn map_template(text: &[u8], name: &str, value: &[u8]) -> Vec<u8> {
    let tag = format!("<@{name}@>").into_bytes();
    let mut out = Vec::with_capacity(text.len() + value.len());
    let mut at = 0;
    while let Some(found) = text[at..]
        .windows(tag.len())
        .position(|w| w == tag.as_slice())
    {
        out.extend_from_slice(&text[at..at + found]);
        out.extend_from_slice(value);
        at += found + tag.len();
    }
    out.extend_from_slice(&text[at..]);
    out
}

// NCBI reference (2.17.0): c++/src/corelib/ncbistr.cpp:5610-5626 (`s_Encode`, rows 0x20-0x7F)
// ```c++
//     "+",   "!",   "%22", "%23", "$",   "%25", "%26", "'",
//     "(",   ")",   "*",   "%2B", ",",   "-",   ".",   "%2F",
//     ...
//     "X",   "Y",   "Z",   "%5B", "%5C", "%5D", "%5E", "_",
//     "%60", "a",   "b",   "c",   "d",   "e",   "f",   "g",
// ```
/// `NStr::URLEncode` with its default `eUrlEnc_SkipMarkChars` (ncbistr.cpp:6058-6131): letters,
/// digits and `! $ ' ( ) * , - . _` stay, a space becomes `+`, every other byte `%XX` (upper-case
/// hex).
pub fn url_encode(text: &[u8], out: &mut Vec<u8>) {
    const HEX: &[u8; 16] = b"0123456789ABCDEF";
    for &byte in text {
        match byte {
            b'A'..=b'Z' | b'a'..=b'z' | b'0'..=b'9' => out.push(byte),
            b'!' | b'$' | b'\'' | b'(' | b')' | b'*' | b',' | b'-' | b'.' | b'_' => out.push(byte),
            b' ' => out.push(b'+'),
            _ => out.extend_from_slice(&[
                b'%',
                HEX[usize::from(byte >> 4)],
                HEX[usize::from(byte & 15)],
            ]),
        }
    }
}

// NCBI reference (2.17.0): c++/src/objtools/align_format/align_format_util.cpp:4979-4996
// ```c++
// bool CAlignFormatUtil::GetTextSeqID(const list<CRef<CSeq_id> > & ids, string *textSeqID)
// {
//     bool hasTextSeqID = false;
//     CConstRef<CSeq_id> seqID = FindTextseq_id(ids);
//     ...
//     if(seqID.Empty()) {
//         seqID = GetSeq_idByType(ids, CSeq_id::e_Pdb);
//     }
//     if(seqID.Empty()) {
//         seqID = GetSeq_idByType(ids, CSeq_id::e_Patent);
//     }
//     if(!seqID.Empty()) {
//         hasTextSeqID = true;
// ```
/// `GetTextSeqID` of an id list: one of its ids is a Textseq-id, a `pdb` or a patent id.
fn has_text_seq_id(ids: &[ParsedSeqId]) -> bool {
    ids.iter().any(|id| {
        id.textseq().is_some() || matches!(id, ParsedSeqId::Pdb(_) | ParsedSeqId::Patent(_))
    })
}

// NCBI reference (2.17.0): c++/src/objtools/align_format/align_format_util.cpp:3673-3691
// ```c++
// static string s_MapCommonUrlParams(string urlTemplate, CAlignFormatUtil::SSeqURLInfo *seqUrlInfo)
// {
//     string db,logstr_moltype;
//     if(seqUrlInfo->isDbNa) {
//         db = "nucleotide";
//         logstr_moltype = "nucl";
//     } else {
//         db = "protein";
//         logstr_moltype ="prot";
//     }
//     string logstr_location = (seqUrlInfo->isAlignLink) ? "align" : "top";
//     string url_link = CAlignFormatUtil::MapTemplate(urlTemplate,"db",db);
//     url_link = CAlignFormatUtil::MapTemplate(url_link,"gi", GI_TO(TIntId, seqUrlInfo->gi));
//     url_link = CAlignFormatUtil::MapTemplate(url_link,"log",logstr_moltype + logstr_location);
//     url_link = CAlignFormatUtil::MapTemplate(url_link,"blast_rank",seqUrlInfo->blast_rank);
//     url_link = CAlignFormatUtil::MapTemplate(url_link,"rid",seqUrlInfo->rid);
//     url_link = CAlignFormatUtil::MapTemplate(url_link,"acc",seqUrlInfo->accession);
//     url_link = CAlignFormatUtil::MapProtocol(url_link);
//     return url_link;
// }
// ```
/// The Entrez link of a Textseq-id, `pdb` or patent id (`GetIDUrlGen` :3788-3798, without
/// templates: `acc` again, then `cssInf` and `target` empty). The ENTREZ template has no `gi`
/// tag, so the GI is not needed.
fn entrez_url(info: &SeqUrlInfo<'_>, config: &IdUrlConfig<'_>) -> Vec<u8> {
    let (db, mol): (&[u8], &str) = if info.nucleotide {
        (b"nucleotide", "nucl")
    } else {
        (b"protein", "prot")
    };
    let place = match info.place {
        LinkPlace::Description => "top",
        LinkPlace::Heading => "align",
    };
    let mut url = map_template(ENTREZ, "db", db);
    url = map_template(&url, "log", format!("{mol}{place}").as_bytes());
    url = map_template(&url, "blast_rank", info.rank.to_string().as_bytes());
    url = map_template(&url, "rid", rid(info.place));
    url = map_template(&url, "acc", info.accession);
    url = map_template(&url, "protocol", config.protocol);
    url = map_template(&url, "acc", info.accession);
    url = map_template(&url, "cssInf", b"");
    map_template(&url, "target", b"")
}

/// The RID of the links: NCBI's empty `m_Rid` of the description table, `0` of the alignments.
fn rid(place: LinkPlace) -> &'static [u8] {
    match place {
        LinkPlace::Description => b"",
        LinkPlace::Heading => b"0",
    }
}

// NCBI reference (2.17.0): c++/src/objtools/align_format/align_format_util.cpp:3766-3840
// ```c++
// string CAlignFormatUtil::GetIDUrlGen(SSeqURLInfo *seqUrlInfo,const CBioseq::TId* ids)
// {
//     string url_link = NcbiEmptyString;
//     CConstRef<CSeq_id> wid = FindBestChoice(*ids, CSeq_id::WorstRank);
//     bool hasTextSeqID = GetTextSeqID(*ids);
//     ...
//     else if (hasTextSeqID) {
//         ... (Entrez, entrez_url)
//     } else {//seqid general, dbtag specified
//         if(wid->Which() == CSeq_id::e_General){
//             const CDbtag& dtg = wid->GetGeneral();
//             const string& dbname = dtg.GetDb();
//             if(NStr::CompareNocase(dbname, "TI") == 0){
//                 string actual_id = CAlignFormatUtil::GetGnlID(dtg);
//                 ...
//                 else {
//                     string traceUrl = CAlignFormatUtil::MapTagToHTML("TRACE");
//                     url_link = CAlignFormatUtil::MapTemplate(traceUrl,"val",actual_id);
//                     ...
//                     url_link = CAlignFormatUtil::MapTemplate(url_link,"cssInf",(seqUrlInfo->addCssInfo) ? temp_class_info.c_str() : "");
//                     url_link = CAlignFormatUtil::MapTemplate(url_link,"rid",seqUrlInfo->rid);
//                 }
//             }
//         } else if (wid->Which() == CSeq_id::e_Local){
//             ...
//             if (m_Reg) {
//                 user_url = (seqUrlInfo->addCssInfo) ? m_Reg->Get("LOCAL_ID","TOOL_URL_ALIGN") : m_Reg->Get("LOCAL_ID","TOOL_URL");
//             }
//             string id_string;
//             wid->GetLabel(&id_string, CSeq_id::eContent);
//             url_link = CAlignFormatUtil::MapTemplate(user_url,"seq_id", NStr::URLEncode(id_string));
//             url_link = CAlignFormatUtil::MapTemplate(url_link,"db_name", NStr::URLEncode(seqUrlInfo->database));
//             url_link = CAlignFormatUtil::MapTemplate(url_link,"taxid", TAX_ID_TO(int, seqUrlInfo->taxid));
//             ...
//             url_link = CAlignFormatUtil::MapTemplate(url_link,"cssInf",(seqUrlInfo->addCssInfo) ? temp_class_info.c_str() : "");
//             url_link = CAlignFormatUtil::MapTemplate(url_link,"title", id_string);
//             url_link = CAlignFormatUtil::MapTemplate(url_link,"target",seqUrlInfo->new_win ? "TARGET=\"EntrezView\"" : "");
//         }
//     }
//     url_link = CAlignFormatUtil::MapProtocol(url_link);
//     seqUrlInfo->seqUrl = url_link;
// 	return url_link;
// }
// ```
/// `GetIDUrlGen` of the id list `ids` (the Bioseq's for a description row, the definition line's
/// for a heading): the Entrez link of a list with a Textseq-id, `pdb` or patent id; the trace
/// link of a `gnl|TI|n` best id; the `[LOCAL_ID] TOOL_URL` template with its tags filled (no
/// `<a`) for a local best id; else nothing (empty). The applications set no templates, so the
/// WGS branch (only with `useTemplates`) is never taken; `MapProtocol` ends every kind.
pub fn id_url(ids: &[ParsedSeqId], info: &SeqUrlInfo<'_>, config: &IdUrlConfig<'_>) -> Vec<u8> {
    let url = if has_text_seq_id(ids) {
        entrez_url(info, config)
    } else {
        match find_best(ids, text_score) {
            Some(best @ ParsedSeqId::General(dbtag)) if dbtag.db.eq_ignore_ascii_case(b"TI") => {
                // `GetGnlID`: the tag (`format_label`'s form of a `general` id).
                let mut val = Vec::new();
                label::format_label(best, false, &mut val);
                let url = map_template(TRACE, "val", &val);
                let url = map_template(&url, "cssInf", b"");
                map_template(&url, "rid", rid(info.place))
            }
            Some(local @ ParsedSeqId::Local(_)) => local_url(local, info, config),
            _ => Vec::new(),
        }
    };
    map_template(&url, "protocol", config.protocol)
}

/// The `[LOCAL_ID] TOOL_URL` template of a local id `local` with its tags filled
/// (`GetIDUrlGen` :3818-3834); empty without the template.
fn local_url(local: &ParsedSeqId, info: &SeqUrlInfo<'_>, config: &IdUrlConfig<'_>) -> Vec<u8> {
    let template = config.local_id_url.unwrap_or_default();
    let mut id_string = Vec::new();
    label::content(local, true, &mut id_string);
    let mut encoded = Vec::new();
    url_encode(&id_string, &mut encoded);
    let mut url = map_template(template, "seq_id", &encoded);
    encoded.clear();
    url_encode(config.db_name, &mut encoded);
    url = map_template(&url, "db_name", &encoded);
    url = map_template(&url, "taxid", info.taxid.to_string().as_bytes());
    url = map_template(&url, "cssInf", b"");
    url = map_template(&url, "title", &id_string);
    map_template(&url, "target", b"")
}

/// `GetLabel(eContent)` of a description row's display id (showdefline.cpp:463-464): the
/// accession of its link and, without a GI, its score label.
pub fn description_accession(display_id: &ParsedSeqId, out: &mut Vec<u8>) {
    label::content(display_id, true, out);
}

/// `CAlignFormatUtil::GetLabel(seqID, CSeq_id::eContent)` of a heading's id (showalign.cpp:2208,
/// 2266; `eContent` converts to `with_version` true): the accession of its link and, without a
/// GI, the id of the "more sequence titles" block.
pub fn heading_label(id: &ParsedSeqId, out: &mut Vec<u8>) {
    label::format_label(id, true, out);
}

// NCBI reference (2.17.0): c++/src/objtools/align_format/showdefline.cpp:465-471
// ```c++
//         sdl->score_url = !useTemplates ? "<a href=#" : "";
//         if (!useTemplates && m_PositionIndex >= 0) {
//             sdl->score_url += "_" + NStr::IntToString(m_PositionIndex) + "_";
//         }
//         sdl->score_url += sdl->gi == ZERO_GI ? accession :
//             NStr::NumericToString(sdl->gi);
//         sdl->score_url += !useTemplates ? ">" : "";
// ```
/// The score link of a description row: `<a href=#` and the GI, or without one the
/// accession (`description_accession`), then `>`.
pub fn score_link(gi: i64, accession: &[u8], out: &mut Vec<u8>) {
    out.extend_from_slice(b"<a href=#");
    if gi == 0 {
        out.extend_from_slice(accession);
    } else {
        out.extend_from_slice(gi.to_string().as_bytes());
    }
    out.push(b'>');
}

// NCBI reference (2.17.0): c++/src/objtools/align_format/showalign.cpp:2357-2359, 2443-2449
// ```c++
//                 if(alnDispParams->gi != ZERO_GI){
//                     out<<"<a name="<<alnDispParams->gi<<"></a>";
//                 } else {
//                     out<<"<a name="<<alnDispParams->seqID->GetSeqIdString()<<"></a>";
//     ...
//                             out<<"<a name="<<alnDispParams->seqID->GetSeqIdString(true)<<"></a>";
// ```
/// The anchor of a heading line: `<a name=` the GI, or without one `GetSeqIdString` of the id,
/// with its version for a database's definition line and without for a `-subject` Bioseq
/// (`with_version`), then `></a>`.
pub fn anchor(gi: i64, id: &ParsedSeqId, with_version: bool, out: &mut Vec<u8>) {
    out.extend_from_slice(b"<a name=");
    if gi == 0 {
        label::content(id, with_version, out);
    } else {
        out.extend_from_slice(gi.to_string().as_bytes());
    }
    out.extend_from_slice(b"></a>");
}

#[cfg(test)]
#[path = "html_ids_tests.rs"]
mod tests;
