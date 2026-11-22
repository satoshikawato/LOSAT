//! The ids of `-html` reports: the link of a sequence's id (NCBI's `CAlignFormatUtil::GetIDUrlGen`),
//! the label of its score link and the anchor of its alignment heading, and the run's and the
//! report's `-html` data that they read (`HtmlLinks`, `HtmlReport`). Reached only with `-html`:
//! the plain report never classifies an id here.

use std::io::{self, Write};
use std::sync::Arc;

use super::formatter_registry::{map_template, InitConfig, ReferenceRegistry};
use crate::report::db_subject::ParsedSubject;
use crate::seqdb::Oid;
use crate::seqid::label;
use crate::seqid::model::{Dbtag, ObjectId, ParsedSeqId};
use crate::seqid::rank::{find_best, text_score};

/// The `-html` links of a run (`PairwiseConfig::html` holds them through `HtmlReport`): NCBI's
/// formatter registry (`CAlignFormatUtil::InitConfig`, `m_Reg`, `m_Protocol`: statics, one per
/// run), read once at the program's setup and only with `-html` and outfmt 0, whose ids link
/// (`blastinput::app::html_links`, so that LOSAT's refusal of the file comes before the first
/// byte; outfmt 18 takes no links); the `-db` text as NCBI keeps
/// it (`CBlastFormat::m_DbName`, `NStr::TruncateSpaces` of `-db`; empty for `-subject`), the
/// `<@db_name@>` of a local id's link; whether the database holds nucleotides (`m_IsDbNa`); and
/// the `./.ncbirc` of the prolog's reference links (`ReferenceRegistry`, read at the same setup;
/// `pairwise::write_prolog`).
#[derive(Debug)]
pub struct HtmlLinks {
    pub(crate) registry: InitConfig,
    pub(crate) db_name: Vec<u8>,
    pub(crate) nucleotide: bool,
    pub(crate) references: ReferenceRegistry,
}

/// The `-html` data of one pairwise report (`PairwiseConfig::html`; `None` is the plain report):
/// the run's links and the OID of each database sequence that the report can show, by subject
/// index, which names a sequence without parsed Seq-ids (`gnl|BL_ORD_ID|<oid>`); empty for a
/// `-subject` report, whose records without parsed Seq-ids are `lcl|Subject_<index + 1>`. The
/// one OID source of the five programs' reports: each program hands its reported sequences'
/// OIDs (`HtmlReport::new`) where it makes the report's `PairwiseConfig`, only with `-html`.
#[derive(Clone, Debug)]
pub struct HtmlReport {
    pub(crate) links: Arc<HtmlLinks>,
    oids: Vec<Oid>,
}

impl HtmlReport {
    /// The report of the run's links `links` whose database sequences have the OIDs `oids` by
    /// subject index (none for `-subject`).
    pub(crate) fn new(links: &Arc<HtmlLinks>, oids: impl IntoIterator<Item = Oid>) -> Self {
        HtmlReport {
            links: Arc::clone(links),
            oids: oids.into_iter().collect(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:459-494
    /// The link (`sdl->id_url`, empty for none) and the score link (`sdl->score_url`) of the
    /// description row of rank `rank` (from 1, NCBI's `blast_rank`) for the subject of index
    /// `s_idx`, `parsed` when it has parsed Seq-ids: the link of its Bioseq's ids (`GetIDUrl`
    /// with the handle of the Seq-align id) with the display id's label (`GetLabel`, `eContent`) as
    /// the accession; without parsed ids a database sequence is `gnl|BL_ORD_ID|<oid>` and a
    /// `-subject` record `lcl|Subject_<s_idx + 1>`. The registry's posts, its `GETURL_DEBUG` line
    /// and its exception go to standard error at the first link of the run (`id_url`).
    pub(crate) fn description_links(
        &self,
        parsed: Option<&ParsedSubject>,
        s_idx: u32,
        rank: usize,
    ) -> io::Result<(Vec<u8>, Vec<u8>)> {
        let own;
        let (ids, display, gi, taxid) = match parsed {
            Some(subject) => match subject.description_ids() {
                Some((display, gi)) => {
                    let taxid = match display {
                        ParsedSeqId::Local(_) => subject.first_taxid(),
                        _ => 0,
                    };
                    (subject.bioseq_ids(), display, gi, taxid)
                }
                None => return Ok((Vec::new(), Vec::new())),
            },
            None => {
                own = [self.unparsed_id(s_idx)];
                (&own[..], &own[0], 0, 0)
            }
        };
        let mut accession = Vec::new();
        description_accession(display, &mut accession);
        let info = SeqUrlInfo {
            accession: &accession,
            rank,
            nucleotide: self.links.nucleotide,
            place: LinkPlace::Description,
            taxid,
        };
        let link = id_url(ids, &info, &self.links, &mut io::stderr())?;
        let mut score = Vec::new();
        score_link(gi, &accession, &mut score);
        Ok((link, score))
    }

    /// The id of the subject of index `s_idx` without parsed Seq-ids: `gnl|BL_ORD_ID|<oid>` of a
    /// database sequence, `lcl|Subject_<s_idx + 1>` of a `-subject` record.
    fn unparsed_id(&self, s_idx: u32) -> ParsedSeqId {
        match self.oids.get(s_idx as usize) {
            Some(oid) => ParsedSeqId::General(Dbtag {
                db: b"BL_ORD_ID".to_vec(),
                tag: ObjectId::Id(i32::try_from(oid.get()).expect("an OID is NCBI's int")),
            }),
            None => ParsedSeqId::Local(ObjectId::Str(
                format!("Subject_{}", u64::from(s_idx) + 1).into_bytes(),
            )),
        }
    }
}

/// The `-html` data of one subject's alignment heading (`x_PrintDefLine`): the report's
/// (`HtmlReport`), the subject's index `s_idx` and the heading's rank `rank` (NCBI's `m_cur_align`,
/// from 1: `QueryGroup::heading_rank`), the `blast_rank` of its links.
#[derive(Clone, Copy, Debug)]
pub(crate) struct HtmlHeading<'a> {
    pub(crate) report: &'a HtmlReport,
    pub(crate) s_idx: u32,
    pub(crate) rank: usize,
}

impl HtmlHeading<'_> {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2413-2449
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2209-2221, 2265-2271
    /// The id of one line of the heading after its `>` or space: the link of the id list `ids` (a
    /// definition line's, or a `-subject` Bioseq's) with `heading_label` of `display` (their
    /// `WorstRank`-best id) as the accession, the line's taxid `taxid` (0 for a `-subject` Bioseq)
    /// and the heading's rank, then the plain id text (`id_text`), `</a>` after a link, and the
    /// anchor (`anchor` of the line's GI `gi`, with the version for a definition line,
    /// `with_version`). Returns the line's id label (`aln_vec_info->id_label`: the GI, else the
    /// accession), the id of the "more sequence titles" block (`more_titles`). The registry's
    /// exception (`id_url`) is the `io::Error`, which NCBI's `DisplaySeqalign` catches for each
    /// alignment (`skips`).
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn line_id(
        &self,
        ids: &[ParsedSeqId],
        display: &ParsedSeqId,
        gi: i64,
        taxid: i32,
        with_version: bool,
        id_text: impl FnOnce(&mut Vec<u8>),
        out: &mut Vec<u8>,
    ) -> io::Result<Vec<u8>> {
        let mut accession = Vec::new();
        heading_label(display, &mut accession);
        let info = SeqUrlInfo {
            accession: &accession,
            rank: self.rank,
            nucleotide: self.report.links.nucleotide,
            place: LinkPlace::Heading,
            taxid,
        };
        let link = id_url(ids, &info, &self.report.links, &mut io::stderr())?;
        out.extend_from_slice(&link);
        id_text(out);
        if !link.is_empty() {
            out.extend_from_slice(b"</a>");
        }
        anchor(gi, display, with_version, out);
        Ok(if gi == 0 {
            accession
        } else {
            gi.to_string().into_bytes()
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2312-2361
    /// The first line of the heading of a subject without parsed Seq-ids up to its title: `>`, the
    /// link of its id (`unparsed_id`: none for `gnl|BL_ORD_ID|<oid>`, the registry's `[LOCAL_ID]
    /// TOOL_URL` for `lcl|Subject_<s_idx + 1>`) and `</a>` after one, the anchor
    /// (`<a name=BL_ORD_ID:<oid>></a>`, `<a name=Subject_<s_idx + 1>></a>`; such an id has no
    /// version, so a definition line's `GetSeqIdString` with versions and a `-subject` Bioseq's
    /// `GetSeqIdString` without them agree) and a space; neither id is shown. The registry's
    /// exception as `line_id`'s.
    pub(crate) fn unparsed_start(&self, out: &mut Vec<u8>) -> io::Result<()> {
        let id = self.report.unparsed_id(self.s_idx);
        out.push(b'>');
        self.line_id(std::slice::from_ref(&id), &id, 0, 0, true, |_| {}, out)?;
        out.push(b' ');
        Ok(())
    }

    /// `GetSeqIdString()` of the subject without parsed Seq-ids (`BL_ORD_ID:<oid>`,
    /// `Subject_` and the subject's number, `s_idx` plus one): the id of NCBI's line in place of
    /// each of its alignments where the registry throws (`skips`).
    pub(crate) fn unparsed_skipped_id(&self) -> Vec<u8> {
        let mut id = Vec::new();
        label::content(&self.report.unparsed_id(self.s_idx), false, &mut id);
        id
    }

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1954-1981
    /// Whether the error `error` of a heading's id link is the registry's exception, which NCBI
    /// catches for each alignment of the subject (its line in place of the alignment, with
    /// `skipped_id`): `Ok` with that id; any other error stays the report's.
    pub(crate) fn skips(
        error: io::Error,
        skipped_id: impl FnOnce() -> Vec<u8>,
    ) -> io::Result<Vec<u8>> {
        match crate::ncbi_exception::NcbiException::of_io_error(&error) {
            Some(_) => Ok(skipped_id()),
            None => Err(error),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2389-2397, 2474-2476
/// The lines a heading shows before it hides the rest (`k_MinDeflinesToShow`), where it has more
/// than `MAX_SHOWN_LINES` (`k_MaxDeflinesToShow`).
pub(crate) const SHOWN_LINES: usize = 3;
pub(crate) const MAX_SHOWN_LINES: usize = 8;

/// The "more sequence titles" link of a heading that hides `hidden` lines, after the space of the
/// first hidden line: its id is the label of the last shown line (`id_label`, `line_id`), and the
/// hidden lines follow in the `info_<id_label>` box, which `</div>` closes before `Length=`.
pub(crate) fn more_titles(id_label: &[u8], hidden: usize, out: &mut Vec<u8>) {
    out.extend_from_slice(b"<a href=\"#\" title=\"Other sequence titles\"  onmouseover=\"showInfo(this)\" class=\"resArrowLinkW mdl hiding\" id=\"");
    out.extend_from_slice(id_label);
    out.extend_from_slice(b"\">");
    out.extend_from_slice(hidden.to_string().as_bytes());
    out.extend_from_slice(b" more sequence titles</a>\n <div id=\"info_");
    out.extend_from_slice(id_label);
    out.extend_from_slice(b"\" class=\"helpbox mdlbox hidden\">");
}

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
    /// `accession`: the label of the description's display id (`GetLabel`, `eContent`), or the
    /// heading's (`CAlignFormatUtil::GetLabel` with versions; `description_accession`,
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

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:181,187
const ENTREZ: &[u8] = b"<a title=\"Show report for <@acc@>\" <@cssInf@>href=\"https://www.ncbi.nlm.nih.gov/<@db@>/<@acc@>?report=genbank&log$=<@log@>&blast_rank=<@blast_rank@>&RID=<@rid@>\" <@target@>>";
const TRACE: &[u8] = b"<a title=\"Show report for <@val@>\" <@cssInf@>href=\"https://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=retrieve&dopt=fasta&val=<@val@>&RID=<@rid@>\">";

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:4979-4996
/// `GetTextSeqID` of an id list: one of its ids is a Textseq-id, a `pdb` or a patent id.
fn has_text_seq_id(ids: &[ParsedSeqId]) -> bool {
    ids.iter().any(|id| {
        id.textseq().is_some() || matches!(id, ParsedSeqId::Pdb(_) | ParsedSeqId::Patent(_))
    })
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3673-3691
/// The Entrez link of a Textseq-id, `pdb` or patent id (`GetIDUrlGen` :3788-3798, without
/// templates: `acc` again, then `cssInf` and `target` empty). The ENTREZ template has no `gi`
/// tag, so the GI is not needed. `s_MapCommonUrlParams` ends with the registry's `MapProtocol`
/// (`InitConfig::map_protocol`), whose first call of the run reads the registry.
fn entrez_url(
    info: &SeqUrlInfo<'_>,
    links: &HtmlLinks,
    diagnostics: &mut dyn Write,
) -> io::Result<Vec<u8>> {
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
    url = links.registry.map_protocol(&url, diagnostics)?;
    url = map_template(&url, "acc", info.accession);
    url = map_template(&url, "cssInf", b"");
    Ok(map_template(&url, "target", b""))
}

/// The RID of the links: NCBI's empty `m_Rid` of the description table, `0` of the alignments.
fn rid(place: LinkPlace) -> &'static [u8] {
    match place {
        LinkPlace::Description => b"",
        LinkPlace::Heading => b"0",
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:3766-3840
/// `GetIDUrlGen` of the id list `ids` (the Bioseq's for a description row, the definition line's
/// for a heading): the Entrez link of a list with a Textseq-id, `pdb` or patent id; the trace
/// link of a `gnl|TI|n` best id; the registry's `[LOCAL_ID] TOOL_URL` with its tags filled
/// (`InitConfig::local_id_url`, no `<a`) for a local best id (`Subject_N` too); else nothing
/// (empty). The applications set no templates, so the WGS branch (only with `useTemplates`) is
/// never taken. `MapProtocol` ends every kind, so the registry is read at the first id of the run
/// whatever its kind: its posts and the `GETURL_DEBUG` line go to `diagnostics` (standard error)
/// then, and a malformed file's exception is the `io::Error` of every call (§1.9 of SURVEY_T4).
pub fn id_url(
    ids: &[ParsedSeqId],
    info: &SeqUrlInfo<'_>,
    links: &HtmlLinks,
    diagnostics: &mut dyn Write,
) -> io::Result<Vec<u8>> {
    let url = if has_text_seq_id(ids) {
        entrez_url(info, links, diagnostics)?
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
            Some(local @ ParsedSeqId::Local(_)) => {
                // The best id's `GetLabel` in its `eContent` form; `local_id_url` ends with
                // `MapProtocol`.
                let mut id_string = Vec::new();
                label::content(local, true, &mut id_string);
                return links.registry.local_id_url(
                    &id_string,
                    &links.db_name,
                    info.taxid,
                    diagnostics,
                );
            }
            _ => Vec::new(),
        }
    };
    links.registry.map_protocol(&url, diagnostics)
}

/// The label of a description row's display id (`GetLabel`, `eContent`; showdefline.cpp:463-464):
/// the accession of its link and, without a GI, its score label.
pub fn description_accession(display_id: &ParsedSeqId, out: &mut Vec<u8>) {
    label::content(display_id, true, out);
}

/// `CAlignFormatUtil::GetLabel` of a heading's id with `CSeq_id::eContent` (showalign.cpp:2208,
/// 2266; `eContent` converts to `with_version` true): the accession of its link and, without a
/// GI, the id of the "more sequence titles" block.
pub fn heading_label(id: &ParsedSeqId, out: &mut Vec<u8>) {
    label::format_label(id, true, out);
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:465-471
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

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2357-2359, 2443-2449
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
