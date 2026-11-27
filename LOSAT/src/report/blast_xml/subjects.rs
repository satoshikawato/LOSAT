//! How a program's final hit lists name and measure their subjects in the XML (`Hit_id`,
//! `Hit_def`, `Hit_accession`, `Hit_len`; DESIGN_OF4 2.3).

use std::sync::Arc;

use anyhow::{Context, Result};

use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::db_subject::{HitDescription, ParsedSubject};

use super::names::{self, Refusal};
use super::records::HitTaxonomy;

/// The subjects of a search, indexed by the subject number of its HSPs.
#[derive(Clone, Copy)]
pub enum XmlSubjects<'a> {
    /// A `-subject` set: the records as searched (`records[s_idx]`), named `Subject_N` with
    /// their titles; `length(s_idx)` is the record's own length (the whole record under
    /// `-subject_loc`). `protein`: whether the records are proteins (the definition line of a
    /// protein record is made differently).
    Records {
        records: &'a [FastaRecord],
        length: &'a dyn Fn(usize) -> usize,
        protein: bool,
    },
    /// A database: the sequences with hits by rank (`ordinals[rank]` the ordinal id, whose OID
    /// `oids` gives, the stored title `titles.get(rank)`, `length(rank)` the sequence length).
    Database {
        ordinals: &'a [u32],
        oids: crate::seqdb::OidMap<'a>,
        titles: XmlTitles<'a>,
        /// The Blast-def-line taxid of each sequence by rank, 0 for none; only XML2 shows
        /// it. Empty where the program has not read them (outfmt 5 alone): every taxid is
        /// then 0.
        taxids: &'a [i32],
        /// What the taxonomy database said about each sequence's taxid by rank
        /// (`HitTaxonomy`); empty where nothing was asked: every one is `Unnamed`.
        taxonomy: &'a [HitTaxonomy],
        length: &'a dyn Fn(usize) -> usize,
        /// The subjects with parsed Seq-ids by rank (empty without them), whose names come from
        /// their definition lines (`SubjectNames::of_parsed`).
        parsed: &'a [Option<Arc<ParsedSubject>>],
        /// What the taxonomy database said about the taxid of each definition line of a subject
        /// with parsed Seq-ids, by rank (empty where nothing was asked;
        /// `algorithm::common::database::read_reported_taxonomy`).
        parsed_taxonomy: &'a [Vec<HitTaxonomy>],
        /// What the reports do at the `Hit` of each sequence by rank (`db_subject::HitLoads`: one
        /// whose `Hit` throws, one that no report loads), for `hit_names`.
        loads: &'a crate::report::db_subject::HitLoads,
    },
}

/// The stored titles of a database's sequences by rank.
#[derive(Clone, Copy)]
pub enum XmlTitles<'a> {
    /// The records of the listed sequences (`records[rank].title`).
    Records(&'a [FastaRecord]),
    /// The titles alone (blastp's `LOSAT_X_BPLISTEDLEAN`: the sequences in no final hit list
    /// share one empty title).
    Shared(&'a [Arc<[u8]>]),
}

impl<'a> XmlTitles<'a> {
    fn get(self, rank: usize) -> Option<&'a [u8]> {
        match self {
            XmlTitles::Records(records) => records.get(rank).map(|record| &record.title[..]),
            XmlTitles::Shared(titles) => titles.get(rank).map(|title| &title[..]),
        }
    }

    fn len(self) -> usize {
        match self {
            XmlTitles::Records(records) => records.len(),
            XmlTitles::Shared(titles) => titles.len(),
        }
    }
}

/// The names of a `Hit`.
pub struct SubjectNames {
    pub id: Vec<u8>,
    /// The title as generated (`-subject`) or stored (`-db`), not yet replaced by
    /// `No definition line` when empty (`names::NO_DEFINITION_LINE`).
    pub title: Vec<u8>,
    pub accession: Vec<u8>,
    /// The `-db` sequence's Blast-def-line taxid, 0 for none and for a `-subject` record.
    pub taxid: i32,
    /// What the taxonomy database said about `taxid` (XML2 and JSON only).
    pub taxonomy: HitTaxonomy,
    /// The XML2 and JSON `HitDescr` list of a subject with parsed Seq-ids, one per definition
    /// line (empty for the others, which have one made of `id`, `title` and `accession`).
    pub descriptions: Vec<HitDescription>,
    /// The taxids that the volume looks up when it fetches a subject with parsed Seq-ids
    /// (`ParsedSubject::volume_taxids`; empty for the others, whose one is `taxid`).
    pub volume_taxids: Vec<i32>,
    /// Outfmt 5's `catch` names the `Hit` (`of_unknown`): `GenerateDefline` threw for the
    /// `-subject` record, or the walk of a database sequence's definition lines threw at one
    /// without a Seq-id (`of_parsed`). Outfmt 5 writes `Hit_len` 0.
    pub defline_thrown: bool,
    /// NCBI's scope gets no Bioseq for the database sequence
    /// (`ParsedSubject::scope_finds_none`): the names are the `catch`'s (`of_unknown`), and so is
    /// each `Hsp` (`writer::write_caught_hsp`).
    pub unloaded: bool,
    /// The exception that NCBI's lookup of the database sequence's Seq-align id throws where its
    /// `Hit` loads it (`of_thrown`, `db_subject::HitLoad::Throws`): the writer of each document
    /// decides what the `Hit` is (`XmlHit::lookup_throws`).
    pub lookup_throws: Option<crate::ncbi_exception::NcbiException>,
}

impl SubjectNames {
    /// The names of a `-subject` record `id` with the FASTA title `title` (`Subject_N` and the
    /// generated definition line), or where `GenerateDefline` throws the catch's
    /// (`SubjectNames::of_unknown` of `lcl|Subject_N`). `protein`: the record is a protein.
    pub fn of_record(id: &[u8], title: &[u8], protein: bool) -> Result<Self, Refusal> {
        let Some((xml_id, title)) = names::xml_subject_names(id, title, protein)? else {
            return Ok(Self::of_unknown(names::fasta_id(id, None)));
        };
        Ok(SubjectNames {
            id: xml_id,
            title,
            accession: id.to_vec(),
            taxid: 0,
            taxonomy: HitTaxonomy::Unnamed,
            descriptions: Vec::new(),
            volume_taxids: Vec::new(),
            defline_thrown: false,
            unloaded: false,
            lookup_throws: None,
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:399-407
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:412-417
    /// The names of a `-subject` record whose `GenerateDefline` threw (`NStr::HtmlDecode` cannot
    /// guess its title's encoding), as NCBI's catch writes them: outfmt 5's `Hit_id`
    /// `kSeqId.AsFastaString()` (`fasta_id`), `Hit_def` and `Hit_accession` `Unknown` and `Hit_len`
    /// 0 (`defline_thrown`); XML2's and JSON's one `HitDescr` of that id alone, with the record's
    /// length.
    fn of_unknown(fasta_id: Vec<u8>) -> Self {
        SubjectNames {
            descriptions: vec![HitDescription {
                id: fasta_id.clone(),
                accession: Vec::new(),
                title: None,
                taxid: 0,
                taxonomy: HitTaxonomy::Unnamed,
            }],
            id: fasta_id,
            title: b"Unknown".to_vec(),
            accession: b"Unknown".to_vec(),
            taxid: 0,
            taxonomy: HitTaxonomy::Unnamed,
            volume_taxids: Vec::new(),
            defline_thrown: true,
            unloaded: false,
            lookup_throws: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:376-407
    /// The names of a database sequence whose lookup of its Seq-align id throws where its `Hit`
    /// loads it (`exception`, `db_subject::HitLoad::Throws`): those of NCBI's `catch` of outfmt 5
    /// (`of_unknown`'s, of the Seq-align id's FASTA form, `ParsedSubject::alignment_fasta`), with
    /// the exception, which the XML2 and JSON documents raise instead (`lookup_throws`).
    pub fn of_thrown(
        subject: &ParsedSubject,
        exception: crate::ncbi_exception::NcbiException,
    ) -> Self {
        SubjectNames {
            defline_thrown: false,
            lookup_throws: Some(exception),
            ..Self::of_unknown(subject.alignment_fasta())
        }
    }

    /// The names of a `-subject` record whose ids `-parse_deflines` read (`ParsedSubject::of_record`):
    /// its id list, its generated title and the accession of its best id; those of `of_record` for a
    /// record shown as one without the option (`ParsedSubject::shown_as_record`). Like a record
    /// without the option it has one `HitDescr` of these (no definition lines).
    pub fn of_believed(
        ids: &crate::blastinput::fasta_reader::believed::BelievedIds,
        title: &[u8],
        protein: bool,
    ) -> Result<Self, Refusal> {
        if ids.subject_local() {
            return Self::of_record(ids.bare_id().as_bytes(), title, protein);
        }
        let subject = ParsedSubject::of_record(ids, title);
        let (id, accession) = subject.record_xml_ids();
        let Ok(title) = subject.bioseq_title(protein, false) else {
            return Ok(Self::of_unknown(names::fasta_id(b"", Some(ids))));
        };
        if title.reads_past_end {
            return Err(Refusal::ReadsPastEnd);
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:393-406
        // The XML2 `HitDescr` of a record whose best id is not local has its accession too.
        let descriptions = if subject.has_local_best_id() {
            Vec::new()
        } else {
            vec![HitDescription {
                id: id.clone(),
                accession: accession.clone(),
                title: Some(title.text.clone()),
                taxid: 0,
                taxonomy: HitTaxonomy::Unnamed,
            }]
        };
        Ok(SubjectNames {
            id,
            title: title.text,
            accession,
            taxid: 0,
            taxonomy: HitTaxonomy::Unnamed,
            descriptions,
            volume_taxids: Vec::new(),
            defline_thrown: false,
            unloaded: false,
            lookup_throws: None,
        })
    }

    /// The names of the database sequence `oid` with the stored title `title`
    /// (`gnl|BL_ORD_ID|n`, the title as stored: showdefline.cpp:265-303 reads the Blast-def-line's
    /// title).
    pub fn of_database_sequence(
        oid: crate::seqdb::Oid,
        title: &[u8],
        taxid: i32,
        taxonomy: HitTaxonomy,
    ) -> Self {
        SubjectNames {
            id: oid.fasta_id().into_bytes(),
            title: title.to_vec(),
            accession: oid.get().to_string().into_bytes(),
            taxid,
            taxonomy,
            descriptions: Vec::new(),
            volume_taxids: Vec::new(),
            defline_thrown: false,
            unloaded: false,
            lookup_throws: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:370-371
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:380-415,422-428
    /// The names of a database sequence with parsed Seq-ids (`ParsedSubject::xml_hit` and
    /// `hit_descriptions`): `Hit_id`, `Hit_def` and `Hit_accession` of outfmt 5, and the
    /// `HitDescr` of each definition line of XML2 and JSON with what the taxonomy database said
    /// about its taxid (`taxonomy[i]` for definition line `i`; empty where nothing was asked).
    /// The walk of the lines throws at a line without a Seq-id (the handle of its null best id,
    /// before its taxid is looked up) and where `GetTaxInfo` does not know a taxid: the
    /// descriptions before that definition line stay, and the Seq-align id alone follows them.
    /// Where outfmt 5's walk throws (`xml_hit` gives `None`), its `Hit` is the `catch`'s
    /// (`of_unknown`: `Hit_len` 0), with the same XML2 descriptions and the real length there.
    /// A `Hit` whose lookup throws has the names of `of_thrown` (`HitNames::of`).
    pub fn of_parsed(subject: &ParsedSubject, taxonomy: &[HitTaxonomy]) -> Self {
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:399-407
        // An empty Bioseq handle throws at the first read of the `Hit`: the `catch`'s names.
        if subject.scope_finds_none() {
            return SubjectNames {
                unloaded: true,
                ..Self::of_unknown(subject.alignment_fasta())
            };
        }
        let mut descriptions = subject.hit_descriptions();
        let lineless = subject
            .deflines()
            .iter()
            .position(|defline| defline.ids.is_empty());
        let mut stop = lineless.map(|at| (at, HitTaxonomy::Unnamed));
        for at in 0..lineless.unwrap_or(descriptions.len()) {
            match taxonomy.get(at) {
                Some(answer @ HitTaxonomy::NotFound(_)) => {
                    stop = Some((at, answer.clone()));
                    break;
                }
                Some(answer) => descriptions[at].taxonomy = answer.clone(),
                None => {}
            }
        }
        if let Some((at, answer)) = stop {
            descriptions.truncate(at);
            descriptions.push(HitDescription {
                id: subject.alignment_fasta(),
                accession: Vec::new(),
                title: None,
                taxid: 0,
                taxonomy: answer,
            });
        }
        let volume_taxids = subject.volume_taxids();
        let Some((id, title, accession)) = subject.xml_hit() else {
            return SubjectNames {
                descriptions,
                volume_taxids,
                ..Self::of_unknown(subject.alignment_fasta())
            };
        };
        SubjectNames {
            id,
            title,
            accession,
            taxid: 0,
            taxonomy: HitTaxonomy::Unnamed,
            descriptions,
            volume_taxids,
            defline_thrown: false,
            unloaded: false,
            lookup_throws: None,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:357-360, 412-417
// (cited at `report::db_subject::ParsedSubject::lookup_throws`: the `Hit` loads its subject by
// the Seq-align id, and a lookup that throws there throws again outside the hit's `catch`)
/// What a `Hit` shows of its subject, made before the report (`SubjectNameCache`, and TBLASTN's
/// and BLASTX's names of the subjects with hits), by what the reports do at it
/// (`db_subject::HitLoads`, the one answer): the subject's names, those of NCBI's `catch` with the
/// exception of the lookup of its Seq-align id where the `Hit` loads it (`SubjectNames::of_thrown`,
/// which each document's writer writes or raises, `XmlHit::lookup_throws`); or nothing for a
/// subject that no report loads (its names are not made, so that nothing NCBI never reads is
/// refused). The names are boxed: a `Hit` that no report loads holds nothing.
pub enum HitNames {
    Names(Box<SubjectNames>),
    Unread,
}

impl HitNames {
    /// The names of the subject `key` of `loads` as its `Hit` shows them: those of NCBI's `catch`
    /// for a `Hit` that throws (`HitLoad::Throws`, `SubjectNames::of_thrown` of `parsed`, the
    /// subject with parsed Seq-ids, which every subject whose lookup throws is), nothing for a
    /// subject that no report loads (`HitLoad::Unread`), else the names `make` makes, whose error
    /// (LOSAT's refusal of a title) is returned.
    pub fn of(
        loads: &crate::report::db_subject::HitLoads,
        key: usize,
        parsed: Option<&ParsedSubject>,
        make: impl FnOnce() -> Result<SubjectNames>,
    ) -> Result<Self> {
        use crate::report::db_subject::HitLoad;
        Ok(match (loads.load(key), loads.thrown(key)) {
            (HitLoad::Throws, Some(exception)) => {
                let subject = parsed.context(
                    "LOSAT internal error: a Hit whose lookup throws has no parsed Seq-ids",
                )?;
                HitNames::Names(Box::new(SubjectNames::of_thrown(
                    subject,
                    exception.clone(),
                )))
            }
            (HitLoad::Unread, _) => HitNames::Unread,
            _ => HitNames::Names(Box::new(make()?)),
        })
    }

    /// The names that the `Hit` shows. The `Hit` of a subject that no report loads is LOSAT's
    /// internal error: the hit loop ends at the `Hit` that throws before it (outfmt 13 to 16,
    /// whose writers raise its exception), as `HitLoads` does.
    pub fn names(&self) -> Result<&SubjectNames> {
        match self {
            HitNames::Names(names) => Ok(names),
            HitNames::Unread => {
                anyhow::bail!("LOSAT: the Hit of a subject that the reports do not load")
            }
        }
    }
}

impl XmlSubjects<'_> {
    /// The names of subject `s_idx`; `program` (`BLASTP`) names the program in the refusal of
    /// a title that NCBI crashes on (`Refusal`).
    pub fn names(&self, s_idx: usize, program: &str) -> Result<SubjectNames> {
        match self {
            XmlSubjects::Records {
                records, protein, ..
            } => {
                let record = records.get(s_idx).context("HSP of an unknown subject")?;
                match record.believed() {
                    Some(ids) => SubjectNames::of_believed(ids, &record.title, *protein),
                    None => {
                        SubjectNames::of_record(record.local_id.as_bytes(), &record.title, *protein)
                    }
                }
                .map_err(|refusal| refusal.error(program, "a subject title"))
            }
            XmlSubjects::Database {
                ordinals,
                oids,
                titles,
                taxids,
                taxonomy,
                parsed,
                parsed_taxonomy,
                ..
            } => {
                if let Some(subject) = parsed.get(s_idx).and_then(Option::as_deref) {
                    let answers = parsed_taxonomy.get(s_idx).map_or(&[][..], Vec::as_slice);
                    return Ok(SubjectNames::of_parsed(subject, answers));
                }
                let ordinal = *ordinals
                    .get(s_idx)
                    .context("HSP of an unknown database sequence")?;
                let title = titles
                    .get(s_idx)
                    .context("HSP of an unknown database sequence")?;
                let taxid = taxids.get(s_idx).copied().unwrap_or(0);
                let taxonomy = taxonomy.get(s_idx).cloned().unwrap_or_default();
                Ok(SubjectNames::of_database_sequence(
                    oids.oid(ordinal),
                    title,
                    taxid,
                    taxonomy,
                ))
            }
        }
    }

    /// The number of subjects that `names` can be asked for.
    pub fn count(&self) -> usize {
        match self {
            XmlSubjects::Records { records, .. } => records.len(),
            XmlSubjects::Database {
                ordinals, titles, ..
            } => ordinals.len().min(titles.len()),
        }
    }

    /// `Hit_len` of subject `s_idx`.
    pub fn length(&self, s_idx: usize) -> i32 {
        match self {
            XmlSubjects::Records { length, .. } | XmlSubjects::Database { length, .. } => {
                length(s_idx) as i32
            }
        }
    }

    /// `names` of subject `s_idx` as its `Hit` shows them (`HitNames::of`, by the database's
    /// `loads`; the records of a subject set are always loaded).
    pub fn hit_names(&self, s_idx: usize, program: &str) -> Result<HitNames> {
        match self {
            XmlSubjects::Records { .. } => self
                .names(s_idx, program)
                .map(|names| HitNames::Names(Box::new(names))),
            XmlSubjects::Database { loads, parsed, .. } => HitNames::of(
                loads,
                s_idx,
                parsed.get(s_idx).and_then(Option::as_deref),
                || self.names(s_idx, program),
            ),
        }
    }

    /// Whether the HSPs of subject `s_idx`, which follows subject `previous` in a hit list, are
    /// shown in `previous`'s `Hit`: copies of one record of a database
    /// (`db_subject::joins_as_copy`).
    pub fn joins(&self, previous: usize, s_idx: usize) -> bool {
        match self {
            XmlSubjects::Records { .. } => false,
            XmlSubjects::Database { parsed, .. } => {
                crate::report::db_subject::rank_joins_as_copy(parsed, previous, s_idx)
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:465-487
// (XML2 and JSON: blastxml2_format.cpp:480-497, the same.)
/// The HSP lists of a query's hit list `lists` that the XML, XML2 and JSON reports show as one
/// `Hit` each, in order: ranges of the lists, each of a subject and the copies that follow it
/// (`joins`, `XmlSubjects::joins`; a list of the same subject as the one before goes on its
/// `Hit`, as an id matches itself). `subject_of` gives the subject of a list, `None` for a list
/// without HSPs, which shows nothing (it can lie inside a range).
pub fn hit_ranges<T>(
    lists: &[T],
    subject_of: impl Fn(&T) -> Option<usize>,
    joins: impl Fn(usize, usize) -> bool,
) -> Vec<std::ops::Range<usize>> {
    let mut ranges: Vec<std::ops::Range<usize>> = Vec::new();
    let mut previous: Option<usize> = None;
    for (index, list) in lists.iter().enumerate() {
        let Some(s_idx) = subject_of(list) else {
            continue;
        };
        match (ranges.last_mut(), previous) {
            (Some(range), Some(previous)) if previous == s_idx || joins(previous, s_idx) => {
                range.end = index + 1
            }
            _ => ranges.push(index..index + 1),
        }
        previous = Some(s_idx);
    }
    ranges
}

/// The names of the subjects of a report, made once per subject (the definition line of a
/// subject is generated when it first has a hit, not again for each hit or each query that hits
/// it).
pub struct SubjectNameCache<'s, 'a> {
    subjects: &'s XmlSubjects<'a>,
    names: Vec<Option<HitNames>>,
}

impl<'s, 'a> SubjectNameCache<'s, 'a> {
    pub fn new(subjects: &'s XmlSubjects<'a>) -> Self {
        let mut names = Vec::new();
        names.resize_with(subjects.count(), || None);
        Self { subjects, names }
    }

    /// `XmlSubjects::hit_names` of subject `s_idx`, made on first use (a refusal is the error;
    /// the exception of a lookup is kept for the `Hit`, `SubjectNames::lookup_throws`).
    pub fn get(&mut self, s_idx: usize, program: &str) -> Result<&HitNames> {
        if self.names.get(s_idx).is_none_or(Option::is_none) {
            let made = self.subjects.hit_names(s_idx, program)?;
            *self
                .names
                .get_mut(s_idx)
                .context("HSP of an unknown subject")? = Some(made);
        }
        Ok(self.names[s_idx].as_ref().expect("made above"))
    }
}

#[cfg(test)]
#[path = "subjects_tests.rs"]
mod tests;
