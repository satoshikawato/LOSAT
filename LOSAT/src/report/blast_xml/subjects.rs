//! How a program's final hit lists name and measure their subjects in the XML (`Hit_id`,
//! `Hit_def`, `Hit_accession`, `Hit_len`; DESIGN_OF4 2.3).

use std::sync::Arc;

use anyhow::{Context, Result};

use crate::blastinput::fasta_reader::FastaRecord;

use super::names::{self, Refusal};
use super::records::HitTaxonomy;

/// The subjects of a search, indexed by the subject number of its HSPs.
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
    /// A database: the sequences with hits by rank (`ordinals[rank]` the ordinal id, the
    /// stored title `titles.get(rank)`, `length(rank)` the sequence length).
    Database {
        ordinals: &'a [u32],
        titles: XmlTitles<'a>,
        /// The Blast-def-line taxid of each sequence by rank, 0 for none; only XML2 shows
        /// it. Empty where the program has not read them (outfmt 5 alone): every taxid is
        /// then 0.
        taxids: &'a [i32],
        /// What the taxonomy database said about each sequence's taxid by rank
        /// (`HitTaxonomy`); empty where nothing was asked: every one is `Unnamed`.
        taxonomy: &'a [HitTaxonomy],
        length: &'a dyn Fn(usize) -> usize,
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
}

impl SubjectNames {
    /// The names of a `-subject` record `id` with the FASTA title `title` (`Subject_N` and the
    /// generated definition line). `protein`: the record is a protein.
    pub fn of_record(id: &[u8], title: &[u8], protein: bool) -> Result<Self, Refusal> {
        let (xml_id, title) = names::xml_subject_names(id, title, protein)?;
        Ok(SubjectNames {
            id: xml_id,
            title,
            accession: id.to_vec(),
            taxid: 0,
            taxonomy: HitTaxonomy::Unnamed,
        })
    }

    /// The names of the database sequence `ordinal` with the stored title `title`
    /// (`gnl|BL_ORD_ID|n`, the title as stored).
    pub fn of_database_sequence(
        ordinal: impl std::fmt::Display,
        title: &[u8],
        taxid: i32,
        taxonomy: HitTaxonomy,
    ) -> Result<Self, Refusal> {
        Ok(SubjectNames {
            id: format!("gnl|BL_ORD_ID|{ordinal}").into_bytes(),
            title: names::xml_db_title(title)?,
            accession: ordinal.to_string().into_bytes(),
            taxid,
            taxonomy,
        })
    }
}

impl XmlSubjects<'_> {
    /// The names of subject `s_idx`; `program` (`BLASTP`) names the program in the refusal of
    /// a title that NCBI cannot write.
    pub fn names(&self, s_idx: usize, program: &str) -> Result<SubjectNames> {
        match self {
            XmlSubjects::Records {
                records, protein, ..
            } => {
                let record = records.get(s_idx).context("HSP of an unknown subject")?;
                SubjectNames::of_record(record.local_id.as_bytes(), &record.title, *protein)
                    .map_err(|refusal| refusal.error(program, "a subject title"))
            }
            XmlSubjects::Database {
                ordinals,
                titles,
                taxids,
                taxonomy,
                ..
            } => {
                let ordinal = *ordinals
                    .get(s_idx)
                    .context("HSP of an unknown database sequence")?;
                let title = titles
                    .get(s_idx)
                    .context("HSP of an unknown database sequence")?;
                let taxid = taxids.get(s_idx).copied().unwrap_or(0);
                let taxonomy = taxonomy.get(s_idx).cloned().unwrap_or_default();
                SubjectNames::of_database_sequence(ordinal, title, taxid, taxonomy)
                    .map_err(|refusal: Refusal| refusal.error(program, "a database title"))
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
}

/// The names of the subjects of a report, made once per subject (the definition line of a
/// subject is generated when it first has a hit, not again for each hit or each query that hits
/// it).
pub struct SubjectNameCache<'s, 'a> {
    subjects: &'s XmlSubjects<'a>,
    names: Vec<Option<SubjectNames>>,
}

impl<'s, 'a> SubjectNameCache<'s, 'a> {
    pub fn new(subjects: &'s XmlSubjects<'a>) -> Self {
        let mut names = Vec::new();
        names.resize_with(subjects.count(), || None);
        Self { subjects, names }
    }

    /// `XmlSubjects::names` of subject `s_idx`, made on first use.
    pub fn get(&mut self, s_idx: usize, program: &str) -> Result<&SubjectNames> {
        if self.names.get(s_idx).is_none_or(Option::is_none) {
            let made = self.subjects.names(s_idx, program)?;
            *self
                .names
                .get_mut(s_idx)
                .context("HSP of an unknown subject")? = Some(made);
        }
        Ok(self.names[s_idx].as_ref().expect("made above"))
    }
}
