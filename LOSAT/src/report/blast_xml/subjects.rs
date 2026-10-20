//! How a program's final hit lists name and measure their subjects in the XML (`Hit_id`,
//! `Hit_def`, `Hit_accession`, `Hit_len`; DESIGN_OF4 2.3).

use std::sync::Arc;

use anyhow::{Context, Result};

use crate::blastinput::fasta_reader::FastaRecord;

use super::names::{self, Refusal};

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
    pub def: Vec<u8>,
    pub accession: Vec<u8>,
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
                let (id, def) =
                    names::xml_subject_names(record.local_id.as_bytes(), &record.title, *protein)
                        .map_err(|refusal| refusal.error(program, "a subject title"))?;
                Ok(SubjectNames {
                    id,
                    def,
                    accession: record.local_id.as_bytes().to_vec(),
                })
            }
            XmlSubjects::Database {
                ordinals, titles, ..
            } => {
                let ordinal = *ordinals
                    .get(s_idx)
                    .context("HSP of an unknown database sequence")?;
                let title = titles
                    .get(s_idx)
                    .context("HSP of an unknown database sequence")?;
                let def = names::xml_db_def(title)
                    .map_err(|refusal: Refusal| refusal.error(program, "a database title"))?;
                Ok(SubjectNames {
                    id: format!("gnl|BL_ORD_ID|{ordinal}").into_bytes(),
                    def,
                    accession: ordinal.to_string().into_bytes(),
                })
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
