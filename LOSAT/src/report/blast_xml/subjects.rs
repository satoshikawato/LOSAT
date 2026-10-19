//! How a program's final hit lists name and measure their subjects in the XML (`Hit_id`,
//! `Hit_def`, `Hit_accession`, `Hit_len`; DESIGN_OF4 2.3).

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
    /// stored title `records[rank].title`, `length(rank)` the sequence length).
    Database {
        ordinals: &'a [u32],
        records: &'a [FastaRecord],
        length: &'a dyn Fn(usize) -> usize,
    },
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
                ordinals, records, ..
            } => {
                let ordinal = *ordinals
                    .get(s_idx)
                    .context("HSP of an unknown database sequence")?;
                let record = records
                    .get(s_idx)
                    .context("HSP of an unknown database sequence")?;
                let def = names::xml_db_def(&record.title)
                    .map_err(|refusal: Refusal| refusal.error(program, "a database title"))?;
                Ok(SubjectNames {
                    id: format!("gnl|BL_ORD_ID|{ordinal}").into_bytes(),
                    def,
                    accession: ordinal.to_string().into_bytes(),
                })
            }
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
