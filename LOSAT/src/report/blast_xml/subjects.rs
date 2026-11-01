//! How a program's final hit lists name and measure their subjects in the XML (`Hit_id`,
//! `Hit_def`, `Hit_accession`, `Hit_len`; DESIGN_OF4 2.3).

use std::sync::Arc;

use anyhow::{Context, Result};

use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::db_subject::{HitDescription, ParsedSubject};

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
        /// The subjects with parsed Seq-ids by rank (empty without them), whose names come from
        /// their definition lines (`SubjectNames::of_parsed`).
        parsed: &'a [Option<Arc<ParsedSubject>>],
        /// What the taxonomy database said about the taxid of each definition line of a subject
        /// with parsed Seq-ids, by rank (empty where nothing was asked;
        /// `algorithm::common::database::read_reported_taxonomy`).
        parsed_taxonomy: &'a [Vec<HitTaxonomy>],
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
            descriptions: Vec::new(),
            volume_taxids: Vec::new(),
        })
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
        let title = subject
            .bioseq_title(protein, false)
            .map_err(|_| Refusal::UnknownEncoding)?;
        if title.reads_past_end {
            return Err(Refusal::ReadsPastEnd);
        }
        names::check_text(&id)?;
        names::check_text(&title.text)?;
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:393-406
        // ```c++
        //         	CRef<CSeq_id> best_id = FindBestChoice(ids, CSeq_id::Score);
        //         	if(!best_id->IsLocal()) {
        //         		CSeq_id_Handle id_handle = CSeq_id_Handle::GetHandle(*best_id);
        //         		string accession = CAlignFormatUtil::GetLabel(id_handle.GetSeqId());
        //         		if(accession != kEmptyStr)
        //         			hit_exp->SetAccession(accession);
        //         	}
        //
        //    			hit_exp->SetTitle(sequence::CDeflineGenerator().GenerateDefline(subj_handle));
        // ```
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
            descriptions: Vec::new(),
            volume_taxids: Vec::new(),
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml2_format.cpp:380-388, 422-428
    // ```c++
    //         		if(defline.IsSetTaxid() && defline.GetTaxid() != ZERO_TAX_ID) {
    //         			TTaxId tax_id = defline.GetTaxid();
    //        				hit_exp->SetTaxid(tax_id);
    //         			if(hasTaxDB) {
    //         				 SSeqDBTaxInfo taxinfo;
    //         				 CSeqDB::GetTaxInfo(tax_id, taxinfo);
    //         				 hit_exp->SetSciname(taxinfo.scientific_name);
    //         			}
    //         		}
    //         		descr_list.push_back(hit_exp);
    //     ...
    //     } catch (const CException&) {
    //     	CRef<blastxml2::CHitDescr> hit_exp(new blastxml2::CHitDescr);
    //     	hit_exp->SetId(kSeqId.AsFastaString());
    //         hit->SetDescription().push_back(hit_exp);
    // ```
    /// The names of a database sequence with parsed Seq-ids (`ParsedSubject::xml_hit` and
    /// `hit_descriptions`): `Hit_id`, `Hit_def` and `Hit_accession` of outfmt 5, and the
    /// `HitDescr` of each definition line of XML2 and JSON with what the taxonomy database said
    /// about its taxid (`taxonomy[i]` for definition line `i`; empty where nothing was asked).
    /// `GetTaxInfo` throws for a taxid that the taxonomy database does not know: the
    /// descriptions before that definition line stay, and the Seq-align id alone follows them.
    pub fn of_parsed(subject: &ParsedSubject, taxonomy: &[HitTaxonomy]) -> Result<Self, Refusal> {
        let (id, title, accession) = subject.xml_hit();
        names::check_text(&id)?;
        names::check_text(&title)?;
        let mut descriptions = subject.hit_descriptions();
        for (at, answer) in taxonomy.iter().enumerate().take(descriptions.len()) {
            if matches!(answer, HitTaxonomy::NotFound(_)) {
                descriptions.truncate(at);
                descriptions.push(HitDescription {
                    id: subject.alignment_fasta(),
                    accession: Vec::new(),
                    title: None,
                    taxid: 0,
                    taxonomy: answer.clone(),
                });
                break;
            }
            descriptions[at].taxonomy = answer.clone();
        }
        Ok(SubjectNames {
            id,
            title,
            accession,
            taxid: 0,
            taxonomy: HitTaxonomy::Unnamed,
            descriptions,
            volume_taxids: subject.volume_taxids(),
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
                titles,
                taxids,
                taxonomy,
                parsed,
                parsed_taxonomy,
                ..
            } => {
                if let Some(subject) = parsed.get(s_idx).and_then(Option::as_deref) {
                    let answers = parsed_taxonomy.get(s_idx).map_or(&[][..], Vec::as_slice);
                    return SubjectNames::of_parsed(subject, answers)
                        .map_err(|refusal: Refusal| refusal.error(program, "a database title"));
                }
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
