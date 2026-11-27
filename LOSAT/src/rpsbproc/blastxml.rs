//! The BlastOutput XML documents of rpsblast and rpstblastn `-outfmt 5`, which rpsbproc reads one
//! after another from its input (`ObjStreamIn` of a `CBlastOutput`, in the format `eSerial_Xml`):
//! the module `NCBI-BlastOutput` read with `report::serial::xml_reader`, keeping what the output
//! prints (the program's version, the database, the matrix and E-value threshold of the parameters)
//! and, for each `Iteration`, the query that `ParseBlastOutput` makes of it (`offline`), built HSP
//! by HSP as the text is read: no document is held as text.
//!
//! NCBI reference (598d8ae6): c++/src/objects/blastxml/blastxml.asn:3-81
//! (The classes that datatool generates from this module read each member as its type:
//! VisibleString, INTEGER (`Int4`), BigInt (`Int8`), REAL; a member of a SEQUENCE type is that
//! class, a SEQUENCE OF a named type is a list of that class. Getters of OPTIONAL members
//! without a default throw when the member is not set; a list can always be got.)

use std::io::Read;

use super::data::CddData;
use super::offline::{Hsp, OfflineQuery};
use super::query::{DomQuery, ProcessError};
use crate::ncbi_str::c_atoi;
use crate::report::serial::asn_text_reader::{AsnSource, StreamSource};
use crate::report::serial::xml_reader::{member, optional, Member, Result, XmlError, XmlReader};

const BLAST_OUTPUT: [Member; 11] = [
    member("program"),
    member("version"),
    member("reference"),
    member("db"),
    member("query-ID"),
    member("query-def"),
    member("query-len"),
    optional("query-seq"),
    member("param"),
    member("iterations"),
    optional("mbstat"),
];

const PARAMETERS: [Member; 10] = [
    optional("matrix"),
    member("expect"),
    optional("include"),
    optional("sc-match"),
    optional("sc-mismatch"),
    member("gap-open"),
    member("gap-extend"),
    optional("filter"),
    optional("pattern"),
    optional("entrez-query"),
];

const STATISTICS: [Member; 7] = [
    member("db-num"),
    member("db-len"),
    member("hsp-len"),
    member("eff-space"),
    member("kappa"),
    member("lambda"),
    member("entropy"),
];

const ITERATION: [Member; 7] = [
    member("iter-num"),
    optional("query-ID"),
    optional("query-def"),
    optional("query-len"),
    optional("hits"),
    optional("stat"),
    optional("message"),
];

const HIT: [Member; 6] = [
    member("num"),
    member("id"),
    member("def"),
    member("accession"),
    member("len"),
    optional("hsps"),
];

const HSP: [Member; 20] = [
    member("num"),
    member("bit-score"),
    member("score"),
    member("evalue"),
    member("query-from"),
    member("query-to"),
    member("hit-from"),
    member("hit-to"),
    optional("pattern-from"),
    optional("pattern-to"),
    optional("query-frame"),
    optional("hit-frame"),
    optional("identity"),
    optional("positive"),
    optional("gaps"),
    optional("align-len"),
    optional("density"),
    member("qseq"),
    member("hseq"),
    optional("midline"),
];

/// What rpsbproc prints of one document: its session line and its queries (each the query of
/// an `Iteration`, or what `ParseBlastOutput` threw there).
#[derive(Debug, Default)]
pub(crate) struct Document {
    /// `BlastOutput_version`.
    pub(crate) version: Vec<u8>,
    /// `BlastOutput_db`.
    pub(crate) db: Vec<u8>,
    /// `Parameters_matrix`, OPTIONAL.
    pub(crate) matrix: Option<Vec<u8>>,
    /// `Parameters_expect`.
    pub(crate) expect: f64,
    pub(crate) queries: Vec<std::result::Result<DomQuery, ProcessError>>,
}

/// The documents of rpsbproc's input, read one after another from one stream.
pub(crate) struct XmlStream<R: Read> {
    reader: XmlReader<StreamSource<R>>,
    done: bool,
}

impl<R: Read> XmlStream<R> {
    /// The documents of `input`, whose first bytes `first` `Init` has read.
    pub(crate) fn new(input: R, first: &[u8]) -> Self {
        XmlStream {
            reader: XmlReader::new(StreamSource::new(input, first)),
            done: false,
        }
    }

    /// The next document (`ObjStreamIn<CBlastOutput>`) with its queries for the E-value cutoff
    /// `evcut` and the data `data`; `Ok(None)` when it does not read, where rpsbproc's loop ends
    /// (and after); `Err` with LOSAT's refusal.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:741-755
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/objutils.hpp:75-84
    // NCBI reference (598d8ae6): c++/src/util/bytesrc.cpp:194-199
    // (After a document the object stream gives back the bytes it did not read and clears the
    // stream's state, so the loop reads the next document from the byte after its last `>`.
    // When the input ends there, `good()` may be false: the loop ends as a read of nothing
    // fails.)
    pub(crate) fn next_document(
        &mut self,
        evcut: f64,
        data: &CddData,
    ) -> std::result::Result<Option<Document>, String> {
        if self.done {
            return Ok(None);
        }
        match read_document(&mut self.reader, evcut, data) {
            Ok(document) => Ok(Some(document)),
            // (`NotBool` is the object stream's only.)
            Err(XmlError::Format(_) | XmlError::NotBool) => {
                self.done = true;
                Ok(None)
            }
            Err(XmlError::Unsupported(message)) => {
                self.done = true;
                Err(message)
            }
        }
    }
}

/// A member of a type that rpsbproc does not read: read for the document's grammar only.
fn skip_string<S: AsnSource>(reader: &mut XmlReader<S>, scratch: &mut Vec<u8>) -> Result<()> {
    reader.string(scratch)
}

/// One document.
fn read_document<S: AsnSource>(
    reader: &mut XmlReader<S>,
    evcut: f64,
    data: &CddData,
) -> Result<Document> {
    reader.begin_document("BlastOutput")?;
    let mut document = Document::default();
    let mut scratch = Vec::new();
    let mut hsp = Hsp::default();
    reader.class("BlastOutput", &BLAST_OUTPUT, |reader, index| match index {
        1 => reader.string(&mut document.version),
        3 => reader.string(&mut document.db),
        6 => reader.int4().map(drop),
        8 => reader.class("Parameters", &PARAMETERS, |reader, index| match index {
            0 => {
                let mut matrix = Vec::new();
                reader.string(&mut matrix)?;
                document.matrix = Some(matrix);
                Ok(())
            }
            1 => {
                document.expect = reader.real()?;
                Ok(())
            }
            2 => reader.real().map(drop),
            3..=6 => reader.int4().map(drop),
            _ => skip_string(reader, &mut scratch),
        }),
        9 => reader.container(|reader| {
            let query = read_iteration(reader, evcut, &mut scratch, &mut hsp)?;
            document.queries.push(query.finish(data));
            Ok(())
        }),
        10 => read_statistics(reader),
        _ => skip_string(reader, &mut scratch),
    })?;
    Ok(document)
}

/// A `Statistics`, which rpsbproc does not read.
fn read_statistics<S: AsnSource>(reader: &mut XmlReader<S>) -> Result<()> {
    reader.class("Statistics", &STATISTICS, |reader, index| match index {
        0 | 2 => reader.int4().map(drop),
        1 => reader.int8().map(drop),
        _ => reader.real().map(drop),
    })
}

/// One `Iteration` and the query that `ParseBlastOutput` makes of it.
fn read_iteration<S: AsnSource>(
    reader: &mut XmlReader<S>,
    evcut: f64,
    scratch: &mut Vec<u8>,
    hsp: &mut Hsp,
) -> Result<OfflineQuery> {
    let mut query = OfflineQuery::default();
    reader.class("Iteration", &ITERATION, |reader, index| match index {
        0 => reader.int4().map(drop),
        // (The id and the title are read into their own strings, not copied from `scratch`: a
        // long one is not held twice.)
        1 => {
            let mut id = Vec::new();
            reader.string(&mut id)?;
            query.set_accession(id);
            Ok(())
        }
        2 => {
            let mut def = Vec::new();
            reader.string(&mut def)?;
            query.set_title(def);
            Ok(())
        }
        3 => {
            query.set_length(reader.int4()?);
            Ok(())
        }
        4 => reader.container(|reader| read_hit(reader, &mut query, evcut, scratch, hsp)),
        5 => read_statistics(reader),
        _ => skip_string(reader, scratch),
    })?;
    Ok(query)
}

/// One `Hit`, whose HSPs go to `query` as they are read.
fn read_hit<S: AsnSource>(
    reader: &mut XmlReader<S>,
    query: &mut OfflineQuery,
    evcut: f64,
    scratch: &mut Vec<u8>,
    hsp: &mut Hsp,
) -> Result<()> {
    let mut pssm_id = 0u32;
    reader.class("Hit", &HIT, |reader, index| match index {
        0 | 4 => reader.int4().map(drop),
        3 => {
            reader.string(scratch)?;
            // The hit's `GetAccession` text goes through `atoi` and is cast to `unsigned int`
            // (offl_cd_align_proc.cpp:155, cited at `OfflineQuery::parse_hsp`).
            pssm_id = c_atoi(scratch) as u32;
            Ok(())
        }
        5 => reader.container(|reader| {
            read_hsp(reader, scratch, hsp)?;
            query.hsp(pssm_id, hsp, evcut);
            Ok(())
        }),
        _ => skip_string(reader, scratch),
    })
}

/// One `Hsp`, into `hsp` (its OPTIONAL members `None` when the text lacks them).
fn read_hsp<S: AsnSource>(
    reader: &mut XmlReader<S>,
    scratch: &mut Vec<u8>,
    hsp: &mut Hsp,
) -> Result<()> {
    hsp.query_frame = None;
    hsp.identity = None;
    hsp.align_len = None;
    reader.class("Hsp", &HSP, |reader, index| {
        match index {
            1 => hsp.bit_score = reader.real()?,
            2 => {
                reader.real()?;
            }
            3 => hsp.evalue = reader.real()?,
            4 => hsp.query_from = reader.int4()?,
            5 => hsp.query_to = reader.int4()?,
            6 => hsp.hit_from = reader.int4()?,
            10 => hsp.query_frame = Some(reader.int4()?),
            12 => hsp.identity = Some(reader.int4()?),
            15 => hsp.align_len = Some(reader.int4()?),
            17 => reader.string(&mut hsp.qseq)?,
            18 => reader.string(&mut hsp.hseq)?,
            19 => reader.string(scratch)?,
            _ => {
                reader.int4()?;
            }
        }
        Ok(())
    })
}

#[cfg(test)]
#[path = "blastxml_tests.rs"]
mod tests;
