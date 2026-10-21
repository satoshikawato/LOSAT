//! The framings of an XML2 / JSON document: where a query's report goes and what comes before and
//! after the reports (DESIGN_OF5 1.1, 1.2).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2332-2356
//! ```c++
//! void CBlastFormat::x_WriteXML2(CCmdLineBlastXML2ReportData & report_data)
//! {
//! 	if(m_FormatType == CFormattingArgs::eXml2_S) {
//! 			BlastXML2_FormatReport(&report_data, &m_Outfile);
//! 	}
//! 	else if (m_FormatType == CFormattingArgs::eJson_S) {
//! 			m_XMLFileCount++;
//! 			if(m_XMLFileCount > 1) {
//! 				m_Outfile << ",\n";
//! 			}
//! 			BlastJSON_FormatReport(&report_data, &m_Outfile);
//! 	}
//! ```
//! There is one document per query and every document is complete (the report header, the
//! parameters, one `Search`). Outfmt 16 writes them one after the other after the XML2 root's
//! start tag, outfmt 15 as the elements of one JSON array, outfmts 14 and 13 each to a file of
//! its own (`<base>_<n>.xml`, `<base>_<n>.json`) with a master file at `-out` (OF5c).
//!
//! `-subject` searches keep the last query's result back until the next query's arrives
//! (`m_IsBl2Seq`, the accumulated results of blast_format.cpp:2363-2404) and write it in the
//! epilog. NCBI's epilog is not reached after a failed query batch, so that held report is
//! never written: the document of a failed run lacks the last formatted query's report.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2364, :2383-2389
//! ```c++
//! 	if((m_IsIterative) || (m_IsBl2Seq)) {
//! 		if(m_AccumulatedResults.empty()) {
//! 			_ASSERT(m_AccumulatedQueries->size() == 0);
//! 			m_AccumulatedResults.push_back(res);
//! ...
//! 				CRef <CBlastSearchQuery> q = m_AccumulatedQueries->GetBlastSearchQuery(0);
//! 				if(m_IsBl2Seq) {
//! 				    CCmdLineBlastXML2ReportData report_data(q, m_AccumulatedResults, m_Options,
//! 				    		   	   	   	   	   	   	    m_Scope, m_SeqInfoSrc);
//! 				    x_WriteXML2(report_data);
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2206-2219
//! ```c++
//!     	if(!m_AccumulatedResults.empty()) {
//!     		CRef <CBlastSearchQuery> q = m_AccumulatedQueries->GetBlastSearchQuery(0);
//!     		if(m_IsBl2Seq) {
//!     			...
//! 		    	x_WriteXML2(report_data);
//! ```

use std::io::{self, Write};

use anyhow::{bail, Result};

use crate::report::serial::xml::XML2_ROOT_ATTRIBUTES;

/// Which of NCBI's four document forms a run writes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Xml2Kind {
    /// Outfmt 16 (`eXml2_S`): all reports on the output stream.
    Single16,
    /// Outfmt 15 (`eJson_S`): all reports in one JSON array on the output stream.
    Single15,
    /// Outfmt 14 (`eXml2`): one `<base>_<n>.xml` file per report and a master file.
    Files14,
    /// Outfmt 13 (`eJson`): one `<base>_<n>.json` file per report and a master file.
    Files13,
}

/// The files of a 13 / 14 run: `create(n)` opens `<base>_<n>.<ext>` for the n-th report.
pub trait SideFiles {
    fn create(&mut self, n: u32) -> io::Result<Box<dyn Write + '_>>;
}

/// The side files of a run that writes none (outfmt 15 and 16): outfmt 13 and 14 are written
/// from OF5c on.
pub struct NoSideFiles;

impl SideFiles for NoSideFiles {
    fn create(&mut self, _n: u32) -> io::Result<Box<dyn Write + '_>> {
        Err(io::Error::other(
            "outfmt 13 and 14 are not written by this build of LOSAT",
        ))
    }
}

/// The start of an outfmt 16 document (`BlastXML2_PrintHeader`, blastxml2_format.cpp:781-802):
/// the declaration and the `BlastXML2` root's start tag with its namespace attributes, which ends
/// the line before `>`.
pub fn header16() -> String {
    format!("<?xml version=\"1.0\" encoding=\"US-ASCII\"?>\n<BlastXML2{XML2_ROOT_ATTRIBUTES}>")
}

/// The start of an outfmt 15 document (blastxml2_format.cpp:834-837).
pub const HEAD15: &[u8] = b"{\n\"BlastOutput2\": [\n";

/// The end of an outfmt 16 document (`PrintEpilog`, blast_format.cpp:2419-2423).
pub const TAIL16: &[u8] = b"</BlastXML2>\n";

/// The end of an outfmt 15 document (blast_format.cpp:2440-2445).
pub const TAIL15: &[u8] = b"]\n}\n";

/// What `PrintProlog` writes for a document of `kind` (nothing for the side-file forms).
pub fn prolog_bytes(kind: Xml2Kind) -> Vec<u8> {
    match kind {
        Xml2Kind::Single16 => header16().into_bytes(),
        Xml2Kind::Single15 => HEAD15.to_vec(),
        Xml2Kind::Files14 | Xml2Kind::Files13 => Vec::new(),
    }
}

/// The document of one run: the reports of the formatted queries and the state that the
/// framing needs (the count of reports, the held report of a `-subject` run).
pub struct Xml2Document<F: SideFiles = NoSideFiles> {
    kind: Xml2Kind,
    /// `m_XMLFileCount`: the reports written so far.
    written: u32,
    /// Whether the document's start has been written (NCBI writes it before the search; LOSAT
    /// at the first write, which gives the same bytes on the output stream).
    started: bool,
    /// A `-subject` run (`m_IsBl2Seq`): a report is written when the next one arrives.
    delayed: bool,
    pending: Vec<u8>,
    has_pending: bool,
    files: F,
}

impl Xml2Document<NoSideFiles> {
    /// A document of kind `kind` on the output stream; `delayed` for a `-subject` search.
    pub fn new(kind: Xml2Kind, delayed: bool) -> Self {
        Self::with_files(kind, delayed, NoSideFiles)
    }
}

impl<F: SideFiles> Xml2Document<F> {
    pub fn with_files(kind: Xml2Kind, delayed: bool, files: F) -> Self {
        Self {
            kind,
            written: 0,
            started: false,
            delayed,
            pending: Vec::new(),
            has_pending: false,
            files,
        }
    }

    pub fn kind(&self) -> Xml2Kind {
        self.kind
    }

    /// The document's start, once (`PrintProlog`, blast_format.cpp:360-367).
    fn prolog(&mut self, out: &mut Vec<u8>) {
        if std::mem::replace(&mut self.started, true) {
            return;
        }
        out.extend_from_slice(&prolog_bytes(self.kind));
    }

    /// Hands one query's complete report (`x_PrintXML2Report`): appends to `out` what the
    /// framing writes now. `report` is left empty. For a `-subject` run the report is held
    /// and the one held before it is written.
    pub fn report(&mut self, report: &mut Vec<u8>, out: &mut Vec<u8>) -> Result<()> {
        self.prolog(out);
        if self.delayed {
            // The report just made is held; the one held before it is written now.
            let had_pending = std::mem::replace(&mut self.has_pending, true);
            std::mem::swap(report, &mut self.pending);
            if !had_pending {
                report.clear();
                return Ok(());
            }
        }
        let result = self.write(report, out);
        report.clear();
        result
    }

    /// `x_WriteXML2`: one report into the stream or its own file.
    fn write(&mut self, report: &[u8], out: &mut Vec<u8>) -> Result<()> {
        self.written += 1;
        match self.kind {
            Xml2Kind::Single16 => out.extend_from_slice(report),
            Xml2Kind::Single15 => {
                if self.written > 1 {
                    out.extend_from_slice(b",\n");
                }
                out.extend_from_slice(report);
            }
            Xml2Kind::Files14 | Xml2Kind::Files13 => {
                let mut file = self.files.create(self.written)?;
                file.write_all(report)?;
                file.flush()?;
            }
        }
        Ok(())
    }

    /// `PrintEpilog` (blast_format.cpp:2203-2230): the held report of a `-subject` run, then
    /// the document's end. A run that stops at a failed batch does not call it.
    pub fn epilog(&mut self, out: &mut Vec<u8>) -> Result<()> {
        self.prolog(out);
        if std::mem::replace(&mut self.has_pending, false) {
            let pending = std::mem::take(&mut self.pending);
            self.write(&pending, out)?;
        }
        match self.kind {
            Xml2Kind::Single16 => out.extend_from_slice(TAIL16),
            Xml2Kind::Single15 => out.extend_from_slice(TAIL15),
            Xml2Kind::Files14 | Xml2Kind::Files13 => {
                bail!("the master file of outfmt 13 and 14 is not written by this build of LOSAT")
            }
        }
        Ok(())
    }
}
