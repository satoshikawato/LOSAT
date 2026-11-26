//! The retrieval of sequences: NCBI's `x_ProcessSearchRequest` (the configuration of
//! `-outfmt`, `-strand`, `-ctrl_a`, `-mask_sequence_with`, `-range`), its dispatch
//! (`-entry all`, `-entry`, `-entry_batch` with and without `-get_dups`), and the lookup of
//! ids: in a version 5 database the LMDB table of accessions and the ISAM of GIs, in a version
//! 4 database `SeqDB_SimplifyAccession` and the ISAMs of its volumes (strings, GIs, PIGs), with
//! `gnl|BL_ORD_ID|n` read without a file, and the answers of a database built without
//! `-parse_seqids`.

use std::cell::{Cell, RefCell};
use std::collections::BTreeSet;
use std::ffi::OsStr;
use std::io::Write;

use crate::algorithm::common::database::TaxdbRun;
use crate::ncbi_app::args::Args;
use crate::ncbi_app::files::{input_lines, read_input_file, Cin, LineEnds, OutputFile};
use crate::ncbi_str::string_to_non_negative_int;
use crate::report::ostream::{with_stream_inserter, CutBuffer, ReportWrite, StreamInserter};
use crate::seqdb::isam::{seqdb_simplify_accession, SimplifyError};
use crate::seqdb::lookup::{IdLookup, ID_OVERFLOWS};
use crate::seqdb::tax_files::NO_ACCESSION_INFO;
use crate::seqdb::taxdb::{SourceNames, TaxCache};
use crate::seqdb::{MaskBlock, MaskRanges, SeqDb, SeqDbError, SeqLocation, SeqType};
use crate::seqid::accguide::parse_flags::{ANY_LOCAL, DEFAULT, PARTIAL_OK, RAW_TEXT};
use crate::seqid::label::{content, fasta};
use crate::seqid::parse::seq_id_from_text;
use crate::seqid::{DeflineSet, ParsedSeqId};

use super::args;
use super::bioseq::{write_asn_text, AsnTextDocument, BioseqFormatter, DbSeqData};
use super::fasta::FastaFormatter;
use super::format::FormatSpec;
use super::range::{parse_sequence_range_open_end, SeqRange, NO_RANGE};
use super::Diag;
use crate::sequence::stdaa_letters::{protein_fasta_letters, protein_string_letters};

/// What ends the retrieval before its last sequence.
#[derive(Debug)]
pub(super) enum Stop {
    /// An exception NCBI's `x_ProcessSearchRequest` catches: its message as a diagnostic
    /// (`Error: [blastdbcmd] <text>`), exit status 1.
    Ncbi(String),
    /// The same, for a message that is not UTF-8 (the unknown specifier read from the NUL
    /// that ends `-outfmt`).
    NcbiBytes(Vec<u8>),
    /// LOSAT's refusal: a diagnostic line, exit status 1.
    Refused(String),
    /// Not an end: NCBI's warning of a UTF8String that its read of a record took for a
    /// VisibleString (the member it names, `seqid::ber::utf8_string_path`), posted at its first
    /// such read (`Source::deflines`), before the record's text. `Run::write` posts it and
    /// writes the record again; the warning is posted once.
    Utf8Post(String),
    /// A write found the pipe closed where SIGPIPE ends NCBI's process (`OutputFile::closed`):
    /// the retrieval ends at once, exit status 0, no message. Any other failed write is
    /// dropped by the `OutputFile` and the retrieval goes on.
    ClosedPipe,
}

impl Stop {
    /// The error of a write to the retrieval's `OutputFile`, which fails only on a closed pipe.
    pub fn write(_: std::io::Error) -> Stop {
        Stop::ClosedPipe
    }

    pub fn refused(reason: String) -> Stop {
        Stop::Refused(super::refusal(&reason))
    }
}

impl From<SeqDbError> for Stop {
    /// A reader error inside the request: NCBI's own (`CSeqDBException`, the `CFileException`
    /// of a lookup file, or another exception such as the `CSerialException` of a header record
    /// that NCBI's read fails on) is caught like any exception there (exit status 1, not the 2 of
    /// an error at opening), with its message (`GetMsg()`, `SeqDbError::ncbi_request_message`);
    /// LOSAT's is its error of a file that cannot be read or a refusal
    /// (`common::database::losat_reader_error`).
    fn from(error: SeqDbError) -> Stop {
        match error.ncbi_request_message() {
            Some(text) => Stop::NcbiBytes(text),
            None => Stop::Refused(crate::algorithm::common::database::losat_reader_error(
                &error,
                super::PROGRAM,
            )),
        }
    }
}

/// `CBlastDB_FormatterConfig`. `m_FmtAlgoId` is in the `%` formatter, one id per `%m`
/// (line BDC DECISIONS 12).
#[derive(Clone, Debug)]
pub(super) struct Config {
    /// `m_SeqRange`.
    pub range: SeqRange,
    /// `m_Strand` is `eNa_strand_minus` (never for a protein database).
    pub strand_minus: bool,
    /// `m_UseCtrlA`.
    pub ctrl_a: bool,
    /// `m_FiltAlgoId`: the filtering algorithm whose masks are in lower case, -1 for none.
    pub filt_algo: i32,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            range: NO_RANGE,
            strand_minus: false,
            ctrl_a: false,
            filt_algo: -1,
        }
    }
}

/// The database as the formatters read it.
pub(super) struct Source<'a> {
    pub db: &'a SeqDb,
    /// Whether the database is one of proteins.
    pub protein: bool,
    /// `m_DbIsProtein`: `protein`, but for a taxid list whether `-dbtype` is `prot`
    /// (`x_InitBlastDB_TaxIdList`). It decides whether `-strand minus` is taken.
    dbtype_protein: bool,
    /// The taxonomy database of the search path (NCBI's one per process) and the lookups of the
    /// volumes' caches, each `GetBioseq` of the FASTA formatter a load.
    taxdb: TaxdbRun,
    /// What `tax_info` looked up, as NCBI's volume keeps it (`seqdb::taxdb::TaxCache`).
    source_names: TaxCache,
    /// The buffers of `masks`, kept from one sequence to the next (`SeqDb::mask_ranges_in`).
    mask_block: RefCell<MaskBlock>,
    /// The definition lines of the last read that `deflines` gave back (`HeaderSet`), whose
    /// buffers the next read writes over (`SeqDb::header_set_in`).
    kept_set: Cell<DeflineSet>,
}

/// The definition lines that [`Source::deflines`] read, given back to the source when dropped:
/// the next read writes over their buffers, so a run allocates the lines, titles and Seq-id lists
/// of its sequences once, not for each sequence.
pub(super) struct HeaderSet<'s> {
    set: DeflineSet,
    home: &'s Cell<DeflineSet>,
}

impl<'s> HeaderSet<'s> {
    /// The set that `home` kept (a new one where a read holds it).
    fn take(home: &'s Cell<DeflineSet>) -> Self {
        HeaderSet {
            set: home.take(),
            home,
        }
    }
}

impl std::ops::Deref for HeaderSet<'_> {
    type Target = DeflineSet;

    fn deref(&self) -> &DeflineSet {
        &self.set
    }
}

impl Drop for HeaderSet<'_> {
    fn drop(&mut self) {
        self.home.set(std::mem::take(&mut self.set));
    }
}

impl<'a> Source<'a> {
    fn new(db: &'a SeqDb, dbtype_protein: bool) -> Self {
        Source {
            db,
            protein: db.seqtype() == SeqType::Protein,
            dbtype_protein,
            taxdb: TaxdbRun::sequential(),
            source_names: TaxCache::default(),
            mask_block: RefCell::default(),
            kept_set: Cell::default(),
        }
    }

    /// NCBI's warning of the UTF8String that its read of the header at `at` takes for a
    /// VisibleString (`Stop::Utf8Post` of the member, `seqid::ber::utf8_string_path`), where
    /// the run has not posted it (`SeqDb::take_utf8_post`: once per process).
    #[cold]
    #[inline(never)]
    fn utf8_post(&self, at: SeqLocation) -> Option<Stop> {
        self.db.take_utf8_post_at(at).map(Stop::Utf8Post)
    }

    /// The stop of a read of the header at `at` that failed with `error`: NCBI's warning first
    /// where its read took a UTF8String before the failure (`utf8_post`), else the error.
    #[cold]
    #[inline(never)]
    fn failed_read(&self, at: SeqLocation, error: crate::seqdb::SeqDbError) -> Stop {
        self.utf8_post(at).unwrap_or_else(|| error.into())
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:411-425
    // Every read of the formatters finds the sequence by NCBI's OID (`FindVol`), not by the
    // OID list: the OID of a PIG that an alias file's range leaves out is read as any other
    // (`-ipg`, probe p4c_c of line RG D1a-2).
    /// The place of the OID `oid` (`SeqDb::locate_oid`; NCBI's `OID not found` past the
    /// volumes, which the lookups never give).
    fn at(&self, oid: u32) -> Result<SeqLocation, Stop> {
        self.db
            .locate_oid(oid)
            .map_err(|_| Stop::Ncbi(OID_NOT_FOUND.to_string()))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:286-298
    /// `GetSeqLength` of `oid`: signed 32-bit `int` for proteins (wrapping `>= 2^31` to negative
    /// as in NCBI's `GetSeqLengthProt`), widened unsigned 32-bit for nucleotides.
    pub fn length(&self, oid: u32) -> Result<i64, Stop> {
        let at = self.at(oid)?;
        if self.protein {
            Ok(i64::from(self.db.protein_length_i32_at(at)?))
        } else {
            let length = self.db.nucleotide_length_at(at)?;
            u32::try_from(length)
                .map(i64::from)
                .map_err(|_| Stop::refused(format!("sequence {oid} is longer than 32 bits")))
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1380-1386
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1784-1804
    /// Validates the protein `.psq` open and range for a sequence whose signed 32-bit length is
    /// below 1, matching `GetBioseq`'s `x_GetSequence` call before it returns an empty Bioseq.
    pub(super) fn check_non_positive_protein_range(&self, oid: u32) -> Result<(), Stop> {
        if self.protein {
            self.db.check_protein_range_at(self.at(oid)?)?;
        }
        Ok(())
    }

    /// `GetHdr` of `oid`: the definition lines with their Seq-ids; `None`, NCBI's null set, for a
    /// stored header of no bytes ([`Source::null_set`]). The first read of a record with a
    /// UTF8String where a VisibleString is expected gives NCBI's warning first (`Stop::Utf8Post`),
    /// also where the read then fails.
    pub fn deflines(&self, oid: u32) -> Result<Option<HeaderSet<'_>>, Stop> {
        let at = self.at(oid)?;
        let mut set = HeaderSet::take(&self.kept_set);
        let read = self.db.header_set_in(at, &mut set.set);
        if !read.map_err(|error| self.failed_read(at, error))? {
            return self.null_set(at);
        }
        self.check_read(oid, at, &set)?;
        Ok(Some(set))
    }

    /// `GetHdr` of `oid` and the header that `GetBioseq` keeps in the Bioseq's user object
    /// (`SeqDb::defline_set_with_header_in`); `None` and NCBI's warning as for
    /// [`Source::deflines`].
    pub fn deflines_with_header(&self, oid: u32) -> Result<Option<(HeaderSet<'_>, Vec<u8>)>, Stop> {
        let at = self.at(oid)?;
        let mut set = HeaderSet::take(&self.kept_set);
        let read = self.db.defline_set_with_header_in(at, &mut set.set);
        let Some(header) = read.map_err(|error| self.failed_read(at, error))? else {
            return self.null_set(at);
        };
        self.check_read(oid, at, &set)?;
        Ok(Some((set, header)))
    }

    /// The checks of a set read from the header at `at`: NCBI's warning of its first UTF8String
    /// read (`utf8_post`), then the refusal of bytes NCBI reads as '#' (`check_read_bytes`).
    fn check_read(&self, oid: u32, at: SeqLocation, set: &DeflineSet) -> Result<(), Stop> {
        if set.utf8 {
            if let Some(post) = self.utf8_post(at) {
                return Err(post);
            }
        }
        check_read_bytes(oid, set)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1997-2003
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2202-2212
    // A database opened with a taxid list (`-taxids`, `-taxidlist` with `-target_only`), or a
    // volume under an alias file's `OID_MASK_TYPE`, filters the set it read, also the null one
    // (probes w1-w4 and O1-O4 of line RG-Dbc); otherwise the readers see no set: the `%` formatter
    // writes nothing, `GetBioseq` gives no Bioseq, and its read of a target's line throws (probes
    // z1-z9, za-zd).
    /// The read of the sequence at `at` whose stored header has no bytes: `None`, or NCBI's null
    /// pointer exception where the read filters its definition lines
    /// ([`SeqDb::filters_header_sets`]).
    fn null_set<T>(&self, at: SeqLocation) -> Result<Option<T>, Stop> {
        if self.db.filters_header_sets(at) {
            return Err(null_pointer());
        }
        Ok(None)
    }

    /// The `Seq-data` of `GetBioseq` of `oid`: the stored NCBIstdaa of a protein, the stored 2-bit
    /// bytes or NCBI4na of a nucleotide (`SequenceBlock::nucleotide_seq_data`).
    pub fn seq_data(&self, oid: u32) -> Result<DbSeqData, Stop> {
        let at = self.at(oid)?;
        if self.protein {
            if self.db.protein_length_i32_at(at)? < 1 {
                self.db.check_protein_range_at(at)?;
                return Ok(DbSeqData::Ncbistdaa(Vec::new()));
            }
            let block = self.db.read_sequence_at(at)?;
            return Ok(DbSeqData::Ncbistdaa(block.protein(0).to_vec()));
        }
        let block = self.db.read_sequence_at(at)?;
        Ok(match block.nucleotide_seq_data(0)? {
            (true, bytes) => DbSeqData::Ncbi2na(bytes),
            (false, bytes) => DbSeqData::Ncbi4na(bytes),
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:1607-1636
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:1153-1158
    /// The letters of `GetSequenceAsString` of `oid`, with `range` when one is given (`-outfmt %s`,
    /// `%h`): iupacaa with the byte 255 for a gap or stop, or the IUPAC letters of a nucleotide.
    pub fn string_letters(&self, oid: u32, range: Option<(u32, u32)>) -> Result<Vec<u8>, Stop> {
        let no_seq = || Stop::Ncbi("Error: could not get sequence or range.".to_string());
        if self.protein {
            let at = self.at(oid)?;
            let base_length = self.db.protein_length_i32_at(at)?;
            if base_length < 1 {
                self.db.check_protein_range_at(at)?;
                if let Some((from, to)) = range {
                    let begin = from as i32;
                    let end = to.wrapping_add(1) as i32;
                    if end > base_length {
                        return Err(Stop::Ncbi(
                            "Error: region beyond sequence range.".to_string(),
                        ));
                    }
                    // Open-ended `-range N-` on a `2^31 + k` protein clips `to` to `length - 1`,
                    // so `end == base_length == -2^31 + k`; 32-bit `end - begin` wraps positive
                    // for `N >= k + 2` and reads only the requested tail slice from `.psq`.
                    let slice_len = end.wrapping_sub(begin);
                    if slice_len >= 1 {
                        let stdaa = self.db.read_protein_slice_at(
                            at,
                            u64::from(from),
                            slice_len as usize,
                        )?;
                        return Ok(protein_string_letters(&stdaa));
                    }
                }
                return Err(no_seq());
            }
            let block = self.db.read_sequence_at(at)?;
            Ok(protein_string_letters(slice(block.protein(0), range)))
        } else {
            self.nucleotide_letters(oid, range)
        }
    }

    /// The letters of the FASTA path (`CSeqVector` in `Ncbieaa` or `Iupacna`).
    pub fn fasta_letters(&self, oid: u32, range: Option<(u32, u32)>) -> Result<Vec<u8>, Stop> {
        if self.protein {
            let at = self.at(oid)?;
            if self.db.protein_length_i32_at(at)? < 1 {
                self.db.check_protein_range_at(at)?;
                return Ok(Vec::new());
            }
            let block = self.db.read_sequence_at(at)?;
            Ok(protein_fasta_letters(slice(block.protein(0), range)))
        } else {
            self.nucleotide_letters(oid, range)
        }
    }

    fn nucleotide_letters(&self, oid: u32, range: Option<(u32, u32)>) -> Result<Vec<u8>, Stop> {
        let at = self.at(oid)?;
        Ok(match range {
            Some((from, to)) => {
                self.db
                    .nucleotide_iupac_range_at(at, from as usize, to as usize + 1)?
            }
            None => self.db.nucleotide_iupac_at(at)?,
        })
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2234-2242
    // `RealToVol` throws for an algorithm the volume does not know, after the cache's key was
    // written, and `GetMaskData` catches it: that sequence has no masks, and the next one of
    // the volume reads the masks of the id cached before, or of an uninitialised one.
    /// `GetMaskData` with `oid` and `algo`: the stored ranges of the sequence, or the refusal of a
    /// volume that has masks but not this algorithm (line NP's M/DECISIONS 9).
    pub fn masks(&self, oid: u32, algo: i32) -> Result<Vec<(u32, u32)>, Stop> {
        let at = self.at(oid)?;
        match self
            .db
            .mask_ranges_in(at, algo, &mut self.mask_block.borrow_mut())?
        {
            MaskRanges::Ranges(ranges) => Ok(ranges),
            MaskRanges::AlgorithmNotInVolume { volume } => Err(Stop::refused(format!(
                "filtering algorithm {algo} in {}, a volume of an alias or a list of names with \
                 masks of other algorithms (NCBI reads the masks of another algorithm there)",
                volume.display()
            ))),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1614-1621
    /// `CSeqDB::GetTaxInfo` of `taxid` (and `GetBioseq`'s lookup of a source,
    /// [`Source::tax_names`]): the names of the taxonomy database on the search path, none where
    /// NCBI's lookup throws (no database, no such taxid), and the `ERR_POST` texts of the lookup
    /// (of opening the database, when it is the run's first lookup, and the offset error), which
    /// NCBI writes where it looks the taxid up. The lookups are kept as NCBI's volume keeps them
    /// (`source_names`), which keeps no lookup that writes a text: a kept answer is that of a
    /// lookup that writes none. Err: LOSAT's refusal of a lookup that reads a record of an empty
    /// `taxdb.btd` (NCBI faults) or whose names NCBI reads past its mapping of `taxdb.btd` without
    /// a defined result (`TaxLookup::into_names`); a super kingdom past the mapping is refused
    /// where the format writes it (`format::DeflineFields`).
    pub fn tax_info(&self, taxid: i32) -> Result<(std::sync::Arc<SourceNames>, Vec<String>), Stop> {
        let mut texts = Vec::new();
        let names = self.source_names.names(taxid, &mut texts, |messages| {
            // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:873-886
            // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/blastdb_dataextract.cpp:902-912
            // The formatter's names catch any exception of the lookup, so a file of the taxonomy
            // database that does not open gives `N/A` (line WQA D18b probe 12, `%S %L %B %K`).
            match self.taxdb.lookup(self.db, taxid, messages) {
                Err(SeqDbError::TaxdbThrows { .. }) => Ok(None),
                lookup => lookup?.into_names(),
            }
        })?;
        Ok((names, texts))
    }

    /// The names of `taxid` that the Bioseq's source holds (a title made from it,
    /// [`Source::tax_info`]): the load of the sequence ([`Source::post_bioseq_texts`], before)
    /// wrote the lookup's texts.
    pub fn tax_names(&self, taxid: i32) -> Result<std::sync::Arc<SourceNames>, Stop> {
        Ok(self.tax_info(taxid)?.0)
    }

    /// The lookups of `CSeqDB::GetBioseq` of the OID `oid`, whose definition lines are `set`, for
    /// the target id `target` (`DeflineSet::taxonomy_lookups`), for a Bioseq that NCBI makes (one
    /// with definition lines and letters): their texts posted to `out` before any byte of the
    /// sequence ([`TaxdbRun::bioseq_texts`]). Err: LOSAT's refusal of a lookup that NCBI faults on
    /// for certain (a record of an empty `taxdb.btd`). Inlined into the formatter, which calls it
    /// for every sequence.
    #[inline]
    pub fn post_bioseq_texts(
        &self,
        oid: u32,
        set: &DeflineSet,
        target: Option<&ParsedSeqId>,
        out: &mut Output<'_, '_>,
    ) -> Result<(), Stop> {
        let texts = self.taxdb.bioseq_texts(self.db, oid, set, target)?;
        for text in texts {
            out.post(&text)?;
        }
        Ok(())
    }
}

/// The refusal of a sequence whose definition lines hold bytes that NCBI reads as '#'.
fn check_read_bytes(oid: u32, set: &DeflineSet) -> Result<(), Stop> {
    if set
        .deflines
        .iter()
        .any(|defline| defline.replaced_title_bytes > 0 || defline.replaced_id_bytes > 0)
    {
        return Err(Stop::refused(format!(
            "sequence {oid}: a definition line with bytes NCBI reads as '#' (it writes an \
             error each time it reads them)"
        )));
    }
    Ok(())
}

/// The retrieval's output and its diagnostics, NCBI's `cout` and the `cerr` tied to it.
pub(super) struct Output<'a, 'b> {
    out: &'a mut dyn ReportWrite,
    diag: &'a mut Diag<'b>,
}

// NCBI's `cerr` is tied to `cout` (the C++ standard's `cerr.tie()`): each line of diagnostics
// comes after NCBI's writes of the output so far (`Diag`), and on a closed pipe SIGPIPE ends
// NCBI's process at that write, before the line.
impl Output<'_, '_> {
    /// `Err(Stop::ClosedPipe)` where a write met a closed pipe before a text.
    fn closed(&self) -> Result<(), Stop> {
        if self.diag.closed() {
            return Err(Stop::ClosedPipe);
        }
        Ok(())
    }

    /// The error that the retrieval posts (`ERR_POST` at severity `Error`, with `message`),
    /// after NCBI's writes of the output; a closed pipe ends the retrieval without the line.
    pub fn post(&mut self, message: impl AsRef<[u8]>) -> Result<(), Stop> {
        self.diag.err_post_bytes(message.as_ref());
        self.closed()
    }

    /// NCBI's warning of a UTF8String read for the VisibleString member `path`
    /// (`Diag::post_utf8_string`), after the output so far, as `post`.
    pub fn post_utf8_string(&mut self, path: &str) -> Result<(), Stop> {
        self.diag.post_utf8_string(path);
        self.closed()
    }

    /// Runs `f` with a direct inserter into the underlying report stream (one lock per record,
    /// zero intermediate `CutBuffer` copy).
    pub fn with_inserter(
        &mut self,
        mut f: impl FnMut(&mut StreamInserter<'_>) -> std::io::Result<()>,
    ) -> Result<(), Stop> {
        with_stream_inserter(&mut *self.out, &mut f).map_err(Stop::write)
    }

    /// NCBI's insertions and flushes `text` (`CutBuffer`), emptied.
    pub fn insertions(&mut self, text: &mut CutBuffer) -> Result<(), Stop> {
        text.write_to(self.out).map_err(Stop::write)
    }

    /// The write of `object` to `m_Out` through `MSerial_AsnText`: the ASN.1 text document
    /// `document`, with NCBI's post of each byte its text stream replaces
    /// (`bioseq::write_asn_text`). The object stream's `Close` writes its buffer and does not
    /// flush `m_Out`.
    pub fn asn_text<D: AsnTextDocument>(&mut self, document: &D) -> Result<(), Stop> {
        let diag = &mut *self.diag;
        write_asn_text(document, &mut *self.out, |replaced| {
            diag.post_replaced(std::slice::from_ref(replaced));
            if diag.closed() {
                return Err(std::io::ErrorKind::BrokenPipe.into());
            }
            Ok(())
        })
    }

    /// The line of `x_Print` whose one field is the ASN.1 text document `document`, which NCBI
    /// writes into the field's string first: the posts of the bytes its text stream replaces, then
    /// the field's insertion, the last separator (empty) and `endl`.
    pub fn asn_text_line<D: AsnTextDocument>(&mut self, document: &D) -> Result<(), Stop> {
        let mut line = Vec::new();
        let mut replaced = Vec::new();
        write_asn_text(document, &mut line, |byte| {
            replaced.push(byte.clone());
            Ok(())
        })?;
        if !replaced.is_empty() {
            self.diag.post_replaced(&replaced);
            self.closed()?;
        }
        let mut text = CutBuffer::default();
        text.insertion(&line);
        text.endl();
        self.insertions(&mut text)
    }
}

impl Write for Output<'_, '_> {
    /// One insertion of NCBI's.
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.out.write_all(buf)?;
        Ok(buf.len())
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        self.out.write_all(buf)
    }

    /// LOSAT's flushes hold (`files::OutputFile`).
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

fn slice(letters: &[u8], range: Option<(u32, u32)>) -> &[u8] {
    match range {
        Some((from, to)) => &letters[from as usize..=to as usize],
        None => letters,
    }
}

/// The formatter of the request (`CBlastDB_Formatter`).
enum Formatter {
    Fasta(FastaFormatter),
    Bioseq(BioseqFormatter),
    Fields(FormatSpec),
}

impl Formatter {
    /// Whether the formatter reads the definition lines of a sequence (the FASTA's `GetBioseq`,
    /// a `%` specifier's `GetHdr`), which build NCBI's OID list (`SeqDb::build_oid_list`).
    fn reads_deflines(&self) -> bool {
        match self {
            Formatter::Fasta(_) | Formatter::Bioseq(_) => true,
            Formatter::Fields(spec) => spec.reads_deflines(),
        }
    }

    /// Whether the formatter writes letters that `-strand` turns (the FASTA, `%s`; the Bioseq
    /// of `%b` holds the stored data).
    fn reads_letters(&self) -> bool {
        match self {
            Formatter::Fasta(_) => true,
            Formatter::Bioseq(_) => false,
            Formatter::Fields(spec) => spec.writes_letters(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:613-628
    /// `Write(-1, config)`: the OID of a PIG that no volume holds. The FASTA and Bioseq
    /// formatters' `GetBioseq` throws `OID not found`; the `%` formatter's first read does
    /// (`FormatSpec::write_missing_oid`).
    fn write_missing(&mut self, out: &mut Output<'_, '_>) -> Result<(), Stop> {
        match self {
            Formatter::Fasta(_) | Formatter::Bioseq(_) => {
                Err(Stop::Ncbi(OID_NOT_FOUND.to_string()))
            }
            Formatter::Fields(spec) => spec.write_missing_oid(out),
        }
    }

    fn write(
        &mut self,
        source: &Source,
        oid: u32,
        config: &Config,
        target: Option<&[u8]>,
        out: &mut Output<'_, '_>,
    ) -> Result<(), Stop> {
        match self {
            Formatter::Fasta(fasta) => fasta.write(source, oid, config, target, out),
            Formatter::Bioseq(bioseq) => bioseq.write(source, oid, target, out),
            Formatter::Fields(spec) => spec.write(source, oid, config, target, out),
        }
    }
}

/// `CSeqDB::kOidNotFound` (seqdb.cpp:53).
pub(super) const OID_NOT_FOUND: &str = "OID not found";

/// The exception of NCBI's `CRef` read of no object (`ncbi_exception::NULL_POINTER`), whose
/// message `x_ProcessSearchRequest` posts: a read of the null set of a stored header of no bytes
/// ([`Source::null_set`], `GetBioseq` of a target).
pub(super) fn null_pointer() -> Stop {
    Stop::Ncbi(crate::ncbi_exception::NULL_POINTER.to_string())
}

/// `NStr::StringToNumeric<Int8>` of `text` with `fConvErr_NoThrow` and `errno` clear: the number of
/// a text that is one (an optional sign and digits that fit 64 bits).
pub(super) fn numeric_target(text: &[u8]) -> Option<i64> {
    args::string_to_int8_message(text).ok()
}

/// `CSeq_id` of `text` with the flags `fParse_PartialOK` and `fParse_Default`: the Seq-id of a
/// `-target_only` target that is not a number; the message of NCBI's exception as the request
/// reports it.
pub(super) fn parse_target(text: &[u8]) -> Result<ParsedSeqId, Stop> {
    parse_target_with(text, PARTIAL_OK | DEFAULT)
}

/// `CSeq_id` of `text` with `flags`, for a `-target_only` target (`parse_target`).
pub(super) fn parse_target_with(text: &[u8], flags: u32) -> Result<ParsedSeqId, Stop> {
    seq_id_from_text(text, flags, &mut |_| {}).map_err(|error| Stop::NcbiBytes(error.message))
}

/// Runs the retrieval of `args` on the open database `db` (`x_ProcessSearchRequest`), the
/// output to `out`, the diagnostics through `diag`; the exit status. `tax_ids` is the list of
/// `-taxids` / `-taxidlist` (`m_TaxIdList`; with `-target_only` the database is restricted to
/// it), `dbtype_protein` NCBI's `m_DbIsProtein`, `stdin` what a batch file `-` reads.
///
/// A failed write leaves NCBI's `ostream` bad and the run goes on (`OutputFile`): the later
/// entries are still looked up and their errors posted, and the exit status is theirs (NCBI
/// 2.17.0 with `-out /dev/full` and an id that is not found after it: `Skipped`, exit status 1;
/// review 2 of line RG item K1, F3). A closed pipe ends the retrieval at once, with exit
/// status 0 and no message, where SIGPIPE ends NCBI's (`Stop::ClosedPipe`): at NCBI's write
/// that meets it, which a diagnostic makes first (`Diag`, `cerr`'s tie to `cout`), so the
/// messages before that write are written and none after.
pub(super) fn process(
    db: &SeqDb,
    args: &Args,
    diag: &mut Diag,
    stdin: &mut Cin<'_>,
    out: &mut OutputFile<'_>,
    tax_ids: Option<&BTreeSet<i32>>,
    dbtype_protein: bool,
) -> i32 {
    let source = Source::new(db, dbtype_protein);
    let result = request(args, &source, tax_ids, diag, stdin, &mut *out);
    if out.closed() || diag.closed() {
        return 0;
    }
    match result {
        Ok(false) => 0,
        Ok(true) => 1,
        Err(Stop::Ncbi(message)) => {
            diag.err_post(&message);
            1
        }
        Err(Stop::NcbiBytes(message)) => {
            diag.err_post_bytes(&message);
            1
        }
        Err(Stop::Refused(reason)) => diag.refuse(&reason),
        Err(Stop::ClosedPipe) => 0,
        // `Run::write` posts it and writes the record again: it ends no retrieval.
        Err(Stop::Utf8Post(path)) => {
            debug_assert!(false, "a post that no write took");
            diag.post_utf8_string(&path);
            0
        }
    }
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:869-897
/// The request: whether an entry was not found (`err_found`), or what stopped it.
fn request(
    args: &Args,
    source: &Source,
    tax_ids: Option<&BTreeSet<i32>>,
    diag: &mut Diag,
    stdin: &mut Cin<'_>,
    out: &mut dyn ReportWrite,
) -> Result<bool, Stop> {
    let (mut config, kind) = init_search_request(args, source, diag)?;
    let line_length = args
        .integer("line_length")
        .and_then(|value| u32::try_from(value).ok())
        .unwrap_or(80);
    let mut formatter = match kind {
        OutputKind::Fasta => {
            Formatter::Fasta(FastaFormatter::new(line_length, args.given("long_seqids")))
        }
        OutputKind::Bioseq => Formatter::Bioseq(BioseqFormatter),
        OutputKind::Fields(spec, mask_algos) => match FormatSpec::parse(&spec, mask_algos) {
            Ok(spec) => Formatter::Fields(spec),
            Err(message) => return Err(Stop::Ncbi(message)),
        },
    };
    let ids = DbIds::of(source, args)?;
    let target_only = args.given("target_only");
    let mut run = Run {
        source,
        ids: &ids,
        formatter: &mut formatter,
        target_only,
        stdin,
        output: Output { out, diag },
    };
    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:819-847
    // A taxid list excludes `-entry` and `-entry_batch` but not `-ipg` and `-ipg_batch`
    // (its exclusion names `pig`, which is no argument): those then go first.
    let (entry, entry_batch) = (args.bytes("entry"), args.os_value("entry_batch"));
    match (entry.as_deref(), entry_batch) {
        (Some(b"all"), _) => {
            run.dump_all(&config)?;
            Ok(false)
        }
        (_, Some(file)) => {
            if args.given("get_dups") {
                run.batch_entry(file, &mut config)
            } else {
                run.batch_entry_no_dup(file, &mut config)
            }
        }
        (Some(entry), None) => run.entry(entry, &config),
        (None, None) => {
            if let Some(pig) = args.value("ipg") {
                // The reader took a value of `int` (`AsInteger`).
                let pig = args.integer("ipg").ok_or_else(|| {
                    Stop::refused(format!("-ipg {pig} is not a number of 32 bits"))
                })?;
                run.ipg(pig, &config)
            } else if let Some(file) = args.os_value("ipg_batch") {
                run.batch_pig(file, &mut config)
            } else if let Some(tax_ids) = tax_ids {
                run.tax_id_list(tax_ids, &config)
            } else {
                Err(Stop::Ncbi(
                    "Must specify query type: one of 'entry', 'entry_batch', or 'pig'".to_string(),
                ))
            }
        }
    }
}

/// What `x_InitSearchRequest` chooses: `m_FASTA`, `m_Asn1Bioseq`, or the format text of the
/// `%` formatter.
enum OutputKind {
    Fasta,
    Bioseq,
    /// The format text without the digits of its `%m`, and the algorithm of each `%m`.
    Fields(Vec<u8>, Vec<i32>),
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:736-817
// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:258-270
/// `s_IsMaskAlgoIdValid` before the retrieval writes anything: whether the database has the
/// algorithm `id` (a negative id is not looked up), with NCBI's error line for one it does not
/// have.
fn is_mask_algo_id_valid(source: &Source, id: i32, diag: &mut Diag) -> Result<bool, Stop> {
    match invalid_mask_algo_id(source, id)? {
        Some(line) => {
            diag.err_post(&line);
            Ok(false)
        }
        None => Ok(true),
    }
}

/// NCBI's error line of `s_IsMaskAlgoIdValid` for the algorithm `id`, or `None` when the
/// database has it (or `id` is negative).
fn invalid_mask_algo_id(source: &Source, id: i32) -> Result<Option<String>, Stop> {
    if id < 0 || source.db.mask_algorithms()?.ids().contains(&id) {
        return Ok(None);
    }
    Ok(Some(format!("Invalid filtering algorithm ID: {id}")))
}

/// The number of the digits at `at` of `outfmt` (`None` without a digit), which it removes.
/// NCBI adds them up in an `int` and overflows; LOSAT refuses a number that does not fit.
fn take_mask_algo_id(outfmt: &mut Vec<u8>, at: usize) -> Result<Option<i32>, Stop> {
    let end = outfmt[at..]
        .iter()
        .position(|byte| !byte.is_ascii_digit())
        .map_or(outfmt.len(), |len| at + len);
    if end == at {
        return Ok(None);
    }
    let id = outfmt[at..end].iter().try_fold(0i32, |id, &digit| {
        id.checked_mul(10)?.checked_add(i32::from(digit - b'0'))
    });
    let Some(id) = id else {
        return Err(Stop::refused(format!(
            "-outfmt %m{} (an algorithm id that overflows NCBI's int)",
            String::from_utf8_lossy(&outfmt[at..end])
        )));
    };
    outfmt.drain(at..end);
    Ok(Some(id))
}

/// Whether the bytes `text` hold `spec` (`std::string::find`).
fn contains(text: &[u8], spec: &[u8]) -> bool {
    text.windows(spec.len()).any(|window| window == spec)
}

/// The algorithm of each `%m` of `outfmt`, whose digits it removes. NCBI reads the first
/// `%m` of the text only, also inside `%%m`, and takes its algorithm for every `%m`; the
/// digits of the others stay in the text. Each `%m<digits>` has its own algorithm here, as
/// in the runs with one `%m` each (line BDC DECISIONS 12,
/// `PD-LOSAT-BLASTDBCMD-DEVIATIONS` item 4); a `%m` without digits has the first one's.
fn read_mask_specifiers(
    outfmt: &mut Vec<u8>,
    source: &Source,
    diag: &mut Diag,
) -> Result<Vec<i32>, Stop> {
    let Some(first_at) = outfmt.windows(2).position(|window| window == b"%m") else {
        return Ok(Vec::new());
    };
    let invalid = || Stop::Ncbi("Invalid filtering algorithm ID for outfmt %m.".to_string());
    let Some(first) = take_mask_algo_id(outfmt, first_at + 2)? else {
        return Err(Stop::Ncbi(
            "The option '-outfmt %m' is not followed by a masking algo ID.".to_string(),
        ));
    };
    if !is_mask_algo_id_valid(source, first, diag)? {
        return Err(invalid());
    }
    let mut algos = Vec::new();
    let mut at = 0;
    while at < outfmt.len() {
        let percent = at;
        if outfmt[percent] != b'%' {
            at += 1;
            continue;
        }
        at += 2;
        if outfmt.get(percent + 1) != Some(&b'm') {
            continue;
        }
        if percent == first_at {
            algos.push(first);
        } else if let Some(id) = take_mask_algo_id(outfmt, at)? {
            if !is_mask_algo_id_valid(source, id, diag)? {
                return Err(invalid());
            }
            algos.push(id);
        } else {
            algos.push(first);
        }
    }
    Ok(algos)
}

/// `x_InitSearchRequest`: the configuration and the kind of output.
fn init_search_request(
    args: &Args,
    source: &Source,
    diag: &mut Diag,
) -> Result<(Config, OutputKind), Stop> {
    let mut config = Config::default();
    let mut kind = OutputKind::Fasta;
    if let Some(given) = args.bytes("outfmt") {
        let mut outfmt = given.into_owned();
        let has = |text: &[u8], spec: &str| contains(text, spec.as_bytes());
        if (has(&outfmt, "%f") && (has(&outfmt, "%b") || has(&outfmt, "%d")))
            || (has(&outfmt, "%b") && has(&outfmt, "%d"))
        {
            return Err(Stop::Ncbi(
                "The %f, %b, %d output format options cannot be specified together.".to_string(),
            ));
        }
        kind = OutputKind::Fields(Vec::new(), Vec::new());
        if has(&outfmt, "%b") {
            outfmt = b"%b".to_vec();
            kind = OutputKind::Bioseq;
        }
        if has(&outfmt, "%f") {
            outfmt = b"%f".to_vec();
            kind = OutputKind::Fasta;
        }
        if has(&outfmt, "%d") {
            outfmt = b"%d".to_vec();
        }
        let mask_algos = read_mask_specifiers(&mut outfmt, source, diag)?;
        if let OutputKind::Fields(text, algos) = &mut kind {
            *text = outfmt;
            *algos = mask_algos;
        }
    }
    if !source.dbtype_protein && args.value("strand") == Some("minus") {
        config.strand_minus = true;
    }
    config.ctrl_a = args.given("ctrl_a");
    if let Some(value) = args.bytes("mask_sequence_with") {
        let id = args::string_to_int8_message(&value)
            .ok()
            .and_then(|id| i32::try_from(id).ok());
        // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:2261-2267
        // An unknown name is the remapper's error "Cannot find string algorithm id in algorithm
        // map."; a database with GI masks gives the names table of its GI-mask reader instead
        // (`MaskAlgorithms::id_of`, which reads the name's bytes).
        let id = match id {
            Some(id) => id,
            None => source
                .db
                .mask_algorithms()?
                .id_of(&value)
                .map_err(Stop::NcbiBytes)?,
        };
        if !is_mask_algo_id_valid(source, id, diag)? {
            return Err(Stop::Ncbi(
                "Invalid filtering algorithm ID for mask_sequence_with.".to_string(),
            ));
        }
        if id < -1 {
            return Err(Stop::refused(format!(
                "-mask_sequence_with {id} (a negative algorithm id, which NCBI accepts \
                 and then looks up) is not served"
            )));
        }
        config.filt_algo = id;
    }
    if let Some(value) = args.bytes("range") {
        if value.is_empty() {
            // `tokens.front()` of an empty vector: NCBI crashed (probe D61).
            return Err(Stop::refused(
                "-range \"\" (NCBI reads past an empty list of range elements) is not served"
                    .to_string(),
            ));
        }
        config.range = parse_sequence_range_open_end(value).map_err(Stop::Ncbi)?;
    }
    Ok((config, kind))
}

/// The lookup tables of the database, when the request looks ids up (the volume lists LOSAT
/// refuses do not stop `-entry all`).
struct DbIds<'a> {
    lookup: Option<IdLookup<'a>>,
    /// `GetBlastDbVersion()` is `eBDB_Version5`.
    version5: bool,
}

impl<'a> DbIds<'a> {
    fn of(source: &Source<'a>, args: &Args) -> Result<Self, Stop> {
        let db = source.db;
        let looks_up = args.given("entry_batch")
            || args.bytes("entry").is_some_and(|entry| *entry != *b"all")
            || args.given("ipg")
            || args.given("ipg_batch");
        Ok(DbIds {
            lookup: if looks_up {
                Some(db.id_lookup()?)
            } else {
                None
            },
            version5: db.is_version5(),
        })
    }
}

/// The answer of `x_GetOids` for one id.
enum Lookup {
    Oids(Vec<u32>),
    /// `Entry not found: <acc>` (the id as NCBI's lookup saw it).
    NotFound(Vec<u8>),
    /// An exception of the lookup with NCBI's text, or LOSAT's refusal.
    Stop(Stop),
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:152-180
/// The key a version 5 database is asked for: the label of an id with `|` or `_`, else the
/// text in upper case.
fn preprocess_accession_for_v5(id: &[u8]) -> Vec<u8> {
    if id.contains(&b'|') || id.contains(&b'_') {
        if let Ok(seqid) = seq_id_from_text(id, RAW_TEXT | ANY_LOCAL | PARTIAL_OK, &mut |_| {}) {
            let mut label = Vec::new();
            match &seqid {
                ParsedSeqId::Pir(_) | ParsedSeqId::Prf(_) => fasta(&seqid, &mut label),
                ParsedSeqId::Pdb(_) => content(&seqid, false, &mut label),
                _ => content(&seqid, true, &mut label),
            }
            return label;
        }
    }
    id.to_ascii_uppercase()
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:183-210
/// `x_GetOids` of `id`: the tables of the database (none for `-entry all`, which looks nothing up).
fn get_oids(ids: &DbIds, id: &[u8]) -> Lookup {
    match &ids.lookup {
        Some(lookup) if ids.version5 => get_oids_v5(lookup, id),
        Some(lookup) => get_oids_v4(lookup, id),
        None => Lookup::NotFound(id.to_vec()),
    }
}

/// `x_GetOids` of `id` of a version 5 database: a GI first when the key is a number, then the
/// accession table.
fn get_oids_v5(lookup: &IdLookup, id: &[u8]) -> Lookup {
    let acc = preprocess_accession_for_v5(id);
    if let Some(gi) = numeric_target(&acc) {
        match lookup.gi_to_oid(gi) {
            Ok(Some(oid)) => return Lookup::Oids(vec![oid]),
            Ok(None) => {}
            Err(error) => return Lookup::Stop(error.into()),
        }
    }
    match lookup.accession_to_oids(&acc) {
        Ok(oids) if oids.is_empty() => Lookup::NotFound(acc),
        Ok(oids) => Lookup::Oids(oids),
        Err(error) => Lookup::Stop(error.into()),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdb.cpp:870-893
/// `x_GetOids` of `id` of a version 4 database: a GI first when the id is a number, then the ISAMs
/// of the volumes. (The GI that `CSeqDB::AccessionToOids` tries when they find nothing is the
/// lookup `x_GetOids` has made before for a number, and finds nothing for any other text.)
fn get_oids_v4(lookup: &IdLookup, id: &[u8]) -> Lookup {
    if let Some(gi) = numeric_target(id) {
        match lookup.gi_to_oid(gi) {
            Ok(Some(oid)) => return Lookup::Oids(vec![oid]),
            Ok(None) => {}
            Err(error) => return Lookup::Stop(error.into()),
        }
    }
    match accession_to_oids_v4(lookup, id) {
        Ok(oids) if oids.is_empty() => Lookup::NotFound(id.to_vec()),
        Ok(oids) => Lookup::Oids(oids),
        Err(stop) => Lookup::Stop(stop),
    }
}

/// `CSeqDBImpl::AccessionToOids` of `id` of a version 4 database, with the exceptions of
/// `SeqDB_SimplifyAccession` and of the volumes as the request reports them.
fn accession_to_oids_v4(lookup: &IdLookup, id: &[u8]) -> Result<Vec<u32>, Stop> {
    let simplified = seqdb_simplify_accession(id).map_err(|error| {
        Stop::NcbiBytes(match error {
            SimplifyError::SeqId(error) => error.message,
            SimplifyError::Exception(exception) => exception.get_msg().to_vec(),
        })
    })?;
    lookup
        .accession_to_oids_v4(id, &simplified)?
        .ok_or_else(|| Stop::Ncbi(ID_OVERFLOWS.to_string()))
}

// NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:382-386
/// The lines of an `-entry_batch` (`-ipg_batch`) file (`args::input_lines`): the id and the
/// rest of the line (`NStr::SplitInTwo` of `line` at a space or a tab, giving `id` and `format`,
/// with the flags `fSplit_MergeDelimiters` and `fSplit_Truncate`), lines without an id left out.
/// A `\r` that the platform's reading leaves in a line stays (Linux: every `\r`), as in NCBI's
/// reading.
fn read_batch(
    argument: &str,
    file: &OsStr,
    stdin: &mut Cin<'_>,
) -> Result<Vec<(Vec<u8>, Vec<u8>)>, Stop> {
    let bytes = read_input_file(argument, file, stdin).map_err(Stop::NcbiBytes)?;
    let delimiter = |c: &u8| *c == b' ' || *c == b'\t';
    let mut lines = Vec::new();
    for line in input_lines(&bytes, LineEnds::HOST) {
        let Some(start) = line.iter().position(|c| !delimiter(c)) else {
            continue;
        };
        let end = line
            .iter()
            .rposition(|c| !delimiter(c))
            .map_or(start, |at| at + 1);
        let trimmed = &line[start..end];
        let (id, format) = match trimmed.iter().position(delimiter) {
            Some(at) => {
                let rest = &trimmed[at..];
                let from = rest
                    .iter()
                    .position(|c| !delimiter(c))
                    .unwrap_or(rest.len());
                (&trimmed[..at], &rest[from..])
            }
            None => (trimmed, &b""[..]),
        };
        lines.push((id.to_vec(), format.to_vec()));
    }
    Ok(lines)
}

/// The retrieval in progress.
struct Run<'a, 'b, 'c> {
    source: &'a Source<'a>,
    ids: &'a DbIds<'a>,
    formatter: &'a mut Formatter,
    target_only: bool,
    /// The standard input, which a batch file `-` reads.
    stdin: &'a mut Cin<'c>,
    output: Output<'a, 'b>,
}

impl Run<'_, '_, '_> {
    /// The error that the retrieval posts (`ERR_POST` at severity `Error`, with `message`)
    /// ([`Output::post`]).
    fn post(&mut self, message: impl AsRef<[u8]>) -> Result<(), Stop> {
        self.output.post(message)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:457-462
    // The minus strand reaches a protein database only through a taxid list read with
    // `-dbtype guess` (`m_DbIsProtein` false): NCBI then writes the protein letters of the
    // FASTA in reverse order (probe J114 of line BDC), and reverse-complements them as
    // nucleotides for `%s`. LOSAT refuses those; the other `%` formats do not read the
    // strand (K08).
    fn write(&mut self, oid: u32, config: &Config, target: Option<&[u8]>) -> Result<(), Stop> {
        if config.strand_minus && self.source.protein && self.formatter.reads_letters() {
            return Err(Stop::refused(
                "-strand minus on a protein database read with a taxid list and -dbtype guess \
                 (NCBI writes the protein letters reversed)"
                    .to_string(),
            ));
        }
        // `oid` is NCBI's OID, which the formatters read as NCBI's do (`Source::at`): one NCBI's
        // check keeps (`SeqDb::check_oid`, `SeqDb::oid_walk`), or, for `-ipg`, the OID of a
        // PIG whose volume holds it (`PigToOid`), kept or not.
        if self.formatter.reads_deflines() {
            self.source.db.build_oid_list()?;
        }
        match self
            .formatter
            .write(self.source, oid, config, target, &mut self.output)
        {
            Err(Stop::Utf8Post(path)) => self.write_after_utf8_post(&path, oid, config, target),
            written => written,
        }
    }

    /// `x_GetOids` of `id` (`get_oids`) with NCBI's warning of a UTF8String that the version check
    /// of its lookup read for a VisibleString (`IdLookup::take_utf8`, `x_CheckVersions`), where
    /// NCBI posts it: after the output so far, before the lookup's own message.
    fn get_oids(&mut self, id: &[u8]) -> Result<Lookup, Stop> {
        let found = get_oids(self.ids, id);
        if let Some(note) = self.ids.lookup.as_ref().and_then(IdLookup::take_utf8) {
            self.post_utf8_read(note.at)?;
        }
        Ok(found)
    }

    /// NCBI's warning of a read of the header at `at` that took a UTF8String for a
    /// VisibleString, where the run has not posted it: after the output so far (`cerr`'s tie).
    #[cold]
    #[inline(never)]
    fn post_utf8_read(&mut self, at: SeqLocation) -> Result<(), Stop> {
        if let Some(path) = self.source.db.take_utf8_post_at(at) {
            self.output.post_utf8_string(&path)?;
        }
        Ok(())
    }

    /// `write` of the sequence `oid` whose first read posts NCBI's warning of a UTF8String
    /// (`Stop::Utf8Post` of the member `path`), before the record's text (the formatters read the
    /// definition lines before they write): after the output so far, as `cerr`'s tie flushes it;
    /// then the record again, which posts nothing (once per process).
    #[cold]
    #[inline(never)]
    fn write_after_utf8_post(
        &mut self,
        path: &str,
        oid: u32,
        config: &Config,
        target: Option<&[u8]>,
    ) -> Result<(), Stop> {
        self.output.post_utf8_string(path)?;
        self.formatter
            .write(self.source, oid, config, target, &mut self.output)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/blastdb_format/seq_formatter.cpp:411-416
    /// `-entry all`: the sequences of NCBI's loop (`SeqDb::oid_walk`) in the order of their
    /// OIDs. Its first check sets NCBI's OID list up: the sequences the database includes,
    /// and, where the iteration range of the only alias file of `-db` leaves listed OIDs out,
    /// also the first listed OID up to the range's end and the OID at the range's end (line RG
    /// D1a-2, survey A, pattern 3). The sequences are read in the order of their OIDs, a window
    /// of each volume file at a time (`seqdb::file::in_order_reads`).
    fn dump_all(&mut self, config: &Config) -> Result<(), Stop> {
        let db = self.source.db;
        let _in_order = crate::seqdb::file::in_order_reads();
        for oid in db.oid_walk()? {
            self.write(oid, config, None)?;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:213-257
    /// `-entry a,b,...`: whether an entry was not found.
    fn entry(&mut self, entry: &[u8], config: &Config) -> Result<bool, Stop> {
        let queries: Vec<&[u8]> = entry.split(|&byte| byte == b',').collect();
        let mut err_found = 0;
        for query in &queries {
            match self.get_oids(query)? {
                Lookup::Oids(oids) => {
                    for oid in oids {
                        self.write(oid, config, self.target_only.then_some(*query))?;
                    }
                }
                Lookup::NotFound(acc) => {
                    self.post([b"Entry not found: ".as_slice(), &acc].concat())?;
                    err_found += 1;
                }
                Lookup::Stop(stop) => return Err(stop),
            }
        }
        if err_found == queries.len() {
            return Err(Stop::Ncbi(
                "Entry or entries not found in BLAST database".to_string(),
            ));
        }
        Ok(err_found > 0)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:375-438
    /// `-entry_batch` without `-get_dups`: one sequence per id, the ids resolved all at once
    /// first (the exception of a version 5 database without `acc2oid` is swallowed there:
    /// the ids keep what the LMDB files before it answered; in a version 4 database the first
    /// OID of each id), a GI after that.
    fn batch_entry_no_dup(&mut self, file: &OsStr, config: &mut Config) -> Result<bool, Stop> {
        let lines = read_batch("entry_batch", file, self.stdin)?;
        let ids: Vec<Vec<u8>> = lines
            .iter()
            .map(|(id, _)| {
                if self.ids.version5 {
                    preprocess_accession_for_v5(id)
                } else {
                    id.clone()
                }
            })
            .collect();
        let mut oids = Vec::with_capacity(ids.len());
        let Some(lookup) = &self.ids.lookup else {
            unreachable!("a batch looks ids up");
        };
        if self.ids.version5 {
            if let Err(error) = lookup.accessions_to_oids(&ids, &mut oids) {
                let swallowed = error.ncbi_message().is_some_and(|message| {
                    message
                        .windows(NO_ACCESSION_INFO.len())
                        .any(|window| window == NO_ACCESSION_INFO.as_bytes())
                });
                if !swallowed {
                    return Err(error.into());
                }
            }
        } else {
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1362-1376
            // The version checks of the lookups read their candidates' headers before anything
            // is written (`IdLookup::take_utf8`: NCBI's warning of a UTF8String there).
            let resolved: Result<(), Stop> = ids.iter().try_for_each(|id| {
                oids.push(accession_to_oids_v4(lookup, id)?.first().copied());
                Ok(())
            });
            if let Some(note) = lookup.take_utf8() {
                self.post_utf8_read(note.at)?;
            }
            resolved?;
        }
        let mut err_found = false;
        for (i, id) in ids.iter().enumerate() {
            let mut oid = oids[i];
            if let (None, Some(gi)) = (oid, numeric_target(id)) {
                oid = lookup.gi_to_oid(gi)?;
            }
            let Some(oid) = oid else {
                err_found = true;
                self.post([b"Skipped ".as_slice(), id].concat())?;
                continue;
            };
            if self.modify_config_for_batch_entry(&lines[i].1, config)? != 0 {
                err_found = true;
                self.post([b"Skipped ".as_slice(), id].concat())?;
                continue;
            }
            self.write(oid, config, self.target_only.then_some(id.as_slice()))?;
        }
        Ok(err_found)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:441-485
    /// `-entry_batch -get_dups`: line by line, every sequence of each id.
    fn batch_entry(&mut self, file: &OsStr, config: &mut Config) -> Result<bool, Stop> {
        let lines = read_batch("entry_batch", file, self.stdin)?;
        let mut err_found = false;
        for (id, format) in &lines {
            if self.modify_config_for_batch_entry(format, config)? != 0 {
                err_found = true;
                self.post([b"Skipped ".as_slice(), id].concat())?;
                continue;
            }
            match self.get_oids(id)? {
                Lookup::Oids(oids) => {
                    for oid in oids {
                        self.write(oid, config, None)?;
                    }
                }
                Lookup::NotFound(acc) => {
                    self.post([b"Entry not found: ".as_slice(), &acc].concat())?;
                    err_found = true;
                    self.post([b"Skipped ".as_slice(), id].concat())?;
                }
                Lookup::Stop(stop) => return Err(stop),
            }
        }
        Ok(err_found)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:1111-1123
    /// The OID lookup of a PIG (`PigToOid`, which asks the volumes, not the OID list: the OID of
    /// a PIG that an alias file's range leaves out is written as any other, probe p4c_c of line
    /// RG D1a-2): `Ok(None)` for NCBI's -1, a refusal where NCBI's OID stays uninitialised (no
    /// volume has PIG files: probes J204, J205, J212, J214 of line BDC, where the stack's value
    /// happened to be no OID).
    fn pig_to_oid(&self, pig: i32) -> Result<Option<u32>, Stop> {
        let lookup = self
            .ids
            .lookup
            .as_ref()
            .expect("an IPG request looks ids up");
        lookup.pig_to_oid(pig)?.ok_or_else(|| {
            Stop::refused(
                "-ipg and -ipg_batch on a database without PIG files (.ppi, .ppd): NCBI \
                 writes an OID it never set"
                    .to_string(),
            )
        })
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:219-222
    // A PIG no volume holds is written as the OID -1 (probes J200-J202).
    /// `-ipg`: the sequence of the PIG `pig`.
    fn ipg(&mut self, pig: i32, config: &Config) -> Result<bool, Stop> {
        match self.pig_to_oid(pig)? {
            Some(oid) => self.write(oid, config, None)?,
            None => self.formatter.write_missing(&mut self.output)?,
        }
        Ok(false)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:488-520
    // An id that is not an `int` is the PIG 0 (`fConvErr_NoThrow`; on Linux a `\r` of a CRLF
    // line stays on the id: J207, O42).
    /// `-ipg_batch`: the sequence of the PIG of each line of `file`, with the line's range,
    /// strand and mask algorithm.
    fn batch_pig(&mut self, file: &OsStr, config: &mut Config) -> Result<bool, Stop> {
        if self.source.db.is_restricted() {
            // NCBI writes the volume OID of each PIG through the database opened with the
            // list, whose header read drops the lines without the list's taxids (probe K38).
            return Err(Stop::refused(
                "-ipg_batch with -target_only and a taxid list".to_string(),
            ));
        }
        let lines = read_batch("ipg_batch", file, self.stdin)?;
        let mut err_found = false;
        for (id, format) in &lines {
            if self.modify_config_for_batch_entry(format, config)? != 0 {
                err_found = true;
                self.post([b"Skipped IPG : ".as_slice(), id].concat())?;
                continue;
            }
            let pig = args::string_to_int8_message(id)
                .ok()
                .and_then(|pig| i32::try_from(pig).ok())
                .unwrap_or(0);
            match self.pig_to_oid(pig)? {
                Some(oid) => self.write(oid, config, None)?,
                None => {
                    err_found = true;
                    self.post([b"Skipped IPG: ".as_slice(), id].concat())?;
                }
            }
        }
        Ok(err_found)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:304-318
    // With `-target_only` the database was opened with the list (`SeqDb::restrict_taxids`):
    // its OIDs are those of the list (the OIDs of its ordinal ids `0..n`, `SeqDb::oid`) and its
    // definition lines those with one of the list's taxids. Without it every definition line
    // of a sequence is written (probes T01-T09, J101-J139, S33-S37 of line BDC).
    /// `-taxids` / `-taxidlist`: the sequences of the taxids `tax_ids`, by OID.
    fn tax_id_list(&mut self, tax_ids: &BTreeSet<i32>, config: &Config) -> Result<bool, Stop> {
        let db = self.source.db;
        let oids: Vec<u32> = if db.is_restricted() {
            (0..db.num_seqs())
                .map(|ordinal| db.oid(ordinal).get())
                .collect()
        } else {
            db.tax_id_lookup()?.tax_ids_to_oids(tax_ids)?
        };
        if oids.is_empty() {
            self.post("No seq found in db for taxonomy list")?;
            return Ok(true);
        }
        for oid in oids {
            self.write(oid, config, None)?;
        }
        Ok(false)
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:271-300
    /// The range, strand and mask algorithm of a batch line (`id [range] [minus] [algo]`):
    /// 1 when the algorithm id is invalid (a word that is not a number is the id -1, which
    /// NCBI does not check).
    fn modify_config_for_batch_entry(
        &mut self,
        format: &[u8],
        config: &mut Config,
    ) -> Result<i32, Stop> {
        let mut status = 0;
        if !self.source.dbtype_protein {
            config.strand_minus = false;
        }
        config.range = NO_RANGE;
        config.filt_algo = -1;
        for token in format
            .split(|&byte| byte == b' ' || byte == b'\t')
            .filter(|token| !token.is_empty())
        {
            if token.contains(&b'-') {
                if let Ok(range) = parse_sequence_range_open_end(token) {
                    config.range = range;
                }
            } else if !self.source.dbtype_protein && token.eq_ignore_ascii_case(b"minus") {
                config.strand_minus = true;
            } else {
                config.filt_algo = string_to_non_negative_int(token);
                // A batch line comes after the output of the lines before it: the line goes
                // through `post` (`cerr`'s tie).
                if let Some(line) = invalid_mask_algo_id(self.source, config.filt_algo)? {
                    self.post(&line)?;
                    status = 1;
                }
            }
        }
        Ok(status)
    }
}

#[cfg(test)]
#[path = "retrieve_tests.rs"]
mod tests;
