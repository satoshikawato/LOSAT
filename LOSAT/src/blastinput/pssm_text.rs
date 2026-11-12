//! The reader of a PSSM in ASN.1 text (`PssmWithParameters`, the file of psiblast's `-out_pssm`
//! and of `-in_pssm`): NCBI's ASN.1 text object stream (`CObjectIStreamAsn`) for the types of
//! `report::serial::pssm`, whose writer gives the same bytes back for a PSSM that psiblast wrote.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1209-1216
//! ```c++
//!     if (args.Exist(kArgPSIInputChkPntFile) && args[kArgPSIInputChkPntFile]) {
//!         CNcbiIstream& in = args[kArgPSIInputChkPntFile].AsInputFile();
//!         _ASSERT(m_Pssm.Empty());
//!         m_Pssm.Reset(new CPssmWithParameters);
//!         try {
//!             switch (CFormatGuess().Format(in)) {
//!             case CFormatGuess::eBinaryASN:
//!                 in >> MSerial_AsnBinary >> *m_Pssm;
//! ```
//!
//! The values are read as NCBI's stream reads them: white space is spaces, tabs and ends of line
//! (`\n`, `\r`, `\r\n` or `\n\r`, each one line), and a comment runs from `--` to the next `--` or
//! the end of the line; the members of a SEQUENCE come in the order of the module; a REAL is
//! `{ mantissa, base, exponent }` or a word (`ReadDouble`, `read_real`); a VisibleString joins its
//! lines (an end of line inside it is dropped) and `""` is one quote; an OCTET STRING is hex digits
//! across lines. The members and variants that LOSAT does not keep are skipped as values
//! (`skip_value`): the descriptors of the query other than its titles, its `annot`, the members
//! of its `Seq-inst` other than `repr`, `mol`, `length` and `seq-data`, and the `constraints` and
//! `annotatedSites` of the parameters; `read_pssm_with_parameters_and_skips` names them, since NCBI
//! writes them back where it writes the PSSM it read (psiblast `-in_pssm` with `-out_pssm`, item
//! PS5b). An empty `descr` is NCBI's error (`CSeq_descr::PostRead`). A string keeps its bytes as
//! read (NCBI replaces a byte outside the visible range by `#` with an error message; a PSSM that
//! LOSAT writes has none, `algorithm::psiblast::pssm`). LOSAT reads the forms below and refuses the others with an error
//! naming its line: a query `Seq-entry` that is a `set`, a `Seq-id` other than `local`, `gi`,
//! `general` and the `Textseq-id` ones, `Seq-data` other than `iupacaa`, `ncbieaa` and
//! `ncbistdaa`, and a REAL mantissa of more than 18 digits.

#[cfg(test)]
use crate::report::serial::asn_text_reader::libm_pow10;
use crate::report::serial::asn_text_reader::{is_first_id_char, Reader, SliceSource};
use crate::report::serial::pssm::{
    FormatRpsDbParameters, ObjectId, Pssm, PssmFinalData, PssmIntermediateData, PssmParameters,
    PssmQuery, PssmSeqId, PssmWithParameters, SeqData, SeqInst, TextseqId, SEQ_INST_MOL,
    SEQ_INST_REPR, TEXTSEQ_VARIANTS,
};

/// A PSSM that the reader cannot read: the line (from 1) where it stopped and why (the error of
/// the shared ASN.1 text reader, `report::serial::asn_text_reader`).
pub type PssmTextError = crate::report::serial::asn_text_reader::AsnTextError;

type Result<T> = std::result::Result<T, PssmTextError>;

/// Reads the `PssmWithParameters` document at the start of `text` (what follows it is not read,
/// as `stream >> object` reads one object).
// NCBI reference (598d8ae6): c++/src/serial/objistrasn.cpp:709-721
// ```c++
// string CObjectIStreamAsn::ReadFileHeader()
// {
//     CTempString id = ReadTypeId(SkipWhiteSpace());
//     string s(id);
//     if ( SkipWhiteSpace() == ':' &&
//          m_Input.PeekCharNoEOF(1) == ':' &&
//          m_Input.PeekCharNoEOF(2) == '=' ) {
//         m_Input.SkipChars(3);
//     }
//     else
//         ThrowError(fFormatError, "'::=' expected");
//     return s;
// }
// ```
#[cfg(test)]
pub fn read_pssm_with_parameters(text: &[u8]) -> Result<PssmWithParameters> {
    read_pssm_with_parameters_and_skips(text).map(|(pssm, _)| pssm)
}

/// `read_pssm_with_parameters`, with the names of the members that the reader skipped (in the
/// order read; `"molinfo"` for a `molinfo` descriptor of the query, `"Seqdesc"` for another
/// descriptor other than a title), which NCBI's object keeps and writes back.
pub fn read_pssm_with_parameters_and_skips(
    text: &[u8],
) -> Result<(PssmWithParameters, Vec<&'static str>)> {
    let mut reader = Reader::new(SliceSource::new(text), "PSSM", "PSSM reader");
    let c = reader.skip_white_space()?;
    let id = reader.scan_end_of_id(is_first_id_char(c));
    let id = reader.src.text(id);
    if id != b"PssmWithParameters" {
        return Err(reader.error(format!(
            "the document is a {}, not a PssmWithParameters",
            String::from_utf8_lossy(id)
        )));
    }
    if reader.skip_white_space()? != b':' || reader.peek_at(1) != b':' || reader.peek_at(2) != b'='
    {
        return Err(reader.error("'::=' expected"));
    }
    reader.pos += 3;
    let pssm = reader.pssm_with_parameters()?;
    Ok((pssm, reader.skipped))
}

// The grammar of `PssmWithParameters` and of the types it holds, over the shared lexical layer.
impl<'a> Reader<SliceSource<'a>> {
    /// `ReadMemberId(SkipWhiteSpace())` (objistrasn.cpp:343-364): a member or variant name.
    fn member_id(&mut self) -> Result<&'a [u8]> {
        let id = self.member_id_range()?;
        Ok(self.src.text(id))
    }

    fn object_id(&mut self) -> Result<ObjectId> {
        let variant = self.member_id()?;
        match variant {
            b"id" => Ok(ObjectId::Id(self.int4()?)),
            b"str" => Ok(ObjectId::Str(self.string()?)),
            _ => Err(self.error(format!(
                "unexpected Object-id variant {}",
                String::from_utf8_lossy(variant)
            ))),
        }
    }

    fn pssm_with_parameters(&mut self) -> Result<PssmWithParameters> {
        let mut members = self.begin_members()?;
        let mut pssm = None;
        let mut params = None;
        while let Some(index) = self.next_member(&mut members, &["pssm", "params"])? {
            match index {
                0 => pssm = Some(self.pssm()?),
                _ => params = Some(self.params()?),
            }
        }
        Ok(PssmWithParameters {
            pssm: pssm.ok_or_else(|| self.error("member pssm missing"))?,
            params,
        })
    }

    fn pssm(&mut self) -> Result<Pssm> {
        const NAMES: [&str; 9] = [
            "isProtein",
            "identifier",
            "numRows",
            "numColumns",
            "rowLabels",
            "byRow",
            "query",
            "intermediateData",
            "finalData",
        ];
        let mut members = self.begin_members()?;
        let mut pssm = Pssm {
            is_protein: None,
            identifier: None,
            num_rows: 0,
            num_columns: 0,
            row_labels: None,
            by_row: None,
            query: None,
            intermediate_data: None,
            final_data: None,
        };
        let (mut rows, mut columns) = (None, None);
        while let Some(index) = self.next_member(&mut members, &NAMES)? {
            match index {
                0 => pssm.is_protein = Some(self.boolean()?),
                1 => pssm.identifier = Some(self.object_id()?),
                2 => rows = Some(self.int4()?),
                3 => columns = Some(self.int4()?),
                4 => pssm.row_labels = Some(self.list(Self::string)?),
                5 => pssm.by_row = Some(self.boolean()?),
                6 => pssm.query = Some(self.seq_entry()?),
                7 => pssm.intermediate_data = Some(self.intermediate_data()?),
                _ => pssm.final_data = Some(self.final_data()?),
            }
        }
        pssm.num_rows = rows.ok_or_else(|| self.error("member numRows missing"))?;
        pssm.num_columns = columns.ok_or_else(|| self.error("member numColumns missing"))?;
        Ok(pssm)
    }

    /// The query's `Seq-entry`: its `seq` variant, a Bioseq.
    fn seq_entry(&mut self) -> Result<PssmQuery> {
        let variant = self.member_id()?;
        if variant != b"seq" {
            return Err(self.error(format!(
                "a query Seq-entry {} is not supported by LOSAT's PSSM reader",
                String::from_utf8_lossy(variant)
            )));
        }
        let mut members = self.begin_members()?;
        let (mut ids, mut titles, mut inst) = (None, None, None);
        while let Some(index) = self.next_member(&mut members, &["id", "descr", "inst", "annot"])? {
            match index {
                0 => ids = Some(self.list(Self::seq_id)?),
                1 => {
                    // NCBI reference (598d8ae6): c++/src/objects/seq/Seq_descr.cpp:62-68
                    // ```c++
                    // void CSeq_descr::PostRead(void) const
                    // {
                    //     if ( !s_SeqDescrAllowEmpty->Get() && Get().empty() ) {
                    //         NCBI_THROW(CSerialException, eInvalidData,
                    //                    "empty Seq-descr is not allowed");
                    // ```
                    // (`OBJECTS_SEQ_DESCR_ALLOW_EMPTY` is false unless set in NCBI's environment or
                    // configuration file, which LOSAT does not read.)
                    let mut descriptors = 0usize;
                    let found: Vec<Vec<u8>> = self
                        .list(|reader| {
                            descriptors += 1;
                            reader.seqdesc_title()
                        })?
                        .into_iter()
                        .flatten()
                        .collect();
                    if descriptors == 0 {
                        return Err(self.error("empty Seq-descr is not allowed"));
                    }
                    titles = (!found.is_empty()).then_some(found);
                }
                2 => inst = Some(self.seq_inst()?),
                _ => self.skip_member("annot")?,
            }
        }
        Ok(PssmQuery {
            ids: ids.ok_or_else(|| self.error("member id missing"))?,
            titles,
            inst: inst.ok_or_else(|| self.error("member inst missing"))?,
        })
    }

    /// One `Seqdesc`: its text when it is a `title`; another variant is skipped.
    fn seqdesc_title(&mut self) -> Result<Option<Vec<u8>>> {
        let variant = self.member_id()?;
        if variant.is_empty() {
            return Err(self.error("choice variant id expected"));
        }
        if variant == b"title" {
            return Ok(Some(self.string()?));
        }
        // A `molinfo` descriptor takes the title out of the query's warning label
        // (`CBlastQuerySourceBioseqSet::GetTitle`, quoted at `algorithm::psiblast::in_pssm`).
        self.skip_member(if variant == b"molinfo" {
            "molinfo"
        } else {
            "Seqdesc"
        })?;
        Ok(None)
    }

    fn seq_id(&mut self) -> Result<PssmSeqId> {
        let variant = self.member_id()?;
        match variant {
            b"local" => Ok(PssmSeqId::Local(self.object_id()?)),
            b"gi" => Ok(PssmSeqId::Gi(self.int8()?)),
            b"general" => {
                let mut members = self.begin_members()?;
                let (mut db, mut tag) = (None, None);
                while let Some(index) = self.next_member(&mut members, &["db", "tag"])? {
                    match index {
                        0 => db = Some(self.string()?),
                        _ => tag = Some(self.object_id()?),
                    }
                }
                Ok(PssmSeqId::General {
                    db: db.ok_or_else(|| self.error("member db missing"))?,
                    tag: tag.ok_or_else(|| self.error("member tag missing"))?,
                })
            }
            _ => {
                let Some(&(name, tag)) = TEXTSEQ_VARIANTS
                    .iter()
                    .find(|(name, _)| name.as_bytes() == variant)
                else {
                    return Err(self.error(format!(
                        "a query Seq-id {} is not supported by LOSAT's PSSM reader",
                        String::from_utf8_lossy(variant)
                    )));
                };
                let mut members = self.begin_members()?;
                let mut id = TextseqId::default();
                while let Some(index) =
                    self.next_member(&mut members, &["name", "accession", "release", "version"])?
                {
                    match index {
                        0 => id.name = Some(self.string()?),
                        1 => id.accession = Some(self.string()?),
                        2 => id.release = Some(self.string()?),
                        _ => id.version = Some(self.int4()?),
                    }
                }
                Ok(PssmSeqId::Textseq {
                    variant: name,
                    tag,
                    id,
                })
            }
        }
    }

    fn seq_inst(&mut self) -> Result<SeqInst> {
        const NAMES: [&str; 9] = [
            "repr", "mol", "length", "fuzz", "topology", "strand", "seq-data", "ext", "hist",
        ];
        let mut members = self.begin_members()?;
        let (mut repr, mut mol, mut length, mut seq_data) = (None, None, None, None);
        while let Some(index) = self.next_member(&mut members, &NAMES)? {
            match index {
                0 => repr = Some(self.enumerated(&SEQ_INST_REPR)?),
                1 => mol = Some(self.enumerated(&SEQ_INST_MOL)?),
                2 => length = Some(self.int4()?),
                6 => seq_data = Some(self.seq_data()?),
                other => self.skip_member(NAMES[other])?,
            }
        }
        Ok(SeqInst {
            repr: repr.ok_or_else(|| self.error("member repr missing"))?,
            mol: mol.ok_or_else(|| self.error("member mol missing"))?,
            length,
            seq_data,
        })
    }

    fn seq_data(&mut self) -> Result<SeqData> {
        let variant = self.member_id()?;
        match variant {
            b"iupacaa" => Ok(SeqData::Iupacaa(self.string()?)),
            b"ncbieaa" => Ok(SeqData::Ncbieaa(self.string()?)),
            b"ncbistdaa" => Ok(SeqData::Ncbistdaa(self.octets()?)),
            _ => Err(self.error(format!(
                "query Seq-data {} is not supported by LOSAT's PSSM reader",
                String::from_utf8_lossy(variant)
            ))),
        }
    }

    fn intermediate_data(&mut self) -> Result<PssmIntermediateData> {
        const NAMES: [&str; 9] = [
            "resFreqsPerPos",
            "weightedResFreqsPerPos",
            "freqRatios",
            "informationContent",
            "gaplessColumnWeights",
            "sigma",
            "intervalSizes",
            "numMatchingSeqs",
            "numIndeptObsr",
        ];
        let mut members = self.begin_members()?;
        let mut data = PssmIntermediateData::default();
        let mut freq_ratios = None;
        while let Some(index) = self.next_member(&mut members, &NAMES)? {
            match index {
                0 => data.res_freqs_per_pos = Some(self.list(Self::int4)?),
                1 => data.weighted_res_freqs_per_pos = Some(self.list(Self::real)?),
                2 => freq_ratios = Some(self.list(Self::real)?),
                3 => data.information_content = Some(self.list(Self::real)?),
                4 => data.gapless_column_weights = Some(self.list(Self::real)?),
                5 => data.sigma = Some(self.list(Self::real)?),
                6 => data.interval_sizes = Some(self.list(Self::int4)?),
                7 => data.num_matching_seqs = Some(self.list(Self::int4)?),
                _ => data.num_indept_obsr = Some(self.list(Self::real)?),
            }
        }
        data.freq_ratios = freq_ratios.ok_or_else(|| self.error("member freqRatios missing"))?;
        Ok(data)
    }

    fn final_data(&mut self) -> Result<PssmFinalData> {
        const NAMES: [&str; 9] = [
            "scores",
            "lambda",
            "kappa",
            "h",
            "scalingFactor",
            "lambdaUngapped",
            "kappaUngapped",
            "hUngapped",
            "wordScoreThreshold",
        ];
        let mut members = self.begin_members()?;
        let (mut scores, mut lambda, mut kappa, mut h) = (None, None, None, None);
        let mut data = PssmFinalData {
            scores: Vec::new(),
            lambda: 0.0,
            kappa: 0.0,
            h: 0.0,
            scaling_factor: None,
            lambda_ungapped: None,
            kappa_ungapped: None,
            h_ungapped: None,
            word_score_threshold: None,
        };
        while let Some(index) = self.next_member(&mut members, &NAMES)? {
            match index {
                0 => scores = Some(self.list(Self::int4)?),
                1 => lambda = Some(self.real()?),
                2 => kappa = Some(self.real()?),
                3 => h = Some(self.real()?),
                4 => data.scaling_factor = Some(self.int4()?),
                5 => data.lambda_ungapped = Some(self.real()?),
                6 => data.kappa_ungapped = Some(self.real()?),
                7 => data.h_ungapped = Some(self.real()?),
                _ => data.word_score_threshold = Some(self.real()?),
            }
        }
        data.scores = scores.ok_or_else(|| self.error("member scores missing"))?;
        data.lambda = lambda.ok_or_else(|| self.error("member lambda missing"))?;
        data.kappa = kappa.ok_or_else(|| self.error("member kappa missing"))?;
        data.h = h.ok_or_else(|| self.error("member h missing"))?;
        Ok(data)
    }

    fn params(&mut self) -> Result<PssmParameters> {
        const NAMES: [&str; 6] = [
            "pseudocount",
            "rpsdbparams",
            "constraints",
            "bitScoreThresh",
            "bitScoreReportingThresh",
            "annotatedSites",
        ];
        let mut members = self.begin_members()?;
        let mut params = PssmParameters::default();
        while let Some(index) = self.next_member(&mut members, &NAMES)? {
            match index {
                0 => params.pseudocount = Some(self.int4()?),
                1 => params.rpsdbparams = Some(self.rpsdbparams()?),
                3 => params.bit_score_thresh = Some(self.real()?),
                4 => params.bit_score_reporting_thresh = Some(self.real()?),
                other => self.skip_member(NAMES[other])?,
            }
        }
        Ok(params)
    }

    fn rpsdbparams(&mut self) -> Result<FormatRpsDbParameters> {
        let mut members = self.begin_members()?;
        let (mut matrix_name, mut gap_open, mut gap_extend) = (None, None, None);
        while let Some(index) =
            self.next_member(&mut members, &["matrixName", "gapOpen", "gapExtend"])?
        {
            match index {
                0 => matrix_name = Some(self.string()?),
                1 => gap_open = Some(self.int4()?),
                _ => gap_extend = Some(self.int4()?),
            }
        }
        Ok(FormatRpsDbParameters {
            matrix_name: matrix_name.ok_or_else(|| self.error("member matrixName missing"))?,
            gap_open,
            gap_extend,
        })
    }
}

#[cfg(test)]
#[path = "pssm_text_tests.rs"]
mod tests;
