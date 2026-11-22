//! The queries of one archive with their domain hits: `CNcbiCdAlignProcessor::ProcessBlastResults`
//! (each query `Bioseq` with the alignments that follow it in the list; a nucleotide query's
//! hits classified per reading frame), the query's label, title, length and residues
//! (`FillSequenceFromBioseq`), and the masks of the results
//! (`ConvertBlastMaskListToSeqLocInfoVector`, which rpsbproc reads but does not print).

use super::archive::{Align, Mask, MaskLocation, Query, QueryRowId, Segs, SeqEntry};
use super::data::{compute_pseudo_pssm_id, CddData};
use super::domains::{
    calculate, parse_align, parse_std_align, sort_reading_frames, AlignIndex, DomAlign,
};
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::seqid::defline_set::{compare, Sic};
use crate::seqid::model::{ObjectId, ParsedSeqId};

/// A failure while the queries of an archive are processed.
#[derive(Debug)]
pub(crate) enum ProcessError {
    /// An exception of NCBI's serial objects thrown while a session is processed
    /// (`CInvalidChoiceSelection::eFail` from a getter of a CHOICE that holds another variant,
    /// `CUnassignedMember::eGet` from a getter of an unset member; rpsbproc: two `Error:` lines,
    /// exit 3).
    Exception(NcbiException),
    /// A `std::exception` that leaves `Run` (`CSimpleException` of `THROW_SIMPLE`): its
    /// `what()`, which NCBI's application layer reports (exit 3).
    Std(String),
    /// What LOSAT does not reproduce (an explicit refusal).
    Unsupported(String),
}

/// The `CInvalidChoiceSelection` message of `C<type>::Get<variant>()` on a CHOICE of
/// `module::type` that holds `current`.
// NCBI reference (598d8ae6): c++/src/serial/exception.cpp:122-150
// NCBI reference (598d8ae6): c++/src/serial/exception.cpp:86-118 (`SPrintIdentifier`: `-` is
// `_`, the first letter and each one after `.` upper-cased)
// (NCBI's rpsbproc printed `CSeq_id::GetGeneral(): Invalid choice selection:
// NCBI-Seqloc::Seq-id.local` for a profile database with `lcl|` ids, PS16a probe 2.)
// (The exception is thrown by the generated `C<type>_Base::ThrowInvalidSelection`, whose place
// LOSAT does not have: `ThrowSite::generated`.)
pub(crate) fn choice_selection(
    module: &str,
    type_name: &str,
    required: &str,
    current: &str,
) -> ProcessError {
    ProcessError::Exception(NcbiException::new(
        ThrowSite::generated("CInvalidChoiceSelection::eFail"),
        format!(
            "C{}::Get{}(): Invalid choice selection: {module}::{type_name}.{current}",
            identifier(type_name),
            identifier(required)
        ),
    ))
}

/// `SPrintIdentifier` of a type or member name without `.`: `-` read as `_`, the first letter
/// upper-cased.
fn identifier(name: &str) -> String {
    let mut out = String::new();
    let mut capitalize = true;
    for c in name.chars() {
        let c = if c == '-' { '_' } else { c };
        if capitalize {
            out.extend(c.to_uppercase());
            capitalize = false;
        } else {
            out.push(c);
        }
    }
    out
}

/// The `CUnassignedMember` exception of `C<type>::Get<member>()` on a SEQUENCE of
/// `module::type` whose `member` is not set.
// NCBI reference (598d8ae6): c++/src/serial/serialobject.cpp:266-306
// (`x_GetVerifyData` is `eSerialVerifyData_Yes` unless `SERIAL_VERIFY_DATA_GET` says otherwise,
// which LOSAT's rpsbproc rejects: `reader_variable`.)
pub(crate) fn unassigned_member(module: &str, type_name: &str, member: &str) -> ProcessError {
    ProcessError::Exception(NcbiException::new(
        UNASSIGNED_MEMBER,
        format!(
            "C{}::Get{}(): Attempt to get unassigned member {module}::{type_name}.{member}",
            identifier(type_name),
            identifier(member)
        ),
    ))
}

// NCBI reference (598d8ae6): c++/src/serial/serialobject.cpp:317-320
// NCBI reference (598d8ae6): c++/src/serial/datatool/classstr.cpp:849-855
// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:805-810
// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:836-842
// The generated getters call `ThrowUnassigned` with the index, so the place is the `NCBI_THROW` of
// serialobject.cpp:306 (`file_name` and `file_line` are 0). Its function is a `const` method,
// whose `__PRETTY_FUNCTION__` ends with `const`, not `)`: `ParseCurrFunctName` names no
// function then.
const UNASSIGNED_MEMBER: ThrowSite = ThrowSite {
    file: "src/serial/serialobject.cpp",
    line: 306,
    kind: "CUnassignedMember::eGet",
    function: None,
};

/// A query of the output (`TDomQuery`) with its hits.
#[derive(Debug)]
pub(crate) struct DomQuery {
    /// `m_strAccession`: the label of the query's best id.
    pub(crate) accession: Vec<u8>,
    /// `m_strTitle`.
    pub(crate) title: Vec<u8>,
    /// `m_uiSeqLen`.
    pub(crate) length: u32,
    /// `m_bIsNa`: `inst.mol` is not `aa`.
    pub(crate) is_na: bool,
    /// `m_strSeqData`: the residues (or bases) as one-letter codes (`Get1LetterSeqData`); empty
    /// when the query has neither `seq-data` nor literals of `ext delta` with data.
    pub(crate) residues: Vec<u8>,
    pub(crate) aligns: Vec<DomAlign>,
    /// `m_dimSplitAligns`: the classification of the hits of each reading frame (index 0 only
    /// for a protein query).
    pub(crate) frames: [AlignIndex; 6],
    /// `_TOfflDomQuery::m_dimTranslated` of a query of the XML input: the residues that its
    /// alignments give each reading frame (`-` elsewhere), which the sites of a nucleotide query
    /// read in place of `TranslateAll` of its bases; `None` for a query of an archive.
    pub(crate) translated: Option<Box<[Vec<u8>; 6]>>,
}

/// `ConvertBlastMaskListToSeqLocInfoVector` as far as it can fail: the first location of each
/// mask that has locations must be a `packed-int`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_blast.cpp:264-300
// The masks are kept with the queries (`m_vecMaskedRegions`), which the text output never prints;
// what is observable is the exception of a location that is not a `packed-int`. A location
// without one id (`GetId` gives no id: `null`, `feat`, `mix`, `equiv`, `bond`, or intervals of
// different ids) is dereferenced by NCBI: refused.
pub(crate) fn check_masks(masks: &[Mask]) -> Result<(), ProcessError> {
    for mask in masks {
        let Some(first) = &mask.first else {
            continue;
        };
        match first {
            MaskLocation::PackedInt(ids) => {
                if ids.is_empty() || ids.iter().any(|id| id != &ids[0]) {
                    return Err(ProcessError::Unsupported(
                        "a mask whose packed-int location is empty or names more than one Seq-id is not supported by LOSAT's rpsbproc".to_string(),
                    ));
                }
            }
            MaskLocation::Other(name @ ("int" | "whole" | "empty" | "pnt" | "packed-pnt")) => {
                return Err(choice_selection(
                    "NCBI-Seqloc",
                    "Seq-loc",
                    "packed-int",
                    name,
                ));
            }
            MaskLocation::Other(name) => {
                return Err(ProcessError::Unsupported(format!(
                    "a mask location of the Seq-loc variant {name} is not supported by LOSAT's rpsbproc"
                )));
            }
        }
    }
    Ok(())
}

/// `CSeq_id::GetStringDescr` of the bioseq with `eFormat_BestWithVersion`, and the clean-up of
/// `FillSequenceFromBioseq`: the label of the best id, the text after its first `:` and then
/// after its first `|`.
// NCBI reference (598d8ae6): c++/src/objects/seqloc/Seq_id.cpp:2292-2345
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:1030-1046
// (The best id's label is never empty when the Bioseq has an id, so the local id's label is
// never used.)
fn query_accession(query: &Query) -> Vec<u8> {
    let is_na = query.mol != 3;
    let score = if is_na {
        crate::seqid::rank::fasta_na_score
    } else {
        crate::seqid::rank::fasta_aa_score
    };
    let mut label = Vec::new();
    if let Some(best) = crate::seqid::rank::find_best(&query.ids, score) {
        crate::seqid::label::both_version(best, &mut label);
    }
    if let Some(at) = label.iter().position(|&b| b == b':') {
        label.drain(..=at);
    }
    if let Some(at) = label.iter().position(|&b| b == b'|') {
        label.drain(..=at);
    }
    label
}

/// `FillSequenceFromBioseq` for the members that rpsbproc prints: the label, the title, the
/// length and the residues (or bases) of a query.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:952-983
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:248-253
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:297-319
// `Get1LetterSeqData` converts the residues, then asks for `inst.length` whether it is set or
// not: `GetLength` of an unset member throws `CUnassignedMember`, which ends rpsbproc as NCBI's
// does (exit 3, after the `SESSION` line). The residues serve the site output (PS16b; a
// nucleotide query's are translated, PS16c). A
// query without a title descriptor but with `annot` takes its title from a feature (refused),
// and a `CFastaReader` DefLine that is not a string is refused (`Query::unread_descriptor`).
// `GetTaxId` (biodata_core.cpp:977) reads the `org` and `source` descriptors but never throws,
// and the taxid is not printed, so those descriptors are only read past.
fn fill_sequence(query: &Query) -> Result<(Vec<u8>, Vec<u8>, u32, Vec<u8>), ProcessError> {
    // `SeqDataTo1LtrSeq` (biodata_core.cpp:149-197) converts the encodings with
    // `CSeqportUtil::Convert`, which LOSAT ports for every encoding but NCBI8aa and NCBIstdaa
    // codes from 28 (refused, before the length, as the conversion comes first; NCBIpaa and
    // NCBIpna throw NCBI's exception), for the `seq-data` and for each literal of `ext delta`
    // (`Get1LetterSeqData` converts each, in order).
    let residues = one_letter_residues(query).map_err(|(what, name)| {
        if name == UNCONVERTIBLE {
            return ProcessError::Std(UNCONVERTIBLE.to_string());
        }
        ProcessError::Unsupported(format!(
            "a query whose {what} is {name} is not supported by LOSAT's rpsbproc"
        ))
    })?;
    let Some(length) = query.length else {
        return Err(unassigned_member("NCBI-Sequence", "Seq-inst", "length"));
    };
    if let Some(what) = query.unread_descriptor {
        let local = query
            .ids
            .iter()
            .any(|id| matches!(id, ParsedSeqId::Local(_)));
        // `ObtainUserDeflineFromBioseq` runs for a query with a local id only.
        if local {
            return Err(ProcessError::Unsupported(format!(
                "a query Bioseq with {what} is not supported by LOSAT's rpsbproc"
            )));
        }
    }
    let title = match &query.title {
        Some(title) => title.clone(),
        None if query.has_annot => {
            return Err(ProcessError::Unsupported(
                "a query Bioseq without a title descriptor but with annot (NCBI takes the title from a feature) is not supported by LOSAT's rpsbproc".to_string(),
            ));
        }
        None => Vec::new(),
    };
    let mut residues = residues;
    // `Get1LetterSeqData` cuts the letters to `inst.length` (`GetLength`).
    residues.truncate(length as usize);
    Ok((query_accession(query), title, length, residues))
}

/// The kind of a `Seq-data` after `SeqDataTo1LtrSeq`: 0 for the protein encodings, 1 for the
/// nucleotide ones.
type LetterType = u8;

/// `SeqDataTo1LtrSeq`: the letters and their kind, or the name of an encoding that LOSAT does
/// not convert. (`gap` data, type -1, is refused as a `seq-data`, as in PS16a: the reader does
/// not check its grammar; the archive reader leaves the literals of `gap` data out.)
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:149-197
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:367-374,651-658,777-784
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1504-1507
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:886-893,950-958
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_tables.cpp:1822-1826
// The tables `C2naToIupacna` and `C4naToIupacna` (sequtil_tables.cpp:498, 1094) spell the
// bases of each byte from its high bits: NCBI2na codes 0 to 3 are `ACGT`, NCBI4na codes 0 to
// 15 are `NACMGRSVTWYHKDBN` (code 0, a gap, is `N`); NCBI8na codes 16 and up are the byte
// 0xff. A whole number of bytes is converted (`AdjustLength` keeps every base of the data), and
// `Get1LetterSeqData` cuts the result to `inst.length`. iupacaa to ncbieaa copies the letters;
// ncbistdaa codes 0 to 27 are the characters `-ABCDEFGHIKLMNPQRSTVWXYZU*OJ` and the other
// codes, 0xff, are refused as in PS16a, and so is NCBI8aa (PS16a). NCBIpaa and NCBIpna have no
// conversion in `CSeqportUtil`: NCBI's `runtime_error` (`UNCONVERTIBLE`).
fn one_letter(data: &(&'static str, Vec<u8>)) -> Result<(LetterType, Vec<u8>), &'static str> {
    const STDAA_TO_EAA: &[u8; 28] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
    const NA2: &[u8; 4] = b"ACGT";
    const NA4: &[u8; 16] = b"NACMGRSVTWYHKDBN";
    let (name, bytes) = data;
    Ok(match *name {
        "ncbieaa" | "iupacaa" => (0, bytes.clone()),
        "ncbistdaa" => {
            if bytes.iter().any(|&code| code >= 28) {
                return Err(name);
            }
            (
                0,
                bytes
                    .iter()
                    .map(|&code| STDAA_TO_EAA[code as usize])
                    .collect(),
            )
        }
        "iupacna" => (1, bytes.clone()),
        "ncbi2na" => (
            1,
            bytes
                .iter()
                .flat_map(|&b| [6u8, 4, 2, 0].map(|shift| NA2[usize::from((b >> shift) & 3)]))
                .collect(),
        ),
        "ncbi4na" => (
            1,
            bytes
                .iter()
                .flat_map(|&b| [NA4[usize::from(b >> 4)], NA4[usize::from(b & 15)]])
                .collect(),
        ),
        "ncbi8na" => (
            1,
            bytes
                .iter()
                .map(|&b| if b < 16 { NA4[usize::from(b)] } else { 0xff })
                .collect(),
        ),
        "ncbipna" | "ncbipaa" => return Err(UNCONVERTIBLE),
        _ => return Err(name),
    })
}

/// `CSeqportUtil`'s exception for an encoding without a conversion (NCBIpna, NCBIpaa): a
/// `std::runtime_error`, reported by `CNcbiApplication` (oracle `pna` of PS16c: `Error: (106.7)
/// Application's execution failed: Requested conversion not implemented`, exit 3).
// NCBI reference (598d8ae6): c++/src/objects/seq/seqport_util.cpp:2029-2035
// NCBI reference (598d8ae6): c++/src/objects/seq/seqport_util.cpp:1947-1959
const UNCONVERTIBLE: &str = "Requested conversion not implemented";

/// `Get1LetterSeqData` before its cut to `inst.length`: the `seq-data`, else the literals of
/// `ext delta` of the kind of the first one with letters, one after another; or what to refuse
/// (the place and the encoding).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:199-246
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/biodata_core.cpp:248-253
fn one_letter_residues(query: &Query) -> Result<Vec<u8>, (&'static str, &'static str)> {
    if let Some(data) = &query.seq_data {
        return Ok(one_letter(data).map_err(|name| ("Seq-data", name))?.1);
    }
    let mut residues = Vec::new();
    let mut literal_type = None;
    for data in &query.delta_literals {
        let (kind, letters) =
            one_letter(data).map_err(|name| ("Seq-data of a delta literal", name))?;
        if *literal_type.get_or_insert(kind) == kind {
            residues.extend_from_slice(&letters);
        }
    }
    Ok(residues)
}

/// `x_GetPssmId`: the PSSM-ID of a subject id (`ProcessCDAlign` keeps the hits whose PSSM-ID
/// is not 0).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_biodata.cpp:175-191
// (The `int` id is returned as `PssmId_t`, an `unsigned int`.)
fn pssm_id(id: &ParsedSeqId) -> Result<u32, ProcessError> {
    let ParsedSeqId::General(dbtag) = id else {
        return Err(choice_selection(
            "NCBI-Seqloc",
            "Seq-id",
            "general",
            id.alternative_name(),
        ));
    };
    Ok(match &dbtag.tag {
        ObjectId::Id(value) => *value as u32,
        ObjectId::Str(text) => compute_pseudo_pssm_id(text),
    })
}

/// `ProcessCDQuery`: the query's hits (`ProcessCDAlign`) and their classification, of all of
/// them for a protein query, of each reading frame's for a nucleotide query.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_biodata.cpp:195-244
// (`ParseAlignSegs` is `domains::parse_align` and `domains::parse_std_align`; it takes the
// hit's id from the second id of a Dense-seg or of the first Std-seg, and reads past `ids` when
// there are fewer: refused. The reading frames of a nucleotide query's hits are the classes of
// `domains::sort_reading_frames`; its exception for a query of length 0 is NCBI's, after the
// hits are read.)
fn process_cd_query(
    aligns: &[Align],
    query: &mut DomQuery,
    data: &CddData,
) -> Result<(), ProcessError> {
    for align in aligns {
        let (dom, hit) = match &align.segs {
            Segs::Denseg(seg) => {
                // (`process_blast_results` refuses a Dense-seg without ids before this.)
                if seg.ids.len() < 2 || seg.starts.len() < 2 * seg.lens.len() {
                    return Err(ProcessError::Unsupported(
                        "a Dense-seg with fewer than 2 ids or fewer starts than 2 per segment (NCBI reads past them) is not supported by LOSAT's rpsbproc".to_string(),
                    ));
                }
                (None, &seg.ids[1])
            }
            Segs::Std(std) => {
                // (`process_blast_results` finds a first id in the segments before this, so
                // there is at least one.)
                let Some(hit) = &std.hit_id else {
                    return Err(ProcessError::Unsupported(
                        "a Seq-align whose first Std-seg has fewer than 2 ids (NCBI reads past them) is not supported by LOSAT's rpsbproc".to_string(),
                    ));
                };
                let dom =
                    parse_std_align(&std.segs, align, 0).map_err(ProcessError::Unsupported)?;
                (Some(dom), hit)
            }
            Segs::Other(_) => unreachable!("the caller refuses other segments"),
        };
        let pssm_id = pssm_id(hit)?;
        if pssm_id == 0 {
            continue;
        }
        if let Some(name) = align.score_of_other_type {
            return Err(ProcessError::Unsupported(format!(
                "a Seq-align whose {name} score has a value of another type is not supported by LOSAT's rpsbproc"
            )));
        }
        // (A NaN E-value or bit score, `NOT-A-NUMBER` in the text, is kept: `domains::calculate`
        // sorts it after every number; the DOMAINS line prints it `-nan`, `ncbi_str::NCBI_NAN`.)
        let dom = match (dom, &align.segs) {
            (Some(mut dom), _) => {
                dom.pssm_id = pssm_id;
                dom
            }
            (None, Segs::Denseg(seg)) => parse_align(seg, align, pssm_id),
            (None, _) => unreachable!("a Dense-seg"),
        };
        // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:1809,1836
        // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_base.cpp:205-206
        // A kept hit whose range starts or ends at 2^31 - 1, of a protein or a nucleotide query:
        // the `int` sums of its DOMAINS and SUPERFAMILIES lines overflow (undefined; PS16c review
        // oracle `plusmax`: both 0.5.1 builds print -2147483648), and so may the overlap sets of
        // `Calculate` that start from the range (`GetTotalResidues`' `to - from + 1`, the `+ 1`
        // of `Cross` and `Merge`): computed as their values (`printed_position`, `segset::Seg`;
        // PD-LOSAT-NCBI-DEFECTS exception 12).
        query.aligns.push(dom);
    }
    if !query.is_na {
        let all: Vec<usize> = (0..query.aligns.len()).collect();
        query.frames[0] = calculate(&mut query.aligns, &all, data);
        return Ok(());
    }
    let frames = sort_reading_frames(&query.aligns, query.length)
        .map_err(|what| ProcessError::Std(what.to_string()))?;
    for (frame, selected) in frames.iter().enumerate() {
        if !selected.is_empty() {
            query.frames[frame] = calculate(&mut query.aligns, selected, data);
        }
    }
    Ok(())
}

/// `CSeq_align::GetSeq_id` for index 0: the id of the query row of an alignment, or what LOSAT refuses
/// (NCBI's exceptions when no segment names one, and the location variants whose id LOSAT does
/// not look up).
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:317-358
// NCBI reference (598d8ae6): c++/src/objects/seqalign/Seq_align.cpp:390-398
// (A `null` location has no id; `empty`, `whole`, `int` and `pnt` have theirs. The exceptions
// end rpsbproc as its other `CException`s of a session do: PS16c oracles `noids2` and
// `packed`, the two `Error:` lines, exit 3. `GetSeq_id` is a `const` method, so the report names
// no function, as `UNASSIGNED_MEMBER`'s. The `packed` segments are read for their grammar only:
// NCBI's `ParseAlignSegs` has a branch for them, which this exception keeps it from reaching.
// `dendiag`, `disc`, `spliced` and `sparse` segments, which `GetSeq_id` reads, are refused.)
fn query_row_id(align: &Align) -> Result<&ParsedSeqId, ProcessError> {
    let seqalign_exception = |line: u32, kind: &'static str, message: &str| {
        ProcessError::Exception(NcbiException::new(
            ThrowSite {
                file: "src/objects/seqalign/Seq_align.cpp",
                line,
                kind,
                function: None,
            },
            message.to_string(),
        ))
    };
    match &align.segs {
        Segs::Denseg(seg) => seg.ids.first().ok_or_else(|| {
            ProcessError::Unsupported(
                "a Dense-seg without ids (NCBI's CSeq_align::GetSeq_id throws) is not supported by LOSAT's rpsbproc".to_string(),
            )
        }),
        Segs::Std(std) => match &std.query_id {
            QueryRowId::Id(id) => Ok(id),
            QueryRowId::Location(name) => Err(ProcessError::Unsupported(format!(
                "a Std-seg without ids whose first location is a {name} is not supported by LOSAT's rpsbproc"
            ))),
            QueryRowId::Missing => Err(seqalign_exception(
                396,
                "CSeqalignException::eInvalidRowNumber",
                "CSeq_align::GetSeq_id(): can not get seq-id for the row requested.",
            )),
        },
        Segs::Other("packed") => Err(seqalign_exception(
            391,
            "CSeqalignException::eUnsupported",
            "CSeq_align::GetSeq_id() currently does not handle this type of alignment.",
        )),
        Segs::Other(name) => Err(ProcessError::Unsupported(format!(
            "a Seq-align with {name} segments is not supported by LOSAT's rpsbproc"
        ))),
    }
}

/// `ProcessBlastResults`: each query with the alignments that follow in the list while their
/// first id matches the query's local id (else its first id).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/cdalignproc_biodata.cpp:308-366
// (`CSeq_id::Match` is true when `Compare` gives `e_YES`, `seqid::defline_set::compare`. The masks of
// each query are not printed. `query_row_id` is `GetSeq_id` for index 0 of a `Dense-seg`, of the
// `Std-seg`s of a translated search and of a `Packed-seg` (NCBI's exception); the other
// segments are refused there.)
pub(crate) fn process_blast_results(
    entries: &[SeqEntry],
    aligns: &[Align],
    data: &CddData,
) -> Result<Vec<DomQuery>, ProcessError> {
    let mut queries = Vec::with_capacity(entries.len());
    let mut next = 0usize;
    for entry in entries {
        let SeqEntry::Seq(query) = entry else {
            return Err(choice_selection("NCBI-Seqset", "Seq-entry", "seq", "set"));
        };
        let (accession, title, length, residues) = fill_sequence(query)?;
        let mut dom_query = DomQuery {
            accession,
            title,
            length,
            is_na: query.mol != 3,
            residues,
            aligns: Vec::new(),
            frames: Default::default(),
            translated: None,
        };
        // NCBI reference (598d8ae6): c++/src/objects/seq/Bioseq.cpp:271-279,327-343
        let current = crate::seqid::rank::find_best(&query.ids, |id| {
            if matches!(id, ParsedSeqId::Local(_)) {
                crate::seqid::rank::best_rank_score(id)
            } else {
                i32::MAX
            }
        })
        .or(query.ids.first());
        let start = next;
        while next < aligns.len() {
            let aligned = query_row_id(&aligns[next])?;
            let Some(current) = current else {
                return Err(ProcessError::Unsupported(
                    "a query Bioseq without ids (NCBI dereferences its missing id) is not supported by LOSAT's rpsbproc".to_string(),
                ));
            };
            if compare(current, aligned) != Sic::Yes {
                break;
            }
            next += 1;
        }
        if next > start {
            process_cd_query(&aligns[start..next], &mut dom_query, data)?;
        }
        queries.push(dom_query);
    }
    Ok(queries)
}

#[cfg(test)]
#[path = "query_tests.rs"]
mod tests;
