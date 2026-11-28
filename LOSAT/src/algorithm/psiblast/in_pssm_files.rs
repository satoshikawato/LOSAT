//! psiblast `-in_pssm` with rounds or with the PSSM file options (series PS, item PS5b): the input
//! PSSM as NCBI's `CPssmWithParameters` holds it when round 1 starts (after `CPsiBlastImpl`'s
//! constructor, which computes the scores of a PSSM of frequency ratios only), which
//! `SavePssmToFile` writes then; the same object after round 1's K/A write-back, which it writes
//! again after a round 1 without alignments with `-save_pssm_after_last_round`; the ASCII PSSM of a
//! file PSSM (`ExtractPssmAncillaryData`: the PSI blocks only, the short form or, with
//! `informationContent`, the long one); the Seq-id that the PSSMs of the later rounds keep; and
//! LOSAT's refusals of what NCBI cannot write.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:385-408
//! With `-in_pssm` the loop's first `SavePssmToFile` of `pssm` writes the input PSSM, the object
//! that `CPsiBlast` searches (the `CPsiBlast` made from `pssm`, `db_adapter` and `psi_opts` holds
//! the same `CRef`); its ASCII PSSM takes the ancillary data of the PSSM itself, which stays for
//! the files written after the loop until the next PSSM is computed (`ComputePssmForNextIteration`
//! replaces `m_AncillaryData`).

use anyhow::{bail, Result};

use crate::core::blast_util::{BLASTAA_SIZE, BLAST_SCORE_MIN};
use crate::core::psi::KarlinBlk;
use crate::report::serial::pssm::{PssmWithParameters, SeqData, SeqEntry};

use super::pssm::{AsciiPssmDiagnostics, PssmRequest};
#[cfg(test)]
use crate::algorithm::pssm_search::in_pssm::with_computed_scores;
use crate::algorithm::pssm_search::in_pssm::InPssmFile;
use crate::report::serial::encoder::ReplacedChar;

/// The query PSSM of psiblast's rounds with `-in_pssm` (`rounds::PsiRoundPlan::input`): when a
/// PSSM file is written, the object that `SavePssmToFile` writes. (The Seq-id that the PSSM of
/// each later round keeps is `PsiSearch::query_seq_loc_id`.)
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct InPssmRounds {
    object: Option<PssmWithParameters>,
    /// LOSAT's refusal of a file NCBI would write back in a way LOSAT does not (`rounds_input`),
    /// made once the database is opened (`refuse_written_back`).
    refusal: Option<String>,
}

/// The input PSSM of a run of psiblast's rounds or PSSM files with `-in_pssm` (`request`: the
/// files), with LOSAT's refusal of a file NCBI would write back in a way LOSAT does not, which the
/// run makes once the database is opened, before any output (`InPssmRounds::refuse_written_back`:
/// NCBI's error of a database it cannot open comes first, and so do the errors that stop NCBI
/// before round 1, `PsiSearch::stop_before_search`):
/// - with `-out_pssm`, a PSSM of frequency ratios only without `rpsdbparams` (NCBI writes part of
///   the file, then stops with a C++ exception: "cannot write unassigned member matrixName",
///   rc 255); every member that the file holds is written back (`report::serial::pssm`);
/// - with `-out_ascii_pssm`, a PSSM with `informationContent` (the long form) whose
///   `weightedResFreqsPerPos`, `informationContent`, `numMatchingSeqs` or the
///   `gaplessColumnWeights` it reads do not cover every column (`PrintAsciiPssm` dereferences
///   the missing matrix or reads past the vectors: undefined behaviour).
///
/// A PSSM that `CPsiBlastValidate::Pssm` refuses gets NCBI's error when the search starts, before
/// any file (the batch stop, `blast_engine::psiblast_stop_at_pssm_query`).
pub(crate) fn rounds_input(file: InPssmFile, request: &PssmRequest) -> InPssmRounds {
    if !file.valid {
        return InPssmRounds {
            object: None,
            refusal: None,
        };
    }
    let refusal = written_back_refusal(&file, request)
        .err()
        .map(|error| error.to_string());
    let files = request.asn_path.is_some() || request.ascii_path.is_some();
    InPssmRounds {
        object: files.then_some(file.object),
        refusal,
    }
}

/// The refusal of a PSSM whose matrix name NCBI cannot write (`InPssmFile::matrix_name_unwritable`),
/// where `writer` names what writes it: `with -out_pssm` here, and the archive of `-outfmt 11`
/// (`algorithm::psiblast::in_pssm::check_in_pssm_options`).
pub(crate) fn matrix_name_refusal(writer: &str) -> String {
    format!(
        "a PSSM of frequency ratios only without params.rpsdbparams {writer} (NCBI BLAST+ cannot write the matrix name it lacks and stops with a C++ exception) is not supported by LOSAT's PSIBLAST"
    )
}

/// The refusals of `rounds_input`.
fn written_back_refusal(file: &InPssmFile, request: &PssmRequest) -> Result<()> {
    if request.asn_path.is_some() && file.matrix_name_unwritable() {
        bail!("{}", matrix_name_refusal("with -out_pssm"));
    }
    if request.ascii_path.is_some() {
        check_long_ascii(&file.object)?;
        if !ascii_letters_defined(&file.object) {
            bail!(
                "an ASCII PSSM (-out_ascii_pssm) of a PSSM query with a residue that NCBIstdaa does not have (NCBI BLAST+ prints its letter from memory it never set) is not supported by LOSAT's PSIBLAST"
            );
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:252-261
// NCBI reference (598d8ae6): c++/src/objects/scoremat/Pssm.cpp:88-93 (`GetQuerySequenceData`)
// NCBI reference (598d8ae6): c++/src/objects/seq/seqport_util.cpp:2045-2054
// NCBI reference (598d8ae6): c++/src/util/sequtil/sequtil_convert_imp.cpp:898-905 and the tables
// of c++/src/util/sequtil/sequtil_tables.cpp:1674-1707 (`CIupacaaToStdaa`: the code of each
// capital for it and for its lower case, 0xff for every other byte) and 1748-1781 (`CEaaToStdaa`:
// the capitals, `*` and `-`, 0xff for every other byte)
/// Whether NCBI's ASCII PSSM prints a letter of NCBIstdaa for every query residue of `object`:
/// the query's data converts through `CSeqConvert`, whose tables give a code to the letters of
/// iupacaa in either case (a lower-case letter prints its capital: M4a-1c LG probe
/// `tk_low_iupac.files`), and to the capitals, `*` and `-` of ncbieaa; any other byte gives 0xff,
/// whose `char` (-1) reads before `NCBISTDAA_TO_AMINOACID`, and an NCBIstdaa code of 28 or more
/// reads past it (NCBI BLAST+ 2.17.0 printed a NUL byte, M4a-1c E probe
/// `old_ps_bad_iupacaa.ps.files`).
fn ascii_letters_defined(object: &PssmWithParameters) -> bool {
    match object
        .pssm
        .query
        .as_ref()
        .and_then(SeqEntry::seq)
        .and_then(|query| query.inst.seq_data.as_ref())
    {
        Some(SeqData::Ncbistdaa(codes)) => {
            codes.iter().all(|&code| usize::from(code) < BLASTAA_SIZE)
        }
        Some(SeqData::Iupacaa(text)) => text.iter().all(u8::is_ascii_alphabetic),
        Some(SeqData::Ncbieaa(text)) => text
            .iter()
            .all(|&c| c.is_ascii_uppercase() || c == b'*' || c == b'-'),
        Some(
            SeqData::Iupacna(_) | SeqData::Ncbi2na(_) | SeqData::Ncbi4na(_) | SeqData::Other(_),
        )
        | None => true,
    }
}

/// The refusal of the long ASCII form of a PSSM file that NCBI cannot print (`rounds_input`).
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:186-188, 278-301 (quoted
// at `pssm::print_ascii_pssm`: the long form reads `weighted_res_freqs` at row `r` and column `i`,
// `info_content[i]`, `num_matching_seqs[i]` and, where that is above 1 and the letter is not X,
// `gapless_col_weights[i]`; `GetWeightedResidueFrequencies` gives NULL without the member)
fn check_long_ascii(object: &PssmWithParameters) -> Result<()> {
    let pssm = &object.pssm;
    let Some(data) = pssm.intermediate_data.as_ref() else {
        return Ok(());
    };
    let info = data.information_content.as_deref().unwrap_or(&[]);
    if info.is_empty() {
        return Ok(());
    }
    let columns = usize::try_from(pssm.num_columns).unwrap_or(0);
    let rows_n = usize::try_from(pssm.num_rows).unwrap_or(0);
    let cells = columns * rows_n;
    let weighted = data.weighted_res_freqs_per_pos.as_deref().unwrap_or(&[]);
    let matching = data.num_matching_seqs.as_deref().unwrap_or(&[]);
    let gapless = data.gapless_column_weights.as_deref().unwrap_or(&[]);
    let letters = query_codes(object);
    let gapless_short = (gapless.len()..columns.min(matching.len())).any(|i| {
        matching[i] > 1 && letters.get(i).copied() != Some(crate::utils::matrix::ncbistdaa::X)
    });
    // `Convert2Matrix` reads the first numRows × numColumns values (a longer list's others are
    // never read; `_ASSERT` only): fewer is a read past the list.
    if weighted.len() < cells || info.len() < columns || matching.len() < columns || gapless_short {
        bail!(
            "an ASCII PSSM (-out_ascii_pssm) of a PSSM file whose informationContent comes without weightedResFreqsPerPos, informationContent, numMatchingSeqs or gaplessColumnWeights for every column (NCBI BLAST+ reads past them) is not supported by LOSAT's PSIBLAST"
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:274-282
    // A weighted frequency whose rounded percentage is not an `int` (or is NaN) is undefined
    // behaviour in that conversion (x86-64 prints -2147483648): LOSAT's refusal, for the cells
    // NCBI prints (the engine's frequencies lie in [0, 1]; only a file reaches this).
    let by_row = pssm.by_row == Some(true);
    let scores = crate::algorithm::pssm_search::in_pssm::convert_2_matrix(
        pssm.final_data
            .as_ref()
            .map_or(&[][..], |data| data.scores.as_slice()),
        rows_n,
        columns,
        by_row,
        BLAST_SCORE_MIN,
    );
    let weighted = crate::algorithm::pssm_search::in_pssm::convert_2_matrix(
        weighted, rows_n, columns, by_row, 0.0,
    );
    let outside_int = scores.iter().zip(&weighted).any(|(scores, weighted)| {
        super::pssm::RESIDUE_ORDER.iter().any(|&r| {
            let value = 100.0 * weighted[r];
            let rounded = value + if value > 0.0 { 0.5 } else { -0.5 };
            scores[r] != BLAST_SCORE_MIN
                && !(rounded > -2_147_483_649.0 && rounded < 2_147_483_648.0)
        })
    });
    if outside_int {
        bail!(
            "an ASCII PSSM (-out_ascii_pssm) of a PSSM file with a weighted residue frequency whose percentage is not an int or is NaN (NCBI BLAST+'s conversion is undefined behaviour) is not supported by LOSAT's PSIBLAST"
        );
    }
    Ok(())
}

/// The query's letters as NCBIstdaa codes (`CPssm::GetQuerySequenceData`, which converts an
/// IUPACaa or NCBIeaa `Seq-data`, an IUPACaa letter in either case: `ascii_letters_defined`;
/// `in_pssm` admits only letters of NCBIstdaa).
fn query_codes(object: &PssmWithParameters) -> Vec<u8> {
    const LETTERS: &[u8; BLASTAA_SIZE] =
        crate::algorithm::pssm_search::in_pssm::NCBISTDAA_TO_AMINOACID;
    let query = object.pssm.query.as_ref().and_then(SeqEntry::seq);
    let mut codes = match query.and_then(|query| query.inst.seq_data.as_ref()) {
        Some(SeqData::Ncbistdaa(codes)) => codes.clone(),
        Some(SeqData::Iupacaa(text)) => text
            .iter()
            .map(|letter| {
                let letter = letter.to_ascii_uppercase();
                LETTERS.iter().position(|&c| c == letter).unwrap_or(0) as u8
            })
            .collect(),
        Some(SeqData::Ncbieaa(text)) => text
            .iter()
            .map(|letter| LETTERS.iter().position(|c| c == letter).unwrap_or(0) as u8)
            .collect(),
        Some(
            SeqData::Iupacna(_) | SeqData::Ncbi2na(_) | SeqData::Ncbi4na(_) | SeqData::Other(_),
        )
        | None => Vec::new(),
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:202, 253-258
    // NCBI reference (598d8ae6): c++/src/objects/scoremat/Pssm.cpp:61-74 (`CPssm::GetQueryLength`)
    // `PrintAsciiPssm` loops up to `GetQueryLength()`: 0 when `Seq-inst.length` is absent, and
    // `length` when the query has one letter more than its Seq-inst's `length` (searched as its
    // first `length` letters, `in_pssm`; NCBI BLAST+ 2.17.0, line PS item PS7c).
    match query.and_then(|query| query.inst.length) {
        None => codes.clear(),
        Some(length)
            if usize::try_from(length)
                .is_ok_and(|length| length >= 1 && length + 1 == codes.len()) =>
        {
            codes.pop();
        }
        Some(_) => {}
    }
    codes
}

/// The files of the input PSSM (`InPssmRounds::files`).
#[derive(Debug)]
pub(crate) struct InputFiles {
    /// Its ASN.1 text (the bytes before NCBI's throw, `thrown`).
    pub(crate) asn: Option<Vec<u8>>,
    /// Its ASCII PSSM (none after NCBI's throw: `SavePssmToFile` writes the ASN.1 file first).
    pub(crate) ascii: Option<Vec<u8>>,
    /// NCBI's posts of the bytes that the ASN.1 text's stream replaces.
    pub(crate) posts: Vec<ReplacedChar>,
    /// NCBI's error where its ASN.1 text stream threw at an element stated as `NULL` in the
    /// query (`pssm::WrittenPssm::thrown`), which ends the run after the file is written.
    pub(crate) thrown: Option<crate::cli::NativeError>,
}

impl InPssmRounds {
    /// LOSAT's refusal of the input PSSM's files (`rounds_input`), made once the database is opened
    /// for a run that searches (`blastp::blast_engine` `search_cli`).
    pub(crate) fn refuse_written_back(&self) -> Result<()> {
        match &self.refusal {
            Some(refusal) => bail!("{refusal}"),
            None => Ok(()),
        }
    }

    /// The files of the input PSSM at the start of round 1 (`SavePssmToFile`): its ASN.1 text and
    /// its ASCII PSSM, as requested.
    ///
    /// With them NCBI's post of each byte that its ASN.1 text stream replaces in the file (a
    /// query title outside the visible ASCII range; `pssm::file_posts`), and NCBI's error where
    /// that stream throws (`InputFiles::thrown`: an element stated as `NULL` in the query, NCBI
    /// BLAST+ 2.17.0 probe D3c-3b `probe4`, `probe5`, `probe6`: the file up to the element, no
    /// ASCII PSSM).
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:392-408
    pub(crate) fn files(&self, request: &PssmRequest) -> Result<InputFiles> {
        let Some(object) = self.object.as_ref() else {
            return Ok(InputFiles {
                asn: None,
                ascii: None,
                posts: Vec::new(),
                thrown: None,
            });
        };
        let (asn, posts, thrown) = match request.asn_path {
            Some(_) => {
                let written = crate::report::serial::pssm::asn_text_written(object, None)
                    .map_err(|error| anyhow::anyhow!("{error}"))?;
                let thrown = written
                    .thrown
                    .map(|thrown| crate::cli::NativeError::thrown("psiblast", &thrown));
                (Some(written.bytes), written.replaced, thrown)
            }
            None => (None, Vec::new(), None),
        };
        let ascii = request
            .ascii_path
            .as_ref()
            .filter(|_| thrown.is_none())
            .map(|_| ascii_pssm(object));
        Ok(InputFiles {
            asn,
            ascii,
            posts,
            thrown,
        })
    }

    /// The ASN.1 text of the input PSSM after round 1's K/A write-back (`ungapped`: `kbp_psi[0]`,
    /// `gapped`: `kbp_gap_psi[0]` of the round's search), which `-save_pssm_after_last_round`
    /// writes after a round 1 without alignments.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:150-168
    /// NCBI reference (598d8ae6): c++/src/objects/scoremat/Pssm.cpp:157-191 (`CPssm::SetLambda`
    /// and the others set the members of `finalData`)
    /// The redo of the composition-based statistics works on copies of the score block
    /// (`s_BlastScoreBlk_Copy` in `Blast_RedoAlignmentCore_MT`, blast_kappa.c:3254-3259), so the
    /// blocks are those that `PsiBlastSetupScoreBlock` set for the search.
    pub(crate) fn written_back(
        &self,
        request: &PssmRequest,
        ungapped: KarlinBlk,
        gapped: KarlinBlk,
    ) -> Result<(Option<Vec<u8>>, Vec<ReplacedChar>)> {
        let (Some(object), Some(_)) = (self.object.as_ref(), request.asn_path.as_ref()) else {
            return Ok((None, Vec::new()));
        };
        // Round 1's write of the object (`files`) read each delay buffer of another format into
        // it (`PssmWithParameters::update_delayed`): this write writes the values read, without
        // the messages of that read again (line RG-F3, PR 2: probe AR4 P).
        let mut object = object.clone();
        object.update_delayed(crate::report::serial::asn_input::Encoding::Text);
        // The six values are written in place of the object's (NCBI sets them in its object).
        let karlin = crate::report::serial::pssm::FinalKarlin {
            lambda: gapped.lambda,
            kappa: gapped.k,
            h: gapped.h,
            lambda_ungapped: ungapped.lambda,
            kappa_ungapped: ungapped.k,
            h_ungapped: ungapped.h,
        };
        crate::report::serial::pssm::asn_text_posting_with_karlin(&object, Some(karlin))
            .map(|(asn, posts)| (Some(asn), posts))
            .map_err(|error| anyhow::anyhow!("{error}"))
    }
}

/// `PrintAsciiPssm` of a file PSSM, with the ancillary data that `ExtractPssmAncillaryData` takes
/// from it: the scores of the file (`GetScores`), the query's letters (`GetQuerySequenceData`), the
/// long form with the file's diagnostics when it has `informationContent`, and the PSI blocks of
/// its `finalData` (a member it lacks is `CPssm::kInvalidStat`, 0).
// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:792-804
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:117-134
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:186-196
fn ascii_pssm(object: &PssmWithParameters) -> Vec<u8> {
    let pssm = &object.pssm;
    let rows_n = usize::try_from(pssm.num_rows).unwrap_or(0);
    let columns = usize::try_from(pssm.num_columns).unwrap_or(0);
    let by_row = pssm.by_row == Some(true);
    let final_data = pssm.final_data.as_ref();
    let scores = crate::algorithm::pssm_search::in_pssm::convert_2_matrix(
        final_data.map_or(&[][..], |data| data.scores.as_slice()),
        rows_n,
        columns,
        by_row,
        BLAST_SCORE_MIN,
    );
    let letters = query_codes(object);
    let data = pssm.intermediate_data.as_ref();
    let info = data.and_then(|data| data.information_content.as_deref());
    let weighted;
    let long = match (data, info) {
        (Some(data), Some(info)) if !info.is_empty() => {
            weighted = crate::algorithm::pssm_search::in_pssm::convert_2_matrix(
                data.weighted_res_freqs_per_pos.as_deref().unwrap_or(&[]),
                rows_n,
                columns,
                by_row,
                0.0,
            );
            Some(AsciiPssmDiagnostics {
                weighted_residue_freqs: &weighted,
                information_content: info,
                num_matching_seqs: data.num_matching_seqs.as_deref().unwrap_or(&[]),
                gapless_column_weights: data.gapless_column_weights.as_deref().unwrap_or(&[]),
            })
        }
        _ => None,
    };
    let stat = |value: Option<f64>| value.unwrap_or(0.0);
    let psi_ungapped = KarlinBlk {
        lambda: stat(final_data.and_then(|data| data.lambda_ungapped)),
        k: stat(final_data.and_then(|data| data.kappa_ungapped)),
        h: stat(final_data.and_then(|data| data.h_ungapped)),
    };
    let psi_gapped = KarlinBlk {
        lambda: stat(final_data.map(|data| data.lambda)),
        k: stat(final_data.map(|data| data.kappa)),
        h: stat(final_data.map(|data| data.h)),
    };
    super::pssm::write_ascii_pssm(
        &scores,
        &letters,
        long,
        &[
            ("PSI Ungapped         ", psi_ungapped),
            ("PSI Gapped           ", psi_gapped),
        ],
    )
}

#[cfg(test)]
#[path = "in_pssm_files_tests.rs"]
mod tests;
