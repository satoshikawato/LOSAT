//! NCBI's Clustal scanner (`CAlnScannerClustal`) and the verification of the base scanner
//! (`CAlnScanner::xVerifyAlignmentData`): the ids and the rows of a Clustal file, or the failure
//! that NCBI's scanner throws (every throw of the scanner is a failed read, `super::read`).

use super::lines::PeekAheadStream;
use crate::ncbi_str::{string_to_int, truncate_spaces};

/// The alignment as the scanner exports it (`SAlignmentFile`): the id string and the row of each
/// sequence, in the order of the first block.
pub(crate) struct Scanned {
    pub(crate) ids: Vec<Vec<u8>>,
    pub(crate) rows: Vec<Vec<u8>>,
}

/// The scanner stopped (a `SShowStopper`, which `CAlnReader::Read` hands to the error reporter
/// that throws it).
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct ScanFailure;

struct BlockInfo {
    seq_count: usize,
    line_length: usize,
    first: bool,
}

/// `sTerminateBlock`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_clustal.cpp:59-85
/// ```c++
///     if (blockInfo.first) {
///         numSeqs = blockInfo.seqCount;
///         blockInfo.first = false;
///     }
///     else
///     if (numSeqs != blockInfo.seqCount) {
///         ...
///         throw SShowStopper(
///             lineCount,
///             EAlnSubcode::eAlnSubcode_BadSequenceCount,
///             description);
///     }
///     blockInfo.seqCount = 0;
///     blockInfo.lineLength = 0;
/// ```
fn terminate_block(num_seqs: &mut usize, block: &mut BlockInfo) -> Result<(), ScanFailure> {
    if block.first {
        *num_seqs = block.seq_count;
        block.first = false;
    } else if *num_seqs != block.seq_count {
        return Err(ScanFailure);
    }
    block.seq_count = 0;
    block.line_length = 0;
    Ok(())
}

/// `sIsConservationLine`: the result of `line.find_first_not_of(":.* \t\r\n")` taken as a `bool`,
/// so true unless that position is 0: a line that starts with `:`, `.` or `*` is a conservation
/// line also when other text follows, and a line of nothing else is one too (`npos` is not 0).
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_clustal.cpp:211-218
/// ```c++
/// bool
/// CAlnScannerClustal::sIsConservationLine(
///     const string& line)
/// {
///     return line.find_first_not_of(":.* \t\r\n");
/// };
/// ```
fn is_conservation_line(line: &[u8]) -> bool {
    line.iter().position(|c| !b":.* \t\r\n".contains(c)) != Some(0)
}

/// `NStr::Split(line, " \t", tokens, NStr::fSplit_Tokenize)`: the runs of other bytes, no empty
/// token.
fn tokenize(line: &[u8]) -> Vec<&[u8]> {
    line.split(|&c| c == b' ' || c == b'\t')
        .filter(|token| !token.is_empty())
        .collect()
}

/// `CAlnScannerClustal::xImportAlignmentData` and `sProcessClustalDataLine`, then the base
/// scanner's `xVerifyAlignmentData`.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_clustal.cpp:90-176
/// ```c++
///     while (iStr.ReadLine(line, lineCount)) {
///         if (lineCount == 1) {
///             if (NStr::StartsWith(line, "clustal", NStr::eNocase)) {
///                 continue;
///             }
///         }
///
///         NStr::TruncateSpacesInPlace(line);
///         if (line.empty()) {
///             if (inBlock) {
///                 sTerminateBlock(lineCount, numSeqs, blockInfo);
///                 inBlock = false;
///             }
///             continue;
///         }
///
///         if (sIsConservationLine(line)) {
///             if (!inBlock) {
///                 ...
///                 throw SShowStopper(...);
///             }
///             sTerminateBlock(lineCount, numSeqs, blockInfo);
///             inBlock = false;
///             continue;
///         }
///         // Search for end of block
///         vector<string> tokens;
///         NStr::Split(line, " \t", tokens, NStr::fSplit_Tokenize);
///         const auto num_tokens = tokens.size();
///         if (num_tokens < 2 || num_tokens > 3) {
///             ...
///             throw SShowStopper(...);
///         }
///
///         int seqLength = 0;
///         if (num_tokens == 3) {
///             seqLength = NStr::StringToInt(tokens[2], NStr::fConvErr_NoThrow);
///             if (!seqLength) {
///                 throw SShowStopper(...);
///             }
///         }
///
///         if (!inBlock) {
///             inBlock = true;
///         }
///
///         sProcessClustalDataLine(
///              tokens, lineCount,
///              blockInfo.seqCount,
///              numSeqs,
///              blockInfo.first,
///              blockInfo.lineLength);
///
///         mSequences[blockInfo.seqCount].push_back({tokens[1], lineCount});
///         ++(blockInfo.seqCount);
///     }
///
///     if (inBlock) {
///         ...
///         throw SShowStopper(...);
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner_clustal.cpp:220-307
/// ```c++
///     auto seqId = tokens[0];
///     if (firstBlock) {
///         TLineInfo existingInfo;
///         auto idComparison = xGetExistingSeqIdInfo(seqId, existingInfo);
///         if (idComparison != ESeqIdComparison::eDifferentChars) {
///             ...
///             throw SShowStopper(...);
///         }
///         mSeqIds.push_back({seqId, lineNum});
///         mSequences.push_back(vector<TLineInfo>());
///     }
///     else {
///         if (seqCount >= numSeqs) {
///             ...
///             throw SShowStopper(...);
///         }
///
///         if (seqId != mSeqIds[seqCount].mData) {
///             ... (every branch throws)
///         }
///     }
///
///     if (seqCount==0) {
///         blockLineLength = tokens[1].size();
///         return;
///     }
///
///     auto currentLineLength = tokens[1].size();
///     if (currentLineLength != blockLineLength) {
///         ...
///         throw SShowStopper(...);
///     }
/// ```
/// The base scanner (`xGetExistingSeqIdInfo`: the same id, or one that differs by case only, is
/// a duplicate):
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner.cpp:196-221
/// ```c++
///     for (const auto& idInfo : mSeqIds) {
///         if (seqId == idInfo.mData) {
///             existingInfo = idInfo;
///             return ESeqIdComparison::eIdentical;
///         }
///     }
///
///     string seqIdLower(seqId);
///     NStr::ToLower(seqIdLower);
///     for (const auto& idInfo : mSeqIds) {
///         string idLower(idInfo.mData);
///         NStr::ToLower(idLower);
///         if (seqIdLower == idLower) {
///             existingInfo = idInfo;
///             return ESeqIdComparison::eDifferByCase;
///         }
///     }
///
///     return ESeqIdComparison::eDifferentChars;
/// ```
pub(crate) fn scan(input: &mut PeekAheadStream<'_>) -> Result<Scanned, ScanFailure> {
    let mut ids: Vec<Vec<u8>> = Vec::new();
    let mut sequences: Vec<Vec<&[u8]>> = Vec::new();
    let mut in_block = false;
    let mut num_seqs = 0usize;
    let mut block = BlockInfo {
        seq_count: 0,
        line_length: 0,
        first: true,
    };
    while let Some((raw, line_count)) = input.read_line() {
        if line_count == 1 && raw.len() >= 7 && raw[..7].eq_ignore_ascii_case(b"clustal") {
            continue;
        }
        let line = truncate_spaces(raw);
        if line.is_empty() {
            if in_block {
                terminate_block(&mut num_seqs, &mut block)?;
                in_block = false;
            }
            continue;
        }
        if is_conservation_line(line) {
            if !in_block {
                return Err(ScanFailure);
            }
            terminate_block(&mut num_seqs, &mut block)?;
            in_block = false;
            continue;
        }
        let tokens = tokenize(line);
        if tokens.len() < 2 || tokens.len() > 3 {
            return Err(ScanFailure);
        }
        if tokens.len() == 3 {
            // `StringToInt(..., fConvErr_NoThrow)`: 0 for a token that is no number.
            if string_to_int(tokens[2], 0).unwrap_or(0) == 0 {
                return Err(ScanFailure);
            }
        }
        in_block = true;
        let id = tokens[0];
        if block.first {
            if ids.iter().any(|existing| existing.eq_ignore_ascii_case(id)) {
                return Err(ScanFailure);
            }
            ids.push(id.to_vec());
            sequences.push(Vec::new());
        } else {
            if block.seq_count >= num_seqs {
                return Err(ScanFailure);
            }
            if id != ids[block.seq_count].as_slice() {
                return Err(ScanFailure);
            }
        }
        if block.seq_count == 0 {
            block.line_length = tokens[1].len();
        } else if tokens[1].len() != block.line_length {
            return Err(ScanFailure);
        }
        sequences[block.seq_count].push(tokens[1]);
        block.seq_count += 1;
    }
    if in_block {
        return Err(ScanFailure);
    }
    let rows: Vec<Vec<u8>> = sequences.iter().map(|parts| parts.concat()).collect();
    verify_rows(&rows)?;
    Ok(Scanned { ids, rows })
}

/// `CAlnScanner::xVerifySingleSequenceData` with the sequence info of `SetClustal(eAlpha_Protein)`
/// (the alphabet of `GetAlphabetLetters(eAlpha_Protein)`, the gap `-` at every place, no
/// missing or match characters): a row may hold the letters of the alphabet and `-` only. NCBI's
/// states (the head allows gaps, the first other byte must be in the alphabet; the body allows the
/// alphabet and gaps; the tail would allow gaps only) end in an error at the first byte that is
/// neither, for the tail starts at one.
///
/// NCBI reference (598d8ae6): c++/src/objtools/readers/aln_scanner.cpp:71-130
/// ```c++
///     const string& alphabet = sequenceInfo.Alphabet();
///     const string legalInHead =
///         sequenceInfo.BeginningGap() + sequenceInfo.Missing();
///     const string legalInBody =
///         alphabet
///         + sequenceInfo.MiddleGap()
///         + sequenceInfo.Missing()
///         + sequenceInfo.Match();
///     const string legalInTail =
///         sequenceInfo.EndGap() + sequenceInfo.Missing();
///     ...
///         if (seqPart == ESeqPart::HEAD) {
///             auto startBody = seqData.find_first_not_of(legalInHead);
///             if (startBody == string::npos) {
///                 continue;
///             }
///             seqPart = ESeqPart::BODY;
///             seqData = seqData.substr(startBody);
///             if (alphabet.find(seqData[0]) == string::npos) {
///                 ... throw SShowStopper(...);
///             }
///         }
///         if (seqPart == ESeqPart::BODY) {
///             auto startTail = seqData.find_first_not_of(legalInBody);
///             if (startTail == string::npos) {
///                 continue;
///             }
///             seqPart = ESeqPart::TAIL;
///             seqData = seqData.substr(startTail);
///         }
///         if (seqPart == ESeqPart::TAIL) {
///             auto startBad = seqData.find_first_not_of(legalInTail);
///             if (startBad == string::npos) {
///                 continue;
///             }
///             ... throw SShowStopper(...);
///         }
/// ```
/// (The pieces of a row are the tokens of its lines: the states run over them in order, which
/// gives the same result as over the whole row.)
fn verify_rows(rows: &[Vec<u8>]) -> Result<(), ScanFailure> {
    if rows.iter().all(|row| {
        row.iter()
            .all(|&c| c == b'-' || c == b'*' || c.is_ascii_alphabetic())
    }) {
        Ok(())
    } else {
        Err(ScanFailure)
    }
}
