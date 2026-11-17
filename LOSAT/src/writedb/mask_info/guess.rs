//! The format guess of `makeblastdb -mask_data`: `CFormatGuess::GuessFormat` with the preferred
//! formats `eBinaryASN` and `eTextASN` and every other format disabled.
//!
//! `GuessFormat` reads a test buffer (the one port of `CFormatGuess`, `utils::format_guess`), then tries
//! the formats in the order of `sm_CheckOrder`; of the two that are enabled `eTextASN` comes
//! first, `eBinaryASN` last. A file that is neither (or that has no test buffer: empty, or
//! made of comments that end before the buffer is full) is `eUnknown`. A file with no test buffer
//! gets one more test, `TestFormatNewick`, whose answer `eNewick` is just as unusable for
//! `s_ReadObject` as `eUnknown`: "Unknown encoding".

use std::io::Read;

use super::MaskFileFormat;
use crate::utils::format_guess::FormatGuess;

// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:453-509
// ```c++
// CFormatGuess::EFormat
// CFormatGuess::GuessFormat( EOnError onerror )
// {
//     ...
//     if (!EnsureTestBuffer()) {
//         ...
//         if(TestFormatNewick(eQuick)) {
//             return CFormatGuess::eNewick;
//         }
//         return CFormatGuess::eUnknown;
//     }
//
//     EMode mode = eQuick;
//
//     // First, try to use hints
//     if ( !m_Hints.IsEmpty() ) {
//         for (size_t f = 0; f < sm_CheckOrder_Size; ++f) {
//             EFormat fmt = EFormat( sm_CheckOrder[f] );
//             if (m_Hints.IsPreferred(fmt)  &&  x_TestFormat(fmt, mode)) {
//                 return fmt;
//             }
//         }
//     }
//
//     // Check other formats, skip the ones that are disabled through hints
//     for (size_t f = 0; f < sm_CheckOrder_Size; ++f) {
//         EFormat fmt = EFormat( sm_CheckOrder[f] );
//         if ( ! m_Hints.IsDisabled(fmt)  &&  x_TestFormat(fmt, mode) ) {
//             return fmt;
//         }
//     }
//     return eUnknown;
// }
// ```
// NCBI reference (598d8ae6): c++/include/util/format_guess.hpp:420-424
// ```c++
// CFormatGuess::CFormatHints::DisableAllNonpreferred(void)
// {
//     m_Disabled = ~m_Preferred;
//     return *this;
// }
// ```
// NCBI reference (598d8ae6): c++/src/util/format_guess.cpp:46-88
// ```c++
// static const CFormatGuess::EFormat sm_CheckOrder[] =
// {
//     ...
//     CFormatGuess::eFasta,
//     CFormatGuess::eTextASN,
//     CFormatGuess::eAlignment,
//     CFormatGuess::eTaxplot,
//     CFormatGuess::eTable,
//     CFormatGuess::eBinaryASN,
//     ...
// };
// ```
/// The format of the mask file `input`, `None` for `eUnknown` (and `eNewick`), with the bytes
/// taken from `input` for the guess (the start of the file; NCBI gives them back to the stream).
pub(super) fn guess_format<R: Read>(input: &mut R) -> Option<(MaskFileFormat, Vec<u8>)> {
    let mut guess = FormatGuess::new(input);
    if !guess.ensure_test_buffer() {
        return None;
    }
    let format = if guess.test_format_text_asn() {
        MaskFileFormat::TextAsn
    } else if guess.test_format_binary_asn() {
        MaskFileFormat::BinaryAsn
    } else {
        return None;
    };
    Some((format, guess.into_bytes()))
}
