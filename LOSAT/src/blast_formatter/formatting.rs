//! What blast_formatter makes of its formatting arguments for one archive object: NCBI's
//! `CFormattingArgs(false, CFormattingArgs::eIsSAM)::ExtractAlgorithmOptions` on the object's
//! options (blast_args.cpp:2870-2990), then blast_formatter's own check of SAM
//! (blast_formatter.cpp:351-355). NCBI runs both in `PrintFormattedOutput`, after the object was
//! read, so `-outfmt` is checked there and not with the other arguments. The readers are the
//! search programs': `-outfmt` is `blastinput::app::parse_formatting_string` and
//! `formatting_handler_check_with` (with `eIsSAM`), the counts and warnings
//! `blastinput::report_counts` (from the object's hit list size, with no prefix on the posts).

use crate::blastinput::app::{
    formatting_handler_check_with, options_error, output_format, parse_formatting_string,
    FormatChoice,
};
use crate::blastinput::report_counts::{
    kept_hitlist_counts, report_counts_from, FormattingArgs, ReportCounts,
};

use super::Options;

/// The formatting of one archive object (`CFormattingArgs` after `ExtractAlgorithmOptions`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Formatting {
    /// The output format and its custom specification (`m_OutputFormat`,
    /// `m_CustomOutputFormatSpec`, `m_CustomDelim`).
    pub(crate) choice: FormatChoice,
    /// The hit list size (`opt.SetHitlistSize`) and the numbers of descriptions and alignments.
    pub(crate) counts: ReportCounts,
    /// `m_LineLength`, `m_HitsSortOption` and `m_HspsSortOption` as the handler keeps them: the
    /// line length and the sort of the hits for the formats up to 4, the sort of the HSPs for
    /// format 0, none when the handler kept the object's hit list size (`kept_hitlist_counts`).
    pub(crate) formatting: FormattingArgs,
    /// `-show_gis` and `-html` (`m_ShowGis`, `m_Html`).
    pub(crate) show_gis: bool,
    pub(crate) html: bool,
    /// The handler's warnings (`ERR_POST(Warning << ...)`), in NCBI's order, without a prefix.
    pub(crate) warnings: Vec<u8>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2870-2894
// ```c++
// void
// CFormattingArgs::ExtractAlgorithmOptions(const CArgs& args,
//                                          CBlastOptions& opt)
// {
//     ParseFormattingString(args, m_OutputFormat, m_CustomOutputFormatSpec,m_CustomDelim);
//     if((m_OutputFormat == eSAM) && !(m_FormatFlags & eIsSAM) ){
//     		NCBI_THROW(CInputException, eInvalidInput,
//     		                        "SAM format is only applicable to blastn" );
//     }
//     ...
//     s_ValidateCustomDelim(m_CustomOutputFormatSpec,m_CustomDelim);
//     m_ShowGis = static_cast<bool>(args[kArgShowGIs]);
//     if(m_IsIgBlast){
//         m_Html = false;
//     } else {
//         m_Html = static_cast<bool>(args[kArgProduceHtml]);
//     }
//     // Default hitlist size 500, value can be changed if import search strategy is used
//     int hitlist_size = opt.GetHitlistSize();
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2938-2945
// ```c++
//     	if (args[kArgLineLength]) {
//     	    m_LineLength = args[kArgLineLength].AsInteger();
//     	}
//         if(args.Exist(kArgSortHits) && args[kArgSortHits])
//         {
//        	    m_HitsSortOption = args[kArgSortHits].AsInteger();
//         }
//     }
// ```
// (the branch of the formats up to 4; the other formats keep the defaults; the sort of the HSPs
// is kept for format 0 alone, quoted at `ReportCounts::warnings_with`.)
/// `ExtractAlgorithmOptions` of the command line `options` for an archive object whose options
/// hold the hit list size `hitlist_size` (`opt.GetHitlistSize()`, the request's; 500 when it has
/// none): NCBI's errors of `-outfmt` (`'x' is not a valid output format`, the delimiter's, the
/// formats of other programs; `Formatting choice is out of range` with exit status 255), or the
/// formatting and its warnings.
pub(crate) fn extract_formatting(
    options: &Options,
    hitlist_size: usize,
) -> anyhow::Result<Formatting> {
    let choice = parse_formatting_string(&options.outfmt)?;
    formatting_handler_check_with(&choice, true)?;
    let outfmt = choice.number;
    let (show_gis, html) = (options.show_gis, options.html);
    if let Some(counts) = kept_hitlist_counts(
        hitlist_size,
        outfmt,
        options.num_descriptions.is_some(),
        options.num_alignments.is_some(),
    ) {
        return Ok(Formatting {
            choice,
            counts,
            formatting: FormattingArgs::default(),
            show_gis,
            html,
            warnings: Vec::new(),
        });
    }
    let counts = report_counts_from(
        hitlist_size,
        outfmt,
        options.max_target_seqs,
        options.num_descriptions,
        options.num_alignments,
    );
    let given = FormattingArgs {
        line_length: options.line_length,
        sort_hits: options.sort_hits,
        sort_hsps: options.sort_hsps,
        hitlist_start: None,
    };
    let warnings = counts.posted_warnings(None, outfmt, options.num_descriptions.is_some(), &given);
    let pairwise_family = outfmt <= output_format::FLAT_QUERY_ANCHORED_NO_IDENTITIES;
    Ok(Formatting {
        choice,
        counts,
        formatting: FormattingArgs {
            line_length: given.line_length.filter(|_| pairwise_family),
            sort_hits: given.sort_hits.filter(|_| pairwise_family),
            sort_hsps: given
                .sort_hsps
                .filter(|_| outfmt == output_format::PAIRWISE),
            hitlist_start: None,
        },
        show_gis,
        html,
        warnings,
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:350-355
// ```c++
//     const EBlastProgramType p = opts.GetProgramType();
//     if((fmt_args.GetFormattedOutputChoice() == CFormattingArgs::eSAM) &&
//        (p != eBlastTypeBlastn )) {
// 		NCBI_THROW(CInputException, eInvalidInput,
// 		                        "SAM format is only applicable to blastn results" );
//     }
// ```
/// blast_formatter's check of SAM (17) for an archive object of a program that is blastn
/// (`blastn`, its tasks megablast and dc-megablast included) or not (`is_blastn`), after
/// `extract_formatting` (probes P2.197 and P2.199: blastp archives).
pub(crate) fn sam_check(choice: &FormatChoice, is_blastn: bool) -> anyhow::Result<()> {
    if choice.number == output_format::SAM && !is_blastn {
        return Err(options_error(
            "SAM format is only applicable to blastn results",
        ));
    }
    Ok(())
}

#[cfg(test)]
#[path = "formatting_tests.rs"]
mod formatting_tests;
