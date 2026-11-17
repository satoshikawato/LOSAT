//! NCBI's BLAST publications and the start of the pairwise report that cites them: the table of
//! `CReference` (api/version.cpp), the version line (`CBlastFormatUtil::BlastPrintVersionInfo`),
//! one reference (`CBlastFormatUtil::BlastPrintReference`, blastfmtutil.cpp) and the references of
//! `CBlastFormat::PrintProlog` (blast_format.cpp), in their plain and their `-html` forms.
//!
//! The pairwise prologs of every program (`pairwise::write_blastp_pairwise_intro` and the
//! others) and the XML reports (`blast_xml::REFERENCE_GAPPED_BLAST`) read the table here.

use std::borrow::Cow;
use std::io::{self, Write};

/// NCBI's BLAST publications (`CReference::EPublication`), in its order; `MaxPublications` is the
/// sentinel, whose text and URL are empty.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/version.hpp:69-79
/// ```c++
///     enum EPublication {
///         eGappedBlast = 0,           ///< 1997 NAR paper
///         ePhiBlast,                  ///< 1998 NAR paper
///         eMegaBlast,                 ///< 2000 J Compt Biol paper
///         eCompBasedStats,            ///< 2001 NAR paper
///         eCompAdjustedMatrices,      ///< submitted for publication
///         eIndexedMegablast,          ///< 2008 Bioinformatics on indexed megablast
///         eDeltaBlast,          	    ///< 2012 Biology Direct on DeltaBLAST
///         eMaxPublications            ///< Used as sentinel value
///     };
/// ```
/// (The whole table is ported: PHI-BLAST, indexed megablast and DELTA-BLAST, which LOSAT does
/// not run, are read only by the tests of NCBI's unit cases.)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[allow(dead_code)]
pub(crate) enum Publication {
    GappedBlast = 0,
    PhiBlast,
    MegaBlast,
    CompBasedStats,
    CompAdjustedMatrices,
    IndexedMegablast,
    DeltaBlast,
    MaxPublications,
}

/// The number of entries of the tables (`eMaxPublications + 1`).
const TABLE_LEN: usize = Publication::MaxPublications as usize + 1;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/version.cpp:46-82
// ```c++
// static const string kReferences[(int)CReference::eMaxPublications+1] = {
//     // eGappedBlast
//     "Stephen F. Altschul, Thomas L. Madden, \
// Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. \
// Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of protein \
// database search programs\", Nucleic Acids Res. 25:3389-3402.",
//     // ePhiBlast
//     "Zheng Zhang, Alejandro A. Sch&auml;ffer, Webb Miller, \
// Thomas L. Madden, David J. Lipman, Eugene V. Koonin, and Stephen F. \
// Altschul (1998), \"Protein sequence similarity searches using patterns \
// as seeds\", Nucleic Acids Res. 26:3986-3990.",
//     // eMegaBlast
//     "Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), \
// \"A greedy algorithm for aligning DNA sequences\", \
// J Comput Biol 2000; 7(1-2):203-14.",
//     // eCompBasedStats
//     "Alejandro A. Sch&auml;ffer, L. Aravind, Thomas L. Madden, Sergei Shavirin, \
// John L. Spouge, Yuri I. Wolf, Eugene V. Koonin, and Stephen F. Altschul \
// (2001), \"Improving the accuracy of PSI-BLAST protein database searches \
// with composition-based statistics and other refinements\", Nucleic Acids \
// Res. 29:2994-3005.",
//     // eCompAdjustedMatrices
//     "Stephen F. Altschul, John C. Wootton, E. Michael Gertz, Richa Agarwala, \
// Aleksandr Morgulis, Alejandro A. Sch&auml;ffer, and Yi-Kuo Yu (2005) \"Protein \
// database searches using compositionally adjusted substitution matrices\", \
// FEBS J. 272:5101-5109.",
//     // eIndexedMegablast
//     "Aleksandr Morgulis, George Coulouris, Yan Raytselis, \
// Thomas L. Madden, Richa Agarwala, Alejandro A. Sch&auml;ffer \
// (2008), \"Database Indexing for Production MegaBLAST Searches\", \
// Bioinformatics 24:1757-1764.",
//     // eDeltaBlast
//     "Grzegorz M. Boratyn, Alejandro A. Schaffer, Richa Agarwala, Stephen F. Altschul, \
// David J. Lipman and Thomas L. Madden (2012) \"Domain enhanced lookup time \
// accelerated BLAST\", Biology Direct 7:12.",
//     // eMaxPublications
//     kEmptyStr
// };
// ```
// (A backslash at a line's end continues the C++ string literal without a newline. The space
// before the newline after `7(1-2):203-14.",` is outside the literal.)
const REFERENCES: [&str; TABLE_LEN] = [
    "Stephen F. Altschul, Thomas L. Madden, Alejandro A. Sch&auml;ffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of protein database search programs\", Nucleic Acids Res. 25:3389-3402.",
    "Zheng Zhang, Alejandro A. Sch&auml;ffer, Webb Miller, Thomas L. Madden, David J. Lipman, Eugene V. Koonin, and Stephen F. Altschul (1998), \"Protein sequence similarity searches using patterns as seeds\", Nucleic Acids Res. 26:3986-3990.",
    "Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb Miller (2000), \"A greedy algorithm for aligning DNA sequences\", J Comput Biol 2000; 7(1-2):203-14.",
    "Alejandro A. Sch&auml;ffer, L. Aravind, Thomas L. Madden, Sergei Shavirin, John L. Spouge, Yuri I. Wolf, Eugene V. Koonin, and Stephen F. Altschul (2001), \"Improving the accuracy of PSI-BLAST protein database searches with composition-based statistics and other refinements\", Nucleic Acids Res. 29:2994-3005.",
    "Stephen F. Altschul, John C. Wootton, E. Michael Gertz, Richa Agarwala, Aleksandr Morgulis, Alejandro A. Sch&auml;ffer, and Yi-Kuo Yu (2005) \"Protein database searches using compositionally adjusted substitution matrices\", FEBS J. 272:5101-5109.",
    "Aleksandr Morgulis, George Coulouris, Yan Raytselis, Thomas L. Madden, Richa Agarwala, Alejandro A. Sch&auml;ffer (2008), \"Database Indexing for Production MegaBLAST Searches\", Bioinformatics 24:1757-1764.",
    "Grzegorz M. Boratyn, Alejandro A. Schaffer, Richa Agarwala, Stephen F. Altschul, David J. Lipman and Thomas L. Madden (2012) \"Domain enhanced lookup time accelerated BLAST\", Biology Direct 7:12.",
    "",
];

// NCBI reference (598d8ae6): c++/src/algo/blast/api/version.cpp:84-106
// ```c++
// static const string kPubMedUrls[(int)CReference::eMaxPublications+1] = {
//     // eGappedBlast
//     "//www.ncbi.nlm.nih.gov/\
// entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=9254694&dopt=Citation",
//     // ePhiBlast
//     "//www.ncbi.nlm.nih.gov/\
// entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=9705509&dopt=Citation",
//     // eMegaBlast
//     "//www.ncbi.nlm.nih.gov/\
// entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=10890397&dopt=Citation",
//     // eCompBasedStats
//     "//www.ncbi.nlm.nih.gov/\
// entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=11452024&dopt=Citation",
//     // eCompAdjustedMatrices
//     "//www.ncbi.nlm.nih.gov/\
// entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=16218944&dopt=Citation",
//     // eIndexedMegablast
//     "//www.ncbi.nlm.nih.gov/pubmed/18567917",
//     // eDeltaBlast
//     "//www.ncbi.nlm.nih.gov/pubmed/22510480",
//     // eMaxPublications
//     kEmptyStr
// };
// ```
const PUBMED_URLS: [&str; TABLE_LEN] = [
    "//www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=9254694&dopt=Citation",
    "//www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=9705509&dopt=Citation",
    "//www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=10890397&dopt=Citation",
    "//www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=11452024&dopt=Citation",
    "//www.ncbi.nlm.nih.gov/entrez/query.fcgi?db=PubMed&cmd=Retrieve&list_uids=16218944&dopt=Citation",
    "//www.ncbi.nlm.nih.gov/pubmed/18567917",
    "//www.ncbi.nlm.nih.gov/pubmed/22510480",
    "",
];

impl Publication {
    /// `CReference::GetString`: the publication's text, with NCBI's HTML entity (`&auml;`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/version.cpp:108-112
    /// ```c++
    /// string
    /// CReference::GetString(EPublication pub)
    /// {
    ///     return kReferences[(int) pub];
    /// }
    /// ```
    pub(crate) const fn string(self) -> &'static str {
        REFERENCES[self as usize]
    }

    /// `CReference::GetHTMLFreeString`: the text with its first `&auml;` written `a` (no text
    /// has two).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/version.cpp:114-123
    /// ```c++
    /// string
    /// CReference::GetHTMLFreeString(EPublication pub)
    /// {
    ///     string pub_string = GetString(pub);
    ///     string::size_type offset = pub_string.find("&auml;");
    ///     if (offset != string::npos)
    ///         pub_string.replace(offset, 6, "a");
    ///
    ///     return pub_string;
    /// }
    /// ```
    pub(crate) fn html_free_string(self) -> Cow<'static, str> {
        let text = self.string();
        match text.find("&auml;") {
            Some(offset) => {
                Cow::Owned([&text[..offset], "a", &text[offset + "&auml;".len()..]].concat())
            }
            None => Cow::Borrowed(text),
        }
    }

    /// `CReference::GetPubmedUrl`: the publication's PubMed URL without its scheme.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/version.cpp:125-129
    /// ```c++
    /// string
    /// CReference::GetPubmedUrl(EPublication pub)
    /// {
    ///     return kPubMedUrls[(int) pub];
    /// }
    /// ```
    pub(crate) const fn pubmed_url(self) -> &'static str {
        PUBMED_URLS[self as usize]
    }
}

/// The scheme of NCBI's URLs where `[BLASTFMTUTIL] PROTOCOL` of `./.ncbirc` does not set one
/// (`httpProt`, blastfmtutil.cpp:100; `report::formatter_registry::protocol`).
pub(crate) const DEFAULT_PROTOCOL: &[u8] = b"https:";

/// The line length of NCBI's prolog (`kFormatLineLength`), at which the references are wrapped.
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/format/blastfmtutil.hpp:57
/// ```c++
/// const size_t kFormatLineLength = 68;
/// ```
const FORMAT_LINE_LENGTH: usize = 68;

/// `CBlastFormatUtil::BlastPrintVersionInfo`: `program` in upper case and `version`
/// (`BlastGetVersion`), in bold under `-html`, and a newline.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:60-73
/// ```c++
/// string CBlastFormatUtil::BlastGetVersion(const string program)
/// {
///     string program_uc = program;
///     return NStr::ToUpper(program_uc) + " " + blast::CBlastVersion().Print();
/// }
///
/// void CBlastFormatUtil::BlastPrintVersionInfo(const string program, bool html,
///                                              CNcbiOstream& out)
/// {
///     if (html)
///         out << "<b>" << BlastGetVersion(program) << "</b>" << "\n";
///     else
///         out << BlastGetVersion(program) << "\n";
/// }
/// ```
pub(crate) fn write_version_line<W: Write + ?Sized>(
    out: &mut W,
    html: bool,
    program: &str,
    version: &str,
) -> io::Result<()> {
    let program = program.to_ascii_uppercase();
    if html {
        writeln!(out, "<b>{program} {version}</b>")
    } else {
        writeln!(out, "{program} {version}")
    }
}

/// `CBlastFormatUtil::BlastPrintReference`: the reference's label and the publication's text,
/// wrapped as one FlatFile string at 68 columns (`x_WrapOutputLine` without `html`), and a blank
/// line. Under `html` the label is a bold link to the publication's PubMed URL with the scheme
/// `protocol` (`[BLASTFMTUTIL] PROTOCOL`, `DEFAULT_PROTOCOL`) and ends the first part of the
/// string with a newline, before the text with its entity; the wrap counts the tags, so it can
/// break the URL. Otherwise the label and `: ` lead the text without its entity.
/// `is_psiblast` names the composition-based statistics reference psiblast's (`starting in round
/// 2`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastfmtutil.cpp:75-121
/// ```c++
/// void
/// CBlastFormatUtil::BlastPrintReference(bool html, size_t line_len,
///                                       CNcbiOstream& out,
///                                       blast::CReference::EPublication pub,
///                                       bool is_psiblast /* = false */)
/// {
///     string reference("Reference");
///     if (pub == blast::CReference::eCompAdjustedMatrices) {
///         reference += " for compositional score matrix adjustment";
///     } else if (pub == blast::CReference::eCompBasedStats) {
///         reference += " for composition-based statistics";
///         if (is_psiblast) {
///             reference += " starting in round 2";
///         }
///     } else if (pub == blast::CReference::eIndexedMegablast) {
///         reference += " for database indexing";
///     } else if (pub == blast::CReference::eDeltaBlast) {
///         reference += " for DELTA-BLAST";
///     }
///
///     ostringstream str;
///     if(html)
///     {
///         CNcbiIfstream  config_file(".ncbirc");
///         CNcbiRegistry config_reg(config_file);
///         string httpProt = "https:";
///         if(!config_reg.Empty()) {
///             if(config_reg.HasEntry("BLASTFMTUTIL","PROTOCOL")) {
///                 httpProt = config_reg.Get("BLASTFMTUTIL","PROTOCOL");
///             }
///         }
///         str << "<b><a href=\""
///             << httpProt
///             << blast::CReference::GetPubmedUrl(pub)
///             << "\">" << reference << "</a>:</b>"
///             << "\n";
///         x_WrapOutputLine(str.str() + blast::CReference::GetString(pub),
///                    line_len, out);
///     }
///     else
///     {
///         str << reference << ": ";
///         x_WrapOutputLine(str.str() + blast::CReference::GetHTMLFreeString(pub),
///                    line_len, out);
///     }
///
///     out << "\n";
/// }
/// ```
/// (The read of `./.ncbirc` is the caller's: `protocol`.)
pub(crate) fn write_reference<W: Write + ?Sized>(
    out: &mut W,
    html: bool,
    protocol: &[u8],
    publication: Publication,
    is_psiblast: bool,
) -> io::Result<()> {
    let mut reference = String::from("Reference");
    match publication {
        Publication::CompAdjustedMatrices => {
            reference.push_str(" for compositional score matrix adjustment");
        }
        Publication::CompBasedStats => {
            reference.push_str(" for composition-based statistics");
            if is_psiblast {
                reference.push_str(" starting in round 2");
            }
        }
        Publication::IndexedMegablast => reference.push_str(" for database indexing"),
        Publication::DeltaBlast => reference.push_str(" for DELTA-BLAST"),
        _ => {}
    }
    let mut text = Vec::new();
    if html {
        text.extend_from_slice(b"<b><a href=\"");
        text.extend_from_slice(protocol);
        text.extend_from_slice(publication.pubmed_url().as_bytes());
        text.extend_from_slice(b"\">");
        text.extend_from_slice(reference.as_bytes());
        text.extend_from_slice(b"</a>:</b>\n");
        text.extend_from_slice(publication.string().as_bytes());
    } else {
        text.extend_from_slice(reference.as_bytes());
        text.extend_from_slice(b": ");
        text.extend_from_slice(publication.html_free_string().as_bytes());
    }
    super::pairwise::write_flatfile_wrapped(out, &text, FORMAT_LINE_LENGTH)?;
    crate::cli::insert(out, b"\n")
}

/// The references of `CBlastFormat::PrintProlog` after the version line and its two newlines, and
/// the two newlines before the database report: the megablast reference for a megablast search
/// (`megablast`), else the Gapped BLAST reference; then psiblast's compositionally adjusted
/// matrices, and the composition-based statistics of psiblast (`starting in round 2`) and blastp,
/// each after a newline. `program` is NCBI's `m_Program` (`blastp`, `psiblast`, `rpsblast`,
/// `blastn`, ...; compared without case, so the callers may pass the version line's name).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:395-437
/// ```c++
///     m_Outfile << NcbiEndl << NcbiEndl;
///     if (m_Program == "deltablast") {
///         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
///                               m_Outfile, CReference::eDeltaBlast);
///         m_Outfile << "\n";
///     }
///
///     if (m_Megablast)
///         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
///                                           m_Outfile, CReference::eMegaBlast);
///     else
///         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
///                                           m_Outfile);
///
///     if (m_Megablast && m_IndexedMegablast)
///     {
///         m_Outfile << "\n";
///         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
///                               m_Outfile, CReference::eIndexedMegablast);
///     }
///
///     if (m_Program == "psiblast" || m_Program == "deltablast") {
///         m_Outfile << "\n";
///         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
///                               m_Outfile, CReference::eCompAdjustedMatrices);
///     }
///     if (m_Program == "psiblast" || m_Program == "blastp") {
///         m_Outfile << "\n";
///         CBlastFormatUtil::BlastPrintReference(m_IsHTML, kFormatLineLength,
///                               m_Outfile, CReference::eCompBasedStats,
///                               (bool)(m_Program == "psiblast"));
///     }
///
///     if (m_Program == "deltablast" || !m_DomainDbInfo.empty()) {
///         m_Outfile << "\n\n";
///         ...
///     }
///     else {
///         m_Outfile << "\n\n";
///     }
/// ```
/// (The `NcbiEndl` pair is the callers': their flushes differ. LOSAT has no deltablast, refuses
/// `-use_index true` with an index (`m_IndexedMegablast`) and has no domain database: their
/// branches are not taken.)
pub(crate) fn write_prolog_references<W: Write + ?Sized>(
    out: &mut W,
    html: bool,
    protocol: &[u8],
    program: &str,
    megablast: bool,
) -> io::Result<()> {
    let psiblast = program.eq_ignore_ascii_case("psiblast");
    let primary = if megablast {
        Publication::MegaBlast
    } else {
        Publication::GappedBlast
    };
    write_reference(out, html, protocol, primary, false)?;
    if psiblast {
        crate::cli::insert(out, b"\n")?;
        write_reference(
            out,
            html,
            protocol,
            Publication::CompAdjustedMatrices,
            false,
        )?;
    }
    if psiblast || program.eq_ignore_ascii_case("blastp") {
        crate::cli::insert(out, b"\n")?;
        write_reference(out, html, protocol, Publication::CompBasedStats, psiblast)?;
    }
    crate::cli::insert(out, b"\n\n")
}

#[cfg(test)]
#[path = "reference_tests.rs"]
mod tests;
