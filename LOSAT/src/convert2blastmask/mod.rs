//! `losat convert2blastmask`: NCBI BLAST+ 2.17.0 convert2blastmask
//! (`CConvert2BlastMaskApplication`, c++/src/app/blastdb/convert2blastmask.cpp), which turns the
//! lower case letters of FASTA input, and its lines `from - to`, into the masking information
//! that makeblastdb's `-mask_data` reads. Its command line is NCBI's (`ncbi_app`); its reader is
//! NCBI's `CFastaReader` with the flags of `CMaskFromFasta`
//! (`writedb::fasta_source::reader::ReaderFlags::mask_from_fasta`) and its writer the masking
//! programs' `CMaskWriterBlastDbMaskInfo` (`seqmasks_io::writer`). The work of its `Run` is in
//! `work`.

mod work;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::ncbi_app::args::{flag, key, mandatory_key, Constraint, Desc, Kind};
use crate::ncbi_app::{no_dependencies, Tool};

/// The name in LOSAT's refusals.
const PROGRAM: &str = "CONVERT2BLASTMASK";

// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:141-145
// ```c++
// const char * const CConvert2BlastMaskApplication::USAGE_LINE
//     = "Convert masking information in lower-case masked FASTA input to file formats suitable for makeblastdb";
//
// void CConvert2BlastMaskApplication::Init(void) {
//     HideStdArgs(fHideLogfile | fHideConffile | fHideFullVersion | fHideXmlHelp | fHideDryRun);
// ```
/// The text of `convert2blastmask -help` (NCBI 2.17.0's, frozen): `-logfile`, `-conffile`,
/// `-version-full`, `-xmlhelp` and `-dryrun` are hidden. Its part before `OPTIONAL ARGUMENTS`
/// is the usage that `-h` and the argument errors print.
const HELP: &str = include_str!("help.txt");

// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:154-186
// ```c++
//     arg_desc->AddDefaultKey("in", "input_file_name",
//                             "Input file name",
//                             CArgDescriptions::eInputFile, "-");
//
//     arg_desc->AddDefaultKey("out", "output_file_name",
//                             "Output file name",
//                             CArgDescriptions::eOutputFile, "-");
//
//     arg_desc->AddDefaultKey("outfmt", "output_format",
//                             "Output file format",
//                             CArgDescriptions::eString, "maskinfo_asn1_text");
//
//     CArgAllow_Strings* strings_allowed = new CArgAllow_Strings();
//     strings_allowed->Allow("maskinfo_asn1_text");
//     strings_allowed->Allow("maskinfo_asn1_bin");
//     strings_allowed->Allow("maskinfo_xml");
//     strings_allowed->Allow("interval");
//     arg_desc->SetConstraint("outfmt", strings_allowed);
//
//     arg_desc->AddFlag      ( "parse_seqids",
//                              "Parse Seq-ids in FASTA input", true );
//
//     arg_desc->AddKey       ("masking_algorithm", "mask_program_name",
//                             "Masking algorithm name (e.g.: dust, seg, "
//                             "windowmasker, repeat). Use 'other' for "
//                             "user-defined type",
//                             CArgDescriptions::eString);
//
//     arg_desc->AddKey     ("masking_options", "mask_program_options",
//                           "Masking algorithm options to create the masked input"
//                           " (free text to describe/include (command line) "
//                           "options used to create the masking)",
//                           CArgDescriptions::eString);
// ```
/// The arguments of `CConvert2BlastMaskApplication::Init` in the order it describes them; the
/// values that `-outfmt` allows are listed in the order of their `std::set`.
const DESCS: &[Desc] = &[
    key("in", Kind::InputFile, Some("-"), None),
    key("out", Kind::OutputFile, Some("-"), None),
    key(
        "outfmt",
        Kind::Text,
        Some("maskinfo_asn1_text"),
        Some(Constraint::Strings(&[
            "interval",
            "maskinfo_asn1_bin",
            "maskinfo_asn1_text",
            "maskinfo_xml",
        ])),
    ),
    flag("parse_seqids"),
    mandatory_key("masking_algorithm", Kind::Text),
    mandatory_key("masking_options", Kind::Text),
];

// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:111-120
// ```c++
//     CConvert2BlastMaskApplication() {
//         CRef<CVersion> version(new CVersion());
//         version->SetVersionInfo(new CBlastVersion());
//         SetFullVersion(version);
//         m_StopWatch.Start();
//         if (m_UsageReport.IsEnabled()) {
//         	m_UsageReport.AddParam(CBlastUsageReport::eVersion, GetVersion().Print());
//         	m_UsageReport.AddParam(CBlastUsageReport::eProgram, (string) "convert2blastmask");
//         }
//     }
// ```
/// convert2blastmask as NCBI's application runs it (`ncbi_app`): BLAST's version, BLAST's
/// usage report (`ncbi_environment`), no dependencies between its arguments and no prefix on
/// its diagnostics.
static TOOL: Tool = Tool {
    program: "convert2blastmask",
    refusal_name: PROGRAM,
    version: crate::NCBI_BLAST_VERSION,
    help: HELP,
    descs: DESCS,
    dependencies: no_dependencies,
    post_prefix: None,
};

// NCBI reference (598d8ae6): c++/src/app/blastdb/convert2blastmask.cpp:259-263
// ```c++
// int main(int argc, const char* argv[])
// {
//     // Execute main application function
//     return CConvert2BlastMaskApplication().AppMain(argc, argv);
// }
// ```
/// Runs `convert2blastmask` with `words` (the command line after `losat convert2blastmask`) on
/// the process's standard streams; the exit status.
pub(crate) fn run(words: &[OsString]) -> i32 {
    crate::ncbi_app::run(&TOOL, words, work::run)
}

/// The run of one command line with `input` as its standard input, its output on `out` and its
/// errors on `err`; the exit status (NCBI's application around convert2blastmask's `Run`:
/// `ncbi_app::run_with`). The library API's `run_convert2blastmask`.
pub(crate) fn run_with(
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
) -> i32 {
    crate::ncbi_app::run_with(&TOOL, words, input, out, err, work::run)
}

#[cfg(test)]
mod tests;
