//! The command line of `losat rpsbproc` (NCBI's rpsbproc 0.5.1): `-i`, `-o`, `-e`, `-m`, `-t`,
//! `-d`, `-f`, `-q` and their long forms `--infile`, `--outfile`, `--evalue`, `--data-mode`,
//! `--target-data`, `--data-path`, `--show-families`, `--quiet`.
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:152-279
//! ```c++
//! 	{
//! 		"i",	//argInFile
//! 		TArgDefinition::eKey,	//enum TArgDefinition::EArgCategory(argument type: eKey, ePos, eOpenPos, eFlag, eNegFlag)
//! 		ncbi::CArgDescriptions::eString,	//ncbi::CArgDescriptions::EType(process type)
//! 		"Input file(s) that holds asn data produced by rpsblast with \"-outfmt 11\" switch. If omitted, default to stdin.",	//string (description)
//! 		true,	//if optional
//! 		TArgDefinition::EMPTYSTR,	//Synopsis (short description)
//! 		"-infile",	//alias (can be used to specify at command line but not indexed in program
//! 	...
//! ```
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/argwrapper.cpp:29-43
//! ```c++
//! 			case TArgDefinition::eKey:
//! 				if (argdef.optional)
//! 				{
//!
//! 					if (TArgDefinition::EMPTYSTR != argdef.def_val || TArgDefinition::EMPTYSTR != argdef.val_env)
//! 						pMyArgs->AddDefaultKey(argdef.name, argdef.synopsis, argdef.descr, argdef.type, argdef.def_val, argdef.flags, argdef.val_env, argdef.def_val);
//! 					else
//!                     {
//! 						pMyArgs->AddOptionalKey(argdef.name, argdef.synopsis, argdef.descr, argdef.type, argdef.flags);
//!                     }
//! 				}
//! 				...
//! 				if (TArgDefinition::EMPTYSTR != argdef.alias)
//! 					pMyArgs->AddAlias(argdef.alias, argdef.name);
//! ```
//! An alias `-infile` of `i` is the argument `--infile` on NCBI's command line (`-` and the
//! alias name). The parser, its errors and `-help` are LOSAT's (`PD-LOSAT-CLI-NONSEARCH-
//! DIFFERENCES`); the values are read as NCBI reads them.

use clap::Args;
use std::ffi::OsString;

use super::domains::DataMode;

/// `TTargetData`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TargetData {
    Doms,
    Feats,
    Both,
}

/// The arguments of `losat rpsbproc`.
#[derive(Args, Debug, Clone)]
pub struct RpsbprocArgs {
    #[arg(
        long = "i",
        value_name = "FILE",
        value_parser = clap::builder::OsStringValueParser::new(),
        help = "Input file that holds asn data produced by rpsblast with \"-outfmt 11\" switch (also --infile). If omitted, default to stdin."
    )]
    pub infile: Option<OsString>,
    #[arg(
        long = "o",
        value_name = "FILE",
        value_parser = clap::builder::OsStringValueParser::new(),
        help = "Output file to write the processed result (also --outfile). If omitted, default to stdout."
    )]
    pub outfile: Option<OsString>,
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:184-199
    // ```c++
    // 		"e",	//argEVCutoff
    // 		TArgDefinition::eKey,
    // 		ncbi::CArgDescriptions::eDouble,
    // 		"EValue cut-off. Program will only process hits with evalues better than this value.",
    // 		...
    // 		NumericDefVal(DEF_EVALCUTOFF),	//Default Value
    // ```
    // (`ParseBlastOutput` skips an HSP whose E-value is above it, `offline::OfflineQuery`; the
    // archive path does not pass it, quoted at `run_sessions`.)
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:1485-1486
    // ```c++
    // constexpr const double PRECALC_EVALUE = 0.01;	//default evalue cut off
    // constexpr const double DEF_EVALCUTOFF = PRECALC_EVALUE;	//default evalue cut off
    // ```
    #[arg(
        long = "e",
        value_name = "REAL",
        default_value = "0.01",
        value_parser = rpsbproc_real,
        help = "E-value cut-off, shown in the header (also --evalue); it drops the HSPs above it in the XML input and selects no hit of an ASN.1 archive."
    )]
    pub evalue: f64,
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:103-104
    // ```c++
    // const char * TDataModes::dimLits[] = {rep, std, full};
    // const char * TDataModes::dimLabels[] = {"Concise", "Standard", "Full"};
    // ```
    #[arg(
        long = "m",
        value_name = "MODE",
        default_value = "rep",
        value_parser = ["rep", "std", "full"],
        help = "Redundancy level of domain hit data: rep (concise), std (standard) or full (all hits) (also --data-mode)"
    )]
    pub data_mode: String,
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:198-199
    // ```c++
    // const char * TTargetData::dimLits[] = {"doms", "feats", "both"};
    // const char * TTargetData::dimLabels[] = {"Domain hits", "Site annotations", "Domain hits and site annotations"};
    // ```
    #[arg(
        long = "t",
        value_name = "DATA",
        default_value = "both",
        value_parser = ["doms", "feats", "both"],
        help = "Target data: doms (domain hits), feats (sites and motifs) or both (also --target-data)"
    )]
    pub target_data: String,
    #[arg(
        long = "d",
        value_name = "DIR",
        value_parser = clap::builder::OsStringValueParser::new(),
        help = "Directory that contains the data files (also --data-path). By default, the 'data' directory in the current directory and where the program is located"
    )]
    pub data_path: Option<OsString>,
    #[arg(
        long = "f",
        help = "Show the superfamilies of the domain hits (also --show-families)"
    )]
    pub show_families: bool,
    #[arg(
        long = "q",
        help = "Quiet mode: no progress on stderr (also --quiet; the program information is still shown)"
    )]
    pub quiet: bool,
}

/// A real argument as NCBI's `CArg_Double` reads it (`ncbi_double`).
fn rpsbproc_real(value: &str) -> Result<f64, String> {
    crate::blastinput::value_parsers::ncbi_double(value)
}

/// The long forms of rpsbproc's arguments (`AddAlias`): `--infile` is `-i`, and so on.
const ALIASES: [(&str, &str); 8] = [
    ("--infile", "-i"),
    ("--outfile", "-o"),
    ("--evalue", "-e"),
    ("--data-mode", "-m"),
    ("--target-data", "-t"),
    ("--data-path", "-d"),
    ("--show-families", "-f"),
    ("--quiet", "-q"),
];

/// The argument word `text` with an alias replaced by its argument (`--infile=X` is `-i=X`).
pub(crate) fn alias_word(text: &str) -> Option<String> {
    let (name, value) = match text.split_once('=') {
        Some((name, value)) => (name, Some(value)),
        None => (text, None),
    };
    let (_, short) = ALIASES.iter().find(|(alias, _)| *alias == name)?;
    Some(match value {
        Some(value) => format!("{short}={value}"),
        None => (*short).to_string(),
    })
}

impl RpsbprocArgs {
    /// `-m`.
    pub(crate) fn mode(&self) -> DataMode {
        match self.data_mode.as_str() {
            "std" => DataMode::Std,
            "full" => DataMode::Full,
            _ => DataMode::Rep,
        }
    }

    /// `-t`.
    pub(crate) fn target(&self) -> TargetData {
        match self.target_data.as_str() {
            "doms" => TargetData::Doms,
            "feats" => TargetData::Feats,
            _ => TargetData::Both,
        }
    }
}
