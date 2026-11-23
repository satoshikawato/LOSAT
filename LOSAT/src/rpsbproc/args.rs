//! The command line of `losat rpsbproc` (NCBI's rpsbproc 0.5.1): `-i`, `-o`, `-e`, `-m`, `-t`,
//! `-d`, `-f`, `-q` and their long forms `--infile`, `--outfile`, `--evalue`, `--data-mode`,
//! `--target-data`, `--data-path`, `--show-families`, `--quiet`.
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:152-279
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/argwrapper.cpp:29-43
//! An alias `-infile` of `i` is the argument `--infile` on NCBI's command line (`-` and the
//! alias name). The parser, its errors and `-help` are LOSAT's (`PD-LOSAT-CLI-NONSEARCH-
//! DIFFERENCES`); the values are read as NCBI reads them.

use crate::blastinput::value_parsers::ncbi_double;
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
    // (`ParseBlastOutput` skips an HSP whose E-value is above it, `offline::OfflineQuery`; the
    // archive path does not pass it, cited at `run_sessions`.)
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:1485-1486
    #[arg(
        long = "e",
        value_name = "REAL",
        default_value = "0.01",
        value_parser = ncbi_double,
        help = "E-value cut-off, shown in the header (also --evalue); it drops the HSPs above it in the XML input and selects no hit of an ASN.1 archive."
    )]
    pub evalue: f64,
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:103-104
    #[arg(
        long = "m",
        value_name = "MODE",
        default_value = "rep",
        help = "Redundancy level of domain hit data: rep (concise), std (standard) or full (all hits) (also --data-mode)"
    )]
    pub data_mode: String,
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:198-199
    #[arg(
        long = "t",
        value_name = "DATA",
        default_value = "both",
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

impl RpsbprocArgs {
    // No NCBI counterpart: NCBI's programs read `CArgs` where their option handles take each value
    // (line RG DECISIONS K4i-3-1).
    /// The options of `args`, rpsbproc's arguments as NCBI's argument reader read them
    /// (`ncbi_app::args`): each key as given, or its default in NCBI's tables where the option
    /// has one, each flag whether it was given, and the run's state, which the run fills in at
    /// its start.
    pub(crate) fn from_args(
        args: &crate::ncbi_app::args::Args,
    ) -> Result<Self, crate::ncbi_app::args::NotConverted> {
        Ok(RpsbprocArgs {
            infile: args.given_option("i")?,
            outfile: args.given_option("o")?,
            evalue: args.option_or_default("e")?,
            data_mode: args.option_or_default("m")?,
            target_data: args.option_or_default("t")?,
            data_path: args.given_option("d")?,
            show_families: args.flag("f"),
            quiet: args.flag("q"),
        })
    }
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
