//! The argument descriptions of the programs whose NCBI `-xmlhelp` LOSAT embeds
//! (`src/ncbi_xmlhelp/<program>.xml`, NCBI BLAST+ 2.17.0's and RpsbProc 0.5.1's XML), written as the
//! Rust tables of `ncbi_app::args` (`PROGRAMS` and the arrays it spans in `$OUT_DIR/ncbi_args.rs`,
//! which `src/ncbi_app/args.rs` includes): one producer of the definitions (line RG DECISIONS
//! K4i-1), read when LOSAT is built, so that a run reads no XML.
//!
//! The tables hold no pointer: every text is a span (`Text`, offset and length) of one string
//! constant, `TEXT`, and every list a span (`Span`) of an array of plain rows. A `&str` or a slice
//! in a static is an address that the dynamic loader relocates when LOSAT starts (one
//! `R_X86_64_RELATIVE` each in the PIE binary, some 2,800 for these tables), in every run of every
//! program; integers are not (`args_tests` checks the generated file for it).
//!
//! The XML is read as NCBI's `CArgDescriptions::PrintUsageXml` writes it: one element or line per
//! line, every text escaped by `NStr::XmlEncode`. An element, a type, a flag or a constraint that this
//! reader does not know stops the build with its text.
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1006-1051
//! (`CArgDesc::PrintUsageXml`: the start tag of a key or a flag with its name, a key's type and
//! `optional="true"` for an optional argument or a flag, then one line for each of its parts,
//! and a `<constraint>` element of the usage text and the constraint's own lines.)
//!
//! Four things the XML cannot carry are taken on trust. Each holds for the 16 XMLs LOSAT embeds,
//! and is to be checked again whenever an XML is replaced: for the 15 BLAST+ XMLs against the
//! 598d8ae6 source of the programs (2026-11-21, line RG K4i review L1), for rpsbproc's against
//! RpsbProc 0.5.1's source (2026-11-22, line RG K4i-2a). (a) and (d) cannot apply to rpsbproc (its
//! XML has no real bound and no integer key with permissible values); (b) and (c) hold for it as
//! said below.
//!
//! (a) A real bound is read back from its text, `%g` with six significant digits: the `<min>` and
//! `<max>` of `<Doubles>` and the numbers of BLAST's usage texts. A bound with more digits would
//! be checked rounded. The bounds of these programs are 0, 0.5, 1, 2, 4 and 100.
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4864-4872
//! NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:2316-2339
//!
//! (b) A constraint whose usage text is empty (an integer range from the least to the greatest
//! `int`, say) has no `<constraint>` element, so its table has no constraint, though NCBI's reader
//! checks it. No argument of these programs has one. rpsbproc constrains `-m` and `-t` alone, each
//! to the names of an enumeration (a set of strings, whose usage text is not empty), and passes no
//! constraint for its other arguments.
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1032-1043
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:211
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:227
//!
//! (c) `<default>` is the display value of a default, which is the default itself unless the
//! program gives another (a display value, or a default read from the environment): the table
//! takes it as the default. No argument of these programs gives one. rpsbproc's wrapper of
//! `AddDefaultKey` passes the default as its display value, and none of the arguments it installs
//! names an environment variable.
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1093-1095
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:1368-1371
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/argwrapper.cpp:33-34
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:152-279
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/rpsbproc_main.cpp:395
//!
//! (d) A usage text of permissible values on an integer key is read as BLAST's set of integers,
//! which compares numbers; BLAST's set of strings on such a key would write the same text and
//! compare the texts (`016` is not `16`). The one such key, blastn's `-template_length`, is a set
//! of integers.
//! NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:214-240

use std::fmt::Write as _;
use std::fs;
use std::path::Path;

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2207-2219
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2253-2350
// The constructor of the description adds the help flags; `AddStdArguments` adds `-logfile`,
// `-conffile`, `-version`, the three `-version-full` flags and `-dryrun` that the program has not,
// unless the program hides them.
/// The arguments that every NCBI application's description has (the help flags of the
/// `CArgDescriptions` constructor, `-logfile` and `-conffile` of `AddDefaultFileArguments`, and the
/// ones `AddStdArguments` adds before the XML is printed): the reader of `ncbi_app::args` knows the
/// ones of `HELP_FLAGS` and `FILE_ARGUMENTS`, and the application reads the others before the
/// reader (`ncbi_environment::preparse`), so the programs' tables leave them out.
const STANDARD: [&str; 11] = [
    "h",
    "help",
    "help-full",
    "xmlhelp",
    "logfile",
    "conffile",
    "version",
    "version-full",
    "version-full-xml",
    "version-full-json",
    "dryrun",
];

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2207-2219
// The constructor of `CArgDescriptions` adds the flags `-h` (with `m_AutoHelp`), `-help`,
// `-help-full` and `-xmlhelp` before any argument of the program. The XML does not list them.
/// The help flags that every reader of `ncbi_app::args` has before the program's arguments
/// (`HELP_FLAGS`). `HideStdArgs` takes effect later (`x_AddDefaultArgs`, after a help or an
/// argument error, before the usage is printed): `-help-full` and `-xmlhelp` are known to the
/// reader (probes C05, C06 of line BDC), and `-version` (which NCBI's application takes before the
/// reader, `ncbi_environment::preparse`) is not.
const HELP_FLAGS: [&str; 4] = ["h", "help", "help-full", "xmlhelp"];

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1213-1218
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2229-2242
// Before it reads the command line, the application adds to each description the optional keys
// `-logfile` (an output file) and `-conffile` (an input file) when the program has not.
/// The arguments that every reader of `ncbi_app::args` has after the program's
/// (`FILE_ARGUMENTS`, `AddDefaultFileArguments`): the name and the `Kind` variant of each. The
/// application opens `-logfile` itself (`ncbi_environment::diag_target`), and LOSAT refuses
/// `-conffile` before the reader (`ncbi_environment::preparse`).
const FILE_ARGUMENTS: [(&str, &str); 2] = [("logfile", "OutputFile"), ("conffile", "InputFile")];

/// How LOSAT refuses one of NCBI's arguments that it does not run (`ncbi_app::args::Refusal`).
enum Refused {
    /// By name, where its reader meets the argument: the refusal names the program and, if any,
    /// says why in parentheses after it.
    Read(Option<&'static str>),
    /// Where NCBI starts the search, whatever the value: NCBI's reader reads the argument, so that
    /// NCBI's own errors about it and the other arguments come first; the refusal says what the
    /// argument is.
    Run(&'static str),
}

/// The arguments of NCBI's programs that LOSAT refuses (AGENTS.md rule 2: explicit unsupported
/// errors): the program, the name and how (`Refused`). The one list of them (line RG DECISIONS
/// K4i-2): the readers' refusals and the run's (`cli::named_refusal`, `cli::run_refusal`) and
/// `-xmlhelp`, which leaves them out of NCBI's XML (`cli::xml_help_text`), read it. Un-refusing an
/// argument is deleting its row. The build checks that each row names an argument of its program's
/// XML, once.
const REFUSED: &[(&str, &str, Refused)] = &[
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastn_args.cpp:63-70
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:63-125
    // `-db` and `-dbsize` are implemented (`blastn/args.rs`), as are the X-drops,
    // `-window_size`, `-off_diagonal_range`, `-no_greedy` (DI-6), `-min_raw_gapped_score`
    // (DI-6b), the database masks (`blastinput::db_mask`, series NP M4a-2), `-filtering_db`
    // (`algorithm::blastn::repeat_filter`, line RG item D4) and `-best_hit_overhang` /
    // `-best_hit_score_edge` (series NP O6a).
    ("blastn", "entrez_query", Refused::Read(None)),
    ("blastn", "gilist", Refused::Read(None)),
    ("blastn", "html", Refused::Read(None)),
    ("blastn", "import_search_strategy", Refused::Read(None)),
    ("blastn", "negative_gilist", Refused::Read(None)),
    ("blastn", "remote", Refused::Read(None)),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastp_args.cpp:44-50
    // `-best_hit_overhang` / `-best_hit_score_edge` are implemented (series NP O6a-2, psiblast's
    // too).
    ("blastp", "entrez_query", Refused::Read(None)),
    ("blastp", "gilist", Refused::Read(None)),
    ("blastp", "html", Refused::Read(None)),
    ("blastp", "import_search_strategy", Refused::Read(None)),
    ("blastp", "ipglist", Refused::Read(None)),
    ("blastp", "negative_gilist", Refused::Read(None)),
    ("blastp", "negative_ipglist", Refused::Read(None)),
    ("blastp", "remote", Refused::Read(None)),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastx_args.cpp:62-66
    // `-best_hit_overhang` / `-best_hit_score_edge` are implemented (series NP O6a-3), as is
    // `-use_sw_tback` (line RG item X4b-2, #1500).
    ("blastx", "html", Refused::Read(None)),
    ("blastx", "gilist", Refused::Read(None)),
    ("blastx", "negative_gilist", Refused::Read(None)),
    ("blastx", "entrez_query", Refused::Read(None)),
    ("blastx", "ipglist", Refused::Read(None)),
    ("blastx", "negative_ipglist", Refused::Read(None)),
    ("blastx", "import_search_strategy", Refused::Read(None)),
    ("blastx", "remote", Refused::Read(None)),
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:246-250
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeblastdb.cpp:279-288
    // `-gi_mask` and `-gi_mask_name` make GI-indexed mask files, which LOSAT's writer does not
    // write yet (line NP, item M4c); `-mask_data` itself works.
    ("makeblastdb", "gi_mask", Refused::Read(Some(GI_MASK))),
    ("makeblastdb", "gi_mask_name", Refused::Read(Some(GI_MASK))),
    // NCBI reference (598d8ae6): c++/src/app/blastdb/makeprofiledb.cpp:398-400,503
    // The PSSM files in ASN.1 binary (psiblast line item PS10b) are refused where the run starts,
    // after makeprofiledb's `Init` (`makeprofiledb::work::run`).
    (
        "makeprofiledb",
        "binary",
        Refused::Run("PSSM files in ASN.1 binary"),
    ),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1156-1164
    // PHI-BLAST (PS8; the ASCII and ASN.1 PSSM files are PS2's and PS3's,
    // `algorithm::psiblast::pssm`, `-in_pssm` PS4's, `algorithm::psiblast::in_pssm`, the MSA
    // inputs PS6's, `algorithm::psiblast::in_msa`), and the arguments that blastp's engine does
    // not run either (blastp's rows; psiblast has no database masks, NCBI's unknown arguments of
    // its parser, and no `-mt_mode`). The search strategy file is written, of a FASTA query and of
    // `-in_pssm` (items F1b-1 and F1b-2 of line RG).
    ("psiblast", "phi_pattern", Refused::Read(None)),
    ("psiblast", "entrez_query", Refused::Read(None)),
    ("psiblast", "gilist", Refused::Read(None)),
    ("psiblast", "html", Refused::Read(None)),
    ("psiblast", "import_search_strategy", Refused::Read(None)),
    ("psiblast", "ipglist", Refused::Read(None)),
    ("psiblast", "negative_gilist", Refused::Read(None)),
    ("psiblast", "negative_ipglist", Refused::Read(None)),
    ("psiblast", "remote", Refused::Read(None)),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:917-918 (cited at
    // blastp's check, `validate_requested_blastp_support`)
    (
        "psiblast",
        "use_sw_tback",
        Refused::Run("the Smith-Waterman traceback"),
    ),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpsblast_args.cpp:85-127 (cited at
    // `algorithm::rpsblast::args`)
    // `-html` (PS15a) and the remote options and `-import_search_strategy` of every program.
    ("rpsblast", "entrez_query", Refused::Read(None)),
    ("rpsblast", "html", Refused::Read(None)),
    ("rpsblast", "import_search_strategy", Refused::Read(None)),
    ("rpsblast", "remote", Refused::Read(None)),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/rpstblastn_args.cpp:44-110 (cited
    // at `algorithm::rpstblastn::args`)
    // `-html` (PS15a) and the remote options and `-import_search_strategy` of every program. The
    // composition-based statistics of a later item are refused once per run
    // (`rpstblastn::search::TranslatedRpsSearch::check_options`).
    ("rpstblastn", "entrez_query", Refused::Read(None)),
    ("rpstblastn", "html", Refused::Read(None)),
    ("rpstblastn", "import_search_strategy", Refused::Read(None)),
    ("rpstblastn", "remote", Refused::Read(None)),
    (
        "rpstblastn",
        "ungapped",
        Refused::Run("NCBI BLAST+ 2.17.0's rpstblastn crashes after an ungapped search"),
    ),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastn_args.cpp:44-50
    // `-best_hit_overhang` / `-best_hit_score_edge` are implemented (series NP O6a-3), as is
    // `-use_sw_tback` (line RG item X4b-1, #1481).
    ("tblastn", "entrez_query", Refused::Read(None)),
    ("tblastn", "gilist", Refused::Read(None)),
    ("tblastn", "html", Refused::Read(None)),
    ("tblastn", "import_search_strategy", Refused::Read(None)),
    ("tblastn", "negative_gilist", Refused::Read(None)),
    ("tblastn", "remote", Refused::Read(None)),
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/tblastx_args.cpp:44-50
    // `-best_hit_overhang` / `-best_hit_score_edge` are implemented (series NP O6a-3).
    ("tblastx", "entrez_query", Refused::Read(None)),
    ("tblastx", "gilist", Refused::Read(None)),
    ("tblastx", "html", Refused::Read(None)),
    ("tblastx", "import_search_strategy", Refused::Read(None)),
    ("tblastx", "negative_gilist", Refused::Read(None)),
    ("tblastx", "remote", Refused::Read(None)),
];

/// Why makeblastdb's `-gi_mask` and `-gi_mask_name` are refused.
const GI_MASK: &str = "it writes GI-indexed mask files, which LOSAT does not write yet; -mask_data, -mask_id and -mask_desc work";

// No NCBI counterpart: PD-TLOSAN-LOCAL-GENCODE-32. NCBI's constraint allows the codes of its usage
// text (c++/src/algo/blast/blastinput/blast_args.cpp:998-1012, cited at `Constraint::GeneticCode`).
/// The keys whose genetic codes LOSAT widens by code 32 (gc.prt ID 32), keeping NCBI's usage text:
/// tblastn's `-db_gencode` (PD-TLOSAN-LOCAL-GENCODE-32, LOSAT's TBLASTN translates with it). The
/// build checks that each is a key of `CArgAllowGeneticCodeInteger`.
const GENETIC_CODE_32: [(&str, &str); 1] = [("tblastn", "db_gencode")];

/// One dependency of the XML: the first name, the `Dep` variant of `ncbi_app::args`, the second
/// name.
type Dependency = (String, &'static str, String);

/// One argument of the XML.
struct Argument {
    name: String,
    /// The `Kind` variant of `ncbi_app::args`.
    kind: &'static str,
    mandatory: bool,
    default: Option<String>,
    constraint: Option<Spec>,
    aliases: Vec<String>,
}

/// A `Constraint` of `ncbi_app::args`: the classes with a list of values, whose values go to the
/// value arrays, and the others, as Rust.
enum Spec {
    Strings(Vec<String>),
    StringSet(Vec<String>),
    IntegerSet(Vec<i32>),
    /// `CArgAllowGeneticCodeInteger`, with code 32 or not (`GENETIC_CODE_32`).
    GeneticCode {
        with_32: bool,
    },
    Rust(String),
}

/// The tables being written: the text that every `Text` spans, and the rows of each array.
#[derive(Default)]
struct Tables {
    text: String,
    descs: Vec<String>,
    aliases: Vec<String>,
    dependencies: Vec<String>,
    string_values: Vec<String>,
    integer_values: Vec<i32>,
    refusals: Vec<String>,
    programs: Vec<String>,
}

impl Tables {
    /// The `Text` of `value`: where it is in `text` (as a whole or inside a longer text), else
    /// added at the end.
    fn text(&mut self, value: &str) -> String {
        let at = self.text.find(value).unwrap_or_else(|| {
            self.text.push_str(value);
            self.text.len() - value.len()
        });
        let (at, len) = (u16::try_from(at), u16::try_from(value.len()));
        let (Ok(at), Ok(len)) = (at, len) else {
            panic!("the text of the argument tables passes 64 KiB")
        };
        format!("Text {{ at: {at}, len: {len} }}")
    }

    /// The `Span` of the rows of an array from `first` to its end (`len`).
    fn span(first: usize, len: usize) -> String {
        let first = u16::try_from(first).expect("fewer than 65,536 rows");
        let count = u16::try_from(len - usize::from(first)).expect("fewer than 65,536 rows");
        format!("Span {{ first: {first}, count: {count} }}")
    }

    /// The Rust of one `Desc`.
    fn desc(
        &mut self,
        name: &str,
        kind: &str,
        default: Option<&str>,
        constraint: Option<&Spec>,
        mandatory: bool,
    ) -> String {
        let name = self.text(name);
        let default = match default {
            Some(default) => format!("Some({})", self.text(default)),
            None => "None".to_string(),
        };
        let constraint = match constraint {
            Some(constraint) => format!("Some({})", self.constraint(constraint)),
            None => "None".to_string(),
        };
        format!(
            "Desc {{ name: {name}, kind: Kind::{kind}, default: {default}, constraint: \
             {constraint}, mandatory: {mandatory} }}"
        )
    }

    /// The Rust of one `Constraint`, its values added to the value arrays.
    fn constraint(&mut self, spec: &Spec) -> String {
        match spec {
            Spec::Strings(values) | Spec::StringSet(values) => {
                let first = self.string_values.len();
                for value in values {
                    let text = self.text(value);
                    self.string_values.push(text);
                }
                let span = Self::span(first, self.string_values.len());
                let class = if matches!(spec, Spec::Strings(_)) {
                    "Strings"
                } else {
                    "StringSet"
                };
                format!("Constraint::{class}({span})")
            }
            Spec::IntegerSet(values) => {
                let first = self.integer_values.len();
                self.integer_values.extend(values);
                let span = Self::span(first, self.integer_values.len());
                format!("Constraint::IntegerSet({span})")
            }
            Spec::GeneticCode { with_32 } => {
                format!("Constraint::GeneticCode {{ with_32: {with_32} }}")
            }
            Spec::Rust(rust) => rust.clone(),
        }
    }

    /// Adds the rows of one program.
    fn program(&mut self, program: &str, arguments: &[Argument], dependencies: &[Dependency]) {
        let first = self.descs.len();
        for argument in arguments {
            let desc = self.desc(
                &argument.name,
                argument.kind,
                argument.default.as_deref(),
                argument.constraint.as_ref(),
                argument.mandatory,
            );
            self.descs.push(desc);
        }
        let descs = Self::span(first, self.descs.len());
        let first = self.aliases.len();
        for argument in arguments {
            for alias in &argument.aliases {
                let row = format!("({}, {})", self.text(alias), self.text(&argument.name));
                self.aliases.push(row);
            }
        }
        let aliases = Self::span(first, self.aliases.len());
        let first = self.dependencies.len();
        for (first_name, dep, second) in dependencies {
            let row = format!(
                "({}, Dep::{dep}, {})",
                self.text(first_name),
                self.text(second)
            );
            self.dependencies.push(row);
        }
        let dependencies = Self::span(first, self.dependencies.len());
        let first = self.refusals.len();
        for (_, name, refused) in REFUSED.iter().filter(|row| row.0 == program) {
            let refusal = match refused {
                Refused::Read(None) => "Refused::Read(None)".to_string(),
                Refused::Read(Some(why)) => format!("Refused::Read(Some({}))", self.text(why)),
                Refused::Run(what) => format!("Refused::Run({})", self.text(what)),
            };
            let row = format!("({}, {refusal})", self.text(name));
            self.refusals.push(row);
        }
        let refusals = Self::span(first, self.refusals.len());
        let program = self.text(program);
        self.programs.push(format!(
            "({program}, Arguments {{ descs: {descs}, aliases: {aliases}, dependencies: \
             {dependencies}, refusals: {refusals} }})"
        ));
    }

    /// The Rust of the tables.
    fn write(mut self) -> String {
        let help_flags: Vec<String> = HELP_FLAGS
            .iter()
            .map(|name| self.desc(name, "Flag", None, None, false))
            .collect();
        let file_arguments: Vec<String> = FILE_ARGUMENTS
            .iter()
            .map(|(name, kind)| self.desc(name, kind, None, None, false))
            .collect();
        let mut out = String::from(
            "// Generated by LOSAT/build.rs (build/ncbi_args.rs) from LOSAT/src/ncbi_xmlhelp/*.xml.\n\
             // No pointer in a static: texts are spans of TEXT, lists spans of the arrays.\n",
        );
        let _ = writeln!(out, "const TEXT: &str = {:?};", self.text);
        let mut array = |name: &str, row: &str, rows: &[String]| {
            let _ = writeln!(out, "static {name}: [{row}; {}] = [", rows.len());
            for row in rows {
                let _ = writeln!(out, "    {row},");
            }
            out.push_str("];\n");
        };
        array("HELP_FLAGS", "Desc", &help_flags);
        array("FILE_ARGUMENTS", "Desc", &file_arguments);
        array("DESCS", "Desc", &self.descs);
        array("ALIASES", "(Text, Text)", &self.aliases);
        array("DEPENDENCIES", "(Text, Dep, Text)", &self.dependencies);
        array("STRING_VALUES", "Text", &self.string_values);
        let integers: Vec<String> = self.integer_values.iter().map(i32::to_string).collect();
        array("INTEGER_VALUES", "i32", &integers);
        array("REFUSALS", "(Text, Refused)", &self.refusals);
        array("PROGRAMS", "(Text, Arguments)", &self.programs);
        out
    }
}

/// Writes `$OUT_DIR/ncbi_args.rs` from the XML files of `manifest_dir/src/ncbi_xmlhelp`.
pub fn write_tables(manifest_dir: &Path, out_dir: &Path) {
    let dir = manifest_dir.join("src/ncbi_xmlhelp");
    println!("cargo:rerun-if-changed={}", dir.display());
    let mut programs: Vec<String> = fs::read_dir(&dir)
        .unwrap_or_else(|error| panic!("{}: {error}", dir.display()))
        .map(|entry| entry.expect("an entry of src/ncbi_xmlhelp").path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "xml"))
        .map(|path| {
            let stem = path.file_stem().expect("a file name").to_str();
            stem.expect("a UTF-8 file name").to_string()
        })
        .collect();
    programs.sort();
    let mut tables = Tables::default();
    for program in &programs {
        let path = dir.join(format!("{program}.xml"));
        let xml =
            fs::read_to_string(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        let (mut arguments, dependencies) =
            read_xml(&xml).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
        check_refused(program, &arguments);
        for (_, name) in GENETIC_CODE_32.iter().filter(|row| row.0 == program) {
            let argument = arguments.iter_mut().find(|argument| argument.name == *name);
            match argument.and_then(|argument| argument.constraint.as_mut()) {
                Some(Spec::GeneticCode { with_32 }) => *with_32 = true,
                _ => panic!("GENETIC_CODE_32: {program} -{name} is no key of a genetic code"),
            }
        }
        tables.program(program, &arguments, &dependencies);
    }
    let rows = REFUSED
        .iter()
        .map(|row| row.0)
        .chain(GENETIC_CODE_32.iter().map(|row| row.0));
    for program in rows {
        assert!(
            programs.iter().any(|known| known == program),
            "a row of REFUSED or GENETIC_CODE_32 names {program}, whose XML LOSAT has not"
        );
    }
    let target = out_dir.join("ncbi_args.rs");
    let out = tables.write();
    fs::write(&target, out).unwrap_or_else(|error| panic!("{}: {error}", target.display()));
}

/// Checks the rows of `REFUSED` of `program`: each names an argument of its XML (`arguments`),
/// once.
fn check_refused(program: &str, arguments: &[Argument]) {
    let rows: Vec<&str> = REFUSED
        .iter()
        .filter(|row| row.0 == program)
        .map(|row| row.1)
        .collect();
    for (at, name) in rows.iter().enumerate() {
        assert!(
            arguments.iter().any(|argument| argument.name == *name),
            "REFUSED: {program} has no argument -{name} in its XML"
        );
        assert!(
            !rows[..at].contains(name),
            "REFUSED: {program} -{name} is refused twice"
        );
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3899-3915
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3925-3943
// The keys are written first, then the flags, each in the order of the description, each with
// its aliases; the dependencies follow the arguments.
/// The arguments (but the standard ones) and the dependencies of one XML, in the XML's order.
fn read_xml(xml: &str) -> Result<(Vec<Argument>, Vec<Dependency>), String> {
    let mut lines = xml.lines();
    lines
        .by_ref()
        .find(|line| *line == "<arguments>")
        .ok_or("no <arguments> line (or a positional mode)")?;
    let mut arguments = Vec::new();
    let mut dependencies = Vec::new();
    while let Some(line) = lines.next() {
        if line == "</arguments>" {
            return Ok((arguments, dependencies));
        }
        if line == "<dependencies>" {
            read_dependencies(&mut lines, &mut dependencies)?;
            continue;
        }
        let argument = read_argument(line, &mut lines)?;
        if !STANDARD.contains(&argument.name.as_str()) {
            arguments.push(argument);
        }
    }
    Err("no </arguments> line".to_string())
}

/// The argument whose element starts at `line` (`<key ...>` or `<flag ...>`), up to its end tag.
fn read_argument<'a>(
    line: &str,
    lines: &mut impl Iterator<Item = &'a str>,
) -> Result<Argument, String> {
    let (tag, rest) = line
        .strip_prefix('<')
        .and_then(|rest| rest.split_once(' '))
        .filter(|(tag, _)| matches!(*tag, "key" | "flag"))
        .ok_or_else(|| format!("an element this reader does not know: {line}"))?;
    let rest = rest.strip_prefix("name=\"").ok_or(line)?;
    let (name, rest) = rest.split_once('"').ok_or(line)?;
    let (kind, rest) = if tag == "flag" {
        ("Flag", rest)
    } else {
        let rest = rest.strip_prefix(" type=\"").ok_or(line)?;
        let (type_name, rest) = rest.split_once('"').ok_or(line)?;
        (kind_of(type_name).ok_or(line)?, rest)
    };
    let optional = match rest {
        " optional=\"true\">" => true,
        ">" if tag == "key" => false,
        _ => return Err(format!("attributes this reader does not know: {line}")),
    };
    let mut argument = Argument {
        name: decode(name)?,
        kind,
        mandatory: !optional,
        default: None,
        constraint: None,
        aliases: Vec::new(),
    };
    let end = format!("</{tag}>");
    loop {
        let line = lines.next().ok_or_else(|| format!("no {end}"))?;
        if line == end {
            break;
        }
        if line == "<constraint>" {
            let description = lines.next().and_then(|line| element(line, "description"));
            let description = decode(description.ok_or(line)?)?;
            argument.constraint = Some(read_constraint(kind, &description, lines)?);
            continue;
        }
        match element_of(line) {
            Some(("description" | "group" | "synopsis", _)) => {}
            Some(("default", value)) => argument.default = Some(decode(value)?),
            Some(("alias", value)) => argument.aliases.push(decode(value)?),
            // `fHidden` (makeblastdb's `-randomize`): the argument is read as any other.
            _ if line == "<flags><hidden/></flags>" => {}
            // `fAppend` of makeblastdb's `-logfile`, a standard argument that the tables leave out.
            _ if line == "<flags><append/></flags>" && argument.name == "logfile" => {}
            _ => return Err(format!("{name}: a line this reader does not know: {line}")),
        }
    }
    if argument.mandatory && argument.default.is_some() {
        return Err(format!("{name}: a mandatory key with a default"));
    }
    Ok(argument)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:2387-2400
/// The `Kind` of a key's type name (`CArgDescriptions::GetTypeName`).
fn kind_of(type_name: &str) -> Option<&'static str> {
    Some(match type_name {
        "String" => "Text",
        "Boolean" => "Boolean",
        "Int8" => "Int8",
        "Integer" => "Integer",
        "Real" => "Double",
        "File_In" => "InputFile",
        "File_Out" => "OutputFile",
        _ => return None,
    })
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4639-4652
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4729-4741
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:4864-4872
// A set of strings writes its case sensitivity and one `<value>` per string; a range of integers
// (written `Integers`, or `Int8s` for 64 bits) or of reals one `<min>` and `<max>` pair per
// interval, the reals with `NStr::DoubleToString`.
/// The `Constraint` of a `<constraint>` element whose usage text is `description`, the lines after
/// it up to `</constraint>`. NCBI's toolkit classes write their values (`<Strings>`, `<Doubles>`,
/// `<Integers>`); BLAST's write only their usage text (`CArgAllow::PrintUsageXml` writes nothing),
/// which names one class in these programs.
fn read_constraint<'a>(
    kind: &str,
    description: &str,
    lines: &mut impl Iterator<Item = &'a str>,
) -> Result<Spec, String> {
    let line = lines.next().ok_or("no </constraint>")?;
    if line == "</constraint>" {
        return blast_constraint(kind, description);
    }
    let constraint = match line {
        "<Strings case_sensitive=\"true\">" => {
            let values = read_values(lines, "</Strings>", &["value"])?;
            Spec::Strings(values.into_iter().map(|(_, value)| value).collect())
        }
        "<Doubles>" | "<Integers>" => {
            let end = if line == "<Doubles>" {
                "</Doubles>"
            } else {
                "</Integers>"
            };
            let values = read_values(lines, end, &["min", "max"])?;
            let [(_, min), (_, max)] = values.as_slice() else {
                return Err(format!("not one interval: {description}"));
            };
            if line == "<Doubles>" {
                let min: f64 = min.parse().map_err(|_| min.clone())?;
                let max: f64 = max.parse().map_err(|_| max.clone())?;
                Spec::Rust(format!(
                    "Constraint::Doubles {{ min: {min:?}, max: {max:?} }}"
                ))
            } else {
                let min: i64 = min.parse().map_err(|_| min.clone())?;
                let max: i64 = max.parse().map_err(|_| max.clone())?;
                Spec::Rust(format!("Constraint::Integers {{ min: {min}, max: {max} }}"))
            }
        }
        _ => return Err(format!("a constraint this reader does not know: {line}")),
    };
    match lines.next() {
        Some("</constraint>") => Ok(constraint),
        other => Err(format!(
            "{description}: {other:?} in place of </constraint>"
        )),
    }
}

/// The `<tag>value</tag>` lines up to `end`, each tag one of `tags` in turn.
fn read_values<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
    end: &str,
    tags: &[&str],
) -> Result<Vec<(String, String)>, String> {
    let mut values = Vec::new();
    loop {
        let line = lines.next().ok_or_else(|| format!("no {end}"))?;
        if line == end {
            return Ok(values);
        }
        let tag = tags[values.len() % tags.len()];
        let value = element(line, tag).ok_or_else(|| format!("{line} in place of <{tag}>"))?;
        values.push((tag.to_string(), decode(value)?));
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:79-185
// NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_input_aux.hpp:214-240
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:998-1013
// BLAST's own constraints write no values of their own, only their usage texts: the ranges at
// least, at most and between, the sets of permissible values and the genetic codes.
/// The class of BLAST's constraint whose `GetUsage` is `usage` (`blast_input_aux.hpp`,
/// `blast_args.cpp`; the `Constraint` of each says how NCBI writes it), on a key of `kind`.
fn blast_constraint(kind: &str, usage: &str) -> Result<Spec, String> {
    let number = |text: &str| -> Result<f64, String> {
        text.parse::<f64>()
            .map_err(|_| format!("a bound that is no number: {usage}"))
    };
    if let Some(rest) = usage.strip_prefix("(>=") {
        let (min, max) = rest
            .strip_suffix(')')
            .and_then(|rest| rest.split_once(" and =<"))
            .ok_or_else(|| format!("a range this reader does not know: {usage}"))?;
        let (min, max) = (number(min)?, number(max)?);
        return Ok(Spec::Rust(format!(
            "Constraint::Between {{ min: {min:?}, max: {max:?}, inclusive: true }}"
        )));
    }
    if let Some(rest) = usage.strip_prefix("(>") {
        let (min, max) = rest
            .strip_suffix(')')
            .and_then(|rest| rest.split_once(" and <"))
            .ok_or_else(|| format!("a range this reader does not know: {usage}"))?;
        let (min, max) = (number(min)?, number(max)?);
        return Ok(Spec::Rust(format!(
            "Constraint::Between {{ min: {min:?}, max: {max:?}, inclusive: false }}"
        )));
    }
    if let Some(min) = usage.strip_prefix(">=") {
        return Ok(Spec::Rust(format!(
            "Constraint::AtLeast({:?})",
            number(min)?
        )));
    }
    if let Some(max) = usage.strip_prefix("<=") {
        return Ok(Spec::Rust(format!(
            "Constraint::AtMost({:?})",
            number(max)?
        )));
    }
    if let Some(rest) = usage.strip_prefix("Permissible values: ") {
        let values: Vec<&str> = rest
            .strip_suffix(' ')
            .map(|rest| rest.split(' ').collect())
            .unwrap_or_default();
        let values: Option<Vec<&str>> = values
            .iter()
            .map(|value| value.strip_prefix('\'')?.strip_suffix('\''))
            .collect();
        let values = values
            .filter(|values| !values.is_empty())
            .ok_or_else(|| format!("a set this reader does not know: {usage}"))?;
        // `CArgAllowIntegerSet` on an integer key (`-template_length`), `CArgAllowStringSet` on a
        // string key (`-task`).
        return match kind {
            "Integer" => {
                let values: Result<Vec<i32>, _> =
                    values.iter().map(|value| value.parse()).collect();
                let values =
                    values.map_err(|_| format!("an integer set of other values: {usage}"))?;
                Ok(Spec::IntegerSet(values))
            }
            "Text" => Ok(Spec::StringSet(
                values.iter().map(|value| value.to_string()).collect(),
            )),
            _ => Err(format!("a set of values on a key of kind {kind}: {usage}")),
        };
    }
    if usage == "values between: 1-6, 9-16, 21-31, 33" {
        return Ok(Spec::GeneticCode { with_32: false });
    }
    if let Some(max) = usage.strip_prefix("file name length < ") {
        let max: u32 = max
            .parse()
            .map_err(|_| format!("a length this reader does not know: {usage}"))?;
        return Ok(Spec::Rust(format!("Constraint::MaxFileNameLength({max})")));
    }
    Err(format!("a constraint that maps to no class: {usage}"))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3925-3943
/// The dependencies of a `<dependencies>` element, up to `</dependencies>`: NCBI writes its
/// `multimap` twice, the `eRequires` entries first, then the `eExcludes` ones, each in the map's
/// order (an exclusion is in the map in both directions).
fn read_dependencies<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
    dependencies: &mut Vec<Dependency>,
) -> Result<(), String> {
    loop {
        let line = lines.next().ok_or("no </dependencies>")?;
        let dep = match line {
            "</dependencies>" => return Ok(()),
            "<first_requires_second>" => "Requires",
            "<first_excludes_second>" => "Excludes",
            _ => return Err(format!("a dependency this reader does not know: {line}")),
        };
        let first = lines.next().and_then(|line| element(line, "arg1"));
        let second = lines.next().and_then(|line| element(line, "arg2"));
        let (Some(first), Some(second)) = (first, second) else {
            return Err(format!("a dependency without two names after {line}"));
        };
        dependencies.push((decode(first)?, dep, decode(second)?));
        let end = format!("</{}", &line[1..]);
        if lines.next() != Some(end.as_str()) {
            return Err(format!("no {end}"));
        }
    }
}

/// The text of `line` = `<tag>text</tag>` (`s_WriteXmlLine`).
fn element<'a>(line: &'a str, tag: &str) -> Option<&'a str> {
    line.strip_prefix('<')?
        .strip_prefix(tag)?
        .strip_prefix('>')?
        .strip_suffix('>')?
        .strip_suffix(tag)?
        .strip_suffix("</")
}

/// The tag and the text of a line `<tag>text</tag>`.
fn element_of(line: &str) -> Option<(&str, &str)> {
    let (tag, _) = line.strip_prefix('<')?.split_once('>')?;
    Some((tag, element(line, tag)?))
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:4041-4058
/// The text that `NStr::XmlEncode` wrote as `text`: its five entities and the character references
/// (`&#xa;`) of the control characters.
fn decode(text: &str) -> Result<String, String> {
    let mut out = String::with_capacity(text.len());
    let mut rest = text;
    while let Some(at) = rest.find('&') {
        out.push_str(&rest[..at]);
        let (entity, after) = rest[at + 1..]
            .split_once(';')
            .ok_or_else(|| format!("an entity without its end: {text}"))?;
        let decoded = match entity {
            "amp" => '&',
            "lt" => '<',
            "gt" => '>',
            "apos" => '\'',
            "quot" => '"',
            _ => entity
                .strip_prefix("#x")
                .and_then(|hex| u32::from_str_radix(hex, 16).ok())
                .and_then(char::from_u32)
                .ok_or_else(|| format!("an entity this reader does not know: &{entity};"))?,
        };
        out.push(decoded);
        rest = after;
    }
    out.push_str(rest);
    Ok(out)
}
