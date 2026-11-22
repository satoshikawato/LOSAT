//! The request of an archive object (`Blast4-request`): what NCBI's archive reader takes from it
//! (`CRemoteBlast::x_GetRequestInfoFromFile` through `CImportStrategy`) and what blast_formatter
//! reads of its queries and subjects. The sequences, locations and PSSM stay the values read
//! (`asn_types::Value`), with typed access to what the readers ask of them.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1717-1741

use super::Stop;
use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::serial::asn_types::table::named;
use crate::report::serial::asn_types::view::V;
use crate::report::serial::asn_types::{member_index, members, Value};
use crate::seqid::rank::{best_rank_score, find_best};
use crate::seqid::ParsedSeqId;

/// The request of a search.
#[derive(Debug)]
pub(crate) struct Request {
    /// `program` and `service` (`GetProgram`, `GetService`), as stored.
    pub(crate) program: Vec<u8>,
    pub(crate) service: Vec<u8>,
    pub(crate) queries: Queries,
    pub(crate) subject: Subject,
    /// The parameters of `algorithm-options`, `program-options` and `format-options` in the
    /// order of the file (`params`).
    params: Vec<Param>,
}

/// The queries of a request (`Blast4-queries`), as read.
#[derive(Debug)]
pub(crate) enum Queries {
    /// `pssm`: the `PssmWithParameters`.
    Pssm(Value),
    /// `seq-loc-list`: the `Seq-loc`s.
    SeqLocList(Vec<Value>),
    /// `bioseq-set`: the `Bioseq-set`.
    BioseqSet(Value),
}

/// The subject of a request (`Blast4-subject`), as read.
#[derive(Debug)]
pub(crate) enum Subject {
    /// `database`: the name as stored (`x_SetDatabase`).
    Database(Vec<u8>),
    /// `sequences`: the `Bioseq`s (`m_SubjectSequences`).
    Sequences(Vec<Value>),
    /// `seq-loc-list`: the `Seq-loc`s. NCBI's archive reader does not read them: its
    /// `SetSequences()` makes the subject an empty list of sequences, in the request it keeps
    /// too, so blast_formatter has no subject.
    SeqLocList,
}

/// The NCBI code that reads a request: the archive reader (`CRemoteBlast::LoadFromArchive`,
/// then `x_GetRequestInfoFromFile`) or the strategy import (`CImportStrategy`).
#[derive(Clone, Copy, PartialEq, Eq)]
enum Reader {
    Archive,
    Strategy,
}

/// The list of options a parameter comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ParamSet {
    /// `algorithm-options` (`GetAlgoOptions`).
    Algorithm,
    /// `program-options` (`GetProgramOptions`).
    Program,
    /// `format-options` (`GetWebFormatOptions`).
    Format,
}

/// A `Blast4-parameter`, as read.
#[derive(Debug)]
pub(crate) struct Param {
    pub(crate) set: ParamSet,
    pub(crate) name: Vec<u8>,
    /// The `Blast4-value`.
    value: Value,
}

impl Param {
    /// The `Blast4-value` (a CHOICE: `choice()` names its variant).
    pub(crate) fn value(&self) -> V<'_> {
        V::new(named::BLAST4_VALUE, &self.value)
    }
}

/// The member `name` of the SEQUENCE `value` of the type `ty`, moved out.
fn take(value: &mut Value, ty: u16, name: &str) -> Option<Value> {
    value.take_member(member_index(ty, name))
}

/// The variant of the CHOICE `value` of the type `ty`: its name and value, moved out.
fn take_variant(value: &mut Value, ty: u16) -> (&'static str, Value) {
    match value {
        Value::Choice(index, variant) => (members(ty)[usize::from(*index)].name(), variant.take()),
        other => unreachable!("a CHOICE of the table, not {other:?}"),
    }
}

/// The elements of a SEQUENCE OF value, moved out.
fn take_items(value: &mut Value) -> Vec<Value> {
    match value {
        Value::List(items) => std::mem::take(items),
        other => unreachable!("a list of values of the table, not {other:?}"),
    }
}

/// The bytes of a string value, moved out.
fn take_str(value: &mut Value) -> Vec<u8> {
    match value {
        Value::Str(bytes) => std::mem::take(bytes),
        other => unreachable!("a string of the table, not {other:?}"),
    }
}

impl Request {
    /// The request of a `Blast4-request` value: `None` for an error archive's
    /// (`IsErrMsgArchive`: a `get-request-info` whose `request-id` is "Error"); a body that is not
    /// a `queue-search` stops `CImportStrategy`, and then a request without `ident`
    /// `GetCreatedBy` (`x_GetRequestInfoFromFile` asks for the program and service, then for it).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2683-2693
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:62-75
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:164-168
    // NCBI reference (598d8ae6): c++/include/objects/blast/Blast4_request_.hpp:191-198
    // (`GetProgram` and `GetService` read mandatory members of the `queue-search`.)
    pub(crate) fn from_value(request: Value) -> Result<Option<Request>, Stop> {
        Request::read(request, Reader::Archive)
    }

    /// The request of a search strategy (`-import_search_strategy`): `CImportStrategy` over the
    /// `Blast4-request` the program read. A body that is not a `queue-search` stops it (an error
    /// archive's `get-request-info` too); its import never asks for `ident` (NCBI's exports
    /// write none; `E.*.ident` of line F2's probes: ignored).
    pub(crate) fn from_strategy(request: Value) -> Result<Request, NcbiException> {
        match Request::read(request, Reader::Strategy) {
            Ok(request) => Ok(request.expect("a strategy's request is a search")),
            Err(Stop::Exception(exception)) => Err(exception),
            Err(stop) => unreachable!("a strategy's import stops only with an exception: {stop:?}"),
        }
    }

    fn read(mut request: Value, reader: Reader) -> Result<Option<Request>, Stop> {
        let ty = named::BLAST4_REQUEST;
        // `ident` (`GetCreatedBy`, the version of the program that wrote the archive) changes no
        // output (probes P8.001-012); NCBI's archive reader stops on a search request without it.
        let has_ident = take(&mut request, ty, "ident").is_some();
        let mut body = take(&mut request, ty, "body").expect("the mandatory body");
        let (variant, mut search) = take_variant(&mut body, named::BLAST4_REQUEST_BODY);
        if reader == Reader::Archive && variant == "get-request-info" {
            let id = V::new(named::BLAST4_GET_REQUEST_INFO_REQUEST, &search).field("request-id");
            if id.str() == b"Error" {
                return Ok(None);
            }
        }
        if variant != "queue-search" {
            return Err(Stop::Exception(NcbiException::new(
                ThrowSite {
                    file: "src/algo/blast/api/search_strategy.cpp",
                    line: 72,
                    kind: "CBlastException::eInvalidArgument",
                    function: Some("BLAST::ncbi::blast::CImportStrategy::CImportStrategy"),
                },
                "No body in CBlast4_request",
            )));
        }
        if reader == Reader::Archive && !has_ident {
            return Err(Stop::Unassigned {
                type_name: "Blast4-request",
                member: "ident",
            });
        }
        let ty = named::BLAST4_QUEUE_SEARCH_REQUEST;
        let mut text = |name: &str| {
            let mut value = take(&mut search, ty, name).expect("a mandatory string");
            take_str(&mut value)
        };
        let (program, service) = (text("program"), text("service"));
        let mut queries = take(&mut search, ty, "queries").expect("the mandatory queries");
        let queries = match take_variant(&mut queries, named::BLAST4_QUERIES) {
            ("pssm", pssm) => Queries::Pssm(pssm),
            ("seq-loc-list", mut locs) => Queries::SeqLocList(take_items(&mut locs)),
            (_, set) => Queries::BioseqSet(set),
        };
        let mut subject = take(&mut search, ty, "subject").expect("the mandatory subject");
        let subject = match take_variant(&mut subject, named::BLAST4_SUBJECT) {
            ("database", mut name) => Subject::Database(take_str(&mut name)),
            ("sequences", mut bioseqs) => Subject::Sequences(take_items(&mut bioseqs)),
            _ => Subject::SeqLocList,
        };
        let mut params = Vec::new();
        for (set, name) in [
            (ParamSet::Algorithm, "algorithm-options"),
            (ParamSet::Program, "program-options"),
            (ParamSet::Format, "format-options"),
        ] {
            let Some(mut list) = take(&mut search, ty, name) else {
                continue;
            };
            params.extend(take_items(&mut list).into_iter().map(|mut param| {
                let ty = named::BLAST4_PARAMETER;
                let name = take_str(&mut take(&mut param, ty, "name").expect("the name"));
                let value = take(&mut param, ty, "value").expect("the value");
                Param { set, name, value }
            }));
        }
        Ok(Some(Request {
            program,
            service,
            queries,
            subject,
            params,
        }))
    }

    /// The parameters of the request: `algorithm-options`, `program-options` and
    /// `format-options` in the order of the file (each list in its order, a name as often as it
    /// is stored). The one accessor of the request's options: NCBI's option builder reads the
    /// first two lists in this order, each parameter in turn (`x_ProcessOptions`), and
    /// `GetParamByName` the first parameter of a name in a list ([`Request::param`]).
    pub(crate) fn params(&self) -> &[Param] {
        &self.params
    }

    /// `GetParamByName(name)` of the list `set`: its first parameter named `name`.
    // NCBI reference (598d8ae6): c++/src/objects/blast/Blast4_parameters.cpp:186-198
    #[cfg(test)]
    pub(crate) fn param(&self, set: ParamSet, name: &[u8]) -> Option<&Param> {
        self.params
            .iter()
            .find(|param| param.set == set && param.name == name)
    }

    /// The size of the hit list (`CBlastOptionsHandle::GetHitlistSize()` of the options that
    /// `CBlastOptionsBuilder::GetSearchOptions` builds): the last `HitlistSize` of the
    /// algorithm and program options, read through [`Request::params`], else NCBI's default,
    /// `BLAST_HITLIST_SIZE` (500, blast_options.h:160, which every options handle sets). The one
    /// option this module reads: the whole option set of a request is line F2's
    /// (`read/options.rs`), which takes this over; until then, only this parameter's own error
    /// is NCBI's (any other parameter that NCBI's builder refuses is not checked here).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:298-306
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:531-545
    // NCBI reference (598d8ae6): c++/src/objects/blast/names.cpp:367-373
    // (`HitlistSize` is an `e_Integer` field, names.cpp; `Which()` is the variant's position from
    // 1. The archive reader builds without `ignore_unsupported_options`. The tests' reader of
    // the one option; the run reads `search_options`.)
    #[cfg(test)]
    pub(crate) fn hitlist_size(&self) -> Result<i32, Stop> {
        self.search_options().map(|options| options.hitlist_size)
    }

    /// Whether the database of the request is protein (`x_SetDatabase`'s residue type).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1537-1546
    pub(crate) fn db_is_protein(&self) -> bool {
        let program = self.program.as_slice();
        program == b"blastp"
            || program == b"blastx"
            || (program == b"tblastn" && self.service == b"rpsblast")
    }

    /// Whether the program translates the query (`Blast_QueryIsTranslated` of
    /// `NetworkProgram2BlastProgramType(program, service)`), which decides the frames of the
    /// masks; NCBI aborts on the `rpsblast` service with a program other than `blastp`,
    /// `tblastn` and `blastx`.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:1933-1962
    // (`BlastProgram2Number`, blast_util.c:278-310, compares the names without case; the
    // programs whose type has `TRANSLATED_QUERY_MASK` are blastx, tblastx and rpstblastn,
    // blast_program.h:72-89. An unknown name is `eBlastTypeUndefined`, 0.)
    pub(crate) fn query_is_translated(&self) -> Result<bool, Stop> {
        let program = self.program.as_slice();
        if self.service == b"rpsblast" {
            return match program {
                b"blastp" => Ok(false),
                b"tblastn" | b"blastx" => Ok(true),
                _ => Err(Stop::Abort(
                    "NetworkProgram2BlastProgramType: the rpsblast service with another program",
                )),
            };
        }
        if self.service == b"psi" {
            return Ok(false);
        }
        Ok(["blastx", "tblastx", "rpstblastn"]
            .iter()
            .any(|name| program.eq_ignore_ascii_case(name.as_bytes())))
    }
}

/// What blast_formatter reads of a request's options: what NCBI's archive reader takes from the
/// options handle that `CRemoteBlast::GetSearchOptions` builds for the archive's first object and
/// keeps for every later one (`m_CBOH`, `m_TaxidList`, `m_DbFilteringAlgorithmId`; decision F3-6
/// (2) of line RG). An interim reading of these four values through the raw accessor
/// ([`Request::params`]), kept to what outfmt 8, 9 and 12 need: the whole option set of a request,
/// with NCBI's locality for blast_formatter (`eRemote`, remote_blast.cpp:2018-2027) and the
/// builder's errors of every other parameter, is line F2's (`read/options.rs`), which takes these
/// values over.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SearchOptions {
    /// The program of the options handle ([`Request::search_program`]).
    pub(crate) program: SearchProgram,
    /// `GetHitlistSize()` ([`Request::hitlist_size`]).
    pub(crate) hitlist_size: i32,
    /// `TaxidList` (`GetTaxidList`): blast_formatter opens the database with it.
    pub(crate) taxids: Option<Vec<i32>>,
    /// `DbFilteringAlgorithmId` other than -1 (`GetDbFilteringAlgorithmId`).
    pub(crate) db_filtering_algorithm: Option<i32>,
}

/// The program of a search's options handle (`CBlastOptionsFactory::Create` of the program that
/// `ComputeProgram` and `AdjustProgram` give), as blast_formatter's formatter names it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum SearchProgram {
    /// A program whose type the formatter names (`Blast_ProgramNameFromType`, then
    /// `ProgramNameToEnum`).
    Named(ProgramType),
    /// `EProgram` of a program whose Seq-annot NCBI's formatter does not name as the others
    /// (PHI-BLAST: `ProgramNameToEnum` of its type's name throws; the mapper of magicblast): its
    /// task name (`EProgramToTaskName`). LOSAT's blast_formatter refuses its archives.
    Other(&'static str),
}

/// The program type of a search as blast_formatter's formatter writes it: the label and value of
/// the Seq-annot's "Blast Type" (`EProgramToTaskName(ProgramNameToEnum(m_Program))` and that
/// `EProgram`), whether its subjects are nucleotide (`Blast_SubjectIsNucleotide`) and whether it
/// is `eBlastTypeBlastn` (the check of SAM).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ProgramType {
    pub(crate) task: &'static str,
    pub(crate) number: i32,
    pub(crate) subject_nucleotide: bool,
    pub(crate) blastn: bool,
}

impl ProgramType {
    const fn new(task: &'static str, number: i32, subject_nucleotide: bool) -> Self {
        ProgramType {
            task,
            number,
            subject_nucleotide,
            blastn: number == 1,
        }
    }
}

/// `EProgram` values (blast_types.hpp:56-76) that the option builder can give.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum EProgram {
    Blastn,
    Blastp,
    Blastx,
    Tblastn,
    Tblastx,
    RpsBlast,
    RpsTblastn,
    Megablast,
    DiscMegablast,
    PsiBlast,
    PsiTblastn,
    PhiBlastp,
    PhiBlastn,
    DeltaBlast,
    VecScreen,
    Mapper,
}

/// The place of an exception of NCBI's option builder.
fn builder_exception(line: u32, kind: &'static str, function: &'static str, msg: String) -> Stop {
    Stop::Exception(NcbiException::new(
        ThrowSite {
            file: "src/algo/blast/api/blast_options_builder.cpp",
            line,
            kind,
            function: Some(function),
        },
        msg,
    ))
}

impl EProgram {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:757-811
    // ```c++
    // EProgram ProgramNameToEnum(const std::string& program_name)
    // {
    //     ...
    //     if (NStr::StartsWith(lowercase_program_name, "blastn")) {
    //         return eBlastn;
    //     // -RMH- support new toolkit program
    //     } else if (NStr::StartsWith(lowercase_program_name, "rmblastn")) {
    //         return eBlastn;
    //     } else if (NStr::StartsWith(lowercase_program_name, "blastp")) {
    //         return eBlastp;
    //     } else if (lowercase_program_name == "blastx") {
    //         return eBlastx;
    //     } else if (lowercase_program_name == "tblastn") {
    //         return eTblastn;
    //     } else if (lowercase_program_name == "tblastx") {
    //         return eTblastx;
    //     } else if (lowercase_program_name == "rpsblast") {
    //         return eRPSBlast;
    //     } else if (lowercase_program_name == "rpstblastn") {
    //         return eRPSTblastn;
    //     } else if (lowercase_program_name == "megablast") {
    //         return eMegablast;
    //     } else if (lowercase_program_name == "psiblast") {
    //         return ePSIBlast;
    //     } else if (lowercase_program_name == "psitblastn") {
    //         return ePSITblastn;
    //     } else if (lowercase_program_name == "dc-megablast") {
    //         return eDiscMegablast;
    //     } else if (lowercase_program_name == "deltablast") {
    //         return eDeltaBlast;
    //     } else if (lowercase_program_name == "vecscreen") {
    //         return eVecScreen;
    //     ...
    //     } else if (lowercase_program_name == "mapper") {
    //         return eMapper;
    //     } else if (lowercase_program_name == "mapr2g") {
    //         return eMapper;
    //     } else if (lowercase_program_name == "mapr2r") {
    //         return eMapper;
    //     } else if (lowercase_program_name == "mapg2g") {
    //         return eMapper;
    //     } else {
    //         NCBI_THROW(CBlastException, eNotSupported,
    //                    "Program type '" + program_name + "' not supported");
    //     }
    // }
    // ```
    /// `ProgramNameToEnum` of a lower-case name: `None` is its `Program type '...' not supported`.
    fn from_name(name: &[u8]) -> Option<EProgram> {
        Some(match name {
            _ if name.starts_with(b"blastn") || name.starts_with(b"rmblastn") => EProgram::Blastn,
            _ if name.starts_with(b"blastp") => EProgram::Blastp,
            b"blastx" => EProgram::Blastx,
            b"tblastn" => EProgram::Tblastn,
            b"tblastx" => EProgram::Tblastx,
            b"rpsblast" => EProgram::RpsBlast,
            b"rpstblastn" => EProgram::RpsTblastn,
            b"megablast" => EProgram::Megablast,
            b"psiblast" => EProgram::PsiBlast,
            b"psitblastn" => EProgram::PsiTblastn,
            b"dc-megablast" => EProgram::DiscMegablast,
            b"deltablast" => EProgram::DeltaBlast,
            b"vecscreen" => EProgram::VecScreen,
            b"mapper" | b"mapr2g" | b"mapr2r" | b"mapg2g" => EProgram::Mapper,
            _ => return None,
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:124-200
    // (`CBlastOptionsFactory::Create`: the nucleotide handles of eBlastn, eMegablast,
    // eDiscMegablast and eVecScreen are `eBlastTypeBlastn`; eRPSBlast `CBlastRPSOptionsHandle`,
    // eRPSTblastn `CRPSTBlastnOptionsHandle`, ePSIBlast `CPSIBlastOptionsHandle`, ePSITblastn the
    // same with `SetPSITblastnDefaults`, eDeltaBlast `CDeltaBlastOptionsHandle`, whose type is
    // `eBlastTypePsiBlast` as `EProgramToEBlastProgramType` gives it, blast_aux.cpp:739-741.)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:319-345 (`BlastNumber2Program`:
    // `eBlastTypeBlastn` "blastn", ... `eBlastTypeRpsBlast` "rpsblast", `eBlastTypeRpsTblastn`
    // "rpstblastn", `eBlastTypePsiBlast` "psiblast", `eBlastTypePsiTblastn` "psitblastn"), which
    // `ProgramNameToEnum` maps back (above) and `EProgramToTaskName` names (blast_aux.cpp:676-694);
    // the numbers are `EProgram`'s (blast_types.hpp:56-76) and the subjects' kind the program
    // type's (blast_program.h:72-89: the PSSM types are blastp's and tblastn's with a PSSM query,
    // blastp's and blastx's with a PSSM subject).
    /// The program of the options handle made for this program, as the formatter names it.
    fn search_program(self) -> SearchProgram {
        SearchProgram::Named(match self {
            EProgram::Blastn
            | EProgram::Megablast
            | EProgram::DiscMegablast
            | EProgram::VecScreen => ProgramType::new("blastn", 1, true),
            EProgram::Blastp => ProgramType::new("blastp", 2, false),
            EProgram::Blastx => ProgramType::new("blastx", 3, false),
            EProgram::Tblastn => ProgramType::new("tblastn", 4, true),
            EProgram::Tblastx => ProgramType::new("tblastx", 5, true),
            EProgram::RpsBlast => ProgramType::new("rpsblast", 6, false),
            EProgram::RpsTblastn => ProgramType::new("rpstblastn", 7, false),
            EProgram::PsiBlast | EProgram::DeltaBlast => ProgramType::new("psiblast", 10, false),
            EProgram::PsiTblastn => ProgramType::new("psitblastn", 11, true),
            EProgram::PhiBlastp => return SearchProgram::Other("phiblastp"),
            EProgram::PhiBlastn => return SearchProgram::Other("phiblastn"),
            EProgram::Mapper => return SearchProgram::Other("mapr2g"),
        })
    }
}

/// `CBlast4Field::Get(<field>).Match(param)` of the field `name` whose value is of the variant
/// `variant` (names.cpp:367-373, quoted at [`Request::hitlist_size`]).
fn typed(param: &Param, name: &[u8], variant: &str) -> bool {
    param.name == name && param_choice(param).0 == variant
}

/// The value of a parameter as NCBI's `CBlast4Field::Match` reads it: its variant's name, and
/// `Which()` (the variant's place from 1) for the builder's error.
fn param_choice(param: &Param) -> (&'static str, V<'_>, usize) {
    let value = param.value();
    let which = match value.value {
        Value::Choice(index, _) => usize::from(*index) + 1,
        _ => 0,
    };
    let (name, inner) = value.choice();
    (name, inner, which)
}

impl Request {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:64-138
    // ```c++
    // CBlastOptionsBuilder::ComputeProgram(const string & program,
    //                                      const string & service)
    // {
    //     string p = program;
    //     string s = service;
    //
    //     NStr::ToLower(p);
    //     NStr::ToLower(s);
    //     ...
    //     bool found = false;
    //
    //     if (p == "blastp") {
    //         if (s == "rpsblast") {
    //             p = "rpsblast";
    //             found = true;
    //         } else if (s == "psi") {
    //             p = "psiblast";
    //             found = true;
    //         } else if (s == "phi") {
    //             // phi is just treated as a blastp here
    //             found = true;
    //         } else if (s == "delta_blast") {
    //             p = "deltablast";
    //             found = true;
    //         }
    //         else if (s == "indexed") {
    //             found = true;
    //         }
    //     } else if (p == "blastn") {
    //         if (s == "megablast") {
    //             p = "megablast";
    //             found = true;
    //         }
    //         if (s == "vecscreen") {
    //             p = "vecscreen";
    //             found = true;
    //         }
    //         if ((s == "sra") || (s == "wgs") || (s=="indexed")) {
    //         	found = true;
    //         }
    //     } else if (p == "tblastn") {
    //         if (s == "rpsblast") {
    //             p = "rpstblastn";
    //             found = true;
    //         } else if (s == "psi") {
    //             p = "psitblastn";
    //             found = true;
    //         }
    //         if ((s == "sra") || (s == "wgs")) {
    //         	found = true;
    //         }
    //     } else if (p == "tblastx") {
    //         found = true;
    //     } else if (p == "blastx") {
    //         if (s == "rpsblast") {
    //             p = "rpstblastn";
    //             found = true;
    //         }
    //     }
    //
    //     if (!found && !(s=="plain" || s=="multi_blast")) {
    //         string msg = "Unsupported combination of program (";
    //         msg += program;
    //         msg += ") and service (";
    //         msg += service;
    //         msg += ").";
    //
    //         NCBI_THROW(CBlastException, eInvalidArgument, msg);
    //     }
    //
    //     return ProgramNameToEnum(p);
    // }
    // ```
    /// `ComputeProgram(program, service)`.
    fn compute_program(&self) -> Result<EProgram, Stop> {
        let p = self.program.to_ascii_lowercase();
        let s = self.service.to_ascii_lowercase();
        let (name, found): (&[u8], bool) = match (p.as_slice(), s.as_slice()) {
            (b"blastp", b"rpsblast") => (b"rpsblast", true),
            (b"blastp", b"psi") => (b"psiblast", true),
            (b"blastp", b"phi" | b"indexed") => (b"blastp", true),
            (b"blastp", b"delta_blast") => (b"deltablast", true),
            (b"blastn", b"megablast") => (b"megablast", true),
            (b"blastn", b"vecscreen") => (b"vecscreen", true),
            (b"blastn", b"sra" | b"wgs" | b"indexed") => (b"blastn", true),
            (b"tblastn", b"rpsblast") => (b"rpstblastn", true),
            (b"tblastn", b"psi") => (b"psitblastn", true),
            (b"tblastn", b"sra" | b"wgs") => (b"tblastn", true),
            (b"tblastx", _) => (b"tblastx", true),
            (b"blastx", b"rpsblast") => (b"rpstblastn", true),
            _ => (p.as_slice(), false),
        };
        if !found && s != b"plain" && s != b"multi_blast" {
            return Err(builder_exception(
                134,
                "CBlastException::eInvalidArgument",
                "BLAST::ncbi::blast::CBlastOptionsBuilder::ComputeProgram",
                format!(
                    "Unsupported combination of program ({}) and service ({}).",
                    String::from_utf8_lossy(&self.program),
                    String::from_utf8_lossy(&self.service)
                ),
            ));
        }
        EProgram::from_name(name).ok_or_else(|| {
            Stop::Exception(NcbiException::new(
                ThrowSite {
                    file: "src/algo/blast/api/blast_aux.cpp",
                    line: 808,
                    kind: "CBlastException::eNotSupported",
                    function: Some("BLAST::ncbi::blast::ProgramNameToEnum"),
                },
                format!(
                    "Program type '{}' not supported",
                    String::from_utf8_lossy(name)
                ),
            ))
        })
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:585-646
    // ```c++
    //     // PHI-BLAST pattern trumps all other options
    //     if (find_if(L->begin(), L->end(),
    //                 SBlast4ParamFinder(eBlastOpt_PHIPattern)) != L->end()) {
    //         switch(program) {
    //         case ePHIBlastn:
    //         case eBlastn:
    //             _TRACE("Adjusting program to phiblastn");
    //             return ePHIBlastn;
    //
    //         case ePHIBlastp:
    //         case eBlastp:
    //             _TRACE("Adjusting program to phiblastp");
    //             return ePHIBlastp;
    //
    //         default:
    //             {
    //                 string msg = "Incorrect combination of option (";
    //                 msg += CBlast4Field::GetName(eBlastOpt_PHIPattern);
    //                 msg += ") and program (";
    //                 msg += program_string;
    //                 msg += ")";
    //
    //                 NCBI_THROW(CRemoteBlastException,
    //                            eServiceNotAvailable,
    //                            msg);
    //             }
    //             break;
    //         }
    //     }
    //
    //     ITERATE(TValueList, iter, *L) {
    //         CBlast4_parameter & p = const_cast<CBlast4_parameter&>(**iter);
    //         const CBlast4_value & v = p.GetValue();
    //
    //         if (CBlast4Field::Get(eBlastOpt_MBTemplateLength).Match(p)) {
    //             if (v.GetInteger() != 0) {
    //                 _TRACE("Adjusting program to discontiguous Megablast");
    //                 return eDiscMegablast;
    //             }
    //         } else if (CBlast4Field::Get(eBlastOpt_Web_StepNumber).Match(p) ||
    //                    CBlast4Field::Get(eBlastOpt_Web_RunPsiBlast).Match(p) ||
    //                    CBlast4Field::Get(eBlastOpt_PseudoCount).Match(p) ||
    //                    CBlast4Field::Get(eBlastOpt_IgnoreMsaMaster).Match(p)
    //                    ) {
    //             // FIXME: should we handle DELTA-BLAST here too?
    //             _TRACE("Adjusting program to psiblast");
    //             return ePSIBlast;
    //         }
    //     }
    //
    //     return program;
    // }
    // ```
    // (`L` is the algorithm, program and format options in that order, `s_MergeCBlast4_parameters`;
    // the fields: `PHIPattern` string, `MBTemplateLength`, `Web_StepNumber` and `PseudoCountWeight`
    // integers, `Web_RunPsiBlast` and `IgnoreMsaMaster` booleans, names.cpp.)
    /// `AdjustProgram` of `program` over every parameter of the request.
    fn adjust_program(&self, program: EProgram) -> Result<EProgram, Stop> {
        if self
            .params
            .iter()
            .any(|p| typed(p, b"PHIPattern", "string"))
        {
            return match program {
                EProgram::PhiBlastn | EProgram::Blastn => Ok(EProgram::PhiBlastn),
                EProgram::PhiBlastp | EProgram::Blastp => Ok(EProgram::PhiBlastp),
                _ => Err(builder_exception(
                    616,
                    "CRemoteBlastException::eServiceNotAvailable",
                    "BLAST::ncbi::blast::CBlastOptionsBuilder::AdjustProgram",
                    format!(
                        "Incorrect combination of option (PHIPattern) and program ({})",
                        String::from_utf8_lossy(&self.program)
                    ),
                )),
            };
        }
        for param in self.params() {
            if typed(param, b"MBTemplateLength", "integer") {
                if param_choice(param).1.int() != 0 {
                    return Ok(EProgram::DiscMegablast);
                }
            } else if typed(param, b"Web_StepNumber", "integer")
                || typed(param, b"Web_RunPsiBlast", "boolean")
                || typed(param, b"PseudoCountWeight", "integer")
                || typed(param, b"IgnoreMsaMaster", "boolean")
            {
                return Ok(EProgram::PsiBlast);
            }
        }
        Ok(program)
    }

    /// The program of the request's options handle (`GetSearchOptions`: `ComputeProgram`, then
    /// `AdjustProgram`), or NCBI's exception there.
    pub(crate) fn search_program(&self) -> Result<SearchProgram, Stop> {
        let program = self.compute_program()?;
        Ok(self.adjust_program(program)?.search_program())
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:205-211
    // ```c++
    //         } else if (CBlast4Field::Get(eBlastOpt_DbFilteringAlgorithmId).Match(p)) {
    //             m_DbFilteringAlgorithmId = v.GetInteger();
    //         } else if (CBlast4Field::Get(eBlastOpt_DbFilteringAlgorithmKey).Match(p)) {
    //             m_DbFilteringAlgorithmKey = v.GetString();
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_builder.cpp:481-495
    // ```c++
    //     case 'T':
    //     	if (CBlast4Field::Get(eBlastOpt_TaxidList).Match(p)) {
    //     ...
    //             m_TaxidList = v.GetInteger_list();
    //     ...
    //     	} else {
    //     		found = false;
    //     	}
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2045-2060
    // ```c++
    //         if(bob.HaveTaxidList()){
    //         	list<TTaxId> list = bob.GetTaxidList();
    //         	m_TaxidList.clear();
    //         	m_TaxidList.insert(list.begin(), list.end());
    //         }
    //         ...
    //         if (bob.HasDbFilteringAlgorithmId() &&
    //             bob.GetDbFilteringAlgorithmId() != -1) {
    //             m_DbFilteringAlgorithmId = bob.GetDbFilteringAlgorithmId();
    //         }
    // ```
    // (`HitlistSize`: [`Request::hitlist_size`]; a parameter of these names with another value
    // type is the builder's `Internal: Error processing option` there, in the order of the
    // parameters. `NegativeTaxidList` is read into a list that blast_formatter does not use:
    // blast_formatter.cpp:378-382 builds its negative list from `GetTaxidList()`, which is empty
    // without `TaxidList`, and an empty list opens the database without one, uniform_search.cpp:
    // 254-291.)
    /// The options blast_formatter reads ([`SearchOptions`]), or NCBI's exception of the program
    /// (`search_program`) or of the first of these options whose value has another type.
    pub(crate) fn search_options(&self) -> Result<SearchOptions, Stop> {
        let program = self.search_program()?;
        let mut options = SearchOptions {
            program,
            hitlist_size: 500,
            taxids: None,
            db_filtering_algorithm: None,
        };
        let mut filtering = None;
        for param in self.params() {
            if param.set == ParamSet::Format {
                continue;
            }
            let (variant, value, which) = param_choice(param);
            let wanted = match param.name.as_slice() {
                b"HitlistSize" | b"DbFilteringAlgorithmId" => "integer",
                b"TaxidList" => "integer-list",
                _ => continue,
            };
            if variant != wanted {
                return Err(option_type_error(&param.name, which));
            }
            match param.name.as_slice() {
                b"HitlistSize" => options.hitlist_size = value.int() as i32,
                b"DbFilteringAlgorithmId" => filtering = Some(value.int() as i32),
                _ => options.taxids = Some(value.ints().map(|id| id as i32).collect()),
            }
        }
        options.db_filtering_algorithm = filtering.filter(|&id| id != -1);
        Ok(options)
    }
}

/// The option builder's exception for the parameter `name` whose value is of the variant
/// `which` (`Which()`), quoted at [`Request::hitlist_size`].
fn option_type_error(name: &[u8], which: usize) -> Stop {
    builder_exception(
        541,
        "CRemoteBlastException::eServiceNotAvailable",
        "BLAST::ncbi::blast::CBlastOptionsBuilder::x_ProcessOneOption",
        format!(
            "Internal: Error processing option [{}] type [{which}].",
            String::from_utf8_lossy(name)
        ),
    )
}

impl Queries {
    /// `CBlast4_queries::GetNumQueries()`: 1 for a PSSM, the locations, the entries at the top
    /// of the Bioseq-set.
    // NCBI reference (598d8ae6): c++/src/objects/blast/Blast4_queries.cpp:57-77
    pub(crate) fn num_queries(&self) -> usize {
        match self {
            Queries::Pssm(_) => 1,
            Queries::SeqLocList(locs) => locs.len(),
            Queries::BioseqSet(set) => V::new(named::BIOSEQ_SET, set).field("seq-set").len(),
        }
    }

    /// The PSSM of a `pssm` query.
    pub(crate) fn pssm(&self) -> Option<V<'_>> {
        match self {
            Queries::Pssm(pssm) => Some(V::new(named::PSSMWITHPARAMETERS, pssm)),
            _ => None,
        }
    }

    /// The locations of a `seq-loc-list` query.
    pub(crate) fn seq_locs(&self) -> impl Iterator<Item = V<'_>> {
        let locs: &[Value] = match self {
            Queries::SeqLocList(locs) => locs,
            _ => &[],
        };
        locs.iter().map(|loc| V::new(named::SEQ_LOC, loc))
    }

    /// The Bioseqs of the queries in the order of NCBI's iterator (`CTypeConstIterator<CBioseq>`
    /// over the Bioseq-set, or over the PSSM's query when it is a set): every Bioseq of the set's
    /// entries, depth first. None for a `seq-loc-list`; a PSSM without a query has none.
    pub(crate) fn bioseqs(&self) -> Vec<Bioseq<'_>> {
        let mut out = Vec::new();
        match self {
            Queries::BioseqSet(set) => set_bioseqs(V::new(named::BIOSEQ_SET, set), &mut out),
            Queries::Pssm(pssm) => {
                if let Some(entry) = V::new(named::PSSMWITHPARAMETERS, pssm)
                    .field("pssm")
                    .get("query")
                {
                    entry_bioseqs(entry, &mut out);
                }
            }
            Queries::SeqLocList(_) => {}
        }
        out
    }

    /// `x_ExtractQueryIds`: the id of each query, in order (`None`: NCBI's null `CRef`, of a
    /// Bioseq whose ids are all ranked worst or a location without one id).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2089-2117
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/remote_blast.cpp:2078-2087
    // (A PSSM without a query: `GetQuery()` of an unset member.)
    pub(crate) fn query_ids(&self) -> Result<Vec<Option<&ParsedSeqId>>, Stop> {
        if let Queries::Pssm(pssm) = self {
            let query = V::new(named::PSSMWITHPARAMETERS, pssm)
                .field("pssm")
                .get("query");
            if query.is_none() {
                return Err(Stop::Unassigned {
                    type_name: "Pssm",
                    member: "query",
                });
            }
        }
        if let Queries::SeqLocList(_) = self {
            return Ok(self.seq_locs().map(super::seq_loc_id).collect());
        }
        Ok(self
            .bioseqs()
            .into_iter()
            .map(|bioseq| bioseq.best_id())
            .collect())
    }
}

impl Subject {
    /// The Bioseqs of a `sequences` subject.
    pub(crate) fn bioseqs(&self) -> Vec<Bioseq<'_>> {
        match self {
            Subject::Sequences(bioseqs) => bioseqs
                .iter()
                .map(|bioseq| Bioseq(V::new(named::BIOSEQ, bioseq)))
                .collect(),
            _ => Vec::new(),
        }
    }
}

/// The Bioseqs of a Bioseq-set's entries, depth first (`CTypeConstIterator<CBioseq>`: in the
/// types of the table, a Bioseq is only the `seq` of a `Seq-entry`, and the entries of a set
/// are its `seq-set`).
fn set_bioseqs<'a>(set: V<'a>, out: &mut Vec<Bioseq<'a>>) {
    let mut todo: Vec<V<'a>> = set.field("seq-set").items().collect();
    todo.reverse();
    while let Some(entry) = todo.pop() {
        match entry.choice() {
            ("seq", bioseq) => out.push(Bioseq(bioseq)),
            (_, set) => {
                let at = todo.len();
                todo.extend(set.field("seq-set").items());
                todo[at..].reverse();
            }
        }
    }
}

/// The Bioseqs of a Seq-entry (`entry_bioseqs` of its set, or its Bioseq).
pub(super) fn entry_bioseqs<'a>(entry: V<'a>, out: &mut Vec<Bioseq<'a>>) {
    match entry.choice() {
        ("seq", bioseq) => out.push(Bioseq(bioseq)),
        (_, set) => set_bioseqs(set, out),
    }
}

/// A Bioseq of a request, as read.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Bioseq<'a>(V<'a>);

impl<'a> Bioseq<'a> {
    /// The Bioseq of a value of type `Bioseq` (a subject of a request).
    pub(super) fn of(value: V<'a>) -> Self {
        Bioseq(value)
    }

    /// The Bioseq value.
    pub(crate) fn value(self) -> V<'a> {
        self.0
    }

    /// Its ids (`GetId()`), in order.
    pub(crate) fn ids(self) -> impl Iterator<Item = &'a ParsedSeqId> {
        self.0.field("id").items().map(|id| id.seq_id())
    }

    /// `FindBestChoice(GetId(), CSeq_id::BestRank)`: the first id of the best rank, `None` when
    /// every id ranks `kMax_Int` or there is none.
    pub(crate) fn best_id(self) -> Option<&'a ParsedSeqId> {
        find_best(self.ids(), best_rank_score)
    }

    /// The descriptors (`descr`) in order.
    pub(super) fn descriptors(self) -> impl Iterator<Item = (&'static str, V<'a>)> {
        self.0
            .get("descr")
            .into_iter()
            .flat_map(|descr| descr.items())
            .map(|desc| desc.choice())
    }

    /// The first `title` descriptor, as stored.
    pub(crate) fn title(self) -> Option<&'a [u8]> {
        self.descriptors()
            .find(|(variant, _)| *variant == "title")
            .map(|(_, title)| title.str())
    }

    /// The `DefLine` field of the first `CFastaReader` user object, as stored (the raw defline
    /// line of the FASTA reader, `NStr::PrintableString` of it).
    #[cfg(test)] // read by the formats of later pull requests (decision C3-3)
    pub(crate) fn fasta_defline(self) -> Option<&'a [u8]> {
        self.descriptors()
            .filter(|(variant, _)| *variant == "user")
            .map(|(_, user)| user)
            .find(|user| matches!(user.field("type").choice(), ("str", kind) if kind.str() == b"CFastaReader"))
            .and_then(|user| {
                user.field("data").items().find_map(|field| {
                    match (field.field("label").choice(), field.field("data").choice()) {
                        (("str", label), ("str", text)) if label.str() == b"DefLine" => {
                            Some(text.str())
                        }
                        _ => None,
                    }
                })
            })
    }

    /// `inst.mol` (`Seq-inst.mol`, ENUMERATED: 0 not-set, 1 dna, 2 rna, 3 aa, 4 na, 255 other).
    pub(crate) fn molecule(self) -> i64 {
        self.0.field("inst").field("mol").int()
    }

    /// `inst.length`, when set.
    pub(crate) fn length(self) -> Option<i64> {
        self.0
            .field("inst")
            .get("length")
            .map(|length| length.int())
    }

    /// `inst.seq-data` as read (a `Seq-data` CHOICE of any coding).
    #[cfg(test)] // read by the formats of later pull requests (decision C3-3)
    pub(crate) fn seq_data(self) -> Option<V<'a>> {
        self.0.field("inst").get("seq-data")
    }
}
