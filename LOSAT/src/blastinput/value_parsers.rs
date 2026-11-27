//! Shared CLI value grammar; engine configuration remains typed.
use crate::utils::seg::SegParams;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:410-420
#[derive(Debug, Clone, PartialEq)]
pub enum DustSpec {
    No,
    Yes,
    Parameters {
        level: u32,
        window: usize,
        linker: usize,
    },
}

impl DustSpec {
    pub fn is_enabled(&self) -> bool {
        !matches!(self, Self::No)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:46-48,62-64:
    pub fn params(&self) -> Option<(u32, usize, usize)> {
        match *self {
            Self::No => None,
            Self::Yes => Some((20, 64, 1)),
            Self::Parameters {
                level,
                window,
                linker,
            } => Some((level, window, linker)),
        }
    }
}

/// A `-dust` value as NCBI reads it: `no`, `yes`, or three signed decimal numbers separated
/// by single spaces. The error is NCBI's message; the masker replaces values out of its
/// ranges with its defaults (`CSymDustMasker`, `utils/dust.rs`), after NCBI's conversion
/// to `Uint4` (`as u32`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:375-384
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:410-426
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:195-198
/// `NStr::Split` keeps the empty fields between adjacent spaces; `NStr::StringToInt` (flags 0)
/// is `ncbi_str::string_to_int`.
pub fn parse_dust_filtering(value: impl AsRef<[u8]>) -> Result<DustSpec, String> {
    let value = value.as_ref();
    match value {
        b"no" => return Ok(DustSpec::No),
        b"yes" => return Ok(DustSpec::Yes),
        _ => {}
    }
    let tokens: Vec<&[u8]> = value.split(|&b| b == b' ').collect();
    if tokens.len() != 3 {
        return Err("Invalid number of arguments to filtering option".into());
    }
    // `NStr::StringToInt` of each token, flags 0.
    let number = |token: &[u8]| {
        crate::ncbi_str::string_to_int(token, 0)
            .map_err(|_| "Invalid input for filtering parameters".to_string())
    };
    let (level, window, linker) = (number(tokens[0])?, number(tokens[1])?, number(tokens[2])?);
    Ok(DustSpec::Parameters {
        level: level as u32,
        window: window as u32 as usize,
        linker: linker as u32 as usize,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:386-407
#[derive(Debug, Clone, PartialEq)]
pub enum SegSpec {
    No,
    Yes,
    WindowLocutHicut { window: i32, locut: f64, hicut: f64 },
}

impl SegSpec {
    #[inline]
    pub fn is_enabled(&self) -> bool {
        !matches!(self, Self::No)
    }

    /// The options as NCBI's filter string prints them (`BlastFilteringOptionsToString`): the
    /// defaults for `yes`, the given window, locut and hicut otherwise (the substitution of
    /// the defaults for a value not above 0 is `params`', after the options are set).
    pub fn filter_options(&self) -> Option<(i32, f64, f64)> {
        match *self {
            Self::No => None,
            Self::Yes => Some((12, 2.2, 2.5)),
            Self::WindowLocutHicut {
                window,
                locut,
                hicut,
            } => Some((window, locut, hicut)),
        }
    }

    /// The SEG parameters of the query filter: a window, locut or hicut that is not above
    /// 0 keeps NCBI's default (`SegParametersNewAa`), then SEG checks the parameters
    /// (`SegParams::new`, blast_seg.c:2247-2258).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1147-1154
    #[inline]
    pub fn params(&self) -> Option<SegParams> {
        match *self {
            Self::No => None,
            Self::Yes => Some(SegParams::default()),
            Self::WindowLocutHicut {
                window,
                locut,
                hicut,
            } => {
                let defaults = SegParams::default();
                Some(SegParams::new(
                    usize::try_from(window)
                        .ok()
                        .filter(|&window| window > 0)
                        .unwrap_or(defaults.window),
                    if locut > 0.0 { locut } else { defaults.locut },
                    if hicut > 0.0 { hicut } else { defaults.hicut },
                ))
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:396-407
    pub fn to_ncbi_cli_string(&self) -> String {
        match self {
            Self::No => "no".to_string(),
            Self::Yes => "yes".to_string(),
            Self::WindowLocutHicut {
                window,
                locut,
                hicut,
            } => format!("{window} {locut} {hicut}"),
        }
    }
}

/// A `-seg` value as NCBI reads it: `no`, `yes`, or a window (an `int`), a locut and a
/// hicut separated by single spaces (`x_TokenizeFilteringArgs`, see `parse_dust_filtering`).
/// The errors are NCBI's messages; LOSAT also rejects a locut or hicut that is not finite.
/// The values are used as NCBI's query filter uses them (`SegSpec::params`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:396-406
#[cfg(any(losat_web_abi_v1, test))]
pub fn parse_seg_filtering(value: &str) -> Result<SegSpec, String> {
    match value {
        "no" => return Ok(SegSpec::No),
        "yes" => return Ok(SegSpec::Yes),
        _ => {}
    }
    let tokens: Vec<&str> = value.split(' ').collect();
    if tokens.len() != 3 {
        return Err("Invalid number of arguments to filtering option".into());
    }
    let invalid = || "Invalid input for filtering parameters".to_string();
    let window = tokens[0].parse::<i32>().map_err(|_| invalid())?;
    let locut = tokens[1].parse::<f64>().map_err(|_| invalid())?;
    let hicut = tokens[2].parse::<f64>().map_err(|_| invalid())?;
    if !locut.is_finite() || !hicut.is_finite() {
        return Err("a SEG locut or hicut that is not finite is not supported by LOSAT".into());
    }
    Ok(SegSpec::WindowLocutHicut {
        window,
        locut,
        hicut,
    })
}

/// The TBLASTX `-outfmt` value of web ABI v1, which is frozen (plan TD-1): it keeps the
/// only format that it accepted before TBLASTX outfmt 0 and 7 were ported (session S08).
#[cfg(losat_web_abi_v1)]
pub fn tblastx_v1_outfmt(value: &str) -> Result<String, String> {
    if value.trim() != "6" {
        return Err(
            "unsupported TBLASTX outfmt: only 6 without custom fields is implemented".into(),
        );
    }
    Ok(value.into())
}

#[cfg(test)]
#[path = "value_parsers_tests.rs"]
mod tests;

/// The `-task` of a program: NCBI's reader allows the names of its set (`CArgAllowStringSet`,
/// `c++/src/algo/blast/blastinput/blast_args.cpp:105`; `ncbi_app::args::choices`); `check` checks
/// a task that a caller of the library set itself.
#[derive(Clone)]
pub struct TaskParser {
    /// The program's name in LOSAT's messages (`BLASTN`).
    program: &'static str,
    /// The program whose set of tasks NCBI's argument tables hold (`blastn`).
    name: &'static str,
}

/// The `-task` of `program` (its name in LOSAT's messages, `BLASTN`; `name` in NCBI's argument
/// tables, `blastn`), `TaskParser`.
pub const fn task_parser(program: &'static str, name: &'static str) -> TaskParser {
    TaskParser { program, name }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:211-255
/// The tasks of NCBI's options factory (`GetTasks` with `eAll`).
const NCBI_TASKS: &[&str] = &[
    "blastn",
    "blastn-short",
    "megablast",
    "dc-megablast",
    "vecscreen",
    "rmblastn",
    "blastp",
    "blastp-short",
    "blastp-fast",
    "psiblast",
    "phiblastp",
    "rpsblast",
    "rpstblastn",
    "blastx",
    "blastx-fast",
    "deltablast",
    "tblastn",
    "tblastn-fast",
    "psitblastn",
    "tblastx",
    "kblastp",
    "mapper",
    "mapr2g",
    "mapr2r",
    "mapg2g",
];

impl TaskParser {
    /// The check of a task that a caller of `run_local` set itself (the command line's parser
    /// makes it at parsing), as NCBI's options factory makes it: the name in lower case (the
    /// engine compares it exactly), then a task of NCBI's set of the program; LOSAT's refusal
    /// for another task of NCBI's factory (another program's, which NCBI's API would set up);
    /// NCBI's error for a name NCBI does not know.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_handle.cpp:336-342
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:662-673
    pub fn check(&self, task: &mut str) -> anyhow::Result<()> {
        task.make_ascii_lowercase();
        let task = &*task;
        let accepted = crate::ncbi_app::args::choices(self.name, "task")
            .expect("a program with -task has NCBI's set of tasks");
        if accepted.iter().any(|name| name == task) {
            return Ok(());
        }
        if NCBI_TASKS.contains(&task) {
            let (last, rest) = accepted.split_last().expect("a program has a task");
            let tasks = if rest.is_empty() {
                last.to_string()
            } else {
                format!("{} or {last}", rest.join(", "))
            };
            anyhow::bail!(
                "the task {task} is not supported by LOSAT's {} (use {tasks})",
                self.program
            );
        }
        anyhow::bail!("'{task}' is not a supported task")
    }
}
