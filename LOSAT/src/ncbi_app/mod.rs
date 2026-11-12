//! NCBI's application (`CNcbiApplication::AppMain`, c++/src/corelib/ncbiapp.cpp) for the programs
//! that LOSAT reads with NCBI's argument reader rather than clap's: `losat blastdbcmd` and `losat
//! dustmasker`. Each program is a `Tool` (its name, version, frozen `-help` text and argument
//! descriptions) and a body, its `Run`; `run_with` does what NCBI's application does around it,
//! in NCBI's order: the diagnostics of `-logfile`, the words taken before the reader (`-version`,
//! `-conffile`, `-dryrun`), the configuration load, the arguments (`args`, NCBI's
//! `CArgDescriptions`), the usage and the report of an argument error, and LOSAT's refusals of
//! the NCBI C++ Toolkit options that it does not serve.

pub(crate) mod args;
pub(crate) mod files;

use std::ffi::OsString;
use std::io::{Read, Write};

use crate::blastinput::ncbi_environment::{
    application_start_of_process, check_ncbi_application_settings, diag_target, log_file, preparse,
    ApplicationSettings, DiagTarget, Preparsed,
};

use args::{ArgException, Args, ArgsError, Dep, Desc};

/// A program run by NCBI's application and argument reader.
pub(crate) struct Tool {
    /// The program's name (`GetProgramBasename`): the start of its `-version` line, its usage,
    /// its `<program>.ini` and its settings (`check_ncbi_application_settings`).
    pub(crate) program: &'static str,
    /// The program's name in LOSAT's refusals (`BLASTDBCMD`).
    pub(crate) refusal_name: &'static str,
    /// The version that `-version` prints (`GetFullVersion`).
    pub(crate) version: &'static str,
    /// NCBI's `-help` text (`CArgDescriptions::PrintUsage` of the descriptions with the details),
    /// frozen. Its part before `OPTIONAL ARGUMENTS` is the usage that `-h` and the argument
    /// errors print (`short_usage`).
    pub(crate) help: &'static str,
    /// The program's arguments in the order its `Init` describes them; the reader adds those
    /// of every `CArgDescriptions` (`args::descs`).
    pub(crate) descs: &'static [Desc],
    /// The `SetDependency` calls of its `Init`, in their order (`args::dependencies` makes
    /// NCBI's `multimap` of them).
    pub(crate) dependencies: fn() -> Vec<(&'static str, Dep, &'static str)>,
    /// The prefix of the program's diagnostics (`SetDiagPostPrefix`), if it sets one.
    pub(crate) post_prefix: Option<&'static str>,
}

impl Tool {
    /// LOSAT's refusal of the option `word`.
    pub(crate) fn unsupported_option(&self, word: &str) -> String {
        format!(
            "the NCBI BLAST+ option {word} is not supported by LOSAT's {}",
            self.refusal_name
        )
    }
}

/// `Tool::dependencies` of a program that sets none.
pub(crate) fn no_dependencies() -> Vec<(&'static str, Dep, &'static str)> {
    Vec::new()
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:3829-3834
// ```c++
//     if (detailed) {
//         x.AddDetails(arr);
//     } else {
//         arr.push_back(kEmptyStr);
//         arr.push_back("Use '-help' to print detailed descriptions of command line arguments");
//     }
// ```
/// The text of `-h`: the usage and the description of `help` (a `Tool::help`), then NCBI's
/// pointer to `-help`.
pub(crate) fn short_usage(help: &str) -> String {
    let end = help.find("OPTIONAL ARGUMENTS").unwrap_or(help.len());
    format!(
        "{}Use '-help' to print detailed descriptions of command line arguments\n",
        &help[..end]
    )
}

/// Where a program's diagnostics go: the target of `-logfile` (`stderr` for `cerr`), and the
/// standard error stream, on which the usage of an argument error is printed.
pub(crate) struct Diag<'a> {
    target: DiagTarget,
    stderr: &'a mut dyn Write,
    tool: &'static Tool,
}

impl Diag<'_> {
    /// Posts one diagnostic line (`LOG_POST`, `ERR_POST`).
    pub(crate) fn post(&mut self, line: &str) {
        self.post_bytes(line.as_bytes());
    }

    /// `post` of a line that is not UTF-8.
    pub(crate) fn post_bytes(&mut self, line: &[u8]) {
        let mut text = line.to_vec();
        text.push(b'\n');
        self.write_text(&text);
    }

    /// Writes `text`, diagnostic lines that end with their line feeds (an NCBI report made by
    /// `crate::ncbi_exception`: `err_post_x`, `NcbiException::report_default_ex`).
    pub(crate) fn write_text(&mut self, text: &[u8]) {
        self.target.write_text(text, self.stderr);
    }

    /// LOSAT's refusal: a diagnostic line, exit status 1.
    pub(crate) fn refuse(&mut self, reason: &str) -> i32 {
        self.post(&format!("Error: {reason}"));
        1
    }

    // NCBI reference (598d8ae6): c++/src/app/blastdb/blastdbcmd.cpp:1175-1177
    // ```c++
    //     SetDiagPostLevel(eDiag_Warning);
    //     SetDiagPostPrefix("blastdbcmd");
    // ```
    /// `ERR_POST(Error << message)` inside `Run`: the severity, the prefix of the application
    /// (`Tool::post_prefix`), the message.
    pub(crate) fn err_post(&mut self, message: &str) {
        self.err_post_bytes(message.as_bytes());
    }

    /// `err_post` of a message that is not UTF-8.
    pub(crate) fn err_post_bytes(&mut self, message: &[u8]) {
        let mut line = b"Error: ".to_vec();
        if let Some(prefix) = self.tool.post_prefix {
            line.extend_from_slice(format!("[{prefix}] ").as_bytes());
        }
        line.extend_from_slice(message);
        self.post_bytes(&line);
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:704-710
    // ```c++
    //         if ( s_HandleExceptions() ) {
    //             try {
    //                 *exit_code = m_DryRun ? DryRun() : Run();
    //             }
    //             catch (CArgException& e) {
    //                 NCBI_RETHROW_SAME(e, "Application's execution failed");
    //             }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:274-285
    // ```c++
    // void CException::AddBacklog(const CDiagCompileInfo& info,
    //                             const string& message,
    //                             EDiagSev severity)
    // {
    //     const CException* prev = m_Predecessor;
    //     m_Predecessor = x_Clone();
    //     ...
    //     x_Init(info, message, 0, severity);
    //     m_MainText = false;
    // }
    // ```
    /// A `CArgException` that the program's `Run` throws (`files::open_input_file`,
    /// `Args::output_file`): the application adds `Application's execution failed` to it and
    /// reports it as an argument error, after the usage; exit status 1.
    pub(crate) fn run_arg_exception(&mut self, error: ArgException) -> i32 {
        self.report_arg_exception(&error, Some("Application's execution failed"))
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:1086-1112
    // ```c++
    //     catch (const CArgException& e) {
    //         // Print USAGE and the exception error message
    //         if ( e.GetErrCode() != CArgException::eNoValue &&  m_ArgDesc.get() ) {
    //             x_AddDefaultArgs();
    //             string str;
    //             if ( !m_ArgDesc->IsSetMiscFlag(CArgDescriptions::fNoUsage) ) {
    //                 m_ArgDesc->PrintUsage(str);
    //                 cerr << str;
    //             }
    //             if ( !m_ArgDesc->IsSetMiscFlag(CArgDescriptions::fDupErrToCerr) ) {
    //                 CStreamDiagHandler* errh =
    //                     dynamic_cast<CStreamDiagHandler*>(GetDiagHandler());
    //                 if (!errh  ||  errh->GetStream() != &cerr) {
    //                     cerr << "Error in command-line arguments. "
    //                         "See error logs for more details." << endl;
    //                 }
    //             }
    //             ...
    //             cerr << string(72, '=') << endl << endl;
    //         }
    //         SetDiagPostAllFlags(eDPF_Severity);
    //         NCBI_REPORT_EXCEPTION_X(18, "", e);
    //         got_exception = true;
    //         exit_code = 1;
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7877-7888
    // ```c++
    //     for (pex = &ex; pex; pex = pex->GetPredecessor()) {
    //         pile.push(pex);
    //         if (!main_pex  &&  pex->HasMainText())
    //             main_pex = pex;
    //     }
    //     if (!main_pex)
    //         main_pex = pile.top();
    //     if ( !IsOssEmpty(*m_Buffer.m_Stream) ) {
    //         *this << "(" << main_pex->GetType() << "::"
    //                      << main_pex->GetErrCodeString() << ") "
    //               << s_GetExceptionText(main_pex);
    //     }
    // ```
    // With `eDPF_Severity` alone each exception of the chain is posted as `Error: <text>`, oldest
    // first, then the title's line, `Error: ` and the empty title, ` (CArgException::<code>) `
    // and the text of the newest exception that has its main text: the one the reader threw (a
    // backlog has none).
    /// NCBI's report of the argument error `error` (the backlog `backlog` added, when `Run`
    /// threw it): the usage, the line of `=` and the posts; exit status 1.
    fn report_arg_exception(&mut self, error: &ArgException, backlog: Option<&str>) -> i32 {
        if error.code != "eNoValue" {
            let _ = self
                .stderr
                .write_all(short_usage(self.tool.help).as_bytes());
            // The diagnostics of a `-logfile` go to a stream that is not `cerr`: NCBI then says
            // where they are.
            if !matches!(self.target, DiagTarget::Stderr) {
                let _ = writeln!(
                    self.stderr,
                    "Error in command-line arguments. See error logs for more details."
                );
            }
            let _ = write!(self.stderr, "{}\n\n", "=".repeat(72));
            let _ = self.stderr.flush();
        }
        for message in error.messages.iter().map(String::as_str).chain(backlog) {
            self.post(&format!("Error: {message}"));
        }
        self.post(&format!(
            "Error:  (CArgException::{}) {}",
            error.code,
            error.message()
        ));
        1
    }
}

/// Runs `tool` with `words` (the command line after `losat <program>`) and its body `body` on
/// the process's standard streams; the exit status.
pub(crate) fn run(tool: &'static Tool, words: &[OsString], body: Body) -> i32 {
    let stdout = std::io::stdout();
    let stderr = std::io::stderr();
    let mut out = stdout.lock();
    let mut err = stderr.lock();
    let status = run_with(tool, words, &mut std::io::stdin(), &mut out, &mut err, body);
    let _ = out.flush();
    let _ = err.flush();
    status
}

/// A program's `Run`: it gets the arguments, the application's settings, the standard input
/// (`files::Cin`, what an `eInputFile` argument `-` reads: `files::open_input_file`), the
/// standard output and the diagnostics, and returns the exit status.
pub(crate) type Body = fn(
    args: &Args,
    settings: &ApplicationSettings,
    input: &mut files::Cin<'_>,
    out: &mut dyn Write,
    diag: &mut Diag,
) -> i32;

// NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:413-424 (`x_TryInit` loads the
// configuration before the arguments are read; quoted at `ApplicationStart`)
// ```c++
// void CNcbiApplicationAPI::x_TryInit(EAppDiagStream diag, const char* conf)
// {
//     // Load registry from the config file
//     if ( conf ) {
//         string x_conf(conf);
//         LoadConfig(*m_Config, &x_conf);
// ```
/// The run of `tool` with `words`, its standard input `input`, its output on `out` and its
/// errors on `err`; the exit status.
///
/// The order is NCBI's: the diagnostic file of `-logfile` is opened first (append; stderr when
/// it cannot be opened), then `-version` is read, then the configuration is loaded (its
/// messages go where the diagnostics go; a `<program>.ini` that NCBI reports as malformed ends
/// the run with exit status 2: g_reg probes `o2`, `inibad_null`), then the arguments; an
/// argument error prints the usage and a line of `=` on stderr and NCBI's report of the
/// exception as diagnostics. LOSAT's refusals come where NCBI would start the work they stand
/// for, before anything is written: `-conffile` without the command line's application start
/// (which loads its file) and the full versions where NCBI reads them,
/// `-xmlhelp` where the reader takes it, and `-dryrun` without that start; with it, `-dryrun`
/// ends the run with exit status 0 after the arguments and the settings. Then `body` runs.
pub(crate) fn run_with(
    tool: &'static Tool,
    words: &[OsString],
    input: &mut dyn Read,
    out: &mut dyn Write,
    err: &mut dyn Write,
    body: Body,
) -> i32 {
    let Some(words) = words
        .iter()
        .map(|word| word.to_str().map(str::to_string))
        .collect::<Option<Vec<String>>>()
    else {
        let _ = writeln!(
            err,
            "Error: a word of the command line is not UTF-8; NCBI BLAST+ reads its bytes, which \
             is not supported by LOSAT's {}",
            tool.refusal_name
        );
        return 1;
    };
    let mut diag = Diag {
        target: diag_target(log_file(&words)),
        stderr: err,
        tool,
    };
    // NCBI's application start (`ncbi_environment::ApplicationStart`), which the command line's
    // entry has made: the file of `-logfile` gets the usage report's messages when NCBI opens it,
    // before `-version`.
    let start = application_start_of_process();
    if let Some(start) = start {
        start.write_log_setup(&mut diag.target);
    }
    let (words, dryrun) = match preparse(&words) {
        Preparsed::Version => {
            let _ = write!(
                out,
                "{}: {}\n Package: LOSAT {}, output compatible with blast {}\n",
                tool.program,
                tool.version,
                env!("CARGO_PKG_VERSION"),
                crate::NCBI_BLAST_RELEASE
            );
            return 0;
        }
        Preparsed::Refused(word) => return diag.refuse(&tool.unsupported_option(word)),
        // The file of `-conffile` is the configuration of the command line's application start
        // (`start_application`); the library API runs none (`cli::toolkit_file_word_refusal`).
        Preparsed::Words {
            conffile: Some(_), ..
        } if start.is_none() => {
            return diag.refuse(&crate::cli::toolkit_file_word_refusal(
                tool.program,
                "-conffile",
            ))
        }
        // So is `-dryrun`, which the command line's entry ends after the arguments.
        Preparsed::Words { dryrun: true, .. } if start.is_none() => {
            return diag.refuse(&crate::cli::toolkit_file_word_refusal(
                tool.program,
                "-dryrun",
            ))
        }
        Preparsed::Words { words, dryrun, .. } => (words, dryrun),
    };
    if let Some(start) = start {
        match start.write_config_load(&mut diag.target, diag.stderr) {
            Ok(None) => {}
            Ok(Some(status)) => return status,
            Err(reason) => return diag.refuse(&reason),
        }
    }
    let args = match args::parse(tool, &words) {
        Ok(args) => args,
        Err(ArgsError::Usage) => {
            let _ = out.write_all(short_usage(tool.help).as_bytes());
            return 0;
        }
        Err(ArgsError::Help) => {
            let _ = out.write_all(tool.help.as_bytes());
            return 0;
        }
        // `eNoValue` (which prints no usage) is not raised by the reader.
        Err(ArgsError::Invalid(error)) => return diag.report_arg_exception(&error, None),
        Err(ArgsError::Refused(reason)) => return diag.refuse(&reason),
    };
    let settings = match check_ncbi_application_settings(tool.program) {
        Ok(settings) => settings,
        Err(reason) => return diag.refuse(&reason),
    };
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiapp.cpp:705-707 (quoted at
    // `cli::main_entry`: blastdbcmd and dustmasker keep the default `DryRun`, which returns 0
    // after the arguments have passed their definitions; LOSAT refuses the settings it refuses
    // in a run first). NCBI BLAST+ 2.17.0 dustmasker: `-in missing.fa -dryrun` and
    // `-infmt blastdb -in missing -dryrun` give exit status 0 and no output.
    if dryrun {
        return 0;
    }
    body(
        &args,
        &settings,
        &mut files::Cin::new(input),
        out,
        &mut diag,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_run_s_arg_exception_is_reported_after_the_usage_with_its_backlog() {
        static TOOL: Tool = Tool {
            program: "p",
            refusal_name: "P",
            version: "1.0.0",
            help: "USAGE\n  p [-h]\n\nOPTIONAL ARGUMENTS\n -h\n",
            descs: &[],
            dependencies: no_dependencies,
            post_prefix: None,
        };
        let mut stderr = Vec::new();
        let mut diag = Diag {
            target: DiagTarget::Stderr,
            stderr: &mut stderr,
            tool: &TOOL,
        };
        let error = ArgException::new("eNoFile", "Argument \"in\". X".to_string());
        assert_eq!(diag.run_arg_exception(error), 1);
        assert_eq!(
            String::from_utf8(stderr).unwrap(),
            format!(
                "USAGE\n  p [-h]\n\nUse '-help' to print detailed descriptions of command line \
                 arguments\n{}\n\nError: Argument \"in\". X\nError: Application's execution \
                 failed\nError:  (CArgException::eNoFile) Argument \"in\". X\n",
                "=".repeat(72)
            )
        );
        // `eNoValue` prints no usage.
        let mut stderr = Vec::new();
        let mut diag = Diag {
            target: DiagTarget::Stderr,
            stderr: &mut stderr,
            tool: &TOOL,
        };
        diag.run_arg_exception(ArgException::new("eNoValue", "m".to_string()));
        assert_eq!(
            String::from_utf8(stderr).unwrap(),
            "Error: m\nError: Application's execution failed\nError:  (CArgException::eNoValue) \
             m\n"
        );
    }
}
