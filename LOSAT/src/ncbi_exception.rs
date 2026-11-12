//! NCBI's exception reports: the text of `CException::what()` that NCBI's applications print
//! after `Error: ` (their last `catch (const CException&)` clause), and the warnings that quote
//! an exception, with one `T0` line per exception of the chain; and the lines that
//! `NCBI_REPORT_EXCEPTION_X` posts with NCBI's default diagnostic flags
//! ([`NcbiException::report_default_ex`]: rpsbproc's application handler, `rpsbproc::report_exception`).
//! This module is the one owner of that text; each place NCBI throws from keeps its
//! [`ThrowSite`] next to the port of the code that throws (`ncbi_str`, `writedb::data_size`,
//! `writedb::build_db`, `utils::winmask`, `seqid::parse`, `blastinput::db_filter`, BLASTN's
//! `-use_index` warning in `blastn::blast_engine::run`, the serial objects' exceptions of
//! rpsbproc in `rpsbproc::query`).
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:370-384
//! ```c
//! string CException::ReportAll(TDiagPostFlags flags) const
//! {
//!     // invert the order
//!     stack<const CException*> pile;
//!     const CException* pex;
//!     for (pex = this; pex; pex = pex->GetPredecessor()) {
//!         pile.push(pex);
//!     }
//!     CNcbiOstrstream os;
//!     os << "NCBI C++ Exception:" << '\n';
//!     for (; !pile.empty(); pile.pop()) {
//!         //indentation
//!         os << "    ";
//!         os << pile.top()->ReportThis(flags) << '\n';
//!     }
//! ```
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:397-406
//! ```c
//! string CException::ReportThis(TDiagPostFlags flags) const
//! {
//!     CNcbiOstrstream os, osex;
//!     ReportStd(os, flags);
//!     ReportExtra(osex);
//!     if ( !IsOssEmpty(osex) ) {
//!         os << " (" << (string)CNcbiOstrstreamToString(osex) << ')';
//!     }
//!     return CNcbiOstrstreamToString(os);
//! }
//! ```
//! NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:5624-5646
//! ```c
//! #ifdef NCBI_THREADS
//!     if (IsSetDiagPostFlag(eDPF_TID, m_Flags)) {
//!         os << 'T' << CThread::GetSelf() << ' ';
//!     }
//! #endif
//!     // "<file>"
//!     bool print_file = (m_File  &&  *m_File  &&
//!                        IsSetDiagPostFlag(eDPF_File, m_Flags));
//!     ...
//!         os << '"' << x_file << '"';
//!     ...
//!         os << (print_file ? ", line " : "line ") << m_Line;
//! ```
//! `ReportStd` writes the message with `SDiagMessage::x_OldWrite` (ncbidiag.cpp:5612-): the
//! thread (`T0`: the report is made in the main thread's `catch`), the source file as
//! `__FILE__` names it in NCBI's build tree, the line, the severity (`Error`), the class and
//! error code (`GetType()::GetErrCodeString()`), the function (`module::class::function() - `,
//! when NCBI knows one) and the message (`GetMsg()`). `CStringException` adds its position
//! (`ReportExtra` of its base `CParseTemplException<CCoreException>`: `m_Pos = <pos>`,
//! ncbistr.hpp:4435-4439, 4499).

/// The build tree that NCBI's reports name: that of the official Linux x64 release of BLAST+
/// 2.17.0 (`ncbi-blast-2.17.0+-x64-linux.tar.gz`). Other builds of 2.17.0 name their own tree
/// and print the same text around it (line NP, M/DECISIONS.md 14).
pub const NCBI_BUILD_ROOT: &str = "/home/coremake/release_build/build/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_1751374548/c++/compilers/unix/../..";

/// Where NCBI throws an exception (`NCBI_THROW`, `NCBI_THROW2`, `NCBI_RETHROW`): what its
/// report prints of the place.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThrowSite {
    /// `__FILE__` below the build tree: `src/corelib/ncbistr.cpp`.
    pub file: &'static str,
    /// `__LINE__` of the throw.
    pub line: u32,
    /// `GetType()::GetErrCodeString()`: `CStringException::eConvert`.
    pub kind: &'static str,
    /// The module, class and function (`NCBI_CURRENT_FUNCTION`, with `NCBI_MODULE` before it
    /// where the library is built with one): `ncbi::NStr::StringToInt8`. `None` where NCBI's
    /// report names none.
    pub function: Option<&'static str>,
}

impl ThrowSite {
    // NCBI reference (598d8ae6): c++/src/serial/datatool/choicestr.cpp:887-893
    // ```c++
    //             "void "<<methodPrefix<<"ThrowInvalidSelection(" STATE_ENUM " index) const\n"
    //             "{\n"
    //             "    throw NCBI_NS_NCBI::CInvalidChoiceSelection(DIAG_COMPILE_INFO";
    //         if ( 1 ) { // add extra argument for better error message
    //             methods << ", this";
    //         }
    //         methods << ", m_choice, index, sm_SelectionNames, sizeof(sm_SelectionNames)/sizeof(sm_SelectionNames[0]));\n"
    // ```
    /// The place of a throw in the code that NCBI's `datatool` generates for a serial type
    /// (`C<type>_Base::ThrowInvalidSelection` of a CHOICE): the generated files are not in NCBI's
    /// source tree, so LOSAT has neither their names and lines nor the function, and leaves
    /// them out (an empty file, line 0, no function). Only the reports that print no place are
    /// made of such an exception ([`NcbiException::report_default_ex`]); its `what()` would not
    /// name the place that NCBI's names.
    pub const fn generated(kind: &'static str) -> ThrowSite {
        ThrowSite {
            file: "",
            line: 0,
            kind,
            function: None,
        }
    }
}

/// An exception of NCBI's C++ Toolkit (`CException`): where it was thrown, its message, the
/// extra text of its class, and the exception it was thrown for (`NCBI_RETHROW`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NcbiException {
    site: ThrowSite,
    message: String,
    extra: Option<String>,
    predecessor: Option<Box<NcbiException>>,
}

impl NcbiException {
    /// `NCBI_THROW(class, code, message)` at `site`.
    pub fn new(site: ThrowSite, message: impl Into<String>) -> Self {
        NcbiException {
            site,
            message: message.into(),
            extra: None,
            predecessor: None,
        }
    }

    // `CStringException` derives from `CParseTemplException<CCoreException>` (ncbistr.hpp:4499),
    // whose `ReportExtra` prints the position.
    // NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:4435-4439
    // ```c++
    //     /// Report error position.
    //     virtual void ReportExtra(ostream& out) const override
    //     {
    //         out << "m_Pos = " << (unsigned long)m_Pos;
    //     }
    // ```
    /// `NCBI_THROW2(CStringException, code, message, pos)` at `site`: the report adds the
    /// position.
    pub fn with_pos(site: ThrowSite, message: impl Into<String>, pos: usize) -> Self {
        NcbiException {
            extra: Some(format!("m_Pos = {pos}")),
            ..Self::new(site, message)
        }
    }

    /// `NCBI_RETHROW(self, class, code, message)` at `site`: a new exception whose predecessor
    /// is `self`.
    pub fn rethrow(self, site: ThrowSite, message: impl Into<String>) -> Self {
        NcbiException {
            predecessor: Some(Box::new(self)),
            ..Self::new(site, message)
        }
    }

    /// Where the exception was thrown (for the tests; the report reads the field in `report_this`).
    #[cfg(test)]
    pub fn site(&self) -> ThrowSite {
        self.site
    }

    /// `GetMsg()`: the message, without the place or the extra text.
    pub fn get_msg(&self) -> &str {
        &self.message
    }

    /// `ReportThis()`: the line of this exception, without its predecessors.
    fn report_this(&self, out: &mut String) {
        let ThrowSite {
            file,
            line,
            kind,
            function,
        } = self.site;
        out.push_str(&format!(
            "T0 \"{NCBI_BUILD_ROOT}/{file}\", line {line}: Error: ({kind}) "
        ));
        if let Some(function) = function {
            out.push_str(function);
            out.push_str("() - ");
        }
        out.push_str(&self.text());
    }

    /// The exceptions of the chain, newest first (`this`, then each `GetPredecessor()`).
    fn chain(&self) -> Vec<&NcbiException> {
        let mut pile = Vec::new();
        let mut pex = Some(self);
        while let Some(e) = pex {
            pile.push(e);
            pex = e.predecessor.as_deref();
        }
        pile
    }

    /// `CException::what()` (`ReportAll`): the report of the chain, oldest exception first,
    /// ending with a line feed.
    pub fn what(&self) -> String {
        let mut os = String::from("NCBI C++ Exception:\n");
        for e in self.chain().iter().rev() {
            os.push_str("    ");
            e.report_this(&mut os);
            os.push('\n');
        }
        os
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7844-7856
    // ```c++
    // static string
    // s_GetExceptionText(const CException* pex)
    // {
    //     string text(pex->GetMsg());
    //     stringstream os;
    //     pex->ReportExtra(os);
    //     string s = os.str();
    //     if ( !s.empty() ) {
    //         text += " (";
    //         text += s;
    //         text += ')';
    //     }
    //     return text;
    // }
    // ```
    /// `s_GetExceptionText`: the message and the extra text of the class (the text of each
    /// exception of a chain that NCBI's diagnostics post, `ncbi_app::args`).
    pub fn text(&self) -> String {
        match &self.extra {
            Some(extra) => format!("{} ({extra})", self.message),
            None => self.message.clone(),
        }
    }

    // NCBI reference (598d8ae6): c++/include/corelib/ncbiexpt.hpp:761-770
    // ```c++
    // #define NCBI_REPORT_EXCEPTION_X(err_subcode, title, ex)                 \
    //     NCBI_REPORT_EXCEPTION_XX(NCBI_USE_ERRCODE_X, err_subcode, title, ex)
    // ...
    // #define NCBI_REPORT_EXCEPTION_XX(err_name, err_subcode, title, ex)   \
    //     NCBI_CHECK_ERR_SUBCODE_X_NAME(err_name, err_subcode);            \
    //     NCBI_NS_NCBI::CExceptionReporter::ReportDefaultEx                \
    //     (NCBI_ERRCODE_X_NAME(err_name), err_subcode,                     \
    //      DIAG_COMPILE_INFO, title, ex, NCBI_NS_NCBI::eDPF_Default)
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:676-700
    // ```c++
    // void CExceptionReporter::ReportDefaultEx(int err_code, int err_subcode,
    //     const CDiagCompileInfo& info, const string& title,const exception& ex,
    //     TDiagPostFlags flags)
    // {
    //     if ( !sm_DefEnabled )
    //         return;
    //     ...
    //     if ( sm_DefHandler ) {
    //         ...
    //     } else {
    //         CNcbiDiag d(info, cex->GetSeverity(), flags);
    //         d.SetOmitStackTrace(true);
    //         d << ErrCode(err_code, err_subcode) << title << " " <<  *cex;
    //     }
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7873-7939
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
    //     ...
    //     for (; !pile.empty(); pile.pop()) {
    //         pex = pile.top();
    //         string text(s_GetExceptionText(pex));
    //         ...
    //         string err_type(pex->GetType());
    //         err_type += "::";
    //         err_type += pex->GetErrCodeString();
    //         ...
    //         SDiagMessage diagmsg
    //             (pex_sev,
    //             text.c_str(),
    //             text.size(),
    //             pex->GetFile().c_str(),
    //             pex->GetLine(),
    //             GetPostFlags(),
    //             NULL,
    //             pex->GetErrCode(),
    //             0,
    //             err_type.c_str(),
    //             ...
    //         m_Buffer.PrintMessage(diagmsg, *this);
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:3984-3987
    // ```c++
    // static const TDiagPostFlags s_DefaultPostFlags =
    //     eDPF_Prefix | eDPF_Severity | eDPF_ErrorID |
    //     eDPF_ErrCodeMessage | eDPF_ErrCodeExplanation |
    //     eDPF_ErrCodeUseSeverity;
    // ```
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:5672-5695
    // ```c++
    //     // <severity>:
    //     if (IsSetDiagPostFlag(eDPF_Severity, m_Flags)  &&
    //         (m_Severity != eDiag_Info || !IsSetDiagPostFlag(eDPF_OmitInfoSev))) {
    //         string sev = CNcbiDiag::SeverityName(m_Severity);
    //         ...
    //             os << sev;
    //         ...
    //         os << ": ";
    //     }
    //
    //     // (<err_code>.<err_subcode>) or (err_text)
    //     if ((m_ErrCode  ||  m_ErrSubCode || m_ErrText)  &&
    //         IsSetDiagPostFlag(eDPF_ErrorID, m_Flags)) {
    //         os << '(';
    //         if (m_ErrText) {
    //             os << m_ErrText;
    //         } else {
    //             os << m_ErrCode << '.' << m_ErrSubCode;
    //         }
    //         os << ") ";
    //     }
    // ```
    // `x_Put` posts a message for each exception of the chain, oldest first, as it runs; the
    // title's message is posted when `d` goes out of scope, after them. With the default flags
    // (an application that sets none, as rpsbproc) `SDiagMessage::x_OldWrite` writes neither the
    // file and line (`eDPF_File`, `eDPF_Line`), the thread (`eDPF_TID`) nor the function
    // (`eDPF_Location`), and no error code description (none without the registry's
    // `[DEBUG] MessageFile`, ncbiapp.cpp:1654-1667); each message ends with a line feed. The
    // severity is the exceptions' `Error` (`sm_DefEnabled` is true and `sm_DefHandler` null
    // unless an application changes them, ncbiexpt.cpp:603-605). The main text is the newest
    // exception's: `rethrow` makes an exception whose predecessor loses it (the constructor's
    // `prev_exception->m_MainText = false`, ncbiexpt.cpp:188-189).
    /// `NCBI_REPORT_EXCEPTION_X(err_subcode, title, e)` in a source file whose
    /// `NCBI_USE_ERRCODE_X` is `err_code`, as an application with NCBI's default diagnostic
    /// flags posts it: `Error: (<class>::<code>) <message>` for each exception of the chain,
    /// oldest first, then `Error: (<err_code>.<err_subcode>) <title> (<class>::<code>)
    /// <message>` of the newest, each line ending with a line feed.
    pub fn report_default_ex(&self, err_code: u32, err_subcode: u32, title: &str) -> String {
        let mut out = String::new();
        for e in self.chain().iter().rev() {
            out.push_str(&post_error(e.site.kind, &e.text()));
        }
        out.push_str(&err_post_x(
            err_code,
            err_subcode,
            &format!("{title} ({}) {}", self.site.kind, self.text()),
        ));
        out
    }
}

/// A message of severity `Error` with the error id `id` (`m_ErrText`, or `<err_code>.<err_subcode>`)
/// as an application with NCBI's default diagnostic flags writes it (`SDiagMessage::x_OldWrite`,
/// quoted at [`NcbiException::report_default_ex`]): `Error: (<id>) <message>` and a line feed.
fn post_error(id: &str, message: &str) -> String {
    format!("Error: ({id}) {message}\n")
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbidiag.hpp:550-551, 564-567
// ```c++
// #define ERR_POST_X(err_subcode, message)                  \
//     ERR_POST_XX(NCBI_USE_ERRCODE_X, err_subcode, message)
// ...
// #define ERR_POST_XX(error_name, err_subcode, message)                   \
//     ERR_POST_EX(NCBI_ERRCODE_X_NAME(error_name),                        \
//                 NCBI_ERR_SUBCODE_X_NAME(error_name, err_subcode),       \
//                 message)
// ```
// NCBI reference (598d8ae6): c++/include/corelib/ncbidiag.hpp:252-256
// ```c++
// #define ERR_POST_EX(err_code, err_subcode, message)          \
//     ( NCBI_NS_NCBI::CNcbiDiag(DIAG_COMPILE_INFO).GetRef()    \
//       << NCBI_NS_NCBI::ErrCode( (err_code), (err_subcode) )  \
//       << message                                             \
//       << NCBI_NS_NCBI::Endm )
// ```
// NCBI reference (598d8ae6): c++/include/corelib/ncbidiag.hpp:934-938
// ```c++
//     NCBI_XNCBI_EXPORT  CNcbiDiag
//     (const CDiagCompileInfo& info,                      ///< File, line, module
//      EDiagSev                sev        = eDiag_Error,  ///< Severity level
//      TDiagPostFlags          post_flags = eDPF_Default  ///< What to post
//      );
// ```
/// `ERR_POST_X(err_subcode, message)` in a source file whose `NCBI_USE_ERRCODE_X` is
/// `err_code`, as an application with NCBI's default diagnostic flags writes it (a `CNcbiDiag`
/// is of severity `Error`): `Error: (<err_code>.<err_subcode>) <message>` and a line feed. The
/// application's report of an exception's title line, and of a `std::exception`
/// (ncbiapp.cpp:675-679: rpsbproc's data loader).
pub fn err_post_x(err_code: u32, err_subcode: u32, message: &str) -> String {
    post_error(&format!("{err_code}.{err_subcode}"), message)
}

#[cfg(test)]
mod tests {
    use super::*;

    // makeblastdb 2.17.0 (the official Linux x64 release) with `MAX_LMDB_TXN_ENTRY=-1x` and two
    // sequences (probe of 2026-11-08, case txn_endptr_2seq): a report with the position.
    #[test]
    fn a_string_exception_reports_its_position() {
        let site = ThrowSite {
            file: "src/corelib/ncbistr.cpp",
            line: 870,
            kind: "CStringException::eConvert",
            function: Some("ncbi::NStr::StringToInt8"),
        };
        let e = NcbiException::with_pos(site, "Cannot convert string '-1x' to Int8", 2);
        assert_eq!(e.get_msg(), "Cannot convert string '-1x' to Int8");
        assert_eq!(e.site(), site);
        assert_eq!(
            e.what(),
            "NCBI C++ Exception:\n    T0 \"/home/coremake/release_build/build/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_1751374548/c++/compilers/unix/../../src/corelib/ncbistr.cpp\", line 870: Error: (CStringException::eConvert) ncbi::NStr::StringToInt8() - Cannot convert string '-1x' to Int8 (m_Pos = 2)\n"
        );
    }

    // blastn 2.17.0 (the official Linux x64 release) `-window_masker_db` with a counts file whose
    // header is cut short (line NP, M/STATE.md "M2 の実機確認"): the reader's exception, without a
    // function, and the factory's rethrow, oldest first.
    #[test]
    fn a_rethrow_reports_the_chain_oldest_first() {
        let reader = NcbiException::new(
            ThrowSite {
                file: "src/algo/winmask/seq_masker_istat_obinary.cpp",
                line: 63,
                kind: "Exception::format error",
                function: None,
            },
            "file too short",
        );
        let factory = reader.rethrow(
            ThrowSite {
                file: "src/algo/winmask/seq_masker_istat_factory.cpp",
                line: 358,
                kind: "Exception::creation failure",
                function: Some("ncbi::CSeqMaskerIstatFactory::create"),
            },
            "could not create a unit counts container",
        );
        assert_eq!(
            factory.get_msg(),
            "could not create a unit counts container"
        );
        assert_eq!(
            factory.what(),
            "NCBI C++ Exception:\n    T0 \"/home/coremake/release_build/build/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_1751374548/c++/compilers/unix/../../src/algo/winmask/seq_masker_istat_obinary.cpp\", line 63: Error: (Exception::format error) file too short\n    T0 \"/home/coremake/release_build/build/PrepareRelease_Linux64-Alma_JSID_01_170042_130.14.18.128_9008__PrepareRelease_Linux64-Alma_1751374548/c++/compilers/unix/../../src/algo/winmask/seq_masker_istat_factory.cpp\", line 358: Error: (Exception::creation failure) ncbi::CSeqMaskerIstatFactory::create() - could not create a unit counts container\n"
        );
    }

    // blastn 2.17.0 (the official Linux x64 release) with a `blastn.ini` of the line `garbage`
    // (env parity g_reg probe `x_ini`, 2026-11-11): `NCBI_REPORT_EXCEPTION_X(15, "Application's
    // initialization failed", e)` of the application layer (error code 106), and of a chain.
    #[test]
    fn a_report_with_a_title_has_a_line_per_exception_then_the_title() {
        let site = ThrowSite {
            file: "src/corelib/ncbireg.cpp",
            line: 722,
            kind: "CRegistryException::eEntry",
            function: Some("ncbi::IRWRegistry::x_Read"),
        };
        let e = NcbiException::with_pos(site, "Invalid registry entry format: 'garbage'", 1);
        assert_eq!(
            e.report_default_ex(106, 15, "Application's initialization failed"),
            "Error: (CRegistryException::eEntry) Invalid registry entry format: 'garbage' (m_Pos = 1)\nError: (106.15) Application's initialization failed (CRegistryException::eEntry) Invalid registry entry format: 'garbage' (m_Pos = 1)\n"
        );
        let outer = NcbiException::new(site, "first").rethrow(
            ThrowSite {
                kind: "CAppException::eSetupDiag",
                ..site
            },
            "second",
        );
        assert_eq!(
            outer.report_default_ex(106, 15, "T"),
            "Error: (CRegistryException::eEntry) first\nError: (CAppException::eSetupDiag) second\nError: (106.15) T (CAppException::eSetupDiag) second\n"
        );
    }

    // NCBI's rpsbproc 0.5.1 (tests/fixtures/rpsbproc/oracle/d_wsline.err and lcl.err): the
    // application's report of an exception with the default diagnostic flags, without the place
    // or the function, the position of a `CStringException` included; and that of an exception
    // of generated code.
    #[test]
    fn the_default_report_names_the_class_and_the_message() {
        let site = ThrowSite {
            file: "src/corelib/ncbistr.cpp",
            line: 895,
            kind: "CStringException::eConvert",
            function: Some("ncbi::NStr::StringToUInt8"),
        };
        let e = NcbiException::with_pos(site, "Cannot convert string '' to Uint8", 0);
        assert_eq!(
            e.report_default_ex(106, 15, "Application's initialization failed"),
            "Error: (CStringException::eConvert) Cannot convert string '' to Uint8 (m_Pos = 0)\nError: (106.15) Application's initialization failed (CStringException::eConvert) Cannot convert string '' to Uint8 (m_Pos = 0)\n"
        );
        let e = NcbiException::new(
            ThrowSite::generated("CInvalidChoiceSelection::eFail"),
            "CSeq_id::GetGeneral(): Invalid choice selection: NCBI-Seqloc::Seq-id.local",
        );
        assert_eq!(
            e.report_default_ex(106, 16, "Application's execution failed"),
            "Error: (CInvalidChoiceSelection::eFail) CSeq_id::GetGeneral(): Invalid choice selection: NCBI-Seqloc::Seq-id.local\nError: (106.16) Application's execution failed (CInvalidChoiceSelection::eFail) CSeq_id::GetGeneral(): Invalid choice selection: NCBI-Seqloc::Seq-id.local\n"
        );
    }

    /// The source files of the crate below `dir`, with their paths.
    fn rust_files(dir: &std::path::Path, out: &mut Vec<std::path::PathBuf>) {
        for entry in std::fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                rust_files(&path, out);
            } else if path.extension().is_some_and(|ext| ext == "rs") {
                out.push(path);
            }
        }
    }

    /// The lines of `text` that are code outside the tests: comment lines and the top-level
    /// `#[cfg(test)] mod name { ... }` blocks are left out (rustfmt puts the closing brace of a
    /// top-level module at the start of a line).
    fn code_lines(text: &str) -> Vec<(usize, &str)> {
        let mut lines = Vec::new();
        let mut in_tests = false;
        let mut after_cfg_test = false;
        for (number, line) in text.lines().enumerate() {
            if in_tests {
                in_tests = line != "}";
                continue;
            }
            if line == "#[cfg(test)]" {
                after_cfg_test = true;
                continue;
            }
            if after_cfg_test {
                if line.starts_with("#[") {
                    continue;
                }
                after_cfg_test = false;
                if line.starts_with("mod ") && line.ends_with('{') {
                    in_tests = true;
                    continue;
                }
            }
            if !line.trim_start().starts_with("//") {
                lines.push((number + 1, line));
            }
        }
        lines
    }

    // The text of an exception report has one owner: a second format string for it (as four
    // modules had until 2026-11-08: winmask, the Seq-id parser, the taxid map and the seqid list
    // filter) fails here. Test files (`tests.rs`, `*_tests.rs`, files under a `tests/` directory)
    // and test modules hold NCBI's bytes as expected values and are not read. The needles are
    // the escaped forms in ordinary string literals: a raw string (`r#"..."#`) or a text put
    // together with `concat!` or from pieces is not caught.
    #[test]
    fn exception_reports_are_formatted_only_here() {
        let src = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
        let mut files = Vec::new();
        rust_files(&src, &mut files);
        assert!(
            files.len() > 100,
            "{} files below {}",
            files.len(),
            src.display()
        );
        let mut offenders = Vec::new();
        for path in files {
            let relative = path.strip_prefix(&src).unwrap();
            let stem = path.file_stem().unwrap().to_string_lossy();
            let test_file = stem == "tests"
                || stem.ends_with("_tests")
                || relative
                    .components()
                    .any(|part| part.as_os_str() == "tests");
            if relative == std::path::Path::new("ncbi_exception.rs") || test_file {
                continue;
            }
            let text = std::fs::read_to_string(&path).unwrap();
            for (number, line) in code_lines(&text) {
                if ["NCBI C++ Exception:\\n", "T0 \\\"", ": Error: ("]
                    .iter()
                    .any(|needle| line.contains(needle))
                {
                    offenders.push(format!("{}:{number}: {}", relative.display(), line.trim()));
                }
            }
        }
        assert!(
            offenders.is_empty(),
            "format an NCBI exception report with ncbi_exception::NcbiException:\n{}",
            offenders.join("\n")
        );
    }

    #[test]
    fn code_lines_leave_out_comments_and_test_modules() {
        let text =
            "fn a() {}\n// \"NCBI C++ Exception:\\n\"\n#[cfg(test)]\n#[path = \"x.rs\"]\nmod x;\n\
                    fn b() {}\n#[cfg(test)]\nmod tests {\n    fn c() {}\n}\nfn d() {}\n";
        assert_eq!(
            code_lines(text),
            [
                (1, "fn a() {}"),
                (5, "mod x;"),
                (6, "fn b() {}"),
                (11, "fn d() {}")
            ]
        );
    }
}
