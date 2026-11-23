//! NCBI's exception reports: the text of `CException::what()` that NCBI's applications print
//! after `Error: ` (their last `catch` clause, which takes a `CException` by reference), and the
//! warnings that quote an exception, with one `T0` line per exception of the chain; and the lines
//! that `NCBI_REPORT_EXCEPTION_X` posts with NCBI's default diagnostic flags
//! ([`NcbiException::report_default_ex`]: rpsbproc's application handler,
//! `rpsbproc::report_exception`). This module is the one owner of that text; each place NCBI throws
//! from keeps its [`ThrowSite`] next to the port of the code that throws (`ncbi_str`,
//! `writedb::data_size`, `writedb::build_db`, `utils::winmask`, `seqid::parse`,
//! `blastinput::db_filter`, BLASTN's `-use_index` warning in `blastn::blast_engine::run`, the
//! serial objects' exceptions of rpsbproc in `rpsbproc::query`).
//!
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:370-384
//! NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:397-406
//! NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:5624-5646
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

// NCBI reference (598d8ae6): c++/include/corelib/ncbiobj.hpp:965-975
// NCBI reference (598d8ae6): c++/src/corelib/ncbiobj.cpp:1002-1014
/// The message of NCBI's exception at a read through a `CRef` of no object
/// (`CObject::ThrowNullPointerException`, `CCoreException::eNullPtr`).
pub const NULL_POINTER: &str = "Attempt to access NULL pointer.";

/// The build tree that NCBI rpsbproc 0.5.1's reports name: that of its Linux x64 binary of the
/// 2025-05-13 package (the oracle of `tests/fixtures/rpsbproc`), built apart from BLAST+ on its own
/// toolkit. Its `ncbireg.cpp` throws at the lines of BLAST+ 2.17.0's (env parity probe I of
/// 2026-11-14 on ws: the seven throws of `IRWRegistry::x_Read`); bioconda's rpsbproc 0.5.1 names
/// another tree (`/opt/conda/conda-bld/rpsbproc_<n>/work/ncbi_cxx/c++`).
pub const RPSBPROC_BUILD_ROOT: &str =
    "/export/home/tomcat/TeamCity/Agent2/work/b895474bf74772ad/struct_tree";

/// The build tree that the reports of NCBI's `program` name: rpsbproc's own
/// ([`RPSBPROC_BUILD_ROOT`]), else BLAST+ 2.17.0's ([`NCBI_BUILD_ROOT`]).
pub fn build_root(program: &str) -> &'static str {
    match program {
        "rpsbproc" => RPSBPROC_BUILD_ROOT,
        _ => NCBI_BUILD_ROOT,
    }
}

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
///
/// The message is bytes, as NCBI's `std::string`: a name or a word of the input that is not UTF-8
/// (a member name of a text ASN.1 file, a file name) keeps its bytes in the report.
///
/// The chain of predecessors is as long as the frames an error leaves (`rethrow_same`: one for
/// each value open in a PSSM nested to any depth), so it is dropped, cloned, compared and
/// shown link by link, never by a recursion as deep as the chain (below).
pub struct NcbiException {
    site: ThrowSite,
    message: Vec<u8>,
    extra: Option<String>,
    predecessor: Option<Box<NcbiException>>,
}

impl NcbiException {
    /// `NCBI_THROW` with a class, a code and a message, at `site`.
    pub fn new(site: ThrowSite, message: impl Into<Vec<u8>>) -> Self {
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
    /// `NCBI_THROW2` with `CStringException`, a code, a message and `pos`, at `site`: the report
    /// adds the position.
    pub fn with_pos(site: ThrowSite, message: impl Into<Vec<u8>>, pos: usize) -> Self {
        let mut exception = Self::new(site, message);
        exception.extra = Some(format!("m_Pos = {pos}"));
        exception
    }

    /// `NCBI_RETHROW(self, class, code, message)` at `site`: a new exception whose predecessor
    /// is `self`.
    pub fn rethrow(self, site: ThrowSite, message: impl Into<Vec<u8>>) -> Self {
        let mut exception = Self::new(site, message);
        exception.predecessor = Some(Box::new(self));
        exception
    }

    // NCBI reference (598d8ae6): c++/include/corelib/ncbiexpt.hpp:749-752
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:274-285
    /// `NCBI_RETHROW_SAME(self, message)` in `function` at `line` of `file`: the same exception
    /// (its class, code and extra text, such as a `CStringException`'s position) with the new
    /// place and message, the old one its predecessor: NCBI's serial object stack adds one for
    /// each frame it leaves (`END_OBJECT_FRAME`: the writer's `report::serial::xml`, the reader's
    /// `report::serial::asn_types::read`).
    pub fn rethrow_same(
        self,
        file: &'static str,
        line: u32,
        function: &'static str,
        message: impl Into<Vec<u8>>,
    ) -> Self {
        NcbiException {
            site: ThrowSite {
                file,
                line,
                kind: self.site.kind,
                function: Some(function),
            },
            message: message.into(),
            extra: self.extra.clone(),
            predecessor: Some(Box::new(self)),
        }
    }

    /// Where the exception was thrown (for the tests; the report reads the field in `report_this`).
    #[cfg(test)]
    pub fn site(&self) -> ThrowSite {
        self.site
    }

    /// `GetType()::GetErrCodeString()` of the exception: `CSeqDBException::eInvalid`, which
    /// tells the class that the programs' `CATCH_ALL` reports apart (`cli::NativeError::exception`).
    pub fn kind(&self) -> &'static str {
        self.site.kind
    }

    /// Whether LOSAT knows the place of every throw of the chain, which `what()` names (a
    /// [`ThrowSite::generated`] names none).
    pub fn places_known(&self) -> bool {
        self.chain().iter().all(|e| !e.site.file.is_empty())
    }

    // NCBI reference (598d8ae6): c++/src/serial/exception.cpp:47-61
    /// `CSerialException::AddFrameInfo` with `frame_info`, which the object stack's `END_OBJECT_FRAME`
    /// calls with `GetFrameName()` of each frame a `CSerialException` leaves: the frame stack is
    /// the extra text of the exception's report (`ReportExtra`, ` at <stack>`; a
    /// `CSerialException` has no other). `NCBI_RETHROW_SAME` (`rethrow_same`) keeps the stack in
    /// the new exception and freezes a copy in its predecessor, as `AddBacklog` clones it.
    pub(crate) fn add_frame_info(&mut self, frame_info: &str) {
        if frame_info.is_empty() {
            return;
        }
        let stack = match self.extra.take() {
            Some(extra) => {
                let stack = extra.strip_prefix(" at ").unwrap_or(&extra);
                format!(" at {frame_info}.{stack}")
            }
            None => format!(" at {frame_info}"),
        };
        self.extra = Some(stack);
    }

    /// `GetMsg()`: the message, without the place or the extra text.
    pub fn get_msg(&self) -> &[u8] {
        &self.message
    }

    /// `ReportThis()`: the line of this exception, without its predecessors, with the build tree
    /// `root`.
    fn report_this(&self, out: &mut Vec<u8>, root: &str) {
        let ThrowSite {
            file,
            line,
            kind,
            function,
        } = self.site;
        out.extend_from_slice(
            format!("T0 \"{root}/{file}\", line {line}: Error: ({kind}) ").as_bytes(),
        );
        if let Some(function) = function {
            out.extend_from_slice(function.as_bytes());
            out.extend_from_slice(b"() - ");
        }
        out.extend_from_slice(&self.text());
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
    /// ending with a line feed, in BLAST+ 2.17.0's build ([`NCBI_BUILD_ROOT`]).
    pub fn what(&self) -> Vec<u8> {
        self.what_in(NCBI_BUILD_ROOT)
    }

    /// [`NcbiException::what`] of the build whose tree is `root` ([`build_root`]), in LOSAT's
    /// words: the messages of the chain, without NCBI's header and places
    /// (`crate::wording::losat_text`, `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 12).
    pub fn what_in(&self, root: &str) -> Vec<u8> {
        crate::wording::losat_bytes(&self.ncbi_what_in(root)).into_owned()
    }

    /// `<text><what()>`, as NCBI posts an exception after a text (an `ERR_POST` of `e.GetMsg()` and
    /// `e.what()`, or of the text "... Q(0/1" and `e.what()`), in LOSAT's words: the joined text is
    /// mapped as a whole, so the report is set off from the text before it as in the mapping of
    /// NCBI's post (`crate::wording::losat_text`).
    pub(crate) fn what_after(&self, text: &[u8]) -> Vec<u8> {
        crate::wording::losat_bytes(&[text, &self.ncbi_what_in(NCBI_BUILD_ROOT)].concat())
            .into_owned()
    }

    /// NCBI's text of [`NcbiException::what`], which LOSAT's maps (`crate::wording::losat_bytes`).
    pub(crate) fn ncbi_what(&self) -> Vec<u8> {
        self.ncbi_what_in(NCBI_BUILD_ROOT)
    }

    /// NCBI's text of [`NcbiException::what_in`], which LOSAT's maps.
    fn ncbi_what_in(&self, root: &str) -> Vec<u8> {
        let mut os = b"NCBI C++ Exception:\n".to_vec();
        for e in self.chain().iter().rev() {
            os.extend_from_slice(b"    ");
            e.report_this(&mut os, root);
            os.push(b'\n');
        }
        os
    }

    /// The exception as the error of a report's writer, where NCBI's formatter throws it while the
    /// report is written (the merge of a query-anchored section,
    /// `report::query_anchored::write_section`): the writers before the program's report pass
    /// `io::Error`s on, and the report finds the exception again ([`NcbiException::of_io_error`],
    /// `cli::NativeError::of_writer`), as NCBI's leaves the formatter for the application's
    /// `CATCH_ALL`.
    pub(crate) fn into_io_error(self) -> std::io::Error {
        std::io::Error::other(self)
    }

    /// The exception that an error of a report's writer carries ([`NcbiException::into_io_error`]).
    pub(crate) fn of_io_error(error: &std::io::Error) -> Option<&NcbiException> {
        error.get_ref()?.downcast_ref()
    }

    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7844-7856
    /// `s_GetExceptionText`: the message and the extra text of the class (the text of each
    /// exception of a chain that NCBI's diagnostics post, `ncbi_app::args`).
    /// In LOSAT's words (`crate::wording::losat_text`: the position as `(position <n>)`).
    pub(crate) fn text(&self) -> Vec<u8> {
        match &self.extra {
            Some(extra) => crate::wording::losat_bytes(
                &[self.message.as_slice(), b" (", extra.as_bytes(), b")"].concat(),
            )
            .into_owned(),
            None => self.message.clone(),
        }
    }

    // NCBI reference (598d8ae6): c++/include/corelib/ncbiexpt.hpp:761-770
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:676-700
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:7873-7939
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:3984-3987
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:5672-5695
    // `x_Put` posts a message for each exception of the chain, oldest first, as it runs; the
    // title's message is posted when `d` goes out of scope, after them. With the default flags
    // (an application that sets none, as rpsbproc) `SDiagMessage::x_OldWrite` writes neither the
    // file and line (`eDPF_File`, `eDPF_Line`), the thread (`eDPF_TID`) nor the function
    // (`eDPF_Location`), and no error code description (none without the registry's
    // `[DEBUG] MessageFile`, ncbiapp.cpp:1654-1667); each message ends with a line feed. The
    // severity is the exceptions' `Error` (`sm_DefEnabled` is true and `sm_DefHandler` null
    // unless an application changes them, ncbiexpt.cpp:603-605). The main text is the newest
    // exception's: `rethrow` makes an exception whose predecessor loses it (the constructor's
    // setting `prev_exception->m_MainText` to false, ncbiexpt.cpp:188-189).
    /// `NCBI_REPORT_EXCEPTION_X` with `err_subcode`, `title` and the exception `e`, in a source
    /// file whose `NCBI_USE_ERRCODE_X` is `err_code`, as an application with NCBI's default
    /// diagnostic flags posts it: `Error: (<class>::<code>) <message>` for each exception of the
    /// chain, oldest first, then `Error: (<err_code>.<err_subcode>) <title> (<class>::<code>)
    /// <message>` of the newest, each line ending with a line feed; in LOSAT's words
    /// (`crate::wording::losat_text`, `PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 12).
    pub fn report_default_ex(&self, err_code: u32, err_subcode: u32, title: &str) -> Vec<u8> {
        crate::wording::losat_bytes(&self.ncbi_report_default_ex(err_code, err_subcode, title))
            .into_owned()
    }

    /// NCBI's text of [`Self::report_default_ex`].
    fn ncbi_report_default_ex(&self, err_code: u32, err_subcode: u32, title: &str) -> Vec<u8> {
        let mut out = Vec::new();
        for e in self.chain().iter().rev() {
            out.extend(post_error_bytes(e.site.kind, &e.text()));
        }
        out.extend(err_post_x_bytes(
            err_code,
            err_subcode,
            &self.titled_text(title),
        ));
        out
    }

    /// `<title> (<class>::<code>) <message>` of the newest exception, the line that ends NCBI's
    /// report of an exception with a title.
    fn titled_text(&self, title: &str) -> Vec<u8> {
        [
            format!("{title} ({}) ", self.site.kind).as_bytes(),
            &self.text(),
        ]
        .concat()
    }

    // NCBI reference (598d8ae6): c++/include/corelib/ncbiexpt.hpp:580-584, 755-757
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:670-674
    /// `NCBI_CATCH` with `title`, of this exception in an application whose diagnostic prefix is
    /// `prefix` (`ReportDefaultEx` without an error code, cited at `report_default_ex`): `Error:
    /// (<class>::<code>) [<prefix>] <message>` for each exception of the chain, oldest first, then
    /// `Error: [<prefix>] <title> (<class>::<code>) <message>` of the newest (NCBI BLAST+ 2.17.0,
    /// line RG D3-1 `om/` cases t08, t10, t76: the Seq-table checks of the object manager).
    pub(crate) fn report_default(&self, prefix: &str, title: &str) -> Vec<u8> {
        let mut out = Vec::new();
        for e in self.chain().iter().rev() {
            out.extend(error_post_line(Some(e.site.kind), Some(prefix), &e.text()));
        }
        out.extend(error_post_line(
            None,
            Some(prefix),
            &self.titled_text(title),
        ));
        out
    }
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbiexpt.cpp:259-265 (each exception deletes its
// predecessor: the chain goes with the newest)
/// The chain is dropped link by link (the drop of each link inside the one that holds it would
/// recurse as deep as the chain).
impl Drop for NcbiException {
    fn drop(&mut self) {
        let mut next = self.predecessor.take();
        while let Some(mut link) = next {
            next = link.predecessor.take();
        }
    }
}

/// A copy of the chain, made from its oldest link.
impl Clone for NcbiException {
    fn clone(&self) -> Self {
        let mut copy: Option<Box<NcbiException>> = None;
        for link in self.chain().into_iter().rev() {
            copy = Some(Box::new(NcbiException {
                site: link.site,
                message: link.message.clone(),
                extra: link.extra.clone(),
                predecessor: copy,
            }));
        }
        *copy.expect("the chain holds the exception itself")
    }
}

/// Two chains are equal when each link is (its place, message and extra text).
impl PartialEq for NcbiException {
    fn eq(&self, other: &Self) -> bool {
        let (mut a, mut b) = (Some(self), Some(other));
        loop {
            match (a, b) {
                (None, None) => return true,
                (Some(x), Some(y)) => {
                    if x.site != y.site || x.message != y.message || x.extra != y.extra {
                        return false;
                    }
                    a = x.predecessor.as_deref();
                    b = y.predecessor.as_deref();
                }
                _ => return false,
            }
        }
    }
}

impl Eq for NcbiException {}

/// The links of the chain, newest first, each with its place, message and extra text.
impl std::fmt::Debug for NcbiException {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        struct Link<'a>(&'a NcbiException);
        impl std::fmt::Debug for Link<'_> {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                f.debug_struct("NcbiException")
                    .field("site", &self.0.site)
                    .field("message", &String::from_utf8_lossy(&self.0.message))
                    .field("extra", &self.0.extra)
                    .finish()
            }
        }
        f.debug_list()
            .entries(self.chain().into_iter().map(Link))
            .finish()
    }
}

/// `what()` ([`NcbiException::what`]) as text, for the places that show an error as text only (a
/// panic, a test's message): a writer of messages finds the exception in an error chain and writes
/// the bytes of `what()` (`cli::NativeError`, `api::error::chain_bytes`).
impl std::fmt::Display for NcbiException {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(&self.what()))
    }
}

/// The error of a report's writer ([`NcbiException::into_io_error`]).
impl std::error::Error for NcbiException {}

/// A message of severity `Error` with the error id `id` (`m_ErrText`, or `<err_code>.<err_subcode>`)
/// as an application with NCBI's default diagnostic flags writes it (`SDiagMessage::x_OldWrite`,
/// cited at [`NcbiException::report_default_ex`]): `Error: (<id>) <message>` and a line feed.
/// The message is bytes: NCBI's diagnostic stream writes them as they are.
fn post_error_bytes(id: &str, message: &[u8]) -> Vec<u8> {
    error_post_line(Some(id), None, message)
}

// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:5685-5694, 5749-5755
/// A post of severity `Error` as `SDiagMessage::x_OldWrite` writes it with NCBI's default flags
/// (cited above and at [`NcbiException::report_default_ex`]): `Error: `, the error id `(<id>) `
/// when the post has one (`<err_code>.<err_subcode>` or `m_ErrText`), the application's prefix
/// `[<prefix>] ` when it set one (`SetDiagPostPrefix`), the message as it is, and a line feed.
/// The one place that makes this line: the exception reports here, the posts of the bytes that
/// NCBI's object streams replace (`report::serial::encoder::ReplacedChar::post`) and the
/// applications' `ERR_POST` of an `Error` (`ncbi_app::Diag::err_post_bytes`).
pub(crate) fn error_post_line(id: Option<&str>, prefix: Option<&str>, message: &[u8]) -> Vec<u8> {
    post_line(Severity::Error, id, prefix, message)
}

/// The severity of a post (`CNcbiDiag::SeverityName`: `Warning`, `Error`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Severity {
    Warning,
    Error,
}

/// A post of severity `severity` from NCBI's main thread (`error_post_line`'s line, of any
/// severity): after the main thread's diagnostic stream failed (`main_diag_failbit`), the header
/// alone.
pub(crate) fn post_line(
    severity: Severity,
    id: Option<&str>,
    prefix: Option<&str>,
    message: &[u8],
) -> Vec<u8> {
    post_line_in(severity, id, prefix, message, main_diag_failbit())
}

/// `post_line` with the state of the stream given: `failed` drops the message (the post's header,
/// `<severity>: (<id>) [<prefix>] `, and its line feed stay; see `main_diag_failbit`).
pub(crate) fn post_line_in(
    severity: Severity,
    id: Option<&str>,
    prefix: Option<&str>,
    message: &[u8],
    failed: bool,
) -> Vec<u8> {
    let mut line = match severity {
        Severity::Warning => b"Warning: ".to_vec(),
        Severity::Error => b"Error: ".to_vec(),
    };
    if let Some(id) = id {
        line.extend_from_slice(format!("({id}) ").as_bytes());
    }
    if let Some(prefix) = prefix {
        line.extend_from_slice(format!("[{prefix}] ").as_bytes());
    }
    if !failed {
        line.extend_from_slice(message);
    }
    line.push(b'\n');
    // LOSAT's words: an exception's class and code (`id`) and NCBI's exception reports are left
    // out (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 12).
    crate::wording::losat_bytes(&line).into_owned()
}

// NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:676-697
// NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:4389-4435 (`CDiagBuffer::Flush`)
// `s` is an output string stream: its buffer gives nothing to read, so inserting `s.rdbuf()` adds no
// character and sets `failbit` on the thread's diagnostic stream (`CDiagBuffer::m_Stream`), which
// `Flush` empties and gives its flags back but never `clear()`s. Every later message of that
// thread is lost: its post prints the header alone, and `CATCH_ALL`'s `LOG_POST` (no header) its
// line feed (NCBI BLAST+ 2.17.0, line RG D3-1 `om/` cases k01-k06 and `om/mute/` e01-e16: the query
// warning line, the taxonomy database warning, the engine error, with 1 to 4 threads).
// The flag is NCBI's for one run (one application, one main thread), and each LOSAT run owns
// its own (`MainDiag`): the command line's process is one run (`COMMAND_LINE`); a library call
// (`api::start::CallStart`) holds a stream of its own, which the call's code reads on the threads
// that run it (`MainDiag::enter` where the call runs, `utils::threading::with_search_pool` for the
// workers of its searches). Two calls in one process, one after the other, at once on several
// threads, or as two `run_batches` iterators that one thread advances in turn, never see each
// other's stream (line RG D3-1 review 4 B1).
static COMMAND_LINE: MainDiag = MainDiag::new();

thread_local! {
    /// The stream of the library call whose code this thread runs (`MainDiag::enter`).
    static CURRENT: std::cell::RefCell<Option<std::sync::Arc<MainDiag>>> =
        const { std::cell::RefCell::new(None) };
}

/// The diagnostic stream of NCBI's main thread of one run: whether it has failed.
#[derive(Debug, Default)]
pub(crate) struct MainDiag {
    failed: std::sync::atomic::AtomicBool,
}

impl MainDiag {
    pub(crate) const fn new() -> MainDiag {
        MainDiag {
            failed: std::sync::atomic::AtomicBool::new(false),
        }
    }

    /// The stream fails.
    pub(crate) fn fail(&self) {
        self.failed
            .store(true, std::sync::atomic::Ordering::Relaxed);
    }

    /// Whether the stream has failed.
    pub(crate) fn failed(&self) -> bool {
        self.failed.load(std::sync::atomic::Ordering::Relaxed)
    }

    /// The stream of the run whose code this thread runs: the library call's that this thread
    /// entered (`enter`), else the command line's.
    fn with_current<R>(read: impl FnOnce(&MainDiag) -> R) -> R {
        CURRENT.with(|current| match current.borrow().as_deref() {
            Some(diag) => read(diag),
            None => read(&COMMAND_LINE),
        })
    }

    /// The library call's stream that this thread runs for (`None`: the command line's), to hand to
    /// the threads that work for the same run (`utils::threading::with_search_pool`, which a build
    /// without threads does not have).
    #[cfg(feature = "parallel")]
    pub(crate) fn current() -> Option<std::sync::Arc<MainDiag>> {
        CURRENT.with(|current| current.borrow().clone())
    }

    /// A thread made for one search of a run (a worker of its pool, `with_search_pool`) takes
    /// that run's stream (`current` on the thread that made it) for the rest of its life.
    #[cfg(feature = "parallel")]
    pub(crate) fn adopt(diag: Option<std::sync::Arc<MainDiag>>) {
        CURRENT.with(|current| *current.borrow_mut() = diag);
    }

    /// This thread runs `diag`'s run (a library call) until the guard is dropped, which gives
    /// the thread back the stream it had.
    pub(crate) fn enter(diag: Option<std::sync::Arc<MainDiag>>) -> MainDiagGuard {
        MainDiagGuard {
            before: CURRENT.with(|current| current.replace(diag)),
        }
    }
}

/// The guard of `MainDiag::enter`.
pub(crate) struct MainDiagGuard {
    before: Option<std::sync::Arc<MainDiag>>,
}

impl Drop for MainDiagGuard {
    fn drop(&mut self) {
        let before = self.before.take();
        CURRENT.with(|current| *current.borrow_mut() = before);
    }
}

/// Whether the diagnostic stream of NCBI's main thread has failed in this run (`ERR_POST_X` of
/// the object manager, subcode 6, with "Failed to parse location of " and `s.rdbuf()`,
/// `algorithm::pssm_search::annot_index`): the run's later posts print their header only
/// (`post_line`) and `CATCH_ALL` its line feeds (`cli::catch_all_text`). False unless a psiblast or
/// tblastn `-in_pssm` run posted it: the command line's, or a library call's (`Inputs::pssm`; each
/// call is a run with a stream of its own, `MainDiag`); LOSAT Web refuses `-in_pssm`.
pub(crate) fn main_diag_failbit() -> bool {
    MainDiag::with_current(MainDiag::failed)
}

/// Marks NCBI's main-thread diagnostic stream failed for the rest of the run
/// (`main_diag_failbit`); called where the run writes the post that fails it.
pub(crate) fn set_main_diag_failbit() {
    MainDiag::with_current(MainDiag::fail);
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbidiag.hpp:550-551, 564-567
// NCBI reference (598d8ae6): c++/include/corelib/ncbidiag.hpp:252-256
// NCBI reference (598d8ae6): c++/include/corelib/ncbidiag.hpp:934-938
/// `ERR_POST_X` with `err_subcode` and `message` in a source file whose `NCBI_USE_ERRCODE_X` is
/// `err_code`, as an application with NCBI's default diagnostic flags writes it (a `CNcbiDiag`
/// is of severity `Error`): `Error: (<err_code>.<err_subcode>) <message>` and a line feed. The
/// application's report of an exception's title line, and of a `std::exception`
/// (ncbiapp.cpp:675-679: rpsbproc's data loader).
pub fn err_post_x(err_code: u32, err_subcode: u32, message: &str) -> String {
    String::from_utf8(err_post_x_bytes(err_code, err_subcode, message.as_bytes()))
        .expect("the text of UTF-8 parts is UTF-8")
}

/// `err_post_x` of a message that is bytes (a `std::exception`'s `what()` that quotes data
/// bytes, rpsbproc's motif parse error).
pub fn err_post_x_bytes(err_code: u32, err_subcode: u32, message: &[u8]) -> Vec<u8> {
    post_error_bytes(&format!("{err_code}.{err_subcode}"), message)
}

#[cfg(test)]
#[path = "ncbi_exception_tests.rs"]
mod tests;
