//! `-window_masker_db` and `-window_masker_taxid`: the query masks of WindowMasker
//! (`utils::winmask`), which replace the DUST and lower-case masks of the queries, and NCBI's
//! texts for its errors.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:216-224
//! ```c++
//!         if (NULL != m_Options->GetWindowMaskerDatabase() ||
//!             0    != m_Options->GetWindowMaskerTaxId()) {
//!
//!             if (m_QueryVector.NotEmpty()) {
//!                 Blast_FindWindowMaskerLoc(*m_QueryVector, m_Options);
//!             } else {
//!                 Blast_FindWindowMaskerLoc(*m_TSeqLocVector, m_Options);
//!             }
//!         }
//! ```
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:280-313
//! ```c++
//! Blast_FindWindowMaskerLoc(CBlastQueryVector & queries, const string & lstat)
//! {
//!     AutoPtr<CSeqMasker> masker(s_BuildSeqMasker(lstat));
//!
//!     for(size_t j = 0; j < queries.Size(); j++) {
//!         ...
//!         CSeqVector psv(*seqloc,
//!                        *queries.GetScope(j),
//!                        CBioseq_Handle::eCoding_Iupac,
//!                        eNa_strand_plus);
//!         ...
//!         AutoPtr<CSeqMasker::TMaskList> pos_masks((*masker)(psv));
//!
//!         TMaskedQueryRegions mqr;
//!
//!         s_BuildMaskedRanges(*pos_masks,
//!                             *seqloc,
//!                             *query_seq_id,
//!                             & mqr,
//!                             0);
//!
//!         query.SetMaskedRegions(mqr);
//!     }
//! }
//! ```
//! The masks run after DUST in `x_CalculateMasks` and `SetMaskedRegions` replaces every mask
//! of the query, its DUST and lower-case masks too (NCBI 2.17.0 gives the same bytes with
//! and without `-dust` and `-lcase_masking`). NCBI reads the counts file for each query batch
//! and query chunk (`s_BuildSeqMasker`); LOSAT reads it once, when the filtering options are
//! read (`BlastnArgs::resolve_window_masker`), and keeps the result for the run.

use std::fmt;
use std::sync::Arc;

use crate::blastinput::ncbi_environment::WindowMaskerRegistry;
use crate::cli::NativeError;
use crate::report::blast_xml::filter::WindowMasker as FilterOption;
use crate::report::diag_log::Posts;
use crate::utils::dust::MaskedInterval;
use crate::utils::winmask::{CountsSource, MaskerParams, SeqMasker, StatType, WinMaskError};

/// The masker of `-window_masker_db` or `-window_masker_taxid`, or the error of its counts
/// file, which NCBI reports when it sets up the first query batch (`batch_error`).
pub struct WindowMasker {
    /// The counts file: the value of `-window_masker_db` as given (BLAST XML's filter string
    /// `W -d <name>;`, the archive's `WindowMaskerDatabase`), or the file of
    /// `-window_masker_taxid` (empty when there is none, as in NCBI).
    pub name: String,
    /// The value of `-window_masker_taxid` (not 0; the filter string's `W -t <taxid>;`, the
    /// archive's `WindowMaskerTaxId`), or `None` for `-window_masker_db`.
    pub taxid: Option<i32>,
    masker: Result<SeqMasker, WinMaskError>,
    /// What NCBI writes each time it loads the registry of `-window_masker_taxid` to find the
    /// file (`ncbi_environment::WindowMaskerRegistry`; `load_posts`); empty for
    /// `-window_masker_db`, which loads no registry.
    registry_posts: Posts,
    /// Whether only the run's first load writes `registry_posts`.
    registry_posts_once: bool,
}

impl fmt::Debug for WindowMasker {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("WindowMasker")
            .field("name", &self.name)
            .field("taxid", &self.taxid)
            .field("masker", &self.masker.as_ref().map(|_| "SeqMasker"))
            .finish()
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:446-458
/// ```c++
///     if (args.Exist(kArgWindowMaskerDatabase) &&
///         args[kArgWindowMaskerDatabase]) {
///         const string& stat_file = args[kArgWindowMaskerDatabase].AsString();
///         const CSeqMaskerIstatFactory::EStatType type =
///             CSeqMaskerIstatFactory::DiscoverStatType(stat_file);
///         if (type != CSeqMaskerIstatFactory::eOBinary &&
///             type != CSeqMaskerIstatFactory::eBinary) {
///             string msg("Only optimized binary windowmasker stat files are supported");
///             NCBI_THROW(CInputException, eInvalidInput, msg);
///         }
///
///         opt.SetWindowMaskerDatabase(stat_file.c_str());
///         filter_dbs++;
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
/// ```c++
///     catch (const CException& e) {                                           \
///         LOG_POST(Error << "Error: " << e.what());                           \
///         exit_code = BLAST_UNKNOWN_ERROR;                                    \
///     }                                                                       \
/// ```
/// A file that cannot be opened is `CATCH_ALL`'s exception (exit 255, its report and an empty
/// line); a format other than the two binary ones is the input error of the filtering options
/// (exit 1). LOSAT's check reads the file on the stream functions NCBI's does; its other
/// errors (the stream could not seek) give the unknown format, as NCBI's stream failure does.
pub(crate) fn open(name: &str, source: &dyn CountsSource) -> Result<WindowMasker, NativeError> {
    match source.discover(name) {
        Ok(found) if matches!(found.stat_type, StatType::OBinary | StatType::Binary) => {}
        Err(error) if error.discover_report().is_some() => {
            return Err(NativeError {
                exit: 255,
                message: format!(
                    "Error: {}\n",
                    error.discover_report().expect("checked above")
                ),
            });
        }
        Ok(_) | Err(_) => {
            return Err(NativeError {
                exit: 1,
                message: "BLAST query/options error: Only optimized binary windowmasker stat files are supported\nPlease refer to the BLAST+ user manual.\n".to_string(),
            });
        }
    }
    Ok(WindowMasker {
        name: name.to_string(),
        taxid: None,
        masker: build(name, source),
        registry_posts: Posts::new(),
        registry_posts_once: false,
    })
}

/// `s_BuildSeqMasker( lstat )`: the masker made from the counts file `name`, or its error.
fn build(name: &str, source: &dyn CountsSource) -> Result<SeqMasker, WinMaskError> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:99-109
    // ```c++
    // CSeqMasker* s_BuildSeqMasker(const string & lstat)
    // {
    //     Uint1 arg_window_size            = 0; // [allow setting of this field?]
    //     Uint4 arg_window_step            = 1;
    //     Uint1 arg_unit_step              = 1;
    //     Uint4 arg_textend                = 0; // [allow setting of this field?]
    //     Uint4 arg_cutoff_score           = 0; // [allow setting of this field?]
    // ```
    let params = MaskerParams::default();
    source
        .counts(name, &params)
        .and_then(|counts| SeqMasker::new(counts, &params))
}

/// The masker of `-window_masker_taxid <taxid>` (not 0) with the counts file
/// `<directory>/<taxid>/wmasker.obinary` (`registry`: the directory of the registry entry
/// `[WINDOW_MASKER] WINDOW_MASKER_PATH`, `None` or empty for the current directory, and what
/// each load of the registry writes, `load_posts`), which NCBI looks for when it makes the masks
/// of each query batch and chunk: a path that is no file is the empty name, whose error
/// (`could not open `) each batch reports (`batch_error`). NCBI checks no format of the
/// file (unlike `-window_masker_db`'s, `open`) and reads the ascii counts too, which LOSAT does not
/// read: it refuses them before the search.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:418-422
/// ```c++
/// Blast_FindWindowMaskerLocTaxId(CBlastQueryVector & queries, int taxid)
/// {
///     string db = WindowMaskerTaxidToDb(taxid);
///     Blast_FindWindowMaskerLoc(queries, db);
/// }
/// ```
/// The directory and the file do not change during a run, so LOSAT looks for the file once.
pub(crate) fn open_taxid(
    taxid: i32,
    registry: &WindowMaskerRegistry,
    source: &dyn CountsSource,
) -> anyhow::Result<WindowMasker> {
    let directory = registry.path.as_deref();
    let file = crate::utils::winmask::taxid_stat_file(directory, taxid).map_err(|error| {
        anyhow::anyhow!(
            "-window_masker_taxid {taxid}: the current directory, where NCBI BLAST+ looks for the counts file, cannot be read ({error}); this is not supported by LOSAT's BLASTN"
        )
    })?;
    let name = match file {
        Some(path) => path.into_os_string().into_string().map_err(|path| {
            anyhow::anyhow!(
                "-window_masker_taxid {taxid}: the path of the counts file {} is not valid UTF-8; this is not supported by LOSAT's BLASTN",
                std::path::Path::new(&path).display()
            )
        })?,
        None => String::new(),
    };
    let masker = build(&name, source);
    if let Err(WinMaskError::UnsupportedFormat(format)) = &masker {
        let format = match format {
            StatType::Ascii => "ascii",
            StatType::OAscii => "optimized ascii",
            StatType::Binary => "binary",
            StatType::OBinary => "optimized binary",
            StatType::Unknown => "unknown",
        };
        anyhow::bail!(
            "-window_masker_taxid {taxid}: the windowmasker unit counts {name} are in the {format} format, which NCBI BLAST+ reads and LOSAT's BLASTN does not; this is not supported by LOSAT's BLASTN"
        );
    }
    Ok(WindowMasker {
        name,
        taxid: Some(taxid),
        masker,
        registry_posts: registry.load_posts.clone(),
        registry_posts_once: registry.posts_once,
    })
}

impl WindowMasker {
    /// What BLAST XML's filter string writes for the masker (`W -t <taxid>;` or `W -d <name>;`).
    pub(crate) fn filter_option(&self) -> FilterOption<'_> {
        match self.taxid {
            Some(taxid) => FilterOption::Taxid(taxid),
            None => FilterOption::Database(&self.name),
        }
    }

    /// The error of the counts file when NCBI could not make the masker from it.
    pub(crate) fn build_error(&self) -> Option<&WinMaskError> {
        self.masker.as_ref().err()
    }

    /// What NCBI writes when it makes the masks of a query factory (`x_CalculateMasks`, once
    /// for a query batch and then once for each chunk of a split batch): with
    /// `-window_masker_taxid` it loads the registry that names the file's directory then
    /// (`s_FindPathToWM`), which may post about `.ncbirc` (`WindowMaskerRegistry::load_posts`).
    /// `first`: the run's first load (the set-up of its first query batch).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:216-224
    /// ```c++
    ///         if (NULL != m_Options->GetWindowMaskerDatabase() ||
    ///             0    != m_Options->GetWindowMaskerTaxId()) {
    ///
    ///             if (m_QueryVector.NotEmpty()) {
    ///                 Blast_FindWindowMaskerLoc(*m_QueryVector, m_Options);
    /// ```
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:411-422
    /// ```c++
    /// string WindowMaskerTaxidToDb(int taxid)
    /// {
    ///     string path = s_FindPathToWM();
    ///     return WindowMaskerTaxidToDb(path, taxid);
    /// }
    ///
    /// void
    /// Blast_FindWindowMaskerLocTaxId(CBlastQueryVector & queries, int taxid)
    /// {
    ///     string db = WindowMaskerTaxidToDb(taxid);
    /// ```
    pub(crate) fn load_posts(&self, first: bool) -> &Posts {
        static NONE: Posts = Posts::new();
        if self.registry_posts_once && !first {
            &NONE
        } else {
            &self.registry_posts
        }
    }

    /// `batch_error` of a batch whose masks failed with `error`, with the registry's posts of
    /// its queries' loads (`first`: the batch is the run's first).
    pub(crate) fn batch_error(
        &self,
        error: &WinMaskError,
        queries: usize,
        first: bool,
    ) -> Option<(Posts, NativeError)> {
        batch_error(
            error,
            queries,
            self.load_posts(first),
            self.load_posts(false),
        )
    }

    /// The masks of a query whose searched letters are `letters` (`psv`): WindowMasker's
    /// intervals of the letters in upper case, or with `-strand minus` of their reverse
    /// complement, moved by nothing (the letters start at the query's location).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:162-189
    /// ```c++
    ///     TSeqPos query_start = seqloc.GetStart(eExtreme_Positional);
    ///     ...
    ///     ITERATE(CSeqMasker::TMaskList, pr, masks) {
    ///         CRef<CSeq_interval> ival(new CSeq_interval);
    ///
    ///         TSeqPos
    ///             start  = pr->first,
    ///             end    = pr->second;
    ///
    ///         ival->SetFrom (query_start + start);
    ///         ival->SetTo   (query_start + end);
    ///         ival->SetId   (query_id);
    ///         ival->SetStrand(eNa_strand_both);
    ///
    ///         if (mqr) {
    ///             CRef<CSeqLocInfo> info_plus
    ///                 (new CSeqLocInfo(&* ival, CSeqLocInfo::eFramePlus1));
    ///             mqr->push_back(info_plus);
    ///
    ///             CRef<CSeqLocInfo> info_minus
    ///                 (new CSeqLocInfo(&* ival, CSeqLocInfo::eFrameMinus1));
    ///             mqr->push_back(info_minus);
    ///         }
    /// ```
    /// `psv` is the vector of the query's Seq-loc on `eNa_strand_plus` relative to the
    /// location: with `-strand minus` the location is on the minus strand and the vector is
    /// the reverse complement of the searched letters. The intervals of that vector are put
    /// on the plus strand of the query unchanged (NCBI 2.17.0: the masks of `-strand minus` are
    /// windowmasker's intervals of the reverse complement), and every interval masks both
    /// strands. The masker sees only A, C, G and T (`CSeqMaskerWindow`), so the other
    /// letters are all one to it.
    pub(crate) fn masks(
        &self,
        letters: &[u8],
        minus: bool,
    ) -> Result<Vec<MaskedInterval>, WinMaskError> {
        let masker = self.masker.as_ref().map_err(Clone::clone)?;
        let data: Vec<u8> = if minus {
            super::lookup::reverse_complement(letters)
        } else {
            letters.to_ascii_uppercase()
        };
        Ok(masker
            .mask(&data)?
            .into_iter()
            .map(|(start, end)| MaskedInterval::new(start as usize, end as usize + 1))
            .collect())
    }

    /// The masks of the queries of a batch or chunk, in order; the first error stops it, as it
    /// stops NCBI's loop over the queries.
    pub(crate) fn batch_masks<'q>(
        &self,
        letters: impl IntoIterator<Item = &'q [u8]>,
        minus: bool,
    ) -> Result<Vec<Vec<MaskedInterval>>, WinMaskError> {
        letters
            .into_iter()
            .map(|letters| self.masks(letters, minus))
            .collect()
    }
}

/// NCBI's stop at the set-up of a query batch whose masks could not be made (also a batch
/// that NCBI splits: it sets up the whole batch before it splits it), with `queries` queries
/// in it: the messages that NCBI posts while
/// it reads the counts file for each query, and the error. Each query's read of the file comes
/// after a load of the registry of `-window_masker_taxid` (`WindowMasker::load_posts`):
/// `first_load` is the batch's first query's post, `later_loads` the others'.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:149-151
/// ```c++
///     catch(CException & e) {
///         NCBI_THROW(CBlastException, eSetup, e.what());
///     }
/// ```
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
/// ```c++
///         } catch (const CException& e) {
///             ...
///             CRef<CSearchMessage> m
///                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
///             messages[index].push_back(m);
///             s_InvalidateQueryContexts(qinfo, index);
///         }
///     ...
///     if (BlastSetup_Validate(qinfo, NULL) != 0 && messages.HasMessages()) {
///         NCBI_THROW(CBlastException, eSetup, messages.ToString());
///     }
/// ```
/// `x_CalculateMasks` is called again for each query (`GetMaskedRegions`) while it has not
/// succeeded, so each query reads the file again and gets the same exception: the query's
/// message is the severity, the exception's text and a space (`TSearchMessages::ToString`),
/// and `CATCH_ALL` writes `BLAST engine error: ` (exit 3). The message of the reader's
/// exception is its report (`what()`); that of a bad index of the hash table, raised while
/// the masker runs, is its own text. `None` for the errors that NCBI does not report this way
/// (LOSAT's own rejections, `rejection`). NCBI BLAST+ 2.17.0 on 2026-11-12 (line NP M3c,
/// `taxid_rc_pre_missing`): a malformed `.ncbirc` and no counts file give the registry's post
/// once for each of the batch's 4 queries, then the error.
pub(crate) fn batch_error(
    error: &WinMaskError,
    queries: usize,
    first_load: &Posts,
    later_loads: &Posts,
) -> Option<(Posts, NativeError)> {
    let text = match error {
        WinMaskError::BadIndex { .. } => error.to_string(),
        _ => error.create_report()?,
    };
    // NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:171-174
    // ```c
    //             else if( !input_stream.read( (char *)cba, cba_size*sizeof( Uint4 ) ) )
    //             {
    //                 ERR_POST( Warning << "file read failed: "
    //                                   << "bit array optimizations are not used." );
    // ```
    let mut warnings = Posts::new();
    let mut message = String::from("BLAST engine error: ");
    for query in 0..queries {
        warnings.append(if query == 0 { first_load } else { later_loads });
        for warning in error.warnings() {
            warnings.lines(format!("Warning: [blastn] {warning}\n").as_bytes());
        }
        message.push_str("Warning: ");
        message.push_str(&text);
        message.push(' ');
    }
    message.push('\n');
    Some((warnings, NativeError { exit: 3, message }))
}

/// LOSAT's rejection of the errors that NCBI does not report as `batch_error` does: a
/// format LOSAT does not read (a counts source that names one after the check), the masker's
/// window check (not reached with blastn's settings), a cache bit array of 0 words that NCBI
/// reads past (`WinMaskError::EmptyCacheBitArray`), a query of 2^32 letters or more, and an
/// error of the masks of a chunk of a split batch whose whole queries were masked (NCBI masks
/// each chunk when it searches it, after the batch's set-up).
pub(crate) fn rejection(name: &str, error: &WinMaskError) -> anyhow::Error {
    anyhow::anyhow!(
        "the windowmasker unit counts {name} give an error that LOSAT's BLASTN does not report as NCBI BLAST+ does ({error}); this is not supported by LOSAT's BLASTN"
    )
}

/// The counts source of a run: the one the caller gave (LOSAT Web), or the file system.
pub(crate) fn source_of(given: Option<&Arc<dyn CountsSource>>) -> Arc<dyn CountsSource> {
    given
        .cloned()
        .unwrap_or_else(|| Arc::new(crate::utils::winmask::FileCounts))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ncbi_exception::NCBI_BUILD_ROOT;
    use crate::utils::winmask::{FileCounts, FormatSite};

    fn counts_path(name: &str) -> String {
        format!("{}/tests/winmask/counts/{name}", env!("CARGO_MANIFEST_DIR"))
    }

    // NCBI BLAST+ 2.17.0 `windowmasker -ustat mg1655_u9.obinary -outfmt interval` of the
    // first record of tests/fixtures/blastn_regression/inputs/wm_query.fa (132 - 146, ...),
    // and of its reverse complement (390 - 402, ...), which are the masks of blastn
    // `-window_masker_db mg1655_u9.obinary` and of `-strand minus` (its archive's LCaseMask).
    #[test]
    fn the_masks_are_windowmaskers_of_the_strand_read() {
        let fasta = std::fs::read_to_string(format!(
            "{}/tests/fixtures/blastn_regression/inputs/wm_query.fa",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let first: String = fasta.split('>').nth(1).unwrap().lines().skip(1).collect();
        let masker = open(&counts_path("mg1655_u9.obinary"), &FileCounts).unwrap();
        let plus = masker.masks(first.as_bytes(), false).unwrap();
        assert_eq!(
            plus[..3],
            [
                MaskedInterval::new(132, 147),
                MaskedInterval::new(414, 429),
                MaskedInterval::new(489, 502)
            ]
        );
        // The lower-case letters of the record are letters like the others.
        assert_eq!(
            masker
                .masks(first.to_ascii_uppercase().as_bytes(), false)
                .unwrap(),
            plus
        );
        let minus = masker.masks(first.as_bytes(), true).unwrap();
        assert_eq!(
            minus[..2],
            [MaskedInterval::new(390, 403), MaskedInterval::new(559, 572)]
        );
    }

    #[test]
    fn the_checks_of_the_file_are_ncbis() {
        let missing = open("nosuch.obinary", &FileCounts).unwrap_err();
        assert_eq!(missing.exit, 255);
        assert_eq!(
            missing.message,
            format!("Error: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 171: Error: (Exception::open failed) ncbi::CSeqMaskerIstatFactory::DiscoverStatType() - could not open nosuch.obinary\n\n")
        );
        // A directory opens and reads nothing on Unix: the unknown format. Windows cannot open a
        // directory as a file, so the open fails first there: the "could not open" exception.
        let directory = open(env!("CARGO_MANIFEST_DIR"), &FileCounts).unwrap_err();
        if cfg!(unix) {
            assert_eq!(directory.exit, 1);
            assert_eq!(
                directory.message,
                "BLAST query/options error: Only optimized binary windowmasker stat files are supported\nPlease refer to the BLAST+ user manual.\n"
            );
        } else {
            assert_eq!(directory.exit, 255);
            assert!(
                directory.message.ends_with(&format!(
                    "could not open {}\n\n",
                    env!("CARGO_MANIFEST_DIR")
                )),
                "{}",
                directory.message
            );
        }
        let binary = open(&counts_path("dm2L_u9.binary"), &FileCounts).unwrap();
        assert!(binary.build_error().is_none());
        assert_eq!(binary.name, counts_path("dm2L_u9.binary"));
    }

    // NCBI BLAST+ 2.17.0 blastn `-window_masker_taxid` (line NP M/STATE.md, "NCBI 2.17.0 の実機"):
    // the file under the directory is used whatever its format is called; a path that is no file
    // is the empty name, which fails each batch with `could not open ` and the factory's error;
    // LOSAT refuses the ascii counts that NCBI reads, before the search.
    #[test]
    fn the_file_of_a_taxid_is_the_masker_or_the_empty_name() {
        let dir = std::env::temp_dir().join(format!(
            "losat-blastn-window-masker-taxid-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        for (taxid, file) in [("562", "mg1655_u9.obinary"), ("7227", "dm2L_u9.binary")] {
            std::fs::create_dir_all(dir.join(taxid)).unwrap();
            std::fs::copy(counts_path(file), dir.join(taxid).join("wmasker.obinary")).unwrap();
        }
        std::fs::create_dir_all(dir.join("10")).unwrap();
        std::fs::write(
            dir.join("10/wmasker.obinary"),
            b"##windowmasker-statistics-format-version:ascii 1.1.0\n9\nAAAAAAAAA 5\n",
        )
        .unwrap();
        let root = WindowMaskerRegistry {
            path: Some(dir.to_str().unwrap().to_string()),
            ..WindowMaskerRegistry::default()
        };
        let root = &root;
        let found = open_taxid(562, root, &FileCounts).unwrap();
        assert_eq!(found.taxid, Some(562));
        assert!(found.name.ends_with("wmasker.obinary"), "{}", found.name);
        assert!(found.build_error().is_none());
        assert_eq!(found.filter_option(), FilterOption::Taxid(562));
        assert!(found.load_posts(true).is_empty() && found.load_posts(false).is_empty());
        // What each load of the registry writes, or only the first.
        let mut posted = Posts::new();
        posted.post(crate::report::diag_log::CRITICAL, "posted\n");
        for (once, later) in [(false, posted.clone()), (true, Posts::new())] {
            let registry = WindowMaskerRegistry {
                load_posts: posted.clone(),
                posts_once: once,
                ..root.clone()
            };
            let masker = open_taxid(562, &registry, &FileCounts).unwrap();
            assert_eq!(masker.load_posts(true), &posted);
            assert_eq!(masker.load_posts(false), &later);
        }
        // A binary counts file under the optimized binary's name.
        assert!(open_taxid(7227, root, &FileCounts)
            .unwrap()
            .build_error()
            .is_none());
        let missing = open_taxid(9606, root, &FileCounts).unwrap();
        assert_eq!(missing.name, "");
        let none = &Posts::new();
        let (warnings, error) = batch_error(missing.build_error().unwrap(), 2, none, none).unwrap();
        assert!(warnings.is_empty());
        // Each query's read of the file comes after a load of the registry.
        let (first, later) = (Posts::of_lines(b"first\n".to_vec()), posted);
        let (posts, _) = batch_error(missing.build_error().unwrap(), 3, &first, &later).unwrap();
        assert_eq!(posts.text(), b"first\nposted\nposted\n");
        let one = format!("Warning: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 171: Error: (Exception::open failed) ncbi::CSeqMaskerIstatFactory::DiscoverStatType() - could not open \n    T0 \"{NCBI_BUILD_ROOT}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 358: Error: (Exception::creation failure) ncbi::CSeqMaskerIstatFactory::create() - could not create a unit counts container\n ");
        assert_eq!(error.exit, 3);
        assert_eq!(error.message, format!("BLAST engine error: {one}{one}\n"));
        let ascii = open_taxid(10, root, &FileCounts).unwrap_err().to_string();
        assert!(ascii.contains("ascii format"), "{ascii}");
        assert!(ascii.contains("not supported by LOSAT"), "{ascii}");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    // NCBI BLAST+ 2.17.0 (the official Linux x64 release) blastn with a counts file whose unit
    // size is 17, and with one whose cache bit array is cut short (line NP M/STATE.md, "M2 の
    // 実機確認"; the bytes are those of tests/blastn_window_masker.rs's fixtures).
    #[test]
    fn a_counts_file_that_fails_stops_the_batch_with_one_message_per_query() {
        let none = &Posts::new();
        let factory = format!("    T0 \"{NCBI_BUILD_ROOT}/src/algo/winmask/seq_masker_istat_factory.cpp\", line 358: Error: (Exception::creation failure) ncbi::CSeqMaskerIstatFactory::create() - could not create a unit counts container\n");
        let unit = format!("Warning: NCBI C++ Exception:\n    T0 \"{NCBI_BUILD_ROOT}/src/algo/winmask/seq_masker_istat_obinary.cpp\", line 114: Error: (Exception::bad parameter) ncbi::CSeqMaskerIstatOBinary::CSeqMaskerIstatOBinary() - unit size must be in [1,16]\n{factory} ");
        let (warnings, error) =
            batch_error(&WinMaskError::format(FormatSite::UnitSize), 2, none, none).unwrap();
        assert!(warnings.is_empty());
        assert_eq!(error.exit, 3);
        assert_eq!(error.message, format!("BLAST engine error: {unit}{unit}\n"));
        let short = WinMaskError::Format {
            site: FormatSite::HashTable,
            warnings: vec!["file read failed: bit array optimizations are not used.".into()],
        };
        let (warnings, _) = batch_error(&short, 2, none, none).unwrap();
        assert_eq!(
            warnings.text(),
            b"Warning: [blastn] file read failed: bit array optimizations are not used.\nWarning: [blastn] file read failed: bit array optimizations are not used.\n"
        );
        let (_, bad_index) = batch_error(
            &WinMaskError::BadIndex {
                key: 8764,
                value: 14779,
            },
            1,
            none,
            none,
        )
        .unwrap();
        assert_eq!(
            bad_index.message,
            "BLAST engine error: Warning: bad index at key 8764 : 14779 \n"
        );
        assert!(
            batch_error(&WinMaskError::EmptyCacheBitArray { unit: 2 }, 1, none, none).is_none()
        );
    }
}
