//! `-window_masker_db` and `-window_masker_taxid`: the query masks of WindowMasker
//! (`utils::winmask`), which replace the DUST and lower-case masks of the queries, and NCBI's
//! texts for its errors.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:216-224
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:280-313
//! The masks run after DUST in `x_CalculateMasks` and `SetMaskedRegions` replaces every mask
//! of the query, its DUST and lower-case masks too (NCBI 2.17.0 gives the same bytes with
//! and without `-dust` and `-lcase_masking`). NCBI reads the counts file for each query batch
//! and query chunk (`s_BuildSeqMasker`); LOSAT reads it once, when the filtering options are
//! read (`BlastnArgs::resolve_window_masker`), and keeps the result for the run.

use crate::wording::{ENGINE_ERROR, QUERY_OPTIONS_ERROR};
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
    /// What NCBI writes at the run's first load of the masker and at each later one
    /// (`load_posts`): the posts of the registry of `-window_masker_taxid`, which it loads to
    /// find the file (`ncbi_environment::WindowMaskerRegistry`; none for `-window_masker_db`,
    /// which loads no registry), then the warnings of the counts reader
    /// (`UnitCounts::warnings`, at every load).
    first_load: Posts,
    later_loads: Posts,
}

/// The posts of a load of the masker: those of the registry (`registry`), then the counts
/// reader's warnings when the masker was made (NCBI 2.17.0, line NP M3d: an ascii file's
/// second, higher `>t_low` writes its warning at each load, after the post of a malformed
/// `.ncbirc`).
fn load_posts(registry: &Posts, masker: &Result<SeqMasker, WinMaskError>) -> Posts {
    let mut posts = registry.clone();
    if let Ok(masker) = masker {
        for warning in masker.warnings() {
            posts.lines(format!("Warning: [blastn] {warning}\n").as_bytes());
        }
    }
    posts
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
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251
/// A file that cannot be opened is `CATCH_ALL`'s exception (exit 255, its report and an empty
/// line); a format other than the two binary ones is the input error of the filtering options
/// (exit 1). LOSAT's check reads the file on the stream functions NCBI's does; its other
/// errors (the stream could not seek) give the unknown format, as NCBI's stream failure does.
pub(crate) fn open(name: &str, source: &dyn CountsSource) -> Result<WindowMasker, NativeError> {
    match source.discover(name) {
        Ok(found) if matches!(found.stat_type, StatType::OBinary | StatType::Binary) => {}
        Err(error) if error.discover_exception().is_some() => {
            return Err(NativeError::exception(
                &error.discover_exception().expect("checked above"),
            ));
        }
        Ok(_) | Err(_) => {
            // See docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md exception 12.
            return Err(NativeError {
                exit: 1,
                message: format!(
                    "{QUERY_OPTIONS_ERROR} Only optimized binary windowmasker stat files are supported\n{}\n",
                    crate::wording::USER_MANUAL_FOOTER
                ).into(),
            });
        }
    }
    let masker = build(name, source);
    let posts = load_posts(&Posts::new(), &masker);
    Ok(WindowMasker {
        name: name.to_string(),
        taxid: None,
        masker,
        first_load: posts.clone(),
        later_loads: posts,
    })
}

/// `s_BuildSeqMasker` of `lstat`: the masker made from the counts file `name`, or its error.
fn build(name: &str, source: &dyn CountsSource) -> Result<SeqMasker, WinMaskError> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:99-109
    let params = MaskerParams::default();
    source
        .counts(name, &params)
        .and_then(|counts| SeqMasker::new(counts, &params))
}

/// The masker of `-window_masker_taxid <taxid>` (not 0) with the counts file
/// `<directory>/<taxid>/wmasker.obinary` (`registry`: the directory of the registry entry
/// `[WINDOW_MASKER] WINDOW_MASKER_PATH`, `None` or empty for the current directory, and what
/// each load of the registry writes, `load_posts`; LOSAT Web's registry of its registered files),
/// which NCBI looks for when it makes the masks of each query batch and chunk (in `source`): a
/// path that is no file is the empty name, whose error
/// (`could not open `) each batch reports (`batch_error`). NCBI checks no format of the
/// file (unlike `-window_masker_db`'s, `open`): the ascii and optimized ascii counts are read
/// too.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:418-422
/// The directory and the file do not change during a run, so LOSAT looks for the file once.
pub(crate) fn open_taxid(
    taxid: i32,
    registry: &WindowMaskerRegistry,
    source: &dyn CountsSource,
) -> anyhow::Result<WindowMasker> {
    let directory = registry.path.as_deref();
    let file = crate::utils::winmask::taxid_stat_file(directory, taxid, source).map_err(|error| {
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
    // The counts that only LOSAT refuses (NCBI's result is undefined there) are refused when
    // the options are read, before any output.
    if let Err(error @ (WinMaskError::UnitSize { .. } | WinMaskError::UnsetHashParameter { .. })) =
        &masker
    {
        return Err(rejection(&name, error));
    }
    let first_load = load_posts(&registry.load_posts, &masker);
    let later_loads = if registry.posts_once {
        load_posts(&Posts::new(), &masker)
    } else {
        first_load.clone()
    };
    Ok(WindowMasker {
        name,
        taxid: Some(taxid),
        masker,
        first_load,
        later_loads,
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
    /// (`s_FindPathToWM`), which may post about `.ncbirc` (`WindowMaskerRegistry::load_posts`),
    /// and then reads the counts file, whose reader may post warnings (`s_BuildSeqMasker`).
    /// `first`: the run's first load (the set-up of its first query batch).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_objmgr_tools.cpp:216-224
    /// NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:411-422
    pub(crate) fn load_posts(&self, first: bool) -> &Posts {
        if first {
            &self.first_load
        } else {
            &self.later_loads
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
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
        WinMaskError::BadIndex { .. } => error.to_string().into_bytes(),
        _ => error.create_report()?,
    };
    // NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_istat_obinary.cpp:171-174
    let mut warnings = Posts::new();
    let mut message = format!("{ENGINE_ERROR} ").into_bytes();
    for query in 0..queries {
        warnings.append(if query == 0 { first_load } else { later_loads });
        for warning in error.warnings() {
            warnings.lines(format!("Warning: [blastn] {warning}\n").as_bytes());
        }
        message.extend_from_slice(b"Warning: ");
        message.extend_from_slice(&text);
        message.push(b' ');
    }
    message.push(b'\n');
    Some((warnings, NativeError { exit: 3, message }))
}

/// LOSAT's rejection of the errors that NCBI does not report as `batch_error` does: an ascii
/// file's unit size above 16 and an optimized ascii file's unset hash parameter (both where
/// NCBI's result is undefined; `open_taxid` refuses them before any output), a cache bit array of 0 words that NCBI reads past (`WinMaskError::EmptyCacheBitArray`), a
/// query of 2^32 letters or more, an error of a counts source of LOSAT Web, and an error of the
/// masks of a chunk of a split batch whose whole queries were masked (NCBI masks each chunk
/// when it searches it, after the batch's set-up).
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
#[path = "window_masker_tests.rs"]
mod tests;
