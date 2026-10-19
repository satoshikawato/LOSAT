//! WindowMasker: the masking part of NCBI's `algo/winmask` (`CSeqMasker`), shared by
//! blastn `-window_masker_db` and the windowmasker tool.
//!
//! Ported: the masker with the settings blastn builds (`s_BuildSeqMasker`,
//! `c++/src/algo/blast/api/winmask_filter.cpp:99-154`) and the integer thresholds of the
//! windowmasker application (`-window`, `-t_extend`, `-t_thres`, `-t_high`, `-t_low`,
//! `-set_t_high`, `-set_t_low`); the optimized binary and binary counts readers; the format
//! check. Not ported: the merge pass, the "min" trigger and discontiguous units (neither caller
//! can turn them on, `win_mask_config.cpp:336-343`), the `-t_*_pct` settings, and the ascii
//! readers.

mod istat;
mod score;
mod taxid;
mod window;

#[cfg(test)]
mod tests;

use std::fmt;
use std::path::Path;

pub use istat::{discover_stat_type, StatFile, StatType, UnitCounts};
pub use taxid::taxid_stat_file;

use score::ScoreMean;
use window::Window;

/// The arguments of `CSeqMasker` that its callers set. `Default` is what blastn passes
/// (all 0: the window is unit + 4 letters and the thresholds come from the counts file).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/winmask_filter.cpp:101-109
/// ```c
///     Uint1 arg_window_size            = 0; // [allow setting of this field?]
///     Uint4 arg_window_step            = 1;
///     Uint1 arg_unit_step              = 1;
///     Uint4 arg_textend                = 0; // [allow setting of this field?]
///     Uint4 arg_cutoff_score           = 0; // [allow setting of this field?]
///     Uint4 arg_max_score              = 0; // [allow setting of this field?]
///     Uint4 arg_min_score              = 0; // [allow setting of this field?]
///     Uint4 arg_set_max_score          = 0; // [allow setting of this field?]
///     Uint4 arg_set_min_score          = 0; // [allow setting of this field?]
/// ```
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct MaskerParams {
    /// `-window`; 0 = unit size + 4.
    pub window_size: u8,
    /// T_extend (`-t_extend`); 0 = the file's value.
    pub textend: u32,
    /// T_threshold (`-t_thres`); 0 = the file's value.
    pub cutoff_score: u32,
    /// T_high (`-t_high`); 0 = the file's value.
    pub max_score: u32,
    /// T_low (`-t_low`); 0 = the file's value.
    pub min_score: u32,
    /// The count used above T_high (`-set_t_high`); 0 = T_high.
    pub set_max_score: u32,
    /// The count used below T_low (`-set_t_low`); 0 = (T_low + 1) / 2.
    pub set_min_score: u32,
}

/// The errors of the counts readers and the masker. The messages are NCBI's exception texts;
/// the caller adds NCBI's framing (the exception chain or the BLAST error line).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WinMaskError {
    /// `CSeqMaskerIstatFactory::eOpen`: "could not open <name>"
    /// (`seq_masker_istat_factory.cpp:171`).
    Open { name: String },
    /// An exception of a reader (`eFormat`, `eBadParam`, `eBadHashParam`); the factory
    /// rethrows it as "could not create a unit counts container"
    /// (`seq_masker_istat_factory.cpp:358-359`). `warnings` are the warnings NCBI posted
    /// before it.
    Format {
        message: String,
        warnings: Vec<String>,
    },
    /// The factory's own `eCreateFail`: "unrecognized unit counts format"
    /// (`seq_masker_istat_factory.cpp:340`) or a `std::exception` it caught (`:361-366`).
    CreateFail { message: String },
    /// A format NCBI reads but LOSAT does not (the ascii and optimized ascii counts).
    UnsupportedFormat(StatType),
    /// `CSeqMaskerException::eValidation` (`seq_masker.cpp:121-126`).
    WindowSize { window: u8, unit: u8 },
    /// A cache bit array of 0 words (units of 1 or 2 letters) that the cache boost would read:
    /// NCBI reads past its `new Uint4[0]` (`seq_masker_cache_boost.cpp:40-48`).
    EmptyCacheBitArray { unit: u8 },
    /// `CSeqMaskerUsetHash::eBadIndex` (`seq_masker_uset_hash.cpp:95-99`).
    BadIndex { key: u32, value: u32 },
}

impl WinMaskError {
    pub(crate) fn format(message: &str) -> Self {
        Self::Format {
            message: message.to_string(),
            warnings: Vec::new(),
        }
    }
}

impl fmt::Display for WinMaskError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Open { name } => write!(f, "could not open {name}"),
            Self::Format { message, .. } => f.write_str(message),
            Self::CreateFail { message } => f.write_str(message),
            Self::UnsupportedFormat(t) => write!(
                f,
                "windowmasker unit counts in the {} format are not supported by LOSAT",
                match t {
                    StatType::Ascii => "ascii",
                    StatType::OAscii => "oascii",
                    StatType::Binary => "binary",
                    StatType::OBinary => "obinary",
                    StatType::Unknown => "unknown",
                }
            ),
            // NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker.cpp:121-125
            // ```c
            //         os << "window size (" << (int)window_size << ") "
            //               "must be greater or equal to unit size (" <<
            //               (int)ustat->UnitSize() << ")";
            // ```
            Self::WindowSize { window, unit } => write!(
                f,
                "window size ({window}) must be greater or equal to unit size ({unit})"
            ),
            Self::EmptyCacheBitArray { unit } => write!(
                f,
                "the cache bit array of units of {unit} letters has no words to read"
            ),
            // NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_uset_hash.cpp:95-99
            // ```c
            //             r << "bad index at key " << hash.first
            //               << " : " << htp[hash.first];
            // ```
            Self::BadIndex { key, value } => write!(f, "bad index at key {key} : {value}"),
        }
    }
}

impl std::error::Error for WinMaskError {}

/// The masker (`CSeqMasker`) with the "mean" trigger and no merge pass.
pub struct SeqMasker {
    ustat: UnitCounts,
    window_size: u8,
}

/// The window advances one letter at a time and every unit of a window is used
/// (`arg_window_step = 1`, `arg_unit_step = 1` in blastn and the application).
const WINDOW_STEP: u32 = 1;
const UNIT_STEP: u8 = 1;

impl SeqMasker {
    /// `s_BuildSeqMasker( lstat )` with `params` in place of its zeros.
    pub fn open(path: &Path, params: &MaskerParams) -> Result<Self, WinMaskError> {
        Self::new(UnitCounts::open(path, params)?, params)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker.cpp:119-129
    /// ```c
    ///     if( window_size == 0 ) window_size = ustat->UnitSize() + 4;
    ///
    ///     if( window_size < ustat->UnitSize() ) {
    ///         std::ostringstream os;
    ///         os << "window size (" << (int)window_size << ") "
    ///               "must be greater or equal to unit size (" <<
    ///               (int)ustat->UnitSize() << ")";
    ///         NCBI_THROW( CSeqMaskerException, eValidation, os.str() );
    ///     }
    ///
    ///     trigger_score = score = new CSeqMaskerScoreMean( ustat );
    /// ```
    /// `ustat` must have been read with the same `params` (the thresholds are set there).
    pub fn new(ustat: UnitCounts, params: &MaskerParams) -> Result<Self, WinMaskError> {
        let unit = ustat.unit_size();
        let window_size = if params.window_size == 0 {
            unit + 4
        } else {
            params.window_size
        };
        if window_size < unit {
            return Err(WinMaskError::WindowSize {
                window: window_size,
                unit,
            });
        }
        Ok(Self { ustat, window_size })
    }

    pub fn unit_counts(&self) -> &UnitCounts {
        &self.ustat
    }

    pub fn window_size(&self) -> u8 {
        self.window_size
    }

    /// Masks `data` (iupacna, upper case, plus strand) and returns the intervals in order, both
    /// ends included (`CSeqMasker::TMaskList`).
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker.cpp:170-262
    /// ```c
    ///     ustat->total_ = 0;
    ///     unique_ptr<TMaskList> mask(new TMaskList);
    ///     Uint4 cutoff_score = ustat->get_threshold();
    ///     Uint4 textend = ustat->get_textend();
    /// ...
    ///     CSeqMaskerWindow & window = *window_ptr;
    ///     score->SetWindow( window );
    /// ...
    ///     Uint4 start = 0, end = 0, cend = 0;
    ///     Uint4 limit = textend;
    ///     const CSeqMaskerIstat::optimization_data * od
    ///         = ustat->get_optimization_data();
    ///
    ///     CSeqMaskerCacheBoost booster( window, od );
    ///
    ///     while( window )
    ///     {
    ///         Uint4 ts = (*trigger_score)();
    ///         Uint4 s = (*score)();
    ///         Uint4 adv = window_step;
    ///
    ///         if( s < limit )
    ///         {
    ///             if( end > start )
    ///             {
    ///                 if( window.Start() > cend )
    ///                 {
    ///                     mask->push_back( TMaskedInterval( start, end ) );
    ///                     start = end = cend = 0;
    ///                 }
    ///             }
    ///
    ///             if( od != 0 && od->cba_ != 0 )
    ///             {
    ///                 adv = window.Start();
    ///
    ///                 if( !booster.Check() )
    ///                     break;
    ///
    ///                 adv = window_step*( 1 + window.Start() - adv );
    ///             }
    ///         }
    ///         else if( ts < cutoff_score )
    ///         {
    ///             if( end  > start )
    ///             {
    ///                 if( window.Start() > cend + 1 )
    ///                 {
    ///                     mask->push_back( TMaskedInterval( start, end ) );
    ///                     start = end = cend = 0;
    ///                 }
    ///                 else cend = window.End();
    ///             }
    ///         }
    ///         else
    ///         {
    ///             if( end > start )
    ///             {
    ///                 if( window.Start() > cend + 1 )
    ///                 {
    ///                     mask->push_back( TMaskedInterval( start, end ) );
    ///                     start = window.Start();
    ///                 }
    ///             }
    ///             else start = window.Start();
    ///
    ///             cend = end = window.End();
    ///         }
    ///
    ///
    ///         if( adv == window_step )
    ///             ++window;
    ///
    ///         score->PostAdvance( adv );
    ///     }
    ///
    ///     if( end > start )
    ///         mask->push_back( TMaskedInterval( start, end ) );
    /// ```
    /// With the "mean" trigger `trigger_score` and `score` are the same object, so `ts == s`.
    /// `merge_pass` is false (`winmask_filter.cpp:110`, `win_mask_config.cpp:336`).
    pub fn mask(&self, data: &[u8]) -> Result<Vec<(u32, u32)>, WinMaskError> {
        assert!(
            u32::try_from(data.len()).is_ok(),
            "WindowMasker positions are 32-bit (TSeqPos)"
        );
        let mut mask = Vec::new();
        let cutoff_score = self.ustat.threshold();
        let limit = self.ustat.textend();
        let mut window = Window::new(
            data,
            self.ustat.unit_size(),
            self.window_size,
            WINDOW_STEP,
            UNIT_STEP,
            0,
            data.len() as u32,
        );
        let mut score = ScoreMean::new(&window, &self.ustat)?;
        let (mut start, mut end, mut cend) = (0u32, 0u32, 0u32);
        let cba = self.ustat.cache_bit_array();
        let mut booster = CacheBoost::new(window.num_units());

        while window.state() {
            let s = score.value();
            let ts = s;
            let mut adv = WINDOW_STEP;

            if s < limit {
                if end > start && window.start() > cend {
                    mask.push((start, end));
                    start = 0;
                    end = 0;
                    cend = 0;
                }

                if let Some(cba) = cba {
                    if cba.is_empty() {
                        return Err(WinMaskError::EmptyCacheBitArray {
                            unit: self.ustat.unit_size(),
                        });
                    }
                    adv = window.start();
                    if !booster.check(&mut window, cba) {
                        break;
                    }
                    adv = WINDOW_STEP
                        .wrapping_mul(1u32.wrapping_add(window.start()).wrapping_sub(adv));
                }
            } else if ts < cutoff_score {
                if end > start {
                    if window.start() > cend.wrapping_add(1) {
                        mask.push((start, end));
                        start = 0;
                        end = 0;
                        cend = 0;
                    } else {
                        cend = window.end();
                    }
                }
            } else {
                if end > start {
                    if window.start() > cend.wrapping_add(1) {
                        mask.push((start, end));
                        start = window.start();
                    }
                } else {
                    start = window.start();
                }
                end = window.end();
                cend = end;
            }

            if adv == WINDOW_STEP {
                window.advance(WINDOW_STEP);
            }

            score.post_advance(adv, &window, &self.ustat)?;
        }

        if end > start {
            mask.push((start, end));
        }
        Ok(mask)
    }
}

/// `CSeqMaskerCacheBoost`: skips the windows whose units all have a 0 bit in the cache bit
/// array (counts below the T_extend the file was written with). With the file's thresholds
/// (blastn) every skipped window scores below T_extend and the result is that of the plain
/// loop; with a smaller `-t_extend` or a larger `-set_t_low` of the application a skipped window
/// can score above it, so the skip is ported as it is, not left out.
struct CacheBoost {
    last_checked: u32,
    nu: usize,
}

impl CacheBoost {
    /// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_cache_boost.hpp:55-58
    /// ```c
    ///         CSeqMaskerCacheBoost( CSeqMaskerWindow & window,
    ///                               const CSeqMaskerIstat::optimization_data * od )
    ///             : window_( window ), od_( od ), last_checked_( 0 )
    ///         { nu_ = window_.NumUnits(); }
    /// ```
    fn new(nu: usize) -> Self {
        Self {
            last_checked: 0,
            nu,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_cache_boost.cpp:40-48
    /// ```c
    /// inline Uint1 CSeqMaskerCacheBoost::bit_at( TUnit pos ) const
    /// {
    ///     pos /= od_->divisor_;
    ///     TSeqPos word = pos/(8*sizeof( Uint4 ));
    ///     TSeqPos bit = pos%(8*sizeof( Uint4 ));
    ///     Uint1 res = (((od_->cba_[word])>>bit)&0x1) == 0 ? 0 : 1;
    ///     return res;
    /// }
    /// ```
    /// `divisor_` is `32/(8*sizeof( Uint4 ))` = 1 (`seq_masker_istat.hpp:67-68`,
    /// `seq_masker_istat_obinary.cpp:185`).
    #[inline]
    fn bit_at(cba: &[u32], pos: u32) -> bool {
        (cba[(pos / 32) as usize] >> (pos % 32)) & 0x1 != 0
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_cache_boost.cpp:50-56
    /// ```c
    ///     for( unsigned int i = 0; i < nu_; ++i )
    ///         if( bit_at( window_[i] ) != 0 )
    ///             return false;
    ///
    ///     return true;
    /// ```
    fn full_check(&self, window: &Window<'_>, cba: &[u32]) -> bool {
        (0..self.nu).all(|i| !Self::bit_at(cba, window.unit(i)))
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_cache_boost.cpp:60-79
    /// ```c
    ///     if( od_ == 0 || od_->cba_ == 0 )
    ///         return true;
    ///
    ///     while( window_ )
    ///     {
    ///         if( last_checked_ + 1 != window_.End() )
    ///         {
    ///             if( !full_check() )
    ///                 break;
    ///         }
    ///         else if( bit_at( window_[nu_-1] ) != 0 )
    ///                 break;
    ///
    ///         last_checked_ = window_.End();
    ///         ++window_;
    ///     }
    ///
    ///     return bool( window_ );
    /// ```
    fn check(&mut self, window: &mut Window<'_>, cba: &[u32]) -> bool {
        while window.state() {
            if self.last_checked.wrapping_add(1) != window.end() {
                if !self.full_check(window, cba) {
                    break;
                }
            } else if Self::bit_at(cba, window.unit(self.nu - 1)) {
                break;
            }
            self.last_checked = window.end();
            window.advance(window.step());
        }
        window.state()
    }
}
