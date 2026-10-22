//! The window of units that slides along the sequence (`CSeqMaskerWindow`).

/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_window.cpp:46, 63-69
/// ```c
/// Uint1 CSeqMaskerWindow::LOOKUP[kMax_UI1];
/// ...
///         LOOKUP[unsigned('A')] = 1;
///         LOOKUP[unsigned('C')] = 2;
///         LOOKUP[unsigned('G')] = 3;
///         LOOKUP[unsigned('T')] = 4;
/// ```
/// Every other letter is 0 (an ambiguity: the window is filled again after it). The data is
/// iupacna, upper case (`CBioseq_Handle::eCoding_Iupac`).
static LOOKUP: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = 1;
    table[b'C' as usize] = 2;
    table[b'G' as usize] = 3;
    table[b'T' as usize] = 4;
    table
};

pub(super) struct Window<'a> {
    data: &'a [u8],
    state: bool,
    unit_size: u8,
    unit_step: u8,
    window_size: u8,
    window_step: u32,
    start: u32,
    end: u32,
    first_unit: usize,
    num_units: usize,
    units: Vec<u32>,
    unit_mask: u32,
    winend: u32,
}

impl<'a> Window<'a> {
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_window.cpp:49-58, 78-86
    /// ```c
    /// CSeqMaskerWindow::CSeqMaskerWindow( const CSeqVector & arg_data,
    ///                                     Uint1 arg_unit_size,
    ///                                     Uint1 arg_window_size,
    ///                                     Uint4 arg_window_step,
    ///                                     Uint1 arg_unit_step,
    ///                                     Uint4 winstart,
    ///                                     Uint4 arg_winend )
    /// ...
    ///     units.resize( NumUnits(), 0 );
    ///     unit_mask = (unit_size == 16) ? 0xFFFFFFFF : (1ULL << (unit_size << 1)) - 1;
    ///
    ///     if( winend == 0 )
    ///         winend = data.size();
    ///
    ///     FillWindow( winstart );
    /// ```
    /// `num_units` is `NumUnits()` (`seq_masker_window.hpp`:
    /// `(window_size - unit_size)/unit_step + 1`), 1 to 255 for a window of at least one unit.
    pub(super) fn new(
        data: &'a [u8],
        unit_size: u8,
        window_size: u8,
        window_step: u32,
        unit_step: u8,
        winstart: u32,
        winend: u32,
    ) -> Self {
        let num_units = units_in_window(unit_size, window_size, unit_step) as usize;
        let mut window = Self {
            data,
            state: false,
            unit_size,
            unit_step,
            window_size,
            window_step,
            start: 0,
            end: 0,
            first_unit: 0,
            num_units,
            units: vec![0; num_units],
            unit_mask: if unit_size == 16 {
                0xFFFF_FFFF
            } else {
                ((1u64 << (unit_size << 1)) - 1) as u32
            },
            winend: if winend == 0 {
                data.len() as u32
            } else {
                winend
            },
        };
        window.fill_window(winstart);
        window
    }

    /// `operator bool()`: false once the end of the sequence was reached.
    #[inline]
    pub(super) fn state(&self) -> bool {
        self.state
    }

    #[inline]
    pub(super) fn start(&self) -> u32 {
        self.start
    }

    #[inline]
    pub(super) fn end(&self) -> u32 {
        self.end
    }

    #[inline]
    pub(super) fn step(&self) -> u32 {
        self.window_step
    }

    #[inline]
    pub(super) fn unit_step(&self) -> u8 {
        self.unit_step
    }

    #[inline]
    pub(super) fn num_units(&self) -> usize {
        self.num_units
    }

    /// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_window.hpp:128-133
    /// ```c
    ///     TUnit operator[]( Uint1 index ) const
    ///     {
    ///         return first_unit + index < NumUnits() ?
    ///             units[first_unit + index] :
    ///             units[first_unit + index - NumUnits()];
    ///     }
    /// ```
    #[inline]
    pub(super) fn unit(&self, index: usize) -> u32 {
        let i = self.first_unit + index;
        if i < self.num_units {
            self.units[i]
        } else {
            self.units[i - self.num_units]
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_window.cpp:93-128
    /// ```c
    ///     if( step >= window_size || unit_step > 1 )
    ///     {
    ///         FillWindow( start + step );
    ///         return;
    ///     }
    ///
    ///     Uint1 num_units = NumUnits();
    ///     Uint1 last_unit = first_unit ? first_unit - 1 : num_units - 1;
    ///     Uint4 unit = units[last_unit];
    ///     Uint4 iter = 0;
    ///
    ///     for( ; ++end < winend && iter < step ; ++iter, ++start )
    ///     {
    ///         Uint1 letter = LOOKUP[unsigned(data[end])];
    ///
    ///         if( !(letter--) )
    ///         {
    ///             FillWindow( end );
    ///             return;
    ///         }
    ///
    ///         unit = ((unit<<2)&unit_mask) + letter;
    ///
    ///         if( ++first_unit == num_units ) first_unit = 0;
    ///
    ///         if( ++last_unit == num_units ) last_unit = 0;
    ///
    ///         units[last_unit] = unit;
    ///     }
    ///
    ///     --end;
    ///
    ///     if( iter != step ) state = false;
    /// ```
    pub(super) fn advance(&mut self, step: u32) {
        if step >= self.window_size as u32 || self.unit_step > 1 {
            self.fill_window(self.start.wrapping_add(step));
            return;
        }

        let num_units = self.num_units;
        let mut last_unit = if self.first_unit != 0 {
            self.first_unit - 1
        } else {
            num_units - 1
        };
        let mut unit = self.units[last_unit];
        let mut iter = 0u32;

        loop {
            // `++end < winend && iter < step`: `++end` runs before every test.
            self.end = self.end.wrapping_add(1);
            if !(self.end < self.winend && iter < step) {
                break;
            }
            let letter = LOOKUP[self.data[self.end as usize] as usize];
            if letter == 0 {
                self.fill_window(self.end);
                return;
            }
            unit = ((unit << 2) & self.unit_mask) + (letter - 1) as u32;
            self.first_unit += 1;
            if self.first_unit == num_units {
                self.first_unit = 0;
            }
            last_unit += 1;
            if last_unit == num_units {
                last_unit = 0;
            }
            self.units[last_unit] = unit;
            iter += 1;
            self.start = self.start.wrapping_add(1);
        }

        self.end = self.end.wrapping_sub(1);

        if iter != step {
            self.state = false;
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_window.cpp:131-157
    /// ```c
    ///     first_unit = 0;
    ///     TUnit unit = 0;
    ///     Int4 iter = 0;
    ///     end = winstart;
    ///
    ///     for( ; iter < window_size && end < data.size(); ++iter, ++end )
    ///     {
    ///         Uint1 letter = LOOKUP[unsigned(data[end])];
    ///
    ///         if( !(letter--) )
    ///         {
    ///             iter = -1;
    ///             continue;
    ///         }
    ///
    ///         unit = ((unit<<2)&unit_mask) + letter;
    ///
    ///         if( iter >= unit_size - 1 )
    ///             if( !((iter + 1 - unit_size)%unit_step) )
    ///                 units[(iter + 1- unit_size)/unit_step] = unit;
    ///     }
    ///
    ///     start = (end--) - window_size;
    ///     state = (iter == window_size);
    /// ```
    fn fill_window(&mut self, winstart: u32) {
        self.first_unit = 0;
        let mut unit: u32 = 0;
        let mut iter: i32 = 0;
        self.end = winstart;
        let window_size = self.window_size as i32;
        let unit_size = self.unit_size as i32;
        let unit_step = self.unit_step as i32;

        while iter < window_size && (self.end as usize) < self.data.len() {
            let letter = LOOKUP[self.data[self.end as usize] as usize];
            if letter == 0 {
                iter = -1;
            } else {
                unit = ((unit << 2) & self.unit_mask) + (letter - 1) as u32;
                if iter >= unit_size - 1 && (iter + 1 - unit_size) % unit_step == 0 {
                    self.units[((iter + 1 - unit_size) / unit_step) as usize] = unit;
                }
            }
            iter += 1;
            self.end = self.end.wrapping_add(1);
        }

        self.start = self.end.wrapping_sub(self.window_size as u32);
        self.end = self.end.wrapping_sub(1);
        self.state = iter == window_size;
    }
}

/// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_window.hpp:176
/// ```c
///     Uint1 NumUnits() const{ return (window_size - unit_size)/unit_step + 1; }
/// ```
fn units_in_window(unit_size: u8, window_size: u8, unit_step: u8) -> u8 {
    (((window_size as i32 - unit_size as i32) / unit_step as i32) + 1) as u8
}
