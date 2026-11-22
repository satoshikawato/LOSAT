//! The unit counts generator of the windowmasker application's `-mk_counts`
//! (`CWinMaskCountsGenerator::operator()`): the unit size, the split of the units into prefixes
//! that `-mem` allows, the two passes over the input (the first finds the counts at the
//! percentages of `-t_*_pct`, the second writes the units from `T_low` with their counts up to
//! `T_high`), and the statistics that the counts file keeps (comments, `##pct:` lines,
//! parameters).
//!
//! The input is read through a `GenomeSource` (the tool's readers, which post their warnings as
//! they read). NCBI reads the whole input once per prefix and pass; where both passes have one
//! prefix the counts of the two passes are the same, so LOSAT counts once and the source repeats
//! the posts of that read for the second pass (`GenomeSource::repeat_posts`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_gen_counts.cpp:226-550

use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::sam::push_g;

use super::istat::reverse_complement;
use super::ostat::CountsWriter;

/// The settings of the generator (`CWinMaskConfig` of `-mk_counts`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GenCountsSettings {
    /// `-unit`, 0 when not given.
    pub unit_size: u8,
    /// `-genome_size` as NCBI's `Uint8` (a negative value wraps), 0 when not given.
    pub genome_size: u64,
    /// `-mem` in megabytes.
    pub mem: u32,
    /// `-t_low`, 0 when not given.
    pub min_count: u32,
    /// `-t_high`, 0 when not given.
    pub max_count: u32,
    /// The percentages of T_low, T_extend, T_threshold and T_high: 90, 99, 99.5, 99.8, each
    /// replaced by its `-t_*_pct` when the counts go to a named file.
    pub th: [f64; 4],
}

/// The default percentages (`CWinMaskConfig`'s `th`, read with `atof`).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_config.cpp:353
pub const DEFAULT_TH: [f64; 4] = [90.0, 99.0, 99.5, 99.8];

/// The input of the generator: every sequence that the run considers, file after file.
pub trait GenomeSource {
    type Error;
    /// Reads the whole input once (`CInputBioseq_CI` over each file), giving the letters of each
    /// considered sequence (IUPACna) to `f` and posting what NCBI's readers post.
    fn read(&mut self, f: &mut dyn FnMut(&[u8])) -> Result<(), Self::Error>;
    /// Posts again what the last `read` posted: a read of the same input that NCBI makes and
    /// LOSAT does not repeat.
    fn repeat_posts(&mut self) -> Result<(), Self::Error>;
    /// `LOG_POST` of `line`.
    fn log(&mut self, line: &str) -> Result<(), Self::Error>;
}

/// Why the generator wrote no counts.
#[derive(Debug)]
pub enum GenError<E> {
    /// An error of the source.
    Source(E),
    /// `GenCountsException::eNullGenome`: the input has no letters.
    EmptyGenome(NcbiException),
    /// What LOSAT does not reproduce: the reason, a noun phrase.
    Refused(String),
}

/// The largest count that the statistics tell apart (`max_count`).
const MAX_COUNT: usize = 500;

/// `CWinMaskCountsGenerator::operator()` after the duplicates check, with `source` as the input
/// and `writer` as the counts file (not yet finalized).
pub fn generate<S: GenomeSource>(
    settings: &GenCountsSettings,
    source: &mut S,
    writer: &mut CountsWriter,
) -> Result<(), GenError<S::Error>> {
    let mut unit_size = settings.unit_size;
    if unit_size == 0 {
        let mut genome = settings.genome_size;
        if genome == 0 {
            source
                .log("computing the genome length")
                .map_err(GenError::Source)?;
            let mut total = 0u64;
            source
                .read(&mut |letters| total = total.wrapping_add(letters.len() as u64))
                .map_err(GenError::Source)?;
            if total == 0 {
                return Err(GenError::EmptyGenome(empty_genome()));
            }
            genome = total;
        }
        unit_size = unit_size_for(genome);
    }
    let prefix_size = unit_size - suffix_size(unit_size, settings.mem);
    writer.set_unit_size(unit_size);
    let mut stats = Stats::new(settings);
    let mut counter = UnitCounter::new(unit_size, prefix_size)?;
    let prefixes = 1u32 << (2 * prefix_size);
    source.log("pass 1").map_err(GenError::Source)?;
    let shared = prefix_size == 0 && !stats.no_extra_pass;
    for prefix in 0..prefixes {
        counter.count(prefix, source).map_err(GenError::Source)?;
        stats.tally(&counter, stats.no_extra_pass, writer);
    }
    if stats.total_ecodes == 0 {
        return Err(GenError::Refused(
            "-mk_counts of an input without a unit of unambiguous letters, whose percentages NCBI BLAST+ prints as the C runtime prints 0/0 (-nan on x86-64 Linux)".to_string(),
        ));
    }
    let index = stats.first_pass_statistics(writer);
    if !stats.no_extra_pass {
        stats.prepare_second_pass(&index);
        source.log("pass 2").map_err(GenError::Source)?;
        for prefix in 0..prefixes {
            if shared {
                source.repeat_posts().map_err(GenError::Source)?;
            } else {
                counter.count(prefix, source).map_err(GenError::Source)?;
            }
            stats.tally(&counter, true, writer);
        }
        stats.second_pass_statistics(writer);
    }
    drop(counter);
    writer.set_comment(b"");
    for (th, at) in stats.th.iter().zip(index) {
        let mut text = b" ".to_vec();
        push_g(&mut text, *th);
        text.extend_from_slice(format!("%% threshold at {at}").as_bytes());
        writer.set_comment(&text);
    }
    for (i, value) in index.into_iter().enumerate() {
        writer.set_param(i, value);
    }
    Ok(())
}

/// `GenCountsException::eNullGenome` with no message.
fn empty_genome() -> NcbiException {
    NcbiException::new(
        ThrowSite {
            file: "src/algo/winmask/win_mask_gen_counts.cpp",
            line: 266,
            kind: "GenCountsException::empty genome",
            function: Some("ncbi::CWinMaskCountsGenerator::operator()"),
        },
        "",
    )
}

/// The unit size of a genome of `genome` letters: one more than the largest size below 16 whose
/// number of units, `genome / 4^size`, is at least 5 (1 when none is).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_gen_counts.cpp:270-277
pub fn unit_size_for(genome: u64) -> u8 {
    let size = (1..=15u8)
        .rev()
        .find(|&size| (genome >> (2 * u32::from(size))) >= 5)
        .unwrap_or(0);
    size + 1
}

/// The suffix size: the largest size, at most `unit_size`, whose `4^size` counters of 4 bytes
/// fit `mem` megabytes (the rest of the unit is the prefix, one pass over the input per prefix).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_gen_counts.cpp:282-293
pub fn suffix_size(unit_size: u8, mem: u32) -> u8 {
    let n_units = u64::from(mem) * 1024 * 1024 / 4;
    let mut suffix = unit_size;
    while suffix > 0 && (1u64 << (2 * u32::from(suffix))) > n_units {
        suffix -= 1;
    }
    suffix
}

/// The code of each IUPACna letter: 0 to 3 for `ACGT` in either case, 4 for any other letter
/// (an ambiguous letter starts the units again).
const LETTER: [u8; 256] = {
    let mut table = [4u8; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

/// The counts of the units of one prefix (`process`'s `counts`).
struct UnitCounter {
    unit_size: u8,
    suffix_size: u8,
    unit_mask: u32,
    prefix_mask: u32,
    suffix_mask: u32,
    /// The prefix counted last, shifted to its place in the unit.
    prefix: u32,
    counts: Vec<u32>,
}

impl UnitCounter {
    /// The counters of units of `unit_size` letters whose first `prefix_size` letters are a
    /// prefix: `4^(unit_size - prefix_size)` of them, as NCBI's `vector< Uint4 >`. An
    /// allocation that fails is NCBI's `std::bad_alloc`, which depends on the memory of the
    /// host: refused.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_gen_counts.cpp:409-426
    fn new<E>(unit_size: u8, prefix_size: u8) -> Result<Self, GenError<E>> {
        let suffix_size = unit_size - prefix_size;
        let size = 1usize << (2 * u32::from(suffix_size));
        let mut counts = Vec::new();
        if counts.try_reserve_exact(size).is_err() {
            return Err(GenError::Refused(format!(
                "-mem large enough for {size} counters of units of {unit_size} letters, which cannot be allocated here (NCBI BLAST+ then fails for want of memory, which depends on the host)"
            )));
        }
        counts.resize(size, 0);
        let full = |bits: u32| {
            if bits >= 32 {
                u32::MAX
            } else {
                (1u32 << bits) - 1
            }
        };
        let (unit_bits, suffix_bits) = (2 * u32::from(unit_size), 2 * u32::from(suffix_size));
        Ok(UnitCounter {
            unit_size,
            suffix_size,
            unit_mask: full(unit_bits),
            prefix_mask: if suffix_bits >= 32 {
                0
            } else {
                full(2 * u32::from(prefix_size)).wrapping_shl(suffix_bits)
            },
            suffix_mask: full(suffix_bits),
            prefix: 0,
            counts,
        })
    }

    /// Counts the units of `prefix` in every sequence of `source`: each unit of `unit_size`
    /// unambiguous letters is counted under the smaller of its value and the value of its
    /// reverse complement (a palindrome twice), the counts stopping at `0xffffffff`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_gen_counts.cpp:452-511
    fn count<S: GenomeSource + ?Sized>(
        &mut self,
        prefix: u32,
        source: &mut S,
    ) -> Result<(), S::Error> {
        let suffix_bits = 2 * u32::from(self.suffix_size);
        self.prefix = prefix.checked_shl(suffix_bits).unwrap_or(0);
        self.counts.fill(0);
        let (unit_size, unit_mask) = (u32::from(self.unit_size), self.unit_mask);
        let (prefix_mask, suffix_mask, prefix) = (self.prefix_mask, self.suffix_mask, self.prefix);
        let top = 2 * (unit_size - 1);
        let counts = &mut self.counts;
        source.read(&mut |letters| {
            let (mut unit, mut runit, mut count) = (0u32, 0u32, 0u32);
            for &c in letters {
                let letter = u32::from(LETTER[usize::from(c)]);
                if letter > 3 {
                    (unit, runit, count) = (0, 0, 0);
                    continue;
                }
                unit = ((unit << 2) & unit_mask).wrapping_add(letter);
                runit = (runit >> 2) | ((3 - letter) << top);
                if count >= unit_size - 1 {
                    if unit <= runit && unit & prefix_mask == prefix {
                        let c = &mut counts[(unit & suffix_mask) as usize];
                        *c = c.saturating_add(1);
                    }
                    if runit <= unit && runit & prefix_mask == prefix {
                        let c = &mut counts[(runit & suffix_mask) as usize];
                        *c = c.saturating_add(1);
                    }
                }
                count += 1;
            }
        })
    }
}

/// The statistics of the passes (`total_ecodes`, `score_counts`, the thresholds).
struct Stats {
    th: [f64; 4],
    min_count: u32,
    t_high: u32,
    has_min_count: bool,
    no_extra_pass: bool,
    total_ecodes: u32,
    score_counts: Vec<u32>,
}

impl Stats {
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_gen_counts.cpp:178-224
    fn new(settings: &GenCountsSettings) -> Self {
        Stats {
            th: settings.th,
            min_count: settings.min_count.max(1),
            t_high: settings.max_count,
            has_min_count: settings.min_count != 0,
            no_extra_pass: settings.min_count != 0 && settings.max_count != 0,
            total_ecodes: 0,
            score_counts: vec![0; MAX_COUNT],
        }
    }

    /// The end of `process` for the counts of one prefix: the number of units seen (a
    /// palindrome once, another unit twice, for itself and its reverse complement), the number
    /// of units of each count from `min_count` (500 and above together), and with `output` the
    /// units from `min_count` with their counts up to `t_high`, in ascending order.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_gen_counts.cpp:527-549
    fn tally(&mut self, counter: &UnitCounter, output: bool, writer: &mut CountsWriter) {
        let unit_size = counter.unit_size;
        for (i, &count) in counter.counts.iter().enumerate() {
            if count == 0 {
                continue;
            }
            let unit = counter.prefix.wrapping_add(i as u32);
            let weight = if unit == reverse_complement(unit, unit_size) {
                1
            } else {
                2
            };
            self.total_ecodes = self.total_ecodes.wrapping_add(weight);
            if count >= self.min_count {
                let slot = &mut self.score_counts[(count as usize).min(MAX_COUNT) - 1];
                *slot = slot.wrapping_add(weight);
                if output {
                    writer.set_unit_count(unit, count.min(self.t_high));
                }
            }
        }
    }

    /// The cumulative numbers of units by count, and the units of counts above 500 added to
    /// each (`offset`).
    fn cumulate(&mut self) -> u32 {
        for i in 1..MAX_COUNT {
            self.score_counts[i] = self.score_counts[i].wrapping_add(self.score_counts[i - 1]);
        }
        self.total_ecodes
            .wrapping_sub(self.score_counts[MAX_COUNT - 1])
    }

    /// The percentage of units with count `i` or less (`current`).
    fn current(&self, i: usize, offset: u32) -> f64 {
        let at = self.score_counts[i - 1].wrapping_add(offset);
        100.0 * (f64::from(at) / f64::from(self.total_ecodes))
    }

    /// The comment of count `i` with the percentage `current`.
    fn count_comment(&self, i: usize, offset: u32, current: f64) -> Vec<u8> {
        let at = self.score_counts[i - 1].wrapping_add(offset);
        let mut text = format!(" {i}\t{at}\t").into_bytes();
        push_g(&mut text, current);
        text
    }

    /// After the first pass: the counts at which the percentages of `th` are reached (`index`,
    /// 0 for one never crossed) and, when the first pass was the only one, the comments and
    /// percentages of its counts.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_gen_counts.cpp:305-338
    fn first_pass_statistics(&mut self, writer: &mut CountsWriter) -> [u32; 4] {
        let offset = self.cumulate();
        let mut index = [0u32; 4];
        let mut previous = 0.0f64;
        if self.no_extra_pass {
            writer.set_comment(format!(" {} ecodes", self.total_ecodes).as_bytes());
        }
        let mut count_map = vec![0.0f64; MAX_COUNT + 1];
        for i in 1..=MAX_COUNT {
            let current = self.current(i, offset);
            if self.no_extra_pass {
                writer.set_comment(&self.count_comment(i, offset, current));
                count_map[i] = current;
            }
            for (slot, &th) in index.iter_mut().zip(&self.th) {
                if previous < th && current >= th {
                    *slot = i as u32;
                }
            }
            previous = current;
        }
        writer.set_count_map(count_map);
        index
    }

    /// The settings of the second pass from the first's `index`.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_gen_counts.cpp:341-361
    fn prepare_second_pass(&mut self, index: &[u32; 4]) {
        self.total_ecodes = 0;
        if !self.has_min_count {
            self.min_count = index[0];
        }
        if self.t_high == 0 {
            self.t_high = index[3];
        }
        if self.min_count == 0 {
            self.min_count = 1;
        }
        self.score_counts.fill(0);
    }

    /// After the second pass: its comments and the percentages of its counts.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_gen_counts.cpp:366-389
    fn second_pass_statistics(&mut self, writer: &mut CountsWriter) {
        let offset = self.cumulate();
        writer.set_comment(format!(" {} ecodes", self.total_ecodes).as_bytes());
        let mut count_map = vec![0.0f64; MAX_COUNT + 1];
        for (i, slot) in count_map.iter_mut().enumerate().skip(1) {
            let current = self.current(i, offset);
            writer.set_comment(&self.count_comment(i, offset, current));
            *slot = current;
        }
        writer.set_count_map(count_map);
    }
}

#[cfg(test)]
#[path = "gen_counts_tests.rs"]
mod tests;
