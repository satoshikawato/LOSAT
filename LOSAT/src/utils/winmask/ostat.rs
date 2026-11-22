//! The writers of the four unit counts formats (`CSeqMaskerOstat` and its ascii, binary,
//! optimized ascii and optimized binary kinds, as `CSeqMaskerOstatFactory::create` makes them),
//! which the windowmasker application's `-mk_counts` and `-convert` write and whose files the
//! readers of `istat` read. The optimized formats choose the parameters of their hash table as
//! `CSeqMaskerOstatOpt::doFinalize` does (`optimize`).
//!
//! A writer collects the unit size, the counts of the units in ascending order, the comments,
//! the four parameters and the percentages, and writes the whole file at the end (`finalize`),
//! as NCBI's writers do.

use std::io::{self, Write};

use crate::ncbi_exception::{NcbiException, ThrowSite};
use crate::report::sam::push_g;

use super::istat::{hash_code, reverse_complement, UnitCounts};
use super::WinMaskError;

/// A counts format; the optimized ones with the memory (in megabytes) that their hash table may
/// take (`size_requested`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CountsFormat {
    Ascii,
    Binary,
    OAscii { size_requested: u16 },
    OBinary { size_requested: u16 },
}

// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_factory.cpp:55-79
// NCBI reference (598d8ae6): c++/include/algo/winmask/win_mask_config.hpp:462-467
// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_ostat_opt.hpp:173
impl CountsFormat {
    /// The format of `-sformat` `sformat` (one of the four names its constraint allows) with
    /// `-smem` `smem`. NCBI's factory takes the format from the start of the string `sformat`
    /// followed by `smem` printed as a `Uint4`, and reads the optimized formats' memory from the
    /// rest with `atoi` into a `Uint2`: the low 16 bits of `smem` (probes of the survey of line
    /// NP item M5: `-smem 65536` is 0, `70000` is 4464, `-1` is 65535).
    pub fn from_sformat(sformat: &str, smem: i32) -> Self {
        let size_requested = smem as u16;
        match sformat {
            "ascii" => CountsFormat::Ascii,
            "binary" => CountsFormat::Binary,
            "oascii" => CountsFormat::OAscii { size_requested },
            "obinary" => CountsFormat::OBinary { size_requested },
            other => unreachable!("-sformat {other} passed its constraint"),
        }
    }

    /// The name in the file's format version line.
    fn name(self) -> &'static str {
        match self {
            CountsFormat::Ascii => "ascii",
            CountsFormat::Binary => "binary",
            CountsFormat::OAscii { .. } => "oascii",
            CountsFormat::OBinary { .. } => "obinary",
        }
    }

    /// What the writer of this format writes to a stream as soon as it is made: the binary
    /// format's stream constructor writes a word 0 before the file (so a binary file written to
    /// the standard output starts with an extra word, which NCBI's reader then rejects); a writer
    /// that opens a file by name writes nothing then.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_bin.cpp:62-67
    pub fn stream_prefix(self) -> &'static [u8] {
        match self {
            CountsFormat::Binary => &[0; 4],
            _ => &[],
        }
    }
}

/// The version of the counts algorithm that a file names in its
/// `##windowmasker-statistics-algorithm-version:` line (`CSeqMaskerVersion` of
/// `CSeqMaskerOstat::STAT_ALGO_COMPONENT_NAME`); NCBI's own is 1.0.0.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat.cpp:41-56
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StatAlgoVersion {
    pub major: i32,
    pub minor: i32,
    pub patch: i32,
}

impl Default for StatAlgoVersion {
    fn default() -> Self {
        StatAlgoVersion {
            major: 1,
            minor: 0,
            patch: 0,
        }
    }
}

/// The name of the counts algorithm's version line (`CSeqMaskerOstat::STAT_ALGO_COMPONENT_NAME`).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat.cpp:45-46
pub(super) const STAT_ALGO_COMPONENT_NAME: &str = "windowmasker-statistics-algorithm-version";

/// The names of the four parameters in their order (`CSeqMaskerOstat::PARAMS`).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_opt.cpp:49
const PARAMS: [&str; 4] = ["t_low", "t_extend", "t_threshold", "t_high"];

/// One megabyte, the unit of the optimized formats' memory.
const MB: u64 = 1024 * 1024;

/// The counts of the units: as given (`Uint4`) for the ascii and binary formats, as the
/// optimized formats keep them (`vector< Uint2 >`, the low 16 bits).
enum Counts {
    Full(Vec<u32>),
    Short(Vec<u16>),
}

/// Why `CountsWriter::finalize` wrote no file.
#[derive(Debug)]
pub enum FinalizeError {
    /// `CSeqMaskerOstatOpt::Exception::eMemory`: no hash table fits the memory of `-smem`.
    Ncbi(NcbiException),
    /// What LOSAT does not reproduce: the reason, a noun phrase.
    Refused(String),
    /// A write of the standard output that found the pipe closed (the only error that the
    /// output streams return: a named file loses what does not reach it, as NCBI's `ofstream`
    /// does).
    ClosedPipe,
}

impl From<io::Error> for FinalizeError {
    fn from(_: io::Error) -> Self {
        FinalizeError::ClosedPipe
    }
}

/// A writer of one counts file (`CSeqMaskerOstat`).
pub struct CountsWriter {
    format: CountsFormat,
    /// The `-meta` text (`##note:`), empty for none.
    metadata: Vec<u8>,
    algo_version: StatAlgoVersion,
    unit_size: u8,
    units: Vec<u32>,
    counts: Counts,
    comments: Vec<Vec<u8>>,
    params: [u32; 4],
    /// The percentage of each count from 0 (`count_map`); empty for none.
    count_map: Vec<f64>,
}

// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat.cpp:69-204
// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_opt.cpp:52-117
// NCBI reference (598d8ae6): c++/include/algo/winmask/seq_masker_ostat.hpp:142-160
impl CountsWriter {
    /// A writer of `format` with the metadata text `metadata` (`-meta`, or the input's note).
    pub fn new(format: CountsFormat, metadata: &[u8]) -> Self {
        let counts = match format {
            CountsFormat::Ascii | CountsFormat::Binary => Counts::Full(Vec::new()),
            CountsFormat::OAscii { .. } | CountsFormat::OBinary { .. } => Counts::Short(Vec::new()),
        };
        CountsWriter {
            format,
            metadata: metadata.to_vec(),
            algo_version: StatAlgoVersion::default(),
            unit_size: 0,
            units: Vec::new(),
            counts,
            comments: Vec::new(),
            params: [0; 4],
            count_map: Vec::new(),
        }
    }

    /// `setUnitSize`.
    pub fn set_unit_size(&mut self, unit_size: u8) {
        self.unit_size = unit_size;
    }

    /// `setUnitCount`: the next unit, above the units before it, and its count. The optimized
    /// formats keep the low 16 bits of the count, and grow their lists by a tenth (at least a
    /// million units) at a time, as NCBI's do.
    pub fn set_unit_count(&mut self, unit: u32, count: u32) {
        match &mut self.counts {
            Counts::Full(counts) => {
                self.units.push(unit);
                counts.push(count);
            }
            Counts::Short(counts) => {
                if self.units.len() == self.units.capacity() {
                    let grow = (self.units.len() / 10).max(1 << 20);
                    self.units.reserve_exact(grow);
                    counts.reserve_exact(grow);
                }
                self.units.push(unit);
                counts.push(count as u16);
            }
        }
    }

    /// `setComment`: a comment line of the ascii format (the others have none).
    pub fn set_comment(&mut self, text: &[u8]) {
        if self.format == CountsFormat::Ascii {
            self.comments.push(text.to_vec());
        }
    }

    /// `setParam` of the parameter `index` (`t_low`, `t_extend`, `t_threshold`, `t_high`).
    pub fn set_param(&mut self, index: usize, value: u32) {
        self.params[index] = value;
    }

    /// `SetStatAlgoVersion`.
    pub fn set_stat_algo_version(&mut self, version: StatAlgoVersion) {
        self.algo_version = version;
    }

    /// `SetMaxCount` with `max_count`, and `SetCount` of every count with its percentage: `count_map` holds the
    /// percentage of each count from 0 to `max_count`.
    pub fn set_count_map(&mut self, count_map: Vec<f64>) {
        self.count_map = count_map;
    }

    /// The format version line without its `##` (`GetStatFmtVersion().Print()`).
    fn format_version(&self) -> String {
        format!(
            "windowmasker-statistics-format-version:{} 1.1.0",
            self.format.name()
        )
    }

    /// The algorithm version line without its `##` (`fmt_gen_algo_ver.Print()`).
    fn algo_version_line(&self) -> String {
        let v = self.algo_version;
        format!(
            "{}:{}.{}.{}",
            STAT_ALGO_COMPONENT_NAME, v.major, v.minor, v.patch
        )
    }

    /// `FormatParameters`: a `t_low` of 0 is printed as 1.
    fn parameters_line(&self) -> String {
        let t_low = if self.params[0] == 0 {
            1
        } else {
            self.params[0]
        };
        format!(
            "##parameters:unit={} t_low={t_low} t_high={}",
            self.unit_size, self.params[3]
        )
    }

    /// The metadata lines of the text formats (`FormatMetaData`).
    fn text_metadata(&self) -> Vec<u8> {
        let mut text = Vec::new();
        for line in [self.format_version(), self.algo_version_line()] {
            text.extend_from_slice(b"##");
            text.extend_from_slice(line.as_bytes());
            text.push(b'\n');
        }
        text.extend_from_slice(self.parameters_line().as_bytes());
        text.push(b'\n');
        if !self.metadata.is_empty() {
            text.extend_from_slice(b"##note:");
            text.extend_from_slice(&self.metadata);
            text.push(b'\n');
        }
        for (i, &pct) in self.count_map.iter().enumerate() {
            text.extend_from_slice(format!("##pct: {i} ").as_bytes());
            push_g(&mut text, pct);
            text.push(b'\n');
        }
        text
    }

    /// The metadata of the binary formats (`WriteBinMetaData`): the length of the strings,
    /// then each string followed by a 0 byte.
    fn binary_metadata(&self) -> Vec<u8> {
        let mut strings: Vec<Vec<u8>> = vec![
            format!("##{}", self.format_version()).into_bytes(),
            format!("##{}", self.algo_version_line()).into_bytes(),
            self.parameters_line().into_bytes(),
        ];
        if !self.metadata.is_empty() {
            let mut note = b"##note:".to_vec();
            note.extend_from_slice(&self.metadata);
            strings.push(note);
        }
        if !self.count_map.is_empty() {
            let mut pct = format!("##pct: {}", self.count_map.len() - 1).into_bytes();
            for &value in &self.count_map {
                pct.push(b' ');
                push_g(&mut pct, value);
            }
            strings.push(pct);
        }
        let size: usize = strings.iter().map(|s| s.len() + 1).sum();
        let mut bytes = (size as u32).to_le_bytes().to_vec();
        for string in strings {
            bytes.extend_from_slice(&string);
            bytes.push(0);
        }
        bytes
    }

    /// `finalize`: writes the file to `out`. The optimized formats first post `optimizing the
    /// data structure` to `log` and choose their hash table (`optimize`).
    pub fn finalize(
        mut self,
        out: &mut dyn Write,
        log: &mut dyn FnMut(&str),
    ) -> Result<(), FinalizeError> {
        let counts = std::mem::replace(&mut self.counts, Counts::Full(Vec::new()));
        match (self.format, counts) {
            (CountsFormat::Ascii, Counts::Full(counts)) => self.write_ascii(&counts, out)?,
            (CountsFormat::Binary, Counts::Full(counts)) => self.write_binary(&counts, out)?,
            (
                CountsFormat::OAscii { size_requested } | CountsFormat::OBinary { size_requested },
                Counts::Short(counts),
            ) => {
                log("optimizing the data structure");
                let binary = matches!(self.format, CountsFormat::OBinary { .. });
                let tables = optimize(
                    self.unit_size,
                    &self.units,
                    &counts,
                    size_requested,
                    binary.then_some(self.params[1]),
                )?;
                if binary {
                    self.write_obinary(&tables, out)?;
                } else {
                    self.write_oascii(&tables, out)?;
                }
            }
            _ => unreachable!("the counts are kept as the format keeps them"),
        }
        out.flush()?;
        Ok(())
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_ascii.cpp:114-134
    /// The text goes out in blocks of 64 KiB: the standard output is line-buffered, so a write
    /// per line would be a system call per unit.
    fn write_ascii(&self, counts: &[u32], out: &mut dyn Write) -> io::Result<()> {
        let mut text = Vec::with_capacity(1 << 16);
        text.extend_from_slice(&self.text_metadata());
        writeln!(text, "{}", self.unit_size)?;
        for (&unit, &count) in self.units.iter().zip(counts) {
            writeln!(text, "{unit:x} {count}")?;
            if text.len() >= 1 << 16 {
                out.write_all(&text)?;
                text.clear();
            }
        }
        text.push(b'\n');
        for comment in &self.comments {
            text.push(b'#');
            text.extend_from_slice(comment);
            text.push(b'\n');
        }
        text.push(b'\n');
        for (name, value) in PARAMS.iter().zip(self.params) {
            writeln!(text, ">{name} {value}")?;
        }
        text.push(b'\n');
        out.write_all(&text)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_bin.cpp:87-103
    fn write_binary(&self, counts: &[u32], out: &mut dyn Write) -> io::Result<()> {
        let mut bytes = 3u32.to_le_bytes().to_vec();
        bytes.extend_from_slice(&self.binary_metadata());
        bytes.extend_from_slice(&0u32.to_le_bytes());
        bytes.extend_from_slice(&u32::from(self.unit_size).to_le_bytes());
        out.write_all(&bytes)?;
        let mut words = Vec::with_capacity(8 * 4096);
        for (&unit, &count) in self.units.iter().zip(counts) {
            words.extend_from_slice(&unit.to_le_bytes());
            words.extend_from_slice(&count.to_le_bytes());
            if words.len() >= 8 * 4096 {
                out.write_all(&words)?;
                words.clear();
            }
        }
        for value in self.params {
            words.extend_from_slice(&value.to_le_bytes());
        }
        out.write_all(&words)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_opt_ascii.cpp:76-95
    fn write_oascii(&self, t: &OptTables, out: &mut dyn Write) -> io::Result<()> {
        let mut text = Vec::with_capacity(1 << 16);
        text.extend_from_slice(&self.text_metadata());
        text.extend_from_slice(b"AAAA\n");
        writeln!(text, "{}", self.unit_size)?;
        writeln!(text, "{} {} {} {}", t.m, t.k, t.roff, t.bc)?;
        for value in self.params {
            writeln!(text, "{value}")?;
        }
        let table_words = t.ht[..1usize << t.k].iter().copied();
        for value in table_words.chain(t.vt.iter().map(|&v| u32::from(v))) {
            writeln!(text, "{value}")?;
            if text.len() >= 1 << 16 {
                out.write_all(&text)?;
                text.clear();
            }
        }
        out.write_all(&text)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_opt_bin.cpp:72-108
    fn write_obinary(&self, t: &OptTables, out: &mut dyn Write) -> io::Result<()> {
        let mut head = 3u32.to_le_bytes().to_vec();
        head.extend_from_slice(&self.binary_metadata());
        let mut words = vec![2, u32::from(self.unit_size), t.m];
        words.extend([u32::from(t.k), u32::from(t.roff), u32::from(t.bc)]);
        words.extend(self.params);
        for word in words {
            head.extend_from_slice(&word.to_le_bytes());
        }
        match &t.cba {
            Some(cba) => {
                head.extend_from_slice(&1u32.to_le_bytes());
                out.write_all(&head)?;
                write_words(out, cba)?;
            }
            None => {
                head.extend_from_slice(&0u32.to_le_bytes());
                out.write_all(&head)?;
            }
        }
        write_words(out, &t.ht[..1usize << t.k])?;
        let mut bytes = Vec::with_capacity(2 * t.vt.len());
        for value in &t.vt {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        out.write_all(&bytes)
    }
}

/// `CWinMaskCountsConverter::operator()` between the writer's creation and its `finalize`: the
/// unit size of `istat`, the count of each unit that is not above its reverse complement and
/// whose stored count (`trueat`) is not 0, in ascending order, the four parameters of the
/// reader, its counts algorithm's version and its percentages (`count_map`, none when empty);
/// `log` gets NCBI's messages. Units of more than 16 letters (an ascii file's) are refused
/// ([`WinMaskError::UnitSize`]): NCBI converts them without an error, but its reverse complement
/// of a unit (`CSeqMaskerUtil::reverse_complement`, called for each unit) then shifts the 32-bit
/// unit by 32 bits or more (undefined behaviour in C++).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/win_mask_counts_converter.cpp:85-136
pub fn convert_counts(
    istat: &UnitCounts,
    count_map: Vec<f64>,
    writer: &mut CountsWriter,
    log: &mut dyn FnMut(&str),
) -> Result<(), WinMaskError> {
    let unit_size = istat.unit_size();
    if unit_size > 16 {
        return Err(WinMaskError::UnitSize { unit: unit_size });
    }
    writer.set_unit_size(unit_size);
    let num_units: u64 = 1 << (2 * u32::from(unit_size));
    log("converting counts...");
    for unit in 0..num_units {
        let unit = unit as u32;
        if unit <= reverse_complement(unit, unit_size) {
            let count = istat.true_count(unit)?;
            if count != 0 {
                writer.set_unit_count(unit, count);
            }
        }
    }
    log("converting parameters...");
    let params = [
        istat.min_count(),
        istat.textend(),
        istat.threshold(),
        istat.max_count(),
    ];
    for (i, value) in params.into_iter().enumerate() {
        writer.set_param(i, value);
    }
    log("final processing...");
    writer.set_stat_algo_version(istat.stat_algo_version());
    if !count_map.is_empty() {
        writer.set_count_map(count_map);
    }
    Ok(())
}

/// Writes `words` as little-endian `Uint4` words.
fn write_words(out: &mut dyn Write, words: &[u32]) -> io::Result<()> {
    let mut bytes = Vec::with_capacity(1 << 16);
    for chunk in words.chunks(1 << 14) {
        bytes.clear();
        for word in chunk {
            bytes.extend_from_slice(&word.to_le_bytes());
        }
        out.write_all(&bytes)?;
    }
    Ok(())
}

/// The hash table of an optimized counts file (`CSeqMaskerOstatOpt::params`): `m` units in the
/// secondary table `vt`, the key of `k` bits at `roff` bits from the right, `bc` bits of the
/// number of collisions in each word of the table `ht`, and the cache bit array of the
/// optimized binary format.
#[derive(Debug)]
struct OptTables {
    m: u32,
    k: u8,
    roff: u8,
    bc: u8,
    ht: Vec<u32>,
    vt: Vec<u16>,
    cba: Option<Vec<u32>>,
}

/// `CSeqMaskerOstatOpt::Exception::eMemory` at `line` of `doFinalize`.
fn memory_exception(line: u32) -> FinalizeError {
    FinalizeError::Ncbi(NcbiException::new(
        ThrowSite {
            file: "src/algo/winmask/seq_masker_ostat_opt.cpp",
            line,
            kind: "Exception::insufficient memory",
            function: Some("ncbi::CSeqMaskerOstatOpt::doFinalize"),
        },
        "Can not find parameters to satisfy memory requirements",
    ))
}

/// `CSeqMaskerOstatOpt::doFinalize` before its `write_out`: the cache bit array (for the
/// optimized binary format: `textend` is its `t_extend`), then the largest key size `k` whose
/// table fits `size_requested` megabytes, the best right offset for it (`find_best_roff`), and
/// the tables.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_opt.cpp:212-332
///
/// The key size is a `Uint1` and its bound, 7 less than the unit's bits, an `int`: with units of 1 to 3 letters and
/// no memory (`-smem` a multiple of 65536), the first loop takes `k` below 0, where it wraps,
/// and stops when the halved estimate reaches 0, at 253; the second loop then allocates 8 GiB
/// and shifts by 253 bits (undefined behaviour in C++): LOSAT refuses that case. With units of 0
/// letters (an ascii file's unit size line of 0, or of a multiple of 256, read as a `Uint1`) the
/// key size starts at 255 and the loop that sizes the estimate counts a `Uint1` up to 257, so it
/// never ends (undefined behaviour in C++ for a loop without side effects): refused too.
fn optimize(
    unit_size: u8,
    units: &[u32],
    counts: &[u16],
    size_requested: u16,
    cba_textend: Option<u32>,
) -> Result<OptTables, FinalizeError> {
    if unit_size == 0 {
        return Err(FinalizeError::Refused(
            "the optimized formats (oascii, obinary) of units of 0 letters, where NCBI BLAST+'s search for the size of its hash table never ends (undefined behaviour in C++)".to_string(),
        ));
    }
    let ubs = 2 * i32::from(unit_size);
    let limit = u64::from(size_requested) * MB;
    if ubs < 7 && limit == 0 {
        return Err(FinalizeError::Refused(format!(
            "-smem of a multiple of 65536 (whose low 16 bits, 0, NCBI BLAST+ keeps) for units of {unit_size} letters, where NCBI BLAST+ sizes its hash table for a key of 253 bits, allocates 8 GiB and shifts past the width of its integers (undefined behaviour in C++)"
        )));
    }
    // NCBI makes the array for the optimized ascii format too, which does not write it: LOSAT
    // only checks that it could be allocated.
    let cba = match cba_textend {
        Some(textend) => Some(cache_bit_array(unit_size, units, counts, textend)?),
        None => {
            cache_bit_array_space(unit_size)?;
            None
        }
    };
    let mut k = ubs - 1;
    let mut emem: u64 = 1 << (k + 2);
    while k >= ubs - 7 && emem > limit {
        emem /= 2;
        k -= 1;
    }
    if k < ubs - 7 {
        return Err(memory_exception(242));
    }
    let size = 1usize << k;
    let mut ht = vec![0u32; size];
    let (mut roff, mut m, mut bc) = (0u8, 0u32, 0u8);
    while k >= ubs - 7 {
        let best = find_best_roff(ubs as u8, k as u8, units, &mut ht);
        roff = best.0;
        m = best.2;
        bc = 0;
        while (1u64 << bc) <= u64::from(best.1) {
            bc += 1;
        }
        if bc <= 7 {
            let mut b0 = 0;
            while (1u64 << b0) <= u64::from(m) {
                b0 += 1;
            }
            if b0 + bc <= 32 && u64::from(2u32.wrapping_mul(m)) + (1u64 << (k + 2)) <= limit {
                break;
            }
        }
        k -= 1;
    }
    if k < ubs - 7 {
        return Err(memory_exception(275));
    }
    let k = k as u8;
    let vt = fill_tables(&mut ht, units, counts, k, roff, bc, m);
    Ok(OptTables {
        m,
        k,
        roff,
        bc,
        ht,
        vt,
        cba,
    })
}

/// `createCacheBitArray`: one bit per unit, set for each unit whose count is at least `textend`
/// and for its reverse complement; an array of no words (units of 1 or 2 letters) has none set.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_opt.cpp:120-160
fn cache_bit_array(
    unit_size: u8,
    units: &[u32],
    counts: &[u16],
    textend: u32,
) -> Result<Vec<u32>, FinalizeError> {
    let (mut cba, size) = cache_bit_array_space(unit_size)?;
    if size == 0 {
        return Ok(cba);
    }
    cba.resize(size, 0u32);
    for (&unit, &count) in units.iter().zip(counts) {
        if u32::from(count) >= textend {
            let rcunit = reverse_complement(unit, unit_size);
            cba[(unit / 32) as usize] |= 1 << (unit % 32);
            cba[(rcunit / 32) as usize] |= 1 << (rcunit % 32);
        }
    }
    Ok(cba)
}

/// Room for the cache bit array of units of `unit_size` letters, and its number of words
/// (`4^unit_size / 32`, none for units of 1 or 2 letters). NCBI posts a warning and writes no
/// array when the array cannot be allocated, which depends on the memory of the host: LOSAT
/// refuses that case.
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_opt.cpp:122-139
fn cache_bit_array_space(unit_size: u8) -> Result<(Vec<u32>, usize), FinalizeError> {
    let total: u64 = 1 << (2 * u32::from(unit_size));
    let size = (total / 32) as usize;
    let mut cba = Vec::new();
    if cba.try_reserve_exact(size).is_err() {
        return Err(FinalizeError::Refused(format!(
            "a cache bit array of {size} words that cannot be allocated here (NCBI BLAST+ then writes the file without one, which depends on the memory of the host)"
        )));
    }
    Ok((cba, size))
}

/// `findBestRoff` of the key size `k` for units of `ubs` bits: for each right offset from 0 to `k` less than `ubs`,
/// the number of units of each key of `ht`; the offset with the lowest average number of units
/// in a key that has more than one (the first such), its largest number of units in a key (as
/// a `Uint1`) and its number of units in keys with more than one (`M`).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_opt.cpp:163-209
fn find_best_roff(ubs: u8, k: u8, units: &[u32], ht: &mut [u32]) -> (u8, u8, u32) {
    let size = 1usize << k;
    let table = &mut ht[..size];
    let mut best: Option<(f64, u8, u8, u32)> = None;
    for roff in 0..=(ubs - k) {
        table.fill(0);
        for &unit in units {
            let slot = &mut table[hash_code(unit, k, roff).0 as usize];
            *slot = slot.wrapping_add(1);
        }
        let mxcoll = table.iter().copied().max().unwrap_or(0) as u8;
        let (mut t, mut tc) = (0u32, 0u32);
        for &n in table.iter() {
            if n > 1 {
                tc = tc.wrapping_add(1);
                t = t.wrapping_add(n);
            }
        }
        let avcoll = if tc > 0 {
            f64::from(t) / f64::from(tc)
        } else {
            0.0
        };
        if best.is_none_or(|(least, ..)| avcoll < least) {
            best = Some((avcoll, roff, mxcoll, t));
        }
    }
    let (_, roff, mxcoll, m) = best.expect("at least one offset");
    (roff, mxcoll, m)
}

/// The tables of the chosen parameters: `ht` holds for each key its number of units in the low
/// `bc` bits, and above them the count and check bits of a key with one unit, or the index in
/// `vt` of the last of its units for a key with more (filled from the last slot down); each
/// `vt` entry is the check bits above 9 bits of count (`Uint2`: larger counts run into the check
/// bits, as in NCBI's files).
///
/// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_ostat_opt.cpp:286-328
fn fill_tables(
    ht: &mut [u32],
    units: &[u32],
    counts: &[u16],
    k: u8,
    roff: u8,
    bc: u8,
    m: u32,
) -> Vec<u16> {
    ht.fill(0);
    for &unit in units {
        let slot = &mut ht[hash_code(unit, k, roff).0 as usize];
        *slot = slot.wrapping_add(1);
    }
    let mut vt = vec![0u16; m as usize];
    let mut vend: u32 = 0;
    let coll_mask: u32 = ((1u64 << bc) - 1) as u32;
    for (&unit, &count) in units.iter().zip(counts) {
        let (key, check) = hash_code(unit, k, roff);
        let word = &mut ht[key as usize];
        let ccount = (*word & coll_mask) as u8;
        if ccount == 0 {
            continue;
        }
        if ccount == 1 {
            *word = word.wrapping_add(u32::from(count) << bc);
            *word = word.wrapping_add(u32::from(check) << 24);
        } else {
            if *word & !coll_mask == 0 {
                vend = vend.wrapping_add(u32::from(ccount));
                *word = word.wrapping_add(vend.wrapping_sub(1) << bc);
            } else {
                *word = word.wrapping_sub(1 << bc);
            }
            vt[(*word >> bc) as usize] =
                (u32::from(check) << 9).wrapping_add(u32::from(count)) as u16;
        }
    }
    vt
}

#[cfg(test)]
#[path = "ostat_tests.rs"]
mod tests;
