//! The CDD annotation data of rpsbproc (`COfflDomClusterData`): the six text files of NCBI's
//! `pub/mmdb/cdd/` directory, read as NCBI's offline loader reads them, and the lookups by
//! PSSM-ID, accession and cluster that the classification makes.
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.hpp:12-26
//!
//! NCBI keeps a `TDomain` (or `TCluster`) per row with every column as a `std::string`; LOSAT
//! keeps the columns that the outputs print or decide by (the PSSM-ID, accession, short name,
//! length, curated and structural flags, cluster, bit score threshold, hierarchy parent and
//! root, and the sites of `cddannot*.dat`: title, motif, positions, type, size) with the texts
//! in one buffer. The definition lines and the site columns no output reads (index, evidence
//! flags) are read as NCBI reads them, for their errors, but not kept.
//! The lookups are NCBI's `CPtrMap`: the rows sorted with NCBI's comparison and searched by
//! NCBI's binary search. The sort is Rust's standard sort with the row's place in the file as
//! the last key (`sort_ties_in_input_order`), so rows of one key stay in the order of the file
//! (NCBI's `std::sort` may leave them in another order; PD-LOSAT-NCBI-DEFECTS 1.67, the order of
//! tied hits).

use std::io::{BufRead, BufReader};

use crate::ncbi_exception::NcbiException;
use crate::ncbi_str::{printable_string, string_to_double, string_to_int, string_to_uint};
use crate::report::serial::asn_input::VECTOR_RESERVE_LENGTH_ERROR;

/// `INVALIDPSSMID`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:26-30
pub(crate) const INVALID_PSSM_ID: u32 = 100_000_000;
const DOMSRCMUL: u32 = 1_000_000;
const INVALID_CLUSTER_ID: i32 = -1;
const SINGLE_MEMBER_CLUSTER: i32 = i32::MAX;

/// A text in `CddData::text`.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Span {
    start: u32,
    len: u32,
}

/// `TDomain`: a model of `cddid.tbl` that is not a superfamily.
#[derive(Clone, Debug)]
pub(crate) struct Domain {
    pub(crate) pssm_id: u32,
    pub(crate) accession: Span,
    pub(crate) short_name: Span,
    /// `m_uiLength`.
    pub(crate) length: u32,
    /// `m_iClusterId`: the superfamily, `INVALIDCLUSTERID` (-1) when none.
    pub(crate) cluster_id: i32,
    /// `m_dMinBitScore` (from `bitscore_specific.txt`).
    pub(crate) min_bit_score: f64,
    /// `m_bCurated`: a `cd` accession.
    pub(crate) curated: bool,
    /// `m_bIsStructDom`: an `sd` accession.
    pub(crate) is_struct_dom: bool,
    /// `m_uiHierarchyParent` (from `cdtrack.txt`, 0 when not set).
    pub(crate) hierarchy_parent: u32,
    /// `m_uiHierarchyRoot` (from `cdtrack.txt`, 0 when not set).
    pub(crate) hierarchy_root: u32,
    /// `m_lstSpecificFeatures`: the range in `CddData::specific_sites`.
    specific: SiteRange,
    /// `m_lstGenericFeatures`: the range in `CddData::generic_sites`.
    generic: SiteRange,
}

/// A run of sites of one domain in a site table (`start`, `len`).
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct SiteRange {
    start: u32,
    len: u32,
}

/// `TDomSite`, with the members the site and motif output reads. Its positions
/// (`CSegSet`, after `SetData`: sorted, merged ranges) are `CddData::site_segs[segs]`.
#[derive(Clone, Debug)]
pub(crate) struct Site {
    /// `m_strTitle`.
    pub(crate) title: Span,
    /// `m_strMotif` (empty: no motif).
    pub(crate) motif: Span,
    /// `m_iType`.
    pub(crate) site_type: i32,
    /// `m_iCompleteSize`: the number of positions of the line for a specific site (duplicates
    /// included), 0 for a generic site (`LoadData_real` sets it for the specific ones only).
    pub(crate) complete_size: i32,
    segs: SiteRange,
    /// The domain of the line (until the sites are grouped by domain).
    domain: u32,
}

/// `TDomSite::eType_StructMotif`.
pub(crate) const TYPE_STRUCT_MOTIF: i32 = 7;

/// `TCluster`: a superfamily row of `cddid.tbl` (a `cl` accession).
#[derive(Clone, Debug)]
pub(crate) struct Cluster {
    pub(crate) pssm_id: u32,
    pub(crate) accession: Span,
    pub(crate) short_name: Span,
    pub(crate) cluster_id: i32,
}

/// The data of one run.
#[derive(Debug, Default)]
pub(crate) struct CddData {
    text: Vec<u8>,
    pub(crate) domains: Vec<Domain>,
    pub(crate) clusters: Vec<Cluster>,
    /// `CDomClusterIndex::m_pssmid2cd`, `m_acxn2cd`, `m_clid2fam` (the doms path never asks
    /// `m_pssmid2fam`, `FindClusterByPssmId`, which is not kept).
    domain_by_pssm_id: Vec<u32>,
    domain_by_accession: Vec<u32>,
    cluster_by_id: Vec<u32>,
    /// The specific and generic sites, grouped by domain in file order (`Domain::specific`,
    /// `Domain::generic`), and the ranges of their positions.
    specific_sites: Vec<Site>,
    generic_sites: Vec<Site>,
    site_segs: Vec<[i32; 2]>,
}

/// Why the data does not load.
#[derive(Debug)]
pub(crate) enum DataError {
    /// `CSimpleException` "Data file `<path>` does not exist" (rpsbproc: "Annotation data loading
    /// unsuccessful: ...", exit 2).
    Missing(String),
    /// A number that NCBI's conversion refuses (rpsbproc: two `Error:` lines, exit 2).
    Conversion(NcbiException),
    /// An exception of the C++ library that leaves `Init` (rpsbproc: one `Error: (106.6)` line
    /// with its `what()`, exit 2).
    Std(String),
    /// A line that LOSAT does not read (NCBI reads outside the line).
    Unsupported(String),
}

/// The paths of the six files (`LoadData_real`'s arguments).
pub(crate) struct DataFiles {
    pub(crate) cddids: String,
    pub(crate) spfeats: String,
    pub(crate) genfeats: String,
    pub(crate) bscores: String,
    pub(crate) cdtrack: String,
    pub(crate) famlinks: String,
}

/// `std::string::npos`.
const NPOS: usize = usize::MAX;

/// One line of a data file with the `std::string` operations NCBI's loader makes on it: `find`
/// from a position (`NPOS` when not found; positions past the end find nothing), `substr`, and
/// `operator[]` (the terminating 0 at the line's length).
struct Line<'a> {
    bytes: &'a [u8],
    file: &'a str,
    /// Whether `operator[]` read the byte before the line (`buf[npos]`, `at`).
    read_before: std::cell::Cell<bool>,
}

impl<'a> Line<'a> {
    fn new(bytes: &'a [u8], file: &'a str) -> Self {
        Line {
            bytes,
            file,
            read_before: std::cell::Cell::new(false),
        }
    }

    fn find(&self, byte: u8, from: usize) -> usize {
        if from >= self.bytes.len() {
            return NPOS;
        }
        self.bytes[from..]
            .iter()
            .position(|&b| b == byte)
            .map_or(NPOS, |at| from + at)
    }

    /// `substr(pos, n)`: `n` clipped to the end; `pos` past the end throws `out_of_range`.
    // The C++ standard's `basic_string::substr` throws `out_of_range` when `pos > size()`; the text
    // is that of NCBI's build, observed: NCBI's rpsbproc printed `Error: (106.6) Application's initialization failed:
    // basic_string::substr: __pos (which is 18446744073709551615) > this->size() (which is 10)`
    // for a `cdtrack.txt` line of an accession alone, PS16a review oracle A4-1 case `ct3`.)
    //
    // The text is that of NCBI's build, whose `size_t` has 64 bits (NCBI BLAST+ 2.17.0 for Linux
    // x64, the only reference): `pos` prints as `%zu` of its 64-bit value on every target, where
    // `npos`, or the difference of two positions that wrapped below 0, is near 2^64 (not
    // 4294967295 as a 32-bit `usize` would print). A line is shorter than 2 GiB, so a position
    // of 2^31 or more is such a wrapped one: it takes its 64-bit value by sign extension.
    fn substr(&self, pos: usize, n: usize) -> Result<&'a [u8], DataError> {
        if pos > self.bytes.len() {
            let pos = pos as isize as i64 as u64;
            return Err(DataError::Std(format!(
                "basic_string::substr: __pos (which is {pos}) > this->size() (which is {})",
                self.bytes.len() as u64
            )));
        }
        let end = pos.saturating_add(n).min(self.bytes.len());
        Ok(&self.bytes[pos..end])
    }

    /// `buf[pos]`: the byte, 0 at the end. Other reads outside the line are undefined in C++:
    /// `buf[npos]`, the byte before the string's data, reads as 0 (the top byte of the string's
    /// length in its own buffer, or of malloc's size word before a longer string's data) and is
    /// noted (`read_before`); a read past the terminating 0 is refused.
    fn at(&self, pos: usize) -> Result<u8, DataError> {
        match pos.cmp(&self.bytes.len()) {
            std::cmp::Ordering::Less => Ok(self.bytes[pos]),
            std::cmp::Ordering::Equal => Ok(0),
            std::cmp::Ordering::Greater if pos == NPOS => {
                self.read_before.set(true);
                Ok(0)
            }
            std::cmp::Ordering::Greater => Err(self.unsupported("a read past its end")),
        }
    }

    fn unsupported(&self, what: &str) -> DataError {
        DataError::Unsupported(format!(
            "a line of {} where NCBI's loader makes {what} is not supported by LOSAT's rpsbproc: '{}'",
            self.file,
            String::from_utf8_lossy(&printable_string(self.bytes))
        ))
    }
}

/// `x_GetTrimPos` and the other trims of rpsbproc's library: a byte is white space when it is at
/// most `' '` as a signed `char` (so the bytes from 0x80 too).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basealgo.cpp:18-28
fn is_trim_space(byte: u8) -> bool {
    (byte as i8) <= b' ' as i8
}

/// `TrimString`.
fn trim(text: &[u8]) -> &[u8] {
    let start = text
        .iter()
        .position(|&b| !is_trim_space(b))
        .unwrap_or(text.len());
    let end = text
        .iter()
        .rposition(|&b| !is_trim_space(b))
        .map_or(start, |at| at + 1);
    &text[start..end.max(start)]
}

/// `LTrimString`.
fn ltrim(text: &[u8]) -> &[u8] {
    let start = text
        .iter()
        .position(|&b| !is_trim_space(b))
        .unwrap_or(text.len());
    &text[start..]
}

// NCBI reference (598d8ae6): c++/include/corelib/ncbistr.hpp:3582-3593
// The loader's `NStr::StringToNumeric<T>(text)` (flags 0, base 10) is `StringToUInt` for the
// unsigned ids and lengths, `StringToInt` for the signed ones and `StringToDouble` for the bit
// scores (ncbistr.hpp:3632-3637): `ncbi_str::string_to_uint`, `string_to_int` and
// `string_to_double`. Their `CStringException` (`eConvert`) is not caught by
// `catch (CSimpleException e)` (rpsbproc_main.cpp:654), so it leaves `Init` to the
// application's handler (`DataError::Conversion`).

/// `TDomSrcCount::ESrcIdx` and the signatures of the source databases.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:117-140
const DOMSRCSIGS: [&[u8]; 12] = [
    b"CD", b"PFAM", b"TIGR", b"COG", b"KOG", b"SMART", b"PRK", b"CHL", b"MTH", b"PHA", b"PLN",
    b"PTZ",
];
/// `TOTALSIGS`.
pub(crate) const TOTALSIGS: usize = 12;

pub(crate) fn dom_acc_type(accession: &[u8]) -> usize {
    DOMSRCSIGS
        .iter()
        .position(|sig| {
            sig.len() < accession.len()
                && sig
                    .iter()
                    .zip(accession)
                    .all(|(&s, &a)| a.to_ascii_uppercase() == s)
        })
        .unwrap_or(TOTALSIGS)
}

/// `TDomSrcCount::ComputePseudoPssmId`: the PSSM-ID of a model known only by accession.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:168-193
// (`PssmId_t` is `unsigned int`: the arithmetic wraps.)
pub(crate) fn compute_pseudo_pssm_id(accession: &[u8]) -> u32 {
    if accession.is_empty() {
        return 0;
    }
    let source = dom_acc_type(accession) as u32;
    let mut pseudo = INVALID_PSSM_ID.wrapping_add(source.wrapping_mul(DOMSRCMUL));
    let start = accession
        .iter()
        .position(u8::is_ascii_digit)
        .unwrap_or(accession.len());
    let mut order: u32 = 0;
    for &digit in accession[start..].iter().take_while(|b| b.is_ascii_digit()) {
        order = order.wrapping_mul(10).wrapping_add(u32::from(digit - b'0'));
    }
    pseudo = pseudo.wrapping_add(order);
    pseudo
}

/// Reads the `\n`-separated lines of `path` (`getline` while the stream is good: a last line
/// without its `\n` is read, and a file that does not open has no lines).
fn for_each_line(
    path: &str,
    mut each: impl FnMut(&[u8]) -> Result<(), DataError>,
) -> Result<(), DataError> {
    let Ok(file) = std::fs::File::open(path) else {
        return Ok(());
    };
    let mut reader = BufReader::with_capacity(1 << 16, file);
    let mut buf = Vec::new();
    loop {
        buf.clear();
        match reader.read_until(b'\n', &mut buf) {
            Ok(0) | Err(_) => return Ok(()),
            Ok(_) => {
                if buf.last() == Some(&b'\n') {
                    buf.pop();
                }
                each(&buf)?;
            }
        }
    }
}

/// `CFile(path).Exists()`: a regular file (links followed).
fn file_exists(path: &str) -> bool {
    std::fs::metadata(path).is_ok_and(|meta| meta.is_file())
}

impl CddData {
    /// The text `span`.
    pub(crate) fn text(&self, span: Span) -> &[u8] {
        &self.text[span.start as usize..(span.start + span.len) as usize]
    }

    fn push_text(&mut self, text: &[u8]) -> Span {
        let start = self.text.len() as u32;
        self.text.extend_from_slice(text);
        Span {
            start,
            len: text.len() as u32,
        }
    }

    /// `CPtrMap::Find` of `m_pssmid2cd`: the domain of `pssm_id`.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/ptrmap.hpp:109-127
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:796-833
    // (The `int` difference is compared without its overflow: `int_order`.)
    pub(crate) fn find_domain(&self, pssm_id: u32) -> Option<usize> {
        let found = binary_search(&self.domain_by_pssm_id, |index| {
            int_order(pssm_id as i32, self.domains[index as usize].pssm_id as i32)
        })?;
        Some(found as usize)
    }

    /// `FindCdInfo(acxn)`: the domain of an accession.
    pub(crate) fn find_domain_by_accession(&self, accession: &[u8]) -> Option<usize> {
        let found = binary_search(&self.domain_by_accession, |index| {
            compare_bytes(accession, self.text(self.domains[index as usize].accession))
        })?;
        Some(found as usize)
    }

    /// `FindClusterInfo(clusterid)`.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:809-820
    pub(crate) fn find_cluster(&self, cluster_id: i32) -> Option<usize> {
        let found = binary_search(&self.cluster_by_id, |index| {
            int_order(cluster_id, self.clusters[index as usize].cluster_id)
        })?;
        Some(found as usize)
    }

    /// `COfflDomClusterData::LoadData_real`.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.cpp:66-83
    // The files are read in this order: `cddid.tbl`, `cdtrack.txt`, `family_superfamily_links`,
    // `cddannot.dat`, `cddannot_generic.dat`, `bitscore_specific.txt`; each one is checked when
    // its turn comes.
    pub(crate) fn load(files: &DataFiles) -> Result<CddData, DataError> {
        let mut data = CddData::default();
        data.load_cddids(&files.cddids)?;
        data.load_cdtrack(&files.cdtrack)?;
        data.load_famlinks(&files.famlinks)?;
        data.load_sites(&files.spfeats, true)?;
        data.load_sites(&files.genfeats, false)?;
        data.load_bscores(&files.bscores)?;
        data.text.shrink_to_fit();
        Ok(data)
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.cpp:100-213
    // (`pos1 + 1` of `npos` is 0: a line with fewer tabs than NCBI expects reads its fields from
    // the start again. The definition lines are not kept.)
    fn load_cddids(&mut self, path: &str) -> Result<(), DataError> {
        if !file_exists(path) {
            return Err(DataError::Missing(path.to_string()));
        }
        const FILE: &str = "cddid.tbl";
        let mut domains = Vec::new();
        let mut clusters = Vec::new();
        for_each_line(path, |bytes| {
            if bytes.is_empty() || bytes[0] == b'#' {
                return Ok(());
            }
            let line = Line::new(bytes, FILE);
            let mut pos0 = 0usize;
            let mut pos1 = line.find(b'\t', pos0);
            let pssm_id_text = ltrim(line.substr(pos0, pos1.wrapping_sub(pos0))?);
            let mut pssm_id = 0u32;
            if !pssm_id_text.is_empty() {
                pssm_id = string_to_uint(pssm_id_text).map_err(DataError::Conversion)?;
            }
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            let accession = trim(line.substr(pos0, pos1.wrapping_sub(pos0))?);
            if pssm_id == 0 {
                pssm_id = compute_pseudo_pssm_id(accession);
            }
            // `tmp[i]` of a `std::string` is its terminating 0 at its length.
            let acc = |i: usize| accession.get(i).copied().unwrap_or(0);
            if acc(0) == b'c' && acc(1) == b'l' && acc(2).is_ascii_digit() {
                let cluster_id =
                    string_to_int(&accession[2..], 0).map_err(DataError::Conversion)?;
                pos0 = pos1.wrapping_add(1);
                pos1 = line.find(b'\t', pos0);
                let short_name = trim(line.substr(pos0, pos1.wrapping_sub(pos0))?);
                pos0 = pos1.wrapping_add(1);
                pos1 = line.find(b'\t', pos0);
                line.substr(pos0, pos1.wrapping_sub(pos0))?;
                let accession = self.push_text(accession);
                let short_name = self.push_text(short_name);
                clusters.push(Cluster {
                    pssm_id,
                    accession,
                    short_name,
                    cluster_id,
                });
                return Ok(());
            }
            // `m_strAccession[2]` of an accession shorter than 2 bytes reads past its end, which is
            // undefined, but no flag can be set then: `[1]` of a 1-byte accession and `[0]` of an
            // empty one are the terminating 0, never `d`, `c` or `s` (`acc` gives 0 there).
            let (mut curated, mut is_struct_dom) = (false, false);
            if acc(2).is_ascii_digit() && acc(1) == b'd' {
                if acc(0) == b'c' {
                    curated = true;
                } else if acc(0) == b's' {
                    is_struct_dom = true;
                }
            }
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            let short_name = trim(line.substr(pos0, pos1.wrapping_sub(pos0))?);
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            line.substr(pos0, pos1.wrapping_sub(pos0))?;
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            while pos1 != NPOS {
                pos0 = pos1.wrapping_add(1);
                pos1 = line.find(b'\t', pos0);
            }
            let length_text = trim(line.substr(pos0, NPOS)?);
            let length = string_to_uint(length_text).map_err(DataError::Conversion)?;
            let accession = self.push_text(accession);
            let short_name = self.push_text(short_name);
            domains.push(Domain {
                pssm_id,
                accession,
                short_name,
                length,
                cluster_id: INVALID_CLUSTER_ID,
                min_bit_score: 0.0,
                curated,
                is_struct_dom,
                hierarchy_parent: 0,
                hierarchy_root: 0,
                specific: SiteRange::default(),
                generic: SiteRange::default(),
            });
            Ok(())
        })?;
        self.domains = domains;
        self.clusters = clusters;
        // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/ptrmap.hpp:167-179
        // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/ptrmap.hpp:63-70
        let domains = &self.domains;
        let clusters = &self.clusters;
        let text = &self.text;
        let span = |s: Span| &text[s.start as usize..(s.start + s.len) as usize];
        // (LOSAT sorts with Rust's standard sort (`sort_ties_in_input_order`) and the comparisons
        // that the lookups search with, `int_order` and `compare_bytes`, which order every pair of
        // keys, then the row's index: rows of one key keep the order of the file, where NCBI's
        // `sort` may leave them in another. No CDD v3.22 file has two rows of one key.
        // PD-LOSAT-NCBI-DEFECTS 1.67, the order of tied hits.)
        let mut by_pssm: Vec<u32> = (0..domains.len() as u32).collect();
        sort_ties_in_input_order(
            &mut by_pssm,
            |&row| row,
            |&a, &b| {
                int_order(
                    domains[a as usize].pssm_id as i32,
                    domains[b as usize].pssm_id as i32,
                )
                .cmp(&0)
            },
        );
        let mut by_acc: Vec<u32> = (0..domains.len() as u32).collect();
        sort_ties_in_input_order(
            &mut by_acc,
            |&row| row,
            |&a, &b| {
                compare_bytes(
                    span(domains[a as usize].accession),
                    span(domains[b as usize].accession),
                )
                .cmp(&0)
            },
        );
        let mut cl_by_id: Vec<u32> = (0..clusters.len() as u32).collect();
        sort_ties_in_input_order(
            &mut cl_by_id,
            |&row| row,
            |&a, &b| {
                int_order(
                    clusters[a as usize].cluster_id,
                    clusters[b as usize].cluster_id,
                )
                .cmp(&0)
            },
        );
        self.domain_by_pssm_id = by_pssm;
        self.domain_by_accession = by_acc;
        self.cluster_by_id = cl_by_id;
        Ok(())
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.cpp:217-308
    // The parent and the root set `m_uiHierarchyParent` and `m_uiHierarchyRoot`, which the site
    // output reads (`pRootCdInfo`; the `-t doms` output does not). The root column that is
    // neither the parent nor the domain itself sets the *parent* (NCBI's quirk, survey 8.4), and
    // the root stays 0. A `find` that fails
    // gives `npos`; the loops `while(buf[pos1] <= ' ')` (a signed `char`) then read `buf[npos]`,
    // the byte before the string's data, which is undefined but 0 in NCBI's build (`Line::at`),
    // and go on from position 0, and a `substr` from `npos` throws `std::out_of_range` (oracle
    // A4-1 case `ct3`, a line of the accession alone: rpsbproc exits 2 with libstdc++'s message).
    // A line that ends without an exception after reading `buf[npos]` (oracle A4-1 case `ct1`,
    // `smart00421 HTH_LUXR 197715 N/A`, a line whose root is missing) is refused: its result
    // rests on the undefined read. NCBI's published `cdtrack.txt` (the v3.23 preview, oracle A4-2:
    // LOSAT's bytes and NCBI's agree on it) has no such line; a read past the terminating 0 is
    // refused too.
    fn load_cdtrack(&mut self, path: &str) -> Result<(), DataError> {
        if !file_exists(path) {
            return Err(DataError::Missing(path.to_string()));
        }
        const FILE: &str = "cdtrack.txt";
        for_each_line(path, |raw| {
            let bytes = ltrim(raw);
            if bytes.is_empty() || bytes[0] == b'#' {
                return Ok(());
            }
            let line = Line::new(bytes, FILE);
            let skip_le_space = |mut pos: usize| -> Result<usize, DataError> {
                while is_trim_space(line.at(pos)?) {
                    pos = pos.wrapping_add(1);
                }
                Ok(pos)
            };
            let mut pos1 = line.find(b' ', 0);
            let current = trim(line.substr(0, pos1)?);
            let Some(dom) = self.find_domain_by_accession(current) else {
                return Ok(());
            };
            pos1 = skip_le_space(pos1)?;
            let mut pos0 = pos1;
            pos1 = line.find(b' ', pos0);
            while line.at(pos1)? == b' ' {
                pos1 = pos1.wrapping_add(1);
            }
            if line.at(pos1)?.is_ascii_digit() {
                pos0 = pos1;
                pos1 = line.find(b' ', pos0);
                pos1 = skip_le_space(pos1)?;
            }
            pos0 = pos1;
            pos1 = line.find(b' ', pos0);
            let parent_acc = trim(line.substr(pos0, pos1.wrapping_sub(pos0))?);
            let own_acc = self.text(self.domains[dom].accession);
            let own_pssm_id = self.domains[dom].pssm_id;
            if parent_acc != b"N/A" {
                if parent_acc == own_acc {
                    self.domains[dom].hierarchy_parent = own_pssm_id;
                } else if let Some(parent) = self.find_domain_by_accession(parent_acc) {
                    self.domains[dom].hierarchy_parent = self.domains[parent].pssm_id;
                }
            }
            pos1 = skip_le_space(pos1)?;
            pos0 = pos1;
            pos1 = line.find(b' ', pos0);
            let root_acc = trim(line.substr(pos0, pos1.wrapping_sub(pos0))?);
            let own_acc = self.text(self.domains[dom].accession);
            if root_acc == parent_acc {
                self.domains[dom].hierarchy_root = self.domains[dom].hierarchy_parent;
            } else if root_acc == own_acc {
                self.domains[dom].hierarchy_root = own_pssm_id;
            } else if let Some(root) = self.find_domain_by_accession(root_acc) {
                self.domains[dom].hierarchy_parent = self.domains[root].pssm_id;
            }
            if line.read_before.get() {
                return Err(line.unsupported("a read of the byte before the line (buf[npos])"));
            }
            Ok(())
        })
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.cpp:310-386
    // (The domain is looked up by accession: the PSSM-ID column is not read.)
    fn load_famlinks(&mut self, path: &str) -> Result<(), DataError> {
        if !file_exists(path) {
            return Err(DataError::Missing(path.to_string()));
        }
        const FILE: &str = "family_superfamily_links";
        for_each_line(path, |raw| {
            let bytes = ltrim(raw);
            if bytes.is_empty() || bytes[0] == b'#' {
                return Ok(());
            }
            let line = Line::new(bytes, FILE);
            let mut pos0 = 0usize;
            let mut pos1 = line.find(b'\t', pos0);
            let Some(dom) =
                self.find_domain_by_accession(line.substr(pos0, pos1.wrapping_sub(pos0))?)
            else {
                return Ok(());
            };
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            let sf = line.substr(pos0, pos1.wrapping_sub(pos0))?;
            let sf_at = |i: usize| sf.get(i).copied().unwrap_or(0);
            if sf_at(0) == b'c' && sf_at(1) == b'l' {
                let number = line.substr(pos0 + 2, pos1.wrapping_sub(pos0).wrapping_sub(2))?;
                let cluster_id = string_to_int(number, 0).map_err(DataError::Conversion)?;
                if self.find_cluster(cluster_id).is_some() {
                    self.domains[dom].cluster_id = cluster_id;
                }
            } else {
                self.domains[dom].cluster_id = SINGLE_MEMBER_CLUSTER;
            }
            Ok(())
        })
    }

    /// The sites of `cddannot.dat` (`specific`) or `cddannot_generic.dat`.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.cpp:388-504
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.cpp:518-609
    // The generic file is read without the trims of the specific one (the line, the accession,
    // the index and the title keep their white space). A coordinates field at the line's end
    // (`pos1` is `npos`) makes NCBI reserve `(npos - pos0) * 3` positions, more than a
    // `vector<int>` can hold whatever `pos0` is, which throws `std::length_error` (the C++
    // standard's `vector::reserve` throws it for a count above `max_size()`; the text is observed:
    // NCBI's rpsbproc printed `Error: (106.6) Application's initialization failed:
    // vector::reserve` and exited 2, PS16a review oracle A4-1 cases `annot10` and `gen10`.)
    fn load_sites(&mut self, path: &str, specific: bool) -> Result<(), DataError> {
        if !file_exists(path) {
            return Err(DataError::Missing(path.to_string()));
        }
        let file = if specific {
            "cddannot.dat"
        } else {
            "cddannot_generic.dat"
        };
        let mut sites: Vec<Site> = Vec::new();
        let mut site_segs = std::mem::take(&mut self.site_segs);
        let mut residues: Vec<i32> = Vec::new();
        for_each_line(path, |raw| {
            let bytes = if specific { ltrim(raw) } else { raw };
            if bytes.is_empty() || bytes[0] == b'#' {
                return Ok(());
            }
            let line = Line::new(bytes, file);
            let mut pos0;
            let mut pos1 = line.find(b'\t', 0);
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            let accession = line.substr(pos0, pos1.wrapping_sub(pos0))?;
            let accession = if specific { trim(accession) } else { accession };
            let Some(domain) = self.find_domain_by_accession(accession) else {
                return Ok(());
            };
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            let index_text = line.substr(pos0, pos1.wrapping_sub(pos0))?;
            let index_text = if specific {
                trim(index_text)
            } else {
                index_text
            };
            string_to_int(index_text, 0).map_err(DataError::Conversion)?;
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            let title = line.substr(pos0, pos1.wrapping_sub(pos0))?;
            let title = if specific { trim(title) } else { title };
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            let mut motif: &[u8] = b"";
            if line.at(pos0)? != b'0' {
                motif = trim(line.substr(pos0, pos1.wrapping_sub(pos0))?);
            }
            for _evidence_flag in 0..3 {
                pos0 = pos1.wrapping_add(1);
                pos1 = line.find(b'\t', pos0);
                line.at(pos0)?;
            }
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            if pos1 == NPOS {
                return Err(DataError::Std(VECTOR_RESERVE_LENGTH_ERROR.to_string()));
            }
            residues.clear();
            let mut pos_c = line.find(b',', pos0);
            while pos_c < pos1 {
                residues.push(
                    string_to_int(line.substr(pos0, pos_c - pos0)?, 0)
                        .map_err(DataError::Conversion)?,
                );
                pos0 = pos_c + 1;
                pos_c = line.find(b',', pos0);
            }
            residues.push(
                string_to_int(line.substr(pos0, pos1.wrapping_sub(pos0))?, 0)
                    .map_err(DataError::Conversion)?,
            );
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            let type_text = if pos1 == NPOS {
                line.substr(pos0, NPOS)?
            } else {
                line.substr(pos0, pos1.wrapping_sub(pos0))?
            };
            let site_type = string_to_int(type_text, 0).map_err(DataError::Conversion)?;
            // `SetData(vecSitePos)` and `m_iCompleteSize = (int)vecSitePos.size()`.
            let segs_start = site_segs.len() as u32;
            let complete_size = if specific { residues.len() as i32 } else { 0 };
            for seg in site_data_segs(&mut residues) {
                site_segs.push(seg);
            }
            let segs = SiteRange {
                start: segs_start,
                len: site_segs.len() as u32 - segs_start,
            };
            let title = self.push_text(title);
            let motif = self.push_text(motif);
            sites.push(Site {
                title,
                motif,
                site_type,
                complete_size,
                segs,
                domain: domain as u32,
            });
            Ok(())
        })?;
        self.site_segs = site_segs;
        // Group the sites by domain, keeping the file order inside a domain (the order of
        // `m_lstSpecificFeatures`).
        sites.sort_by_key(|site| site.domain);
        let mut first = 0usize;
        while first < sites.len() {
            let domain = sites[first].domain;
            let mut end = first;
            while end < sites.len() && sites[end].domain == domain {
                end += 1;
            }
            let range = SiteRange {
                start: first as u32,
                len: (end - first) as u32,
            };
            if specific {
                self.domains[domain as usize].specific = range;
            } else {
                self.domains[domain as usize].generic = range;
            }
            first = end;
        }
        sites.shrink_to_fit();
        if specific {
            self.specific_sites = sites;
        } else {
            self.generic_sites = sites;
        }
        Ok(())
    }

    /// `m_lstSpecificFeatures` (`specific`) or `m_lstGenericFeatures` of a domain.
    pub(crate) fn sites(&self, domain: usize, specific: bool) -> &[Site] {
        let (range, table) = if specific {
            (self.domains[domain].specific, &self.specific_sites)
        } else {
            (self.domains[domain].generic, &self.generic_sites)
        };
        &table[range.start as usize..(range.start + range.len) as usize]
    }

    /// The ranges `from..=to` of a site's positions (`CSegSet` of the `TDomSite`, `ori_from`
    /// equal to `from`).
    pub(crate) fn site_segs(&self, site: &Site) -> &[[i32; 2]] {
        &self.site_segs[site.segs.start as usize..(site.segs.start + site.segs.len) as usize]
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.cpp:617-657
    fn load_bscores(&mut self, path: &str) -> Result<(), DataError> {
        if !file_exists(path) {
            return Err(DataError::Missing(path.to_string()));
        }
        const FILE: &str = "bitscore_specific.txt";
        for_each_line(path, |bytes| {
            if bytes.is_empty() || bytes[0] == b'#' {
                return Ok(());
            }
            let line = Line::new(bytes, FILE);
            let mut pos0;
            let mut pos1 = line.find(b'\t', 0);
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            let Some(dom) =
                self.find_domain_by_accession(line.substr(pos0, pos1.wrapping_sub(pos0))?)
            else {
                return Ok(());
            };
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            let text = if pos1 == NPOS {
                line.substr(pos0, NPOS)?
            } else {
                line.substr(pos0, pos1.wrapping_sub(pos0))?
            };
            self.domains[dom].min_bit_score =
                string_to_double(text).map_err(DataError::Conversion)?;
            Ok(())
        })
    }
}

/// `CSegSetTmpl::SetData(residues)`: the sorted, merged ranges of a list of positions.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/impl/segset_impl.hpp:25-55
// (`AddSeg` of a range that is more than one position from the last one appends it: the
// ranges come in order.)
fn site_data_segs(residues: &mut [i32]) -> Vec<[i32; 2]> {
    let mut segs: Vec<[i32; 2]> = Vec::new();
    if residues.is_empty() {
        return segs;
    }
    residues.sort_unstable();
    let (mut from, mut to) = (residues[0], residues[0]);
    for &residue in residues.iter() {
        if i64::from(residue) > i64::from(to) + 1 {
            segs.push([from, to]);
            from = residue;
            to = residue;
        } else if residue > to {
            to = residue;
        }
    }
    segs.push([from, to]);
    segs
}

/// The sign of NCBI's `int` difference `a - b` of two ids (`TSortByPssmId`, `TSortByClusterId`),
/// computed without its overflow. PD-LOSAT-NCBI-DEFECTS exception 9 (version 1.13): two ids whose
/// `int` values are 2^31 or more apart overflow NCBI's difference (undefined in C++); NCBI's
/// binary sorts its maps by the signed values (GCC folds `a - b < 0` to `a < b`) but its search
/// reads the wrapped difference, and can miss such an id (PS16a review oracle A2, run 2). LOSAT
/// compares the signed values in both, the order NCBI's maps have.
fn int_order(a: i32, b: i32) -> i32 {
    (i64::from(a) - i64::from(b)).signum() as i32
}

/// Sorts `items` by `order`, and items that `order` leaves equal by `position`, their place in the
/// input, which must increase strictly over `items` as given: the order of a stable sort. LOSAT's
/// sorts of rpsbproc's hits and data rows (PD-LOSAT-NCBI-DEFECTS 1.67, the order of tied hits). It
/// uses Rust's `sort_unstable_by`, which sorts in place, because `slice::sort_by` allocates a
/// buffer of up to the slice's length beyond its short slices; with `position` as the last key no
/// two items are equal, so the two give the same order.
pub(crate) fn sort_ties_in_input_order<T, P: Ord>(
    items: &mut [T],
    position: impl Fn(&T) -> P,
    mut order: impl FnMut(&T, &T) -> std::cmp::Ordering,
) {
    debug_assert!(items
        .windows(2)
        .all(|pair| position(&pair[0]) < position(&pair[1])));
    items.sort_unstable_by(|a, b| order(a, b).then_with(|| position(a).cmp(&position(b))));
}

/// `std::string::compare` as a sign: bytes compared as unsigned, a prefix first.
fn compare_bytes(a: &[u8], b: &[u8]) -> i32 {
    match a.cmp(b) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Greater => 1,
    }
}

/// `CPtrMap::Find` over `sorted` with `ord(row)` the comparison of the key with a row.
fn binary_search(sorted: &[u32], mut ord: impl FnMut(u32) -> i32) -> Option<u32> {
    let (mut begin, mut end) = (0usize, sorted.len());
    while end > begin {
        let mid = begin + ((end - begin) >> 1);
        let order = ord(sorted[mid]);
        if order == 0 {
            return Some(sorted[mid]);
        } else if order > 0 {
            begin = mid + 1;
        } else {
            end = mid;
        }
    }
    None
}

#[cfg(test)]
#[path = "data_tests.rs"]
mod tests;
