//! The CDD annotation data of rpsbproc (`COfflDomClusterData`): the six text files of NCBI's
//! `pub/mmdb/cdd/` directory, read as NCBI's offline loader reads them, and the lookups by
//! PSSM-ID, accession and cluster that the classification makes.
//!
//! NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.hpp:12-26
//! ```c++
//! class COfflDomClusterData: public CDomClusterIndex
//! {
//! 	...
//! 	static constexpr const char * __cdd_ids_file = "cddid.tbl";
//! 	static constexpr const char * __spec_feats_file = "cddannot.dat";
//! 	static constexpr const char * __gen_feats_file = "cddannot_generic.dat";
//! 	static constexpr const char * __bs_thrlds_file = "bitscore_specific.txt";
//! 	static constexpr const char * __cdtrack_file = "cdtrack.txt";
//! 	static constexpr const char * __domfam_link_file = "family_superfamily_links";
//! ```
//!
//! NCBI keeps a `TDomain` (or `TCluster`) per row with every column as a `std::string`; LOSAT
//! keeps the columns that `-t doms` prints or decides by (the PSSM-ID, accession, short name,
//! length, curated and structural flags, cluster and bit score threshold) with the texts in one
//! buffer. The definition lines, the hierarchy (`cdtrack.txt`) and the sites (`cddannot*.dat`)
//! are read as NCBI reads them, for their errors, but not kept: no byte of `-t doms` depends on
//! them (item PS16b keeps what its outputs read).
//! The lookups are NCBI's `CPtrMap`: the rows sorted by libstdc++'s `std::sort` with NCBI's
//! comparison and searched by NCBI's binary search, so that a key held by two rows finds the row
//! that NCBI finds.

use std::io::{BufRead, BufReader};

use crate::ncbi_exception::NcbiException;
use crate::ncbi_str::{printable_string, string_to_double, string_to_int, string_to_uint};
use crate::utils::std_sort::std_sort;

/// `INVALIDPSSMID`.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:26-30
// ```c++
// constexpr const PssmId_t INVALIDPSSMID = 100000000;
// constexpr const PssmId_t DOMSRCMUL = 1000000;
// constexpr const ArchId_t MINSPECARCHID = 10000000;
// constexpr const ClusterId_t INVALIDCLUSTERID = -1;
// constexpr const ClusterId_t SINGLEMEMBERCLUSTER = k_INT_MAX;
// ```
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
}

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
}

/// Why the data does not load.
#[derive(Debug)]
pub(crate) enum DataError {
    /// `CSimpleException` "Data file <path> does not exist" (rpsbproc: "Annotation data loading
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
    // libstdc++ (GCC 13): bits/basic_string.h:385-393
    // ```c++
    //       size_type
    //       _M_check(size_type __pos, const char* __s) const
    //       {
    // 	if (__pos > this->size())
    // 	  __throw_out_of_range_fmt(__N("%s: __pos (which is %zu) > "
    // 				       "this->size() (which is %zu)"),
    // 				   __s, __pos, this->size());
    // 	return __pos;
    //       }
    // ```
    // (`substr` calls `_M_check(__pos, "basic_string::substr")`, bits/basic_string.h:3163-3165;
    // NCBI's rpsbproc printed `Error: (106.6) Application's initialization failed:
    // basic_string::substr: __pos (which is 18446744073709551615) > this->size() (which is 10)`
    // for a `cdtrack.txt` line of an accession alone, PS16a review oracle A4-1 case `ct3`.)
    fn substr(&self, pos: usize, n: usize) -> Result<&'a [u8], DataError> {
        if pos > self.bytes.len() {
            return Err(DataError::Std(format!(
                "basic_string::substr: __pos (which is {pos}) > this->size() (which is {})",
                self.bytes.len()
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
            printable_string(self.bytes)
        ))
    }
}

/// `x_GetTrimPos` and the other trims of rpsbproc's library: a byte is white space when it is at
/// most `' '` as a signed `char` (so the bytes from 0x80 too).
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basealgo.cpp:18-28
// ```c++
// void x_GetTrimPos(const string & ori, size_t &pos0, size_t &pos1)
// {
// 	pos0 = 0;
// 	pos1 = ori.size();;
//
// 	while (pos0 < pos1 && ori[pos0] <= ' ')
// 		++pos0;
//
// 	while (pos1 > pos0 && ori[pos1 - 1] <= ' ')
// 		--pos1;
//
// }
// ```
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
// ```c++
//     template <typename TNumeric>
//     static typename enable_if< is_integral<TNumeric>::value && is_signed<TNumeric>::value && (sizeof(TNumeric) == sizeof(int) && !is_same<TNumeric, long>::value), TNumeric>::type
//     x_StringToNumeric(const CTempString str, TStringToNumFlags flags, int base)
//     {
//         return StringToInt(str, flags, base);
//     }
//     template <typename TNumeric>
//     static typename enable_if< is_integral<TNumeric>::value && is_unsigned<TNumeric>::value && (sizeof(TNumeric) == sizeof(unsigned int) && !is_same<TNumeric, unsigned long>::value), TNumeric>::type
//     x_StringToNumeric(const CTempString str, TStringToNumFlags flags, int base)
//     {
//         return StringToUInt(str, flags, base);
//     }
// ```
// The loader's `NStr::StringToNumeric<T>(text)` (flags 0, base 10) is `StringToUInt` for the
// unsigned ids and lengths, `StringToInt` for the signed ones and `StringToDouble` for the bit
// scores (ncbistr.hpp:3632-3637): `ncbi_str::string_to_uint`, `string_to_int` and
// `string_to_double`. Their `CStringException` (`eConvert`) is not caught by
// `catch (CSimpleException e)` (rpsbproc_main.cpp:654), so it leaves `Init` to the
// application's handler (`DataError::Conversion`).

/// `TDomSrcCount::ESrcIdx` and the signatures of the source databases.
// NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.cpp:117-140
// ```c++
// array<const char * , TDomSrcCount::TOTALSIGS > TDomSrcCount::DOMSRCSIGS = {"CD", "PFAM", "TIGR", "COG", "KOG", "SMART", "PRK", "CHL", "MTH", "PHA", "PLN", "PTZ"};
//
// TDomSrcCount::ESrcIdx TDomSrcCount::DomAccType(const string &acxn)
// {
// 	size_t sig = 0, acxn_len = acxn.size();
// 	while (sig < TOTALSIGS)
// 	{
//
// 		size_t siglen = strlen(TDomSrcCount::DOMSRCSIGS[sig]);
// 		if (siglen < acxn_len)
// 		{
// 			for (size_t cidx = 0; cidx < siglen; ++cidx)
// 			{
// 				if (toupper((unsigned char)acxn[cidx]) != (unsigned char)TDomSrcCount::DOMSRCSIGS[sig][cidx])
// 					goto labelNext;
// 			}
// 			return (ESrcIdx)sig;
// 		}
// 	labelNext:
// 		++sig;
// 	}
// 	return TOTALSIGS;
// }
// ```
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
// ```c++
// PssmId_t TDomSrcCount::ComputePseudoPssmId(const string &acxn)
// {
// 	PssmId_t pseudo_pssmid = 0;
// 	if (!acxn.empty())
// 	{
// 		ESrcIdx srctype = DomAccType(acxn);
// 		pseudo_pssmid = INVALIDPSSMID + (PssmId_t)(srctype * DOMSRCMUL);
//
//
// 		size_t idx = 0, idxEnd = acxn.size();
//
// 		while (idx < idxEnd && !isdigit(acxn[idx]))
// 			++idx;
//
// 		PssmId_t odr = 0;
// 		while (idx < idxEnd && isdigit(acxn[idx]))
// 		{
// 			odr = odr * 10 + (PssmId_t)(acxn[idx] - '0');
// 			++idx;
// 		}
//
// 		pseudo_pssmid += odr;
//
// 	}
// 	return pseudo_pssmid;
// }
// ```
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
    // ```c++
    // CPtrMap<TKey, TVal, TComp> :: Find(const TKey &k) const
    // {
    // 	std::size_t rbegin = 0, rend = m_ptr_vec.size(), dist = rend - rbegin;
    // 	while (dist > 0)
    // 	{
    // 		std::size_t rmid = rbegin + (dist >> 1);
    // 		int ord = m_comp(k, *m_ptr_vec[rmid]);
    // 		if (0 == ord)
    // 			return m_ptr_vec[rmid];	//found
    // 		else if (ord > 0) rbegin = rmid + 1;
    // 		else rend = rmid;
    // 		dist = rend - rbegin;
    // 	}
    // 	return nullptr;
    // }
    // ```
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/basedata.hpp:796-833
    // ```c++
    // 	struct TSortByPssmId
    // 	{
    // 		int operator () (const TCluster & v1, const TCluster & v2) const
    // 		{
    // 			return (int)v1.m_uiPssmId - (int)v2.m_uiPssmId;
    // 		}
    //
    // 		int operator () (PssmId_t pssmid, const TCluster & v) const
    // 		{
    // 			return (int)pssmid - (int)v.m_uiPssmId;
    // 		}
    // 	};
    // 	...
    // 	struct TSortByAcxn
    // 	{
    // 		...
    // 		int operator () (const std::string & acxn, const TCluster & v) const
    // 		{
    // 			return acxn.compare(v.m_strAccession);
    // 		}
    // 	};
    // ```
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
    // ```c++
    // 	struct TSortByClusterId
    // 	{
    // 		int operator () (const TCluster & v1, const TCluster & v2) const
    // 		{
    // 			return v1.m_iClusterId - v2.m_iClusterId;
    // 		}
    //
    // 		int operator () (ClusterId_t clstid, const TCluster & v) const
    // 		{
    // 			return (int)clstid - (int)v.m_iClusterId;
    // 		}
    // 	};
    // ```
    pub(crate) fn find_cluster(&self, cluster_id: i32) -> Option<usize> {
        let found = binary_search(&self.cluster_by_id, |index| {
            int_order(cluster_id, self.clusters[index as usize].cluster_id)
        })?;
        Some(found as usize)
    }

    /// `COfflDomClusterData::LoadData_real`.
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.cpp:66-83
    // ```c++
    // void COfflDomClusterData::LoadData_real(const string &cddids, const string &spfeats, const string &genfeats, const string &bscores, const string &cdtrack, const string &famlinks)
    // {
    // 	...
    // 	if (m_CDStore.IsEmpty())	//not loaded yet, load it
    // 	{
    // 		ifstream dfstream;
    // 		string buf;
    //
    // 		CFile srcfile;
    // 		srcfile.Reset(cddids);
    // 		...
    // 		if (!srcfile.Exists())
    // 			throw CSimpleException(__FILE__, __LINE__, "Data file " + srcfile.GetPath() + " does not exist");
    // ```
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
    // ```c++
    // 			while (dfstream.good())
    // 			{
    // 				buf.clear();
    // 				getline(dfstream, buf);
    // 				//LTrimString(buf);
    // 				if (!buf.empty() && '#' != buf[0])
    // 				{
    //
    // 					size_t pos0 = 0, pos1 = buf.find('\t', pos0);
    //
    // 					string pssmid_lit = buf.substr(pos0, pos1 - pos0);
    // 					LTrimString(pssmid_lit);
    //
    // 					PssmId_t pssmid = 0;
    //
    // 					if (!pssmid_lit.empty())
    // 						pssmid = NStr::StringToNumeric<PssmId_t> (pssmid_lit);
    //
    // 					pos0 = pos1 + 1;
    // 					pos1 = buf.find('\t', pos0);
    // 					string tmp(buf.substr(pos0, pos1 - pos0));
    // 					TrimString(tmp);
    //
    // 					if (0 == pssmid)	//make pseudo pssmid
    // 						pssmid = TDomSrcCount::ComputePseudoPssmId(tmp);
    //
    // 					if ('c' == tmp[0] && 'l' == tmp[1] && (tmp[2] >= '0' && tmp[2] <= '9'))	//cluster
    // 					{
    //
    // 						// is cluster
    // 						TCluster * p_clst = m_ClusterStore.Append();
    // 						p_clst->m_uiPssmId = pssmid;
    // 						p_clst->m_strAccession = move(tmp);
    // 						p_clst->m_iClusterId = NStr::StringToNumeric<ClusterId_t> (p_clst->m_strAccession.substr(2));
    //
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);
    // 						p_clst->m_strShortName = buf.substr(pos0, pos1 - pos0);
    // 						TrimString(p_clst->m_strShortName);
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);
    // 						p_clst->m_strDefline = buf.substr(pos0, pos1 - pos0);
    // 						TrimString(p_clst->m_strDefline);
    // 						...
    // 						clst_ptrs.emplace_back(p_clst);
    // 					}
    // 					else	//cd
    // 					{
    // 						TDomain * p_cd = m_CDStore.Append();
    // 						p_cd->m_uiPssmId = pssmid;
    // 						p_cd->m_strAccession = buf.substr(pos0, pos1 - pos0);
    // 						TrimString(p_cd->m_strAccession);
    // 						if (isdigit(p_cd->m_strAccession[2]))	//cdxxx, sdxxx or clxxx
    // 						{
    // 							if ('d' == p_cd->m_strAccession[1])	//cd or sd
    // 							{
    // 								if ('c' == p_cd->m_strAccession[0])	//curated
    // 								{
    // 									p_cd->m_bCurated = true;
    // 									//cd_ptrs.emplace_back(p_cd);
    // 								}
    //
    // 								else if ('s' == p_cd->m_strAccession[0])	//curated
    // 									p_cd->m_bIsStructDom = true;
    // 							}
    // 						}
    //
    //
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);
    // 						p_cd->m_strShortName = buf.substr(pos0, pos1 - pos0);
    // 						TrimString(p_cd->m_strShortName);
    //
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);
    // 						p_cd->m_strDefline = buf.substr(pos0, pos1 - pos0);
    // 						TrimString(p_cd->m_strDefline);
    //
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);
    //
    // 						// -- assume length is the last field
    // 						while (string::npos != pos1)
    // 						{
    // 							pos0 = pos1 + 1;
    // 							pos1 = buf.find('\t', pos0);
    // 						}
    //
    // 						string tmp1(buf.substr(pos0));
    // 						TrimString(tmp1);
    //
    // 						p_cd->m_uiLength = NStr::StringToNumeric<SeqLen_t> (tmp1);
    // 						dom_ptrs.emplace_back(p_cd);
    // 					}
    // 				}
    // 			}
    // 			dfstream.close();
    //
    // 			// -- now create maps,
    // 			clst_ptrs.shrink_to_fit();	//cluster
    // 			//cd_ptrs.shrink_to_fit();	//curated
    // 			dom_ptrs.shrink_to_fit();	//all domains
    //
    // 			SetPssmId2CdPtrs(dom_ptrs);
    // 			SetPssmId2ClusterPtrs(clst_ptrs);
    //
    // 			SetAcxn2CdPtrs(move(dom_ptrs));
    // 			SetClusterId2ClusterPtrs(move(clst_ptrs));
    // ```
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
            });
            Ok(())
        })?;
        self.domains = domains;
        self.clusters = clusters;
        // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/ptrmap.hpp:167-179
        // ```c++
        // void CPtrMap<TKey, TVal, TComp> :: Reset(const TSORTED_PTRS & ptrs)
        // {
        // 	m_ptr_vec = ptrs;
        // 	if (m_ptr_vec.size() > 1) sort(m_ptr_vec.begin(), m_ptr_vec.end(), m_sort_cmp);
        // }
        // ```
        // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/ptrmap.hpp:63-70
        // ```c++
        // 	struct TSortComp
        // 	{
        // 		const TComp &m_comp;
        // 		bool operator() (CONST_VAL_PTR p1, CONST_VAL_PTR p2)
        // 		{
        // 			return m_comp(*p1, *p2) < 0;
        // 		}
        // ```
        let domains = &self.domains;
        let clusters = &self.clusters;
        let text = &self.text;
        let span = |s: Span| &text[s.start as usize..(s.start + s.len) as usize];
        let mut by_pssm: Vec<u32> = (0..domains.len() as u32).collect();
        std_sort(&mut by_pssm, &mut |&a: &u32, &b: &u32| {
            int_order(
                domains[a as usize].pssm_id as i32,
                domains[b as usize].pssm_id as i32,
            ) < 0
        });
        let mut by_acc: Vec<u32> = (0..domains.len() as u32).collect();
        std_sort(&mut by_acc, &mut |&a: &u32, &b: &u32| {
            compare_bytes(
                span(domains[a as usize].accession),
                span(domains[b as usize].accession),
            ) < 0
        });
        let mut cl_by_id: Vec<u32> = (0..clusters.len() as u32).collect();
        std_sort(&mut cl_by_id, &mut |&a: &u32, &b: &u32| {
            int_order(
                clusters[a as usize].cluster_id,
                clusters[b as usize].cluster_id,
            ) < 0
        });
        self.domain_by_pssm_id = by_pssm;
        self.domain_by_accession = by_acc;
        self.cluster_by_id = cl_by_id;
        Ok(())
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.cpp:217-308
    // ```c++
    // 		srcfile.Reset(cdtrack);
    // 		if (!srcfile.Exists())
    // 			throw CSimpleException(__FILE__, __LINE__, "Data file " + srcfile.GetPath() + " does not exist");
    // 		else
    // 		{
    // 			dfstream.open(srcfile.GetPath().c_str());
    //
    // 			while (dfstream.good())
    // 			{
    // 				buf.clear();
    // 				getline(dfstream, buf);
    // 				LTrimString(buf);
    //
    // 				if (!buf.empty() && '#' != buf[0])
    // 				{
    // 					// -- accession
    // 					size_t pos0 = 0, pos1 = buf.find(' ', pos0);
    //
    // 					string currAcc = buf.substr(pos0, pos1 - pos0);
    //
    // 					TrimString(currAcc);
    // 					TDomain * p_dom = const_cast< TDomain* > (FindCdInfo(currAcc));
    // 					...
    // 					if (nullptr != p_dom)
    // 					{
    // 						// -- shortname - not used
    // 						while(buf[pos1] <= ' ') ++pos1;
    // 						pos0 = pos1;
    // 						pos1 = buf.find(' ', pos0);
    //
    // 						// -- // -- skip short-name. to go pssmid
    // 						while(' ' == buf[pos1]) ++pos1;
    //
    //
    // 						if (isdigit(buf[pos1]))	//is pssmid
    // 						{
    // 							pos0 = pos1;
    // 							pos1 = buf.find(' ', pos0);
    // 							// -- parent acc
    // 							while(buf[pos1] <= ' ') ++pos1;
    // 							//pssmid = NStr::StringToNumeric<PssmId_t> (buf.substr(pos0, pos1 - pos0));
    // 						}
    // 						//	otherwise this should already goes to next field
    // 						pos0 = pos1;
    // 						pos1 = buf.find(' ', pos0);
    //
    // 						string parentAcc = buf.substr(pos0, pos1 - pos0);
    //
    // 						TrimString(parentAcc);
    // 						if ("N/A" != parentAcc)
    // 						{
    // 							if (parentAcc == p_dom->m_strAccession)
    // 								p_dom->m_uiHierarchyParent = p_dom->m_uiPssmId;
    // 							else
    // 							{
    // 								const TDomain * ppcd = FindCdInfo(parentAcc);
    // 								if (nullptr != ppcd)
    // 									p_dom->m_uiHierarchyParent = ppcd->m_uiPssmId;
    // 							}
    // 						}
    //
    // 						// -- root acc
    // 						while(buf[pos1] <= ' ') ++pos1;
    // 						pos0 = pos1;
    // 						pos1 = buf.find(' ', pos0);
    // 						string rootAcc = buf.substr(pos0, pos1 - pos0);
    // 						TrimString(rootAcc);
    //
    // 						if (rootAcc == parentAcc)
    // 							p_dom->m_uiHierarchyRoot = p_dom->m_uiHierarchyParent;
    // 						else if (rootAcc == p_dom->m_strAccession)
    // 							p_dom->m_uiHierarchyRoot = p_dom->m_uiPssmId;
    // 						else
    // 						{
    // 							const TDomain * prcd = FindCdInfo(rootAcc);
    //
    // 							if (nullptr != prcd)
    // 								p_dom->m_uiHierarchyParent = prcd->m_uiPssmId;
    //
    // 						}
    // 					}
    // 				}
    // 			}	//dfstream
    // ```
    // The parent and the root set only `m_uiHierarchyParent` and `m_uiHierarchyRoot`, which no
    // output of `-t doms` reads (the doms path reads the length, threshold, flags, cluster,
    // accession and short name of a domain), so LOSAT does not keep them; what this loader does
    // for `-t doms` is what its `std::string` operations do on each line. A `find` that fails
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
            if self.find_domain_by_accession(current).is_none() {
                return Ok(());
            }
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
            // `parentAcc`
            line.substr(pos0, pos1.wrapping_sub(pos0))?;
            pos1 = skip_le_space(pos1)?;
            pos0 = pos1;
            pos1 = line.find(b' ', pos0);
            // `rootAcc`
            line.substr(pos0, pos1.wrapping_sub(pos0))?;
            if line.read_before.get() {
                return Err(line.unsupported("a read of the byte before the line (buf[npos])"));
            }
            Ok(())
        })
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.cpp:310-386
    // ```c++
    // 		srcfile.Reset(famlinks);
    //
    // 		if (!srcfile.Exists())
    // 			throw CSimpleException(__FILE__, __LINE__, "Data file " + srcfile.GetPath() + " does not exist");
    // 		else
    // 		{
    // 			dfstream.open(srcfile.GetPath().c_str());
    //
    // 			while (dfstream.good())
    // 			{
    // 				buf.clear();
    // 				getline(dfstream, buf);
    // 				LTrimString(buf);
    // 				if (!buf.empty() && '#' != buf[0])
    // 				{
    // 					size_t pos0 = 0, pos1 = buf.find('\t', pos0);
    // 					// -- acxn
    // 					TDomain * p_dom = const_cast< TDomain* > (FindCdInfo(buf.substr(pos0, pos1 - pos0)));
    //
    // 					if (nullptr != p_dom)
    // 					{
    // 						// -- pssmid, unused
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);
    // 						...
    // 						//superfamily accession
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);
    //
    // 						string sf_acxn = buf.substr(pos0, pos1 - pos0);
    //
    //
    // 						if ('c' ==sf_acxn[0] && 'l' == sf_acxn[1])	//superfamily
    // 						{
    // 							ClusterId_t clstid = NStr::StringToNumeric<ClusterId_t> (buf.substr(pos0+2, pos1-pos0-2));
    //
    // 							const TCluster *p_clst = FindClusterInfo(clstid);
    //
    // 							if (nullptr != p_clst)
    // 								p_dom->m_iClusterId = clstid;
    // 						}
    // 						else
    // 						{
    // 							p_dom->m_iClusterId = SINGLEMEMBERCLUSTER;
    // 						}
    // 						...
    // 						p_dom->m_bMultiDom = false;
    // 					}
    // ```
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
    // ```c++
    // 		TDomSite featValue;
    // 		// -- specific site annotations
    // 		srcfile.Reset(spfeats);
    // 		...
    // 				buf.clear();
    // 				getline(dfstream, buf);
    // 				LTrimString(buf);
    // 				if (!buf.empty() && '#' != buf[0])
    // 				{
    // 					size_t pos0 = 0, pos1 = buf.find('\t', pos0);	//pssmid, skip
    //
    //
    // 					pos0 = pos1 + 1;
    // 					pos1 = buf.find('\t', pos0);	//accession
    //
    //
    // 					string tmp(buf.substr(pos0, pos1 - pos0));
    // 					TrimString(tmp);
    // 					TDomain * p_dom = const_cast< TDomain* > (FindCdInfo(tmp));
    // 					...
    // 					if (nullptr != p_dom)
    // 					{
    // 						p_dom->m_lstSpecificFeatures.emplace_back(featValue);
    // 						TDomSite &tgt =  *(p_dom->m_lstSpecificFeatures.rbegin());
    //
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);	//shortname
    //
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);	//index
    //
    // 						tmp = buf.substr(pos0, pos1 - pos0);
    // 						TrimString(tmp);
    // 						tgt.m_iIndex = NStr::StringToNumeric<int>(tmp);
    //
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);	//title
    //
    // 						tgt.m_strTitle = buf.substr(pos0, pos1 - pos0);
    // 						TrimString(tgt.m_strTitle);
    //
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);	//Motif
    //
    // 						if ('0' != buf[pos0])
    // 						{
    // 							//tgt.SetMotifStr(buf.substr(pos0, pos1 - pos0).c_str());
    // 							tgt.m_strMotif = buf.substr(pos0, pos1 - pos0);
    // 							TrimString(tgt.m_strMotif);
    // 						}
    //
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);	//structure based evidence present
    //
    // 						if ('1' == buf[pos0])
    // 							tgt.m_flags |= TDomSite::STRUCTURE_BASED_EVIDENCE;
    // 						...
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);	//site coordinates
    //
    // 						vector<TSignedSeqPos> vecSitePos;	//to store coordinates
    //
    // 						vecSitePos.reserve((pos1 - pos0) * 3);	//estimated
    //
    // 						size_t pos_c = buf.find(',', pos0);
    //
    // 						while (pos_c < pos1)
    // 						{
    // 							vecSitePos.emplace_back(NStr::StringToNumeric<SeqPos_t>(buf.substr(pos0, pos_c - pos0)));
    // 							pos0 = pos_c + 1;
    // 							pos_c = buf.find(',', pos0);
    // 						}
    //
    // 						// -- last one
    // 						vecSitePos.emplace_back(NStr::StringToNumeric<SeqPos_t>(buf.substr(pos0, pos1 - pos0)));
    //
    // 						tgt.SetData(vecSitePos);
    //
    // 						tgt.m_iCompleteSize = (int)vecSitePos.size();
    //
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);	//type
    //
    // 						if (string::npos == pos1)	//no other field
    // 							tgt.m_iType = NStr::StringToNumeric<int>(buf.substr(pos0));
    // 						else	//assume has motif string added
    // 							tgt.m_iType = NStr::StringToNumeric<int>(buf.substr(pos0, pos1 - pos0));
    // 					}
    // ```
    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.cpp:518-609
    // ```c++
    // 			while (dfstream.good())
    // 			{
    // 				buf.clear();
    // 				getline(dfstream, buf);
    // 				if (!buf.empty() && '#' != buf[0])
    // 				{
    // 					size_t pos0 = 0, pos1 = buf.find('\t', pos0);	//pssmid -- skip
    //
    // 					pos0 = pos1 + 1;
    // 					pos1 = buf.find('\t', pos0);	//accession
    //
    // 					//PssmId_t pssmid = NStr::StringToNumeric<PssmId_t> (buf.substr(pos0, pos1 - pos0));
    // 					TDomain * p_dom = const_cast< TDomain* > (FindCdInfo(buf.substr(pos0, pos1-pos0)));
    // 					if (nullptr != p_dom)
    // 					{
    // 						...
    // 						tgt.m_iIndex = NStr::StringToNumeric<int>(buf.substr(pos0, pos1 - pos0));
    //
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);	//title
    //
    // 						tgt.m_strTitle = buf.substr(pos0, pos1 - pos0);
    // 						...
    // ```
    // The generic file is read without the trims of the specific one (the line, the accession,
    // the index and the title keep their white space). A coordinates field at the line's end
    // (`pos1` is `npos`) makes NCBI reserve `(npos - pos0) * 3` positions, more than a
    // `vector<int>` can hold whatever `pos0` is, which throws `std::length_error`:
    // libstdc++ (GCC 13): bits/vector.tcc:64-71
    // ```c++
    //   template<typename _Tp, typename _Alloc>
    //     _GLIBCXX20_CONSTEXPR
    //     void
    //     vector<_Tp, _Alloc>::
    //     reserve(size_type __n)
    //     {
    //       if (__n > this->max_size())
    // 	__throw_length_error(__N("vector::reserve"));
    // ```
    // (NCBI's rpsbproc printed `Error: (106.6) Application's initialization failed:
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
            if self.find_domain_by_accession(accession).is_none() {
                return Ok(());
            }
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
            line.substr(pos0, pos1.wrapping_sub(pos0))?;
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            if line.at(pos0)? != b'0' {
                line.substr(pos0, pos1.wrapping_sub(pos0))?;
            }
            for _evidence_flag in 0..3 {
                pos0 = pos1.wrapping_add(1);
                pos1 = line.find(b'\t', pos0);
                line.at(pos0)?;
            }
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            if pos1 == NPOS {
                return Err(DataError::Std("vector::reserve".to_string()));
            }
            let mut pos_c = line.find(b',', pos0);
            while pos_c < pos1 {
                string_to_int(line.substr(pos0, pos_c - pos0)?, 0)
                    .map_err(DataError::Conversion)?;
                pos0 = pos_c + 1;
                pos_c = line.find(b',', pos0);
            }
            string_to_int(line.substr(pos0, pos1.wrapping_sub(pos0))?, 0)
                .map_err(DataError::Conversion)?;
            pos0 = pos1.wrapping_add(1);
            pos1 = line.find(b'\t', pos0);
            let type_text = if pos1 == NPOS {
                line.substr(pos0, NPOS)?
            } else {
                line.substr(pos0, pos1.wrapping_sub(pos0))?
            };
            string_to_int(type_text, 0).map_err(DataError::Conversion)?;
            Ok(())
        })
    }

    // NCBI reference (RpsbProc 0.5.1, 2025-05-13): RpsbProc/src/common/offl_sparcle_data.cpp:617-657
    // ```c++
    // 		// -- read bitscore thresholds
    // 		srcfile.Reset(bscores);
    //
    // 		if (!srcfile.Exists())
    // 			throw CSimpleException(__FILE__, __LINE__, "Data file " + srcfile.GetPath() + " does not exist");
    // 		else
    // 		{
    //
    // 			dfstream.open(srcfile.GetPath().c_str());
    // 			while (dfstream.good())
    // 			{
    // 				buf.clear();
    // 				getline(dfstream, buf);
    // 				if (!buf.empty() && '#' != buf[0])
    // 				{
    // 					size_t pos0 = 0, pos1 = buf.find('\t', pos0);	//pssmid - skip
    // 					// -- accession
    // 					pos0 = pos1 + 1;
    // 					pos1 = buf.find('\t', pos0);	//accession
    // 					...
    // 					TDomain * p_dom = const_cast< TDomain* > (FindCdInfo(buf.substr(pos0, pos1 - pos0)));
    //
    // 					if (nullptr != p_dom)
    // 					{
    // 						// -- score
    // 						pos0 = pos1 + 1;
    // 						pos1 = buf.find('\t', pos0);
    //
    // 						if (string::npos == pos1)
    // 							p_dom->m_dMinBitScore = NStr::StringToNumeric<double>(buf.substr(pos0));
    // 						else
    // 							p_dom->m_dMinBitScore = NStr::StringToNumeric<double>(buf.substr(pos0, pos1 - pos0));
    //
    //
    // 					}
    // ```
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

/// The sign of NCBI's `int` difference `a - b` of two ids (`TSortByPssmId`, `TSortByClusterId`),
/// computed without its overflow. PD-LOSAT-NCBI-DEFECTS exception 9 (version 1.13): two ids whose
/// `int` values are 2^31 or more apart overflow NCBI's difference (undefined in C++); NCBI's
/// binary sorts its maps by the signed values (GCC folds `a - b < 0` to `a < b`) but its search
/// reads the wrapped difference, and can miss such an id (PS16a review oracle A2, run 2). LOSAT
/// compares the signed values in both, the order NCBI's maps have.
fn int_order(a: i32, b: i32) -> i32 {
    (i64::from(a) - i64::from(b)).signum() as i32
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
