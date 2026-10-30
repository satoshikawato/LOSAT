//! Read-only reader of NCBI's `taxonomy4blast.sqlite3` for `-taxids` expansion.
//!
//! NCBI BLAST+ expands every user taxid to the taxid and all of its descendants
//! (`s_GetTaxIDList`, c++/src/algo/blast/blastinput/blast_args.cpp) with
//! `CTaxonomy4BlastSQLite::GetLeafNodeTaxids`
//! (c++/src/objtools/blast/seqdb_reader/tax4blastsqlite.cpp), which runs one
//! recursive SQL query per taxid on the SQLite file `taxonomy4blast.sqlite3`.
//! This module answers the same question without SQLite: it reads the file
//! format itself (<https://www.sqlite.org/fileformat2.html>), by positioned
//! reads only, and holds no page in memory between calls.
//!
//! ## What the real file holds (NCBI 2.17.0 swissprot download)
//!
//! ```text
//! CREATE TABLE TaxidInfo (taxid INTEGER PRIMARY KEY, parent INTEGER)
//! CREATE INDEX TaxidInfoCompositeIdx_parent ON TaxidInfo(parent,taxid)
//! ```
//!
//! 3,469,310 rows, page size 4096, UTF-8, journal mode `delete`. An index entry is
//! the record `(parent, taxid, rowid)` (`rowid` equals `taxid`, SQLite appends it
//! even to an index that names the rowid alias). Taxid 1 (root) has parent 1,
//! see [`TaxonomySqlite::descendants`] for what that means.
//!
//! ## Which failure happens where (NCBI)
//!
//! * [`TaxonomySqlite::open`] is the constructor: `sqlite3_open_v2` (lazy, the file
//!   is not read) and `x_SanityCheck`, whose first statement reads page 1 and
//!   the schema (a non-database file fails here with `[26] file is not a
//!   database`) and scans `sqlite_master`. The caller turns every `CException`
//!   of the constructor into the warning "The -taxids command line option
//!   requires additional data files. ..." and goes on without expansion.
//! * [`TaxonomySqlite::descendants`] is `GetLeafNodeTaxids`: the query is
//!   prepared at the first call and stepped to the end; a failure there
//!   escapes the per-id `try` of `s_GetTaxIDList` as "Invalid taxidlist file ".
//!
//! The file name resolution (`SeqDB_ResolveDbPath`, "Database '...' not found")
//! belongs to the caller; the reader takes the bytes.
//!
//! ## Not reproduced
//!
//! * SQLite parses the whole schema when it prepares the first statement; a
//!   file whose other `CREATE` statements are malformed fails in NCBI and works
//!   here.
//! * A database in WAL mode (header bytes 18/19 equal 2) is refused: the
//!   write-ahead log is a second file this reader does not see.
//! * Only the NCBI index definition `(parent,taxid)` is read (checked on the
//!   text of its `CREATE INDEX`); anything else is refused by `descendants`.
//! * A taxid column that is not an integer is refused (NCBI's file has integers only).

use std::cmp::Ordering;
use std::collections::HashSet;
use std::sync::{Arc, Mutex};

use super::lmdb_read::ReadAt;

/// NCBI's default file name (`CTaxonomy4BlastSQLite::kDefaultName`), used in the
/// texts of the sanity-check errors by [`TaxonomySqlite::open`].
pub const DEFAULT_NAME: &str = "taxonomy4blast.sqlite3";

const TABLE_NAME: &str = "TaxidInfo";
const INDEX_NAME: &str = "TaxidInfoCompositeIdx_parent";
/// Longest chain of b-tree levels followed; SQLite's own limit is 20.
const MAX_DEPTH: u32 = 24;

/// One column of a record (SQLite file format, section 2.1 "Record Format").
#[derive(Debug, Clone, PartialEq)]
enum Val {
    Null,
    Int(i64),
    Real(f64),
    Text(Vec<u8>),
    Blob(Vec<u8>),
    /// A text (`true`) or blob (`false`) whose body is not in the bytes read;
    /// enough to order it (after every number) but not to read it.
    Cut(bool),
}

/// Result of decoding a record prefix.
enum Rec {
    Cols(Vec<Val>),
    /// The bytes given end inside the header or a needed body.
    Short,
}

/// Text encoding of the database (file header offset 56).
#[derive(Clone, Copy, Debug, PartialEq)]
enum Enc {
    Utf8,
    Utf16Le,
    Utf16Be,
}

/// A `taxonomy4blast.sqlite3` opened for reading.
///
/// Reads go through a small direct-mapped cache of whole pages (at most 1 MiB),
/// so that the upper levels of the index are not read again for every node.
pub struct TaxonomySqlite<'r, R: ReadAt + ?Sized> {
    src: &'r R,
    page_size: usize,
    /// Page size minus the reserved bytes: the `U` of the cell formulas.
    usable: usize,
    page_count: u64,
    enc: Enc,
    index_root: u32,
    /// `Err` when the `CREATE INDEX` text is not NCBI's: reported by `descendants`.
    index_def: Result<(), String>,
    cache: Mutex<Vec<Option<(u32, Arc<Vec<u8>>)>>>,
}

/// A b-tree page and its header (file format section 1.6).
struct Page {
    data: Arc<Vec<u8>>,
    kind: u8,
    ncells: usize,
    /// Rightmost child of an interior page.
    right: u32,
    /// Offset in the page of the cell pointer array.
    ptrs: usize,
    /// Offset in the page before which no cell may start.
    min_cell: usize,
}

/// A decoded index cell: the left child (interior pages only), the first column and
/// the second column (`None` when not asked for, or when the first column ended the decoding).
struct IdxCell {
    left: u32,
    key0: Val,
    key1: Option<Val>,
}

// ---------------------------------------------------------------------------
// Primitives
// ---------------------------------------------------------------------------

fn corrupt(what: &str) -> String {
    format!("[11] database disk image is malformed: {}", what)
}

/// SQLite file format, section 1.8 "Variable-Length Integers": big-endian, 7 bits
/// per byte, a set high bit means more bytes follow, at most 9 bytes (the ninth
/// carries 8 bits). Returns the value and the length, `None` when `b` ends early.
fn varint(b: &[u8]) -> Option<(u64, usize)> {
    let mut v: u64 = 0;
    for i in 0..8 {
        let c = *b.get(i)?;
        v = (v << 7) | u64::from(c & 0x7f);
        if c & 0x80 == 0 {
            return Some((v, i + 1));
        }
    }
    let c = *b.get(8)?;
    Some(((v << 8) | u64::from(c), 9))
}

fn be_uint(b: &[u8]) -> u64 {
    b.iter().fold(0u64, |a, &c| (a << 8) | u64::from(c))
}

fn be_int(b: &[u8]) -> i64 {
    // sign-extend from the first byte
    let mut v: i64 = if b.first().is_some_and(|&c| c & 0x80 != 0) {
        -1
    } else {
        0
    };
    for &c in b {
        v = (v << 8) | i64::from(c);
    }
    v
}

/// Body length of a serial type (section 2.1), `Err` for the reserved types 10 and 11.
fn serial_len(t: u64) -> Result<u64, String> {
    Ok(match t {
        0 | 8 | 9 => 0,
        1 => 1,
        2 => 2,
        3 => 3,
        4 => 4,
        5 => 6,
        6 | 7 => 8,
        10 | 11 => return Err(corrupt("reserved serial type")),
        n if n % 2 == 0 => (n - 12) / 2,
        n => (n - 13) / 2,
    })
}

/// Decode a record (file format section 2.1) from its first bytes `buf` of `total`.
/// At most `max_cols` columns are decoded. A text or blob body that is not in `buf`
/// becomes `Val::Cut` when `cut_ok` (and ends the decoding); any other missing
/// byte gives `Rec::Short`.
fn parse_record(buf: &[u8], total: u64, max_cols: usize, cut_ok: bool) -> Result<Rec, String> {
    let (hsize, mut hp) = match varint(buf) {
        Some(x) => x,
        None => return Ok(Rec::Short),
    };
    if hsize < hp as u64 || hsize > total {
        return Err(corrupt("record header size"));
    }
    let mut types: Vec<u64> = Vec::new();
    while (hp as u64) < hsize && types.len() < max_cols {
        match varint(&buf[hp.min(buf.len())..]) {
            Some((t, n)) if (hp + n) as u64 <= hsize => {
                types.push(t);
                hp += n;
            }
            Some(_) => return Err(corrupt("record header column")),
            None => return Ok(Rec::Short),
        }
    }
    let mut cols = Vec::with_capacity(types.len());
    let mut pos = hsize;
    for t in types {
        let len = serial_len(t)?;
        let end = pos
            .checked_add(len)
            .ok_or_else(|| corrupt("record body size"))?;
        if end > total {
            return Err(corrupt("record body past the payload"));
        }
        if end > buf.len() as u64 {
            if cut_ok && t >= 12 {
                cols.push(Val::Cut(t % 2 == 1));
                break;
            }
            return Ok(Rec::Short);
        }
        let body = &buf[pos as usize..end as usize];
        cols.push(match t {
            0 => Val::Null,
            1..=6 => Val::Int(be_int(body)),
            7 => Val::Real(f64::from_bits(be_uint(body))),
            8 => Val::Int(0),
            9 => Val::Int(1),
            n if n % 2 == 0 => Val::Blob(body.to_vec()),
            _ => Val::Text(body.to_vec()),
        });
        pos = end;
    }
    Ok(Rec::Cols(cols))
}

/// Order of a record's first column against the integer `x`, as SQLite's index
/// comparison does (file format section 2.1 and the sort order of storage classes
/// in <https://www.sqlite.org/datatype3.html>, section 4.1): NULL, then numbers by
/// value, then text, then blobs.
fn cmp_val(v: &Val, x: i64) -> Ordering {
    match v {
        Val::Null => Ordering::Less,
        Val::Int(i) => i.cmp(&x),
        Val::Real(f) => {
            // exact comparison of a double with an integer (no rounding of x)
            if f.is_nan() || *f < -9.223372036854775808e18 {
                Ordering::Less
            } else if *f >= 9.223372036854775808e18 {
                Ordering::Greater
            } else {
                let fi = f.trunc();
                match (fi as i64).cmp(&x) {
                    Ordering::Equal => {
                        let frac = f - fi;
                        if frac > 0.0 {
                            Ordering::Greater
                        } else if frac < 0.0 {
                            Ordering::Less
                        } else {
                            Ordering::Equal
                        }
                    }
                    o => o,
                }
            }
        }
        Val::Text(_) | Val::Blob(_) | Val::Cut(_) => Ordering::Greater,
    }
}

/// UTF-8 text of a column in the database encoding (lossy on invalid bytes).
fn text_of(enc: Enc, b: &[u8]) -> String {
    match enc {
        Enc::Utf8 => String::from_utf8_lossy(b).into_owned(),
        Enc::Utf16Le | Enc::Utf16Be => {
            let units: Vec<u16> = b
                .chunks_exact(2)
                .map(|c| {
                    if enc == Enc::Utf16Le {
                        u16::from_le_bytes([c[0], c[1]])
                    } else {
                        u16::from_be_bytes([c[0], c[1]])
                    }
                })
                .collect();
            String::from_utf16_lossy(&units)
        }
    }
}

impl<'r, R: ReadAt + ?Sized> TaxonomySqlite<'r, R> {
    // -----------------------------------------------------------------------
    // Construction
    // -----------------------------------------------------------------------

    /// Open with the default file name in error texts, see [`Self::open_named`].
    pub fn open(src: &'r R) -> Result<Self, String> {
        Self::open_named(src, DEFAULT_NAME)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/tax4blastsqlite.cpp:77-106
    // ```c++
    // void CTaxonomy4BlastSQLite::x_SanityCheck()
    // {
    //     const map<string, string> kRequiredElements {
    //         { "table", "TaxidInfo" },
    //         { "index", "TaxidInfoCompositeIdx_parent" }
    //     };
    //     const string kSqlStmt =
    //         "SELECT COUNT(*) FROM sqlite_master WHERE type=? and name=?;";
    //     for (auto& [type, name] : kRequiredElements) {
    //         unique_ptr<CSQLITE_Statement> s(new CSQLITE_Statement(&*m_DbConn, kSqlStmt));
    //         s->Bind(1, type);
    //         s->Bind(2, name);
    //         s->Execute();
    //         if (s->Step()) {
    //             if (s->GetInt(0) != 1) {
    //                 CNcbiOstrstream oss;
    //                 oss << "Database '" << m_DbName << "' does not have " << type;
    //                 oss << " " << name << ". Please run the following command or ";
    //                 oss << "contact your system administrator to install it:" << endl;
    //                 oss << "update_blastdb.pl --decompress taxdb";
    //                 NCBI_THROW(CSeqDBException, eArgErr, CNcbiOstrstreamToString(oss));
    //             }
    //         } else {
    //             CNcbiOstrstream oss;
    //             oss << "Failed to check for " << type << " " << name << " in '";
    //             oss << m_DbName << "'";
    //             NCBI_THROW(CSeqDBException, eArgErr, CNcbiOstrstreamToString(oss));
    //         }
    //     }
    // }
    // ```
    // The `map` iterates in key order, so the index is checked before the table.
    // `SELECT COUNT(*)` always returns a row, so the "Failed to check" branch cannot
    // happen. `db_name` is the resolved path NCBI puts in the message. The header
    // checks are SQLite's for the first statement (file format section 1.3,
    // "The Database Header"); a failure there reads "[26] file is not a database"
    // as in the wrapper's `Error preparing statement for "<sql>": [26] ...`.
    /// Open the database bytes `src`: read the header and make NCBI's sanity check
    /// (exactly one `sqlite_master` row for the table `TaxidInfo` and one for the
    /// index `TaxidInfoCompositeIdx_parent`). `db_name` only appears in error texts.
    pub fn open_named(src: &'r R, db_name: &str) -> Result<Self, String> {
        const SQL: &str = "SELECT COUNT(*) FROM sqlite_master WHERE type=? and name=?;";
        let notdb = || {
            format!(
                "Error preparing statement for \"{}\": [26] file is not a database",
                SQL
            )
        };
        let len = src.len();
        let mut me = TaxonomySqlite {
            src,
            page_size: 4096,
            usable: 4096,
            page_count: 0,
            enc: Enc::Utf8,
            index_root: 0,
            index_def: Ok(()),
            cache: Mutex::new(Vec::new()),
        };
        // A file of size 0 is an empty database for SQLite: no sqlite_master rows.
        if len > 0 {
            if len < 100 {
                return Err(notdb());
            }
            let mut h = [0u8; 100];
            src.read_exact_at(&mut h, 0)?;
            if &h[0..16] != b"SQLite format 3\0" {
                return Err(notdb());
            }
            let ps = u16::from_be_bytes([h[16], h[17]]) as usize;
            let ps = if ps == 1 { 65536 } else { ps };
            if !(512..=65536).contains(&ps) || !ps.is_power_of_two() {
                return Err(notdb());
            }
            // bytes 18/19: file format write/read version, 1 = rollback journal, 2 = WAL
            if h[18] == 2 || h[19] == 2 {
                return Err("database is in WAL mode: the write-ahead log (taxonomy4blast.sqlite3-wal) is not available to this reader, copy the file after a checkpoint with journal_mode=DELETE".to_string());
            }
            if h[19] > 2 {
                return Err(notdb());
            }
            let usable = ps - h[20] as usize;
            if usable < 480 || h[21] != 64 || h[22] != 32 || h[23] != 32 {
                return Err(notdb());
            }
            me.enc = match be_uint(&h[56..60]) {
                1 => Enc::Utf8,
                2 => Enc::Utf16Le,
                3 => Enc::Utf16Be,
                _ => return Err(notdb()),
            };
            if be_uint(&h[44..48]) > 4 {
                return Err(notdb());
            }
            // header size in pages is valid only when the change counter matches
            let counter = be_uint(&h[24..28]);
            let dbsize = be_uint(&h[28..32]);
            let valid_for = be_uint(&h[92..96]);
            me.page_size = ps;
            me.usable = usable;
            me.page_count = if dbsize > 0 && counter == valid_for {
                dbsize
            } else {
                len / ps as u64
            };
            if me.page_count == 0 {
                return Err(notdb());
            }
        }
        let slots = (1_048_576 / me.page_size).clamp(8, 256);
        me.cache = Mutex::new(vec![None; slots]);
        // sqlite_master: type, name, tbl_name, rootpage, sql
        let (mut n_table, mut n_index) = (0u32, 0u32);
        let mut index_root = 0u32;
        let mut index_sql: Option<String> = None;
        if len > 0 {
            let mut budget = me.page_count + 64;
            let enc = me.enc;
            me.walk_table(1, 0, &mut budget, &mut |row: Vec<Val>| {
                let text = |i: usize| match row.get(i) {
                    Some(Val::Text(b)) => Some(text_of(enc, b)),
                    _ => None,
                };
                let (ty, name) = (text(0), text(1));
                if ty.as_deref() == Some("table") && name.as_deref() == Some(TABLE_NAME) {
                    n_table += 1;
                }
                if ty.as_deref() == Some("index") && name.as_deref() == Some(INDEX_NAME) {
                    n_index += 1;
                    if let Some(Val::Int(r)) = row.get(3) {
                        index_root = u32::try_from(*r).unwrap_or(0);
                    }
                    index_sql = text(4);
                }
                Ok(())
            })?;
        }
        for (ty, name, n) in [
            ("index", INDEX_NAME, n_index),
            ("table", TABLE_NAME, n_table),
        ] {
            if n != 1 {
                return Err(format!(
                    "Database '{}' does not have {} {}. Please run the following command or contact your system administrator to install it:\nupdate_blastdb.pl --decompress taxdb",
                    db_name, ty, name
                ));
            }
        }
        me.index_root = index_root;
        me.index_def = match index_sql {
            Some(s) => {
                // CREATE INDEX name ON TaxidInfo(parent,taxid): no sort order, collation or WHERE
                let n: String = s
                    .chars()
                    .filter(|c| !c.is_whitespace() && !"\"`[]'".contains(*c))
                    .collect::<String>()
                    .to_ascii_lowercase();
                if n.ends_with("ontaxidinfo(parent,taxid)") {
                    Ok(())
                } else {
                    Err(format!("unsupported definition of index {}: {:?} (this reader needs ON TaxidInfo(parent,taxid))", INDEX_NAME, s))
                }
            }
            None => Err(format!("index {} has no sql", INDEX_NAME)),
        };
        if me.index_root < 2 || u64::from(me.index_root) > me.page_count {
            me.index_def = Err(corrupt("root page of the index"));
        }
        Ok(me)
    }

    // -----------------------------------------------------------------------
    // Pages and cells
    // -----------------------------------------------------------------------

    /// Page `pgno` (1-based), from the cache or by one positioned read.
    fn read_page(&self, pgno: u32) -> Result<Arc<Vec<u8>>, String> {
        if pgno == 0 || u64::from(pgno) > self.page_count {
            return Err(corrupt("page number out of range"));
        }
        let slot = pgno as usize;
        {
            let c = self.cache.lock().unwrap_or_else(|e| e.into_inner());
            if let Some((n, p)) = &c[slot % c.len()] {
                if *n == pgno {
                    return Ok(Arc::clone(p));
                }
            }
        }
        let mut buf = vec![0u8; self.page_size];
        self.src
            .read_exact_at(&mut buf, (u64::from(pgno) - 1) * self.page_size as u64)?;
        let p = Arc::new(buf);
        let mut c = self.cache.lock().unwrap_or_else(|e| e.into_inner());
        let n = c.len();
        c[slot % n] = Some((pgno, Arc::clone(&p)));
        Ok(p)
    }

    /// SQLite file format, section 1.6 "B-tree Pages" (header of 8 bytes for a leaf,
    /// 12 for an interior page; page 1 starts its b-tree header after the 100-byte
    /// database header, but its cell offsets are still counted from the page start).
    fn page_hdr(&self, pgno: u32, want: &[u8]) -> Result<Page, String> {
        let data = self.read_page(pgno)?;
        let off = if pgno == 1 { 100usize } else { 0 };
        let b = &data[off..off + 12];
        let kind = b[0];
        if !want.contains(&kind) {
            return Err(corrupt("unexpected b-tree page type"));
        }
        let interior = kind == 0x02 || kind == 0x05;
        let hs = if interior { 12 } else { 8 };
        let ncells = u16::from_be_bytes([b[3], b[4]]) as usize;
        let min_cell = off + hs + 2 * ncells;
        if min_cell > self.usable {
            return Err(corrupt("cell pointer array past the page"));
        }
        let right = if interior {
            be_uint(&b[8..12]) as u32
        } else {
            0
        };
        Ok(Page {
            data,
            kind,
            ncells,
            right,
            ptrs: off + hs,
            min_cell,
        })
    }

    /// Offset in the page of cell `i`.
    fn cell_off(&self, p: &Page, i: usize) -> Result<usize, String> {
        let at = p.ptrs + 2 * i;
        let c = u16::from_be_bytes([p.data[at], p.data[at + 1]]) as usize;
        if c < p.min_cell || c >= self.usable {
            return Err(corrupt("cell offset outside the cell content area"));
        }
        Ok(c)
    }

    /// Bytes of a payload held on the b-tree page, file format section 1.6 "Cell
    /// Format": `U` usable size, `P` payload size; table leaf `X = U-35`, index
    /// `X = ((U-12)*64/255)-23`, `M = ((U-12)*32/255)-23`; if `P <= X` all is local,
    /// else `K = M+((P-M)%(U-4))` and local is `K` when `K <= X`, otherwise `M`.
    fn local_size(&self, p: u64, table_leaf: bool) -> u64 {
        let u = self.usable as u64;
        let x = if table_leaf {
            u - 35
        } else {
            (u - 12) * 64 / 255 - 23
        };
        let m = (u - 12) * 32 / 255 - 23;
        if p <= x {
            return p;
        }
        let k = m + (p - m) % (u - 4);
        if k <= x {
            k
        } else {
            m
        }
    }

    /// The first `want` bytes (at most `p`) of the payload that starts at offset
    /// `start` of page `pg`, of which `local` bytes are on the page, following the
    /// overflow chain (file format section 1.7 "Overflow Pages": a 4-byte next page
    /// number, then `U-4` bytes of payload).
    fn payload_prefix(
        &self,
        pg: &Page,
        start: usize,
        p: u64,
        local: u64,
        want: u64,
    ) -> Result<Vec<u8>, String> {
        let want = want.min(p);
        let local_end = (start as u64).saturating_add(local);
        if local_end > self.usable as u64 {
            return Err(corrupt("cell payload past the page"));
        }
        let local_end = local_end as usize;
        if want <= local {
            return Ok(pg.data[start..start + want as usize].to_vec());
        }
        if p > self.src.len() || local_end + 4 > self.usable {
            return Err(corrupt("payload larger than the file"));
        }
        let mut out = Vec::with_capacity(want as usize);
        out.extend_from_slice(&pg.data[start..local_end]);
        // the number of the first overflow page follows the local bytes
        let mut next = be_uint(&pg.data[local_end..local_end + 4]);
        let chunk = (self.usable - 4) as u64;
        let mut b4 = [0u8; 4];
        while (out.len() as u64) < want {
            if next == 0 || next > self.page_count {
                return Err(corrupt("overflow chain ends early"));
            }
            let at = (next - 1) * self.page_size as u64;
            self.src.read_exact_at(&mut b4, at)?;
            let take = (want - out.len() as u64).min(chunk) as usize;
            let old = out.len();
            out.resize(old + take, 0);
            self.src.read_exact_at(&mut out[old..], at + 4)?;
            next = be_uint(&b4);
        }
        Ok(out)
    }

    // -----------------------------------------------------------------------
    // sqlite_master: table b-tree
    // -----------------------------------------------------------------------

    /// In-order walk of a table b-tree (file format section 1.6: interior pages
    /// 0x05 hold `(left child, rowid)` cells and a rightmost child, leaf pages 0x0d
    /// hold `(payload size, rowid, payload)` cells); `f` gets every row decoded
    /// to its first five columns.
    fn walk_table(
        &self,
        pgno: u32,
        depth: u32,
        budget: &mut u64,
        f: &mut dyn FnMut(Vec<Val>) -> Result<(), String>,
    ) -> Result<(), String> {
        if depth > MAX_DEPTH || *budget == 0 {
            return Err(corrupt("b-tree too deep or cyclic"));
        }
        *budget -= 1;
        let p = self.page_hdr(pgno, &[0x05, 0x0d])?;
        for i in 0..p.ncells {
            let c = self.cell_off(&p, i)?;
            if p.kind == 0x05 {
                let b4 = p
                    .data
                    .get(c..c + 4)
                    .ok_or_else(|| corrupt("interior cell"))?;
                self.walk_table(be_uint(b4) as u32, depth + 1, budget, f)?;
            } else {
                let hb = &p.data[c..(c + 18).min(self.usable)];
                let (psize, n1) = varint(hb).ok_or_else(|| corrupt("cell header"))?;
                let (_rowid, n2) = varint(&hb[n1..]).ok_or_else(|| corrupt("cell header"))?;
                let local = self.local_size(psize, true);
                let data = self.payload_prefix(&p, c + n1 + n2, psize, local, psize)?;
                match parse_record(&data, psize, 5, false)? {
                    Rec::Cols(cols) => f(cols)?,
                    Rec::Short => return Err(corrupt("short record")),
                }
            }
        }
        if p.kind == 0x05 {
            self.walk_table(p.right, depth + 1, budget, f)?;
        }
        Ok(())
    }

    // -----------------------------------------------------------------------
    // The index b-tree
    // -----------------------------------------------------------------------

    /// Decode index cell `i` (file format section 1.6: interior pages 0x02 hold
    /// `(left child, payload size, payload)`, leaf pages 0x0a hold `(payload size,
    /// payload)`; the payload is the record `(parent, taxid, rowid)`). The second
    /// column is decoded when `want1`.
    fn index_cell(&self, p: &Page, i: usize, want1: bool) -> Result<IdxCell, String> {
        let c = self.cell_off(p, i)?;
        let b = &p.data[c..self.usable];
        let mut pos = 0usize;
        let mut left = 0u32;
        if p.kind == 0x02 {
            let l = b.get(0..4).ok_or_else(|| corrupt("index cell"))?;
            left = be_uint(l) as u32;
            pos = 4;
        }
        let (psize, n) = varint(&b[pos..]).ok_or_else(|| corrupt("index cell"))?;
        pos += n;
        let local = self.local_size(psize, false);
        let max_cols = if want1 { 2 } else { 1 };
        // fast path: the record prefix is on the page; else more bytes through the overflow chain
        let mut owned: Option<Vec<u8>> = None;
        let mut want = (b.len() - pos).min(local.min(usize::MAX as u64) as usize) as u64;
        loop {
            let cur: &[u8] = match &owned {
                Some(v) => v,
                None => &b[pos..pos + want as usize],
            };
            match parse_record(cur, psize, max_cols, true)? {
                Rec::Cols(cols) => {
                    let mut it = cols.into_iter();
                    let key0 = it.next().ok_or_else(|| corrupt("empty index record"))?;
                    return Ok(IdxCell {
                        left,
                        key0,
                        key1: it.next(),
                    });
                }
                Rec::Short => {
                    if want >= psize {
                        return Err(corrupt("short index record"));
                    }
                    want = (want * 4 + 64).min(psize);
                    owned = Some(self.payload_prefix(p, c + pos, psize, local, want)?);
                }
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/tax4blastsqlite.cpp:40-48
    // ```c++
    // const char* kSQLQuery = R"(
    // WITH RECURSIVE descendants_from(taxid, parent) AS
    //     (SELECT ?1, NULL
    //      UNION ALL
    //      SELECT T.taxid, T.parent
    //        FROM descendants_from D, TaxidInfo T
    //        WHERE D.taxid = T.parent)
    // SELECT taxid from descendants_from;
    // )";
    // ```
    // The recursive step `T.parent = D.taxid` is answered by the index
    // `TaxidInfoCompositeIdx_parent (parent, taxid)`: all entries whose first
    // column equals `x` are adjacent, ordered by taxid. The scan seeks the first
    // entry not below `(x)` and walks in order while the first column equals `x`;
    // in an index b-tree the cells of interior pages are entries too (file format
    // section 1.6), so they are visited between their left and right subtrees.
    // Returns `true` when an entry above `x` was met (the scan is over).
    fn index_scan(
        &self,
        pgno: u32,
        x: i64,
        depth: u32,
        budget: &mut u64,
        out: &mut Vec<i64>,
    ) -> Result<bool, String> {
        if depth > MAX_DEPTH || *budget == 0 {
            return Err(corrupt("b-tree too deep or cyclic"));
        }
        *budget -= 1;
        let p = self.page_hdr(pgno, &[0x02, 0x0a])?;
        // first cell whose first column is not below x
        let (mut lo, mut hi) = (0usize, p.ncells);
        while lo < hi {
            let mid = (lo + hi) / 2;
            if cmp_val(&self.index_cell(&p, mid, false)?.key0, x) == Ordering::Less {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        for j in lo..p.ncells {
            let cell = self.index_cell(&p, j, true)?;
            if p.kind == 0x02 && self.index_scan(cell.left, x, depth + 1, budget, out)? {
                return Ok(true);
            }
            match cmp_val(&cell.key0, x) {
                Ordering::Less => {}
                Ordering::Greater => return Ok(true),
                Ordering::Equal => match cell.key1 {
                    Some(Val::Int(t)) => out.push(t),
                    _ => return Err(corrupt("taxid column of the index is not an integer")),
                },
            }
        }
        if p.kind == 0x02 {
            return self.index_scan(p.right, x, depth + 1, budget, out);
        }
        Ok(false)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/tax4blastsqlite.cpp:108-126
    // ```c++
    // void CTaxonomy4BlastSQLite::GetLeafNodeTaxids(const int taxid, vector<int>& descendants)
    // {
    //     descendants.clear();
    //     if (taxid <= 0) return;
    //
    //     if (!m_Statement) {
    //         m_Statement.reset(new CSQLITE_Statement(&*m_DbConn, kSQLQuery));
    //     }
    //     m_Statement->Reset();
    //     m_Statement->ClearBindings();
    //     m_Statement->Bind(1, taxid);
    //     m_Statement->Execute();
    //
    //     while (m_Statement->Step()) {
    //         auto desc_taxid = m_Statement->GetInt(0);
    //         if (desc_taxid != taxid)
    //             descendants.push_back(desc_taxid);
    //     }
    // }
    // ```
    /// All proper descendants of `taxid` (children, their children, ...) in the order
    /// SQLite's recursive CTE emits them, which is breadth first: SQLite keeps the
    /// rows of the recursive table in a queue, outputs the row at the head and
    /// appends the rows its recursive step finds (children of a node in index
    /// order, i.e. ascending taxid). The order does not matter to NCBI: the caller
    /// (`s_GetTaxIDList`) puts the result into a `std::set`; tests compare sorted.
    ///
    /// `taxid <= 0` gives an empty list without a query; a taxid that is not in
    /// the file gives an empty list. A result taxid is `sqlite3_column_int`,
    /// i.e. the low 32 bits of the stored 64-bit value, and every result equal to
    /// `taxid` after that is dropped, as in NCBI.
    ///
    /// Taxid 1 is its own parent in NCBI's file. `UNION ALL` does not drop
    /// duplicates, so under SQLite the query for taxid 1 returns taxid 1 again
    /// and again for ever (verified: `sqlite3` streams the whole tree repeatedly);
    /// NCBI's `-taxids 1` therefore never finishes. LOSAT skips a taxid it has
    /// already visited (DI-8 DECISIONS 10, an approved difference of the
    /// PD-LOSAT-NCBI-DEFECTS kind), so the walk ends: for a start taxid in a parent
    /// cycle the result is every taxid reachable from it, without the start taxid
    /// itself, as for any other taxid. `descendants(1)` is the whole tree of the
    /// real file. `Err` stays for a damaged file (corrupt pages, depth limit, page
    /// budget, more rows than the file can contain).
    pub fn descendants(&self, taxid: i32) -> Result<Vec<i32>, String> {
        if taxid <= 0 {
            return Ok(Vec::new());
        }
        self.index_def.clone()?;
        let x = i64::from(taxid);
        // The queue and the output are the same list: `nodes[head..]` still to expand.
        let mut nodes: Vec<i64> = vec![x];
        // NCBI's `WITH RECURSIVE ... UNION ALL` (tax4blastsqlite.cpp) does not end for a start
        // taxid that lies in a parent cycle. LOSAT skips taxids it has already visited
        // (DI-8 DECISIONS 10, an approved difference of the PD-LOSAT-NCBI-DEFECTS kind).
        let mut seen: HashSet<i64> = HashSet::new();
        seen.insert(x);
        let mut kids: Vec<i64> = Vec::new();
        // an index entry takes at least 6 bytes of a file: more rows than that are a cycle
        let max_rows = self.src.len() / 6 + 16;
        let mut head = 0usize;
        while head < nodes.len() {
            let cur = nodes[head];
            head += 1;
            kids.clear();
            let mut budget = self.page_count + 64;
            self.index_scan(self.index_root, cur, 0, &mut budget, &mut kids)?;
            for &k in &kids {
                if seen.insert(k) {
                    nodes.push(k);
                }
            }
            if nodes.len() as u64 > max_rows {
                return Err(format!("the recursive query for taxid {} returns more rows than the file holds (damaged taxonomy table)", taxid));
            }
        }
        Ok(nodes
            .iter()
            .map(|&t| t as i32)
            .filter(|&t| t != taxid)
            .collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use std::path::PathBuf;

    /// Directory of the fixtures (small databases made with Python's `sqlite3` from the
    /// schema text of NCBI's real file) and of their `.ref` files (NCBI's query run by SQLite).
    fn fx_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/seqdb/taxsqlite")
    }

    fn fnv(ids: &[i32]) -> u64 {
        let mut h: u64 = 0xcbf29ce484222325;
        for i in ids {
            for b in format!("{},", i).bytes() {
                h = (h ^ u64::from(b)).wrapping_mul(0x100000001b3);
            }
        }
        h
    }

    fn load(name: &str) -> Vec<u8> {
        std::fs::read(fx_dir().join(name)).unwrap_or_else(|e| panic!("fixture {}: {}", name, e))
    }

    /// Compare every node of a fixture with the reference written by gen_fx.py.
    fn check_fixture(name: &str) {
        let data = load(&format!("{}.sqlite3", name));
        let db = TaxonomySqlite::open(&data).unwrap_or_else(|e| panic!("{}: open: {}", name, e));
        let refs = std::fs::read_to_string(fx_dir().join(format!("{}.ref", name))).unwrap();
        let mut n = 0;
        let mut full: HashMap<i64, Vec<i32>> = HashMap::new();
        for line in refs.lines() {
            let f: Vec<&str> = line.split(' ').collect();
            if f[0] == "wide" {
                continue;
            }
            if f[0] == "full" {
                full.insert(
                    f[1].parse().unwrap(),
                    f[2..]
                        .iter()
                        .filter(|s| !s.is_empty())
                        .map(|s| s.parse().unwrap())
                        .collect(),
                );
                continue;
            }
            let t: i64 = f[0].parse().unwrap();
            // taxids beyond i32 cannot be asked through the i32 API
            if t > i64::from(i32::MAX) {
                continue;
            }
            if f[1] == "LOOP" {
                // taxid 1 is its own parent: SQLite's query never ends; see `root_in_a_parent_cycle`
                continue;
            }
            let mut got = db
                .descendants(t as i32)
                .unwrap_or_else(|e| panic!("{} taxid {}: {}", name, t, e));
            got.sort_unstable();
            let (cnt, h): (usize, u64) = (f[1].parse().unwrap(), f[2].parse().unwrap());
            assert_eq!((got.len(), fnv(&got)), (cnt, h), "{} taxid {}", name, t);
            if let Some(want) = full.get(&t) {
                assert_eq!(&got, want, "{} taxid {}", name, t);
            }
            n += 1;
        }
        assert!(n > 50, "{}: only {} nodes compared", name, n);
        // taxid <= 0 and absent taxids
        assert_eq!(db.descendants(0).unwrap(), Vec::<i32>::new());
        assert_eq!(db.descendants(-5).unwrap(), Vec::<i32>::new());
        assert_eq!(db.descendants(i32::MAX).unwrap(), Vec::<i32>::new());
        println!("{}: {} nodes equal to sqlite3", name, n);
    }

    #[test]
    fn fixture_page_512() {
        check_fixture("p512");
    }

    #[test]
    fn fixture_utf16() {
        check_fixture("utf16");
    }

    /// Taxid 1 is its own parent (as in NCBI's file). NCBI's query for it never ends; LOSAT
    /// skips visited taxids (DI-8 DECISIONS 10), so the result is every other node reachable from 1.
    #[test]
    fn root_in_a_parent_cycle() {
        for name in ["p512", "utf16"] {
            let data = load(&format!("{name}.sqlite3"));
            let db = TaxonomySqlite::open(&data).unwrap();
            let refs = std::fs::read_to_string(fx_dir().join(format!("{name}.ref"))).unwrap();
            // The reference lists every node of the fixture (the "1 LOOP" line is the root); the
            // nodes that 1 does not reach are 2 (parent NULL), 3 (parent 0) and the two absent ids.
            let mut want: Vec<i32> = Vec::new();
            let mut listed = 0;
            for l in refs.lines() {
                let f: Vec<&str> = l.split(' ').collect();
                if f[0] == "wide" || f[0] == "full" {
                    continue;
                }
                listed += 1;
                let t: i64 = f[0].parse().unwrap();
                if [1, 2, 3, 12345, 999_999_999].contains(&t) {
                    continue;
                }
                // sqlite3_column_int keeps the low 32 bits
                want.push(t as i32);
            }
            want.sort_unstable();
            let mut got = db.descendants(1).unwrap();
            got.sort_unstable();
            assert_eq!(got.len(), want.len(), "{name}");
            assert_eq!(got, want, "{name}");
            assert_eq!(listed, want.len() + 5, "{name}");
            assert!(!got.contains(&1));
        }
    }

    /// Two nodes that are each other's parent (20 and 21), and the root: the walk ends and
    /// returns what NCBI's query would return before it started to repeat itself.
    #[test]
    fn parent_cycles_end() {
        let data = load("cycle.sqlite3");
        let db = TaxonomySqlite::open(&data).unwrap();
        let d = |t: i32| {
            let mut v = db.descendants(t).unwrap();
            v.sort_unstable();
            v
        };
        assert_eq!(d(1), vec![10, 11]);
        assert_eq!(d(10), vec![11]);
        assert_eq!(d(20), vec![21, 22]);
        assert_eq!(d(21), vec![20, 22]);
        assert_eq!(d(22), Vec::<i32>::new());
        assert_eq!(d(30), Vec::<i32>::new());
        assert_eq!(d(31), Vec::<i32>::new());
    }

    #[test]
    fn wide_taxid_truncates_like_getint() {
        // gen_fx.py stores taxid 5_000_000_000 as a child of the node named in the
        // `wide` line of the reference, and 5_000_000_001 as its child. GetInt keeps
        // the low 32 bits: 5_000_000_000 - 2^32 = 705_032_704. (The reference equality
        // of check_fixture already covers it; this names the values.)
        let data = load("p512.sqlite3");
        let db = TaxonomySqlite::open(&data).unwrap();
        let refs = std::fs::read_to_string(fx_dir().join("p512.ref")).unwrap();
        let wide: i32 = refs
            .lines()
            .find_map(|l| l.strip_prefix("wide "))
            .unwrap()
            .parse()
            .unwrap();
        let got = db.descendants(wide).unwrap();
        assert!(
            got.contains(&705_032_704) && got.contains(&705_032_705),
            "{:?}",
            got
        );
    }

    #[test]
    fn sanity_check_failures() {
        let data = load("noindex.sqlite3");
        let e = TaxonomySqlite::open_named(&data, "/db/taxonomy4blast.sqlite3")
            .err()
            .unwrap();
        assert_eq!(
            e,
            "Database '/db/taxonomy4blast.sqlite3' does not have index TaxidInfoCompositeIdx_parent. Please run the following command or contact your system administrator to install it:\nupdate_blastdb.pl --decompress taxdb"
        );
        let e = TaxonomySqlite::open(&data).err().unwrap();
        assert!(
            e.starts_with("Database 'taxonomy4blast.sqlite3' does not have index "),
            "{}",
            e
        );
        // an empty file is an empty database: both missing, the index is reported first
        let empty = load("empty.bin");
        let e = TaxonomySqlite::open(&empty).err().unwrap();
        assert!(
            e.contains("does not have index TaxidInfoCompositeIdx_parent"),
            "{}",
            e
        );
        // the table row removed from sqlite_master, the index there
        let d = load("notable.sqlite3");
        let e = TaxonomySqlite::open(&d).err().unwrap();
        assert!(
            e.contains("does not have table TaxidInfo. Please run"),
            "{}",
            e
        );
    }

    #[test]
    fn not_a_database() {
        let d = load("notdb.bin");
        let e = TaxonomySqlite::open(&d).err().unwrap();
        assert!(e.ends_with("[26] file is not a database"), "{}", e);
        let e = TaxonomySqlite::open(&b"SQLite format 3\0"[..])
            .err()
            .unwrap();
        assert!(e.contains("not a database"), "{}", e);
    }

    #[test]
    fn wal_refused() {
        let mut d = load("p512.sqlite3");
        d[18] = 2;
        d[19] = 2;
        let e = TaxonomySqlite::open(&d).err().unwrap();
        assert!(e.contains("WAL"), "{}", e);
    }

    #[test]
    fn page_in_memory_and_vec_sources() {
        // the two impls of ReadAt answer the same
        let v = load("p512.sqlite3");
        let a = TaxonomySqlite::open(&v).unwrap();
        let b = TaxonomySqlite::open(v.as_slice()).unwrap();
        let first: i32 = {
            let refs = std::fs::read_to_string(fx_dir().join("p512.ref")).unwrap();
            refs.lines()
                .filter(|l| {
                    !l.starts_with("full")
                        && !l.starts_with("wide")
                        && !l.contains("LOOP")
                        && !l.contains(" 0 ")
                })
                .map(|l| l.split(' ').next().unwrap().parse::<i64>().unwrap())
                .find(|&t| t < i64::from(i32::MAX))
                .unwrap() as i32
        };
        assert_eq!(a.descendants(first).unwrap(), b.descendants(first).unwrap());
        assert!(!a.descendants(first).unwrap().is_empty());
    }

    /// Random byte-damaged copies per fixture (the source test used 1500; kept short for a debug build).
    const DAMAGE_ROUNDS: usize = 200;

    /// Damaged files must give `Err` (or a result), never a panic or a hang.
    #[test]
    fn damaged_files_never_panic() {
        let mut seed: u64 = 0x1234_5678_9abc_def0;
        let mut rnd = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        for name in ["p512.sqlite3", "utf16.sqlite3"] {
            let orig = load(name);
            // the root, a node without children and the largest subtrees of the fixture
            let refs =
                std::fs::read_to_string(fx_dir().join(name.replace(".sqlite3", ".ref"))).unwrap();
            let mut big: Vec<(usize, i32)> = refs
                .lines()
                .filter(|l| !l.starts_with("wide") && !l.starts_with("full"))
                .filter_map(|l| {
                    let f: Vec<&str> = l.split(' ').collect();
                    Some((
                        f.get(1)?.parse().ok()?,
                        i32::try_from(f[0].parse::<i64>().ok()?).ok()?,
                    ))
                })
                .collect();
            big.sort_unstable_by(|a, b| b.cmp(a));
            let mut ids: Vec<i32> = vec![1, 2, 705_032_704];
            ids.extend(big.iter().take(3).map(|b| b.1));
            // truncations
            for cut in [
                0usize,
                50,
                99,
                100,
                101,
                511,
                1000,
                orig.len() / 3,
                orig.len() / 2,
                orig.len() - 1,
            ] {
                let d = &orig[..cut.min(orig.len())];
                if let Ok(db) = TaxonomySqlite::open(d) {
                    for &t in &ids {
                        let _ = db.descendants(t);
                    }
                }
            }
            // random byte damage, mostly in the first pages and the cell pointers
            let mut opens = 0;
            for round in 0..DAMAGE_ROUNDS {
                let mut d = orig.clone();
                for _ in 0..(1 + round % 6) {
                    let at = if rnd() % 3 == 0 {
                        (rnd() % 160) as usize
                    } else {
                        (rnd() % d.len() as u64) as usize
                    };
                    d[at] = rnd() as u8;
                }
                if let Ok(db) = TaxonomySqlite::open(&d) {
                    opens += 1;
                    for &t in &ids {
                        let _ = db.descendants(t);
                    }
                }
            }
            assert!(opens > 100, "{}: only {} damaged files opened", name, opens);
        }
    }

    // --- the real file ------------------------------------------------------------

    /// Positioned reads of a real file (pread), the way the integration will use it.
    struct FileSrc(std::fs::File, u64);

    impl ReadAt for FileSrc {
        fn len(&self) -> u64 {
            self.1
        }
        fn read_exact_at(&self, buf: &mut [u8], offset: u64) -> Result<(), String> {
            ReadAt::read_exact_at(&self.0, buf, offset)
        }
    }

    /// Needs the real NCBI file: `LOSAT_TEST_TAX_SQLITE=/path/to/taxonomy4blast.sqlite3` (3.4 million
    /// rows) and `LOSAT_TEST_TAX_SQLITE_EXPECTED`, the file that `tests/seqdb/tax_sqlite_ref.py`
    /// writes from it with Python's `sqlite3`: a line `T: sorted descendants` for 212 taxids (NCBI's
    /// query run by SQLite) and for taxid 1, whose NCBI query never ends (a walk of the whole tree).
    #[test]
    #[ignore = "needs LOSAT_TEST_TAX_SQLITE (the real taxonomy4blast.sqlite3) and LOSAT_TEST_TAX_SQLITE_EXPECTED"]
    fn real_file_equals_sqlite3() {
        let path = std::env::var("LOSAT_TEST_TAX_SQLITE").expect("LOSAT_TEST_TAX_SQLITE");
        let expected = std::env::var("LOSAT_TEST_TAX_SQLITE_EXPECTED")
            .expect("LOSAT_TEST_TAX_SQLITE_EXPECTED");
        let t0 = std::time::Instant::now();
        let f = std::fs::File::open(&path).unwrap();
        let file = FileSrc(f.try_clone().unwrap(), f.metadata().unwrap().len());
        let db = TaxonomySqlite::open_named(&file, &path).unwrap();
        println!("real file: {} bytes, open {:?}", file.1, t0.elapsed());
        let refs = std::fs::read_to_string(&expected).unwrap();
        let mut total = 0usize;
        let mut count = 0usize;
        for line in refs.lines() {
            let (t, rest) = line.split_once(':').unwrap();
            let t: i32 = t.parse().unwrap();
            let want: Vec<i32> = rest
                .split_whitespace()
                .map(|s| s.parse().unwrap())
                .collect();
            let t1 = std::time::Instant::now();
            // root: SQLite never finishes; LOSAT returns the whole tree
            let mut got = db
                .descendants(t)
                .unwrap_or_else(|e| panic!("taxid {}: {}", t, e));
            let el = t1.elapsed();
            got.sort_unstable();
            assert_eq!(got, want, "taxid {}", t);
            total += got.len();
            count += 1;
            if got.len() > 20000 || [1, 9606, 561, 10090].contains(&t) {
                println!("taxid {:>8}: {:>8} descendants in {:?}", t, got.len(), el);
            }
        }
        assert_eq!(count, 213, "{expected}: lines");
        println!(
            "{} taxids equal to sqlite3 ({} descendants in all)",
            count, total
        );
    }
}
