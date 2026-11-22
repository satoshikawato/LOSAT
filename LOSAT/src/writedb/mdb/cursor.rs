// Derived from LMDB 0.9.24 (mdb.c, lmdb.h) as bundled in NCBI BLAST+ 2.17.0,
// c++/src/util/lmdb/, which carries this notice:
//
// Copyright 2011-2019 Howard Chu, Symas Corp.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted only as authorized by the OpenLDAP
// Public License.
//
// A copy of this license is available in the file LICENSE in the
// top-level directory of the distribution or, alternatively, at
// <http://www.OpenLDAP.org/license.html>.
//
// This code is derived from btree.c written by Martin Hedenfalk.
//
// Copyright (c) 2009, 2010 Martin Hedenfalk <martin@bzero.se>
//
// Permission to use, copy, modify, and distribute this software for any
// purpose with or without fee is hereby granted, provided that the above
// copyright notice and this permission notice appear in all copies.
//
// THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
// WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
// MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
// ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
// WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
// ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
// OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
//
// (In LOSAT, the OpenLDAP Public License 2.8 is in THIRD_PARTY_NOTICES, section 2 (LMDB).)

//! Cursors, searches, puts, deletes, page splits and the freelist of the LMDB port
//! (`c++/src/util/lmdb/mdb.c`). The cursors of a call live in `Env::cursors` (a stack: a
//! cursor that mdb.c keeps on the C stack is pushed and truncated in the same order); the
//! cursors linked in `txn->mt_cursors[dbi]` are in `Env::tracked`.

use super::midl::*;
use super::page::*;
use super::*;

impl Env {
    // ---- cursor records --------------------------------------------------------------------

    /// `*mc->mc_db`.
    fn db(&self, c: usize) -> &MdbDb {
        match self.cursors[c].db {
            DbLoc::Txn(d) => &self.t_dbs[d],
            DbLoc::X(x) => &self.cursors[x].xdb,
        }
    }

    fn db_mut(&mut self, c: usize) -> &mut MdbDb {
        match self.cursors[c].db {
            DbLoc::Txn(d) => &mut self.t_dbs[d],
            DbLoc::X(x) => &mut self.cursors[x].xdb,
        }
    }

    /// `*mc->mc_dbflag`.
    fn dbflag(&self, c: usize) -> u8 {
        match self.cursors[c].dbflag {
            DbLoc::Txn(d) => self.t_dbflags[d],
            DbLoc::X(x) => self.cursors[x].xdbflag,
        }
    }

    fn dbflag_or(&mut self, c: usize, f: u8) {
        match self.cursors[c].dbflag {
            DbLoc::Txn(d) => self.t_dbflags[d] |= f,
            DbLoc::X(x) => self.cursors[x].xdbflag |= f,
        }
    }

    fn dbflag_clear(&mut self, c: usize, f: u8) {
        match self.cursors[c].dbflag {
            DbLoc::Txn(d) => self.t_dbflags[d] &= !f,
            DbLoc::X(x) => self.cursors[x].xdbflag &= !f,
        }
    }

    fn top_page(&self, c: usize) -> PageRef {
        let k = &self.cursors[c];
        k.pg[k.top]
    }

    fn top_ki(&self, c: usize) -> usize {
        let k = &self.cursors[c];
        k.ki[k.top] as usize
    }

    fn set_top_ki(&mut self, c: usize, v: usize) {
        let k = &mut self.cursors[c];
        k.ki[k.top] = v as u16;
    }

    /// The `MDB_page *` of the data of a node (`NODEDATA` of the node taken as a sub-page).
    fn sub_ref(&self, page: PageRef, off: usize) -> PageRef {
        match page {
            PageRef::Pg(h) => PageRef::Sub(h, off as u16),
            PageRef::Sub(h, o) => PageRef::Sub(h, o + off as u16),
            _ => panic!("LMDB writer: sub-page in a scratch page"),
        }
    }

    fn slot_bytes(&self, s: Slot) -> &[u8] {
        &self.page(s.page)[s.off..s.off + s.len]
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:7596-7624 (mdb_cursor_init)
    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:7505-7522 (mdb_xcursor_init0)
    /// `mdb_cursor_init` with an xcursor for a DUPSORT database: the index of the cursor.
    pub(super) fn cursor_init(&mut self, dbi: usize) -> usize {
        let c = self.cursors.len();
        self.cursors.push(Cursor {
            dbi,
            db: DbLoc::Txn(dbi),
            dbflag: DbLoc::Txn(dbi),
            cmp: self.dbx_cmp[dbi],
            dcmp: self.dbx_dcmp[dbi],
            flags: 0,
            snum: 0,
            top: 0,
            pg: [PageRef::Null; CURSOR_STACK],
            ki: [0; CURSOR_STACK],
            xc: None,
            xdb: MdbDb::default(),
            xdbflag: 0,
        });
        if self.t_dbs[dbi].flags as u32 & MDB_DUPSORT != 0 {
            let x = c + 1;
            self.cursors.push(Cursor {
                dbi,
                db: DbLoc::X(x),
                dbflag: DbLoc::X(x),
                cmp: self.dbx_dcmp[dbi],
                dcmp: Cmp::None,
                flags: C_SUB,
                snum: 0,
                top: 0,
                pg: [PageRef::Null; CURSOR_STACK],
                ki: [0; CURSOR_STACK],
                xc: None,
                xdb: MdbDb::default(),
                xdbflag: 0,
            });
            self.cursors[c].xc = Some(x);
        }
        if self.dbflag(c) & DB_STALE != 0 {
            self.page_search(c, None, MDB_PS_ROOTONLY);
        }
        c
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:8214-8234 (mdb_cursor_copy)
    /// `mdb_cursor_copy` into a new cursor without an xcursor (`mn.mc_xcursor` is set to `NULL`).
    fn cursor_copy(&mut self, src: usize) -> usize {
        let mut k = self.cursors[src].clone();
        k.xc = None;
        let c = self.cursors.len();
        self.cursors.push(k);
        c
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:7524-7563 (mdb_xcursor_init1)
    fn xcursor_init1(&mut self, c: usize, page: PageRef, node: usize) {
        let x = self.cursors[c].xc.expect("LMDB writer: no xcursor");
        let mflags = self.db(c).flags as u32;
        let p = self.page(page);
        if node_flags(p, node) & F_SUBDATA != 0 {
            let d = node_data(p, node);
            let db = MdbDb::from_bytes(&p[d..d + MDB_DB_SIZE]);
            let k = &mut self.cursors[x];
            k.xdb = db;
            k.pg[0] = PageRef::Null;
            k.snum = 0;
            k.top = 0;
            k.flags = C_SUB;
        } else {
            let d = node_data(p, node);
            let fp = self.sub_ref(page, d);
            let fpb = self.page(fp);
            let entries = numkeys(fpb) as u64;
            let root = pgno(fpb);
            let fpad = pad(fpb);
            let k = &mut self.cursors[x];
            k.xdb = MdbDb {
                pad: 0,
                flags: 0,
                depth: 1,
                branch_pages: 0,
                leaf_pages: 1,
                overflow_pages: 0,
                entries,
                root,
            };
            k.snum = 1;
            k.top = 0;
            k.flags = C_INITIALIZED | C_SUB;
            k.pg[0] = fp;
            k.ki[0] = 0;
            if mflags & MDB_DUPFIXED != 0 {
                k.xdb.flags = MDB_DUPFIXED as u16;
                k.xdb.pad = fpad as u32;
                if mflags & MDB_INTEGERDUP != 0 {
                    k.xdb.flags |= MDB_INTEGERKEY as u16;
                }
            }
        }
        self.cursors[x].xdbflag = DB_VALID | DB_USRVALID | DB_DUPDATA;
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:7565-7594 (mdb_xcursor_init2)
    fn xcursor_init2(&mut self, m2: usize, src_x: usize, new_dupdata: bool) {
        let x = match self.cursors[m2].xc {
            Some(x) => x,
            None => return,
        };
        if new_dupdata {
            let k = &mut self.cursors[x];
            k.snum = 1;
            k.top = 0;
            k.flags |= C_INITIALIZED;
            k.ki[0] = 0;
            k.xdbflag = DB_VALID | DB_USRVALID | DB_DUPDATA;
        } else if self.cursors[x].flags & C_INITIALIZED == 0 {
            return;
        }
        let (db, pg0) = (self.cursors[src_x].xdb, self.cursors[src_x].pg[0]);
        self.cursors[x].xdb = db;
        self.cursors[x].pg[0] = pg0;
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1243-1250 (XCURSOR_REFRESH)
    fn xcursor_refresh(&mut self, m: usize, top: usize, mp: PageRef) {
        let x = match self.cursors[m].xc {
            Some(x) => x,
            None => return,
        };
        if self.cursors[x].flags & C_INITIALIZED == 0 {
            return;
        }
        let ki = self.cursors[m].ki[top] as usize;
        let p = self.page(mp);
        if ki >= numkeys(p) {
            return;
        }
        let n = nodeptr(p, ki);
        if node_flags(p, n) & (F_DUPDATA | F_SUBDATA) == F_DUPDATA {
            let d = node_data(p, n);
            let r = self.sub_ref(mp, d);
            self.cursors[x].pg[0] = r;
        }
    }

    /// Entry `k` (from the head) of `txn->mt_cursors[dbi]`, as the `m2` and `m3` of a fixup
    /// loop: `m2`'s flags, `m2` (None for the dummy of `WITH_CURSOR_TRACKING`), and
    /// `m3`, which is `mx_cursor` of `m2->mc_xcursor` when `mc->mc_flags` has `C_SUB` set, else
    /// `m2`.
    fn tracked_at(&self, dbi: usize, k: usize, sub: bool) -> Option<(u32, Option<usize>, usize)> {
        let list = &self.tracked[dbi];
        match list[list.len() - 1 - k] {
            Track::Cursor(id) => {
                let m3 = if sub { self.cursors[id].xc? } else { id };
                Some((self.cursors[id].flags, Some(id), m3))
            }
            Track::Dummy(mn) => {
                if sub {
                    Some((C_INITIALIZED, None, mn))
                } else {
                    None
                }
            }
        }
    }

    // ---- searching -------------------------------------------------------------------------

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5268-5365 (mdb_node_search)
    // The integer comparisons of branch pages (`mdb_cmp_int`/`mdb_cmp_long` for
    // `mdb_cmp_cint`) order keys of one size as `Cmp::Int` does.
    /// `mdb_node_search`: the node (offset; 0 for a LEAF2 page) and the exact flag.
    fn node_search(&mut self, c: usize, key: &[u8]) -> (Option<usize>, bool) {
        let (mp, cmp) = {
            let k = &self.cursors[c];
            (k.pg[k.top], k.cmp)
        };
        let ks = self.db(c).pad as usize;
        let p = self.page(mp);
        let nkeys = numkeys(p);
        let mut low: isize = if is_leaf(p) { 0 } else { 1 };
        let mut high: isize = nkeys as isize - 1;
        let mut i: usize = 0;
        let mut rc = 0;
        let leaf2 = is_leaf2(p);
        while low <= high {
            i = ((low + high) >> 1) as usize;
            if leaf2 {
                let o = leaf2key(i, ks);
                rc = cmp.cmp(key, &p[o..o + ks]);
            } else {
                let n = nodeptr(p, i);
                let o = node_key(n);
                rc = cmp.cmp(key, &p[o..o + node_ksize(p, n)]);
            }
            if rc == 0 {
                break;
            }
            if rc > 0 {
                low = i as isize + 1;
            } else {
                high = i as isize - 1;
            }
        }
        if rc > 0 {
            i += 1;
        }
        let exact = rc == 0 && nkeys > 0;
        let node = if i < nkeys {
            Some(if leaf2 { 0 } else { nodeptr(p, i) })
        } else {
            None
        };
        self.set_top_ki(c, i);
        (node, exact)
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5381-5417 (mdb_cursor_pop,
    // mdb_cursor_push)
    fn cursor_pop(&mut self, c: usize) {
        let k = &mut self.cursors[c];
        if k.snum > 0 {
            k.snum -= 1;
            if k.snum > 0 {
                k.top -= 1;
            } else {
                k.flags &= !C_INITIALIZED;
            }
        }
    }

    fn cursor_push(&mut self, c: usize, mp: PageRef) -> i32 {
        if self.cursors[c].snum >= CURSOR_STACK {
            self.t_flags |= MDB_TXN_ERROR;
            return MDB_CURSOR_FULL;
        }
        let k = &mut self.cursors[c];
        k.top = k.snum;
        k.snum += 1;
        k.pg[k.top] = mp;
        k.ki[k.top] = 0;
        MDB_SUCCESS
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5480-5550 (mdb_page_search_root)
    fn page_search_root(&mut self, c: usize, key: Option<&[u8]>, flags: u32) -> i32 {
        let mut mp = self.top_page(c);
        while is_branch(self.page(mp)) {
            let nk = numkeys(self.page(mp));
            let mut i: usize;
            let mut ready = false;
            if flags & (MDB_PS_FIRST | MDB_PS_LAST) != 0 {
                i = 0;
                if flags & MDB_PS_LAST != 0 {
                    i = nk - 1;
                    let k = &self.cursors[c];
                    if k.flags & C_INITIALIZED != 0 && k.ki[k.top] as usize == i {
                        let k = &mut self.cursors[c];
                        k.top = k.snum;
                        k.snum += 1;
                        mp = k.pg[k.top];
                        ready = true;
                    }
                }
            } else {
                let (node, exact) = self.node_search(c, key.unwrap());
                if node.is_none() {
                    i = nk - 1;
                } else {
                    i = self.top_ki(c);
                    if !exact {
                        i -= 1;
                    }
                }
            }
            if !ready {
                let p = self.page(mp);
                let child = node_pgno(p, nodeptr(p, i));
                let (np, _) = match self.page_get(child) {
                    Ok(r) => r,
                    Err(rc) => return rc,
                };
                self.set_top_ki(c, i);
                let rc = self.cursor_push(c, np);
                if rc != MDB_SUCCESS {
                    return rc;
                }
                mp = np;
            }
            if flags & MDB_PS_MODIFY != 0 {
                let rc = self.page_touch(c);
                if rc != MDB_SUCCESS {
                    return rc;
                }
                mp = self.top_page(c);
            }
        }
        if !is_leaf(self.page(mp)) {
            self.t_flags |= MDB_TXN_ERROR;
            return MDB_CORRUPTED;
        }
        let k = &mut self.cursors[c];
        k.flags |= C_INITIALIZED;
        k.flags &= !C_EOF;
        MDB_SUCCESS
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5578-5671 (mdb_page_search)
    pub(super) fn page_search(&mut self, c: usize, key: Option<&[u8]>, flags: u32) -> i32 {
        if self.t_flags & MDB_TXN_BLOCKED != 0 {
            return MDB_BAD_TXN;
        }
        if self.dbflag(c) & DB_STALE != 0 {
            let name = self.dbx_name[self.cursors[c].dbi].clone();
            let mc2 = self.cursor_init(MAIN_DBI);
            let rc = self.page_search(mc2, Some(&name), 0);
            if rc != MDB_SUCCESS {
                return rc;
            }
            let (leaf, exact) = self.node_search(mc2, &name);
            if !exact {
                return MDB_NOTFOUND;
            }
            let leaf = leaf.unwrap();
            let lp = self.top_page(mc2);
            if node_flags(self.page(lp), leaf) & (F_DUPDATA | F_SUBDATA) != F_SUBDATA {
                return MDB_INCOMPATIBLE;
            }
            let data = match self.node_read(lp, leaf) {
                Ok(d) => d,
                Err(rc) => return rc,
            };
            let db = MdbDb::from_bytes(self.slot_bytes(data));
            self.cursors.truncate(mc2);
            if self.db(c).flags & PERSISTENT_FLAGS != db.flags {
                return MDB_INCOMPATIBLE;
            }
            *self.db_mut(c) = db;
            self.dbflag_clear(c, DB_STALE);
        }
        let root = self.db(c).root;
        if root == P_INVALID {
            return MDB_NOTFOUND;
        }
        let pg0 = self.cursors[c].pg[0];
        if pg0 == PageRef::Null || pgno(self.page(pg0)) != root {
            match self.page_get(root) {
                Ok((r, _)) => self.cursors[c].pg[0] = r,
                Err(rc) => return rc,
            }
        }
        self.cursors[c].snum = 1;
        self.cursors[c].top = 0;
        if flags & MDB_PS_MODIFY != 0 {
            let rc = self.page_touch(c);
            if rc != MDB_SUCCESS {
                return rc;
            }
        }
        if flags & MDB_PS_ROOTONLY != 0 {
            return MDB_SUCCESS;
        }
        self.page_search_root(c, key, flags)
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5749-5779 (mdb_node_read)
    pub(super) fn node_read(&mut self, page: PageRef, leaf: usize) -> Result<Slot, i32> {
        let p = self.page(page);
        let len = node_dsz(p, leaf);
        let d = node_data(p, leaf);
        if node_flags(p, leaf) & F_BIGDATA == 0 {
            return Ok(Slot { page, off: d, len });
        }
        let pg = get_u64(p, d);
        let (omp, _) = self.page_get(pg)?;
        Ok(Slot {
            page: omp,
            off: PAGEHDRSZ,
            len,
        })
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5805-5858 (mdb_cursor_sibling)
    fn cursor_sibling(&mut self, c: usize, move_right: bool) -> i32 {
        if self.cursors[c].snum < 2 {
            return MDB_NOTFOUND;
        }
        self.cursor_pop(c);
        let ki = self.top_ki(c);
        let at_end = if move_right {
            ki + 1 >= numkeys(self.page(self.top_page(c)))
        } else {
            ki == 0
        };
        if at_end {
            let rc = self.cursor_sibling(c, move_right);
            if rc != MDB_SUCCESS {
                let k = &mut self.cursors[c];
                k.top += 1;
                k.snum += 1;
                return rc;
            }
        } else if move_right {
            self.set_top_ki(c, ki + 1);
        } else {
            self.set_top_ki(c, ki - 1);
        }
        let tp = self.top_page(c);
        let p = self.page(tp);
        let child = node_pgno(p, nodeptr(p, self.top_ki(c)));
        let mp = match self.page_get(child) {
            Ok((r, _)) => r,
            Err(rc) => {
                self.cursors[c].flags &= !(C_INITIALIZED | C_EOF);
                return rc;
            }
        };
        self.cursor_push(c, mp);
        if !move_right {
            let n = numkeys(self.page(mp));
            self.set_top_ki(c, n.wrapping_sub(1));
        }
        MDB_SUCCESS
    }

    /// The key at the cursor (`MDB_GET_KEY`, or `LEAF2KEY` on a LEAF2 page).
    fn key_slot(&self, c: usize, leaf: usize) -> Slot {
        let page = self.top_page(c);
        let p = self.page(page);
        if is_leaf2(p) {
            let ks = self.db(c).pad as usize;
            return Slot {
                page,
                off: leaf2key(self.top_ki(c), ks),
                len: ks,
            };
        }
        Slot {
            page,
            off: node_key(leaf),
            len: node_ksize(p, leaf),
        }
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5860-5939 (mdb_cursor_next)
    /// `mdb_cursor_next` of `mc`, `key` and `data` with `MDB_NEXT`, on a database without dups.
    pub(super) fn cursor_next(
        &mut self,
        c: usize,
        want_data: bool,
    ) -> (i32, Option<Slot>, Option<Slot>) {
        if self.cursors[c].flags & C_INITIALIZED == 0 {
            return self.cursor_first(c, want_data);
        }
        let mut mp = self.top_page(c);
        if self.cursors[c].flags & C_EOF != 0 {
            if self.top_ki(c) >= numkeys(self.page(mp)).wrapping_sub(1) {
                return (MDB_NOTFOUND, None, None);
            }
            self.cursors[c].flags ^= C_EOF;
        }
        if self.db(c).flags as u32 & MDB_DUPSORT != 0 {
            return (
                self.unreachable("mdb_cursor_next on a DUPSORT database"),
                None,
                None,
            );
        }
        if self.cursors[c].flags & C_DEL != 0 {
            self.cursors[c].flags ^= C_DEL;
        } else if self.top_ki(c) + 1 >= numkeys(self.page(mp)) {
            let rc = self.cursor_sibling(c, true);
            if rc != MDB_SUCCESS {
                self.cursors[c].flags |= C_EOF;
                return (rc, None, None);
            }
            mp = self.top_page(c);
        } else {
            let ki = self.top_ki(c);
            self.set_top_ki(c, ki + 1);
        }
        if is_leaf2(self.page(mp)) {
            return (MDB_SUCCESS, Some(self.key_slot(c, 0)), None);
        }
        let leaf = nodeptr(self.page(mp), self.top_ki(c));
        if node_flags(self.page(mp), leaf) & F_DUPDATA != 0 {
            self.xcursor_init1(c, mp, leaf);
        }
        let mut data = None;
        if want_data {
            match self.node_read(mp, leaf) {
                Ok(d) => data = Some(d),
                Err(rc) => return (rc, None, None),
            }
        }
        (MDB_SUCCESS, Some(self.key_slot(c, leaf)), data)
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6030-6216 (mdb_cursor_set)
    /// `mdb_cursor_set`: the return code, the key found (`MDB_SET_RANGE`) and the data.
    pub(super) fn cursor_set(
        &mut self,
        c: usize,
        key: &[u8],
        want_data: bool,
        op: CursorOp,
        mut exactp: Option<&mut bool>,
    ) -> (i32, Option<Slot>, Option<Slot>) {
        if key.is_empty() {
            return (MDB_BAD_VALSIZE, None, None);
        }
        if let Some(x) = self.cursors[c].xc {
            self.cursors[x].flags &= !(C_INITIALIZED | C_EOF);
        }
        let cmp = self.cursors[c].cmp;
        let ks = self.db(c).pad as usize;
        let mut leaf: Option<usize> = None;
        let mut mp;
        let mut rc;
        // 0: page_search, 1: set2, 2: set1
        let mut stage = 0;
        if self.cursors[c].flags & C_INITIALIZED != 0 {
            mp = self.top_page(c);
            let p = self.page(mp);
            let nkeys = numkeys(p);
            if nkeys == 0 {
                self.set_top_ki(c, 0);
                return (MDB_NOTFOUND, None, None);
            }
            let leaf2 = is_leaf2(p);
            let key_at = |p: &[u8], i: usize| -> (usize, usize, usize) {
                if leaf2 {
                    (0, leaf2key(i, ks), ks)
                } else {
                    let n = nodeptr(p, i);
                    (n, node_key(n), node_ksize(p, n))
                }
            };
            let (n0, o, l) = key_at(p, 0);
            if !leaf2 {
                leaf = Some(n0);
            }
            rc = cmp.cmp(key, &p[o..o + l]);
            if rc == 0 {
                self.set_top_ki(c, 0);
                if let Some(e) = exactp.as_deref_mut() {
                    *e = true;
                }
                stage = 2;
            } else {
                if rc > 0 {
                    if nkeys > 1 {
                        let p = self.page(mp);
                        let (n1, o, l) = key_at(p, nkeys - 1);
                        if !leaf2 {
                            leaf = Some(n1);
                        }
                        rc = cmp.cmp(key, &p[o..o + l]);
                        if rc == 0 {
                            self.set_top_ki(c, nkeys - 1);
                            if let Some(e) = exactp.as_deref_mut() {
                                *e = true;
                            }
                            stage = 2;
                        } else if rc < 0 {
                            let ki = self.top_ki(c);
                            if ki < nkeys {
                                let p = self.page(mp);
                                let (nk, o, l) = key_at(p, ki);
                                if !leaf2 {
                                    leaf = Some(nk);
                                }
                                rc = cmp.cmp(key, &p[o..o + l]);
                                if rc == 0 {
                                    if let Some(e) = exactp.as_deref_mut() {
                                        *e = true;
                                    }
                                    stage = 2;
                                }
                            }
                            if stage != 2 {
                                self.cursors[c].flags &= !C_EOF;
                                stage = 1;
                            }
                        }
                    }
                    if stage == 0 {
                        let top = self.cursors[c].top;
                        let mut i = 0;
                        while i < top {
                            let k = &self.cursors[c];
                            if (k.ki[i] as usize) < numkeys(self.page(k.pg[i])).wrapping_sub(1) {
                                break;
                            }
                            i += 1;
                        }
                        if i == top {
                            self.set_top_ki(c, nkeys);
                            return (MDB_NOTFOUND, None, None);
                        }
                    }
                }
                if stage == 0 && self.cursors[c].top == 0 {
                    self.set_top_ki(c, 0);
                    if op == CursorOp::SetRange && exactp.is_none() {
                        stage = 2;
                    } else {
                        return (MDB_NOTFOUND, None, None);
                    }
                }
            }
        } else {
            self.cursors[c].pg[0] = PageRef::Null;
            mp = PageRef::Null;
        }
        if stage == 0 {
            rc = self.page_search(c, Some(key), 0);
            if rc != MDB_SUCCESS {
                return (rc, None, None);
            }
            mp = self.top_page(c);
            stage = 1;
        }
        if stage == 1 {
            let (l, ex) = self.node_search(c, key);
            leaf = l;
            if let Some(e) = exactp {
                *e = ex;
                if !ex {
                    return (MDB_NOTFOUND, None, None);
                }
            }
            if leaf.is_none() {
                let rc = self.cursor_sibling(c, true);
                if rc != MDB_SUCCESS {
                    self.cursors[c].flags |= C_EOF;
                    return (rc, None, None);
                }
                mp = self.top_page(c);
                leaf = Some(nodeptr(self.page(mp), 0));
            }
        }
        // set1:
        self.cursors[c].flags |= C_INITIALIZED;
        self.cursors[c].flags &= !C_EOF;
        let mut key_out = None;
        if is_leaf2(self.page(mp)) {
            if op == CursorOp::SetRange {
                key_out = Some(self.key_slot(c, 0));
            }
            return (MDB_SUCCESS, key_out, None);
        }
        let leaf = leaf.unwrap();
        let dup = node_flags(self.page(mp), leaf) & F_DUPDATA != 0;
        if dup {
            self.xcursor_init1(c, mp, leaf);
        }
        let mut rc = MDB_SUCCESS;
        let mut data = None;
        if want_data {
            if dup {
                let x = self.cursors[c].xc.unwrap();
                let (r, k, _) = self.cursor_first(x, false);
                rc = r;
                data = k;
            } else {
                if let Some(x) = self.cursors[c].xc {
                    self.cursors[x].flags &= !(C_INITIALIZED | C_EOF);
                }
                match self.node_read(mp, leaf) {
                    Ok(d) => data = Some(d),
                    Err(r) => return (r, None, None),
                }
            }
        }
        if op == CursorOp::SetRange {
            key_out = Some(self.key_slot(c, leaf));
        }
        (rc, key_out, data)
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6218-6256 (mdb_cursor_first)
    pub(super) fn cursor_first(
        &mut self,
        c: usize,
        want_data: bool,
    ) -> (i32, Option<Slot>, Option<Slot>) {
        self.cursor_end(c, want_data, false)
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6258-6296 (mdb_cursor_last)
    fn cursor_last(&mut self, c: usize, want_data: bool) -> (i32, Option<Slot>, Option<Slot>) {
        self.cursor_end(c, want_data, true)
    }

    fn cursor_end(
        &mut self,
        c: usize,
        want_data: bool,
        last: bool,
    ) -> (i32, Option<Slot>, Option<Slot>) {
        if let Some(x) = self.cursors[c].xc {
            self.cursors[x].flags &= !(C_INITIALIZED | C_EOF);
        }
        if self.cursors[c].flags & C_INITIALIZED == 0 || self.cursors[c].top != 0 {
            let rc = self.page_search(c, None, if last { MDB_PS_LAST } else { MDB_PS_FIRST });
            if rc != MDB_SUCCESS {
                return (rc, None, None);
            }
        }
        let mp = self.top_page(c);
        if last {
            let n = numkeys(self.page(mp));
            self.set_top_ki(c, n.wrapping_sub(1));
            self.cursors[c].flags |= C_INITIALIZED | C_EOF;
        } else {
            self.cursors[c].flags |= C_INITIALIZED;
            self.cursors[c].flags &= !C_EOF;
            self.set_top_ki(c, 0);
        }
        if is_leaf2(self.page(mp)) {
            return (MDB_SUCCESS, Some(self.key_slot(c, 0)), None);
        }
        let leaf = nodeptr(self.page(mp), self.top_ki(c));
        let mut data = None;
        if want_data {
            if node_flags(self.page(mp), leaf) & F_DUPDATA != 0 {
                self.xcursor_init1(c, mp, leaf);
                let x = self.cursors[c].xc.unwrap();
                let (rc, k, _) = self.cursor_end(x, false, last);
                if rc != MDB_SUCCESS {
                    return (rc, None, None);
                }
                data = k;
            } else {
                match self.node_read(mp, leaf) {
                    Ok(d) => data = Some(d),
                    Err(rc) => return (rc, None, None),
                }
            }
        }
        (MDB_SUCCESS, Some(self.key_slot(c, leaf)), data)
    }

    // ---- touching pages --------------------------------------------------------------------

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6505-6534 (mdb_cursor_touch)
    fn cursor_touch(&mut self, c: usize) -> i32 {
        let mut rc = MDB_SUCCESS;
        let dbi = self.cursors[c].dbi;
        if dbi >= CORE_DBS && self.dbflag(c) & (DB_DIRTY | DB_DUPDATA) == 0 {
            let name = self.dbx_name[dbi].clone();
            let mc2 = self.cursor_init(MAIN_DBI);
            rc = self.page_search(mc2, Some(&name), MDB_PS_MODIFY);
            self.cursors.truncate(mc2);
            if rc != MDB_SUCCESS {
                return rc;
            }
            self.dbflag_or(c, DB_DIRTY);
        }
        self.cursors[c].top = 0;
        if self.cursors[c].snum > 0 {
            loop {
                rc = self.page_touch(c);
                if rc != MDB_SUCCESS {
                    break;
                }
                self.cursors[c].top += 1;
                if self.cursors[c].top >= self.cursors[c].snum {
                    break;
                }
            }
            self.cursors[c].top = self.cursors[c].snum - 1;
        }
        rc
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:2410-2510 (mdb_page_touch)
    fn page_touch(&mut self, c: usize) -> i32 {
        let (mp, top) = {
            let k = &self.cursors[c];
            (k.pg[k.top], k.top)
        };
        if flags(self.page(mp)) & P_DIRTY != 0 {
            return MDB_SUCCESS;
        }
        let new = match self.page_alloc(1) {
            Ok(p) => p,
            Err(rc) => {
                self.t_flags |= MDB_TXN_ERROR;
                return rc;
            }
        };
        let old = pgno(self.page(mp));
        midl_xappend(&mut self.t_free_pgs, old);
        if top > 0 {
            let (parent, pki) = {
                let k = &self.cursors[c];
                (k.pg[top - 1], k.ki[top - 1] as usize)
            };
            let n = nodeptr(self.page(parent), pki);
            set_node_pgno(self.page_mut(parent), n, new);
        } else {
            self.db_mut(c).root = new;
        }
        {
            let src = self
                .clean
                .get(&old)
                .expect("LMDB writer: clean page not read");
            let i = self.t_dirty.binary_search_by(|e| e.0.cmp(&new)).unwrap();
            Self::page_copy(&mut self.t_dirty[i].1, src, self.psize);
            let np = &mut self.t_dirty[i].1;
            set_pgno(np, new);
            let f = flags(np) | P_DIRTY;
            set_flags(np, f);
        }
        let np = PageRef::Pg(new);
        self.cursors[c].pg[top] = np;
        let dbi = self.cursors[c].dbi;
        let sub = self.cursors[c].flags & C_SUB != 0;
        let snum = self.cursors[c].snum;
        for k in 0..self.tracked[dbi].len() {
            let Some((_, m2, m3)) = self.tracked_at(dbi, k, sub) else {
                continue;
            };
            if sub {
                if self.cursors[m3].snum < snum {
                    continue;
                }
                if self.cursors[m3].pg[top] == mp {
                    self.cursors[m3].pg[top] = np;
                }
            } else {
                let m2 = m2.unwrap();
                if self.cursors[m2].snum < snum || m2 == c {
                    continue;
                }
                if self.cursors[m2].pg[top] == mp {
                    self.cursors[m2].pg[top] = np;
                    if is_leaf(self.page(np)) {
                        self.xcursor_refresh(m2, top, np);
                    }
                }
            }
        }
        MDB_SUCCESS
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:1997-2024 (mdb_page_spill)
    // A transaction of makeblastdb dirties far fewer than 131071 pages; spilling is not
    // ported.
    fn page_spill(&mut self, c: usize, kd: Option<(usize, usize)>) -> i32 {
        if self.cursors[c].flags & C_SUB != 0 {
            return MDB_SUCCESS;
        }
        let mut i = self.db(c).depth as usize;
        if self.cursors[c].dbi >= CORE_DBS {
            i += self.t_dbs[MAIN_DBI].depth as usize;
        }
        if let Some((k, d)) = kd {
            i += (NODESIZE + k + d + self.psize) / self.psize;
        }
        i += i;
        if self.t_dirty_room > i {
            return MDB_SUCCESS;
        }
        self.unreachable("mdb_page_spill (dirty list nearly full)")
    }

    // ---- put -------------------------------------------------------------------------------

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6540-7066 (mdb_cursor_put); the
    // parts are cited where they are ported below.
    /// `mdb_cursor_put`: the return code and, for `MDB_RESERVE`, where the data goes.
    pub(super) fn cursor_put(
        &mut self,
        mc: usize,
        key: &[u8],
        data: Val,
        flags: u32,
    ) -> (i32, Option<Loc>) {
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6571-6581
        if flags & MDB_MULTIPLE != 0 {
            return (self.unreachable("MDB_MULTIPLE"), None);
        }
        let nospill = flags & MDB_NOSPILL;
        let mut flags = flags & !MDB_NOSPILL;
        if self.t_flags & MDB_TXN_BLOCKED != 0 {
            return (MDB_BAD_TXN, None);
        }
        if key.len().wrapping_sub(1) >= MDB_MAXKEYSIZE {
            return (MDB_BAD_VALSIZE, None);
        }
        let dbflags = self.db(mc).flags as u32;
        let dlimit = if dbflags & MDB_DUPSORT != 0 {
            MDB_MAXKEYSIZE
        } else {
            MAXDATASIZE
        };
        if data.size() > dlimit {
            return (MDB_BAD_VALSIZE, None);
        }
        let mut dkey: Vec<u8> = Vec::new();
        let rc: i32;
        let mut no_root = false;
        if flags == MDB_CURRENT {
            if self.cursors[mc].flags & C_INITIALIZED == 0 {
                return (EINVAL, None);
            }
            rc = MDB_SUCCESS;
        } else if self.db(mc).root == P_INVALID {
            let k = &mut self.cursors[mc];
            k.snum = 0;
            k.top = 0;
            k.flags &= !C_INITIALIZED;
            rc = MDB_NOTFOUND;
            no_root = true;
        } else {
            let mut exact = false;
            let mut r;
            if flags & MDB_APPEND != 0 {
                let (r0, k2, _) = self.cursor_last(mc, true);
                r = r0;
                if r == MDB_SUCCESS {
                    let cmp = self.cursors[mc].cmp;
                    r = cmp.cmp(key, self.slot_bytes(k2.unwrap()));
                    if r > 0 {
                        r = MDB_NOTFOUND;
                        let ki = self.top_ki(mc);
                        self.set_top_ki(mc, ki + 1);
                    } else {
                        r = MDB_KEYEXIST;
                    }
                }
            } else {
                let (r0, _, _) = self.cursor_set(mc, key, true, CursorOp::Set, Some(&mut exact));
                r = r0;
            }
            if flags & MDB_NOOVERWRITE != 0 && r == MDB_SUCCESS {
                return (MDB_KEYEXIST, None);
            }
            if r != MDB_SUCCESS && r != MDB_NOTFOUND {
                return (r, None);
            }
            rc = r;
        }
        if self.cursors[mc].flags & C_DEL != 0 {
            self.cursors[mc].flags ^= C_DEL;
        }
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6634-6668
        if nospill == 0 {
            let rc2 = self.page_spill(mc, Some((key.len(), data.size())));
            if rc2 != MDB_SUCCESS {
                return (rc2, None);
            }
        }
        if no_root {
            let np = match self.page_new(mc, P_LEAF, 1) {
                Ok(p) => p,
                Err(rc2) => return (rc2, None),
            };
            self.cursor_push(mc, PageRef::Pg(np));
            let db = self.db_mut(mc);
            db.root = np;
            db.depth += 1;
            self.dbflag_or(mc, DB_DIRTY);
            if self.db(mc).flags as u32 & (MDB_DUPSORT | MDB_DUPFIXED) == MDB_DUPFIXED {
                let p = self.page_mut(PageRef::Pg(np));
                let f = flags_of(p) | P_LEAF2;
                set_flags(p, f);
            }
            self.cursors[mc].flags |= C_INITIALIZED;
        } else {
            let rc2 = self.cursor_touch(mc);
            if rc2 != MDB_SUCCESS {
                return (rc2, None);
            }
        }
        let insert_key = rc != MDB_SUCCESS;
        let mut insert_data = rc != MDB_SUCCESS;
        let mut rc = rc;

        let mut do_sub = false;
        let mut sub_root: Option<u64> = None;
        // The node data of a new sub-page or sub-DB record (`xdata`), when not `data`.
        let mut xdata: Option<Vec<u8>> = None;
        let mut put_sub = false;
        let mut reserve_loc: Option<Loc> = None;
        if insert_key {
            // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6670-6683
            // A key of at most 511 bytes and a dup of at most 511 bytes always fit a node of a
            // page of 4096 bytes or more.
            if dbflags & MDB_DUPSORT != 0 && NODESIZE + key.len() + data.size() > self.nodemax {
                return (self.unreachable("new DUPSORT key too big for a node"), None);
            }
        } else {
            let tp = self.top_page(mc);
            if is_leaf2(self.page(tp)) {
                return (self.unreachable("overwrite of an existing LEAF2 key"), None);
            }
            let leaf = nodeptr(self.page(tp), self.top_ki(mc));
            let (leaf_flags, old_off, mut old_size) = {
                let p = self.page(tp);
                (node_flags(p, leaf), node_data(p, leaf), node_dsz(p, leaf))
            };
            let mut go_current = true;
            if dbflags & MDB_DUPSORT != 0 {
                // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6722-6731
                go_current = false;
                let host = pgno(self.page(tp));
                set_pgno(&mut self.pbuf, host);
                // fp, mp: None = me_pbuf.
                let mut fp: Option<PageRef> = None;
                let mut mp: Option<PageRef> = None;
                let mut offset: usize = 0;
                let mut xsize: usize = 0;
                let mut stage_new_sub = true;
                if leaf_flags & F_DUPDATA == 0 {
                    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6733-6769
                    // One `MDB_page` past `fp` is `me_pbuf` plus the `sizeof` of `MDB_page`, 24
                    // bytes on 64-bit hosts.
                    let dbytes = match data {
                        Val::Bytes(b) => b,
                        Val::Reserve(_) => {
                            return (self.unreachable("MDB_RESERVE into DUPSORT"), None)
                        }
                    };
                    let dcmp = self.cursors[mc].dcmp;
                    let same = dcmp.cmp(dbytes, &self.page(tp)[old_off..old_off + old_size]) == 0;
                    if flags == MDB_CURRENT || same {
                        if flags != MDB_CURRENT && flags & (MDB_NODUPDATA | MDB_APPENDDUP) != 0 {
                            return (MDB_KEYEXIST, None);
                        }
                        go_current = true;
                        stage_new_sub = false;
                    } else {
                        dkey = self.page(tp)[old_off..old_off + old_size].to_vec();
                        self.pbuf[24..24 + dkey.len()].copy_from_slice(&dkey);
                        set_flags(&mut self.pbuf, P_LEAF | P_DIRTY | P_SUBP);
                        set_lower(&mut self.pbuf, PAGEHDRSZ as u16);
                        xsize = PAGEHDRSZ + dkey.len() + dbytes.len();
                        if dbflags & MDB_DUPFIXED != 0 {
                            let f = flags_of(&self.pbuf) | P_LEAF2;
                            set_flags(&mut self.pbuf, f);
                            set_pad(&mut self.pbuf, dbytes.len() as u16);
                            xsize += 2 * dbytes.len();
                        } else {
                            return (
                                self.unreachable("sub-page of a non-DUPFIXED database"),
                                None,
                            );
                        }
                        set_upper(&mut self.pbuf, xsize as u16);
                        old_size = xsize;
                    }
                } else if leaf_flags & F_SUBDATA != 0 {
                    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6770-6771
                    flags |= (F_DUPDATA | F_SUBDATA) as u32;
                    put_sub = true;
                    stage_new_sub = false;
                } else {
                    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6775-6797
                    let f = self.sub_ref(tp, old_off);
                    fp = Some(f);
                    let mut fallthrough = flags == MDB_CURRENT;
                    if !fallthrough {
                        if dbflags & MDB_DUPFIXED == 0 {
                            return (
                                self.unreachable("sub-page of a non-DUPFIXED database"),
                                None,
                            );
                        }
                        let fpb = self.page(f);
                        offset = pad(fpb) as usize;
                        if (sizeleft(fpb) as usize) < offset {
                            offset *= 4;
                        } else {
                            fallthrough = true;
                        }
                    }
                    if fallthrough {
                        let fpb = self.page_mut(f);
                        let fl = flags_of(fpb) | P_DIRTY;
                        set_flags(fpb, fl);
                        set_pgno(fpb, host);
                        let x = self.cursors[mc].xc.unwrap();
                        self.cursors[x].pg[0] = f;
                        flags |= F_DUPDATA as u32;
                        put_sub = true;
                        stage_new_sub = false;
                    } else {
                        xsize = old_size + offset;
                    }
                }
                if stage_new_sub {
                    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6800-6850
                    let fpb: &[u8] = match fp {
                        None => &self.pbuf,
                        Some(f) => self.page(f),
                    };
                    let mut fp_flags = flags_of(fpb);
                    let fp_pad = pad(fpb);
                    let fp_lower = lower(fpb);
                    let fp_upper = upper(fpb);
                    let fp_nkeys = numkeys(fpb);
                    let ksz = node_ksize(self.page(tp), leaf);
                    let mut dummy = MdbDb::default();
                    let to_subdb = NODESIZE + ksz + xsize > self.nodemax;
                    if to_subdb {
                        fp_flags &= !P_SUBP;
                        if dbflags & MDB_DUPFIXED != 0 {
                            fp_flags |= P_LEAF2;
                            dummy.pad = fp_pad as u32;
                            dummy.flags = MDB_DUPFIXED as u16;
                            if dbflags & MDB_INTEGERDUP != 0 {
                                dummy.flags |= MDB_INTEGERKEY as u16;
                            }
                        }
                        dummy.depth = 1;
                        dummy.branch_pages = 0;
                        dummy.leaf_pages = 1;
                        dummy.overflow_pages = 0;
                        dummy.entries = fp_nkeys as u64;
                        let np = match self.page_alloc(1) {
                            Ok(p) => p,
                            Err(rc) => return (rc, None),
                        };
                        offset = self.psize - old_size;
                        flags |= (F_DUPDATA | F_SUBDATA) as u32;
                        dummy.root = np;
                        sub_root = Some(np);
                        mp = Some(PageRef::Pg(np));
                    }
                    if mp != fp {
                        let keys: Vec<u8> = {
                            let fpb: &[u8] = match fp {
                                None => &self.pbuf,
                                Some(f) => self.page(f),
                            };
                            if fp_flags & P_LEAF2 != 0 {
                                fpb[PAGEHDRSZ..PAGEHDRSZ + fp_nkeys * fp_pad as usize].to_vec()
                            } else {
                                return (
                                    self.unreachable("sub-page of a non-DUPFIXED database"),
                                    None,
                                );
                            }
                        };
                        let mpb: &mut [u8] = match mp {
                            None => &mut self.pbuf,
                            Some(m) => self.page_mut(m),
                        };
                        set_flags(mpb, fp_flags | P_DIRTY);
                        set_pad(mpb, fp_pad);
                        set_lower(mpb, fp_lower);
                        set_upper(mpb, fp_upper.wrapping_add(offset as u16));
                        mpb[PAGEHDRSZ..PAGEHDRSZ + keys.len()].copy_from_slice(&keys);
                    }
                    xdata = Some(if to_subdb {
                        dummy.to_bytes().to_vec()
                    } else {
                        self.pbuf[..xsize].to_vec()
                    });
                    flags |= F_DUPDATA as u32;
                    do_sub = true;
                    if !insert_key {
                        self.node_del(mc, 0);
                    }
                }
            }
            if go_current {
                // A dirty overflow page is in this transaction's dirty list (`level` 1), a clean
                // one in the file (`level` 0): the unspill branch is never taken.
                if (leaf_flags ^ flags as u16) & F_SUBDATA != 0 {
                    return (MDB_INCOMPATIBLE, None);
                }
                if leaf_flags & F_BIGDATA != 0 {
                    let pg = get_u64(self.page(tp), old_off);
                    let dpages = (PAGEHDRSZ - 1 + data.size()) / self.psize + 1;
                    let (omp, _) = match self.page_get(pg) {
                        Ok(r) => r,
                        Err(rc2) => return (rc2, None),
                    };
                    let ovpages = pages(self.page(omp)) as usize;
                    if ovpages >= dpages && flags_of(self.page(omp)) & P_DIRTY != 0 {
                        set_node_dsz(self.page_mut(tp), leaf, data.size());
                        return match data {
                            Val::Reserve(_) => (
                                MDB_SUCCESS,
                                Some(Loc {
                                    page: omp,
                                    off: PAGEHDRSZ,
                                }),
                            ),
                            Val::Bytes(b) => {
                                self.page_mut(omp)[PAGEHDRSZ..PAGEHDRSZ + b.len()]
                                    .copy_from_slice(b);
                                (MDB_SUCCESS, None)
                            }
                        };
                    }
                    let rc2 = self.ovpage_free(mc, pg);
                    if rc2 != MDB_SUCCESS {
                        return (rc2, None);
                    }
                } else if data.size() == old_size {
                    if flags & MDB_RESERVE != 0 {
                        return (
                            MDB_SUCCESS,
                            Some(Loc {
                                page: tp,
                                off: old_off,
                            }),
                        );
                    } else if self.cursors[mc].flags & C_SUB == 0 {
                        if let Val::Bytes(b) = data {
                            self.page_mut(tp)[old_off..old_off + b.len()].copy_from_slice(b);
                        }
                    } else {
                        return (
                            self.unreachable("same-size overwrite in a sub-cursor"),
                            None,
                        );
                    }
                    return (MDB_SUCCESS, None);
                }
                self.node_del(mc, 0);
            }
        }
        if !put_sub {
            // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6942-6970
            let nflags = flags & (F_DUPDATA as u32 | F_SUBDATA as u32 | MDB_RESERVE | MDB_APPEND);
            let rdata: Val = match &xdata {
                Some(x) => Val::Bytes(x),
                None => data,
            };
            let tp = self.top_page(mc);
            let nsize = if is_leaf2(self.page(tp)) {
                key.len()
            } else {
                self.leaf_size(key.len(), rdata.size())
            };
            if (sizeleft(self.page(tp)) as usize) < nsize {
                let mut nflags = nflags;
                if flags & (F_DUPDATA | F_SUBDATA) as u32 == F_DUPDATA as u32 {
                    nflags &= !MDB_APPEND;
                }
                if !insert_key {
                    nflags |= MDB_SPLIT_REPLACE;
                }
                let (r, loc) = self.page_split(mc, key, Some(rdata), P_INVALID, nflags);
                rc = r;
                reserve_loc = loc;
            } else {
                let ki = self.top_ki(mc);
                let (dsize, dbytes): (usize, &[u8]) = match rdata {
                    Val::Bytes(b) => (b.len(), b),
                    Val::Reserve(n) => (n, &[]),
                };
                let (r, loc) = self.node_add(mc, ki, Some(key), dsize, dbytes, 0, nflags);
                rc = r;
                reserve_loc = loc;
                if rc == MDB_SUCCESS {
                    let dbi = self.cursors[mc].dbi;
                    let sub = self.cursors[mc].flags & C_SUB != 0;
                    let i = self.cursors[mc].top;
                    let mp = self.cursors[mc].pg[i];
                    let (snum, mki) = (self.cursors[mc].snum, self.cursors[mc].ki[i]);
                    for k in 0..self.tracked[dbi].len() {
                        let Some((_, _, m3)) = self.tracked_at(dbi, k, sub) else {
                            continue;
                        };
                        let c3 = &self.cursors[m3];
                        if m3 == mc || c3.snum < snum || c3.pg[i] != mp {
                            continue;
                        }
                        if c3.ki[i] >= mki && insert_key {
                            self.cursors[m3].ki[i] += 1;
                        }
                        self.xcursor_refresh(m3, i, mp);
                    }
                }
            }
            if rc != MDB_SUCCESS {
                self.t_flags |= MDB_TXN_ERROR;
                return (rc, None);
            }
        }
        if do_sub || put_sub {
            // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:6984-7033
            let dbytes = match data {
                Val::Bytes(b) => b,
                Val::Reserve(_) => return (self.unreachable("MDB_RESERVE into DUPSORT"), None),
            };
            let x = self.cursors[mc].xc.unwrap();
            let tp = self.top_page(mc);
            let leaf = nodeptr(self.page(tp), self.top_ki(mc));
            let mut xflags;
            if flags & MDB_CURRENT != 0 {
                xflags = MDB_CURRENT | MDB_NOSPILL;
            } else {
                self.xcursor_init1(mc, tp, leaf);
                xflags = if flags & MDB_NODUPDATA != 0 {
                    MDB_NOOVERWRITE | MDB_NOSPILL
                } else {
                    MDB_NOSPILL
                };
            }
            if let Some(r) = sub_root {
                self.cursors[x].pg[0] = PageRef::Pg(r);
            }
            let new_dupdata = !dkey.is_empty();
            if !dkey.is_empty() {
                let dk = std::mem::take(&mut dkey);
                let (r, _) = self.cursor_put(x, &dk, Val::Bytes(&[]), xflags);
                rc = r;
                if rc != MDB_SUCCESS {
                    return self.bad_sub(rc);
                }
            }
            let leaf_subdata = node_flags(self.page(tp), leaf) & F_SUBDATA != 0;
            if !leaf_subdata || sub_root.is_some() {
                let dbi = self.cursors[mc].dbi;
                let i = self.cursors[mc].top;
                let mp = self.cursors[mc].pg[i];
                let (snum, mki) = (self.cursors[mc].snum, self.cursors[mc].ki[i]);
                for k in 0..self.tracked[dbi].len() {
                    let Some((_, Some(m2), _)) = self.tracked_at(dbi, k, false) else {
                        continue;
                    };
                    let c2 = &self.cursors[m2];
                    if m2 == mc || c2.snum < snum || c2.flags & C_INITIALIZED == 0 {
                        continue;
                    }
                    if c2.pg[i] == mp {
                        if c2.ki[i] == mki {
                            self.xcursor_init2(m2, x, new_dupdata);
                        } else if !insert_key {
                            self.xcursor_refresh(m2, i, mp);
                        }
                    }
                }
            }
            let ecount = self.cursors[x].xdb.entries;
            if flags & MDB_APPENDDUP != 0 {
                xflags |= MDB_APPEND;
            }
            let (r, _) = self.cursor_put(x, dbytes, Val::Bytes(&[]), xflags);
            rc = r;
            if flags & F_SUBDATA as u32 != 0 {
                let db = self.cursors[x].xdb.to_bytes();
                let d = node_data(self.page(tp), leaf);
                self.page_mut(tp)[d..d + MDB_DB_SIZE].copy_from_slice(&db);
            }
            insert_data = self.cursors[x].xdb.entries != ecount;
        }
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:7035-7062
        if insert_data {
            self.db_mut(mc).entries += 1;
        }
        if insert_key {
            if rc != MDB_SUCCESS {
                return self.bad_sub(rc);
            }
            self.cursors[mc].flags |= C_INITIALIZED;
        }
        (rc, reserve_loc)
    }

    fn bad_sub(&mut self, rc: i32) -> (i32, Option<Loc>) {
        let rc = if rc == MDB_KEYEXIST {
            MDB_CORRUPTED
        } else {
            rc
        };
        self.t_flags |= MDB_TXN_ERROR;
        (rc, None)
    }

    // ---- delete ----------------------------------------------------------------------------

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:7068-7167 (mdb_cursor_del)
    fn cursor_del(&mut self, c: usize, flags: u32) -> i32 {
        if self.t_flags & MDB_TXN_BLOCKED != 0 {
            return MDB_BAD_TXN;
        }
        if self.cursors[c].flags & C_INITIALIZED == 0 {
            return EINVAL;
        }
        if self.top_ki(c) >= numkeys(self.page(self.top_page(c))) {
            return MDB_NOTFOUND;
        }
        if flags & MDB_NOSPILL == 0 {
            let rc = self.page_spill(c, None);
            if rc != MDB_SUCCESS {
                return rc;
            }
        }
        let rc = self.cursor_touch(c);
        if rc != MDB_SUCCESS {
            return rc;
        }
        let mp = self.top_page(c);
        if !is_leaf2(self.page(mp)) {
            let leaf = nodeptr(self.page(mp), self.top_ki(c));
            let lf = node_flags(self.page(mp), leaf);
            if lf & F_DUPDATA != 0 {
                return self.unreachable("mdb_cursor_del of a dup item");
            } else if (lf ^ flags as u16) & F_SUBDATA != 0 {
                self.t_flags |= MDB_TXN_ERROR;
                return MDB_INCOMPATIBLE;
            }
            if lf & F_BIGDATA != 0 {
                let pg = get_u64(self.page(mp), node_data(self.page(mp), leaf));
                if let Err(rc) = self.page_get(pg) {
                    self.t_flags |= MDB_TXN_ERROR;
                    return rc;
                }
                let rc = self.ovpage_free(c, pg);
                if rc != MDB_SUCCESS {
                    self.t_flags |= MDB_TXN_ERROR;
                    return rc;
                }
            }
        }
        self.cursor_del0(c)
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:8424-8520 (mdb_cursor_del0)
    // The only cursor that deletes is the freelist cursor of `mdb_freelist_save`, which is not
    // in `mt_cursors`, and FREE_DBI has no other cursors: the fixup loops find no cursor.
    fn cursor_del0(&mut self, c: usize) -> i32 {
        let pad = self.db(c).pad as usize;
        self.node_del(c, pad);
        self.db_mut(c).entries -= 1;
        let dbi = self.cursors[c].dbi;
        if !self.tracked[dbi].is_empty() {
            return self.unreachable("mdb_cursor_del0 with other cursors");
        }
        let rc = self.rebalance(c);
        if rc == MDB_SUCCESS {
            if self.cursors[c].snum == 0 {
                return rc;
            }
            self.cursors[c].flags |= C_DEL;
        }
        if rc != MDB_SUCCESS {
            self.t_flags |= MDB_TXN_ERROR;
        }
        rc
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:8236-8350
    // (mdb_rebalance, root cases)
    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:858-859,863
    // FREE_DBI stays one leaf (and its overflow pages) in makeblastdb's sequence; a root
    // branch with one key and the rebalance of a non-root page are not ported.
    fn rebalance(&mut self, c: usize) -> i32 {
        let tp = self.top_page(c);
        let p = self.page(tp);
        let (minkeys, thresh) = if is_branch(p) { (2, 1) } else { (1, 250) };
        let room = (self.psize - PAGEHDRSZ) as i64;
        let fill = 1000 * (room - sizeleft(p) as i64) / room;
        if fill >= thresh && numkeys(p) >= minkeys {
            return MDB_SUCCESS;
        }
        if self.cursors[c].snum < 2 {
            let mp = self.cursors[c].pg[0];
            let p = self.page(mp);
            if is_subp(p) {
                return MDB_SUCCESS;
            }
            if numkeys(p) == 0 {
                let pg = pgno(p);
                let db = self.db_mut(c);
                db.root = P_INVALID;
                db.depth = 0;
                db.leaf_pages = 0;
                midl_append(&mut self.t_free_pgs, pg);
                let k = &mut self.cursors[c];
                k.snum = 0;
                k.top = 0;
                k.flags &= !C_INITIALIZED;
            } else if is_branch(p) && numkeys(p) == 1 {
                return self.unreachable("mdb_rebalance of a root branch with one key");
            }
            return MDB_SUCCESS;
        }
        self.unreachable("mdb_rebalance of a non-root page")
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:5673-5747 (mdb_ovpage_free)
    fn ovpage_free(&mut self, c: usize, pg: u64) -> i32 {
        let p = self.page_by_pgno(pg);
        let ovpages = pages(p) as usize;
        let dirty = flags_of(p) & P_DIRTY != 0;
        if self.pghead.is_some() && dirty {
            let x = self
                .dirty_index(pg)
                .expect("LMDB writer: dirty overflow page");
            let (_, dp) = self.t_dirty.remove(x);
            self.t_dirty_room += 1;
            self.dpage_free(dp);
            let mop = self.pghead.as_mut().unwrap();
            let n = mop[0] as usize;
            if mop.len() < n + ovpages + 1 {
                mop.resize(n + ovpages + 1, 0);
            }
            let mut j = n + ovpages;
            let mut i = n;
            while i > 0 && mop[i] < pg {
                mop[j] = mop[i];
                j -= 1;
                i -= 1;
            }
            let mut pgv = pg;
            while j > i {
                mop[j] = pgv;
                j -= 1;
                pgv += 1;
            }
            mop[0] += ovpages as u64;
        } else {
            midl_append_range(&mut self.t_free_pgs, pg, ovpages);
        }
        self.db_mut(c).overflow_pages -= ovpages as u64;
        MDB_SUCCESS
    }

    // ---- nodes and pages -------------------------------------------------------------------

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:7169-7205 (mdb_page_new)
    // `mp_pad` is not set: a recycled buffer keeps its old value.
    fn page_new(&mut self, c: usize, flags: u16, num: usize) -> Result<u64, i32> {
        let np = self.page_alloc(num)?;
        let psize = self.psize;
        let p = self.page_mut(PageRef::Pg(np));
        set_flags(p, flags | P_DIRTY);
        set_lower(p, PAGEHDRSZ as u16);
        set_upper(p, psize as u16);
        let (br, lf, ov) = (is_branch(p), is_leaf(p), is_overflow(p));
        if ov && !br && !lf {
            set_pages(p, num as u32);
        }
        let db = self.db_mut(c);
        if br {
            db.branch_pages += 1;
        } else if lf {
            db.leaf_pages += 1;
        } else if ov {
            db.overflow_pages += num as u64;
        }
        Ok(np)
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:7207-7230 (mdb_leaf_size)
    fn leaf_size(&self, ksize: usize, dsize: usize) -> usize {
        let mut sz = NODESIZE + ksize + dsize;
        if sz > self.nodemax {
            sz -= dsize - 8;
        }
        even(sz + 2)
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:7232-7253 (mdb_branch_size)
    fn branch_size(ksize: usize) -> usize {
        NODESIZE + ksize + 2
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:7255-7391 (mdb_node_add)
    /// `mdb_node_add` on the cursor's top page; `dbytes` is the data (the 8-byte page number
    /// for `F_BIGDATA`, unused for `MDB_RESERVE`). For `MDB_RESERVE`, where the data goes.
    #[allow(clippy::too_many_arguments)]
    fn node_add(
        &mut self,
        c: usize,
        indx: usize,
        key: Option<&[u8]>,
        dsize: usize,
        dbytes: &[u8],
        pgno_v: u64,
        flags: u32,
    ) -> (i32, Option<Loc>) {
        let mut flags = flags;
        let mp = self.top_page(c);
        let ks = self.db(c).pad as usize;
        if is_leaf2(self.page(mp)) {
            let key = key.unwrap();
            let p = self.page_mut(mp);
            let ptr = leaf2key(indx, ks);
            let dif = numkeys(p) as isize - indx as isize;
            if dif > 0 {
                p.copy_within(ptr..ptr + dif as usize * ks, ptr + ks);
            }
            p[ptr..ptr + ks].copy_from_slice(&key[..ks]);
            let lo = lower(p) + 2;
            set_lower(p, lo);
            let up = upper(p).wrapping_sub((ks - 2) as u16);
            set_upper(p, up);
            return (MDB_SUCCESS, None);
        }
        let (room, leaf) = {
            let p = self.page(mp);
            (sizeleft(p) as isize - 2, is_leaf(p))
        };
        let mut node_size = NODESIZE + key.map_or(0, |k| k.len());
        let mut ofp: Option<u64> = None;
        let mut evened = false;
        if leaf {
            if flags & F_BIGDATA as u32 != 0 {
                node_size += 8;
            } else if node_size + dsize > self.nodemax {
                let ovpages = (PAGEHDRSZ - 1 + dsize) / self.psize + 1;
                node_size = even(node_size + 8);
                if node_size as isize > room {
                    self.t_flags |= MDB_TXN_ERROR;
                    return (MDB_PAGE_FULL, None);
                }
                match self.page_new(c, P_OVERFLOW, ovpages) {
                    Ok(p) => ofp = Some(p),
                    Err(rc) => return (rc, None),
                }
                flags |= F_BIGDATA as u32;
                evened = true;
            } else {
                node_size += dsize;
            }
        }
        if !evened {
            node_size = even(node_size);
            if node_size as isize > room {
                self.t_flags |= MDB_TXN_ERROR;
                return (MDB_PAGE_FULL, None);
            }
        }
        let p = self.page_mut(mp);
        let n = numkeys(p);
        let mut i = n;
        while i > indx {
            let v = ptr(p, i - 1);
            set_ptr(p, i, v);
            i -= 1;
        }
        let ofs = upper(p) - node_size as u16;
        set_ptr(p, indx, ofs);
        set_upper(p, ofs);
        let lo = lower(p) + 2;
        set_lower(p, lo);
        let node = ofs as usize;
        let ksize = key.map_or(0, |k| k.len());
        put_u16(p, node + 6, ksize as u16);
        set_node_flags(p, node, flags as u16);
        if leaf {
            set_node_dsz(p, node, dsize);
        } else {
            set_node_pgno(p, node, pgno_v);
        }
        if let Some(k) = key {
            p[node_key(node)..node_key(node) + ksize].copy_from_slice(k);
        }
        let mut loc = None;
        if leaf {
            let ndata = node_key(node) + ksize;
            match ofp {
                None => {
                    if flags & F_BIGDATA as u32 != 0 {
                        p[ndata..ndata + 8].copy_from_slice(&dbytes[..8]);
                    } else if flags & MDB_RESERVE != 0 {
                        loc = Some(Loc {
                            page: mp,
                            off: ndata,
                        });
                    } else {
                        p[ndata..ndata + dsize].copy_from_slice(dbytes);
                    }
                }
                Some(o) => {
                    p[ndata..ndata + 8].copy_from_slice(&o.to_le_bytes());
                    let op = PageRef::Pg(o);
                    if flags & MDB_RESERVE != 0 {
                        loc = Some(Loc {
                            page: op,
                            off: PAGEHDRSZ,
                        });
                    } else {
                        self.page_mut(op)[PAGEHDRSZ..PAGEHDRSZ + dsize].copy_from_slice(dbytes);
                    }
                }
            }
        }
        (MDB_SUCCESS, loc)
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:7393-7448 (mdb_node_del)
    fn node_del(&mut self, c: usize, ksize: usize) {
        let mp = self.top_page(c);
        let indx = self.top_ki(c);
        let p = self.page_mut(mp);
        let numkeys = numkeys(p);
        if is_leaf2(p) {
            let x = numkeys - 1 - indx;
            let base = leaf2key(indx, ksize);
            if x > 0 {
                p.copy_within(base + ksize..base + ksize + x * ksize, base);
            }
            let lo = lower(p) - 2;
            set_lower(p, lo);
            let up = upper(p) + (ksize - 2) as u16;
            set_upper(p, up);
            return;
        }
        let node = nodeptr(p, indx);
        let mut sz = NODESIZE + node_ksize(p, node);
        if is_leaf(p) {
            if node_flags(p, node) & F_BIGDATA != 0 {
                sz += 8;
            } else {
                sz += node_dsz(p, node);
            }
        }
        let sz = even(sz);
        let pv = ptr(p, indx);
        let mut j = 0;
        for i in 0..numkeys {
            if i != indx {
                let v = ptr(p, i);
                set_ptr(p, j, if v < pv { v + sz as u16 } else { v });
                j += 1;
            }
        }
        let up = upper(p) as usize;
        p.copy_within(up..pv as usize, up + sz);
        let lo = lower(p) - 2;
        set_lower(p, lo);
        set_upper(p, (up + sz) as u16);
    }

    // ---- split -----------------------------------------------------------------------------

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:8582-9006
    // (mdb_page_split); the parts are cited where they are ported below.
    /// `mdb_page_split`; for `MDB_RESERVE`, where the data goes.
    fn page_split(
        &mut self,
        mc: usize,
        newkey: &[u8],
        newdata: Option<Val>,
        newpgno: u64,
        nflags: u32,
    ) -> (i32, Option<Loc>) {
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:8613-8665
        let psize = self.psize;
        let mp = self.top_page(mc);
        let mp_pgno = pgno(self.page(mp));
        let newindx = self.top_ki(mc);
        let mut nkeys = numkeys(self.page(mp));
        let mp_flags = flags_of(self.page(mp));
        let mp_pad = pad(self.page(mp));
        let mp_leaf = is_leaf(self.page(mp));
        let mut copy: Option<u16> = None;
        let mut new_root: Option<usize> = None;
        let mut did_split = false;
        let mut reserve: Option<Loc> = None;
        let rp_pgno = match self.page_new(mc, mp_flags, 1) {
            Ok(p) => p,
            Err(rc) => return (rc, None),
        };
        let rp = PageRef::Pg(rp_pgno);
        set_pad(self.page_mut(rp), mp_pad);
        let mut ptop;
        let rc: i32 = 'done: {
            if self.cursors[mc].top < 1 {
                let pp = match self.page_new(mc, P_BRANCH, 1) {
                    Ok(p) => p,
                    Err(rc) => break 'done rc,
                };
                let k = &mut self.cursors[mc];
                let mut i = k.snum;
                while i > 0 {
                    k.pg[i] = k.pg[i - 1];
                    k.ki[i] = k.ki[i - 1];
                    i -= 1;
                }
                k.pg[0] = PageRef::Pg(pp);
                k.ki[0] = 0;
                let db = self.db_mut(mc);
                db.root = pp;
                new_root = Some(db.depth as usize);
                db.depth += 1;
                let (rc, _) = self.node_add(mc, 0, None, 0, &[], mp_pgno, 0);
                if rc != MDB_SUCCESS {
                    let k = &mut self.cursors[mc];
                    k.pg[0] = k.pg[1];
                    k.ki[0] = k.ki[1];
                    let db = self.db_mut(mc);
                    db.root = mp_pgno;
                    db.depth -= 1;
                    break 'done rc;
                }
                let k = &mut self.cursors[mc];
                k.snum += 1;
                k.top += 1;
                ptop = 0;
            } else {
                ptop = self.cursors[mc].top - 1;
            }
            let mn = self.cursor_copy(mc);
            {
                let mki = self.cursors[mc].ki[ptop];
                let k = &mut self.cursors[mn];
                k.pg[k.top] = rp;
                k.ki[ptop] = mki + 1;
            }
            // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:8715-8795
            let split_indx: usize;
            let sepkey: Vec<u8>;
            if nflags & MDB_APPEND != 0 {
                let k = &mut self.cursors[mn];
                k.ki[k.top] = 0;
                sepkey = newkey.to_vec();
                split_indx = newindx;
                nkeys = 0;
            } else {
                // split_indx = (nkeys+1) / 2
                let mut si = nkeys.div_ceil(2);
                if is_leaf2(self.page(rp)) {
                    break 'done self.unreachable("split of a LEAF2 page without MDB_APPEND");
                }
                let pmax = (psize - PAGEHDRSZ) as isize;
                let nsize = even(if mp_leaf {
                    self.leaf_size(newkey.len(), newdata.map_or(0, |d| d.size()))
                } else {
                    Self::branch_size(newkey.len())
                }) as isize;
                let mut cbuf = self.page_malloc(1);
                set_pgno(&mut cbuf, mp_pgno);
                set_flags(&mut cbuf, mp_flags);
                set_lower(&mut cbuf, PAGEHDRSZ as u16);
                set_upper(&mut cbuf, psize as u16);
                {
                    let src = self.page(mp);
                    let mut j = 0;
                    for i in 0..nkeys {
                        if i == newindx {
                            set_ptr(&mut cbuf, j, 0);
                            j += 1;
                        }
                        set_ptr(&mut cbuf, j, ptr(src, i));
                        j += 1;
                    }
                }
                self.scratch.push(cbuf);
                let ck = (self.scratch.len() - 1) as u16;
                copy = Some(ck);
                let cp = PageRef::Scratch(ck);
                if nkeys < 32 || nsize > pmax / 16 || newindx >= nkeys {
                    let mut ps: isize = 0;
                    let (mut i, j, k): (isize, isize, isize) = if newindx <= si || newindx >= nkeys
                    {
                        (
                            0,
                            1,
                            if newindx >= nkeys {
                                nkeys as isize
                            } else {
                                si as isize + 1 + mp_leaf as isize
                            },
                        )
                    } else {
                        (nkeys as isize, -1, si as isize - 1)
                    };
                    while i != k {
                        if i as usize == newindx {
                            ps += nsize;
                        } else {
                            let src = self.page(mp);
                            let node = ptr(self.page(cp), i as usize) as usize;
                            ps += (NODESIZE + node_ksize(src, node) + 2) as isize;
                            if mp_leaf {
                                if node_flags(src, node) & F_BIGDATA != 0 {
                                    ps += 8;
                                } else {
                                    ps += node_dsz(src, node) as isize;
                                }
                            }
                            ps = even(ps as usize) as isize;
                        }
                        if ps > pmax || i == k - j {
                            si = (i + (j < 0) as isize) as usize;
                            break;
                        }
                        i += j;
                    }
                }
                split_indx = si;
                if split_indx == newindx {
                    sepkey = newkey.to_vec();
                } else {
                    let src = self.page(mp);
                    let node = ptr(self.page(cp), split_indx) as usize;
                    let o = node_key(node);
                    sepkey = src[o..o + node_ksize(src, node)].to_vec();
                }
            }
            // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:8802-8842
            let pp = self.cursors[mn].pg[ptop];
            if (sizeleft(self.page(pp)) as usize) < Self::branch_size(sepkey.len()) {
                let snum = self.cursors[mc].snum;
                {
                    let k = &mut self.cursors[mn];
                    k.snum -= 1;
                    k.top -= 1;
                }
                did_split = true;
                let dbi = self.cursors[mn].dbi;
                let tr = if self.cursors[mn].flags & C_SUB != 0 {
                    Track::Dummy(mn)
                } else {
                    Track::Cursor(mn)
                };
                self.tracked[dbi].push(tr);
                let (rc, _) = self.page_split(mn, &sepkey, None, rp_pgno, 0);
                self.tracked[dbi].pop();
                if rc != MDB_SUCCESS {
                    break 'done rc;
                }
                if self.cursors[mc].snum > snum {
                    ptop += 1;
                }
                let (mnp, mcp, mcki) = (
                    self.cursors[mn].pg[ptop],
                    self.cursors[mc].pg[ptop],
                    self.cursors[mc].ki[ptop] as usize,
                );
                if mnp != mcp && mcki >= numkeys(self.page(mcp)) {
                    for i in 0..ptop {
                        let (pg, ki) = (self.cursors[mn].pg[i], self.cursors[mn].ki[i]);
                        self.cursors[mc].pg[i] = pg;
                        self.cursors[mc].ki[i] = ki;
                    }
                    self.cursors[mc].pg[ptop] = mnp;
                    let mnki = self.cursors[mn].ki[ptop];
                    if mnki != 0 {
                        self.cursors[mc].ki[ptop] = mnki - 1;
                    } else {
                        self.cursors[mc].ki[ptop] = mnki;
                        self.cursor_sibling(mc, false);
                    }
                }
            } else {
                self.cursors[mn].top -= 1;
                let ki = self.cursors[mn].ki[ptop] as usize;
                let (rc, _) = self.node_add(mn, ki, Some(&sepkey), 0, &[], rp_pgno, 0);
                self.cursors[mn].top += 1;
                if rc != MDB_SUCCESS {
                    break 'done rc;
                }
            }
            let mctop = self.cursors[mc].top;
            if nflags & MDB_APPEND != 0 {
                self.cursors[mc].pg[mctop] = rp;
                self.cursors[mc].ki[mctop] = 0;
                let (dsize, dbytes): (usize, &[u8]) = match newdata {
                    Some(Val::Bytes(b)) => (b.len(), b),
                    Some(Val::Reserve(n)) => (n, &[]),
                    None => (0, &[]),
                };
                let (rc, loc) = self.node_add(mc, 0, Some(newkey), dsize, dbytes, newpgno, nflags);
                if rc != MDB_SUCCESS {
                    break 'done rc;
                }
                reserve = loc;
                for i in 0..mctop {
                    let ki = self.cursors[mn].ki[i];
                    self.cursors[mc].ki[i] = ki;
                }
            } else {
                // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:8853-8927
                // The source nodes are read from `mp`, which changes only after the loop; each
                // is copied out before `mdb_node_add` writes to `rp` or `copy`.
                let cp = PageRef::Scratch(copy.unwrap());
                self.cursors[mc].pg[mctop] = rp;
                let mut i = split_indx;
                let mut j = 0usize;
                let mut kbuf: Vec<u8> = Vec::new();
                let mut dbuf: Vec<u8> = Vec::new();
                loop {
                    let (flags, pg, dsize, rkey_len);
                    let use_new = i == newindx;
                    if use_new {
                        kbuf.clear();
                        kbuf.extend_from_slice(newkey);
                        dbuf.clear();
                        let mut ds = 0;
                        if mp_leaf {
                            match newdata {
                                Some(Val::Bytes(b)) => {
                                    dbuf.extend_from_slice(b);
                                    ds = b.len();
                                }
                                Some(Val::Reserve(n)) => ds = n,
                                None => {}
                            }
                        }
                        dsize = ds;
                        pg = newpgno;
                        flags = nflags;
                        self.cursors[mc].ki[mctop] = j as u16;
                    } else {
                        let node = ptr(self.page(cp), i) as usize;
                        let src = self.page(mp);
                        let ks = node_ksize(src, node);
                        kbuf.clear();
                        kbuf.extend_from_slice(&src[node_key(node)..node_key(node) + ks]);
                        dbuf.clear();
                        let fl = node_flags(src, node);
                        if mp_leaf {
                            let d = node_data(src, node);
                            let ds = node_dsz(src, node);
                            let stored = if fl & F_BIGDATA != 0 { 8 } else { ds };
                            dbuf.extend_from_slice(&src[d..d + stored]);
                            dsize = ds;
                            pg = 0;
                        } else {
                            dsize = 0;
                            pg = node_pgno(src, node);
                        }
                        flags = fl as u32;
                    }
                    rkey_len = if !mp_leaf && j == 0 { 0 } else { kbuf.len() };
                    let (rc, loc) =
                        self.node_add(mc, j, Some(&kbuf[..rkey_len]), dsize, &dbuf, pg, flags);
                    if rc != MDB_SUCCESS {
                        break 'done rc;
                    }
                    if use_new {
                        reserve = loc;
                    }
                    if i == nkeys {
                        i = 0;
                        j = 0;
                        self.cursors[mc].pg[mctop] = cp;
                    } else {
                        i += 1;
                        j += 1;
                    }
                    if i == split_indx {
                        break;
                    }
                }
                {
                    let ck = copy.unwrap() as usize;
                    let cbuf = std::mem::take(&mut self.scratch[ck]);
                    let cn = numkeys(&cbuf);
                    let cup = upper(&cbuf) as usize;
                    let p = self.page_mut(mp);
                    for i in 0..cn {
                        set_ptr(p, i, ptr(&cbuf, i));
                    }
                    set_lower(p, lower(&cbuf));
                    set_upper(p, upper(&cbuf));
                    let at = if cn > 0 {
                        ptr(p, cn - 1) as usize
                    } else {
                        upper(p) as usize
                    };
                    let n = psize - cup;
                    let from = if cn > 0 {
                        ptr(&cbuf, cn - 1) as usize
                    } else {
                        cup
                    };
                    p[at..at + n].copy_from_slice(&cbuf[from..from + n]);
                    self.scratch[ck] = cbuf;
                }
                if newindx < split_indx {
                    self.cursors[mc].pg[mctop] = mp;
                } else {
                    self.cursors[mc].pg[mctop] = rp;
                    self.cursors[mc].ki[ptop] += 1;
                    let (mnp, mcp, mcki) = (
                        self.cursors[mn].pg[ptop],
                        self.cursors[mc].pg[ptop],
                        self.cursors[mc].ki[ptop] as usize,
                    );
                    if mnp != mcp && mcki >= numkeys(self.page(mcp)) {
                        for i in 0..=ptop {
                            let (pg, ki) = (self.cursors[mn].pg[i], self.cursors[mn].ki[i]);
                            self.cursors[mc].pg[i] = pg;
                            self.cursors[mc].ki[i] = ki;
                        }
                    }
                }
                if nflags & MDB_RESERVE != 0 {
                    let tp = self.top_page(mc);
                    let p = self.page(tp);
                    let node = nodeptr(p, self.top_ki(mc));
                    if node_flags(p, node) & F_BIGDATA == 0 {
                        reserve = Some(Loc {
                            page: tp,
                            off: node_data(p, node),
                        });
                    }
                }
            }
            // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:8956-8995
            // `new_root` is the depth before a root split (`md_depth` of `mc->mc_db` is incremented
            // and `new_root` takes the old value), so at least 1 when the root was split; `None`
            // stands for the C value 0.
            {
                let dbi = self.cursors[mc].dbi;
                let sub = self.cursors[mc].flags & C_SUB != 0;
                let nk = numkeys(self.page(mp)) as u16;
                let mctop = self.cursors[mc].top;
                for k in 0..self.tracked[dbi].len() {
                    let Some((m2flags, _, m3)) = self.tracked_at(dbi, k, sub) else {
                        continue;
                    };
                    if m3 == mc {
                        continue;
                    }
                    if m2flags & self.cursors[m3].flags & C_INITIALIZED == 0 {
                        continue;
                    }
                    if let Some(nr) = new_root.filter(|&d| d != 0) {
                        if self.cursors[m3].pg[0] != mp {
                            continue;
                        }
                        let k3 = &mut self.cursors[m3];
                        let mut kk = nr as isize;
                        while kk >= 0 {
                            let u = kk as usize;
                            k3.ki[u + 1] = k3.ki[u];
                            k3.pg[u + 1] = k3.pg[u];
                            kk -= 1;
                        }
                        k3.ki[0] = if k3.ki[0] >= nk { 1 } else { 0 };
                        let root = self.cursors[mc].pg[0];
                        let k3 = &mut self.cursors[m3];
                        k3.pg[0] = root;
                        k3.snum += 1;
                        k3.top += 1;
                    }
                    let (m3top, m3pg_top) = (self.cursors[m3].top, self.cursors[m3].pg[mctop]);
                    if m3top >= mctop && m3pg_top == mp {
                        let k3 = &mut self.cursors[m3];
                        if k3.ki[mctop] as usize >= newindx && nflags & MDB_SPLIT_REPLACE == 0 {
                            k3.ki[mctop] += 1;
                        }
                        if k3.ki[mctop] >= nk {
                            k3.pg[mctop] = rp;
                            k3.ki[mctop] -= nk;
                            for i in 0..mctop {
                                let (ki, pg) = (self.cursors[mn].ki[i], self.cursors[mn].pg[i]);
                                self.cursors[m3].ki[i] = ki;
                                self.cursors[m3].pg[i] = pg;
                            }
                        }
                    } else if !did_split
                        && m3top >= ptop
                        && self.cursors[m3].pg[ptop] == self.cursors[mc].pg[ptop]
                        && self.cursors[m3].ki[ptop] >= self.cursors[mc].ki[ptop]
                    {
                        self.cursors[m3].ki[ptop] += 1;
                    }
                    if mp_leaf {
                        let pg = self.cursors[m3].pg[mctop];
                        self.xcursor_refresh(m3, mctop, pg);
                    }
                }
            }
            MDB_SUCCESS
        };
        if let Some(ck) = copy {
            debug_assert_eq!(ck as usize, self.scratch.len() - 1);
            let cbuf = self.scratch.pop().unwrap();
            self.page_free(cbuf);
        }
        if rc != MDB_SUCCESS {
            self.t_flags |= MDB_TXN_ERROR;
            return (rc, None);
        }
        (rc, reserve)
    }

    // ---- freelist --------------------------------------------------------------------------

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3073-3293 (mdb_freelist_save); the
    // parts are cited where they are ported below. There are no loose pages, so the parts of
    // `mt_loose_pgs` (3098-3138, 3246-3261) never run.
    pub(super) fn freelist_save(&mut self) -> i32 {
        let maxfree_1pg = self.maxfree_1pg;
        let mut more: i64 = 1;
        let mut pglast: u64 = 0;
        let mut head_id: u64 = 0;
        let mut freecnt: u64 = 0;
        let mut head_room: i64 = 0;
        let mut total_room: i64 = 0;
        let mut mop_len: i64;
        let clean_limit = maxfree_1pg;
        let mc = self.cursor_init(FREE_DBI);
        if self.pghead.is_some() {
            let rc = self.page_search(mc, None, MDB_PS_FIRST | MDB_PS_MODIFY);
            if rc != MDB_SUCCESS && rc != MDB_NOTFOUND {
                return rc;
            }
        }
        loop {
            // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3152-3162
            while pglast < self.pglast {
                let (rc, key, _) = self.cursor_first(mc, false);
                if rc != MDB_SUCCESS {
                    return rc;
                }
                pglast = get_u64(self.page(key.unwrap().page), key.unwrap().off);
                head_id = pglast;
                total_room = 0;
                head_room = 0;
                let rc = self.cursor_del(mc, 0);
                if rc != MDB_SUCCESS {
                    return rc;
                }
            }
            // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3167-3189
            if freecnt < self.t_free_pgs[0] {
                if freecnt == 0 {
                    let rc = self.page_search(mc, None, MDB_PS_LAST | MDB_PS_MODIFY);
                    if rc != MDB_SUCCESS && rc != MDB_NOTFOUND {
                        return rc;
                    }
                }
                let key = self.t_txnid.to_le_bytes();
                let mut size;
                let mut loc;
                loop {
                    freecnt = self.t_free_pgs[0];
                    size = midl_sizeof(&self.t_free_pgs);
                    let (rc, l) = self.cursor_put(mc, &key, Val::Reserve(size), MDB_RESERVE);
                    if rc != MDB_SUCCESS {
                        return rc;
                    }
                    loc = l.expect("LMDB writer: MDB_RESERVE location");
                    if freecnt >= self.t_free_pgs[0] {
                        break;
                    }
                }
                midl_sort(&mut self.t_free_pgs);
                let n = size / 8;
                let ids: Vec<u64> = self.t_free_pgs[..n].to_vec();
                let dst = &mut self.page_mut(loc.page)[loc.off..loc.off + size];
                for (k, v) in ids.iter().enumerate() {
                    put_u64(dst, 8 * k, *v);
                }
                continue;
            }
            // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3203-3240
            mop_len = self.pghead.as_ref().map_or(0, |m| m[0] as i64);
            if total_room >= mop_len {
                if total_room == mop_len {
                    break;
                }
                more -= 1;
                if more < 0 {
                    break;
                }
            } else if head_room >= maxfree_1pg && head_id > 1 {
                head_id -= 1;
                head_room = 0;
            }
            total_room -= head_room;
            head_room = mop_len - total_room;
            if head_room > maxfree_1pg && head_id > 1 {
                head_room = (head_room as u64 / head_id) as i64;
                head_room += maxfree_1pg - head_room % (maxfree_1pg + 1);
            } else if head_room < 0 {
                head_room = 0;
            }
            let key = head_id.to_le_bytes();
            let size = (head_room as usize + 1) * 8;
            let (rc, l) = self.cursor_put(mc, &key, Val::Reserve(size), MDB_RESERVE);
            if rc != MDB_SUCCESS {
                return rc;
            }
            let loc = l.expect("LMDB writer: MDB_RESERVE location");
            let jmax = if head_room > clean_limit {
                head_room as usize
            } else {
                0
            };
            let dst = &mut self.page_mut(loc.page)[loc.off..];
            for j in 0..=jmax {
                put_u64(dst, 8 * j, 0);
            }
            total_room += head_room;
        }
        // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:3265-3288
        // The puts of MDB_CURRENT allocate no page (the records were reserved above), so
        // `me_pghead` does not change while it is read.
        let mut rc = MDB_SUCCESS;
        if mop_len != 0 {
            let mut pos = mop_len as usize;
            let (r, mut key, mut data) = self.cursor_first(mc, true);
            rc = r;
            while rc == MDB_SUCCESS {
                let k = key.unwrap();
                let id = get_u64(self.page(k.page), k.off);
                let mut len = (data.unwrap().len / 8) as i64 - 1;
                if len > mop_len {
                    len = mop_len;
                }
                pos -= len as usize;
                let mut bytes = Vec::with_capacity((len as usize + 1) * 8);
                bytes.extend_from_slice(&(len as u64).to_le_bytes());
                {
                    let mop = self.pghead.as_ref().unwrap();
                    for v in &mop[pos + 1..=pos + len as usize] {
                        bytes.extend_from_slice(&v.to_le_bytes());
                    }
                }
                let (r, _) =
                    self.cursor_put(mc, &id.to_le_bytes(), Val::Bytes(&bytes), MDB_CURRENT);
                rc = r;
                if rc != MDB_SUCCESS {
                    break;
                }
                mop_len -= len;
                if mop_len == 0 {
                    break;
                }
                let (r, k2, d2) = self.cursor_next(mc, true);
                rc = r;
                key = k2;
                data = d2;
            }
        }
        self.cursors.truncate(mc);
        rc
    }

    // ---- databases -------------------------------------------------------------------------

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:9712-9728 (mdb_default_cmp)
    fn default_cmp(&mut self, dbi: usize) -> i32 {
        let f = self.t_dbs[dbi].flags as u32;
        if f & (MDB_REVERSEKEY | MDB_REVERSEDUP) != 0 {
            return self.unreachable("MDB_REVERSEKEY/MDB_REVERSEDUP");
        }
        self.dbx_cmp[dbi] = if f & MDB_INTEGERKEY != 0 {
            Cmp::Int
        } else {
            Cmp::Memn
        };
        self.dbx_dcmp[dbi] = if f & MDB_DUPSORT == 0 {
            Cmp::None
        } else if f & MDB_INTEGERDUP != 0 {
            Cmp::Int
        } else {
            Cmp::Memn
        };
        MDB_SUCCESS
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:9730-9847 (mdb_dbi_open, named DB)
    /// `mdb_dbi_open` of `txn`, `name` and `flags` (lmdb++ `dbi::open`).
    pub(crate) fn dbi_open(&mut self, name: &str, flags: u32) -> Result<Dbi, LmdbError> {
        let valid = MDB_REVERSEKEY
            | MDB_DUPSORT
            | MDB_INTEGERKEY
            | MDB_DUPFIXED
            | MDB_INTEGERDUP
            | MDB_REVERSEDUP
            | MDB_CREATE;
        if flags & !valid != 0 {
            return Err(self.error("mdb_dbi_open", EINVAL));
        }
        if self.t_flags & MDB_TXN_BLOCKED != 0 {
            return Err(self.error("mdb_dbi_open", MDB_BAD_TXN));
        }
        if self.dbx_cmp[MAIN_DBI] == Cmp::None {
            let rc = self.default_cmp(MAIN_DBI);
            if rc != MDB_SUCCESS {
                return Err(self.error("mdb_dbi_open", rc));
            }
        }
        let name = name.as_bytes();
        let mut unused = 0;
        for i in CORE_DBS..self.t_numdbs {
            if self.dbx_name[i].is_empty() {
                if unused == 0 {
                    unused = i;
                }
                continue;
            }
            if self.dbx_name[i] == name {
                return Ok(Dbi(i));
            }
        }
        if unused == 0 && self.t_numdbs >= self.maxdbs {
            return Err(self.error("mdb_dbi_open", MDB_DBS_FULL));
        }
        if self.t_dbs[MAIN_DBI].flags as u32 & (MDB_DUPSORT | MDB_INTEGERKEY) != 0 {
            let rc = if flags & MDB_CREATE != 0 {
                MDB_INCOMPATIBLE
            } else {
                MDB_NOTFOUND
            };
            return Err(self.error("mdb_dbi_open", rc));
        }
        let mut dbflag = DB_NEW | DB_VALID | DB_USRVALID;
        let mut exact = false;
        let mc = self.cursor_init(MAIN_DBI);
        let (mut rc, _, data) = self.cursor_set(mc, name, true, CursorOp::Set, Some(&mut exact));
        let db: MdbDb;
        if rc == MDB_SUCCESS {
            let tp = self.top_page(mc);
            let p = self.page(tp);
            let node = nodeptr(p, self.top_ki(mc));
            if node_flags(p, node) & (F_DUPDATA | F_SUBDATA) != F_SUBDATA {
                return Err(self.error("mdb_dbi_open", MDB_INCOMPATIBLE));
            }
            db = MdbDb::from_bytes(self.slot_bytes(data.unwrap()));
        } else {
            if rc != MDB_NOTFOUND || flags & MDB_CREATE == 0 {
                return Err(self.error("mdb_dbi_open", rc));
            }
            let dummy = MdbDb {
                root: P_INVALID,
                flags: (flags & PERSISTENT_FLAGS as u32) as u16,
                ..MdbDb::default()
            };
            self.tracked[MAIN_DBI].push(Track::Cursor(mc));
            let (r, _) = self.cursor_put(mc, name, Val::Bytes(&dummy.to_bytes()), F_SUBDATA as u32);
            self.tracked[MAIN_DBI].pop();
            rc = r;
            dbflag |= DB_DIRTY;
            db = dummy;
        }
        self.cursors.truncate(mc);
        if rc != MDB_SUCCESS {
            return Err(self.error("mdb_dbi_open", rc));
        }
        let slot = if unused != 0 { unused } else { self.t_numdbs };
        self.dbx_name[slot] = name.to_vec();
        self.t_dbflags[slot] = dbflag;
        self.t_dbs[slot] = db;
        let rc = self.default_cmp(slot);
        if rc != MDB_SUCCESS {
            return Err(self.error("mdb_dbi_open", rc));
        }
        if unused == 0 {
            self.t_numdbs += 1;
        }
        Ok(Dbi(slot))
    }

    // LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/mdb.c:9008-9031 (mdb_put)
    // lmdb++ (NCBI 598d8ae6 copy): c++/include/util/lmdbxx/lmdb++.h:850-861 (lmdb::dbi_put)
    /// `lmdb::dbi_put`: `false` when the key/data pair exists (`MDB_KEYEXIST`).
    pub(crate) fn put(
        &mut self,
        dbi: Dbi,
        key: &[u8],
        data: &[u8],
        flags: u32,
    ) -> Result<bool, LmdbError> {
        let dbi = dbi.0;
        if dbi >= self.t_numdbs || self.t_dbflags[dbi] & DB_USRVALID == 0 {
            return Err(self.error("mdb_put", EINVAL));
        }
        if flags & !(MDB_NOOVERWRITE | MDB_NODUPDATA | MDB_RESERVE | MDB_APPEND | MDB_APPENDDUP)
            != 0
        {
            return Err(self.error("mdb_put", EINVAL));
        }
        if self.t_flags & MDB_TXN_BLOCKED != 0 {
            return Err(self.error("mdb_put", MDB_BAD_TXN));
        }
        let mc = self.cursor_init(dbi);
        self.tracked[dbi].push(Track::Cursor(mc));
        let (rc, _) = self.cursor_put(mc, key, Val::Bytes(data), flags);
        self.tracked[dbi].pop();
        self.cursors.truncate(mc);
        if rc != MDB_SUCCESS && rc != MDB_KEYEXIST {
            return Err(self.error("mdb_put", rc));
        }
        Ok(rc == MDB_SUCCESS)
    }
}

/// `mp_flags` (named apart from the local variables called `flags`).
#[inline]
fn flags_of(p: &[u8]) -> u16 {
    super::page::flags(p)
}
