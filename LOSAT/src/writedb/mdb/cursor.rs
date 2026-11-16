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

    /// The `MDB_page *` of the data of a node (`NODEDATA(node)` taken as a sub-page).
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:7596-7624 (mdb_cursor_init)
    // ```c
    // 	mc->mc_next = NULL;
    // 	mc->mc_backup = NULL;
    // 	mc->mc_dbi = dbi;
    // 	mc->mc_txn = txn;
    // 	mc->mc_db = &txn->mt_dbs[dbi];
    // 	mc->mc_dbx = &txn->mt_dbxs[dbi];
    // 	mc->mc_dbflag = &txn->mt_dbflags[dbi];
    // 	mc->mc_snum = 0;
    // 	mc->mc_top = 0;
    // 	mc->mc_pg[0] = 0;
    // 	mc->mc_ki[0] = 0;
    // 	mc->mc_flags = 0;
    // 	if (txn->mt_dbs[dbi].md_flags & MDB_DUPSORT) {
    // 		mdb_tassert(txn, mx != NULL);
    // 		mc->mc_xcursor = mx;
    // 		mdb_xcursor_init0(mc);
    // 	} else {
    // 		mc->mc_xcursor = NULL;
    // 	}
    // 	if (*mc->mc_dbflag & DB_STALE) {
    // 		mdb_page_search(mc, NULL, MDB_PS_ROOTONLY);
    // 	}
    // ```
    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:7505-7522 (mdb_xcursor_init0)
    // ```c
    // 	mx->mx_cursor.mc_xcursor = NULL;
    // 	mx->mx_cursor.mc_txn = mc->mc_txn;
    // 	mx->mx_cursor.mc_db = &mx->mx_db;
    // 	mx->mx_cursor.mc_dbx = &mx->mx_dbx;
    // 	mx->mx_cursor.mc_dbi = mc->mc_dbi;
    // 	mx->mx_cursor.mc_dbflag = &mx->mx_dbflag;
    // 	mx->mx_cursor.mc_snum = 0;
    // 	mx->mx_cursor.mc_top = 0;
    // 	mx->mx_cursor.mc_flags = C_SUB;
    // 	mx->mx_dbx.md_name.mv_size = 0;
    // 	mx->mx_dbx.md_name.mv_data = NULL;
    // 	mx->mx_dbx.md_cmp = mc->mc_dbx->md_dcmp;
    // 	mx->mx_dbx.md_dcmp = NULL;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:8214-8234 (mdb_cursor_copy)
    // ```c
    // 	cdst->mc_txn = csrc->mc_txn;
    // 	cdst->mc_dbi = csrc->mc_dbi;
    // 	cdst->mc_db  = csrc->mc_db;
    // 	cdst->mc_dbx = csrc->mc_dbx;
    // 	cdst->mc_snum = csrc->mc_snum;
    // 	cdst->mc_top = csrc->mc_top;
    // 	cdst->mc_flags = csrc->mc_flags;
    // 	for (i=0; i<csrc->mc_snum; i++) {
    // 		cdst->mc_pg[i] = csrc->mc_pg[i];
    // 		cdst->mc_ki[i] = csrc->mc_ki[i];
    // 	}
    // ```
    /// `mdb_cursor_copy` into a new cursor without an xcursor (`mn.mc_xcursor = NULL`).
    fn cursor_copy(&mut self, src: usize) -> usize {
        let mut k = self.cursors[src].clone();
        k.xc = None;
        let c = self.cursors.len();
        self.cursors.push(k);
        c
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:7524-7563 (mdb_xcursor_init1)
    // ```c
    // 	if (node->mn_flags & F_SUBDATA) {
    // 		memcpy(&mx->mx_db, NODEDATA(node), sizeof(MDB_db));
    // 		mx->mx_cursor.mc_pg[0] = 0;
    // 		mx->mx_cursor.mc_snum = 0;
    // 		mx->mx_cursor.mc_top = 0;
    // 		mx->mx_cursor.mc_flags = C_SUB;
    // 	} else {
    // 		MDB_page *fp = NODEDATA(node);
    // 		mx->mx_db.md_pad = 0;
    // 		mx->mx_db.md_flags = 0;
    // 		mx->mx_db.md_depth = 1;
    // 		mx->mx_db.md_branch_pages = 0;
    // 		mx->mx_db.md_leaf_pages = 1;
    // 		mx->mx_db.md_overflow_pages = 0;
    // 		mx->mx_db.md_entries = NUMKEYS(fp);
    // 		COPY_PGNO(mx->mx_db.md_root, fp->mp_pgno);
    // 		mx->mx_cursor.mc_snum = 1;
    // 		mx->mx_cursor.mc_top = 0;
    // 		mx->mx_cursor.mc_flags = C_INITIALIZED|C_SUB;
    // 		mx->mx_cursor.mc_pg[0] = fp;
    // 		mx->mx_cursor.mc_ki[0] = 0;
    // 		if (mc->mc_db->md_flags & MDB_DUPFIXED) {
    // 			mx->mx_db.md_flags = MDB_DUPFIXED;
    // 			mx->mx_db.md_pad = fp->mp_pad;
    // 			if (mc->mc_db->md_flags & MDB_INTEGERDUP)
    // 				mx->mx_db.md_flags |= MDB_INTEGERKEY;
    // 		}
    // 	}
    // 	mx->mx_dbflag = DB_VALID|DB_USRVALID|DB_DUPDATA;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:7565-7594 (mdb_xcursor_init2)
    // ```c
    // 	if (new_dupdata) {
    // 		mx->mx_cursor.mc_snum = 1;
    // 		mx->mx_cursor.mc_top = 0;
    // 		mx->mx_cursor.mc_flags |= C_INITIALIZED;
    // 		mx->mx_cursor.mc_ki[0] = 0;
    // 		mx->mx_dbflag = DB_VALID|DB_USRVALID|DB_DUPDATA;
    // 	} else if (!(mx->mx_cursor.mc_flags & C_INITIALIZED)) {
    // 		return;
    // 	}
    // 	mx->mx_db = src_mx->mx_db;
    // 	mx->mx_cursor.mc_pg[0] = src_mx->mx_cursor.mc_pg[0];
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1243-1250 (XCURSOR_REFRESH)
    // ```c
    // #define XCURSOR_REFRESH(mc, top, mp) do { \
    // 	MDB_page *xr_pg = (mp); \
    // 	MDB_node *xr_node; \
    // 	if (!XCURSOR_INITED(mc) || (mc)->mc_ki[top] >= NUMKEYS(xr_pg)) break; \
    // 	xr_node = NODEPTR(xr_pg, (mc)->mc_ki[top]); \
    // 	if ((xr_node->mn_flags & (F_DUPDATA|F_SUBDATA)) == F_DUPDATA) \
    // 		(mc)->mc_xcursor->mx_cursor.mc_pg[0] = NODEDATA(xr_node); \
    // } while (0)
    // ```
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
    /// `m3 = (mc->mc_flags & C_SUB) ? &m2->mc_xcursor->mx_cursor : m2`.
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5268-5365 (mdb_node_search)
    // ```c
    // 	nkeys = NUMKEYS(mp);
    // 	low = IS_LEAF(mp) ? 0 : 1;
    // 	high = nkeys - 1;
    // 	cmp = mc->mc_dbx->md_cmp;
    // 	...
    // 	if (IS_LEAF2(mp)) {
    // 		nodekey.mv_size = mc->mc_db->md_pad;
    // 		node = NODEPTR(mp, 0);	/* fake */
    // 		while (low <= high) {
    // 			i = (low + high) >> 1;
    // 			nodekey.mv_data = LEAF2KEY(mp, i, nodekey.mv_size);
    // 			rc = cmp(key, &nodekey);
    // 			if (rc == 0)
    // 				break;
    // 			if (rc > 0)
    // 				low = i + 1;
    // 			else
    // 				high = i - 1;
    // 		}
    // 	} else {
    // 		while (low <= high) {
    // 			i = (low + high) >> 1;
    // 			node = NODEPTR(mp, i);
    // 			nodekey.mv_size = NODEKSZ(node);
    // 			nodekey.mv_data = NODEKEY(node);
    // 			rc = cmp(key, &nodekey);
    // 			...
    // 		}
    // 	}
    // 	if (rc > 0) {	/* Found entry is less than the key. */
    // 		i++;	/* Skip to get the smallest entry larger than key. */
    // 		if (!IS_LEAF2(mp))
    // 			node = NODEPTR(mp, i);
    // 	}
    // 	if (exactp)
    // 		*exactp = (rc == 0 && nkeys > 0);
    // 	/* store the key index */
    // 	mc->mc_ki[mc->mc_top] = i;
    // 	if (i >= nkeys)
    // 		/* There is no entry larger or equal to the key. */
    // 		return NULL;
    // 	/* nodeptr is fake for LEAF2 */
    // 	return node;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5381-5417 (mdb_cursor_pop,
    // mdb_cursor_push)
    // ```c
    // 	if (mc->mc_snum) {
    // 		mc->mc_snum--;
    // 		if (mc->mc_snum) {
    // 			mc->mc_top--;
    // 		} else {
    // 			mc->mc_flags &= ~C_INITIALIZED;
    // 		}
    // 	}
    // ...
    // 	if (mc->mc_snum >= CURSOR_STACK) {
    // 		mc->mc_txn->mt_flags |= MDB_TXN_ERROR;
    // 		return MDB_CURSOR_FULL;
    // 	}
    // 	mc->mc_top = mc->mc_snum++;
    // 	mc->mc_pg[mc->mc_top] = mp;
    // 	mc->mc_ki[mc->mc_top] = 0;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5480-5550 (mdb_page_search_root)
    // ```c
    // 	while (IS_BRANCH(mp)) {
    // 		MDB_node	*node;
    // 		indx_t		i;
    // 		...
    // 		if (flags & (MDB_PS_FIRST|MDB_PS_LAST)) {
    // 			i = 0;
    // 			if (flags & MDB_PS_LAST) {
    // 				i = NUMKEYS(mp) - 1;
    // 				/* if already init'd, see if we're already in right place */
    // 				if (mc->mc_flags & C_INITIALIZED) {
    // 					if (mc->mc_ki[mc->mc_top] == i) {
    // 						mc->mc_top = mc->mc_snum++;
    // 						mp = mc->mc_pg[mc->mc_top];
    // 						goto ready;
    // 					}
    // 				}
    // 			}
    // 		} else {
    // 			int	 exact;
    // 			node = mdb_node_search(mc, key, &exact);
    // 			if (node == NULL)
    // 				i = NUMKEYS(mp) - 1;
    // 			else {
    // 				i = mc->mc_ki[mc->mc_top];
    // 				if (!exact) {
    // 					mdb_cassert(mc, i > 0);
    // 					i--;
    // 				}
    // 			}
    // 		}
    // 		mdb_cassert(mc, i < NUMKEYS(mp));
    // 		node = NODEPTR(mp, i);
    // 		if ((rc = mdb_page_get(mc, NODEPGNO(node), &mp, NULL)) != 0)
    // 			return rc;
    // 		mc->mc_ki[mc->mc_top] = i;
    // 		if ((rc = mdb_cursor_push(mc, mp)))
    // 			return rc;
    // ready:
    // 		if (flags & MDB_PS_MODIFY) {
    // 			if ((rc = mdb_page_touch(mc)) != 0)
    // 				return rc;
    // 			mp = mc->mc_pg[mc->mc_top];
    // 		}
    // 	}
    // 	if (!IS_LEAF(mp)) {
    // 		mc->mc_txn->mt_flags |= MDB_TXN_ERROR;
    // 		return MDB_CORRUPTED;
    // 	}
    // 	mc->mc_flags |= C_INITIALIZED;
    // 	mc->mc_flags &= ~C_EOF;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5578-5671 (mdb_page_search)
    // ```c
    // 	if (mc->mc_txn->mt_flags & MDB_TXN_BLOCKED) {
    // 		return MDB_BAD_TXN;
    // 	} else {
    // 		/* Make sure we're using an up-to-date root */
    // 		if (*mc->mc_dbflag & DB_STALE) {
    // 				MDB_cursor mc2;
    // 				...
    // 				mdb_cursor_init(&mc2, mc->mc_txn, MAIN_DBI, NULL);
    // 				rc = mdb_page_search(&mc2, &mc->mc_dbx->md_name, 0);
    // 				if (rc)
    // 					return rc;
    // 				{
    // 					MDB_val data;
    // 					int exact = 0;
    // 					uint16_t flags;
    // 					MDB_node *leaf = mdb_node_search(&mc2,
    // 						&mc->mc_dbx->md_name, &exact);
    // 					if (!exact)
    // 						return MDB_NOTFOUND;
    // 					if ((leaf->mn_flags & (F_DUPDATA|F_SUBDATA)) != F_SUBDATA)
    // 						return MDB_INCOMPATIBLE; /* not a named DB */
    // 					rc = mdb_node_read(&mc2, leaf, &data);
    // 					if (rc)
    // 						return rc;
    // 					memcpy(&flags, ((char *) data.mv_data + offsetof(MDB_db, md_flags)),
    // 						sizeof(uint16_t));
    // 					if ((mc->mc_db->md_flags & PERSISTENT_FLAGS) != flags)
    // 						return MDB_INCOMPATIBLE;
    // 					memcpy(mc->mc_db, data.mv_data, sizeof(MDB_db));
    // 				}
    // 				*mc->mc_dbflag &= ~DB_STALE;
    // 		}
    // 		root = mc->mc_db->md_root;
    // 		if (root == P_INVALID) {		/* Tree is empty. */
    // 			return MDB_NOTFOUND;
    // 		}
    // 	}
    // 	if (!mc->mc_pg[0] || mc->mc_pg[0]->mp_pgno != root)
    // 		if ((rc = mdb_page_get(mc, root, &mc->mc_pg[0], NULL)) != 0)
    // 			return rc;
    // 	mc->mc_snum = 1;
    // 	mc->mc_top = 0;
    // 	if (flags & MDB_PS_MODIFY) {
    // 		if ((rc = mdb_page_touch(mc)))
    // 			return rc;
    // 	}
    // 	if (flags & MDB_PS_ROOTONLY)
    // 		return MDB_SUCCESS;
    // 	return mdb_page_search_root(mc, key, flags);
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5749-5779 (mdb_node_read)
    // ```c
    // 	if (!F_ISSET(leaf->mn_flags, F_BIGDATA)) {
    // 		data->mv_size = NODEDSZ(leaf);
    // 		data->mv_data = NODEDATA(leaf);
    // 		return MDB_SUCCESS;
    // 	}
    // 	/* Read overflow data.
    // 	 */
    // 	data->mv_size = NODEDSZ(leaf);
    // 	memcpy(&pgno, NODEDATA(leaf), sizeof(pgno));
    // 	if ((rc = mdb_page_get(mc, pgno, &omp, NULL)) != 0) {
    // 		return rc;
    // 	}
    // 	data->mv_data = METADATA(omp);
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5805-5858 (mdb_cursor_sibling)
    // ```c
    // 	if (mc->mc_snum < 2) {
    // 		return MDB_NOTFOUND;		/* root has no siblings */
    // 	}
    // 	mdb_cursor_pop(mc);
    // 	if (move_right ? (mc->mc_ki[mc->mc_top] + 1u >= NUMKEYS(mc->mc_pg[mc->mc_top]))
    // 		       : (mc->mc_ki[mc->mc_top] == 0)) {
    // 		if ((rc = mdb_cursor_sibling(mc, move_right)) != MDB_SUCCESS) {
    // 			/* undo cursor_pop before returning */
    // 			mc->mc_top++;
    // 			mc->mc_snum++;
    // 			return rc;
    // 		}
    // 	} else {
    // 		if (move_right)
    // 			mc->mc_ki[mc->mc_top]++;
    // 		else
    // 			mc->mc_ki[mc->mc_top]--;
    // 	}
    // 	mdb_cassert(mc, IS_BRANCH(mc->mc_pg[mc->mc_top]));
    // 	indx = NODEPTR(mc->mc_pg[mc->mc_top], mc->mc_ki[mc->mc_top]);
    // 	if ((rc = mdb_page_get(mc, NODEPGNO(indx), &mp, NULL)) != 0) {
    // 		/* mc will be inconsistent if caller does mc_snum++ as above */
    // 		mc->mc_flags &= ~(C_INITIALIZED|C_EOF);
    // 		return rc;
    // 	}
    // 	mdb_cursor_push(mc, mp);
    // 	if (!move_right)
    // 		mc->mc_ki[mc->mc_top] = NUMKEYS(mp)-1;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5860-5939 (mdb_cursor_next)
    // ```c
    // 	if (!(mc->mc_flags & C_INITIALIZED))
    // 		return mdb_cursor_first(mc, key, data);
    // 	mp = mc->mc_pg[mc->mc_top];
    // 	if (mc->mc_flags & C_EOF) {
    // 		if (mc->mc_ki[mc->mc_top] >= NUMKEYS(mp)-1)
    // 			return MDB_NOTFOUND;
    // 		mc->mc_flags ^= C_EOF;
    // 	}
    // 	if (mc->mc_db->md_flags & MDB_DUPSORT) {
    // 		...
    // 	}
    // 	if (mc->mc_flags & C_DEL) {
    // 		mc->mc_flags ^= C_DEL;
    // 		goto skip;
    // 	}
    // 	if (mc->mc_ki[mc->mc_top] + 1u >= NUMKEYS(mp)) {
    // 		if ((rc = mdb_cursor_sibling(mc, 1)) != MDB_SUCCESS) {
    // 			mc->mc_flags |= C_EOF;
    // 			return rc;
    // 		}
    // 		mp = mc->mc_pg[mc->mc_top];
    // 	} else
    // 		mc->mc_ki[mc->mc_top]++;
    // skip:
    // 	if (IS_LEAF2(mp)) {
    // 		key->mv_size = mc->mc_db->md_pad;
    // 		key->mv_data = LEAF2KEY(mp, mc->mc_ki[mc->mc_top], key->mv_size);
    // 		return MDB_SUCCESS;
    // 	}
    // 	mdb_cassert(mc, IS_LEAF(mp));
    // 	leaf = NODEPTR(mp, mc->mc_ki[mc->mc_top]);
    // 	if (F_ISSET(leaf->mn_flags, F_DUPDATA)) {
    // 		mdb_xcursor_init1(mc, leaf);
    // 	}
    // 	if (data) {
    // 		if ((rc = mdb_node_read(mc, leaf, data)) != MDB_SUCCESS)
    // 			return rc;
    // 		...
    // 	}
    // 	MDB_GET_KEY(leaf, key);
    // ```
    /// `mdb_cursor_next(mc, key, data, MDB_NEXT)` on a database without dups.
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:6030-6216 (mdb_cursor_set)
    // ```c
    // 	if (key->mv_size == 0)
    // 		return MDB_BAD_VALSIZE;
    // 	if (mc->mc_xcursor)
    // 		mc->mc_xcursor->mx_cursor.mc_flags &= ~(C_INITIALIZED|C_EOF);
    // 	/* See if we're already on the right page */
    // 	if (mc->mc_flags & C_INITIALIZED) {
    // 		MDB_val nodekey;
    // 		mp = mc->mc_pg[mc->mc_top];
    // 		if (!NUMKEYS(mp)) {
    // 			mc->mc_ki[mc->mc_top] = 0;
    // 			return MDB_NOTFOUND;
    // 		}
    // 		if (mp->mp_flags & P_LEAF2) {
    // 			nodekey.mv_size = mc->mc_db->md_pad;
    // 			nodekey.mv_data = LEAF2KEY(mp, 0, nodekey.mv_size);
    // 		} else {
    // 			leaf = NODEPTR(mp, 0);
    // 			MDB_GET_KEY2(leaf, nodekey);
    // 		}
    // 		rc = mc->mc_dbx->md_cmp(key, &nodekey);
    // 		if (rc == 0) {
    // 			mc->mc_ki[mc->mc_top] = 0;
    // 			if (exactp)
    // 				*exactp = 1;
    // 			goto set1;
    // 		}
    // 		if (rc > 0) {
    // 			unsigned int i;
    // 			unsigned int nkeys = NUMKEYS(mp);
    // 			if (nkeys > 1) {
    // 				(last node of the page: rc == 0 -> set1; rc < 0 -> the current node, then set2)
    // 			}
    // 			/* If any parents have right-sibs, search.
    // 			 * Otherwise, there's nothing further.
    // 			 */
    // 			for (i=0; i<mc->mc_top; i++)
    // 				if (mc->mc_ki[i] <
    // 					NUMKEYS(mc->mc_pg[i])-1)
    // 					break;
    // 			if (i == mc->mc_top) {
    // 				/* There are no other pages */
    // 				mc->mc_ki[mc->mc_top] = nkeys;
    // 				return MDB_NOTFOUND;
    // 			}
    // 		}
    // 		if (!mc->mc_top) {
    // 			/* There are no other pages */
    // 			mc->mc_ki[mc->mc_top] = 0;
    // 			if (op == MDB_SET_RANGE && !exactp) {
    // 				rc = 0;
    // 				goto set1;
    // 			} else
    // 				return MDB_NOTFOUND;
    // 		}
    // 	} else {
    // 		mc->mc_pg[0] = 0;
    // 	}
    // 	rc = mdb_page_search(mc, key, 0);
    // 	if (rc != MDB_SUCCESS)
    // 		return rc;
    // 	mp = mc->mc_pg[mc->mc_top];
    // set2:
    // 	leaf = mdb_node_search(mc, key, exactp);
    // 	if (exactp != NULL && !*exactp) {
    // 		/* MDB_SET specified and not an exact match. */
    // 		return MDB_NOTFOUND;
    // 	}
    // 	if (leaf == NULL) {
    // 		if ((rc = mdb_cursor_sibling(mc, 1)) != MDB_SUCCESS) {
    // 			mc->mc_flags |= C_EOF;
    // 			return rc;		/* no entries matched */
    // 		}
    // 		mp = mc->mc_pg[mc->mc_top];
    // 		leaf = NODEPTR(mp, 0);
    // 	}
    // set1:
    // 	mc->mc_flags |= C_INITIALIZED;
    // 	mc->mc_flags &= ~C_EOF;
    // 	if (IS_LEAF2(mp)) {
    // 		if (op == MDB_SET_RANGE || op == MDB_SET_KEY) {
    // 			key->mv_size = mc->mc_db->md_pad;
    // 			key->mv_data = LEAF2KEY(mp, mc->mc_ki[mc->mc_top], key->mv_size);
    // 		}
    // 		return MDB_SUCCESS;
    // 	}
    // 	if (F_ISSET(leaf->mn_flags, F_DUPDATA)) {
    // 		mdb_xcursor_init1(mc, leaf);
    // 	}
    // 	if (data) {
    // 		if (F_ISSET(leaf->mn_flags, F_DUPDATA)) {
    // 			if (op == MDB_SET || op == MDB_SET_KEY || op == MDB_SET_RANGE) {
    // 				rc = mdb_cursor_first(&mc->mc_xcursor->mx_cursor, data, NULL);
    // 			} else { ... }
    // 		} else if (op == MDB_GET_BOTH || op == MDB_GET_BOTH_RANGE) {
    // 			...
    // 		} else {
    // 			if (mc->mc_xcursor)
    // 				mc->mc_xcursor->mx_cursor.mc_flags &= ~(C_INITIALIZED|C_EOF);
    // 			if ((rc = mdb_node_read(mc, leaf, data)) != MDB_SUCCESS)
    // 				return rc;
    // 		}
    // 	}
    // 	/* The key already matches in all other cases */
    // 	if (op == MDB_SET_RANGE || op == MDB_SET_KEY)
    // 		MDB_GET_KEY(leaf, key);
    // 	return rc;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:6218-6256 (mdb_cursor_first)
    // ```c
    // 	if (mc->mc_xcursor)
    // 		mc->mc_xcursor->mx_cursor.mc_flags &= ~(C_INITIALIZED|C_EOF);
    // 	if (!(mc->mc_flags & C_INITIALIZED) || mc->mc_top) {
    // 		rc = mdb_page_search(mc, NULL, MDB_PS_FIRST);
    // 		if (rc != MDB_SUCCESS)
    // 			return rc;
    // 	}
    // 	mdb_cassert(mc, IS_LEAF(mc->mc_pg[mc->mc_top]));
    // 	leaf = NODEPTR(mc->mc_pg[mc->mc_top], 0);
    // 	mc->mc_flags |= C_INITIALIZED;
    // 	mc->mc_flags &= ~C_EOF;
    // 	mc->mc_ki[mc->mc_top] = 0;
    // 	if (IS_LEAF2(mc->mc_pg[mc->mc_top])) {
    // 		key->mv_size = mc->mc_db->md_pad;
    // 		key->mv_data = LEAF2KEY(mc->mc_pg[mc->mc_top], 0, key->mv_size);
    // 		return MDB_SUCCESS;
    // 	}
    // 	if (data) {
    // 		if (F_ISSET(leaf->mn_flags, F_DUPDATA)) {
    // 			mdb_xcursor_init1(mc, leaf);
    // 			rc = mdb_cursor_first(&mc->mc_xcursor->mx_cursor, data, NULL);
    // 			if (rc)
    // 				return rc;
    // 		} else {
    // 			if ((rc = mdb_node_read(mc, leaf, data)) != MDB_SUCCESS)
    // 				return rc;
    // 		}
    // 	}
    // 	MDB_GET_KEY(leaf, key);
    // ```
    pub(super) fn cursor_first(
        &mut self,
        c: usize,
        want_data: bool,
    ) -> (i32, Option<Slot>, Option<Slot>) {
        self.cursor_end(c, want_data, false)
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:6258-6296 (mdb_cursor_last)
    // ```c
    // 	if (mc->mc_xcursor)
    // 		mc->mc_xcursor->mx_cursor.mc_flags &= ~(C_INITIALIZED|C_EOF);
    // 	if (!(mc->mc_flags & C_INITIALIZED) || mc->mc_top) {
    // 		rc = mdb_page_search(mc, NULL, MDB_PS_LAST);
    // 		if (rc != MDB_SUCCESS)
    // 			return rc;
    // 	}
    // 	mdb_cassert(mc, IS_LEAF(mc->mc_pg[mc->mc_top]));
    // 	mc->mc_ki[mc->mc_top] = NUMKEYS(mc->mc_pg[mc->mc_top]) - 1;
    // 	mc->mc_flags |= C_INITIALIZED|C_EOF;
    // 	leaf = NODEPTR(mc->mc_pg[mc->mc_top], mc->mc_ki[mc->mc_top]);
    // 	if (IS_LEAF2(mc->mc_pg[mc->mc_top])) {
    // 		key->mv_size = mc->mc_db->md_pad;
    // 		key->mv_data = LEAF2KEY(mc->mc_pg[mc->mc_top], mc->mc_ki[mc->mc_top], key->mv_size);
    // 		return MDB_SUCCESS;
    // 	}
    // 	if (data) {
    // 		if (F_ISSET(leaf->mn_flags, F_DUPDATA)) {
    // 			mdb_xcursor_init1(mc, leaf);
    // 			rc = mdb_cursor_last(&mc->mc_xcursor->mx_cursor, data, NULL);
    // 			if (rc)
    // 				return rc;
    // 		} else {
    // 			if ((rc = mdb_node_read(mc, leaf, data)) != MDB_SUCCESS)
    // 				return rc;
    // 		}
    // 	}
    // 	MDB_GET_KEY(leaf, key);
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:6505-6534 (mdb_cursor_touch)
    // ```c
    // 	if (mc->mc_dbi >= CORE_DBS && !(*mc->mc_dbflag & (DB_DIRTY|DB_DUPDATA))) {
    // 		/* Touch DB record of named DB */
    // 		MDB_cursor mc2;
    // 		MDB_xcursor mcx;
    // 		if (TXN_DBI_CHANGED(mc->mc_txn, mc->mc_dbi))
    // 			return MDB_BAD_DBI;
    // 		mdb_cursor_init(&mc2, mc->mc_txn, MAIN_DBI, &mcx);
    // 		rc = mdb_page_search(&mc2, &mc->mc_dbx->md_name, MDB_PS_MODIFY);
    // 		if (rc)
    // 			 return rc;
    // 		*mc->mc_dbflag |= DB_DIRTY;
    // 	}
    // 	mc->mc_top = 0;
    // 	if (mc->mc_snum) {
    // 		do {
    // 			rc = mdb_page_touch(mc);
    // 		} while (!rc && ++(mc->mc_top) < mc->mc_snum);
    // 		mc->mc_top = mc->mc_snum-1;
    // 	}
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:2410-2510 (mdb_page_touch)
    // ```c
    // 	if (!F_ISSET(mp->mp_flags, P_DIRTY)) {
    // 		if (txn->mt_flags & MDB_TXN_SPILLS) { ... }
    // 		if ((rc = mdb_midl_need(&txn->mt_free_pgs, 1)) ||
    // 			(rc = mdb_page_alloc(mc, 1, &np)))
    // 			goto fail;
    // 		pgno = np->mp_pgno;
    // 		mdb_cassert(mc, mp->mp_pgno != pgno);
    // 		mdb_midl_xappend(txn->mt_free_pgs, mp->mp_pgno);
    // 		/* Update the parent page, if any, to point to the new page */
    // 		if (mc->mc_top) {
    // 			MDB_page *parent = mc->mc_pg[mc->mc_top-1];
    // 			MDB_node *node = NODEPTR(parent, mc->mc_ki[mc->mc_top-1]);
    // 			SETPGNO(node, pgno);
    // 		} else {
    // 			mc->mc_db->md_root = pgno;
    // 		}
    // 	} else if (txn->mt_parent && !IS_SUBP(mp)) {
    // 		...
    // 	} else {
    // 		return 0;
    // 	}
    // 	mdb_page_copy(np, mp, txn->mt_env->me_psize);
    // 	np->mp_pgno = pgno;
    // 	np->mp_flags |= P_DIRTY;
    // done:
    // 	/* Adjust cursors pointing to mp */
    // 	mc->mc_pg[mc->mc_top] = np;
    // 	m2 = txn->mt_cursors[mc->mc_dbi];
    // 	if (mc->mc_flags & C_SUB) {
    // 		for (; m2; m2=m2->mc_next) {
    // 			m3 = &m2->mc_xcursor->mx_cursor;
    // 			if (m3->mc_snum < mc->mc_snum) continue;
    // 			if (m3->mc_pg[mc->mc_top] == mp)
    // 				m3->mc_pg[mc->mc_top] = np;
    // 		}
    // 	} else {
    // 		for (; m2; m2=m2->mc_next) {
    // 			if (m2->mc_snum < mc->mc_snum) continue;
    // 			if (m2 == mc) continue;
    // 			if (m2->mc_pg[mc->mc_top] == mp) {
    // 				m2->mc_pg[mc->mc_top] = np;
    // 				if (IS_LEAF(np))
    // 					XCURSOR_REFRESH(m2, mc->mc_top, np);
    // 			}
    // 		}
    // 	}
    // 	return 0;
    // fail:
    // 	txn->mt_flags |= MDB_TXN_ERROR;
    // 	return rc;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:1997-2024 (mdb_page_spill)
    // ```c
    // 	if (m0->mc_flags & C_SUB)
    // 		return MDB_SUCCESS;
    // 	/* Estimate how much space this op will take */
    // 	i = m0->mc_db->md_depth;
    // 	/* Named DBs also dirty the main DB */
    // 	if (m0->mc_dbi >= CORE_DBS)
    // 		i += txn->mt_dbs[MAIN_DBI].md_depth;
    // 	/* For puts, roughly factor in the key+data size */
    // 	if (key)
    // 		i += (LEAFSIZE(key, data) + txn->mt_env->me_psize) / txn->mt_env->me_psize;
    // 	i += i;	/* double it for good measure */
    // 	need = i;
    // 	if (txn->mt_dirty_room > i)
    // 		return MDB_SUCCESS;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:6540-7066 (mdb_cursor_put); the
    // parts are quoted where they are ported below.
    /// `mdb_cursor_put`: the return code and, for `MDB_RESERVE`, where the data goes.
    pub(super) fn cursor_put(
        &mut self,
        mc: usize,
        key: &[u8],
        data: Val,
        flags: u32,
    ) -> (i32, Option<Loc>) {
        // ```c
        // 	if (flags & MDB_MULTIPLE) { ... }
        // 	nospill = flags & MDB_NOSPILL;
        // 	flags &= ~MDB_NOSPILL;
        // 	if (mc->mc_txn->mt_flags & (MDB_TXN_RDONLY|MDB_TXN_BLOCKED))
        // 		return (mc->mc_txn->mt_flags & MDB_TXN_RDONLY) ? EACCES : MDB_BAD_TXN;
        // 	if (key->mv_size-1 >= ENV_MAXKEY(env))
        // 		return MDB_BAD_VALSIZE;
        // #if SIZE_MAX > MAXDATASIZE
        // 	if (data->mv_size > ((mc->mc_db->md_flags & MDB_DUPSORT) ? ENV_MAXKEY(env) : MAXDATASIZE))
        // 		return MDB_BAD_VALSIZE;
        // ```
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
        // ```c
        // 	dkey.mv_size = 0;
        // 	if (flags == MDB_CURRENT) {
        // 		if (!(mc->mc_flags & C_INITIALIZED))
        // 			return EINVAL;
        // 		rc = MDB_SUCCESS;
        // 	} else if (mc->mc_db->md_root == P_INVALID) {
        // 		/* new database, cursor has nothing to point to */
        // 		mc->mc_snum = 0;
        // 		mc->mc_top = 0;
        // 		mc->mc_flags &= ~C_INITIALIZED;
        // 		rc = MDB_NO_ROOT;
        // 	} else {
        // 		int exact = 0;
        // 		MDB_val d2;
        // 		if (flags & MDB_APPEND) {
        // 			MDB_val k2;
        // 			rc = mdb_cursor_last(mc, &k2, &d2);
        // 			if (rc == 0) {
        // 				rc = mc->mc_dbx->md_cmp(key, &k2);
        // 				if (rc > 0) {
        // 					rc = MDB_NOTFOUND;
        // 					mc->mc_ki[mc->mc_top]++;
        // 				} else {
        // 					/* new key is <= last key */
        // 					rc = MDB_KEYEXIST;
        // 				}
        // 			}
        // 		} else {
        // 			rc = mdb_cursor_set(mc, key, &d2, MDB_SET, &exact);
        // 		}
        // 		if ((flags & MDB_NOOVERWRITE) && rc == 0) {
        // 			*data = d2;
        // 			return MDB_KEYEXIST;
        // 		}
        // 		if (rc && rc != MDB_NOTFOUND)
        // 			return rc;
        // 	}
        // 	if (mc->mc_flags & C_DEL)
        // 		mc->mc_flags ^= C_DEL;
        // ```
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
        // ```c
        // 	/* Cursor is positioned, check for room in the dirty list */
        // 	if (!nospill) {
        // 		...
        // 			rdata = data;
        // 		if ((rc2 = mdb_page_spill(mc, key, rdata)))
        // 			return rc2;
        // 	}
        // 	if (rc == MDB_NO_ROOT) {
        // 		MDB_page *np;
        // 		/* new database, write a root leaf page */
        // 		if ((rc2 = mdb_page_new(mc, P_LEAF, 1, &np))) {
        // 			return rc2;
        // 		}
        // 		mdb_cursor_push(mc, np);
        // 		mc->mc_db->md_root = np->mp_pgno;
        // 		mc->mc_db->md_depth++;
        // 		*mc->mc_dbflag |= DB_DIRTY;
        // 		if ((mc->mc_db->md_flags & (MDB_DUPSORT|MDB_DUPFIXED))
        // 			== MDB_DUPFIXED)
        // 			np->mp_flags |= P_LEAF2;
        // 		mc->mc_flags |= C_INITIALIZED;
        // 	} else {
        // 		/* make sure all cursor pages are writable */
        // 		rc2 = mdb_cursor_touch(mc);
        // 		if (rc2)
        // 			return rc2;
        // 	}
        // 	insert_key = insert_data = rc;
        // ```
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
            // ```c
            // 	if (insert_key) {
            // 		/* The key does not exist */
            // 		if ((mc->mc_db->md_flags & MDB_DUPSORT) &&
            // 			LEAFSIZE(key, data) > env->me_nodemax)
            // 		{
            // 			/* Too big for a node, insert in sub-DB. ... */
            // 			goto prep_subDB;
            // 		}
            // ```
            // A key of at most 511 bytes and a dup of at most 511 bytes always fit a node of a
            // page of 4096 bytes or more.
            if dbflags & MDB_DUPSORT != 0 && NODESIZE + key.len() + data.size() > self.nodemax {
                return (self.unreachable("new DUPSORT key too big for a node"), None);
            }
        } else {
            // ```c
            // 	} else {
            // 		/* there's only a key anyway, so this is a no-op */
            // 		if (IS_LEAF2(mc->mc_pg[mc->mc_top])) {
            // 			... (overwrite a LEAF2 key, fix_parent)
            // 		}
            // more:
            // 		leaf = NODEPTR(mc->mc_pg[mc->mc_top], mc->mc_ki[mc->mc_top]);
            // 		olddata.mv_size = NODEDSZ(leaf);
            // 		olddata.mv_data = NODEDATA(leaf);
            // ```
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
                // ```c
                // 		/* DB has dups? */
                // 		if (F_ISSET(mc->mc_db->md_flags, MDB_DUPSORT)) {
                // 			unsigned	i, offset = 0;
                // 			mp = fp = xdata.mv_data = env->me_pbuf;
                // 			mp->mp_pgno = mc->mc_pg[mc->mc_top]->mp_pgno;
                // ```
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
                    // ```c
                    // 			/* Was a single item before, must convert now */
                    // 			if (!F_ISSET(leaf->mn_flags, F_DUPDATA)) {
                    // 				MDB_cmp_func *dcmp;
                    // 				/* Just overwrite the current item */
                    // 				if (flags == MDB_CURRENT)
                    // 					goto current;
                    // 				dcmp = mc->mc_dbx->md_dcmp;
                    // 				/* does data match? */
                    // 				if (!dcmp(data, &olddata)) {
                    // 					if (flags & (MDB_NODUPDATA|MDB_APPENDDUP))
                    // 						return MDB_KEYEXIST;
                    // 					/* overwrite it */
                    // 					goto current;
                    // 				}
                    // 				/* Back up original data item */
                    // 				dkey.mv_size = olddata.mv_size;
                    // 				dkey.mv_data = memcpy(fp+1, olddata.mv_data, olddata.mv_size);
                    // 				/* Make sub-page header for the dup items, with dummy body */
                    // 				fp->mp_flags = P_LEAF|P_DIRTY|P_SUBP;
                    // 				fp->mp_lower = (PAGEHDRSZ-PAGEBASE);
                    // 				xdata.mv_size = PAGEHDRSZ + dkey.mv_size + data->mv_size;
                    // 				if (mc->mc_db->md_flags & MDB_DUPFIXED) {
                    // 					fp->mp_flags |= P_LEAF2;
                    // 					fp->mp_pad = data->mv_size;
                    // 					xdata.mv_size += 2 * data->mv_size;	/* leave space for 2 more */
                    // 				} else {
                    // 					...
                    // 				}
                    // 				fp->mp_upper = xdata.mv_size - PAGEBASE;
                    // 				olddata.mv_size = xdata.mv_size; /* pretend olddata is fp */
                    // ```
                    // `fp+1` is `me_pbuf + sizeof(MDB_page)`, 24 bytes on 64-bit hosts.
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
                    // ```c
                    // 			} else if (leaf->mn_flags & F_SUBDATA) {
                    // 				/* Data is on sub-DB, just store it */
                    // 				flags |= F_DUPDATA|F_SUBDATA;
                    // 				goto put_sub;
                    // ```
                    flags |= (F_DUPDATA | F_SUBDATA) as u32;
                    put_sub = true;
                    stage_new_sub = false;
                } else {
                    // ```c
                    // 			} else {
                    // 				/* Data is on sub-page */
                    // 				fp = olddata.mv_data;
                    // 				switch (flags) {
                    // 				default:
                    // 					if (!(mc->mc_db->md_flags & MDB_DUPFIXED)) {
                    // 						offset = EVEN(NODESIZE + sizeof(indx_t) +
                    // 							data->mv_size);
                    // 						break;
                    // 					}
                    // 					offset = fp->mp_pad;
                    // 					if (SIZELEFT(fp) < offset) {
                    // 						offset *= 4; /* space for 4 more */
                    // 						break;
                    // 					}
                    // 					/* FALLTHRU: Big enough MDB_DUPFIXED sub-page */
                    // 				case MDB_CURRENT:
                    // 					fp->mp_flags |= P_DIRTY;
                    // 					COPY_PGNO(fp->mp_pgno, mp->mp_pgno);
                    // 					mc->mc_xcursor->mx_cursor.mc_pg[0] = fp;
                    // 					flags |= F_DUPDATA;
                    // 					goto put_sub;
                    // 				}
                    // 				xdata.mv_size = olddata.mv_size + offset;
                    // 			}
                    // ```
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
                    // ```c
                    // 			fp_flags = fp->mp_flags;
                    // 			if (NODESIZE + NODEKSZ(leaf) + xdata.mv_size > env->me_nodemax) {
                    // 					/* Too big for a sub-page, convert to sub-DB */
                    // 					fp_flags &= ~P_SUBP;
                    // prep_subDB:
                    // 					if (mc->mc_db->md_flags & MDB_DUPFIXED) {
                    // 						fp_flags |= P_LEAF2;
                    // 						dummy.md_pad = fp->mp_pad;
                    // 						dummy.md_flags = MDB_DUPFIXED;
                    // 						if (mc->mc_db->md_flags & MDB_INTEGERDUP)
                    // 							dummy.md_flags |= MDB_INTEGERKEY;
                    // 					} else {
                    // 						dummy.md_pad = 0;
                    // 						dummy.md_flags = 0;
                    // 					}
                    // 					dummy.md_depth = 1;
                    // 					dummy.md_branch_pages = 0;
                    // 					dummy.md_leaf_pages = 1;
                    // 					dummy.md_overflow_pages = 0;
                    // 					dummy.md_entries = NUMKEYS(fp);
                    // 					xdata.mv_size = sizeof(MDB_db);
                    // 					xdata.mv_data = &dummy;
                    // 					if ((rc = mdb_page_alloc(mc, 1, &mp)))
                    // 						return rc;
                    // 					offset = env->me_psize - olddata.mv_size;
                    // 					flags |= F_DUPDATA|F_SUBDATA;
                    // 					dummy.md_root = mp->mp_pgno;
                    // 					sub_root = mp;
                    // 			}
                    // 			if (mp != fp) {
                    // 				mp->mp_flags = fp_flags | P_DIRTY;
                    // 				mp->mp_pad   = fp->mp_pad;
                    // 				mp->mp_lower = fp->mp_lower;
                    // 				mp->mp_upper = fp->mp_upper + offset;
                    // 				if (fp_flags & P_LEAF2) {
                    // 					memcpy(METADATA(mp), METADATA(fp), NUMKEYS(fp) * fp->mp_pad);
                    // 				} else {
                    // 					...
                    // 				}
                    // 			}
                    // 			rdata = &xdata;
                    // 			flags |= F_DUPDATA;
                    // 			do_sub = 1;
                    // 			if (!insert_key)
                    // 				mdb_node_del(mc, 0);
                    // 			goto new_sub;
                    // 		}
                    // ```
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
                // ```c
                // current:
                // 		/* LMDB passes F_SUBDATA in 'flags' to write a DB record */
                // 		if ((leaf->mn_flags ^ flags) & F_SUBDATA)
                // 			return MDB_INCOMPATIBLE;
                // 		/* overflow page overwrites need special handling */
                // 		if (F_ISSET(leaf->mn_flags, F_BIGDATA)) {
                // 			MDB_page *omp;
                // 			pgno_t pg;
                // 			int level, ovpages, dpages = OVPAGES(data->mv_size, env->me_psize);
                // 			memcpy(&pg, olddata.mv_data, sizeof(pg));
                // 			if ((rc2 = mdb_page_get(mc, pg, &omp, &level)) != 0)
                // 				return rc2;
                // 			ovpages = omp->mp_pages;
                // 			/* Is the ov page large enough? */
                // 			if (ovpages >= dpages) {
                // 			  if (!(omp->mp_flags & P_DIRTY) &&
                // 				  (level || (env->me_flags & MDB_WRITEMAP)))
                // 			  { ... unspill ... }
                // 			  /* Is it dirty? */
                // 			  if (omp->mp_flags & P_DIRTY) {
                // 				...
                // 				SETDSZ(leaf, data->mv_size);
                // 				if (F_ISSET(flags, MDB_RESERVE))
                // 					data->mv_data = METADATA(omp);
                // 				else
                // 					memcpy(METADATA(omp), data->mv_data, data->mv_size);
                // 				return MDB_SUCCESS;
                // 			  }
                // 			}
                // 			if ((rc2 = mdb_ovpage_free(mc, omp)) != MDB_SUCCESS)
                // 				return rc2;
                // 		} else if (data->mv_size == olddata.mv_size) {
                // 			/* same size, just replace it. ... */
                // 			if (F_ISSET(flags, MDB_RESERVE))
                // 				data->mv_data = olddata.mv_data;
                // 			else if (!(mc->mc_flags & C_SUB))
                // 				memcpy(olddata.mv_data, data->mv_data, data->mv_size);
                // 			else {
                // 				memcpy(NODEKEY(leaf), key->mv_data, key->mv_size);
                // 				goto fix_parent;
                // 			}
                // 			return MDB_SUCCESS;
                // 		}
                // 		mdb_node_del(mc, 0);
                // 	}
                // 	rdata = data;
                // ```
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
            // ```c
            // new_sub:
            // 	nflags = flags & NODE_ADD_FLAGS;
            // 	nsize = IS_LEAF2(mc->mc_pg[mc->mc_top]) ? key->mv_size : mdb_leaf_size(env, key, rdata);
            // 	if (SIZELEFT(mc->mc_pg[mc->mc_top]) < nsize) {
            // 		if (( flags & (F_DUPDATA|F_SUBDATA)) == F_DUPDATA )
            // 			nflags &= ~MDB_APPEND; /* sub-page may need room to grow */
            // 		if (!insert_key)
            // 			nflags |= MDB_SPLIT_REPLACE;
            // 		rc = mdb_page_split(mc, key, rdata, P_INVALID, nflags);
            // 	} else {
            // 		/* There is room already in this leaf page. */
            // 		rc = mdb_node_add(mc, mc->mc_ki[mc->mc_top], key, rdata, 0, nflags);
            // 		if (rc == 0) {
            // 			/* Adjust other cursors pointing to mp */
            // 			MDB_cursor *m2, *m3;
            // 			MDB_dbi dbi = mc->mc_dbi;
            // 			unsigned i = mc->mc_top;
            // 			MDB_page *mp = mc->mc_pg[i];
            // 			for (m2 = mc->mc_txn->mt_cursors[dbi]; m2; m2=m2->mc_next) {
            // 				if (mc->mc_flags & C_SUB)
            // 					m3 = &m2->mc_xcursor->mx_cursor;
            // 				else
            // 					m3 = m2;
            // 				if (m3 == mc || m3->mc_snum < mc->mc_snum || m3->mc_pg[i] != mp) continue;
            // 				if (m3->mc_ki[i] >= mc->mc_ki[i] && insert_key) {
            // 					m3->mc_ki[i]++;
            // 				}
            // 				XCURSOR_REFRESH(m3, i, mp);
            // 			}
            // 		}
            // 	}
            // ```
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
            // ```c
            // put_sub:
            // 			xdata.mv_size = 0;
            // 			xdata.mv_data = "";
            // 			leaf = NODEPTR(mc->mc_pg[mc->mc_top], mc->mc_ki[mc->mc_top]);
            // 			if (flags & MDB_CURRENT) {
            // 				xflags = MDB_CURRENT|MDB_NOSPILL;
            // 			} else {
            // 				mdb_xcursor_init1(mc, leaf);
            // 				xflags = (flags & MDB_NODUPDATA) ?
            // 					MDB_NOOVERWRITE|MDB_NOSPILL : MDB_NOSPILL;
            // 			}
            // 			if (sub_root)
            // 				mc->mc_xcursor->mx_cursor.mc_pg[0] = sub_root;
            // 			new_dupdata = (int)dkey.mv_size;
            // 			/* converted, write the original data first */
            // 			if (dkey.mv_size) {
            // 				rc = mdb_cursor_put(&mc->mc_xcursor->mx_cursor, &dkey, &xdata, xflags);
            // 				if (rc)
            // 					goto bad_sub;
            // 				/* we've done our job */
            // 				dkey.mv_size = 0;
            // 			}
            // 			if (!(leaf->mn_flags & F_SUBDATA) || sub_root) {
            // 				/* Adjust other cursors pointing to mp */
            // 				MDB_cursor *m2;
            // 				MDB_xcursor *mx = mc->mc_xcursor;
            // 				unsigned i = mc->mc_top;
            // 				MDB_page *mp = mc->mc_pg[i];
            // 				for (m2 = mc->mc_txn->mt_cursors[mc->mc_dbi]; m2; m2=m2->mc_next) {
            // 					if (m2 == mc || m2->mc_snum < mc->mc_snum) continue;
            // 					if (!(m2->mc_flags & C_INITIALIZED)) continue;
            // 					if (m2->mc_pg[i] == mp) {
            // 						if (m2->mc_ki[i] == mc->mc_ki[i]) {
            // 							mdb_xcursor_init2(m2, mx, new_dupdata);
            // 						} else if (!insert_key) {
            // 							XCURSOR_REFRESH(m2, i, mp);
            // 						}
            // 					}
            // 				}
            // 			}
            // 			ecount = mc->mc_xcursor->mx_db.md_entries;
            // 			if (flags & MDB_APPENDDUP)
            // 				xflags |= MDB_APPEND;
            // 			rc = mdb_cursor_put(&mc->mc_xcursor->mx_cursor, data, &xdata, xflags);
            // 			if (flags & F_SUBDATA) {
            // 				void *db = NODEDATA(leaf);
            // 				memcpy(db, &mc->mc_xcursor->mx_db, sizeof(MDB_db));
            // 			}
            // 			insert_data = mc->mc_xcursor->mx_db.md_entries - ecount;
            // 		}
            // ```
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
        // ```c
        // 		/* Increment count unless we just replaced an existing item. */
        // 		if (insert_data)
        // 			mc->mc_db->md_entries++;
        // 		if (insert_key) {
        // 			/* Invalidate txn if we created an empty sub-DB */
        // 			if (rc)
        // 				goto bad_sub;
        // 			/* If we succeeded and the key didn't exist before,
        // 			 * make sure the cursor is marked valid.
        // 			 */
        // 			mc->mc_flags |= C_INITIALIZED;
        // 		}
        // 		...
        // 		return rc;
        // bad_sub:
        // 		if (rc == MDB_KEYEXIST)	/* should not happen, we deleted that item */
        // 			rc = MDB_CORRUPTED;
        // 	}
        // 	mc->mc_txn->mt_flags |= MDB_TXN_ERROR;
        // 	return rc;
        // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:7068-7167 (mdb_cursor_del)
    // ```c
    // 	if (!(mc->mc_flags & C_INITIALIZED))
    // 		return EINVAL;
    // 	if (mc->mc_ki[mc->mc_top] >= NUMKEYS(mc->mc_pg[mc->mc_top]))
    // 		return MDB_NOTFOUND;
    // 	if (!(flags & MDB_NOSPILL) && (rc = mdb_page_spill(mc, NULL, NULL)))
    // 		return rc;
    // 	rc = mdb_cursor_touch(mc);
    // 	if (rc)
    // 		return rc;
    // 	mp = mc->mc_pg[mc->mc_top];
    // 	if (IS_LEAF2(mp))
    // 		goto del_key;
    // 	leaf = NODEPTR(mp, mc->mc_ki[mc->mc_top]);
    // 	if (F_ISSET(leaf->mn_flags, F_DUPDATA)) {
    // 		...
    // 	}
    // 	/* LMDB passes F_SUBDATA in 'flags' to delete a DB record */
    // 	else if ((leaf->mn_flags ^ flags) & F_SUBDATA) {
    // 		rc = MDB_INCOMPATIBLE;
    // 		goto fail;
    // 	}
    // 	/* add overflow pages to free list */
    // 	if (F_ISSET(leaf->mn_flags, F_BIGDATA)) {
    // 		MDB_page *omp;
    // 		pgno_t pg;
    // 		memcpy(&pg, NODEDATA(leaf), sizeof(pg));
    // 		if ((rc = mdb_page_get(mc, pg, &omp, NULL)) ||
    // 			(rc = mdb_ovpage_free(mc, omp)))
    // 			goto fail;
    // 	}
    // del_key:
    // 	return mdb_cursor_del0(mc);
    // fail:
    // 	mc->mc_txn->mt_flags |= MDB_TXN_ERROR;
    // 	return rc;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:8424-8520 (mdb_cursor_del0)
    // ```c
    // 	ki = mc->mc_ki[mc->mc_top];
    // 	mp = mc->mc_pg[mc->mc_top];
    // 	mdb_node_del(mc, mc->mc_db->md_pad);
    // 	mc->mc_db->md_entries--;
    // 	{
    // 		/* Adjust other cursors pointing to mp */
    // 		for (m2 = mc->mc_txn->mt_cursors[dbi]; m2; m2=m2->mc_next) {
    // 			m3 = (mc->mc_flags & C_SUB) ? &m2->mc_xcursor->mx_cursor : m2;
    // 			if (! (m2->mc_flags & m3->mc_flags & C_INITIALIZED))
    // 				continue;
    // 			if (m3 == mc || m3->mc_snum < mc->mc_snum)
    // 				continue;
    // 			if (m3->mc_pg[mc->mc_top] == mp) {
    // 				if (m3->mc_ki[mc->mc_top] == ki) {
    // 					m3->mc_flags |= C_DEL;
    // 					if (mc->mc_db->md_flags & MDB_DUPSORT) {
    // 						/* Sub-cursor referred into dataset which is gone */
    // 						m3->mc_xcursor->mx_cursor.mc_flags &= ~(C_INITIALIZED|C_EOF);
    // 					}
    // 					continue;
    // 				} else if (m3->mc_ki[mc->mc_top] > ki) {
    // 					m3->mc_ki[mc->mc_top]--;
    // 				}
    // 				XCURSOR_REFRESH(m3, mc->mc_top, mp);
    // 			}
    // 		}
    // 	}
    // 	rc = mdb_rebalance(mc);
    // 	if (rc == MDB_SUCCESS) {
    // 		/* DB is totally empty now, just bail out.
    // 		 * Other cursors adjustments were already done
    // 		 * by mdb_rebalance and aren't needed here.
    // 		 */
    // 		if (!mc->mc_snum)
    // 			return rc;
    // 		mp = mc->mc_pg[mc->mc_top];
    // 		nkeys = NUMKEYS(mp);
    // 		/* Adjust other cursors pointing to mp */
    // 		for (m2 = mc->mc_txn->mt_cursors[dbi]; !rc && m2; m2=m2->mc_next) {
    // 			... (cursors past the end move to the next sibling; DUPSORT xcursor fixes)
    // 		}
    // 		mc->mc_flags |= C_DEL;
    // 	}
    // 	if (rc)
    // 		mc->mc_txn->mt_flags |= MDB_TXN_ERROR;
    // 	return rc;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:8236-8350 (mdb_rebalance, root cases)
    // ```c
    // 	if (IS_BRANCH(mc->mc_pg[mc->mc_top])) {
    // 		minkeys = 2;
    // 		thresh = 1;
    // 	} else {
    // 		minkeys = 1;
    // 		thresh = FILL_THRESHOLD;
    // 	}
    // 	if (PAGEFILL(mc->mc_txn->mt_env, mc->mc_pg[mc->mc_top]) >= thresh &&
    // 		NUMKEYS(mc->mc_pg[mc->mc_top]) >= minkeys) {
    // 		return MDB_SUCCESS;
    // 	}
    // 	if (mc->mc_snum < 2) {
    // 		MDB_page *mp = mc->mc_pg[0];
    // 		if (IS_SUBP(mp)) {
    // 			return MDB_SUCCESS;
    // 		}
    // 		if (NUMKEYS(mp) == 0) {
    // 			mc->mc_db->md_root = P_INVALID;
    // 			mc->mc_db->md_depth = 0;
    // 			mc->mc_db->md_leaf_pages = 0;
    // 			rc = mdb_midl_append(&mc->mc_txn->mt_free_pgs, mp->mp_pgno);
    // 			if (rc)
    // 				return rc;
    // 			/* Adjust cursors pointing to mp */
    // 			mc->mc_snum = 0;
    // 			mc->mc_top = 0;
    // 			mc->mc_flags &= ~C_INITIALIZED;
    // 			... (other cursors)
    // 		} else if (IS_BRANCH(mp) && NUMKEYS(mp) == 1) {
    // 			...
    // 		} else
    // 			DPUTS("root page doesn't need rebalancing");
    // 		return MDB_SUCCESS;
    // 	}
    // ```
    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:858-859,863
    // ```c
    // #define PAGEFILL(env, p) (1000L * ((env)->me_psize - PAGEHDRSZ - SIZELEFT(p)) / \
    //				((env)->me_psize - PAGEHDRSZ))
    // #define FILL_THRESHOLD	 250
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:5673-5747 (mdb_ovpage_free)
    // ```c
    // 	pgno_t pg = mp->mp_pgno;
    // 	unsigned x = 0, ovpages = mp->mp_pages;
    // 	...
    // 	if (env->me_pghead &&
    // 		!txn->mt_parent &&
    // 		((mp->mp_flags & P_DIRTY) ||
    // 		 (sl && (x = mdb_midl_search(sl, pn)) <= sl[0] && sl[x] == pn)))
    // 	{
    // 		unsigned i, j;
    // 		pgno_t *mop;
    // 		MDB_ID2 *dl, ix, iy;
    // 		rc = mdb_midl_need(&env->me_pghead, ovpages);
    // 		if (rc)
    // 			return rc;
    // 		...
    // 		/* Remove from dirty list */
    // 		dl = txn->mt_u.dirty_list;
    // 		x = dl[0].mid--;
    // 		for (ix = dl[x]; ix.mptr != mp; ix = iy) {
    // 			if (x > 1) {
    // 				x--;
    // 				iy = dl[x];
    // 				dl[x] = ix;
    // 			} else { ... MDB_CORRUPTED }
    // 		}
    // 		txn->mt_dirty_room++;
    // 		if (!(env->me_flags & MDB_WRITEMAP))
    // 			mdb_dpage_free(env, mp);
    // release:
    // 		/* Insert in me_pghead */
    // 		mop = env->me_pghead;
    // 		j = mop[0] + ovpages;
    // 		for (i = mop[0]; i && mop[i] < pg; i--)
    // 			mop[j--] = mop[i];
    // 		while (j>i)
    // 			mop[j--] = pg++;
    // 		mop[0] += ovpages;
    // 	} else {
    // 		rc = mdb_midl_append_range(&txn->mt_free_pgs, pg, ovpages);
    // 		if (rc)
    // 			return rc;
    // 	}
    // 	mc->mc_db->md_overflow_pages -= ovpages;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:7169-7205 (mdb_page_new)
    // ```c
    // 	if ((rc = mdb_page_alloc(mc, num, &np)))
    // 		return rc;
    // 	np->mp_flags = flags | P_DIRTY;
    // 	np->mp_lower = (PAGEHDRSZ-PAGEBASE);
    // 	np->mp_upper = mc->mc_txn->mt_env->me_psize - PAGEBASE;
    // 	if (IS_BRANCH(np))
    // 		mc->mc_db->md_branch_pages++;
    // 	else if (IS_LEAF(np))
    // 		mc->mc_db->md_leaf_pages++;
    // 	else if (IS_OVERFLOW(np)) {
    // 		mc->mc_db->md_overflow_pages += num;
    // 		np->mp_pages = num;
    // 	}
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:7207-7230 (mdb_leaf_size)
    // ```c
    // 	sz = LEAFSIZE(key, data);
    // 	if (sz > env->me_nodemax) {
    // 		/* put on overflow page */
    // 		sz -= data->mv_size - sizeof(pgno_t);
    // 	}
    // 	return EVEN(sz + sizeof(indx_t));
    // ```
    fn leaf_size(&self, ksize: usize, dsize: usize) -> usize {
        let mut sz = NODESIZE + ksize + dsize;
        if sz > self.nodemax {
            sz -= dsize - 8;
        }
        even(sz + 2)
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:7232-7253 (mdb_branch_size)
    // ```c
    // 	sz = INDXSIZE(key);
    // 	...
    // 	return sz + sizeof(indx_t);
    // ```
    fn branch_size(ksize: usize) -> usize {
        NODESIZE + ksize + 2
    }

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:7255-7391 (mdb_node_add)
    // ```c
    // 	if (IS_LEAF2(mp)) {
    // 		/* Move higher keys up one slot. */
    // 		int ksize = mc->mc_db->md_pad, dif;
    // 		char *ptr = LEAF2KEY(mp, indx, ksize);
    // 		dif = NUMKEYS(mp) - indx;
    // 		if (dif > 0)
    // 			memmove(ptr+ksize, ptr, dif*ksize);
    // 		/* insert new key */
    // 		memcpy(ptr, key->mv_data, ksize);
    // 		/* Just using these for counting */
    // 		mp->mp_lower += sizeof(indx_t);
    // 		mp->mp_upper -= ksize - sizeof(indx_t);
    // 		return MDB_SUCCESS;
    // 	}
    // 	room = (ssize_t)SIZELEFT(mp) - (ssize_t)sizeof(indx_t);
    // 	if (key != NULL)
    // 		node_size += key->mv_size;
    // 	if (IS_LEAF(mp)) {
    // 		mdb_cassert(mc, key && data);
    // 		if (F_ISSET(flags, F_BIGDATA)) {
    // 			/* Data already on overflow page. */
    // 			node_size += sizeof(pgno_t);
    // 		} else if (node_size + data->mv_size > mc->mc_txn->mt_env->me_nodemax) {
    // 			int ovpages = OVPAGES(data->mv_size, mc->mc_txn->mt_env->me_psize);
    // 			int rc;
    // 			/* Put data on overflow page. */
    // 			node_size = EVEN(node_size + sizeof(pgno_t));
    // 			if ((ssize_t)node_size > room)
    // 				goto full;
    // 			if ((rc = mdb_page_new(mc, P_OVERFLOW, ovpages, &ofp)))
    // 				return rc;
    // 			flags |= F_BIGDATA;
    // 			goto update;
    // 		} else {
    // 			node_size += data->mv_size;
    // 		}
    // 	}
    // 	node_size = EVEN(node_size);
    // 	if ((ssize_t)node_size > room)
    // 		goto full;
    // update:
    // 	/* Move higher pointers up one slot. */
    // 	for (i = NUMKEYS(mp); i > indx; i--)
    // 		mp->mp_ptrs[i] = mp->mp_ptrs[i - 1];
    // 	/* Adjust free space offsets. */
    // 	ofs = mp->mp_upper - node_size;
    // 	mdb_cassert(mc, ofs >= mp->mp_lower + sizeof(indx_t));
    // 	mp->mp_ptrs[indx] = ofs;
    // 	mp->mp_upper = ofs;
    // 	mp->mp_lower += sizeof(indx_t);
    // 	/* Write the node data. */
    // 	node = NODEPTR(mp, indx);
    // 	node->mn_ksize = (key == NULL) ? 0 : key->mv_size;
    // 	node->mn_flags = flags;
    // 	if (IS_LEAF(mp))
    // 		SETDSZ(node,data->mv_size);
    // 	else
    // 		SETPGNO(node,pgno);
    // 	if (key)
    // 		memcpy(NODEKEY(node), key->mv_data, key->mv_size);
    // 	if (IS_LEAF(mp)) {
    // 		ndata = NODEDATA(node);
    // 		if (ofp == NULL) {
    // 			if (F_ISSET(flags, F_BIGDATA))
    // 				memcpy(ndata, data->mv_data, sizeof(pgno_t));
    // 			else if (F_ISSET(flags, MDB_RESERVE))
    // 				data->mv_data = ndata;
    // 			else
    // 				memcpy(ndata, data->mv_data, data->mv_size);
    // 		} else {
    // 			memcpy(ndata, &ofp->mp_pgno, sizeof(pgno_t));
    // 			ndata = METADATA(ofp);
    // 			if (F_ISSET(flags, MDB_RESERVE))
    // 				data->mv_data = ndata;
    // 			else
    // 				memcpy(ndata, data->mv_data, data->mv_size);
    // 		}
    // 	}
    // 	return MDB_SUCCESS;
    // full:
    // 	mc->mc_txn->mt_flags |= MDB_TXN_ERROR;
    // 	return MDB_PAGE_FULL;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:7393-7448 (mdb_node_del)
    // ```c
    // 	numkeys = NUMKEYS(mp);
    // 	if (IS_LEAF2(mp)) {
    // 		int x = numkeys - 1 - indx;
    // 		base = LEAF2KEY(mp, indx, ksize);
    // 		if (x)
    // 			memmove(base, base + ksize, x * ksize);
    // 		mp->mp_lower -= sizeof(indx_t);
    // 		mp->mp_upper += ksize - sizeof(indx_t);
    // 		return;
    // 	}
    // 	node = NODEPTR(mp, indx);
    // 	sz = NODESIZE + node->mn_ksize;
    // 	if (IS_LEAF(mp)) {
    // 		if (F_ISSET(node->mn_flags, F_BIGDATA))
    // 			sz += sizeof(pgno_t);
    // 		else
    // 			sz += NODEDSZ(node);
    // 	}
    // 	sz = EVEN(sz);
    // 	ptr = mp->mp_ptrs[indx];
    // 	for (i = j = 0; i < numkeys; i++) {
    // 		if (i != indx) {
    // 			mp->mp_ptrs[j] = mp->mp_ptrs[i];
    // 			if (mp->mp_ptrs[i] < ptr)
    // 				mp->mp_ptrs[j] += sz;
    // 			j++;
    // 		}
    // 	}
    // 	base = (char *)mp + mp->mp_upper + PAGEBASE;
    // 	memmove(base + sz, base, ptr - mp->mp_upper);
    // 	mp->mp_lower -= sizeof(indx_t);
    // 	mp->mp_upper += sz;
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:8582-9006 (mdb_page_split); the parts
    // are quoted where they are ported below.
    /// `mdb_page_split`; for `MDB_RESERVE`, where the data goes.
    fn page_split(
        &mut self,
        mc: usize,
        newkey: &[u8],
        newdata: Option<Val>,
        newpgno: u64,
        nflags: u32,
    ) -> (i32, Option<Loc>) {
        // ```c
        // 	mp = mc->mc_pg[mc->mc_top];
        // 	newindx = mc->mc_ki[mc->mc_top];
        // 	nkeys = NUMKEYS(mp);
        // 	/* Create a right sibling. */
        // 	if ((rc = mdb_page_new(mc, mp->mp_flags, 1, &rp)))
        // 		return rc;
        // 	rp->mp_pad = mp->mp_pad;
        // 	if (mc->mc_top < 1) {
        // 		if ((rc = mdb_page_new(mc, P_BRANCH, 1, &pp)))
        // 			goto done;
        // 		/* shift current top to make room for new parent */
        // 		for (i=mc->mc_snum; i>0; i--) {
        // 			mc->mc_pg[i] = mc->mc_pg[i-1];
        // 			mc->mc_ki[i] = mc->mc_ki[i-1];
        // 		}
        // 		mc->mc_pg[0] = pp;
        // 		mc->mc_ki[0] = 0;
        // 		mc->mc_db->md_root = pp->mp_pgno;
        // 		new_root = mc->mc_db->md_depth++;
        // 		/* Add left (implicit) pointer. */
        // 		if ((rc = mdb_node_add(mc, 0, NULL, NULL, mp->mp_pgno, 0)) != MDB_SUCCESS) {
        // 			/* undo the pre-push */
        // 			...
        // 			goto done;
        // 		}
        // 		mc->mc_snum++;
        // 		mc->mc_top++;
        // 		ptop = 0;
        // 	} else {
        // 		ptop = mc->mc_top-1;
        // 	}
        // 	mdb_cursor_copy(mc, &mn);
        // 	mn.mc_xcursor = NULL;
        // 	mn.mc_pg[mn.mc_top] = rp;
        // 	mn.mc_ki[ptop] = mc->mc_ki[ptop]+1;
        // ```
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
            // ```c
            // 	if (nflags & MDB_APPEND) {
            // 		mn.mc_ki[mn.mc_top] = 0;
            // 		sepkey = *newkey;
            // 		split_indx = newindx;
            // 		nkeys = 0;
            // 	} else {
            // 		split_indx = (nkeys+1) / 2;
            // 		if (IS_LEAF2(rp)) {
            // 			... (move half of the keys of a LEAF2 page)
            // 		} else {
            // 			int psize, nsize, k;
            // 			/* Maximum free space in an empty page */
            // 			pmax = env->me_psize - PAGEHDRSZ;
            // 			if (IS_LEAF(mp))
            // 				nsize = mdb_leaf_size(env, newkey, newdata);
            // 			else
            // 				nsize = mdb_branch_size(env, newkey);
            // 			nsize = EVEN(nsize);
            // 			/* grab a page to hold a temporary copy */
            // 			copy = mdb_page_malloc(mc->mc_txn, 1);
            // 			if (copy == NULL) {
            // 				rc = ENOMEM;
            // 				goto done;
            // 			}
            // 			copy->mp_pgno  = mp->mp_pgno;
            // 			copy->mp_flags = mp->mp_flags;
            // 			copy->mp_lower = (PAGEHDRSZ-PAGEBASE);
            // 			copy->mp_upper = env->me_psize - PAGEBASE;
            // 			/* prepare to insert */
            // 			for (i=0, j=0; i<nkeys; i++) {
            // 				if (i == newindx) {
            // 					copy->mp_ptrs[j++] = 0;
            // 				}
            // 				copy->mp_ptrs[j++] = mp->mp_ptrs[i];
            // 			}
            // 			if (nkeys < 32 || nsize > pmax/16 || newindx >= nkeys) {
            // 				/* Find split point */
            // 				psize = 0;
            // 				if (newindx <= split_indx || newindx >= nkeys) {
            // 					i = 0; j = 1;
            // 					k = newindx >= nkeys ? nkeys : split_indx+1+IS_LEAF(mp);
            // 				} else {
            // 					i = nkeys; j = -1;
            // 					k = split_indx-1;
            // 				}
            // 				for (; i!=k; i+=j) {
            // 					if (i == newindx) {
            // 						psize += nsize;
            // 						node = NULL;
            // 					} else {
            // 						node = (MDB_node *)((char *)mp + copy->mp_ptrs[i] + PAGEBASE);
            // 						psize += NODESIZE + NODEKSZ(node) + sizeof(indx_t);
            // 						if (IS_LEAF(mp)) {
            // 							if (F_ISSET(node->mn_flags, F_BIGDATA))
            // 								psize += sizeof(pgno_t);
            // 							else
            // 								psize += NODEDSZ(node);
            // 						}
            // 						psize = EVEN(psize);
            // 					}
            // 					if (psize > pmax || i == k-j) {
            // 						split_indx = i + (j<0);
            // 						break;
            // 					}
            // 				}
            // 			}
            // 			if (split_indx == newindx) {
            // 				sepkey.mv_size = newkey->mv_size;
            // 				sepkey.mv_data = newkey->mv_data;
            // 			} else {
            // 				node = (MDB_node *)((char *)mp + copy->mp_ptrs[split_indx] + PAGEBASE);
            // 				sepkey.mv_size = node->mn_ksize;
            // 				sepkey.mv_data = NODEKEY(node);
            // 			}
            // 		}
            // 	}
            // ```
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
            // ```c
            // 	/* Copy separator key to the parent.
            // 	 */
            // 	if (SIZELEFT(mn.mc_pg[ptop]) < mdb_branch_size(env, &sepkey)) {
            // 		int snum = mc->mc_snum;
            // 		mn.mc_snum--;
            // 		mn.mc_top--;
            // 		did_split = 1;
            // 		/* We want other splits to find mn when doing fixups */
            // 		WITH_CURSOR_TRACKING(mn,
            // 			rc = mdb_page_split(&mn, &sepkey, NULL, rp->mp_pgno, 0));
            // 		if (rc)
            // 			goto done;
            // 		/* root split? */
            // 		if (mc->mc_snum > snum) {
            // 			ptop++;
            // 		}
            // 		/* Right page might now have changed parent.
            // 		 * Check if left page also changed parent.
            // 		 */
            // 		if (mn.mc_pg[ptop] != mc->mc_pg[ptop] &&
            // 		    mc->mc_ki[ptop] >= NUMKEYS(mc->mc_pg[ptop])) {
            // 			for (i=0; i<ptop; i++) {
            // 				mc->mc_pg[i] = mn.mc_pg[i];
            // 				mc->mc_ki[i] = mn.mc_ki[i];
            // 			}
            // 			mc->mc_pg[ptop] = mn.mc_pg[ptop];
            // 			if (mn.mc_ki[ptop]) {
            // 				mc->mc_ki[ptop] = mn.mc_ki[ptop] - 1;
            // 			} else {
            // 				/* find right page's left sibling */
            // 				mc->mc_ki[ptop] = mn.mc_ki[ptop];
            // 				mdb_cursor_sibling(mc, 0);
            // 			}
            // 		}
            // 	} else {
            // 		mn.mc_top--;
            // 		rc = mdb_node_add(&mn, mn.mc_ki[ptop], &sepkey, NULL, rp->mp_pgno, 0);
            // 		mn.mc_top++;
            // 	}
            // 	if (rc != MDB_SUCCESS) {
            // 		goto done;
            // 	}
            // ```
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
                // ```c
                // 	if (nflags & MDB_APPEND) {
                // 		mc->mc_pg[mc->mc_top] = rp;
                // 		mc->mc_ki[mc->mc_top] = 0;
                // 		rc = mdb_node_add(mc, 0, newkey, newdata, newpgno, nflags);
                // 		if (rc)
                // 			goto done;
                // 		for (i=0; i<mc->mc_top; i++)
                // 			mc->mc_ki[i] = mn.mc_ki[i];
                // ```
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
                // ```c
                // 	} else if (!IS_LEAF2(mp)) {
                // 		/* Move nodes */
                // 		mc->mc_pg[mc->mc_top] = rp;
                // 		i = split_indx;
                // 		j = 0;
                // 		do {
                // 			if (i == newindx) {
                // 				rkey.mv_data = newkey->mv_data;
                // 				rkey.mv_size = newkey->mv_size;
                // 				if (IS_LEAF(mp)) {
                // 					rdata = newdata;
                // 				} else
                // 					pgno = newpgno;
                // 				flags = nflags;
                // 				/* Update index for the new key. */
                // 				mc->mc_ki[mc->mc_top] = j;
                // 			} else {
                // 				node = (MDB_node *)((char *)mp + copy->mp_ptrs[i] + PAGEBASE);
                // 				rkey.mv_data = NODEKEY(node);
                // 				rkey.mv_size = node->mn_ksize;
                // 				if (IS_LEAF(mp)) {
                // 					xdata.mv_data = NODEDATA(node);
                // 					xdata.mv_size = NODEDSZ(node);
                // 					rdata = &xdata;
                // 				} else
                // 					pgno = NODEPGNO(node);
                // 				flags = node->mn_flags;
                // 			}
                // 			if (!IS_LEAF(mp) && j == 0) {
                // 				/* First branch index doesn't need key data. */
                // 				rkey.mv_size = 0;
                // 			}
                // 			rc = mdb_node_add(mc, j, &rkey, rdata, pgno, flags);
                // 			if (rc)
                // 				goto done;
                // 			if (i == nkeys) {
                // 				i = 0;
                // 				j = 0;
                // 				mc->mc_pg[mc->mc_top] = copy;
                // 			} else {
                // 				i++;
                // 				j++;
                // 			}
                // 		} while (i != split_indx);
                // 		nkeys = NUMKEYS(copy);
                // 		for (i=0; i<nkeys; i++)
                // 			mp->mp_ptrs[i] = copy->mp_ptrs[i];
                // 		mp->mp_lower = copy->mp_lower;
                // 		mp->mp_upper = copy->mp_upper;
                // 		memcpy(NODEPTR(mp, nkeys-1), NODEPTR(copy, nkeys-1),
                // 			env->me_psize - copy->mp_upper - PAGEBASE);
                // 		/* reset back to original page */
                // 		if (newindx < split_indx) {
                // 			mc->mc_pg[mc->mc_top] = mp;
                // 		} else {
                // 			mc->mc_pg[mc->mc_top] = rp;
                // 			mc->mc_ki[ptop]++;
                // 			/* Make sure mc_ki is still valid.
                // 			 */
                // 			if (mn.mc_pg[ptop] != mc->mc_pg[ptop] &&
                // 				mc->mc_ki[ptop] >= NUMKEYS(mc->mc_pg[ptop])) {
                // 				for (i=0; i<=ptop; i++) {
                // 					mc->mc_pg[i] = mn.mc_pg[i];
                // 					mc->mc_ki[i] = mn.mc_ki[i];
                // 				}
                // 			}
                // 		}
                // 		if (nflags & MDB_RESERVE) {
                // 			node = NODEPTR(mc->mc_pg[mc->mc_top], mc->mc_ki[mc->mc_top]);
                // 			if (!(node->mn_flags & F_BIGDATA))
                // 				newdata->mv_data = NODEDATA(node);
                // 		}
                // ```
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
            // ```c
            // 	{
            // 		/* Adjust other cursors pointing to mp */
            // 		MDB_cursor *m2, *m3;
            // 		MDB_dbi dbi = mc->mc_dbi;
            // 		nkeys = NUMKEYS(mp);
            // 		for (m2 = mc->mc_txn->mt_cursors[dbi]; m2; m2=m2->mc_next) {
            // 			if (mc->mc_flags & C_SUB)
            // 				m3 = &m2->mc_xcursor->mx_cursor;
            // 			else
            // 				m3 = m2;
            // 			if (m3 == mc)
            // 				continue;
            // 			if (!(m2->mc_flags & m3->mc_flags & C_INITIALIZED))
            // 				continue;
            // 			if (new_root) {
            // 				int k;
            // 				/* sub cursors may be on different DB */
            // 				if (m3->mc_pg[0] != mp)
            // 					continue;
            // 				/* root split */
            // 				for (k=new_root; k>=0; k--) {
            // 					m3->mc_ki[k+1] = m3->mc_ki[k];
            // 					m3->mc_pg[k+1] = m3->mc_pg[k];
            // 				}
            // 				if (m3->mc_ki[0] >= nkeys) {
            // 					m3->mc_ki[0] = 1;
            // 				} else {
            // 					m3->mc_ki[0] = 0;
            // 				}
            // 				m3->mc_pg[0] = mc->mc_pg[0];
            // 				m3->mc_snum++;
            // 				m3->mc_top++;
            // 			}
            // 			if (m3->mc_top >= mc->mc_top && m3->mc_pg[mc->mc_top] == mp) {
            // 				if (m3->mc_ki[mc->mc_top] >= newindx && !(nflags & MDB_SPLIT_REPLACE))
            // 					m3->mc_ki[mc->mc_top]++;
            // 				if (m3->mc_ki[mc->mc_top] >= nkeys) {
            // 					m3->mc_pg[mc->mc_top] = rp;
            // 					m3->mc_ki[mc->mc_top] -= nkeys;
            // 					for (i=0; i<mc->mc_top; i++) {
            // 						m3->mc_ki[i] = mn.mc_ki[i];
            // 						m3->mc_pg[i] = mn.mc_pg[i];
            // 					}
            // 				}
            // 			} else if (!did_split && m3->mc_top >= ptop && m3->mc_pg[ptop] == mc->mc_pg[ptop] &&
            // 				m3->mc_ki[ptop] >= mc->mc_ki[ptop]) {
            // 				m3->mc_ki[ptop]++;
            // 			}
            // 			if (IS_LEAF(mp))
            // 				XCURSOR_REFRESH(m3, mc->mc_top, m3->mc_pg[mc->mc_top]);
            // 		}
            // 	}
            // ```
            // `new_root` is the depth before a root split (`new_root = mc->mc_db->md_depth++`),
            // so at least 1 when the root was split; `None` stands for the C value 0.
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
        // ```c
        // done:
        // 	if (copy)					/* tmp page */
        // 		mdb_page_free(env, copy);
        // 	if (rc)
        // 		mc->mc_txn->mt_flags |= MDB_TXN_ERROR;
        // 	return rc;
        // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:3073-3293 (mdb_freelist_save); the
    // parts are quoted where they are ported below. There are no loose pages, so the parts of
    // `mt_loose_pgs` (3098-3138, 3246-3261) never run.
    pub(super) fn freelist_save(&mut self) -> i32 {
        // ```c
        // 	MDB_cursor mc;
        // 	MDB_env	*env = txn->mt_env;
        // 	int rc, maxfree_1pg = env->me_maxfree_1pg, more = 1;
        // 	txnid_t	pglast = 0, head_id = 0;
        // 	pgno_t	freecnt = 0, *free_pgs, *mop;
        // 	ssize_t	head_room = 0, total_room = 0, mop_len, clean_limit;
        // 	mdb_cursor_init(&mc, txn, FREE_DBI, NULL);
        // 	if (env->me_pghead) {
        // 		/* Make sure first page of freeDB is touched and on freelist */
        // 		rc = mdb_page_search(&mc, NULL, MDB_PS_FIRST|MDB_PS_MODIFY);
        // 		if (rc && rc != MDB_NOTFOUND)
        // 			return rc;
        // 	}
        // 	...
        // 	/* MDB_RESERVE cancels meminit in ovpage malloc (when no WRITEMAP) */
        // 	clean_limit = (env->me_flags & (MDB_NOMEMINIT|MDB_WRITEMAP))
        // 		? SSIZE_MAX : maxfree_1pg;
        // ```
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
            // ```c
            // 		/* If using records from freeDB which we have not yet
            // 		 * deleted, delete them and any we reserved for me_pghead.
            // 		 */
            // 		while (pglast < env->me_pglast) {
            // 			rc = mdb_cursor_first(&mc, &key, NULL);
            // 			if (rc)
            // 				return rc;
            // 			pglast = head_id = *(txnid_t *)key.mv_data;
            // 			total_room = head_room = 0;
            // 			mdb_tassert(txn, pglast <= env->me_pglast);
            // 			rc = mdb_cursor_del(&mc, 0);
            // 			if (rc)
            // 				return rc;
            // 		}
            // ```
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
            // ```c
            // 		/* Save the IDL of pages freed by this txn, to a single record */
            // 		if (freecnt < txn->mt_free_pgs[0]) {
            // 			if (!freecnt) {
            // 				/* Make sure last page of freeDB is touched and on freelist */
            // 				rc = mdb_page_search(&mc, NULL, MDB_PS_LAST|MDB_PS_MODIFY);
            // 				if (rc && rc != MDB_NOTFOUND)
            // 					return rc;
            // 			}
            // 			free_pgs = txn->mt_free_pgs;
            // 			/* Write to last page of freeDB */
            // 			key.mv_size = sizeof(txn->mt_txnid);
            // 			key.mv_data = &txn->mt_txnid;
            // 			do {
            // 				freecnt = free_pgs[0];
            // 				data.mv_size = MDB_IDL_SIZEOF(free_pgs);
            // 				rc = mdb_cursor_put(&mc, &key, &data, MDB_RESERVE);
            // 				if (rc)
            // 					return rc;
            // 				/* Retry if mt_free_pgs[] grew during the Put() */
            // 				free_pgs = txn->mt_free_pgs;
            // 			} while (freecnt < free_pgs[0]);
            // 			mdb_midl_sort(free_pgs);
            // 			memcpy(data.mv_data, free_pgs, data.mv_size);
            // 			continue;
            // 		}
            // ```
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
            // ```c
            // 		mop = env->me_pghead;
            // 		mop_len = (mop ? mop[0] : 0) + txn->mt_loose_count;
            // 		/* Reserve records for me_pghead[]. Split it if multi-page,
            // 		 * to avoid searching freeDB for a page range. Use keys in
            // 		 * range [1,me_pglast]: Smaller than txnid of oldest reader.
            // 		 */
            // 		if (total_room >= mop_len) {
            // 			if (total_room == mop_len || --more < 0)
            // 				break;
            // 		} else if (head_room >= maxfree_1pg && head_id > 1) {
            // 			/* Keep current record (overflow page), add a new one */
            // 			head_id--;
            // 			head_room = 0;
            // 		}
            // 		/* (Re)write {key = head_id, IDL length = head_room} */
            // 		total_room -= head_room;
            // 		head_room = mop_len - total_room;
            // 		if (head_room > maxfree_1pg && head_id > 1) {
            // 			/* Overflow multi-page for part of me_pghead */
            // 			head_room /= head_id; /* amortize page sizes */
            // 			head_room += maxfree_1pg - head_room % (maxfree_1pg + 1);
            // 		} else if (head_room < 0) {
            // 			/* Rare case, not bothering to delete this record */
            // 			head_room = 0;
            // 		}
            // 		key.mv_size = sizeof(head_id);
            // 		key.mv_data = &head_id;
            // 		data.mv_size = (head_room + 1) * sizeof(pgno_t);
            // 		rc = mdb_cursor_put(&mc, &key, &data, MDB_RESERVE);
            // 		if (rc)
            // 			return rc;
            // 		/* IDL is initially empty, zero out at least the length */
            // 		pgs = (pgno_t *)data.mv_data;
            // 		j = head_room > clean_limit ? head_room : 0;
            // 		do {
            // 			pgs[j] = 0;
            // 		} while (--j >= 0);
            // 		total_room += head_room;
            // 	}
            // ```
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
        // ```c
        // 	/* Fill in the reserved me_pghead records */
        // 	rc = MDB_SUCCESS;
        // 	if (mop_len) {
        // 		MDB_val key, data;
        // 		mop += mop_len;
        // 		rc = mdb_cursor_first(&mc, &key, &data);
        // 		for (; !rc; rc = mdb_cursor_next(&mc, &key, &data, MDB_NEXT)) {
        // 			txnid_t id = *(txnid_t *)key.mv_data;
        // 			ssize_t	len = (ssize_t)(data.mv_size / sizeof(MDB_ID)) - 1;
        // 			MDB_ID save;
        // 			mdb_tassert(txn, len >= 0 && id <= env->me_pglast);
        // 			key.mv_data = &id;
        // 			if (len > mop_len) {
        // 				len = mop_len;
        // 				data.mv_size = (len + 1) * sizeof(MDB_ID);
        // 			}
        // 			data.mv_data = mop -= len;
        // 			save = mop[0];
        // 			mop[0] = len;
        // 			rc = mdb_cursor_put(&mc, &key, &data, MDB_CURRENT);
        // 			mop[0] = save;
        // 			if (rc || !(mop_len -= len))
        // 				break;
        // 		}
        // 	}
        // 	return rc;
        // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:9712-9728 (mdb_default_cmp)
    // ```c
    // 	uint16_t f = txn->mt_dbs[dbi].md_flags;
    // 	txn->mt_dbxs[dbi].md_cmp =
    // 		(f & MDB_REVERSEKEY) ? mdb_cmp_memnr :
    // 		(f & MDB_INTEGERKEY) ? mdb_cmp_cint  : mdb_cmp_memn;
    // 	txn->mt_dbxs[dbi].md_dcmp =
    // 		!(f & MDB_DUPSORT) ? 0 :
    // 		((f & MDB_INTEGERDUP)
    // 		 ? ((f & MDB_DUPFIXED)   ? mdb_cmp_int   : mdb_cmp_cint)
    // 		 : ((f & MDB_REVERSEDUP) ? mdb_cmp_memnr : mdb_cmp_memn));
    // ```
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:9730-9847 (mdb_dbi_open, named DB)
    // ```c
    // 	if (flags & ~VALID_FLAGS)
    // 		return EINVAL;
    // 	if (txn->mt_flags & MDB_TXN_BLOCKED)
    // 		return MDB_BAD_TXN;
    // 	...
    // 	if (txn->mt_dbxs[MAIN_DBI].md_cmp == NULL) {
    // 		mdb_default_cmp(txn, MAIN_DBI);
    // 	}
    // 	/* Is the DB already open? */
    // 	len = strlen(name);
    // 	for (i=CORE_DBS; i<txn->mt_numdbs; i++) {
    // 		if (!txn->mt_dbxs[i].md_name.mv_size) {
    // 			/* Remember this free slot */
    // 			if (!unused) unused = i;
    // 			continue;
    // 		}
    // 		if (len == txn->mt_dbxs[i].md_name.mv_size &&
    // 			!strncmp(name, txn->mt_dbxs[i].md_name.mv_data, len)) {
    // 			*dbi = i;
    // 			return MDB_SUCCESS;
    // 		}
    // 	}
    // 	/* If no free slot and max hit, fail */
    // 	if (!unused && txn->mt_numdbs >= txn->mt_env->me_maxdbs)
    // 		return MDB_DBS_FULL;
    // 	/* Cannot mix named databases with some mainDB flags */
    // 	if (txn->mt_dbs[MAIN_DBI].md_flags & (MDB_DUPSORT|MDB_INTEGERKEY))
    // 		return (flags & MDB_CREATE) ? MDB_INCOMPATIBLE : MDB_NOTFOUND;
    // 	/* Find the DB info */
    // 	dbflag = DB_NEW|DB_VALID|DB_USRVALID;
    // 	exact = 0;
    // 	key.mv_size = len;
    // 	key.mv_data = (void *)name;
    // 	mdb_cursor_init(&mc, txn, MAIN_DBI, NULL);
    // 	rc = mdb_cursor_set(&mc, &key, &data, MDB_SET, &exact);
    // 	if (rc == MDB_SUCCESS) {
    // 		/* make sure this is actually a DB */
    // 		MDB_node *node = NODEPTR(mc.mc_pg[mc.mc_top], mc.mc_ki[mc.mc_top]);
    // 		if ((node->mn_flags & (F_DUPDATA|F_SUBDATA)) != F_SUBDATA)
    // 			return MDB_INCOMPATIBLE;
    // 	} else {
    // 		if (rc != MDB_NOTFOUND || !(flags & MDB_CREATE))
    // 			return rc;
    // 		...
    // 	}
    // 	...
    // 	if (rc) {
    // 		/* MDB_NOTFOUND and MDB_CREATE: Create new DB */
    // 		data.mv_size = sizeof(MDB_db);
    // 		data.mv_data = &dummy;
    // 		memset(&dummy, 0, sizeof(dummy));
    // 		dummy.md_root = P_INVALID;
    // 		dummy.md_flags = flags & PERSISTENT_FLAGS;
    // 		WITH_CURSOR_TRACKING(mc,
    // 			rc = mdb_cursor_put(&mc, &key, &data, F_SUBDATA));
    // 		dbflag |= DB_DIRTY;
    // 	}
    // 	if (rc) {
    // 		free(namedup);
    // 	} else {
    // 		/* Got info, register DBI in this txn */
    // 		unsigned int slot = unused ? unused : txn->mt_numdbs;
    // 		txn->mt_dbxs[slot].md_name.mv_data = namedup;
    // 		txn->mt_dbxs[slot].md_name.mv_size = len;
    // 		txn->mt_dbxs[slot].md_rel = NULL;
    // 		txn->mt_dbflags[slot] = dbflag;
    // 		...
    // 		memcpy(&txn->mt_dbs[slot], data.mv_data, sizeof(MDB_db));
    // 		*dbi = slot;
    // 		mdb_default_cmp(txn, slot);
    // 		if (!unused) {
    // 			txn->mt_numdbs++;
    // 		}
    // 	}
    // ```
    /// `mdb_dbi_open(txn, name, flags)` (lmdb++ `dbi::open`).
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

    // NCBI reference (598d8ae6): c++/src/util/lmdb/mdb.c:9008-9031 (mdb_put)
    // ```c
    // 	if (!key || !data || !TXN_DBI_EXIST(txn, dbi, DB_USRVALID))
    // 		return EINVAL;
    // 	if (flags & ~(MDB_NOOVERWRITE|MDB_NODUPDATA|MDB_RESERVE|MDB_APPEND|MDB_APPENDDUP))
    // 		return EINVAL;
    // 	if (txn->mt_flags & (MDB_TXN_RDONLY|MDB_TXN_BLOCKED))
    // 		return (txn->mt_flags & MDB_TXN_RDONLY) ? EACCES : MDB_BAD_TXN;
    // 	mdb_cursor_init(&mc, txn, dbi, &mx);
    // 	mc.mc_next = txn->mt_cursors[dbi];
    // 	txn->mt_cursors[dbi] = &mc;
    // 	rc = mdb_cursor_put(&mc, key, data, flags);
    // 	txn->mt_cursors[dbi] = mc.mc_next;
    // 	return rc;
    // ```
    // NCBI reference (598d8ae6): c++/include/util/lmdbxx/lmdb++.h:850-861 (lmdb::dbi_put)
    // ```c++
    //   const int rc = ::mdb_put(txn, dbi, const_cast<MDB_val*>(key), const_cast<MDB_val*>(data), flags);
    //   if (rc != MDB_SUCCESS && rc != MDB_KEYEXIST) {
    //     error::raise("mdb_put", rc);
    //   }
    //   return (rc == MDB_SUCCESS);
    // ```
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
