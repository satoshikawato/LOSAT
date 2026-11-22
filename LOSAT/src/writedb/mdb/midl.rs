// Derived from midl.c and midl.h of LMDB 0.9.24 as bundled in NCBI BLAST+ 2.17.0,
// c++/src/util/lmdb/, which carry this notice:
//
// This work is part of OpenLDAP Software <http://www.openldap.org/>.
//
// Copyright 2000-2019 The OpenLDAP Foundation.
// Portions Copyright 2001-2018 Howard Chu, Symas Corp.
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
// (In LOSAT, the OpenLDAP Public License 2.8 is in THIRD_PARTY_NOTICES, section 2 (LMDB).)

//! ID lists of LMDB (`c++/src/util/lmdb/midl.c`): an IDL is a list of page numbers whose
//! element 0 is the number of IDs, kept in descending order. Here an IDL is a `Vec<u64>` whose
//! length is at least `ids[0] + 1`; the C code's allocation sizes (`mdb_midl_need`,
//! `mdb_midl_grow`, `mdb_midl_shrink`) only decide when memory is reallocated, so `Vec` growth
//! stands in for them.

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/midl.h:50-70
/// `MDB_IDL_UM_MAX`.
pub(super) const MDB_IDL_UM_MAX: usize = (1 << 17) - 1;

/// A new IDL with no IDs.
pub(super) fn midl_new() -> Vec<u64> {
    vec![0]
}

/// Make room for `ids[0] + 1 + extra` elements.
fn room(ids: &mut Vec<u64>, extra: usize) {
    let need = ids[0] as usize + 1 + extra;
    if ids.len() < need {
        ids.resize(need, 0);
    }
}

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/midl.c:31-61 (mdb_midl_search)
#[allow(dead_code)]
pub(super) fn midl_search(ids: &[u64], id: u64) -> usize {
    let mut base = 0usize;
    let mut cursor = 1usize;
    let mut val = 0i32;
    let mut n = ids[0] as usize;
    while n > 0 {
        let pivot = n >> 1;
        cursor = base + pivot + 1;
        val = cmp_u64(ids[cursor], id);
        if val < 0 {
            n = pivot;
        } else if val > 0 {
            base = cursor;
            n -= pivot + 1;
        } else {
            return cursor;
        }
    }
    if val > 0 {
        cursor += 1;
    }
    cursor
}

/// `CMP` of `x` and `y`: -1 when `x` is less than `y`, 1 when it is greater, else 0.
#[inline]
pub(super) fn cmp_u64(x: u64, y: u64) -> i32 {
    if x < y {
        -1
    } else {
        (x > y) as i32
    }
}

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/midl.h:76-79 (mdb_midl_xappend)
/// `mdb_midl_xappend` (the caller made room with `mdb_midl_need`; here the list grows).
pub(super) fn midl_xappend(ids: &mut Vec<u64>, id: u64) {
    room(ids, 1);
    ids[0] += 1;
    let n = ids[0] as usize;
    ids[n] = id;
}

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/midl.c:157-169 (mdb_midl_append)
pub(super) fn midl_append(ids: &mut Vec<u64>, id: u64) {
    midl_xappend(ids, id);
}

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/midl.c:172-187 (mdb_midl_append_range)
pub(super) fn midl_append_range(ids: &mut Vec<u64>, mut id: u64, n: usize) {
    room(ids, n);
    let len = ids[0] as usize;
    ids[0] = (len + n) as u64;
    let mut k = n;
    while k > 0 {
        ids[len + k] = id;
        k -= 1;
        id += 1;
    }
}

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/midl.c:201-213 (mdb_midl_xmerge)
/// `mdb_midl_xmerge` (the caller made room with `mdb_midl_need`; here the list grows).
pub(super) fn midl_xmerge(idl: &mut Vec<u64>, merge: &[u64]) {
    let mut i = merge[0] as usize;
    room(idl, i);
    let mut j = idl[0] as usize;
    let mut k = i + j;
    let total = k;
    idl[0] = u64::MAX;
    let mut old_id = idl[j];
    while i > 0 {
        let merge_id = merge[i];
        i -= 1;
        while old_id < merge_id {
            idl[k] = old_id;
            k -= 1;
            j -= 1;
            old_id = idl[j];
        }
        idl[k] = merge_id;
        k = k.wrapping_sub(1);
    }
    idl[0] = total as u64;
}

// LMDB 0.9.24 (NCBI 598d8ae6 copy): c++/src/util/lmdb/midl.c:205-264 (mdb_midl_sort)
/// `mdb_midl_sort`: sorts `ids[1..=ids[0]]` in descending order. The IDs of an IDL are
/// distinct page numbers, so any correct sort gives the order NCBI's quicksort gives.
pub(super) fn midl_sort(ids: &mut [u64]) {
    let n = ids[0] as usize;
    ids[1..=n].sort_unstable_by(|a, b| b.cmp(a));
}

/// `MDB_IDL_SIZEOF` of `ids`: the `sizeof` of `MDB_ID` times one more than `ids[0]`.
pub(super) fn midl_sizeof(ids: &[u64]) -> usize {
    (ids[0] as usize + 1) * 8
}

#[cfg(test)]
#[path = "midl_tests.rs"]
mod tests;
