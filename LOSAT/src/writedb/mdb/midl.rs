//! ID lists of LMDB (`c++/src/util/lmdb/midl.c`): an IDL is a list of page numbers whose
//! element 0 is the number of IDs, kept in descending order. Here an IDL is a `Vec<u64>` whose
//! length is at least `ids[0] + 1`; the C code's allocation sizes (`mdb_midl_need`,
//! `mdb_midl_grow`, `mdb_midl_shrink`) only decide when memory is reallocated, so `Vec` growth
//! stands in for them.

// NCBI reference (598d8ae6): c++/src/util/lmdb/midl.h:50-70
// ```c
// #define MDB_IDL_LOGN	16	/* DB_SIZE is 2^16, UM_SIZE is 2^17 */
// #define MDB_IDL_DB_SIZE		(1<<MDB_IDL_LOGN)
// #define MDB_IDL_UM_SIZE		(1<<(MDB_IDL_LOGN+1))
// #define MDB_IDL_UM_MAX		(MDB_IDL_UM_SIZE-1)
// #define MDB_IDL_SIZEOF(ids)		(((ids)[0]+1) * sizeof(MDB_ID))
// ```
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

// NCBI reference (598d8ae6): c++/src/util/lmdb/midl.c:31-61 (mdb_midl_search)
// ```c
// unsigned mdb_midl_search( MDB_IDL ids, MDB_ID id )
// {
//     unsigned base = 0; unsigned cursor = 1; int val = 0; unsigned n = ids[0];
//     while( 0 < n ) {
//         unsigned pivot = n >> 1;
//         cursor = base + pivot + 1;
//         val = CMP( ids[cursor], id );
//         if( val < 0 ) { n = pivot; }
//         else if ( val > 0 ) { base = cursor; n -= pivot + 1; }
//         else { return cursor; }
//     }
//     if( val > 0 ) { ++cursor; }
//     return cursor;
// }
// ```
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

/// `CMP(x,y)`: `( (x) < (y) ? -1 : (x) > (y) )`.
#[inline]
pub(super) fn cmp_u64(x: u64, y: u64) -> i32 {
    if x < y {
        -1
    } else {
        (x > y) as i32
    }
}

// NCBI reference (598d8ae6): c++/src/util/lmdb/midl.h:156-161 (mdb_midl_xappend)
// ```c
// #define mdb_midl_xappend(idl, id) do { \
//		MDB_ID *xidl = (idl), xlen = ++(xidl[0]); \
//		xidl[xlen] = (id); \
//	} while (0)
// ```
/// `mdb_midl_xappend` (the caller made room with `mdb_midl_need`; here the list grows).
pub(super) fn midl_xappend(ids: &mut Vec<u64>, id: u64) {
    room(ids, 1);
    ids[0] += 1;
    let n = ids[0] as usize;
    ids[n] = id;
}

// NCBI reference (598d8ae6): c++/src/util/lmdb/midl.c:142-155 (mdb_midl_append)
// ```c
// int mdb_midl_append( MDB_IDL *idp, MDB_ID id )
// {
//     MDB_IDL ids = *idp;
//     /* Too big? */
//     if (ids[0] >= ids[-1]) { if (mdb_midl_grow(idp, MDB_IDL_UM_MAX)) return ENOMEM; ids = *idp; }
//     ids[0]++;
//     ids[ids[0]] = id;
//     return 0;
// }
// ```
pub(super) fn midl_append(ids: &mut Vec<u64>, id: u64) {
    midl_xappend(ids, id);
}

// NCBI reference (598d8ae6): c++/src/util/lmdb/midl.c:172-187 (mdb_midl_append_range)
// ```c
// int mdb_midl_append_range( MDB_IDL *idp, MDB_ID id, unsigned n )
// {
//     MDB_ID *ids = *idp, len = ids[0];
//     /* Too big? */
//     if (len + n > ids[-1]) { if (mdb_midl_grow(idp, n | MDB_IDL_UM_MAX)) return ENOMEM; ids = *idp; }
//     ids[0] = len + n;
//     ids += len;
//     while (n)
//         ids[n--] = id++;
//     return 0;
// }
// ```
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

// NCBI reference (598d8ae6): c++/src/util/lmdb/midl.c:189-203 (mdb_midl_xmerge)
// ```c
// void mdb_midl_xmerge( MDB_IDL idl, MDB_IDL merge )
// {
//     MDB_ID old_id, merge_id, i = merge[0], j = idl[0], k = i+j, total = k;
//     idl[0] = (MDB_ID)-1;		/* delimiter for idl scan below */
//     old_id = idl[j];
//     while (i) {
//         merge_id = merge[i--];
//         for (; old_id < merge_id; old_id = idl[--j])
//             idl[k--] = old_id;
//         idl[k--] = merge_id;
//     }
//     idl[0] = total;
// }
// ```
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

// NCBI reference (598d8ae6): c++/src/util/lmdb/midl.c:205-264 (mdb_midl_sort)
// ```c
// /* Quicksort + Insertion sort for small arrays */
// void mdb_midl_sort( MDB_IDL ids )
// {
//     ...                 /* descending: "if (ids[i] >= a) break;" */
// }
// ```
/// `mdb_midl_sort`: sorts `ids[1..=ids[0]]` in descending order. The IDs of an IDL are
/// distinct page numbers, so any correct sort gives the order NCBI's quicksort gives.
pub(super) fn midl_sort(ids: &mut [u64]) {
    let n = ids[0] as usize;
    ids[1..=n].sort_unstable_by(|a, b| b.cmp(a));
}

/// `MDB_IDL_SIZEOF(ids)`: `(ids[0] + 1) * sizeof(MDB_ID)`.
pub(super) fn midl_sizeof(ids: &[u64]) -> usize {
    (ids[0] as usize + 1) * 8
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn xmerge_keeps_descending_order() {
        let mut a = vec![3, 9, 5, 2];
        midl_xmerge(&mut a, &[3, 8, 4, 1]);
        assert_eq!(&a[..7], &[6, 9, 8, 5, 4, 2, 1]);
        let mut e = midl_new();
        midl_xmerge(&mut e, &[2, 7, 3]);
        assert_eq!(&e[..3], &[2, 7, 3]);
    }

    #[test]
    fn append_range_is_descending() {
        let mut a = vec![1, 50];
        midl_append_range(&mut a, 10, 3);
        assert_eq!(&a[..5], &[4, 50, 12, 11, 10]);
    }

    #[test]
    fn search_finds_insert_position_in_descending_list() {
        let a = [4, 40, 30, 20, 10];
        assert_eq!(midl_search(&a, 30), 2);
        assert_eq!(midl_search(&a, 25), 3);
        assert_eq!(midl_search(&a, 50), 1);
        assert_eq!(midl_search(&a, 5), 5);
    }
}
