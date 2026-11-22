use super::*;

#[test]
fn capacity_is_mapped_from_the_bound_on_glibc_only() {
    assert_eq!(x_mapped_capacity(X_MAPPED_FROM - 1, 1), X_MAPPED_FROM - 1);
    let at = x_mapped_capacity(X_MAPPED_FROM, 1);
    let rows = x_mapped_capacity(X_MAPPED_FROM / 12 + 1, 12);
    if X_MAPPED_ALLOCATOR {
        assert_eq!(at, X_MAPPED_BYTES);
        assert!(rows * 12 >= X_MAPPED_BYTES && (rows - 1) * 12 < X_MAPPED_BYTES);
    } else {
        assert_eq!(at, X_MAPPED_FROM);
        assert_eq!(rows, X_MAPPED_FROM / 12 + 1);
    }
    assert_eq!(x_mapped_capacity(40 << 20, 1), 40 << 20);
}

#[test]
fn reserve_doubles_as_vec_below_the_bound() {
    let mut v: Vec<u8> = Vec::with_capacity(1000);
    v.resize(1000, 7);
    x_reserve_mapped(&mut v, 1);
    assert_eq!(v.capacity(), 2000);
    assert_eq!(v.len(), 1000);
    x_reserve_mapped(&mut v, 10);
    assert_eq!(v.capacity(), 2000);
    x_reserve_mapped(&mut v, X_MAPPED_FROM);
    assert!(v.capacity() >= 1000 + X_MAPPED_FROM);
    assert!(v.iter().all(|&b| b == 7));
}

#[test]
fn sort_by_mapped_is_stable_and_equals_sort_by() {
    // (key, position): sorted by key only, so equal keys must keep their order.
    let mut state = 0x2545_f491_4f6c_dd1du64;
    for &n in &[0usize, 1, 5, 7, 9, 70_000, 200_003] {
        let v: Vec<(u32, u32)> = (0..n as u32)
            .map(|i| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                ((state % 1000) as u32, i)
            })
            .collect();
        let mut a = v.clone();
        a.sort_by(|x, y| y.0.cmp(&x.0));
        let mut b = v.clone();
        static POOL: XPool<(u32, u32)> = XPool::new(|| 64 << 20);
        x_sort_by_mapped(&mut b, |x, y| y.0.cmp(&x.0), &POOL);
        assert_eq!(a, b, "n = {n}");
        // The merges alone, with runs of 4 (every merge path, odd tails).
        let mut c = v.clone();
        for run in c.chunks_mut(4) {
            run.sort_by(|x, y| y.0.cmp(&x.0));
        }
        x_merge_runs(
            &mut c,
            &mut |x: &(u32, u32), y: &(u32, u32)| y.0.cmp(&x.0),
            4,
            &POOL,
        );
        assert_eq!(a, c, "n = {n}, runs of 4");
    }
}

#[test]
fn pool_lends_idle_buffers_and_keeps_small_ones_out() {
    let pool: XPool<u8> = XPool::new(|| 64 << 20);
    let mut small = vec![1u8; 100];
    pool.give(&mut small);
    assert_eq!(small.len(), 100);
    assert!(pool.take(1).is_none());
    let mut v: Vec<u8> = vec![5; 10];
    pool.reserve(&mut v, X_MAPPED_FROM);
    assert_eq!(&v[..], &[5u8; 10][..]);
    assert!(v.capacity() >= 10 + X_MAPPED_FROM);
    let addr = v.as_ptr();
    v.resize(X_MAPPED_FROM + 10, 0);
    pool.give(&mut v);
    // Elsewhere than glibc `give` keeps nothing and leaves `v` as it is.
    assert_eq!(v.capacity() == 0, X_MAPPED_ALLOCATOR);
    if X_MAPPED_ALLOCATOR {
        let mut w: Vec<u8> = vec![9; 3];
        pool.reserve(&mut w, X_MAPPED_FROM);
        assert_eq!(w.as_ptr(), addr);
        assert_eq!(&w[..], &[9u8; 3][..]);
    }
    // Only a buffer that holds what is asked is lent (a smaller one would be moved).
    let mut big: Vec<u8> = Vec::new();
    x_reserve_mapped(&mut big, X_MAPPED_FROM);
    let cap = big.capacity();
    pool.give(&mut big);
    assert!(pool.take(cap + 1).is_none());
    if X_MAPPED_ALLOCATOR {
        assert_eq!(pool.take(cap).map(|v| v.capacity()), Some(cap));
    }
}

#[test]
fn pool_frees_beyond_its_bound() {
    // A buffer that alone touched more than the bound is not kept.
    let none: XPool<u8> = XPool::new(|| 0);
    let mut v: Vec<u8> = Vec::new();
    x_reserve_mapped(&mut v, X_MAPPED_FROM);
    v.resize(X_MAPPED_FROM, 1);
    none.give(&mut v);
    assert_eq!(v.capacity() == 0, X_MAPPED_ALLOCATOR);
    assert!(none.take(1).is_none());
    // With room for one touched buffer, the older one is freed when a second comes back.
    let one: XPool<u8> = XPool::new(|| X_MAPPED_FROM + X_MAPPED_FROM / 2);
    let (mut a, mut b): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
    x_reserve_mapped(&mut a, X_MAPPED_FROM);
    x_reserve_mapped(&mut b, X_MAPPED_FROM);
    a.resize(X_MAPPED_FROM, 1);
    b.resize(X_MAPPED_FROM, 2);
    let b_addr = b.as_ptr();
    one.give(&mut a);
    one.give(&mut b);
    if X_MAPPED_ALLOCATOR {
        let kept = one.take(1).expect("the newer buffer is kept");
        assert_eq!(kept.as_ptr(), b_addr);
        assert!(one.take(1).is_none());
    } else {
        assert!(one.take(1).is_none());
    }
}

// `x_filled` (the database's block reads, C4-2 I): a buffer of the pool keeps the bytes its
// last user left and only the bytes past them are filled; `take` still lends an empty buffer; a
// buffer below `X_MAPPED_FROM` is a new one, filled.
#[test]
fn pool_fills_only_past_the_kept_bytes() {
    let pool: XPool<u8> = XPool::new(|| 64 << 20);
    let mut v: Vec<u8> = Vec::new();
    x_reserve_mapped(&mut v, 3 * X_MAPPED_FROM);
    v.resize(2 * X_MAPPED_FROM, 7);
    let addr = v.as_ptr();
    pool.give(&mut v);
    let longer = pool.x_filled(2 * X_MAPPED_FROM + 5, 0);
    assert_eq!(longer.len(), 2 * X_MAPPED_FROM + 5);
    if X_MAPPED_ALLOCATOR {
        assert_eq!(longer.as_ptr(), addr);
        assert!(longer[..2 * X_MAPPED_FROM].iter().all(|&b| b == 7));
        assert!(longer[2 * X_MAPPED_FROM..].iter().all(|&b| b == 0));
    } else {
        assert!(longer.iter().all(|&b| b == 0));
    }
    let mut longer = longer;
    pool.give(&mut longer);
    let shorter = pool.x_filled(X_MAPPED_FROM, 0);
    assert_eq!(shorter.len(), X_MAPPED_FROM);
    let kept = if X_MAPPED_ALLOCATOR { 7 } else { 0 };
    assert!(shorter.iter().all(|&b| b == kept));
    let mut shorter = shorter;
    pool.give(&mut shorter);
    let empty = pool.take(1);
    assert_eq!(empty.is_some(), X_MAPPED_ALLOCATOR);
    assert!(empty.is_none_or(|empty| empty.is_empty()));
    assert_eq!(pool.x_filled(10, 3), vec![3u8; 10]);
}

// `x_filled_kept` and `x_keep` (a worker's database blocks, C4-3): below `X_MAPPED_FROM` the
// worker's buffer is read again, filled only past its length, grown to exactly the length
// asked; from `X_MAPPED_FROM` on the pool's buffer, which goes back to the pool, and the
// worker's buffer stays.
#[test]
fn a_kept_buffer_is_filled_only_past_its_bytes_and_outliers_go_to_the_pool() {
    let pool: XPool<u8> = XPool::new(|| 64 << 20);
    let mut kept: Vec<u8> = Vec::new();
    let mut v = pool.x_filled_kept(&mut kept, 100, 0);
    assert_eq!((v.len(), v.capacity()), (100, 100));
    assert!(v.iter().all(|&b| b == 0));
    v.fill(7);
    pool.x_keep(&mut kept, v);
    assert_eq!(kept.len(), 100);
    let mut v = pool.x_filled_kept(&mut kept, 150, 0);
    assert!(kept.is_empty() && kept.capacity() == 0);
    assert_eq!((v.len(), v.capacity()), (150, 150));
    assert!(v[..100].iter().all(|&b| b == 7) && v[100..].iter().all(|&b| b == 0));
    v.fill(5);
    pool.x_keep(&mut kept, v);
    let v = pool.x_filled_kept(&mut kept, 40, 0);
    assert_eq!((v.len(), v.capacity()), (40, 150));
    assert!(v.iter().all(|&b| b == 5));
    pool.x_keep(&mut kept, v);
    let small = kept.as_ptr();
    // One long sequence: the pool's buffer; the worker's stays.
    let mut long = pool.x_filled_kept(&mut kept, X_MAPPED_FROM, 0);
    assert_eq!((kept.as_ptr(), kept.len()), (small, 40));
    assert_eq!(long.len(), X_MAPPED_FROM);
    assert!(long.iter().all(|&b| b == 0));
    long.fill(9);
    let addr = long.as_ptr();
    pool.x_keep(&mut kept, long);
    assert_eq!((kept.as_ptr(), kept.len()), (small, 40));
    let again = pool.x_filled_kept(&mut kept, X_MAPPED_FROM + 3, 0);
    assert_eq!(kept.as_ptr(), small);
    if X_MAPPED_ALLOCATOR {
        assert_eq!(again.as_ptr(), addr);
        assert!(again[..X_MAPPED_FROM].iter().all(|&b| b == 9));
        assert!(again[X_MAPPED_FROM..].iter().all(|&b| b == 0));
    } else {
        assert!(again.iter().all(|&b| b == 0));
    }
    let mut again = again;
    pool.give(&mut again);
}

// `x_take_near`, `x_disown` and `x_clear_idle`.
#[test]
fn pool_lends_near_buffers_forgets_kept_ones_and_clears() {
    let pool: XPool<u8> = XPool::new(|| 64 << 20);
    let touched = |pool: &XPool<u8>, addr: usize| {
        let st = pool.state.lock().unwrap();
        st.lent.iter().find(|&&(a, _)| a == addr).map(|&(_, t)| t)
    };
    let mut v: Vec<u8> = Vec::new();
    x_reserve_mapped(&mut v, 4 * X_MAPPED_FROM);
    v.resize(4 * X_MAPPED_FROM, 1);
    pool.give(&mut v);
    // A buffer that touched more than twice what is asked is not lent near.
    assert!(pool.x_take_near(X_MAPPED_FROM).is_none());
    let near = pool.x_take_near(2 * X_MAPPED_FROM);
    assert_eq!(near.is_some(), X_MAPPED_ALLOCATOR);
    if let Some(near) = near {
        let addr = near.as_ptr() as usize;
        assert_eq!(touched(&pool, addr), Some(4 * X_MAPPED_FROM));
        // Kept by its user: the pool forgets it.
        pool.x_disown(&near);
        assert_eq!(touched(&pool, addr), None);
    }
    let mut w: Vec<u8> = Vec::new();
    x_reserve_mapped(&mut w, X_MAPPED_FROM);
    w.resize(X_MAPPED_FROM, 2);
    pool.give(&mut w);
    pool.x_clear_idle();
    assert!(pool.take(1).is_none());
    assert_eq!(pool.state.lock().unwrap().idle_touched, 0);
}
