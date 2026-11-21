use super::*;
use std::cell::{Cell, RefCell};

/// `x_trimmed` on the one table `memo`.
impl Memo {
    fn trimmed(
        &mut self,
        segment: &[u8],
        maxtrim: usize,
        trim: impl FnOnce() -> (usize, usize),
    ) -> (usize, usize) {
        let memo = RefCell::new(std::mem::replace(self, Memo::new()));
        let pair = x_trimmed_in(|| memo.borrow_mut(), segment, maxtrim, trim);
        *self = memo.into_inner();
        pair
    }
}

#[test]
fn memo_returns_the_stored_pair_for_equal_bytes_and_maxtrim_only() {
    let mut memo = Memo::new();
    let calls = Cell::new(0usize);
    let trim = |pair: (usize, usize)| {
        let calls = &calls;
        move || {
            calls.set(calls.get() + 1);
            pair
        }
    };
    let a = [1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];
    let b = [1u8, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12];
    assert_eq!(memo.trimmed(&a, 50, trim((1, 9))), (1, 9));
    // A hit returns the stored pair without calling the trim.
    assert_eq!(memo.trimmed(&a, 50, trim((0, 0))), (1, 9));
    assert_eq!(calls.get(), 1);
    // One letter differs (in the tail word): evaluated.
    assert_eq!(memo.trimmed(&b, 50, trim((2, 8))), (2, 8));
    assert_eq!(calls.get(), 2);
    // Same letters, other maxtrim: evaluated.
    assert_eq!(memo.trimmed(&a, 49, trim((3, 7))), (3, 7));
    assert_eq!(calls.get(), 3);
    // A prefix is another segment.
    assert_eq!(memo.trimmed(&a[..10], 50, trim((0, 9))), (0, 9));
    assert_eq!(calls.get(), 4);
    // Too long for the memo: evaluated every time.
    let long = [7u8; X_MAX_SEGMENT + 1];
    assert_eq!(memo.trimmed(&long, 50, trim((5, 6))), (5, 6));
    assert_eq!(memo.trimmed(&long, 50, trim((5, 6))), (5, 6));
    assert_eq!(calls.get(), 6);
}

#[test]
fn trim_runs_without_the_table_and_a_stored_key_is_not_stored_twice() {
    let memo = RefCell::new(Memo::new());
    let a = [3u8, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5];
    // The trim takes the table itself (a held table would make this borrow panic) and stores
    // the same key first, as another thread that missed it at the same moment would.
    let pair = x_trimmed_in(
        || memo.borrow_mut(),
        &a,
        50,
        || {
            memo.borrow_mut().store(hash_of(&a, 50), &a, 50, (2, 7));
            (2, 7)
        },
    );
    assert_eq!(pair, (2, 7));
    assert_eq!(memo.borrow().entries, 1);
    assert_eq!(
        x_trimmed_in(|| memo.borrow_mut(), &a, 50, || unreachable!()),
        (2, 7)
    );
}

#[test]
fn memo_keeps_every_entry_through_growth_and_clears_at_its_entry_bound() {
    let mut memo = Memo::new();
    // Distinct segments: the index in base 251 over three letters, plus one more.
    let segment = |n: usize| {
        vec![
            (n % 251) as u8,
            (n / 251 % 251) as u8,
            (n / 251 / 251) as u8,
            (n % 7) as u8,
        ]
    };
    let pair = |n: usize| (n % 3, n % 5 + 3);
    let count = 3 * X_FIRST_SLOTS;
    for n in 0..count {
        assert_eq!(memo.trimmed(&segment(n), 50, || pair(n)), pair(n));
    }
    assert_eq!(memo.entries, count);
    assert!(memo.slots.len() >= 2 * count);
    // All of them are still there after the table doubled.
    for n in 0..count {
        assert_eq!(memo.trimmed(&segment(n), 50, || unreachable!()), pair(n));
    }
    // Filling past the bound clears the table once; the results stay the trim's.
    for n in count..X_MAX_ENTRIES + count + 10 {
        assert_eq!(memo.trimmed(&segment(n), 50, || pair(n)), pair(n));
    }
    assert_eq!(memo.entries, count + 10);
    assert!(memo.slots.len() <= 2 * X_MAX_ENTRIES);
    let mut evaluated = 0usize;
    for n in 0..count {
        assert_eq!(
            memo.trimmed(&segment(n), 50, || {
                evaluated += 1;
                pair(n)
            }),
            pair(n)
        );
    }
    assert_eq!(evaluated, count, "the first entries were cleared");
}

#[test]
fn memo_clears_at_its_letter_bound() {
    let mut memo = Memo::new();
    // Distinct segments of X_MAX_SEGMENT letters: the arena fills before the entry bound.
    let segment = |n: usize| {
        let mut s = vec![(n % 251) as u8; X_MAX_SEGMENT];
        s[0] = (n / 251) as u8;
        s
    };
    let fit = X_MAX_ARENA / X_MAX_SEGMENT;
    for n in 0..fit {
        memo.trimmed(&segment(n), 50, || (0, 1));
    }
    assert_eq!(memo.entries, fit);
    assert_eq!(memo.trimmed(&segment(0), 50, || unreachable!()), (0, 1));
    // One more segment does not fit in the arena: the table is cleared first.
    memo.trimmed(&segment(fit), 50, || (0, 2));
    assert_eq!(memo.entries, 1);
    assert_eq!(memo.arena.len(), X_MAX_SEGMENT);
    assert_eq!(memo.trimmed(&segment(fit), 50, || unreachable!()), (0, 2));
}
