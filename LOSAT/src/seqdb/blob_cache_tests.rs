//! The model of NCBI's cache of released blobs (`BlobCache`, `SeqDb::batch_blobs`), against the
//! batches of the probes of line RG-KB 9b-2 (`K/kb/i9b2/probe/probe3.py`, NCBI BLAST+ 2.17.0 under
//! several allocator settings: a subject that a partly emptied release holds or not changes with
//! the allocator, one of a whole release or of none does not).

use super::{BlobCache, BlobHeld, RELEASED_BLOBS};
use BlobHeld::{Held, Out, Unsure};

/// The batch of subjects `ordinals` (every title posts).
fn batch(cache: &mut BlobCache, ordinals: std::ops::Range<u32>) -> Vec<BlobHeld> {
    let subjects: Vec<(u32, bool)> = ordinals.map(|ordinal| (ordinal, true)).collect();
    cache.batch(&subjects, false)
}

#[test]
fn a_release_that_fits_is_held_whole() {
    let mut cache = BlobCache::default();
    assert_eq!(batch(&mut cache, 0..6), vec![Out; 6]);
    assert_eq!(batch(&mut cache, 0..6), vec![Held; 6]);
    // Held again, the subjects are the last release's alone.
    assert_eq!(batch(&mut cache, 6..10), vec![Out; 4]);
    assert_eq!(batch(&mut cache, 0..6), vec![Held; 6]);
}

#[test]
fn a_release_partly_pushed_out_is_unsure() {
    // probe3 `aba`: families A (6), B (5), A: the release of B pushes one of A out.
    let mut cache = BlobCache::default();
    batch(&mut cache, 0..6);
    assert_eq!(batch(&mut cache, 6..11), vec![Out; 5]);
    assert_eq!(batch(&mut cache, 0..6), vec![Unsure; 6]);
    // `bac` (B, A, C) and `acb` (A, C, B): the third batch's subjects were never loaded.
    let mut cache = BlobCache::default();
    batch(&mut cache, 6..11);
    batch(&mut cache, 0..6);
    assert_eq!(batch(&mut cache, 11..14), vec![Out; 3]);
    let mut cache = BlobCache::default();
    batch(&mut cache, 0..6);
    batch(&mut cache, 11..14);
    assert_eq!(batch(&mut cache, 6..11), vec![Out; 5]);
    // `manymb`: a batch of more than the cache keeps leaves some of itself.
    let mut cache = BlobCache::default();
    let big = u32::try_from(RELEASED_BLOBS).unwrap() + 12;
    batch(&mut cache, 0..big);
    assert_eq!(
        batch(&mut cache, big - 2..big + 2),
        vec![Unsure, Unsure, Out, Out]
    );
}

#[test]
fn a_release_pushed_out_whole_is_out() {
    let mut cache = BlobCache::default();
    batch(&mut cache, 0..6);
    batch(&mut cache, 6..10);
    // A (6) and B (4) fill the cache; C (10) pushes them all out.
    assert_eq!(batch(&mut cache, 10..20), vec![Out; 10]);
    assert_eq!(batch(&mut cache, 0..10), vec![Out; 10]);
    // A partly pushed out leaves at the release that pushes out what is left of it.
    let mut cache = BlobCache::default();
    batch(&mut cache, 0..6);
    batch(&mut cache, 6..11);
    batch(&mut cache, 20..24);
    assert_eq!(batch(&mut cache, 0..1), vec![Unsure]);
    let mut cache = BlobCache::default();
    batch(&mut cache, 0..6);
    batch(&mut cache, 6..11);
    batch(&mut cache, 20..25);
    assert_eq!(batch(&mut cache, 0..1), vec![Out]);
}

#[test]
fn batches_before_the_model_leave_any_title_that_does_not_post_unsure() {
    let mut cache = BlobCache::default();
    assert_eq!(
        cache.batch(&[(0, true), (1, false)], true),
        vec![Out, Unsure]
    );
    // Ten blobs released after them push out what the batches before held.
    let subjects: Vec<(u32, bool)> = (10..18).map(|ordinal| (ordinal, false)).collect();
    cache.batch(&subjects, true);
    assert_eq!(cache.batch(&[(2, false)], true), vec![Out]);
}
