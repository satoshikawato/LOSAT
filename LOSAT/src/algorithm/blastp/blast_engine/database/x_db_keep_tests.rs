use super::*;

struct NoSink;
impl PrelimSink for NoSink {
    fn add(&mut self, _result: &mut BlastpSubjectResult) {}
}

fn window(sink: &mut NoSink, width: usize, floor: usize, max: usize) -> XWindow<'_, NoSink> {
    XWindow {
        next_ordinal: 0,
        diag_offset: 0,
        taken: 0,
        merged: 0,
        slots: std::iter::repeat_with(|| None).take(width).collect(),
        floor,
        waiting_bytes: 0,
        max_waiting_bytes: max,
        sink,
        error: None,
        stopped: false,
    }
}

// Past the width nothing is taken; past the floor only while few bytes wait; within the floor
// always, so the oldest block not merged always runs.
#[test]
fn the_window_narrows_to_the_floor_while_many_bytes_wait() {
    let mut sink = NoSink;
    let mut w = window(&mut sink, 64, 8, 1000);
    w.merged = 5;
    for (past, bytes, full) in [
        (0, 5000, false),
        (7, 5000, false),
        (8, 5000, true),
        (8, 1000, false),
        (63, 1000, false),
        (63, 1001, true),
        (64, 0, true),
    ] {
        w.taken = w.merged + past;
        w.waiting_bytes = bytes;
        assert_eq!(w.full(), full, "past {past} bytes {bytes}");
    }
}

#[test]
fn the_size_of_the_window_without_and_with_the_switch() {
    let narrow = XWindowSize::of(4, 8, false);
    assert_eq!((narrow.floor, narrow.width), (32, 32));
    assert_eq!(narrow.max_waiting_bytes, usize::MAX);
    let wide = XWindowSize::of(4, 8, true);
    assert_eq!((wide.floor, wide.width), (32, 256));
    assert_eq!(wide.max_waiting_bytes, 4 << 20);
}

#[test]
fn the_bytes_of_a_block_count_its_vectors() {
    use std::mem::size_of;
    let mut result = BlastpSubjectResult::new(0);
    let list = BlastpHspList {
        oid: 0,
        query_index: 0,
        hsps: Vec::with_capacity(10),
        best_evalue: 0.0,
    };
    let hsps = list.hsps.capacity() * size_of::<XPrelimHsp>();
    result.x_compact_hsp_lists = Vec::with_capacity(2);
    let lists =
        result.x_compact_hsp_lists.capacity() * size_of::<(usize, BlastpHspList<XPrelimHsp>)>();
    result.x_compact_hsp_lists.push((0, list));
    let block = vec![result];
    let blocks = block.capacity() * size_of::<BlastpSubjectResult>();
    assert_eq!(x_kept_bytes(&Ok((block, false))), hsps + lists + blocks);
    assert_eq!(x_kept_bytes(&Err(anyhow::anyhow!("failed"))), 0);
}
