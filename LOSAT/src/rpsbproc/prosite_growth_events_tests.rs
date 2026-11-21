use super::super::RECORD_BYTES;
use std::cell::RefCell;

thread_local! {
    static LOG: RefCell<Vec<String>> = const { RefCell::new(Vec::new()) };
    static PAUSED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Stops (`true`) or resumes recording: LOSAT's second run of a match (for the other bytes
/// past the text) is not an allocation of NCBI's binary.
pub(crate) fn pause(on: bool) {
    PAUSED.with(|paused| paused.set(on));
}

/// The stack of NCBI's `Match` (`std::stack` of a `std::deque`) makes its map of 8 pointers
/// (64 bytes) and its first node of 512 bytes when it is constructed.
pub(crate) fn match_starts() {
    push("m64".to_string());
    push("m512".to_string());
}

fn push(event: String) {
    if !PAUSED.with(std::cell::Cell::get) {
        LOG.with(|log| log.borrow_mut().push(event));
    }
}

pub(crate) fn take() -> Vec<String> {
    LOG.with(|log| std::mem::take(&mut *log.borrow_mut()))
}

pub(super) fn string_grows(_old: u64, new: u64) {
    push(format!("m{new}"));
}

pub(super) fn records_grow(_old: u64, new: u64) {
    push(format!("r{}", new * RECORD_BYTES));
}

pub(super) fn records_copy(len: u64) {
    push(format!("r{}", len * RECORD_BYTES));
}
