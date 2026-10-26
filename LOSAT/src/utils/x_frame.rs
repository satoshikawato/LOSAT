//! No NCBI counterpart (memory placement only): the stack frames of the 256-bit kernels.

/// No NCBI counterpart (memory placement only): called first in a 256-bit kernel, it makes the
/// compiler realign that kernel's stack frame to 64 bytes. With the ABI's 16-byte frame, LLVM
/// spills 256-bit registers with unaligned moves at fixed offsets from `rsp`, so whether a spill
/// crosses a cache line (a split load or store, slower) depends on `rsp mod 64`: on the callers'
/// frames, and for the main thread on the stack offset that address-space randomisation sets per
/// run. A change anywhere up the call chain moved blastp's X-drop DP by +-2 % with the same
/// instructions. Inlined, this function's 64-byte-aligned local is a local of the kernel, so the
/// kernel's frame is realigned and its spill slots lie on 32-byte boundaries; the volatile store
/// of one byte (never read) keeps the local. It does not change any value NCBI computes.
///
/// The realigned frame takes `rbp` as its frame pointer, one general register less in the kernel.
/// In align8's traceback kernel that meant stack loads and stores of general registers in the row
/// code (1.65 % more instructions on dc-megablast, 1.003 of the cycles at one thread), and in
/// megablast's greedy rows and blastp's restricted rows it cost more time than the splits saved,
/// so those kernels keep the ABI's frame.
#[cfg(target_arch = "x86_64")]
#[inline(always)]
pub(crate) fn x_frame_align() {
    #[repr(align(64))]
    struct XFrameAlign(u8);
    let mut x_frame = XFrameAlign(0);
    // SAFETY: `&mut x_frame.0` is a valid, aligned place of a local.
    unsafe { std::ptr::write_volatile(&mut x_frame.0, 0) };
}
