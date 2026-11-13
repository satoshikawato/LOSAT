// The shared memory of the threaded module (docs/web/abi_v2.md §3).
//
// The module declares the largest memory it can use: 65,535 pages, 4 GiB less one page
// (16,384 pages, 1 GiB, in the modules built before plan DW-28). An engine reserves that
// maximum as address space when the memory is created. Where it cannot, JavaScriptCore
// throws RangeError, and V8 without guard regions (on Android) keeps a smaller reservation
// without an error, at worst the initial pages alone, so that the memory cannot grow. A
// memory may have a smaller maximum than the module declares, so the host tries smaller
// ones in turn, down to the 1 GiB that the modules declared before. Below that the
// threaded engine does not start, as before, and the search runs on the serial module
// (plan §4.7): its memory is not reserved in advance and grows as far as the device lets it.
import type { MemoryLimits } from './artifact';

const PAGE_BYTES = 65_536;

/** The maximums tried after the declared one, in pages: 3 GiB, 2 GiB, 1.5 GiB, 1 GiB. */
export const STEP_DOWN_PAGES: readonly number[] = Object.freeze([49_152, 32_768, 24_576, 16_384]);

export interface SharedMemory {
  readonly memory: WebAssembly.Memory;
  /** The maximum of `memory` in pages: the declared one or a smaller one. */
  readonly maximumPages: number;
}

/** The maximums to try for a module's limits, largest first. */
export function candidateMaximums(limits: MemoryLimits): number[] {
  const declared = limits.maximum;
  if (declared === undefined) throw new Error('the threaded engine module declares no maximum for its shared memory');
  return [declared, ...STEP_DOWN_PAGES.filter((pages) => pages < declared && pages >= limits.initial)];
}

/**
 * Creates the shared memory with the largest maximum that this browser provides. Each
 * candidate is created one page short of the initial size and then grown by that page: a
 * memory that was created but cannot grow is not used.
 */
export function createSharedMemory(
  limits: MemoryLimits,
  create: (descriptor: WebAssembly.MemoryDescriptor) => WebAssembly.Memory = (descriptor) => new WebAssembly.Memory(descriptor),
): SharedMemory {
  const candidates = candidateMaximums(limits);
  let last: unknown;
  for (const maximumPages of candidates) {
    try {
      const memory = create({ initial: Math.max(0, limits.initial - 1), maximum: maximumPages, shared: true });
      memory.grow(1);
      return { memory, maximumPages };
    } catch (error) {
      last = error;
    }
  }
  const smallest = (candidates[candidates.length - 1]! * PAGE_BYTES) / 2 ** 30;
  const reason = last instanceof Error ? last.message : String(last);
  throw new Error(`the browser could not reserve ${formatGib(smallest)} GiB of shared memory (${reason})`, { cause: last });
}

function formatGib(value: number): string {
  return Number.isInteger(value) ? String(value) : value.toFixed(2);
}

/**
 * The cells of the waits of the threaded module (the shared-memory guard,
 * LOSAT/tests/wasi_shared_memory.js: WAITS_ENDING, WAIT_TABLE, WAIT_SLOTS). A wait records its
 * address plus 1 in the table while it waits (0 is a free cell); once WAITS_ENDING is set, waits wait in slices with a
 * call of the host between them. WebKit (JavaScriptCore, webkit-2311) never ends a worker
 * whose thread waits without a timeout, and keeps its instance and the shared memory until
 * the page closes: a page can reserve about 15 memories of 4 GiB, so the cancels of threaded
 * searches used them up (losat-baselines/web-flaky-20261103/cancel-timeout/ct3).
 */
export const WAITS_ENDING_ADDRESS = 20;
export const WAIT_TABLE_ADDRESS = 64;
export const WAIT_TABLE_SLOTS = 128;
/** How long endWaits goes on without any change in the table or any thread woken. */
const END_WAITS_MS = 1_000;

/**
 * Lets the workers of the threaded module's shared memory end: the memory keeper (and
 * ThreadHost.terminate) call it before the workers whose threads may be in the module end. It sets
 * WAITS_ENDING, then wakes every wait recorded in the table, again each millisecond, until the
 * table is empty: a woken thread waits again, now in slices, and its worker ends in the host call
 * between two of them. A thread that took its cell but has not begun its wait is woken by a later
 * round. It stops early only when nothing changed in the table and no thread was woken for
 * `deadlineMs` (a cell that stays set without a waiter, after an overflow of the main thread's
 * stack into the cells). The web build serves the threaded module only with the guard's waits
 * (web/app/build/reactors.ts); a memory whose buffer is not shared or too small for the cells is left as
 * it is.
 */
export function endWaits(memory: WebAssembly.Memory, deadlineMs = END_WAITS_MS): void {
  let until = performance.now() + deadlineMs;
  let last = '';
  const pause = new Int32Array(new SharedArrayBuffer(4));
  for (;;) {
    const buffer = memory.buffer as ArrayBufferLike;
    if (!(buffer instanceof SharedArrayBuffer) || buffer.byteLength < WAIT_TABLE_ADDRESS + 4 * WAIT_TABLE_SLOTS) return;
    // A view of the current size: the memory may have grown since the last round.
    const view = new Int32Array(buffer);
    Atomics.store(view, WAITS_ENDING_ADDRESS / 4, 1);
    const cells: number[] = [];
    let woken = 0;
    for (let slot = 0; slot < WAIT_TABLE_SLOTS; slot++) {
      // The address of the wait, plus 1 (0 is a free cell).
      const cell = Atomics.load(view, WAIT_TABLE_ADDRESS / 4 + slot);
      if (cell === 0) continue;
      cells.push(slot, cell);
      const index = ((cell >>> 0) - 1) >>> 2;
      if (index < view.length) woken += Atomics.notify(view, index);
    }
    if (cells.length === 0) return;
    const now = performance.now();
    const table = cells.join(',');
    if (woken > 0 || table !== last) until = now + deadlineMs;
    else if (now >= until) return;
    last = table;
    Atomics.wait(pause, 0, 0, 1);
  }
}

/**
 * Grows a shared memory to its maximum, so that no thread can grow it any more: the keeper
 * of the memory does so when the Gateway asks, before the Gateway ends the Engine worker and
 * its thread workers while the threaded instance is in use (memory-keeper.ts,
 * engine-worker/instance-use.ts). WebKit
 * (JavaScriptCore, webkit-2311) crashes the page when a thread grows a shared memory while
 * the instances of workers that end are taken off it (`SharedArrayBufferContents::grow`
 * reads a null instance); a memory at its maximum is not grown, so that cannot happen
 * (losat-baselines/web-flaky-20261103/crash). The threads may grow it meanwhile, and a grow
 * past the maximum throws, so the size is read again; a browser that cannot give the
 * maximum leaves the memory as it is.
 */
export function growToMaximum({ memory, maximumPages }: SharedMemory): void {
  for (let attempt = 0; attempt < 8; attempt++) {
    const pages = memory.buffer.byteLength / PAGE_BYTES;
    if (pages >= maximumPages) return;
    try {
      memory.grow(maximumPages - pages);
      return;
    } catch {
      // Grown by a thread meanwhile (read the size again), or not available.
    }
  }
}
