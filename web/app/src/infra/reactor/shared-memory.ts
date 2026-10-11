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
