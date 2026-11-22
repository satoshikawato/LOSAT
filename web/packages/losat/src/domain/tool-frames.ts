// The frames of stream 12 of a tool run (docs/web/abi_v2.md §5, version 13; line NP WW3a, D-WW3-6, D-WW3-9, 2026-11-22):
// the files that the run writes, each as it is written. Every `emit` on the stream is one frame: the file's
// number in the run (0, 1, ... in the order the run opens its files) and the frame's kind, each an unsigned
// 32-bit little-endian integer, then a payload. A file's frames come as its beginning, its bytes (each with the
// offset at which they are written), its truncations and its end; the frames of two open files may alternate.
// A run that opens one file twice (`-out X -logfile X`) writes ONE file with POSIX semantics: each handle has its own
// position, so bytes may come at an offset before the file's end, and the run empties the file when it opens it again
// for writing (a truncate to 0).

/** Kind 0: the beginning of a file; the payload is its name as the argv writes it, UTF-8. The file begins empty. */
export const FRAME_BEGIN = 0;
/** Kind 1: bytes of a file; the payload is the offset (u64 LE) at which they are written, then the bytes (at most 1 MiB less 16 bytes). */
export const FRAME_BYTES = 1;
/** Kind 2: the end of a file, which comes when the run ends, after it closed the file's last handle; no payload. */
export const FRAME_END = 2;
/** Kind 3: the file's length becomes the payload's u64 LE (cut, or zero bytes added when it is longer than the file). */
export const FRAME_TRUNCATE = 3;

/** Bytes of a frame before its payload: the file's number and the kind. */
export const FRAME_HEADER_BYTES = 8;

/** Bytes of a u64 in a payload. */
const U64_BYTES = 8;

/**
 * One frame of stream 12, read (`parseToolFrame`). `offset` and `length` are numbers: a value above
 * `Number.MAX_SAFE_INTEGER` is an error.
 */
export type ToolFrame =
  | { readonly number: number; readonly kind: 'begin'; readonly name: string }
  | { readonly number: number; readonly kind: 'bytes'; readonly offset: number; readonly bytes: Uint8Array }
  | { readonly number: number; readonly kind: 'truncate'; readonly length: number }
  | { readonly number: number; readonly kind: 'end' };

/**
 * Reads the frame that an `emit` on stream 12 sent. `bytes` of a frame of bytes is a view of `emitted`, which is
 * valid for the call only: `ownedFrame` copies it. A frame that breaks the ABI is an error: fewer than 8 bytes, a kind
 * other than 0 to 3, a name that is not UTF-8, an end with a payload, a frame of bytes with fewer than 8 bytes of
 * payload (its offset), a truncate whose payload is not 8 bytes, an offset or a length above `Number.MAX_SAFE_INTEGER`.
 */
export function parseToolFrame(emitted: Uint8Array): ToolFrame {
  if (emitted.length < FRAME_HEADER_BYTES) throw new Error(`a frame of stream 12 has ${emitted.length} bytes, fewer than its 8 of header`);
  const view = new DataView(emitted.buffer, emitted.byteOffset, emitted.byteLength);
  const number = view.getUint32(0, true);
  const kind = view.getUint32(4, true);
  const payload = emitted.subarray(FRAME_HEADER_BYTES);
  if (kind === FRAME_BEGIN) {
    try {
      return { number, kind: 'begin', name: new TextDecoder('utf-8', { fatal: true }).decode(payload) };
    } catch {
      throw new Error(`the name of file ${number} of a tool run is not UTF-8`);
    }
  }
  if (kind === FRAME_BYTES) {
    if (payload.length < U64_BYTES) throw new Error(`the bytes of file ${number} of a tool run have ${payload.length} bytes of payload, fewer than their 8 of offset`);
    return { number, kind: 'bytes', offset: safeInteger(view, FRAME_HEADER_BYTES, `the offset of the bytes of file ${number}`), bytes: payload.subarray(U64_BYTES) };
  }
  if (kind === FRAME_END) {
    if (payload.length !== 0) throw new Error(`the end of file ${number} of a tool run has ${payload.length} bytes of payload`);
    return { number, kind: 'end' };
  }
  if (kind === FRAME_TRUNCATE) {
    if (payload.length !== U64_BYTES) throw new Error(`the truncate of file ${number} of a tool run has ${payload.length} bytes of payload, not 8`);
    return { number, kind: 'truncate', length: safeInteger(view, FRAME_HEADER_BYTES, `the length of the truncate of file ${number}`) };
  }
  throw new Error(`a frame of stream 12 has the kind ${kind}, not 0, 1, 2 or 3`);
}

/** The u64 LE at `at`, as a number; an error for one that a number cannot hold exactly. */
function safeInteger(view: DataView, at: number, what: string): number {
  const value = view.getBigUint64(at, true);
  if (value > BigInt(Number.MAX_SAFE_INTEGER)) throw new Error(`${what} of a tool run is ${value}, more than ${Number.MAX_SAFE_INTEGER}`);
  return Number(value);
}

/** The frame with its bytes copied, so that it outlives the `emit` that sent it. */
export function ownedFrame(frame: ToolFrame): ToolFrame {
  return frame.kind === 'bytes' ? { ...frame, bytes: frame.bytes.slice() } : frame;
}
