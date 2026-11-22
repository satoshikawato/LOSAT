// The frames of stream 12 of a tool run (docs/web/abi_v2.md §5, version 13; line NP WW3a, D-WW3-6, 2026-11-22):
// the files that the run writes, each as it is written. Every `emit` on the stream is one frame: the file's
// number in the run (0, 1, ... in the order the run opens its files) and the frame's kind, each an unsigned
// 32-bit little-endian integer, then a payload. A file's frames come as its beginning, its bytes, its end; the
// frames of two open files may alternate, and two openings of one name are two files.

/** Kind 0: the beginning of a file; the payload is its name as the argv writes it, UTF-8. */
export const FRAME_BEGIN = 0;
/** Kind 1: bytes of a file (at most 1 MiB less 8 bytes in a frame). */
export const FRAME_BYTES = 1;
/** Kind 2: the end of a file, which comes when the run closes it; no payload. */
export const FRAME_END = 2;

/** Bytes of a frame before its payload: the file's number and the kind. */
export const FRAME_HEADER_BYTES = 8;

/** One frame of stream 12, read (`parseToolFrame`). */
export type ToolFrame =
  | { readonly number: number; readonly kind: 'begin'; readonly name: string }
  | { readonly number: number; readonly kind: 'bytes'; readonly bytes: Uint8Array }
  | { readonly number: number; readonly kind: 'end' };

/**
 * Reads the frame that an `emit` on stream 12 sent. `bytes` of a frame of bytes is a view of `emitted`, which is
 * valid for the call only: `ownedFrame` copies it. A frame that breaks the ABI (fewer than 8 bytes, a kind other
 * than 0, 1 and 2, a name that is not UTF-8, an end with a payload) is an error.
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
  if (kind === FRAME_BYTES) return { number, kind: 'bytes', bytes: payload };
  if (kind === FRAME_END) {
    if (payload.length !== 0) throw new Error(`the end of file ${number} of a tool run has ${payload.length} bytes of payload`);
    return { number, kind: 'end' };
  }
  throw new Error(`a frame of stream 12 has the kind ${kind}, not 0, 1 or 2`);
}

/** The frame with its bytes copied, so that it outlives the `emit` that sent it. */
export function ownedFrame(frame: ToolFrame): ToolFrame {
  return frame.kind === 'bytes' ? { ...frame, bytes: frame.bytes.slice() } : frame;
}
