// Whether an output holds data, for a failed run's kept output (line WQA's sweep, D15 round 2, 2026-11-21): the bytes of
// an output are data, those of an archive (outfmt 13 and 14: a ZIP of stored entries, docs/web/abi_v2.md §5) are the
// bytes of its entries, not the records that the archive has even when its files are empty.

/** Reads the bytes [start, end) of the output. */
type Read = (start: number, end: number) => Promise<Uint8Array>;

const LOCAL_FILE = 0x04034b50;
const END_OF_DIRECTORY = 0x06054b50;
const DIRECTORY_ENTRY = 0x02014b50;
const END_RECORD = 22;

/**
 * Whether an output of `length` bytes holds data: a text or binary output with bytes does; a ZIP archive does where one
 * of its entries has bytes. An archive that cannot be read as a classic ZIP (ZIP64, a comment, damaged records) counts
 * as holding data: it is kept, as an output that the engine wrote.
 */
export async function hasEntryBytes(read: Read, length: number): Promise<boolean> {
  if (length === 0) return false;
  if (length < END_RECORD) return true;
  const view = (bytes: Uint8Array) => new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const first = view(await read(0, 4));
  if (first.getUint32(0, true) !== LOCAL_FILE && first.getUint32(0, true) !== END_OF_DIRECTORY) return true;
  const end = view(await read(length - END_RECORD, length));
  if (end.getUint32(0, true) !== END_OF_DIRECTORY) return true;
  const entries = end.getUint16(10, true);
  const size = end.getUint32(12, true);
  const offset = end.getUint32(16, true);
  if (entries === 0xffff || size === 0xffffffff || offset === 0xffffffff || offset + size + END_RECORD !== length) return true;
  if (entries === 0) return false;
  const directory = view(await read(offset, offset + size));
  let at = 0;
  for (let i = 0; i < entries; i++) {
    if (at + 46 > size || directory.getUint32(at, true) !== DIRECTORY_ENTRY) return true;
    if (directory.getUint32(at + 24, true) > 0) return true;
    at += 46 + directory.getUint16(at + 28, true) + directory.getUint16(at + 30, true) + directory.getUint16(at + 32, true);
  }
  return false;
}
