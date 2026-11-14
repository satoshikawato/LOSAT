// The files of NCBI's application start and end that are not inputs or outputs of a search (line RG
// item E1w-3b; docs/web/abi_v2.md §4, §5): the registry file of `-conffile`, which a host gives the
// engine (`file_begin` kind 2), and the files that a run makes besides its outputs (stream 8: the file
// of `-logfile NAME`), which the engine sends after the run as a ZIP archive of stored entries.

/**
 * The registry file of `-conffile`: the file, under the name that the argv gives (the engine's
 * *start*). A Blob crosses to a worker without a copy of its bytes; the engine reads it in chunks.
 */
export interface RegistryFile {
  readonly name: string;
  readonly file: Blob;
}

/** A file that a run made besides its outputs (stream 8): the file of `-logfile NAME`, named NAME. */
export interface RunFile {
  readonly name: string;
  /** The bytes as the engine wrote them. */
  readonly bytes: Uint8Array;
}

const LOCAL_HEADER = 0x04034b50;
const CENTRAL_HEADER = 0x02014b50;
const END_OF_CENTRAL_DIRECTORY = 0x06054b50;
/** General purpose flags that a stored entry of the engine's archive may have: bit 11, a UTF-8 name. */
const UTF8_FLAG = 0x0800;

let crcTable: Uint32Array | undefined;

/** The CRC-32 of ZIP (ISO-HDLC: the reflected polynomial 0xEDB88320). */
function crc32(bytes: Uint8Array): number {
  if (crcTable === undefined) {
    crcTable = new Uint32Array(256);
    for (let n = 0; n < 256; n++) {
      let c = n;
      for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
      crcTable[n] = c >>> 0;
    }
  }
  let crc = ~0;
  for (const byte of bytes) crc = crcTable[(crc ^ byte) & 0xff]! ^ (crc >>> 8);
  return ~crc >>> 0;
}

/**
 * The files of the engine's archive of a run's files (stream 8), in its order. The archive is the
 * adapter's (web/adapter/src/zip.rs): stored entries (method 0), each with its CRC-32 and sizes in
 * its local header (no data descriptor) and a UTF-8 name, then a central directory that lists the
 * same entries, then the end record. Anything else, and an entry whose bytes do not have its CRC-32,
 * is an error of the engine's archive.
 */
export function readStoredZip(archive: Uint8Array): RunFile[] {
  const view = new DataView(archive.buffer, archive.byteOffset, archive.byteLength);
  const fail = (reason: string): never => {
    throw new Error(`the engine's archive of the run's files (ABI v2 stream 8) cannot be read: ${reason}`);
  };
  const u16 = (at: number) => (at + 2 <= archive.length ? view.getUint16(at, true) : fail(`it ends at byte ${archive.length}`));
  const u32 = (at: number) => (at + 4 <= archive.length ? view.getUint32(at, true) : fail(`it ends at byte ${archive.length}`));
  const decoder = new TextDecoder('utf-8', { fatal: true });
  const files: RunFile[] = [];
  const entries: Array<{ readonly offset: number; readonly crc: number; readonly size: number; readonly name: Uint8Array }> = [];
  let at = 0;
  while (u32(at) === LOCAL_HEADER) {
    const n = entries.length + 1;
    const flags = u16(at + 6);
    if ((flags & ~UTF8_FLAG) !== 0) fail(`entry ${n} has the flags 0x${flags.toString(16)} (a data descriptor, or encryption)`);
    const method = u16(at + 8);
    if (method !== 0) fail(`entry ${n} is not stored (method ${method})`);
    const crc = u32(at + 14);
    const size = u32(at + 18);
    if (u32(at + 22) !== size) fail(`entry ${n} has two sizes`);
    const nameLength = u16(at + 26);
    const data = at + 30 + nameLength + u16(at + 28);
    if (data + size > archive.length) fail(`entry ${n} ends past the archive`);
    const name = archive.subarray(at + 30, at + 30 + nameLength);
    const bytes = archive.slice(data, data + size);
    if (crc32(bytes) !== crc) fail(`the bytes of entry ${n} do not have its CRC-32`);
    let text: string;
    try {
      text = decoder.decode(name);
    } catch {
      return fail(`the name of entry ${n} is not UTF-8`);
    }
    entries.push({ offset: at, crc, size, name });
    files.push({ name: text, bytes });
    at = data + size;
  }
  // The central directory lists the same entries, and the end record ends the archive.
  const start = at;
  for (const [index, entry] of entries.entries()) {
    const n = index + 1;
    if (u32(at) !== CENTRAL_HEADER) fail(`the central directory has no header for entry ${n}`);
    const nameLength = u16(at + 28);
    const name = archive.subarray(at + 46, at + 46 + nameLength);
    const same =
      u16(at + 10) === 0 &&
      u32(at + 16) === entry.crc &&
      u32(at + 20) === entry.size &&
      u32(at + 24) === entry.size &&
      u32(at + 42) === entry.offset &&
      name.length === entry.name.length &&
      name.every((byte, i) => byte === entry.name[i]);
    if (!same) fail(`the central directory differs from entry ${n}`);
    at += 46 + nameLength + u16(at + 30) + u16(at + 32);
  }
  if (u32(at) !== END_OF_CENTRAL_DIRECTORY) fail(`no end record after ${entries.length} entries`);
  const count = u16(at + 10);
  if (u16(at + 8) !== count || count !== entries.length || u32(at + 12) !== at - start || u32(at + 16) !== start) {
    fail('the end record does not describe the central directory');
  }
  if (at + 22 + u16(at + 20) !== archive.length) fail('bytes follow the end record');
  return files;
}
