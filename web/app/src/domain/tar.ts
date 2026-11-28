// A tar archive of files, as one Blob made of the files' own Blobs (line WW, WW4b2: "Download all" of
// a makeblastdb run, the form in which NCBI hands out its BLAST databases): no byte of a file is read or
// copied into the page; only the 512-byte headers are made here. POSIX ustar headers; a path longer
// than ustar holds has a pax header (`path`) before it; a file of 8 GiB or more has its size in
// base-256 (GNU tar, read by every tar of today).

const BLOCK = 512;
const encoder = new TextEncoder();

export interface TarEntry {
  /** The path in the archive (`db.nin`, `sub/dir/db.nin`). */
  readonly name: string;
  readonly file: Blob;
}

/** The archive of `entries`, each file modified at `mtimeMs` (the end of the run that made them). */
export function tarArchive(entries: readonly TarEntry[], mtimeMs: number): Blob {
  const mtime = Math.max(0, Math.floor(mtimeMs / 1000));
  const parts: BlobPart[] = [];
  for (const entry of entries) {
    const path = encoder.encode(entry.name);
    const split = ustarPath(path);
    if (split === undefined) {
      const record = paxRecord('path', path);
      parts.push(header({ name: encoder.encode('PaxHeader'), prefix: new Uint8Array(0), size: record.length, mtime, type: 'x' }));
      parts.push(record, padding(record.length));
      parts.push(header({ name: path.subarray(Math.max(0, path.length - 100)), prefix: new Uint8Array(0), size: entry.file.size, mtime, type: '0' }));
    } else {
      parts.push(header({ ...split, size: entry.file.size, mtime, type: '0' }));
    }
    parts.push(entry.file, padding(entry.file.size));
  }
  parts.push(new Uint8Array(BLOCK * 2));
  return new Blob(parts, { type: 'application/x-tar' });
}

/** The path as ustar's `name` (100 bytes) and `prefix` (155 bytes, before a `/`), or undefined when it does not fit. */
function ustarPath(path: Uint8Array): { readonly name: Uint8Array; readonly prefix: Uint8Array } | undefined {
  if (path.length <= 100) return { name: path, prefix: new Uint8Array(0) };
  for (let i = path.length - 1; i > 0; i--) {
    if (path[i] !== 0x2f) continue;
    if (path.length - i - 1 > 100) return undefined;
    if (i <= 155) return { name: path.subarray(i + 1), prefix: path.subarray(0, i) };
  }
  return undefined;
}

/** A pax extended header record `<length> <key>=<value>\n`, whose length counts its own digits. */
function paxRecord(key: string, value: Uint8Array): Uint8Array<ArrayBuffer> {
  const body = value.length + key.length + 3;
  let length = body + String(body).length;
  if (String(length).length !== String(body).length) length = body + String(length).length;
  const bytes = new Uint8Array(length);
  bytes.set(encoder.encode(`${length} ${key}=`));
  bytes.set(value, length - value.length - 1);
  bytes[length - 1] = 0x0a;
  return bytes;
}

function padding(size: number): Uint8Array<ArrayBuffer> {
  return new Uint8Array((BLOCK - (size % BLOCK)) % BLOCK);
}

function header(fields: {
  readonly name: Uint8Array;
  readonly prefix: Uint8Array;
  readonly size: number;
  readonly mtime: number;
  readonly type: '0' | 'x';
}): Uint8Array<ArrayBuffer> {
  const block = new Uint8Array(BLOCK);
  block.set(fields.name.subarray(0, 100), 0);
  octal(block, 100, 8, 0o644);
  octal(block, 108, 8, 0);
  octal(block, 116, 8, 0);
  size(block, 124, fields.size);
  octal(block, 136, 12, fields.mtime);
  block.fill(0x20, 148, 156);
  block[156] = fields.type.charCodeAt(0);
  block.set(encoder.encode('ustar\u000000'), 257);
  block.set(fields.prefix.subarray(0, 155), 345);
  let sum = 0;
  for (const byte of block) sum += byte;
  octal(block, 148, 7, sum);
  return block;
}

/** Writes `value` as zero-padded octal digits ending with NUL in the `width` bytes at `at`. */
function octal(block: Uint8Array, at: number, width: number, value: number): void {
  block.set(encoder.encode(`${value.toString(8).padStart(width - 1, '0')}\u0000`), at);
}

/** The size field: octal below 8 GiB, else base-256 (its first byte 0x80, the size big-endian in the rest). */
function size(block: Uint8Array, at: number, value: number): void {
  if (value < 8 * 1024 ** 3) {
    octal(block, at, 12, value);
    return;
  }
  block[at] = 0x80;
  let rest = value;
  for (let i = at + 11; i > at; i--) {
    block[i] = rest % 256;
    rest = Math.floor(rest / 256);
  }
}
