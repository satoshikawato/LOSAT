// The run's working directory of the browser host (src/infra/reactor/run-directory.ts) through its WASI functions, as
// the engine calls them: the preopen `.`, opening a file by its path, reading it in place, a descriptor that one
// thread opened read on another, and a search's directory that writes nothing (line WQA D12a, 2026-11-20); a tool's
// writable directory, the overlay of the run's files on the user's (line NP WW4b1, D-WW4b-1, 2026-11-27): each call
// that makeblastdb's `std::fs` makes, answered as Linux answers it, and the files that the run leaves.
import { describe, expect, it } from 'vitest';
import { CHUNK_BYTES, FILE_FD_BASE, pathNames, PREOPEN_FD, RUN_DIRECTORY, RunDirectory, TABLE_BYTES, type BlobReader } from '../src/infra/reactor/run-directory.ts';

const ERRNO = {
  SUCCESS: 0,
  BADF: 8,
  BUSY: 10,
  EXIST: 20,
  INVAL: 28,
  IO: 29,
  ISDIR: 31,
  NAMETOOLONG: 37,
  NOENT: 44,
  NOTDIR: 54,
  NOTEMPTY: 55,
  ROFS: 69,
  NOTCAPABLE: 76,
};
const RIGHT_FD_READ = 1n << 1n;
const RIGHT_FD_WRITE = 1n << 6n;
const OFLAGS_CREAT = 1;
const OFLAGS_DIRECTORY = 2;
const OFLAGS_TRUNC = 8;
const FILETYPE_DIRECTORY = 3;
const FILETYPE_REGULAR_FILE = 4;

type Call = (...args: Array<number | bigint>) => number;

/** The bytes of a Blob, as a dedicated worker's FileReaderSync gives them; the reads that it was asked for. */
function syncReader(contents: Map<Blob, Uint8Array>): { read: BlobReader; reads: Array<[number, number]> } {
  const reads: Array<[number, number]> = [];
  return {
    reads,
    read: (blob, start, end) => {
      reads.push([start, end]);
      return contents.get(blob)!.slice(start, end);
    },
  };
}

/** A guest of 64 KiB, the directory's WASI functions in front of a shim that answers EBADF, and helpers to call them. */
function guest(directory: RunDirectory) {
  const memory = new WebAssembly.Memory({ initial: 1 });
  const base: Record<string, Call> = {};
  for (const name of [
    'fd_prestat_get',
    'fd_prestat_dir_name',
    'path_open',
    'path_filestat_get',
    'path_readlink',
    'path_create_directory',
    'path_unlink_file',
    'path_remove_directory',
    'fd_filestat_get',
    'fd_fdstat_get',
    'fd_read',
    'fd_pread',
    'fd_seek',
    'fd_tell',
    'fd_write',
    'fd_close',
    'fd_advise',
    'fd_readdir',
    'fd_fdstat_set_flags',
    'fd_sync',
    'fd_pwrite',
    'fd_filestat_set_size',
    'fd_filestat_set_times',
    'fd_allocate',
    'path_rename',
    'path_filestat_set_times',
  ]) {
    base[name] = () => ERRNO.BADF;
  }
  const wasi = directory.imports(base as unknown as WebAssembly.ModuleImports, () => memory) as unknown as Record<string, Call>;
  const view = () => new DataView(memory.buffer);
  const PATH = 1024;
  const OUT = 2048;
  const IOVS = 3072;
  const DATA = 4096;
  const encoder = new TextEncoder();
  const path = (text: string) => {
    const bytes = encoder.encode(text);
    new Uint8Array(memory.buffer).set(bytes, PATH);
    return [PATH, bytes.length] as const;
  };
  return {
    wasi,
    view,
    path,
    open(name: string, oflags = 0, rights = RIGHT_FD_READ): { errno: number; fd: number } {
      const errno = wasi.path_open!(PREOPEN_FD, 0, ...path(name), oflags, rights, 0n, 0, OUT);
      return { errno, fd: view().getUint32(OUT, true) };
    },
    /** `path_open` of `name` relative to the open descriptor `dirfd`. */
    openFrom(dirfd: number, name: string): { errno: number; fd: number } {
      const errno = wasi.path_open!(dirfd, 0, ...path(name), 0, RIGHT_FD_READ, 0n, 0, OUT);
      return { errno, fd: view().getUint32(OUT, true) };
    },
    stat(name: string): { errno: number; type: number; size: bigint } {
      const errno = wasi.path_filestat_get!(PREOPEN_FD, 0, ...path(name), OUT);
      return { errno, type: view().getUint8(OUT + 16), size: view().getBigUint64(OUT + 32, true) };
    },
    /** fd_read into two vectors of `first` and `second` bytes; the bytes read. */
    read(fd: number, first: number, second = 0): { errno: number; bytes: string } {
      view().setUint32(IOVS, DATA, true);
      view().setUint32(IOVS + 4, first, true);
      view().setUint32(IOVS + 8, DATA + first, true);
      view().setUint32(IOVS + 12, second, true);
      const errno = wasi.fd_read!(fd, IOVS, 2, OUT);
      const count = view().getUint32(OUT, true);
      return { errno, bytes: new TextDecoder().decode(new Uint8Array(memory.buffer, DATA, errno === 0 ? count : 0)) };
    },
    pread(fd: number, length: number, offset: bigint): { errno: number; bytes: string } {
      view().setUint32(IOVS, DATA, true);
      view().setUint32(IOVS + 4, length, true);
      const errno = wasi.fd_pread!(fd, IOVS, 1, offset, OUT);
      const count = view().getUint32(OUT, true);
      return { errno, bytes: new TextDecoder().decode(new Uint8Array(memory.buffer, DATA, errno === 0 ? count : 0)) };
    },
    seek(fd: number, offset: bigint, whence: number): { errno: number; at: bigint } {
      const errno = wasi.fd_seek!(fd, offset, whence, OUT);
      return { errno, at: view().getBigUint64(OUT, true) };
    },
  };
}

function directoryOf(files: Record<string, string>, table?: SharedArrayBuffer) {
  const contents = new Map<Blob, Uint8Array>();
  const list = Object.entries(files).map(([name, text]) => {
    const file = new Blob([text]);
    contents.set(file, new TextEncoder().encode(text));
    return { name, file };
  });
  const reader = syncReader(contents);
  const directory = new RunDirectory(table, reader.read);
  directory.open(list);
  return { directory, list, reader, contents };
}

describe("the run's working directory", () => {
  it('is the preopen `/work`, the first directory that WASI asks about', () => {
    const { directory } = directoryOf({});
    const { wasi, view } = guest(directory);
    expect(wasi.fd_prestat_get!(PREOPEN_FD, 0)).toBe(ERRNO.SUCCESS);
    expect(view().getUint8(0)).toBe(0);
    expect(view().getUint32(4, true)).toBe(RUN_DIRECTORY.length);
    expect(wasi.fd_prestat_dir_name!(PREOPEN_FD, 16, RUN_DIRECTORY.length - 1)).toBe(ERRNO.NAMETOOLONG);
    expect(wasi.fd_prestat_dir_name!(PREOPEN_FD, 16, RUN_DIRECTORY.length)).toBe(ERRNO.SUCCESS);
    expect(new TextDecoder().decode(new Uint8Array(view().buffer, 16, RUN_DIRECTORY.length))).toBe('/work');
    // Descriptor 4 is not a preopen: the shim's answer.
    expect(wasi.fd_prestat_get!(4, 0)).toBe(ERRNO.BADF);
  });

  it('opens a file by its path, and reads it in place', () => {
    const { directory, reader } = directoryOf({ 'db/nt.nin': 'header-bytes', 'ids.txt': 'P1\nP2\n' });
    const g = guest(directory);
    for (const name of ['db/nt.nin', './db/nt.nin', 'db//nt.nin', 'db/../db/nt.nin', 'db/./nt.nin']) {
      const opened = g.open(name);
      expect(opened.errno, name).toBe(ERRNO.SUCCESS);
      expect(opened.fd).toBeGreaterThanOrEqual(FILE_FD_BASE);
      expect(g.wasi.fd_close!(opened.fd)).toBe(ERRNO.SUCCESS);
    }
    const { fd } = g.open('db/nt.nin');
    expect(g.read(fd, 6, 3)).toEqual({ errno: 0, bytes: 'header-by' });
    expect(g.read(fd, 100)).toEqual({ errno: 0, bytes: 'tes' });
    // At the end: 0 bytes, as a file system's read.
    expect(g.read(fd, 100)).toEqual({ errno: 0, bytes: '' });
    // Only the bytes that each read asked for were taken from the Blob.
    expect(reader.reads).toEqual([
      [0, 9],
      [9, 12],
    ]);
    expect(g.pread(fd, 5, 7n)).toEqual({ errno: 0, bytes: 'bytes' });
    // fd_pread does not move the position.
    expect(g.seek(fd, 0n, 1)).toEqual({ errno: 0, at: 12n });
    expect(g.seek(fd, -5n, 2)).toEqual({ errno: 0, at: 7n });
    expect(g.read(fd, 3)).toEqual({ errno: 0, bytes: 'byt' });
    expect(g.seek(fd, -1n, 0).errno).toBe(ERRNO.INVAL);
    expect(g.wasi.fd_tell!(fd, 2048)).toBe(ERRNO.SUCCESS);
    expect(g.view().getBigUint64(2048, true)).toBe(10n);
    expect(g.stat('db/nt.nin')).toEqual({ errno: 0, type: 4, size: 12n });
    expect(g.wasi.fd_filestat_get!(fd, 2048)).toBe(ERRNO.SUCCESS);
    expect(g.view().getBigUint64(2048 + 32, true)).toBe(12n);
    expect(g.wasi.fd_close!(fd)).toBe(ERRNO.SUCCESS);
    expect(g.wasi.fd_close!(fd)).toBe(ERRNO.BADF);
    expect(g.read(fd, 1).errno).toBe(ERRNO.BADF);
  });

  it("answers for its folders and for what it does not hold as a file system does", () => {
    const { directory } = directoryOf({ 'db/nt.nin': 'x' });
    const g = guest(directory);
    expect(g.stat('db')).toMatchObject({ errno: 0, type: 3 });
    expect(g.stat('.')).toMatchObject({ errno: 0, type: 3 });
    expect(g.stat('db/nt.pin').errno).toBe(ERRNO.NOENT);
    expect(g.open('nt.nal').errno).toBe(ERRNO.NOENT);
    // Out of the directory: WASI's capability, as the command line's sandbox would refuse.
    expect(g.open('../etc/passwd').errno).toBe(ERRNO.NOTCAPABLE);
    // A path is looked up one name after another, as the OS looks it up: `..` after a name that is not there, or
    // any name after a file, fails (review of D12a, low 4).
    expect(g.open('x/../db/nt.nin').errno).toBe(ERRNO.NOENT);
    expect(g.stat('x/..').errno).toBe(ERRNO.NOENT);
    expect(g.open('db/nt.nin/../nt.nin').errno).toBe(ERRNO.NOTDIR);
    expect(g.stat('db/nt.nin/x').errno).toBe(ERRNO.NOTDIR);
    expect(g.stat('db/nt.nin/').errno).toBe(ERRNO.NOTDIR);
    const { fd } = g.open('db/nt.nin');
    expect(g.wasi.path_open!(fd, 0, 1024, 1, 0, RIGHT_FD_READ, 0n, 0, 2048)).toBe(ERRNO.NOTDIR);
    // No file of the directory is a link.
    g.view().setUint16(1024, 0x6264, true); // "db"
    expect(g.wasi.path_readlink!(PREOPEN_FD, 1024, 2, 0, 0, 2048)).toBe(ERRNO.INVAL);
  });

  // Measured 2026-11-23 with the command line: `windowmasker -ustat counts/u4.ascii -in fa` with `fa` a folder prints "basic_filebuf::xsgetn
  // error reading the file: Is a directory" (exit 3): on Linux the open of a folder for reading succeeds and its read fails with EISDIR.
  it('opens a folder for reading as Linux does, and fails its reads with EISDIR', () => {
    const { directory, reader } = directoryOf({ 'fa/dm2L_3k.fa': 'ACGT', 'fa/sub/x': 'x', 'ids.txt': 'P1\n' });
    const g = guest(directory);
    for (const name of ['fa', 'fa/', './fa', 'fa/sub/..', '.']) {
      const opened = g.open(name);
      expect(opened.errno, name).toBe(ERRNO.SUCCESS);
      expect(opened.fd).toBeGreaterThanOrEqual(FILE_FD_BASE);
      expect(g.wasi.fd_close!(opened.fd)).toBe(ERRNO.SUCCESS);
    }
    const { fd } = g.open('fa');
    expect(g.read(fd, 10).errno).toBe(ERRNO.ISDIR);
    expect(g.pread(fd, 10, 0n).errno).toBe(ERRNO.ISDIR);
    // No byte of any file was read from a Blob.
    expect(reader.reads).toEqual([]);
    // The stat calls say directory, and the descriptor is a directory's.
    expect(g.wasi.fd_filestat_get!(fd, 2048)).toBe(ERRNO.SUCCESS);
    expect(g.view().getUint8(2048 + 16)).toBe(FILETYPE_DIRECTORY);
    expect(g.wasi.fd_fdstat_get!(fd, 2048)).toBe(ERRNO.SUCCESS);
    expect(g.view().getUint8(2048)).toBe(FILETYPE_DIRECTORY);
    // Writes are not allowed on a descriptor that is open for reading; seeks and syncs are.
    expect(g.wasi.fd_write!(fd, 0, 0, 0)).toBe(ERRNO.BADF);
    expect(g.seek(fd, 0n, 0)).toEqual({ errno: 0, at: 0n });
    expect(g.wasi.fd_sync!(fd)).toBe(ERRNO.SUCCESS);
    // Paths are looked up from it, as the OS looks up a path from an open folder: and a file of it is opened.
    const inside = g.openFrom(fd, 'dm2L_3k.fa');
    expect(inside.errno).toBe(ERRNO.SUCCESS);
    expect(g.read(inside.fd, 4)).toEqual({ errno: 0, bytes: 'ACGT' });
    expect(g.openFrom(fd, 'sub/x').errno).toBe(ERRNO.SUCCESS);
    expect(g.openFrom(fd, '../ids.txt').errno).toBe(ERRNO.SUCCESS);
    expect(g.openFrom(fd, '../../x').errno).toBe(ERRNO.NOTCAPABLE);
    expect(g.openFrom(fd, 'nope').errno).toBe(ERRNO.NOENT);
    // The descriptor ends with the close, and with the run.
    expect(g.wasi.fd_close!(fd)).toBe(ERRNO.SUCCESS);
    expect(g.read(fd, 1).errno).toBe(ERRNO.BADF);
    expect(g.wasi.fd_close!(fd)).toBe(ERRNO.BADF);
    const again = g.open('fa').fd;
    directory.open([]);
    expect(g.read(again, 1).errno).toBe(ERRNO.BADF);
  });

  it('still refuses to open a folder for writing, with O_CREAT or O_TRUNC, or as a file', () => {
    const { directory } = directoryOf({ 'db/nt.nin': 'x' });
    const g = guest(directory);
    expect(g.open('db', 0, RIGHT_FD_READ | RIGHT_FD_WRITE).errno).toBe(ERRNO.ISDIR);
    expect(g.open('db', OFLAGS_TRUNC).errno).toBe(ERRNO.ISDIR);
    expect(g.open('db', OFLAGS_CREAT, RIGHT_FD_WRITE).errno).toBe(ERRNO.ISDIR);
    expect(g.open('db', OFLAGS_CREAT | 4).errno).toBe(ERRNO.EXIST);
    // O_DIRECTORY on a folder opens it; on a file, ENOTDIR.
    expect(g.open('db', OFLAGS_DIRECTORY).errno).toBe(ERRNO.SUCCESS);
    expect(g.open('db/nt.nin', OFLAGS_DIRECTORY).errno).toBe(ERRNO.NOTDIR);
    // The stat of a folder by its path is unchanged: directory.
    expect(g.stat('db')).toMatchObject({ errno: 0, type: FILETYPE_DIRECTORY });
    expect(g.stat('db/nt.nin')).toMatchObject({ errno: 0, type: FILETYPE_REGULAR_FILE });
  });

  it('lists the entries of a folder descriptor: `.`, `..`, then its files and folders', () => {
    const { directory } = directoryOf({ 'fa/b.fa': 'x', 'fa/a.fa': 'x', 'fa/sub/x': 'x', 'top': 'x' });
    const g = guest(directory);
    const names = (fd: number, room: number, cookie = 0n) => {
      expect(g.wasi.fd_readdir!(fd, 8192, room, cookie, 2048)).toBe(ERRNO.SUCCESS);
      const used = g.view().getUint32(2048, true);
      const found: Array<{ name: string; type: number; next: bigint }> = [];
      for (let at = 8192; at + 24 <= 8192 + used; ) {
        const length = g.view().getUint32(at + 16, true);
        if (at + 24 + length > 8192 + used) break;
        found.push({
          name: new TextDecoder().decode(new Uint8Array(g.view().buffer, at + 24, length)),
          type: g.view().getUint8(at + 20),
          next: g.view().getBigUint64(at, true),
        });
        at += 24 + length;
      }
      return { used, found };
    };
    const folder = g.open('fa').fd;
    const all = names(folder, 4096);
    expect(all.found.map((entry) => [entry.name, entry.type])).toEqual([
      ['.', FILETYPE_DIRECTORY],
      ['..', FILETYPE_DIRECTORY],
      ['a.fa', FILETYPE_REGULAR_FILE],
      ['b.fa', FILETYPE_REGULAR_FILE],
      ['sub', FILETYPE_DIRECTORY],
    ]);
    // Less than the buffer is used: that is the end. A buffer that is too small is filled, and the cookie goes on.
    expect(all.used).toBeLessThan(4096);
    const small = names(folder, 60);
    expect(small.used).toBe(60);
    expect(small.found.map((entry) => entry.name)).toEqual(['.', '..']);
    expect(names(folder, 4096, small.found.at(-1)!.next).found.map((entry) => entry.name)).toEqual(['a.fa', 'b.fa', 'sub']);
    expect(names(g.open('.').fd, 4096).found.map((entry) => entry.name)).toEqual(['.', '..', 'fa', 'top']);
    // A file is no directory.
    expect(g.wasi.fd_readdir!(g.open('top').fd, 8192, 100, 0n, 2048)).toBe(ERRNO.NOTDIR);
  });

  it('is read-only', () => {
    const { directory } = directoryOf({ 'db/nt.nin': 'x' });
    const g = guest(directory);
    expect(g.open('db/nt.nin', 0, RIGHT_FD_READ | RIGHT_FD_WRITE).errno).toBe(ERRNO.ROFS);
    expect(g.open('db/nt.nin', OFLAGS_TRUNC).errno).toBe(ERRNO.ROFS);
    expect(g.open('out.txt', OFLAGS_CREAT, RIGHT_FD_WRITE).errno).toBe(ERRNO.ROFS);
    expect(g.wasi.path_create_directory!(PREOPEN_FD, 1024, 1)).toBe(ERRNO.ROFS);
    expect(g.wasi.path_unlink_file!(PREOPEN_FD, 1024, 1)).toBe(ERRNO.ROFS);
    const { fd } = g.open('db/nt.nin');
    expect(g.wasi.fd_write!(fd, 0, 0, 0)).toBe(ERRNO.BADF);
  });

  it('gives a failed read of a Blob as EIO', () => {
    const file = new Blob(['abc']);
    const directory = new RunDirectory(undefined, () => {
      throw new DOMException('the file changed', 'NotReadableError');
    });
    directory.open([{ name: 'a', file }]);
    const g = guest(directory);
    const { fd } = g.open('a');
    expect(g.read(fd, 3).errno).toBe(ERRNO.IO);
    // A thread that cannot read synchronously (no FileReaderSync) fails the same way.
    const without = new RunDirectory(undefined, undefined);
    without.open([{ name: 'a', file }]);
    const h = guest(without);
    expect(h.read(h.open('a').fd, 3).errno).toBe(ERRNO.IO);
  });

  it("shares its descriptors with the instance's threads: one thread opens, another reads", () => {
    const table = new SharedArrayBuffer(TABLE_BYTES);
    const { list, contents } = directoryOf({ 'nt.00.nsq': '0123456789' }, table);
    // The Engine worker's directory and a thread worker's, on the same table; the thread takes the run's files.
    const main = new RunDirectory(table, syncReader(contents).read);
    main.open(list);
    const thread = new RunDirectory(table, syncReader(contents).read);
    thread.take(list);
    const a = guest(main);
    const b = guest(thread);
    const { fd } = a.open('nt.00.nsq');
    expect(a.read(fd, 4).bytes).toBe('0123');
    // The position is the descriptor's, whichever thread reads.
    expect(b.read(fd, 3).bytes).toBe('456');
    expect(a.read(fd, 10).bytes).toBe('789');
    expect(b.wasi.fd_close!(fd)).toBe(ERRNO.SUCCESS);
    expect(a.read(fd, 1).errno).toBe(ERRNO.BADF);
  });

  it("ends an earlier run's descriptors when a run opens its files", () => {
    const { directory, list } = directoryOf({ a: 'one' });
    const g = guest(directory);
    const { fd } = g.open('a');
    directory.open(list);
    expect(g.read(fd, 1).errno).toBe(ERRNO.BADF);
    directory.open([]);
    expect(g.open('a').errno).toBe(ERRNO.NOENT);
  });

  it('refuses names that are not one path each', () => {
    const file = new Blob(['x']);
    const directory = new RunDirectory(undefined, () => new Uint8Array(0));
    expect(() => directory.open([{ name: '', file }])).toThrow(/not a path/);
    expect(() => directory.open([{ name: '../a', file }])).toThrow(/not a path/);
    expect(() => directory.open([{ name: 'a', file }, { name: './a', file }])).toThrow(/two files are named a/);
    expect(() => directory.open([{ name: 'db', file }, { name: 'db/a', file }])).toThrow(/both a file and a folder/);
    expect(pathNames('a/./b//c')).toEqual(['a', 'b', 'c']);
    expect(pathNames('a/../../b')).toBeUndefined();
  });
});

/** The files that a run left, by name, with their bytes (`takeWritten`, each piece at most a chunk). */
function writtenOf(directory: RunDirectory): Array<[string, Uint8Array]> {
  return directory.takeWritten().map((file) => {
    const pieces = [...file.pieces()];
    for (const piece of pieces) expect(piece.length).toBeLessThanOrEqual(CHUNK_BYTES);
    const bytes = new Uint8Array(pieces.reduce((sum, piece) => sum + piece.length, 0));
    let at = 0;
    for (const piece of pieces) {
      bytes.set(piece, at);
      at += piece.length;
    }
    expect(bytes.length).toBe(file.size);
    return [file.name, bytes];
  });
}

const textOf = (bytes: Uint8Array) => new TextDecoder().decode(bytes);

/** A tool run's directory over `files` (writable), with the calls that Rust's `std::fs` makes on WASI. */
function writableOf(files: Record<string, string | Uint8Array>) {
  const contents = new Map<Blob, Uint8Array>();
  const list = Object.entries(files).map(([name, value]) => {
    const data = typeof value === 'string' ? new TextEncoder().encode(value) : value;
    const file = new Blob([data as BlobPart]);
    contents.set(file, data);
    return { name, file };
  });
  const reader = syncReader(contents);
  const directory = new RunDirectory(undefined, reader.read);
  directory.open(list, { writable: true });
  const g = guest(directory);
  const OUT = 2048;
  const IOVS = 3072;
  const DATA = 4096;
  const SECOND_PATH = 1536;
  const encoder = new TextEncoder();
  const memory8 = () => new Uint8Array(g.view().buffer);
  /** Puts the parts into the guest's memory as vectors; returns their count. */
  const vectors = (parts: Array<string | Uint8Array>) => {
    let at = DATA;
    parts.forEach((part, i) => {
      const data = typeof part === 'string' ? encoder.encode(part) : part;
      memory8().set(data, at);
      g.view().setUint32(IOVS + 8 * i, at, true);
      g.view().setUint32(IOVS + 8 * i + 4, data.length, true);
      at += data.length;
    });
    return parts.length;
  };
  const twoPaths = (from: string, to: string) => {
    const [fromAt, fromLength] = g.path(from);
    const data = encoder.encode(to);
    memory8().set(data, SECOND_PATH);
    return [fromAt, fromLength, SECOND_PATH, data.length] as const;
  };
  return {
    ...g,
    directory,
    files: list,
    reader,
    /** `path_open` as Rust's `OpenOptions` asks: the rights of reading and writing, O_APPEND, O_CREAT, O_EXCL, O_TRUNC. */
    openWith(name: string, how: { read?: boolean; write?: boolean; append?: boolean; create?: boolean; truncate?: boolean; createNew?: boolean }) {
      const rights = (how.read ? RIGHT_FD_READ : 0n) | (how.write || how.append ? RIGHT_FD_WRITE : 0n);
      const oflags = (how.create || how.createNew ? OFLAGS_CREAT : 0) | (how.createNew ? 4 : 0) | (how.truncate ? OFLAGS_TRUNC : 0);
      const errno = g.wasi.path_open!(PREOPEN_FD, 0, ...g.path(name), oflags, rights, 0n, how.append ? 1 : 0, OUT);
      return { errno, fd: g.view().getUint32(OUT, true) };
    },
    write(fd: number, ...parts: Array<string | Uint8Array>) {
      const errno = g.wasi.fd_write!(fd, IOVS, vectors(parts), OUT);
      return { errno, count: errno === 0 ? g.view().getUint32(OUT, true) : 0 };
    },
    pwrite(fd: number, part: string, offset: bigint) {
      const errno = g.wasi.fd_pwrite!(fd, IOVS, vectors([part]), offset, OUT);
      return { errno, count: errno === 0 ? g.view().getUint32(OUT, true) : 0 };
    },
    tell(fd: number) {
      expect(g.wasi.fd_tell!(fd, OUT)).toBe(ERRNO.SUCCESS);
      return g.view().getBigUint64(OUT, true);
    },
    /** The whole file `name`, read through a descriptor of its own. */
    contents(name: string): string {
      const { errno, fd } = this.openWith(name, { read: true });
      expect(errno, name).toBe(ERRNO.SUCCESS);
      let text = '';
      for (;;) {
        const part = g.read(fd, 1000);
        expect(part.errno).toBe(ERRNO.SUCCESS);
        if (part.bytes === '') break;
        text += part.bytes;
      }
      expect(g.wasi.fd_close!(fd)).toBe(ERRNO.SUCCESS);
      return text;
    },
    fstat(fd: number) {
      expect(g.wasi.fd_filestat_get!(fd, OUT)).toBe(ERRNO.SUCCESS);
      return { type: g.view().getUint8(OUT + 16), nlink: g.view().getBigUint64(OUT + 24, true), size: g.view().getBigUint64(OUT + 32, true), mtime: g.view().getBigUint64(OUT + 48, true) };
    },
    fdstat(fd: number) {
      expect(g.wasi.fd_fdstat_get!(fd, OUT)).toBe(ERRNO.SUCCESS);
      return { flags: g.view().getUint16(OUT + 2, true), rights: g.view().getBigUint64(OUT + 8, true) };
    },
    mkdir: (name: string) => g.wasi.path_create_directory!(PREOPEN_FD, ...g.path(name)),
    unlink: (name: string) => g.wasi.path_unlink_file!(PREOPEN_FD, ...g.path(name)),
    rmdir: (name: string) => g.wasi.path_remove_directory!(PREOPEN_FD, ...g.path(name)),
    rename(from: string, to: string) {
      const [fromAt, fromLength, toAt, toLength] = twoPaths(from, to);
      return g.wasi.path_rename!(PREOPEN_FD, fromAt, fromLength, PREOPEN_FD, toAt, toLength);
    },
    /** The names of a folder's listing from `cookie` in a buffer of `room` bytes, and the cookie after the last. */
    list(fd: number, room: number, cookie = 0n) {
      const BUFFER = 8192;
      expect(g.wasi.fd_readdir!(fd, BUFFER, room, cookie, OUT)).toBe(ERRNO.SUCCESS);
      const used = g.view().getUint32(OUT, true);
      const names: string[] = [];
      let next = cookie;
      for (let at = BUFFER; at + 24 <= BUFFER + used; ) {
        const length = g.view().getUint32(at + 16, true);
        if (at + 24 + length > BUFFER + used) break;
        names.push(new TextDecoder().decode(new Uint8Array(g.view().buffer, at + 24, length)));
        next = g.view().getBigUint64(at, true);
        at += 24 + length;
      }
      return { names, next };
    },
  };
}

describe("a tool run's writable directory", () => {
  it('makes, writes and reads back a file, as `File::create`, `write_vectored` and `File::open` do', () => {
    const w = writableOf({});
    const made = w.openWith('db.pin', { write: true, create: true, truncate: true });
    expect(made.errno).toBe(ERRNO.SUCCESS);
    expect(w.write(made.fd, 'hello ')).toEqual({ errno: 0, count: 6 });
    expect(w.write(made.fd, 'wor', 'ld')).toEqual({ errno: 0, count: 5 });
    expect(w.stat('db.pin')).toEqual({ errno: 0, type: FILETYPE_REGULAR_FILE, size: 11n });
    expect(w.fstat(made.fd)).toMatchObject({ type: FILETYPE_REGULAR_FILE, nlink: 1n, size: 11n });
    expect(w.fstat(made.fd).mtime).toBeGreaterThan(0n);
    // Open for writing only: a read is EBADF, as Linux's read of an O_WRONLY descriptor; the rights say so.
    expect(w.read(made.fd, 4).errno).toBe(ERRNO.BADF);
    expect(w.fdstat(made.fd).rights & RIGHT_FD_WRITE).toBe(RIGHT_FD_WRITE);
    expect(w.fdstat(made.fd).rights & RIGHT_FD_READ).toBe(0n);
    expect(w.contents('db.pin')).toBe('hello world');
    const reading = w.openWith('db.pin', { read: true });
    expect(w.fdstat(reading.fd).rights & RIGHT_FD_WRITE).toBe(0n);
    // Open for reading: a write is EBADF and a cut EINVAL (Linux's ftruncate).
    expect(w.write(reading.fd, 'x').errno).toBe(ERRNO.BADF);
    expect(w.wasi.fd_filestat_set_size!(reading.fd, 0n)).toBe(ERRNO.INVAL);
    expect(w.wasi.fd_allocate!(reading.fd, 0n, 1n)).toBe(ERRNO.BADF);
    // `create_new` of a name that is there is EEXIST; a missing name ENOENT without O_CREAT and in a missing folder;
    // `x/` with O_CREAT EISDIR; a name after a file ENOTDIR.
    expect(w.openWith('db.pin', { write: true, createNew: true }).errno).toBe(ERRNO.EXIST);
    expect(w.openWith('none', { read: true }).errno).toBe(ERRNO.NOENT);
    expect(w.openWith('no/x', { write: true, create: true }).errno).toBe(ERRNO.NOENT);
    expect(w.openWith('new/', { write: true, create: true }).errno).toBe(ERRNO.ISDIR);
    expect(w.openWith('db.pin/x', { write: true, create: true }).errno).toBe(ERRNO.NOTDIR);
    expect(w.wasi.path_open!(PREOPEN_FD, 0, ...w.path('dir'), OFLAGS_CREAT | OFLAGS_DIRECTORY, RIGHT_FD_READ, 0n, 0, 2048)).toBe(ERRNO.INVAL);
    // O_TRUNC empties it.
    const again = w.openWith('db.pin', { write: true, truncate: true });
    expect(w.stat('db.pin').size).toBe(0n);
    expect(w.write(again.fd, 'new').count).toBe(3);
    // Not a link; syncs succeed.
    expect(w.wasi.path_readlink!(PREOPEN_FD, ...w.path('db.pin'), 0, 0, 2048)).toBe(ERRNO.INVAL);
    expect(w.wasi.fd_sync!(again.fd)).toBe(ERRNO.SUCCESS);
    expect(writtenOf(w.directory).map(([name, bytes]) => [name, textOf(bytes)])).toEqual([['db.pin', 'new']]);
  });

  it('writes at the position, past the end with a hole of zeros, and at the end for O_APPEND, pwrite too (Linux)', () => {
    const w = writableOf({});
    const { fd } = w.openWith('f', { read: true, write: true, create: true });
    w.write(fd, 'abc');
    expect(w.seek(fd, 5n, 0)).toEqual({ errno: 0, at: 5n });
    w.write(fd, 'Z');
    expect(w.pread(fd, 10, 0n).bytes).toBe('abc\0\0Z');
    w.seek(fd, 1n, 0);
    w.write(fd, 'B');
    // pwrite writes at its offset and leaves the position.
    expect(w.pwrite(fd, 'C', 2n)).toEqual({ errno: 0, count: 1 });
    expect(w.tell(fd)).toBe(2n);
    expect(w.seek(fd, 0n, 2)).toEqual({ errno: 0, at: 6n });
    // A descriptor that appends writes at the end wherever it is, also with pwrite, as Linux's O_APPEND does.
    const log = w.openWith('f', { append: true });
    expect(w.fdstat(log.fd).flags).toBe(1);
    w.write(log.fd, '!');
    w.seek(log.fd, 0n, 0);
    w.write(log.fd, '?');
    expect(w.pwrite(log.fd, '#', 0n).count).toBe(1);
    expect(w.contents('f')).toBe('aBC\0\0Z!?#');
    // fcntl(F_SETFL) sets or clears O_APPEND on a descriptor open for writing, and nothing else.
    expect(w.wasi.fd_fdstat_set_flags!(log.fd, 0)).toBe(ERRNO.SUCCESS);
    w.seek(log.fd, 0n, 0);
    w.write(log.fd, 'a');
    expect(w.contents('f')).toBe('aBC\0\0Z!?#');
    expect(w.wasi.fd_fdstat_set_flags!(log.fd, 4)).toBe(ERRNO.NOTCAPABLE);
    expect(w.wasi.fd_fdstat_set_flags!(w.openWith('f', { read: true }).fd, 1)).toBe(ERRNO.NOTCAPABLE);
    expect(w.pwrite(fd, 'x', -1n).errno).toBe(ERRNO.INVAL);
  });

  it("cuts and extends a file (`set_len`, `fd_allocate`); a user's bytes past a cut do not come back, and its Blob is not changed", async () => {
    const w = writableOf({ u: '0123456789' });
    const { fd } = w.openWith('u', { read: true, write: true });
    expect(w.wasi.fd_filestat_set_size!(fd, 4n)).toBe(ERRNO.SUCCESS);
    expect(w.contents('u')).toBe('0123');
    expect(w.wasi.fd_filestat_set_size!(fd, 8n)).toBe(ERRNO.SUCCESS);
    expect(w.contents('u')).toBe('0123\0\0\0\0');
    expect(w.wasi.fd_allocate!(fd, 6n, 6n)).toBe(ERRNO.SUCCESS);
    expect(w.stat('u').size).toBe(12n);
    expect(w.wasi.fd_allocate!(fd, 0n, 2n)).toBe(ERRNO.SUCCESS);
    expect(w.stat('u').size).toBe(12n);
    expect(await w.files[0]!.file.text()).toBe('0123456789');
    expect(writtenOf(w.directory).map(([name, bytes]) => [name, textOf(bytes)])).toEqual([['u', '0123\0\0\0\0\0\0\0\0']]);
  });

  it("writes into a large user's file in its view: the chunks that a write touches take the user's bytes of their range, the rest is read in place", () => {
    const size = 3 * CHUNK_BYTES;
    const big = new Uint8Array(size).map((_, i) => i % 251);
    const w = writableOf({ 'big.fa': big });
    const { fd } = w.openWith('big.fa', { read: true, write: true });
    w.seek(fd, BigInt(CHUNK_BYTES - 2), 0);
    expect(w.write(fd, 'WXYZ').count).toBe(4);
    // The two chunks around the write, each read from the Blob once.
    expect(w.reader.reads).toEqual([
      [0, CHUNK_BYTES],
      [CHUNK_BYTES, 2 * CHUNK_BYTES],
    ]);
    const expected = big.slice();
    expected.set(new TextEncoder().encode('WXYZ'), CHUNK_BYTES - 2);
    w.view().setUint32(3072, 4096, true);
    w.view().setUint32(3076, 8, true);
    expect(w.wasi.fd_pread!(fd, 3072, 1, BigInt(CHUNK_BYTES - 4), 2048)).toBe(ERRNO.SUCCESS);
    expect([...new Uint8Array(w.view().buffer, 4096, 8)]).toEqual([...expected.subarray(CHUNK_BYTES - 4, CHUNK_BYTES + 4)]);
    // A read of the third chunk takes the user's bytes in place.
    w.reader.reads.length = 0;
    expect(w.pread(fd, 3, BigInt(2 * CHUNK_BYTES)).errno).toBe(ERRNO.SUCCESS);
    expect(w.reader.reads).toEqual([[2 * CHUNK_BYTES, 2 * CHUNK_BYTES + 3]]);
    expect(w.directory.writtenBytes).toBe(2 * CHUNK_BYTES);
    const [[name, bytes]] = writtenOf(w.directory) as [[string, Uint8Array]];
    expect(name).toBe('big.fa');
    expect(Buffer.from(bytes).equals(Buffer.from(expected))).toBe(true);
    // The chunks were given up as the file went out.
    expect(w.directory.writtenBytes).toBe(0);
  });

  it('makes and removes folders, renames over a name, and removes files, answering as Linux does', () => {
    const w = writableOf({ 'old/db.pin': 'old', 'keep.txt': 'k' });
    expect(w.mkdir('out')).toBe(ERRNO.SUCCESS);
    expect(w.stat('out').type).toBe(FILETYPE_DIRECTORY);
    expect(w.mkdir('out')).toBe(ERRNO.EXIST);
    expect(w.mkdir('out/')).toBe(ERRNO.EXIST);
    expect(w.mkdir('no/x')).toBe(ERRNO.NOENT);
    expect(w.mkdir('keep.txt/x')).toBe(ERRNO.NOTDIR);
    expect(w.mkdir('../x')).toBe(ERRNO.NOTCAPABLE);
    // A volume written under its numbered name, then renamed to the database's name (makeblastdb's last step).
    const volume = w.openWith('out/db.00.pin', { write: true, create: true, truncate: true });
    w.write(volume.fd, 'v');
    expect(w.rename('out/db.00.pin', 'out/db.pin')).toBe(ERRNO.SUCCESS);
    expect(w.stat('out/db.00.pin').errno).toBe(ERRNO.NOENT);
    expect(w.contents('out/db.pin')).toBe('v');
    // A rename over a file replaces it.
    w.write(w.openWith('out/x', { write: true, create: true }).fd, 'X');
    expect(w.rename('out/x', 'out/db.pin')).toBe(ERRNO.SUCCESS);
    expect(w.contents('out/db.pin')).toBe('X');
    expect(w.rename('out/db.pin', 'out/db.pin')).toBe(ERRNO.SUCCESS);
    // A file onto a folder is EISDIR, a folder onto a file ENOTDIR, onto a folder that holds files ENOTEMPTY, into
    // itself EINVAL; `.` and the run directory EBUSY; a missing name ENOENT; a file to `x/` ENOTDIR.
    expect(w.rename('keep.txt', 'out')).toBe(ERRNO.ISDIR);
    expect(w.rename('old', 'keep.txt')).toBe(ERRNO.NOTDIR);
    expect(w.rename('old', 'out')).toBe(ERRNO.NOTEMPTY);
    expect(w.rename('old', 'old/sub')).toBe(ERRNO.INVAL);
    expect(w.rename('.', 'y')).toBe(ERRNO.BUSY);
    expect(w.rename('none', 'y')).toBe(ERRNO.NOENT);
    expect(w.rename('keep.txt', 'y/')).toBe(ERRNO.NOTDIR);
    expect(w.rename('keep.txt/', 'y')).toBe(ERRNO.NOTDIR);
    // unlink: a folder is EISDIR, a missing name ENOENT, `file/` ENOTDIR.
    expect(w.unlink('out')).toBe(ERRNO.ISDIR);
    expect(w.unlink('none')).toBe(ERRNO.NOENT);
    expect(w.unlink('keep.txt/')).toBe(ERRNO.NOTDIR);
    // rmdir: a folder that holds a file is ENOTEMPTY, a file ENOTDIR, `.` EINVAL.
    expect(w.rmdir('out')).toBe(ERRNO.NOTEMPTY);
    expect(w.rmdir('keep.txt')).toBe(ERRNO.NOTDIR);
    expect(w.rmdir('.')).toBe(ERRNO.INVAL);
    expect(w.rmdir('none')).toBe(ERRNO.NOENT);
    // The user's old database goes, and its folder.
    expect(w.unlink('old/db.pin')).toBe(ERRNO.SUCCESS);
    expect(w.rmdir('old/')).toBe(ERRNO.SUCCESS);
    expect(w.stat('old').errno).toBe(ERRNO.NOENT);
    // A folder that is renamed takes its files with it; an empty one onto an empty one replaces it.
    expect(w.mkdir('tmp')).toBe(ERRNO.SUCCESS);
    w.write(w.openWith('tmp/a', { write: true, create: true }).fd, 'A');
    expect(w.mkdir('final')).toBe(ERRNO.SUCCESS);
    expect(w.rename('tmp', 'final')).toBe(ERRNO.SUCCESS);
    expect(w.contents('final/a')).toBe('A');
    expect(w.stat('tmp').errno).toBe(ERRNO.NOENT);
    // The files of the run that are there, by name; the user's file that the run did not change is not one.
    expect(writtenOf(w.directory).map(([name, bytes]) => [name, textOf(bytes)])).toEqual([
      ['final/a', 'A'],
      ['out/db.pin', 'X'],
    ]);
  });

  it('keeps a removed file readable and writable through its descriptors, and leaves it out of what the run leaves', () => {
    const w = writableOf({});
    const lock = w.openWith('db.pdb-lock', { read: true, write: true, create: true });
    w.write(lock.fd, 'L');
    const reading = w.openWith('db.pdb-lock', { read: true });
    expect(w.unlink('db.pdb-lock')).toBe(ERRNO.SUCCESS);
    expect(w.stat('db.pdb-lock').errno).toBe(ERRNO.NOENT);
    expect(w.read(reading.fd, 10).bytes).toBe('L');
    expect(w.fstat(reading.fd).nlink).toBe(0n);
    expect(w.write(lock.fd, 'M').count).toBe(1);
    expect(w.pread(reading.fd, 10, 0n).bytes).toBe('LM');
    expect(w.directory.writtenBytes).toBeGreaterThan(0);
    w.wasi.fd_close!(lock.fd);
    w.wasi.fd_close!(reading.fd);
    // Its bytes go with its last descriptor.
    expect(w.directory.writtenBytes).toBe(0);
    expect(writtenOf(w.directory)).toEqual([]);
  });

  it('lists a folder as it was when the listing began: names that the run removes meanwhile come once', () => {
    const w = writableOf({ 'db/a': 'a', 'db/b': 'b', 'db/c': 'c' });
    w.write(w.openWith('db/d', { write: true, create: true }).fd, 'd');
    const folder = w.openWith('db', { read: true }).fd;
    const first = w.list(folder, 60);
    expect(first.names).toEqual(['.', '..']);
    expect(w.unlink('db/a')).toBe(ERRNO.SUCCESS);
    expect(w.list(folder, 4096, first.next).names).toEqual(['a', 'b', 'c', 'd']);
    expect(w.list(folder, 4096).names).toEqual(['.', '..', 'b', 'c', 'd']);
    expect(w.list(w.openWith('.', { read: true }).fd, 4096).names).toEqual(['.', '..', 'db']);
  });

  it("gives the files that the run made or changed by their names' bytes, a user's file renamed at its new name", () => {
    const w = writableOf({ B: 'b', a: 'a', z: 'zz' });
    expect(w.mkdir('d')).toBe(ERRNO.SUCCESS);
    w.write(w.openWith('d/x', { write: true, create: true }).fd, 'x');
    w.write(w.openWith('a.b', { write: true, create: true }).fd, 'ab');
    expect(w.rename('z', 'é')).toBe(ERRNO.SUCCESS);
    w.write(w.openWith('B', { write: true, append: true }).fd, '+');
    // An empty folder that the run made is not a file of it.
    expect(w.mkdir('empty')).toBe(ERRNO.SUCCESS);
    expect(writtenOf(w.directory).map(([name, bytes]) => [name, textOf(bytes)])).toEqual([
      ['B', 'b+'],
      ['a.b', 'ab'],
      ['d/x', 'x'],
      ['é', 'zz'],
    ]);
    // The next run's files replace what this run made.
    w.directory.open([], { writable: true });
    expect(w.stat('d').errno).toBe(ERRNO.NOENT);
    expect(w.directory.takeWritten()).toEqual([]);
  });

  it("is writable only on an instance of its own table; a search's directory refuses every change with EROFS", () => {
    expect(() => new RunDirectory(new SharedArrayBuffer(TABLE_BYTES), () => new Uint8Array(0)).open([], { writable: true })).toThrow(/serial instance/);
    const { directory } = directoryOf({ 'db/nt.nin': 'x', a: 'a' });
    const g = guest(directory);
    const [from, fromLength] = g.path('a');
    expect(g.wasi.path_rename!(PREOPEN_FD, from, fromLength, PREOPEN_FD, from, fromLength)).toBe(ERRNO.ROFS);
    expect(g.wasi.path_remove_directory!(PREOPEN_FD, ...g.path('db'))).toBe(ERRNO.ROFS);
    expect(g.wasi.path_filestat_set_times!(PREOPEN_FD, 0, ...g.path('a'), 0n, 0n, 8)).toBe(ERRNO.ROFS);
    const { fd } = g.open('a');
    expect(g.wasi.fd_filestat_set_times!(fd, 0n, 0n, 8)).toBe(ERRNO.ROFS);
    expect(g.wasi.fd_pwrite!(fd, 0, 0, 0n, 2048)).toBe(ERRNO.BADF);
    expect(g.wasi.fd_filestat_set_size!(fd, 0n)).toBe(ERRNO.INVAL);
    expect(g.wasi.fd_allocate!(fd, 0n, 1n)).toBe(ERRNO.BADF);
    expect(directory.takeWritten()).toEqual([]);
  });
});
