// The run's working directory of the browser host (src/infra/reactor/run-directory.ts) through its WASI functions, as
// the engine calls them: the preopen `.`, opening a file by its path, reading it in place, a descriptor that one
// thread opened read on another, and a directory that writes nothing (line WQA D12a, 2026-11-20).
import { describe, expect, it } from 'vitest';
import { FILE_FD_BASE, pathNames, PREOPEN_FD, RUN_DIRECTORY, RunDirectory, TABLE_BYTES, type BlobReader } from '../src/infra/reactor/run-directory.ts';

const ERRNO = { SUCCESS: 0, BADF: 8, INVAL: 28, IO: 29, ISDIR: 31, NAMETOOLONG: 37, NOENT: 44, NOTDIR: 54, ROFS: 69, NOTCAPABLE: 76 };
const RIGHT_FD_READ = 1n << 1n;
const RIGHT_FD_WRITE = 1n << 6n;
const OFLAGS_CREAT = 1;
const OFLAGS_TRUNC = 8;

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
    open(name: string, oflags = 0, rights = RIGHT_FD_READ): { errno: number; fd: number } {
      const errno = wasi.path_open!(PREOPEN_FD, 0, ...path(name), oflags, rights, 0n, 0, OUT);
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
    expect(g.open('db').errno).toBe(ERRNO.ISDIR);
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
