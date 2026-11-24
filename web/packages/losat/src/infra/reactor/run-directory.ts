// The run's working directory as the engine sees it (docs/web/abi_v2.md §3; line WQA D12a, 2026-11-20): WASI's
// preopened directory `/work` (descriptor 3, `RUN_DIRECTORY`), read-only, that holds the files that the argv names,
// as the user gave them (`DirectoryFile`: a BLAST database of `-db` with its volumes, `.nal`/`.pal`, `taxdb.*`, the
// lists of the filters). The engine makes it its current directory and opens the files by name, as the command line
// opens them in its working directory: `db/nt` is the file `db/nt` of the run, and an absolute name such as `/db/nt`
// names a file outside it, where there is none (WASI's C library gives a path to the preopen whose name begins it,
// so under a preopen `.` the names `db/nt` and `/db/nt` were one file; review of D12a, low 4). A path is looked up
// one name after another, as the OS looks it up: a name after a file is ENOTDIR and a name after a missing one
// ENOENT, `..` included. A folder of the directory opens for reading as it does on Linux, where the open of a
// directory succeeds and its reads fail: the descriptor is a directory's (`fd_read` and `fd_pread` give EISDIR, the stat
// calls say directory, `fd_readdir` lists it, and a path opened from it is looked up in it); an open of a folder
// for writing is EISDIR (line NP WW3a, 2026-11-23: the command line's `windowmasker -in <folder>` says "Is a directory"
// from the read). The engine reads the files in place: each read takes the bytes that it asks for from the
// file's Blob (`File.slice`) with `FileReaderSync`, so a file is never copied whole into the Wasm memory or into the
// browser's storage. The Node hosts of LOSAT's tests give the engine the same view from the OS
// (LOSAT/tests/wasi_thread_worker.js, `preopens: { "/work": directory }`).
//
// LOSAT opens a database file when a search needs it (once for a run; the repeat database of BLASTN's -filtering_db at
// each query batch) (LOSAT/src/seqdb/file.rs: on WASI, a seek and a read under a lock) and
// reads it from the search's threads. A descriptor that one thread opened must therefore be read on another: each
// thread worker has a WASI shim of its own, so the descriptors live in a table in memory that every thread of the instance shares (`table`, a
// SharedArrayBuffer for the threaded module), each entry the index of the file in the run's list (for a
// folder, its place in the run's list of folders) and its position.
// The list itself is the same on every thread: the Engine worker gives it to each thread that the run starts
// (thread-host.ts, `start`), and a File crosses to a worker without its bytes.
import type { DirectoryFile } from '../../domain/run-files.js';

/** The descriptor of the preopened run directory, the first that WASI's C library asks about. */
export const PREOPEN_FD = 3;
/** The run directory's name: the preopen's, which the engine makes its current directory (abi_v2.md §3). */
export const RUN_DIRECTORY = '/work';
/** The first descriptor of an open file of the directory: far from the shim's own (0-2 and the few it opens). */
export const FILE_FD_BASE = 0x40000000;
/** How many files of the directory may be open at once. */
export const OPEN_FILES = 1024;
/**
 * The bytes of the table of descriptors (`RunDirectory.table`): per slot the file (Int32) and the position (BigInt64). A
 * slot is 0 when free, the file's index plus 1 for a file, and minus the folder's place plus 1 for a folder.
 */
export const TABLE_BYTES = 12 * OPEN_FILES;
/** The most bytes that one read takes from a Blob, so that a large read does not hold a second copy of its bytes. */
const READ_LIMIT = 16 * 1024 * 1024;

// WASI preview1's values (wasi_snapshot_preview1.witx).
const ERRNO_SUCCESS = 0;
const ERRNO_BADF = 8;
const ERRNO_EXIST = 20;
const ERRNO_INVAL = 28;
const ERRNO_IO = 29;
const ERRNO_ISDIR = 31;
const ERRNO_MFILE = 33;
const ERRNO_NAMETOOLONG = 37;
const ERRNO_NOENT = 44;
const ERRNO_NOTDIR = 54;
const ERRNO_ROFS = 69;
const ERRNO_NOTCAPABLE = 76;
const FILETYPE_DIRECTORY = 3;
const FILETYPE_REGULAR_FILE = 4;
const OFLAGS_CREAT = 1;
const OFLAGS_DIRECTORY = 2;
const OFLAGS_EXCL = 4;
const OFLAGS_TRUNC = 8;
const FDFLAGS_APPEND = 1;
const RIGHT_FD_WRITE = 1n << 6n;
/** The rights of a file of the directory: every right but those that write (the directory is read-only). */
const FILE_RIGHTS =
  (1n << 1n) | // fd_read
  (1n << 2n) | // fd_seek
  (1n << 3n) | // fd_fdstat_set_flags
  (1n << 4n) | // fd_sync
  (1n << 5n) | // fd_tell
  (1n << 7n) | // fd_advise
  (1n << 21n); // fd_filestat_get
/** The rights of the directory: to look up, open and stat its files, and those of the files it opens. */
const DIRECTORY_RIGHTS =
  (1n << 13n) | // path_open
  (1n << 14n) | // fd_readdir
  (1n << 15n) | // path_readlink
  (1n << 18n) | // path_filestat_get
  (1n << 21n); // fd_filestat_get

/** The order of two names as bytes. */
function compareBytes(a: Uint8Array, b: Uint8Array): number {
  for (let i = 0; i < Math.min(a.length, b.length); i++) if (a[i] !== b[i]) return a[i]! - b[i]!;
  return a.length - b.length;
}

/** Reads the bytes `start` to `end` of a Blob, synchronously (on a worker: `FileReaderSync`). */
export type BlobReader = (blob: Blob, start: number, end: number) => Uint8Array;

interface FileReaderSyncLike {
  readAsArrayBuffer(blob: Blob): ArrayBuffer;
}

/**
 * The browser's synchronous read of a Blob, which dedicated workers have (the Engine worker and the thread workers
 * of every browser); undefined elsewhere (the main thread, Node).
 */
export function fileReaderSync(): BlobReader | undefined {
  const Reader = (globalThis as { FileReaderSync?: new () => FileReaderSyncLike }).FileReaderSync;
  if (Reader === undefined) return undefined;
  const reader = new Reader();
  return (blob, start, end) => new Uint8Array(reader.readAsArrayBuffer(blob.slice(start, end)));
}

/**
 * The path `path` of the directory as its names (`a/./b//c` is `a`, `b`, `c`; `..` takes the name before it), or
 * undefined where `..` leaves the directory, which WASI does not allow (ENOTCAPABLE).
 */
export function pathNames(path: string): string[] | undefined {
  const names: string[] = [];
  for (const name of path.split('/')) {
    if (name === '' || name === '.') continue;
    if (name === '..') {
      if (names.pop() === undefined) return undefined;
      continue;
    }
    names.push(name);
  }
  return names;
}

/** What the WASI functions of the directory need: the guest's memory, as this thread sees it. */
export type GuestMemory = () => WebAssembly.Memory;

type WasiFunction = (...args: Array<number | bigint>) => unknown;
type WasiImports = Record<string, WasiFunction>;

/**
 * The run's working directory of one instance of the engine on one thread. The Engine worker opens the run's files
 * in it (`open`), which also ends the descriptors of an earlier run; a thread worker shares the table and takes the
 * same list for each thread that it runs (`take`). Its WASI functions (`imports`) stand in front of the shim's.
 */
export class RunDirectory {
  private files: readonly DirectoryFile[] = [];
  /** The files by their path's names joined with `/`. */
  private byPath = new Map<string, number>();
  /** The directories that hold a file: every leading part of a file's path. */
  private directories = new Set<string>(['']);
  /** The folders in the order of `directories` (the run directory first): a folder's descriptor holds its place here. */
  private folders: readonly string[] = [''];
  private readonly slots: Int32Array;
  private readonly positions: BigInt64Array;

  /**
   * @param table the table of descriptors (`TABLE_BYTES`), a SharedArrayBuffer that every thread of the threaded
   *   instance shares, or an ArrayBuffer of the serial instance's own
   * @param read the read of a Blob (`fileReaderSync` in a worker); undefined where there is none, and a read then
   *   fails with EIO
   */
  constructor(
    readonly table: ArrayBuffer | SharedArrayBuffer = new ArrayBuffer(TABLE_BYTES),
    private readonly read: BlobReader | undefined = fileReaderSync(),
  ) {
    this.slots = new Int32Array(table, 0, OPEN_FILES);
    this.positions = new BigInt64Array(table, 4 * OPEN_FILES, OPEN_FILES);
  }

  /** The files that the directory holds. */
  get list(): readonly DirectoryFile[] {
    return this.files;
  }

  /**
   * Gives the directory the files of a run, and ends the descriptors that an earlier run left open (the Engine
   * worker, before a run's first call into the instance, while no thread of the instance runs). Two files of one
   * path are an error.
   */
  open(files: readonly DirectoryFile[]): void {
    this.take(files);
    this.slots.fill(0);
    this.positions.fill(0n);
  }

  /** Gives the directory the files of a run that the Engine worker opened (`open`): a thread worker's. */
  take(files: readonly DirectoryFile[]): void {
    const byPath = new Map<string, number>();
    const directories = new Set<string>(['']);
    files.forEach((file, index) => {
      const names = pathNames(file.name);
      if (names === undefined || names.length === 0) throw new Error(`the file name ${JSON.stringify(file.name)} is not a path in the run's directory`);
      const path = names.join('/');
      if (byPath.has(path)) throw new Error(`two files are named ${path}`);
      byPath.set(path, index);
      for (let i = 1; i < names.length; i++) directories.add(names.slice(0, i).join('/'));
    });
    for (const path of byPath.keys()) {
      if (directories.has(path)) throw new Error(`${path} is both a file and a folder`);
    }
    this.files = files;
    this.byPath = byPath;
    this.directories = directories;
    this.folders = [...directories];
  }

  /**
   * The WASI functions of the directory in front of `base` (the shim's, which take every other descriptor): the
   * preopened run directory, and the files that it opens.
   */
  imports(base: WebAssembly.ModuleImports, memory: GuestMemory): WebAssembly.ModuleImports {
    const fallback = base as WasiImports;
    const own: WasiImports = {};
    const bytes = (pointer: number, length: number): Uint8Array => {
      const buffer = memory().buffer;
      // A shared memory that another thread grew: grow(0) gives this thread its new length.
      if (pointer + length > buffer.byteLength) memory().grow(0);
      return new Uint8Array(memory().buffer, pointer, length);
    };
    const view = (pointer: number, length: number) => {
      const target = bytes(pointer, length);
      return new DataView(target.buffer, target.byteOffset, target.byteLength);
    };
    const pathAt = (pointer: number, length: number) => new TextDecoder().decode(bytes(pointer, length).slice());
    /**
     * The file (its index) or directory (-1) of a path of the directory, or an errno. Each name is looked up in the
     * folder that the names before it give, as the OS looks a path up: a name (`..` too) after a file is ENOTDIR,
     * after a name that is not there ENOENT; `..` past the directory is ENOTCAPABLE, as WASI gives it.
     */
    const lookup = (path: string, base = ''): { readonly file: number; readonly folder?: string } | { readonly errno: number } => {
      const names: string[] = base === '' ? [] : base.split('/');
      for (const name of path.split('/')) {
        const here = names.join('/');
        if (names.length > 0 && !this.directories.has(here)) return { errno: this.byPath.has(here) ? ERRNO_NOTDIR : ERRNO_NOENT };
        if (name === '' || name === '.') continue;
        if (name !== '..') names.push(name);
        else if (names.pop() === undefined) return { errno: ERRNO_NOTCAPABLE };
      }
      const joined = names.join('/');
      const file = this.byPath.get(joined);
      if (file !== undefined) return { file };
      if (this.directories.has(joined)) return { file: -1, folder: joined };
      return { errno: ERRNO_NOENT };
    };
    const slotOf = (fd: number) => fd - FILE_FD_BASE;
    const isFile = (fd: number) => slotOf(fd) >= 0 && slotOf(fd) < OPEN_FILES;
    /** The file of an open descriptor, or undefined (a free slot, and a folder's). */
    const openFile = (fd: number): DirectoryFile | undefined => {
      const index = Atomics.load(this.slots, slotOf(fd)) - 1;
      return index < 0 ? undefined : this.files[index];
    };
    /** The folder (its path, '' for the run directory) that an open descriptor is, or undefined. */
    const openFolder = (fd: number): string | undefined => {
      const state = Atomics.load(this.slots, slotOf(fd));
      return state < 0 ? this.folders[-state - 1] : undefined;
    };
    /** The folder that a path is looked up from: the preopen's, or an open folder's; undefined for any other descriptor. */
    const baseOf = (fd: number): string | undefined => (fd === PREOPEN_FD ? '' : isFile(fd) ? openFolder(fd) : undefined);
    const writeFilestat = (pointer: number, file: number) => {
      const out = view(pointer, 64);
      const blob = file < 0 ? undefined : this.files[file]?.file;
      const time = BigInt(Math.max(0, Math.trunc((blob as File | undefined)?.lastModified ?? 0))) * 1_000_000n;
      out.setBigUint64(0, 0n, true); // dev
      out.setBigUint64(8, BigInt(file + 2), true); // ino: the directory 1, the files 2 on
      out.setBigUint64(16, 0n, true);
      out.setUint8(16, file < 0 ? FILETYPE_DIRECTORY : FILETYPE_REGULAR_FILE);
      out.setBigUint64(24, 1n, true); // nlink
      out.setBigUint64(32, BigInt(blob?.size ?? 0), true);
      for (const at of [40, 48, 56]) out.setBigUint64(at, time, true);
    };
    /** Reads up to the length of the vectors of `iovs` from `file` at `offset` into them; returns the count. */
    const readInto = (file: DirectoryFile, iovs: number, count: number, offset: bigint): number => {
      const size = BigInt(file.file.size);
      if (offset >= size) return 0;
      const vectors = view(iovs, 8 * count);
      let wanted = 0;
      for (let i = 0; i < count; i++) wanted += vectors.getUint32(8 * i + 4, true);
      const start = Number(offset);
      const length = Math.min(wanted, Number(size - offset), READ_LIMIT);
      if (length === 0) return 0;
      if (this.read === undefined) throw new Error('this thread cannot read a file synchronously');
      const data = this.read(file.file, start, start + length);
      let done = 0;
      for (let i = 0; i < count && done < data.length; i++) {
        const pointer = vectors.getUint32(8 * i, true);
        const part = Math.min(vectors.getUint32(8 * i + 4, true), data.length - done);
        bytes(pointer, part).set(data.subarray(done, done + part));
        done += part;
      }
      return done;
    };
    /**
     * Calls `body` for a descriptor of the directory's files; errors of a read are EIO. A descriptor of a folder gets
     * `onFolder`: an errno (what Linux answers for that call on a directory descriptor) or a function of its own.
     */
    const onFile = (
      name: string,
      body: (fd: number, file: DirectoryFile, ...args: Array<number | bigint>) => number,
      onFolder: number | ((fd: number, folder: string, ...args: Array<number | bigint>) => number) = ERRNO_BADF,
    ) => {
      own[name] = (fd, ...args) => {
        const descriptor = Number(fd);
        if (!isFile(descriptor)) return fallback[name]!(fd, ...args);
        const folder = openFolder(descriptor);
        if (folder !== undefined) return typeof onFolder === 'number' ? onFolder : onFolder(descriptor, folder, ...args);
        const file = openFile(descriptor);
        if (file === undefined) return ERRNO_BADF;
        try {
          return body(descriptor, file, ...args);
        } catch {
          return ERRNO_IO;
        }
      };
    };

    own.fd_prestat_get = (fd, pointer) => {
      if (Number(fd) !== PREOPEN_FD) return fallback.fd_prestat_get!(fd, pointer);
      const out = view(Number(pointer), 8);
      out.setUint8(0, 0); // a directory
      out.setUint32(4, RUN_DIRECTORY.length, true); // the length of its name
      return ERRNO_SUCCESS;
    };
    own.fd_prestat_dir_name = (fd, pointer, length) => {
      if (Number(fd) !== PREOPEN_FD) return fallback.fd_prestat_dir_name!(fd, pointer, length);
      if (Number(length) < RUN_DIRECTORY.length) return ERRNO_NAMETOOLONG;
      bytes(Number(pointer), RUN_DIRECTORY.length).set(new TextEncoder().encode(RUN_DIRECTORY));
      return ERRNO_SUCCESS;
    };
    own.path_open = (dirfd, dirflags, path, pathLength, oflags, rightsBase, rightsInheriting, fdflags, opened) => {
      const base = baseOf(Number(dirfd));
      if (base === undefined) {
        if (isFile(Number(dirfd))) return Atomics.load(this.slots, slotOf(Number(dirfd))) === 0 ? ERRNO_BADF : ERRNO_NOTDIR;
        return fallback.path_open!(dirfd, dirflags, path, pathLength, oflags, rightsBase, rightsInheriting, fdflags, opened);
      }
      const found = lookup(pathAt(Number(path), Number(pathLength)), base);
      const flags = Number(oflags);
      const writes = (BigInt(rightsBase) & RIGHT_FD_WRITE) !== 0n || (Number(fdflags) & FDFLAGS_APPEND) !== 0;
      if ('errno' in found) return found.errno === ERRNO_NOENT && flags & OFLAGS_CREAT ? ERRNO_ROFS : found.errno;
      if (flags & OFLAGS_CREAT && flags & OFLAGS_EXCL) return ERRNO_EXIST;
      let state: number;
      if (found.file < 0) {
        // As Linux: a folder opens for reading (its reads fail with EISDIR); for writing, with O_CREAT or O_TRUNC it is EISDIR.
        if (writes || flags & (OFLAGS_CREAT | OFLAGS_TRUNC)) return ERRNO_ISDIR;
        state = -(this.folders.indexOf(found.folder!) + 1);
      } else {
        if (flags & OFLAGS_DIRECTORY) return ERRNO_NOTDIR;
        if (writes || flags & OFLAGS_TRUNC) return ERRNO_ROFS;
        state = found.file + 1;
      }
      for (let slot = 0; slot < OPEN_FILES; slot++) {
        if (Atomics.compareExchange(this.slots, slot, 0, state) !== 0) continue;
        Atomics.store(this.positions, slot, 0n);
        view(Number(opened), 4).setUint32(0, FILE_FD_BASE + slot, true);
        return ERRNO_SUCCESS;
      }
      return ERRNO_MFILE;
    };
    own.path_filestat_get = (dirfd, flags, path, pathLength, pointer) => {
      const base = baseOf(Number(dirfd));
      if (base === undefined) {
        if (isFile(Number(dirfd))) return Atomics.load(this.slots, slotOf(Number(dirfd))) === 0 ? ERRNO_BADF : ERRNO_NOTDIR;
        return fallback.path_filestat_get!(dirfd, flags, path, pathLength, pointer);
      }
      const found = lookup(pathAt(Number(path), Number(pathLength)), base);
      if ('errno' in found) return found.errno;
      writeFilestat(Number(pointer), found.file);
      return ERRNO_SUCCESS;
    };
    own.path_readlink = (dirfd, path, pathLength, buffer, length, used) => {
      const base = baseOf(Number(dirfd));
      if (base === undefined) return fallback.path_readlink!(dirfd, path, pathLength, buffer, length, used);
      const found = lookup(pathAt(Number(path), Number(pathLength)), base);
      // No file of the directory is a symbolic link.
      return 'errno' in found ? found.errno : ERRNO_INVAL;
    };
    for (const name of ['path_create_directory', 'path_unlink_file', 'path_remove_directory']) {
      own[name] = (dirfd, ...args) => (baseOf(Number(dirfd)) !== undefined ? ERRNO_ROFS : fallback[name]!(dirfd, ...args));
    }
    own.fd_filestat_get = (fd, pointer) => {
      if (Number(fd) === PREOPEN_FD) {
        writeFilestat(Number(pointer), -1);
        return ERRNO_SUCCESS;
      }
      if (!isFile(Number(fd))) return fallback.fd_filestat_get!(fd, pointer);
      const state = Atomics.load(this.slots, slotOf(Number(fd)));
      if (state === 0) return ERRNO_BADF;
      writeFilestat(Number(pointer), state < 0 ? -1 : state - 1);
      return ERRNO_SUCCESS;
    };
    own.fd_fdstat_get = (fd, pointer) => {
      const folder = Number(fd) === PREOPEN_FD || (isFile(Number(fd)) && openFolder(Number(fd)) !== undefined);
      if (!folder && !isFile(Number(fd))) return fallback.fd_fdstat_get!(fd, pointer);
      if (!folder && openFile(Number(fd)) === undefined) return ERRNO_BADF;
      const out = view(Number(pointer), 24);
      out.setBigUint64(0, 0n, true);
      out.setUint8(0, folder ? FILETYPE_DIRECTORY : FILETYPE_REGULAR_FILE);
      out.setBigUint64(8, folder ? DIRECTORY_RIGHTS : FILE_RIGHTS, true);
      out.setBigUint64(16, folder ? FILE_RIGHTS | DIRECTORY_RIGHTS : 0n, true);
      return ERRNO_SUCCESS;
    };
    // A read of a folder's descriptor is EISDIR, as on Linux, where the open succeeded.
    onFile(
      'fd_read',
      (fd, file, iovs, count, nread) => {
        const slot = slotOf(fd);
        const got = readInto(file, Number(iovs), Number(count), Atomics.load(this.positions, slot));
        Atomics.add(this.positions, slot, BigInt(got));
        view(Number(nread), 4).setUint32(0, got, true);
        return ERRNO_SUCCESS;
      },
      ERRNO_ISDIR,
    );
    onFile(
      'fd_pread',
      (_fd, file, iovs, count, offset, nread) => {
        const got = readInto(file, Number(iovs), Number(count), BigInt.asUintN(64, BigInt(offset!)));
        view(Number(nread), 4).setUint32(0, got, true);
        return ERRNO_SUCCESS;
      },
      ERRNO_ISDIR,
    );
    /** The seek of a descriptor whose end is `size` bytes (a folder's is 0). */
    const seek = (fd: number, size: number, offset: number | bigint | undefined, whence: number | bigint | undefined, result: number | bigint | undefined) => {
      const slot = slotOf(fd);
      const from = [0n, Atomics.load(this.positions, slot), BigInt(size)][Number(whence)];
      const to = from === undefined ? -1n : from + BigInt.asIntN(64, BigInt(offset!));
      if (to < 0n) return ERRNO_INVAL;
      Atomics.store(this.positions, slot, to);
      view(Number(result), 8).setBigUint64(0, to, true);
      return ERRNO_SUCCESS;
    };
    onFile('fd_seek', (fd, file, offset, whence, result) => seek(fd, file.file.size, offset, whence, result), (fd, _folder, offset, whence, result) =>
      seek(fd, 0, offset, whence, result),
    );
    const tell = (fd: number, result: number | bigint | undefined) => {
      view(Number(result), 8).setBigUint64(0, Atomics.load(this.positions, slotOf(fd)), true);
      return ERRNO_SUCCESS;
    };
    onFile('fd_tell', (fd, _file, result) => tell(fd, result), (fd, _folder, result) => tell(fd, result));
    for (const name of ['fd_advise', 'fd_sync', 'fd_datasync']) onFile(name, () => ERRNO_SUCCESS, ERRNO_SUCCESS);
    const setFlags = (flags: number | bigint | undefined) => (Number(flags) === 0 ? ERRNO_SUCCESS : ERRNO_NOTCAPABLE);
    onFile('fd_fdstat_set_flags', (_fd, _file, flags) => setFlags(flags), (_fd, _folder, flags) => setFlags(flags));
    for (const name of ['fd_write', 'fd_pwrite', 'fd_allocate', 'fd_filestat_set_size', 'fd_filestat_set_times']) {
      // A file or folder of the read-only directory is not open for writing.
      onFile(name, () => ERRNO_BADF);
    }
    onFile('fd_readdir', () => ERRNO_NOTDIR, (_fd, folder, buffer, length, cookie, used) => {
      // The folder's entries, as the OS lists them: `.` and `..`, then the names in the folder, in byte order.
      const prefix = folder === '' ? '' : `${folder}/`;
      const children = new Map<string, number>();
      for (const path of this.byPath.keys()) {
        if (!path.startsWith(prefix)) continue;
        const name = path.slice(prefix.length);
        if (!name.includes('/')) children.set(name, FILETYPE_REGULAR_FILE);
      }
      for (const path of this.folders) {
        if (path === '' || !path.startsWith(prefix)) continue;
        const name = path.slice(prefix.length);
        if (!name.includes('/')) children.set(name, FILETYPE_DIRECTORY);
      }
      const encoder = new TextEncoder();
      const names = [...children.keys()].map((name) => ({ name, bytes: encoder.encode(name) })).sort((a, b) => compareBytes(a.bytes, b.bytes));
      const entries = [
        { bytes: encoder.encode('.'), type: FILETYPE_DIRECTORY },
        { bytes: encoder.encode('..'), type: FILETYPE_DIRECTORY },
        ...names.map(({ name, bytes }) => ({ bytes, type: children.get(name)! })),
      ];
      const room = Number(length);
      const target = bytes(Number(buffer), room);
      let written = 0;
      for (let next = Number(BigInt.asUintN(32, BigInt(cookie!))); next < entries.length && written < room; next++) {
        const entry = entries[next]!;
        const record = new Uint8Array(24 + entry.bytes.length);
        const header = new DataView(record.buffer);
        header.setBigUint64(0, BigInt(next + 1), true); // d_next: the cookie of the entry after this
        header.setBigUint64(8, 1n, true); // d_ino
        header.setUint32(16, entry.bytes.length, true);
        header.setUint8(20, entry.type);
        record.set(entry.bytes, 24);
        const taken = Math.min(record.length, room - written);
        target.set(record.subarray(0, taken), written);
        written += taken;
      }
      view(Number(used), 4).setUint32(0, written, true);
      return ERRNO_SUCCESS;
    });
    own.fd_close = (fd) => {
      const descriptor = Number(fd);
      if (!isFile(descriptor)) return fallback.fd_close!(fd);
      const slot = slotOf(descriptor);
      const state = Atomics.load(this.slots, slot);
      if (state === 0 || Atomics.compareExchange(this.slots, slot, state, 0) !== state) return ERRNO_BADF;
      return ERRNO_SUCCESS;
    };
    const merged: WasiImports = { ...fallback };
    // Only the functions that the module imports from the shim are replaced; the others stay out.
    for (const [name, call] of Object.entries(own)) if (name in fallback) merged[name] = call;
    return merged as WebAssembly.ModuleImports;
  }
}
