// The run's working directory as the engine sees it (docs/web/abi_v2.md §3; line WQA D12a, 2026-11-20): WASI's
// preopened directory `/work` (descriptor 3, `RUN_DIRECTORY`), that holds the files that the argv names, as the user
// gave them (`DirectoryFile`: a BLAST database of `-db` with its volumes, `.nal`/`.pal`, `taxdb.*`, the lists of the
// filters, a tool's `-in`). The engine makes it its current directory and opens the files by name, as the command line
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
// A search's directory is read-only: a call that would make, change or remove a name is EROFS, and a file opens for
// reading only (a search writes no named file on the Web: its outputs are streams). A tool's directory is writable
// (`open(files, { writable: true })`; line NP WW4b1, D-WW4b-1, 2026-11-27): makeblastdb writes its database there as the
// command line writes it into its working directory, with `std::fs`. The directory is then an overlay: the user's files
// below (their Blobs, never changed), and above them what the run makes and writes, in this worker's memory, in chunks
// of at most `CHUNK_BYTES` (a chunk holds the bytes of its range; a range of no chunk past the user's bytes is a hole
// and reads as zeros). The run creates, truncates, appends to, reads back, seeks in, cuts and extends its files, makes,
// renames (over another name) and removes files and folders, and lists them, and each call answers as Linux answers
// it. A file that the run removed stays readable through the descriptors that it has open, as on Linux. A rename or a
// write of a user's file changes the run's view only. After the run, the files that the run made or changed and that
// are still there are its output (`takeWritten`, by name in byte order): what the command line leaves in its folder.
// A file that it removed (a lock, a volume's name before its rename, what a failed run cleans up) is not there. Writes
// are possible only on the serial instance, whose table is its own: tool runs are serial (runtime.ts).
//
// LOSAT opens a database file when a search needs it (once for a run; the repeat database of BLASTN's -filtering_db at
// each query batch) (LOSAT/src/seqdb/file.rs: on WASI, a seek and a read under a lock) and
// reads it from the search's threads. A descriptor that one thread opened must therefore be read on another: each
// thread worker has a WASI shim of its own, so the descriptors live in a table in memory that every thread of the instance shares (`table`, a
// SharedArrayBuffer for the threaded module), each entry the number of the file in the run's list (for a
// folder, its place in the run's list of folders) and its position.
// The list itself is the same on every thread: the Engine worker gives it to each thread that the run starts
// (thread-host.ts, `start`), and a File crosses to a worker without its bytes; each thread builds the same files and
// folders from it, in the same order.
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
 * slot is 0 when free, the file's number plus 1 for a file, and minus the folder's number plus 1 for a folder.
 */
export const TABLE_BYTES = 12 * OPEN_FILES;
/** The most bytes that one read takes from a Blob, so that a large read does not hold a second copy of its bytes. */
const READ_LIMIT = 16 * 1024 * 1024;
/** The bytes of a chunk of what a run writes (a frame of a tool's file holds one chunk at most). */
export const CHUNK_BYTES = 1024 * 1024;
/** The smallest chunk: a small file does not take a whole chunk. */
const MIN_CHUNK_BYTES = 4096;

// WASI preview1's values (wasi_snapshot_preview1.witx).
const ERRNO_SUCCESS = 0;
const ERRNO_BADF = 8;
const ERRNO_BUSY = 10;
const ERRNO_EXIST = 20;
const ERRNO_FBIG = 22;
const ERRNO_INVAL = 28;
const ERRNO_IO = 29;
const ERRNO_ISDIR = 31;
const ERRNO_MFILE = 33;
const ERRNO_NAMETOOLONG = 37;
const ERRNO_NOENT = 44;
const ERRNO_NOSPC = 51;
const ERRNO_NOTDIR = 54;
const ERRNO_NOTEMPTY = 55;
const ERRNO_ROFS = 69;
const ERRNO_XDEV = 75;
const ERRNO_NOTCAPABLE = 76;
const FILETYPE_DIRECTORY = 3;
const FILETYPE_REGULAR_FILE = 4;
const OFLAGS_CREAT = 1;
const OFLAGS_DIRECTORY = 2;
const OFLAGS_EXCL = 4;
const OFLAGS_TRUNC = 8;
const FDFLAGS_APPEND = 1;
const FSTFLAGS_MTIM = 4;
const FSTFLAGS_MTIM_NOW = 8;
const RIGHT_FD_READ = 1n << 1n;
const RIGHT_FD_WRITE = 1n << 6n;
/** The rights of a file of the directory open for reading: every right but those that write. */
const FILE_RIGHTS =
  RIGHT_FD_READ |
  (1n << 2n) | // fd_seek
  (1n << 3n) | // fd_fdstat_set_flags
  (1n << 4n) | // fd_sync
  (1n << 5n) | // fd_tell
  (1n << 7n) | // fd_advise
  (1n << 21n); // fd_filestat_get
/** The rights that a file open for writing adds. */
const FILE_WRITE_RIGHTS =
  (1n << 0n) | // fd_datasync
  RIGHT_FD_WRITE |
  (1n << 8n) | // fd_allocate
  (1n << 22n) | // fd_filestat_set_size
  (1n << 23n); // fd_filestat_set_times
/** The rights of the directory: to look up, open and stat its files, and those of the files it opens. */
const DIRECTORY_RIGHTS =
  (1n << 13n) | // path_open
  (1n << 14n) | // fd_readdir
  (1n << 15n) | // path_readlink
  (1n << 18n) | // path_filestat_get
  (1n << 21n); // fd_filestat_get
/** The rights that a writable directory adds: to make, rename and remove its names. */
const DIRECTORY_WRITE_RIGHTS =
  (1n << 9n) | // path_create_directory
  (1n << 10n) | // path_create_file
  (1n << 16n) | // path_rename_source
  (1n << 17n) | // path_rename_target
  (1n << 19n) | // path_filestat_set_size
  (1n << 20n) | // path_filestat_set_times
  (1n << 25n) | // path_remove_directory
  (1n << 26n); // path_unlink_file
/** What a descriptor of a writable directory may do besides reading (`modes`). */
const MODE_WRITE = 1;
const MODE_APPEND = 2;
const MODE_NO_READ = 4;

/** The order of two names as bytes. */
function compareBytes(a: Uint8Array, b: Uint8Array): number {
  for (let i = 0; i < Math.min(a.length, b.length); i++) if (a[i] !== b[i]) return a[i]! - b[i]!;
  return a.length - b.length;
}

/** Now, in nanoseconds since the epoch: the time that a change gives a file. */
const nowNs = () => BigInt(Date.now()) * 1_000_000n;

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

/** A file of the directory: a user's (its bytes below the run's) or one that the run made. */
interface FileNode {
  readonly kind: 'file';
  /** Its number in `fileNodes`, which a descriptor's slot holds (plus 1). */
  readonly id: number;
  readonly ino: bigint;
  /** The user's bytes, if it was a user's file: those below `baseSize` that no chunk covers show. */
  readonly base: Blob | undefined;
  baseSize: number;
  /** What the run wrote, by chunk (index × CHUNK_BYTES); a chunk also holds the user's bytes of its range. */
  readonly chunks: Map<number, Uint8Array>;
  size: number;
  mtime: bigint;
  /** The path of the user's file that it was, if any. */
  readonly origin: string | undefined;
  /** The run made it, wrote it, cut or extended it. */
  changed: boolean;
  /** The folder that holds it under `name`; undefined once it is removed. */
  parent: DirNode | undefined;
  name: string;
}

/** A folder of the directory: the run directory (number 0), a folder of the user's paths, or one that the run made. */
interface DirNode {
  readonly kind: 'dir';
  /** Its number in `dirNodes`, which a descriptor's slot holds (minus, plus 1). */
  readonly id: number;
  readonly ino: bigint;
  readonly children: Map<string, FileNode | DirNode>;
  /** The folder that holds it under `name`; undefined for the run directory and once it is removed. */
  parent: DirNode | undefined;
  name: string;
  removed: boolean;
  mtime: bigint;
}

type DirectoryNode = FileNode | DirNode;

/** A path looked up: what it names, or the folder and name where it would be (a final `/`: `slash`), or an errno. */
type Found =
  | { readonly node: DirectoryNode; readonly name: string }
  | { readonly node: undefined; readonly parent: DirNode; readonly name: string; readonly slash: boolean }
  | { readonly errno: number };

/** An entry of a folder's listing (`fd_readdir`). */
interface ListEntry {
  readonly bytes: Uint8Array;
  readonly type: number;
  readonly ino: bigint;
}

/** A file that a tool run left in its directory (`RunDirectory.takeWritten`). */
export interface WrittenFile {
  /** Its path in the run's directory (`db.pin`, `sub/db.pin`). */
  readonly name: string;
  readonly size: number;
  /**
   * Its bytes from the start, in pieces of at most CHUNK_BYTES, each the caller's own: a chunk that the run wrote is
   * given up as it is given (so the worker's memory falls as the file goes), the user's bytes are read from the Blob,
   * and a hole is zeros.
   */
  pieces(): Generator<Uint8Array>;
}

/**
 * The run's working directory of one instance of the engine on one thread. The Engine worker opens the run's files
 * in it (`open`), which also ends the descriptors of an earlier run; a thread worker shares the table and takes the
 * same list for each thread that it runs (`take`). Its WASI functions (`imports`) stand in front of the shim's.
 */
export class RunDirectory {
  private files: readonly DirectoryFile[] = [];
  /** Whether a run may make, change and remove names (a tool's run; `open`). */
  private writable = false;
  /** The files, those of the list first (in its order), then those that the run made. */
  private fileNodes: FileNode[] = [];
  /** The folders: the run directory first, then those of the list's paths (in the order they come), then the run's. */
  private dirNodes: DirNode[] = [];
  private nextIno = 2n;
  private readonly slots: Int32Array;
  private readonly positions: BigInt64Array;
  /** Per slot of a writable directory, what its descriptor may do (MODE_*); 0 for a read-only directory's. */
  private readonly modes = new Uint8Array(OPEN_FILES);
  /** Per slot of a folder that is being listed, the entries as the listing began (`fd_readdir`). */
  private readonly listings = new Map<number, readonly ListEntry[]>();

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
    this.take([]);
  }

  /** The files that the user gave the directory. */
  get list(): readonly DirectoryFile[] {
    return this.files;
  }

  /** The bytes that the run's writes hold in this worker's memory (its chunks, also those of removed files still open). */
  get writtenBytes(): number {
    let total = 0;
    for (const node of this.fileNodes) for (const chunk of node.chunks.values()) total += chunk.length;
    return total;
  }

  /**
   * Gives the directory the files of a run, and ends the descriptors and drops the files that an earlier run left (the
   * Engine worker, before a run's first call into the instance, while no thread of the instance runs). Two files of
   * one path are an error. `writable` makes it a tool run's directory, which only an instance of its own table (the
   * serial one) may have.
   */
  open(files: readonly DirectoryFile[], options: { readonly writable?: boolean } = {}): void {
    const writable = options.writable === true;
    if (writable && typeof SharedArrayBuffer === 'function' && this.table instanceof SharedArrayBuffer) {
      throw new Error("a directory whose descriptors the instance's threads share is read-only: a tool runs on the serial instance");
    }
    this.take(files);
    this.writable = writable;
    this.slots.fill(0);
    this.positions.fill(0n);
    this.modes.fill(0);
    this.listings.clear();
  }

  /** Gives the directory the files of a run that the Engine worker opened (`open`): a thread worker's. */
  take(files: readonly DirectoryFile[]): void {
    const root: DirNode = { kind: 'dir', id: 0, ino: 1n, children: new Map(), parent: undefined, name: '', removed: false, mtime: 0n };
    const fileNodes: FileNode[] = [];
    // The folders are numbered in the order that the list's paths give them (the same on every thread), and their
    // inode numbers come after the files'.
    const dirNodes: DirNode[] = [root];
    files.forEach((file, id) => {
      const names = pathNames(file.name);
      if (names === undefined || names.length === 0) throw new Error(`the file name ${JSON.stringify(file.name)} is not a path in the run's directory`);
      const path = names.join('/');
      let folder = root;
      for (let i = 0; i < names.length - 1; i++) {
        const child = folder.children.get(names[i]!);
        if (child?.kind === 'file') throw new Error(`${names.slice(0, i + 1).join('/')} is both a file and a folder`);
        if (child !== undefined) {
          folder = child;
          continue;
        }
        const number = dirNodes.length;
        const made: DirNode = { kind: 'dir', id: number, ino: BigInt(files.length + 1 + number), children: new Map(), parent: folder, name: names[i]!, removed: false, mtime: 0n };
        folder.children.set(made.name, made);
        dirNodes.push(made);
        folder = made;
      }
      const name = names.at(-1)!;
      const there = folder.children.get(name);
      if (there?.kind === 'file') throw new Error(`two files are named ${path}`);
      if (there !== undefined) throw new Error(`${path} is both a file and a folder`);
      const size = file.file.size;
      const mtime = BigInt(Math.max(0, Math.trunc((file.file as File).lastModified ?? 0))) * 1_000_000n;
      const node: FileNode = { kind: 'file', id, ino: BigInt(id + 2), base: file.file, baseSize: size, chunks: new Map(), size, mtime, origin: path, changed: false, parent: folder, name };
      folder.children.set(name, node);
      fileNodes.push(node);
    });
    this.files = files;
    this.fileNodes = fileNodes;
    this.dirNodes = dirNodes;
    this.nextIno = BigInt(files.length + 1 + dirNodes.length);
  }

  /**
   * The files that the run made or changed (written, cut, extended, or renamed from the user's path) and that are still
   * in the directory, by their path in byte order; the folders hold no output of their own. A file's chunks are given up
   * as its pieces are taken, and the directory holds no file of the run once `open` gives it the next run's files.
   */
  takeWritten(): WrittenFile[] {
    const found: Array<{ readonly node: FileNode; readonly name: string; readonly bytes: Uint8Array }> = [];
    const encoder = new TextEncoder();
    const walk = (folder: DirNode, prefix: string) => {
      for (const [name, child] of folder.children) {
        const path = `${prefix}${name}`;
        if (child.kind === 'dir') walk(child, `${path}/`);
        else if (child.changed || child.origin !== path) found.push({ node: child, name: path, bytes: encoder.encode(path) });
      }
    };
    walk(this.dirNodes[0]!, '');
    found.sort((a, b) => compareBytes(a.bytes, b.bytes));
    return found.map(({ node, name }) => ({ name, size: node.size, pieces: () => this.pieces(node) }));
  }

  /** The bytes of `node` from the start in pieces of at most CHUNK_BYTES (`WrittenFile.pieces`). */
  private *pieces(node: FileNode): Generator<Uint8Array> {
    for (let start = 0; start < node.size; start += CHUNK_BYTES) {
      const index = start / CHUNK_BYTES;
      const length = Math.min(CHUNK_BYTES, node.size - start);
      const chunk = node.chunks.get(index);
      if (chunk !== undefined) {
        node.chunks.delete(index);
        if (chunk.length === length) yield chunk;
        else if (chunk.length > length) yield chunk.slice(0, length);
        else {
          // A hole after the bytes that the chunk holds.
          const whole = new Uint8Array(length);
          whole.set(chunk);
          yield whole;
        }
        continue;
      }
      const piece = new Uint8Array(length);
      const fromBase = Math.min(length, Math.max(0, node.baseSize - start));
      if (fromBase > 0) piece.set(this.readBase(node, start, start + fromBase));
      yield piece;
    }
  }

  /** The user's bytes `start` to `end` of `node`'s Blob. */
  private readBase(node: FileNode, start: number, end: number): Uint8Array {
    if (this.read === undefined) throw new Error('this thread cannot read a file synchronously');
    const data = this.read(node.base!, start, end);
    if (data.length !== end - start) throw new Error(`a read of ${end - start} bytes of a file gave ${data.length}`);
    return data;
  }

  /** `length` bytes of `node` from `start` (within its size): its chunks, the user's bytes below them, zeros in a hole. */
  private readNode(node: FileNode, start: number, length: number): Uint8Array {
    const end = start + length;
    if (node.chunks.size === 0 && end <= node.baseSize) return this.readBase(node, start, end);
    const out = new Uint8Array(length);
    let at = start;
    while (at < end) {
      const index = Math.floor(at / CHUNK_BYTES);
      const chunkStart = index * CHUNK_BYTES;
      const chunk = node.chunks.get(index);
      if (chunk !== undefined) {
        const pieceEnd = Math.min(end, chunkStart + CHUNK_BYTES);
        const held = Math.min(pieceEnd, chunkStart + chunk.length);
        if (held > at) out.set(chunk.subarray(at - chunkStart, held - chunkStart), at - start);
        at = pieceEnd;
        continue;
      }
      // The user's bytes up to the next chunk that the run wrote, in one read.
      let runEnd = Math.min(end, chunkStart + CHUNK_BYTES);
      while (runEnd < end && !node.chunks.has(runEnd / CHUNK_BYTES)) runEnd = Math.min(end, runEnd + CHUNK_BYTES);
      const baseEnd = Math.min(runEnd, node.baseSize);
      if (baseEnd > at) out.set(this.readBase(node, at, baseEnd), at - start);
      at = runEnd;
    }
    return out;
  }

  /** The chunk `index` of `node` with room for its first `needed` bytes: made from the user's bytes of its range, or grown. */
  private chunkFor(node: FileNode, index: number, needed: number): Uint8Array {
    const room = (bytes: number) => Math.min(CHUNK_BYTES, Math.max(MIN_CHUNK_BYTES, 2 ** Math.ceil(Math.log2(Math.max(1, bytes)))));
    const held = node.chunks.get(index);
    if (held === undefined) {
      const chunkStart = index * CHUNK_BYTES;
      const fromBase = Math.max(0, Math.min(node.baseSize, chunkStart + CHUNK_BYTES) - chunkStart);
      const chunk = new Uint8Array(room(Math.max(needed, fromBase)));
      if (fromBase > 0) chunk.set(this.readBase(node, chunkStart, chunkStart + fromBase));
      node.chunks.set(index, chunk);
      return chunk;
    }
    if (held.length >= needed) return held;
    const grown = new Uint8Array(room(needed));
    grown.set(held);
    node.chunks.set(index, grown);
    return grown;
  }

  /** Writes `data` into `node` at `offset`, past its end too (a hole between reads as zeros). */
  private writeNode(node: FileNode, offset: number, data: Uint8Array): void {
    const end = offset + data.length;
    for (let at = offset; at < end; ) {
      const index = Math.floor(at / CHUNK_BYTES);
      const chunkStart = index * CHUNK_BYTES;
      const pieceEnd = Math.min(end, chunkStart + CHUNK_BYTES);
      this.chunkFor(node, index, pieceEnd - chunkStart).set(data.subarray(at - offset, pieceEnd - offset), at - chunkStart);
      at = pieceEnd;
    }
    if (end > node.size) node.size = end;
    node.changed = true;
    node.mtime = nowNs();
  }

  /** Makes `node` `size` bytes long: cut (the user's bytes and the chunks past it go), or extended with zeros. */
  private resize(node: FileNode, size: number): void {
    if (size < node.size) {
      node.baseSize = Math.min(node.baseSize, size);
      for (const [index, chunk] of node.chunks) {
        const chunkStart = index * CHUNK_BYTES;
        if (chunkStart >= size) node.chunks.delete(index);
        else if (chunkStart + chunk.length > size) chunk.fill(0, size - chunkStart);
      }
    }
    node.size = size;
    node.changed = true;
    node.mtime = nowNs();
  }

  /** The folder `folder` and every folder in it are no longer in the directory. */
  private static removeTree(folder: DirNode): void {
    folder.removed = true;
    folder.parent = undefined;
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
     * What a path of the directory names, looked up from `from`. Each name is looked up in the folder that the names
     * before it give, as the OS looks a path up: a name (`..` too) after a file is ENOTDIR, after a name that is not
     * there ENOENT (a final `/` after it asks for a folder: `slash`); `..` past the directory is ENOTCAPABLE, as WASI
     * gives it; a folder that the run removed holds nothing.
     */
    const lookup = (path: string, from: DirNode): Found => {
      let node: DirectoryNode | undefined = from;
      let parent: DirNode = from;
      let name = '.';
      const parts = path.split('/');
      for (let i = 0; i < parts.length; i++) {
        const part = parts[i]!;
        if (node === undefined) {
          if (parts.slice(i).every((rest) => rest === '')) return { node: undefined, parent, name, slash: true };
          return { errno: ERRNO_NOENT };
        }
        if (node.kind === 'file') return { errno: ERRNO_NOTDIR };
        if (part === '') continue;
        if (part === '.') {
          name = part;
          continue;
        }
        if (node.removed) return { errno: ERRNO_NOENT };
        if (part === '..') {
          if (node.parent === undefined) return node.id === 0 ? { errno: ERRNO_NOTCAPABLE } : { errno: ERRNO_NOENT };
          node = node.parent;
          name = part;
          continue;
        }
        parent = node;
        name = part;
        node = node.children.get(part);
      }
      return node === undefined ? { node: undefined, parent, name, slash: false } : { node, name };
    };
    const slotOf = (fd: number) => fd - FILE_FD_BASE;
    const isFile = (fd: number) => slotOf(fd) >= 0 && slotOf(fd) < OPEN_FILES;
    /** The file of an open descriptor, or undefined (a free slot, and a folder's). */
    const openFile = (fd: number): FileNode | undefined => {
      const id = Atomics.load(this.slots, slotOf(fd)) - 1;
      return id < 0 ? undefined : this.fileNodes[id];
    };
    /** The folder that an open descriptor is, or undefined. */
    const openFolder = (fd: number): DirNode | undefined => {
      const state = Atomics.load(this.slots, slotOf(fd));
      return state < 0 ? this.dirNodes[-state - 1] : undefined;
    };
    /** The folder that a path is looked up from: the preopen's, or an open folder's; undefined for any other descriptor. */
    const baseOf = (fd: number): DirNode | undefined => (fd === PREOPEN_FD ? this.dirNodes[0] : isFile(fd) ? openFolder(fd) : undefined);
    /** The errno of a path call on a descriptor of the directory that is no folder (a free slot's, a file's). */
    const notFolder = (fd: number) => (Atomics.load(this.slots, slotOf(fd)) === 0 ? ERRNO_BADF : ERRNO_NOTDIR);
    const writeFilestat = (pointer: number, node: DirectoryNode) => {
      const out = view(pointer, 64);
      out.setBigUint64(0, 0n, true); // dev
      out.setBigUint64(8, node.ino, true); // ino: the run directory 1, the user's files 2 on, then the folders and the run's files
      out.setBigUint64(16, 0n, true);
      out.setUint8(16, node.kind === 'dir' ? FILETYPE_DIRECTORY : FILETYPE_REGULAR_FILE);
      const linked = node.kind === 'dir' ? !node.removed : node.parent !== undefined;
      out.setBigUint64(24, linked ? 1n : 0n, true); // nlink
      out.setBigUint64(32, BigInt(node.kind === 'file' ? node.size : 0), true);
      for (const at of [40, 48, 56]) out.setBigUint64(at, node.mtime, true);
    };
    /** Reads up to the length of the vectors of `iovs` from `node` at `offset` into them; returns the count. */
    const readInto = (node: FileNode, iovs: number, count: number, offset: bigint): number => {
      const size = BigInt(node.size);
      if (offset >= size) return 0;
      const vectors = view(iovs, 8 * count);
      let wanted = 0;
      for (let i = 0; i < count; i++) wanted += vectors.getUint32(8 * i + 4, true);
      const length = Math.min(wanted, Number(size - offset), READ_LIMIT);
      if (length === 0) return 0;
      const data = this.readNode(node, Number(offset), length);
      let done = 0;
      for (let i = 0; i < count && done < data.length; i++) {
        const pointer = vectors.getUint32(8 * i, true);
        const part = Math.min(vectors.getUint32(8 * i + 4, true), data.length - done);
        bytes(pointer, part).set(data.subarray(done, done + part));
        done += part;
      }
      return done;
    };
    /** Writes the vectors of `iovs` into `node` from `offset`; returns the count, or an errno past the largest file. */
    const writeFrom = (node: FileNode, iovs: number, count: number, offset: number): { readonly count: number } | { readonly errno: number } => {
      const vectors = view(iovs, 8 * count);
      let total = 0;
      for (let i = 0; i < count; i++) total += vectors.getUint32(8 * i + 4, true);
      if (offset + total > Number.MAX_SAFE_INTEGER) return { errno: ERRNO_FBIG };
      let at = offset;
      for (let i = 0; i < count; i++) {
        const length = vectors.getUint32(8 * i + 4, true);
        if (length === 0) continue;
        this.writeNode(node, at, bytes(vectors.getUint32(8 * i, true), length));
        at += length;
      }
      return { count: total };
    };
    /** Whether a descriptor's slot is open for writing (a file of a writable directory opened so). */
    const writes = (fd: number) => (this.modes[slotOf(fd)]! & MODE_WRITE) !== 0;
    /** Ends the bytes of a removed file once no descriptor has it open. */
    const dropIfGone = (node: FileNode) => {
      if (node.parent !== undefined) return;
      for (let slot = 0; slot < OPEN_FILES; slot++) if (Atomics.load(this.slots, slot) === node.id + 1) return;
      node.chunks.clear();
    };
    /** Takes `node` out of its folder (the folder changes; a file stays readable through its open descriptors). */
    const unlink = (node: DirectoryNode) => {
      const parent = node.parent;
      if (parent === undefined) return;
      parent.children.delete(node.name);
      parent.mtime = nowNs();
      if (node.kind === 'dir') RunDirectory.removeTree(node);
      else {
        node.parent = undefined;
        dropIfGone(node);
      }
    };
    /**
     * Calls `body` for a descriptor of the directory's files; errors of a read are EIO, and a write that finds no
     * memory is ENOSPC (as a full disk). A descriptor of a folder gets `onFolder`: an errno (what Linux answers for that
     * call on a directory descriptor) or a function of its own.
     */
    const onFile = (
      name: string,
      body: (fd: number, file: FileNode, ...args: Array<number | bigint>) => number,
      onFolder: number | ((fd: number, folder: DirNode, ...args: Array<number | bigint>) => number) = ERRNO_BADF,
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
        } catch (error) {
          return error instanceof RangeError ? ERRNO_NOSPC : ERRNO_IO;
        }
      };
    };
    /** Calls `body` for a path of the directory looked up from `dirfd`; other descriptors go to the shim. */
    const onPath = (name: string, body: (from: DirNode, ...args: Array<number | bigint>) => number) => {
      own[name] = (dirfd, ...args) => {
        const descriptor = Number(dirfd);
        const from = baseOf(descriptor);
        if (from !== undefined) return body(from, ...args);
        if (isFile(descriptor)) return notFolder(descriptor);
        return fallback[name]!(dirfd, ...args);
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
    onPath('path_open', (from, _dirflags, path, pathLength, oflags, rightsBase, _rightsInheriting, fdflags, opened) => {
      const found = lookup(pathAt(Number(path), Number(pathLength)), from);
      if ('errno' in found) return found.errno;
      const flags = Number(oflags);
      const appends = (Number(fdflags) & FDFLAGS_APPEND) !== 0;
      const forWriting = (BigInt(rightsBase!) & RIGHT_FD_WRITE) !== 0n || appends;
      let node: DirectoryNode;
      if (found.node === undefined) {
        if (!(flags & OFLAGS_CREAT)) return ERRNO_NOENT;
        if (!this.writable) return ERRNO_ROFS;
        // As Linux: `open("x/", O_CREAT)` is EISDIR, and O_CREAT with O_DIRECTORY is EINVAL.
        if (found.slash) return ERRNO_ISDIR;
        if (flags & OFLAGS_DIRECTORY) return ERRNO_INVAL;
        const made: FileNode = {
          kind: 'file',
          id: this.fileNodes.length,
          ino: this.nextIno++,
          base: undefined,
          baseSize: 0,
          chunks: new Map(),
          size: 0,
          mtime: nowNs(),
          origin: undefined,
          changed: true,
          parent: found.parent,
          name: found.name,
        };
        this.fileNodes.push(made);
        found.parent.children.set(made.name, made);
        found.parent.mtime = made.mtime;
        node = made;
      } else {
        node = found.node;
        if (flags & OFLAGS_CREAT && flags & OFLAGS_EXCL) return ERRNO_EXIST;
        if (node.kind === 'dir') {
          // As Linux: a folder opens for reading (its reads fail with EISDIR); for writing, with O_CREAT or O_TRUNC it is EISDIR.
          if (forWriting || flags & (OFLAGS_CREAT | OFLAGS_TRUNC)) return ERRNO_ISDIR;
        } else {
          if (flags & OFLAGS_DIRECTORY) return ERRNO_NOTDIR;
          if (!this.writable && (forWriting || flags & OFLAGS_TRUNC)) return ERRNO_ROFS;
          if (flags & OFLAGS_TRUNC && node.size > 0) this.resize(node, 0);
        }
      }
      const state = node.kind === 'dir' ? -(node.id + 1) : node.id + 1;
      for (let slot = 0; slot < OPEN_FILES; slot++) {
        if (Atomics.compareExchange(this.slots, slot, 0, state) !== 0) continue;
        Atomics.store(this.positions, slot, 0n);
        this.modes[slot] = !this.writable || node.kind === 'dir' ? 0 : (forWriting ? MODE_WRITE : 0) | (appends ? MODE_APPEND : 0) | ((BigInt(rightsBase!) & RIGHT_FD_READ) === 0n ? MODE_NO_READ : 0);
        this.listings.delete(slot);
        view(Number(opened), 4).setUint32(0, FILE_FD_BASE + slot, true);
        return ERRNO_SUCCESS;
      }
      return ERRNO_MFILE;
    });
    onPath('path_filestat_get', (from, _flags, path, pathLength, pointer) => {
      const found = lookup(pathAt(Number(path), Number(pathLength)), from);
      if ('errno' in found) return found.errno;
      if (found.node === undefined) return ERRNO_NOENT;
      writeFilestat(Number(pointer), found.node);
      return ERRNO_SUCCESS;
    });
    onPath('path_readlink', (from, path, pathLength) => {
      const found = lookup(pathAt(Number(path), Number(pathLength)), from);
      if ('errno' in found) return found.errno;
      // No file of the directory is a symbolic link.
      return found.node === undefined ? ERRNO_NOENT : ERRNO_INVAL;
    });
    onPath('path_create_directory', (from, path, pathLength) => {
      if (!this.writable) return ERRNO_ROFS;
      const found = lookup(pathAt(Number(path), Number(pathLength)), from);
      if ('errno' in found) return found.errno;
      if (found.node !== undefined) return ERRNO_EXIST;
      const made: DirNode = { kind: 'dir', id: this.dirNodes.length, ino: this.nextIno++, children: new Map(), parent: found.parent, name: found.name, removed: false, mtime: nowNs() };
      this.dirNodes.push(made);
      found.parent.children.set(made.name, made);
      found.parent.mtime = made.mtime;
      return ERRNO_SUCCESS;
    });
    onPath('path_unlink_file', (from, path, pathLength) => {
      if (!this.writable) return ERRNO_ROFS;
      const found = lookup(pathAt(Number(path), Number(pathLength)), from);
      if ('errno' in found) return found.errno;
      if (found.node === undefined) return ERRNO_NOENT;
      // As Linux: unlink of a folder is EISDIR; `x/` of a file was ENOTDIR above.
      if (found.node.kind === 'dir') return ERRNO_ISDIR;
      unlink(found.node);
      return ERRNO_SUCCESS;
    });
    onPath('path_remove_directory', (from, path, pathLength) => {
      if (!this.writable) return ERRNO_ROFS;
      const found = lookup(pathAt(Number(path), Number(pathLength)), from);
      if ('errno' in found) return found.errno;
      if (found.node === undefined) return ERRNO_NOENT;
      if (found.node.kind === 'file') return ERRNO_NOTDIR;
      // As Linux: `rmdir(".")` is EINVAL, `rmdir("..")` ENOTEMPTY, the run directory itself EBUSY.
      if (found.name === '.') return ERRNO_INVAL;
      if (found.name === '..') return ERRNO_NOTEMPTY;
      if (found.node.parent === undefined) return ERRNO_BUSY;
      if (found.node.children.size > 0) return ERRNO_NOTEMPTY;
      unlink(found.node);
      return ERRNO_SUCCESS;
    });
    own.path_rename = (oldfd, oldPath, oldLength, newfd, newPath, newLength) => {
      const oldFrom = baseOf(Number(oldfd));
      const newFrom = baseOf(Number(newfd));
      if (oldFrom === undefined && newFrom === undefined) {
        if (isFile(Number(oldfd))) return notFolder(Number(oldfd));
        if (isFile(Number(newfd))) return notFolder(Number(newfd));
        return fallback.path_rename!(oldfd, oldPath, oldLength, newfd, newPath, newLength);
      }
      // One side outside the directory: another file system (EXDEV), or a descriptor that is no folder.
      if (oldFrom === undefined) return isFile(Number(oldfd)) ? notFolder(Number(oldfd)) : ERRNO_XDEV;
      if (newFrom === undefined) return isFile(Number(newfd)) ? notFolder(Number(newfd)) : ERRNO_XDEV;
      if (!this.writable) return ERRNO_ROFS;
      const source = lookup(pathAt(Number(oldPath), Number(oldLength)), oldFrom);
      if ('errno' in source) return source.errno;
      if (source.node === undefined) return ERRNO_NOENT;
      const target = lookup(pathAt(Number(newPath), Number(newLength)), newFrom);
      if ('errno' in target) return target.errno;
      const node = source.node;
      // As Linux: the names `.` and `..` and the run directory itself do not move (EBUSY).
      if (source.name === '.' || source.name === '..' || node.parent === undefined) return ERRNO_BUSY;
      if (target.node === undefined) {
        // `new/` names a folder: a file does not go there (ENOTDIR).
        if (target.slash && node.kind === 'file') return ERRNO_NOTDIR;
      } else {
        if (target.node === node) return ERRNO_SUCCESS;
        if (target.name === '.' || target.name === '..' || target.node.parent === undefined) return ERRNO_BUSY;
        if (node.kind === 'file' && target.node.kind === 'dir') return ERRNO_ISDIR;
        if (node.kind === 'dir' && target.node.kind === 'file') return ERRNO_NOTDIR;
        if (target.node.kind === 'dir' && target.node.children.size > 0) return ERRNO_NOTEMPTY;
      }
      const parent = target.node === undefined ? target.parent : target.node.parent!;
      // A folder does not go into itself or into a folder of its own (EINVAL).
      if (node.kind === 'dir') for (let up: DirNode | undefined = parent; up !== undefined; up = up.parent) if (up === node) return ERRNO_INVAL;
      if (target.node !== undefined) unlink(target.node);
      node.parent!.children.delete(node.name);
      node.parent!.mtime = nowNs();
      node.parent = parent;
      node.name = target.name;
      parent.children.set(node.name, node);
      parent.mtime = nowNs();
      return ERRNO_SUCCESS;
    };
    own.fd_filestat_get = (fd, pointer) => {
      if (Number(fd) === PREOPEN_FD) {
        writeFilestat(Number(pointer), this.dirNodes[0]!);
        return ERRNO_SUCCESS;
      }
      if (!isFile(Number(fd))) return fallback.fd_filestat_get!(fd, pointer);
      const node = openFolder(Number(fd)) ?? openFile(Number(fd));
      if (node === undefined) return ERRNO_BADF;
      writeFilestat(Number(pointer), node);
      return ERRNO_SUCCESS;
    };
    own.fd_fdstat_get = (fd, pointer) => {
      const folder = Number(fd) === PREOPEN_FD || (isFile(Number(fd)) && openFolder(Number(fd)) !== undefined);
      if (!folder && !isFile(Number(fd))) return fallback.fd_fdstat_get!(fd, pointer);
      if (!folder && openFile(Number(fd)) === undefined) return ERRNO_BADF;
      const mode = folder ? 0 : this.modes[slotOf(Number(fd))]!;
      const fileRights = (mode & MODE_NO_READ ? FILE_RIGHTS & ~RIGHT_FD_READ : FILE_RIGHTS) | (mode & MODE_WRITE ? FILE_WRITE_RIGHTS : 0n);
      const folderRights = DIRECTORY_RIGHTS | (this.writable ? DIRECTORY_WRITE_RIGHTS : 0n);
      const out = view(Number(pointer), 24);
      out.setBigUint64(0, 0n, true);
      out.setUint8(0, folder ? FILETYPE_DIRECTORY : FILETYPE_REGULAR_FILE);
      out.setUint16(2, mode & MODE_APPEND ? FDFLAGS_APPEND : 0, true);
      out.setBigUint64(8, folder ? folderRights : fileRights, true);
      out.setBigUint64(16, folder ? FILE_RIGHTS | (this.writable ? FILE_WRITE_RIGHTS : 0n) | folderRights : 0n, true);
      return ERRNO_SUCCESS;
    };
    // A read of a folder's descriptor is EISDIR, as on Linux, where the open succeeded; of a file open for writing only, EBADF.
    onFile(
      'fd_read',
      (fd, file, iovs, count, nread) => {
        const slot = slotOf(fd);
        if (this.modes[slot]! & MODE_NO_READ) return ERRNO_BADF;
        const got = readInto(file, Number(iovs), Number(count), Atomics.load(this.positions, slot));
        Atomics.add(this.positions, slot, BigInt(got));
        view(Number(nread), 4).setUint32(0, got, true);
        return ERRNO_SUCCESS;
      },
      ERRNO_ISDIR,
    );
    onFile(
      'fd_pread',
      (fd, file, iovs, count, offset, nread) => {
        if (this.modes[slotOf(fd)]! & MODE_NO_READ) return ERRNO_BADF;
        const got = readInto(file, Number(iovs), Number(count), BigInt.asUintN(64, BigInt(offset!)));
        view(Number(nread), 4).setUint32(0, got, true);
        return ERRNO_SUCCESS;
      },
      ERRNO_ISDIR,
    );
    // A write goes at the descriptor's position, or at the end of a file opened to append (as Linux's O_APPEND, also for
    // `pwrite`); a descriptor that is not open for writing (a folder's, a file's opened to read) is EBADF.
    onFile('fd_write', (fd, file, iovs, count, written) => {
      const slot = slotOf(fd);
      if (!writes(fd)) return ERRNO_BADF;
      const at = this.modes[slot]! & MODE_APPEND ? file.size : Number(Atomics.load(this.positions, slot));
      const done = writeFrom(file, Number(iovs), Number(count), at);
      if ('errno' in done) return done.errno;
      Atomics.store(this.positions, slot, BigInt(at + done.count));
      view(Number(written), 4).setUint32(0, done.count, true);
      return ERRNO_SUCCESS;
    });
    onFile('fd_pwrite', (fd, file, iovs, count, offset, written) => {
      if (!writes(fd)) return ERRNO_BADF;
      const asked = BigInt.asIntN(64, BigInt(offset!));
      if (asked < 0n) return ERRNO_INVAL;
      if (asked > BigInt(Number.MAX_SAFE_INTEGER)) return ERRNO_FBIG;
      const at = this.modes[slotOf(fd)]! & MODE_APPEND ? file.size : Number(asked);
      const done = writeFrom(file, Number(iovs), Number(count), at);
      if ('errno' in done) return done.errno;
      view(Number(written), 4).setUint32(0, done.count, true);
      return ERRNO_SUCCESS;
    });
    // As Linux's ftruncate: a file not open for writing (or a folder) is EINVAL.
    onFile(
      'fd_filestat_set_size',
      (fd, file, size) => {
        if (!writes(fd)) return ERRNO_INVAL;
        const length = BigInt.asUintN(64, BigInt(size!));
        if (length > BigInt(Number.MAX_SAFE_INTEGER)) return ERRNO_FBIG;
        this.resize(file, Number(length));
        return ERRNO_SUCCESS;
      },
      ERRNO_INVAL,
    );
    // As Linux's posix_fallocate: the file is at least `offset + length` long; a descriptor not open for writing is EBADF.
    onFile('fd_allocate', (fd, file, offset, length) => {
      if (!writes(fd)) return ERRNO_BADF;
      const end = BigInt.asUintN(64, BigInt(offset!)) + BigInt.asUintN(64, BigInt(length!));
      if (end > BigInt(Number.MAX_SAFE_INTEGER)) return ERRNO_FBIG;
      if (Number(end) > file.size) this.resize(file, Number(end));
      return ERRNO_SUCCESS;
    });
    const setTimes = (node: DirectoryNode, mtime: number | bigint | undefined, flags: number | bigint | undefined) => {
      if (!this.writable) return ERRNO_ROFS;
      const which = Number(flags);
      if (which & FSTFLAGS_MTIM_NOW) node.mtime = nowNs();
      else if (which & FSTFLAGS_MTIM) node.mtime = BigInt.asUintN(64, BigInt(mtime!));
      return ERRNO_SUCCESS;
    };
    onFile(
      'fd_filestat_set_times',
      (_fd, file, _atim, mtim, flags) => setTimes(file, mtim, flags),
      (_fd, folder, _atim, mtim, flags) => setTimes(folder, mtim, flags),
    );
    onPath('path_filestat_set_times', (from, _lookupFlags, path, pathLength, _atim, mtim, flags) => {
      const found = lookup(pathAt(Number(path), Number(pathLength)), from);
      if ('errno' in found) return found.errno;
      if (found.node === undefined) return ERRNO_NOENT;
      return setTimes(found.node, mtim, flags);
    });
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
    onFile('fd_seek', (fd, file, offset, whence, result) => seek(fd, file.size, offset, whence, result), (fd, _folder, offset, whence, result) =>
      seek(fd, 0, offset, whence, result),
    );
    const tell = (fd: number, result: number | bigint | undefined) => {
      view(Number(result), 8).setBigUint64(0, Atomics.load(this.positions, slotOf(fd)), true);
      return ERRNO_SUCCESS;
    };
    onFile('fd_tell', (fd, _file, result) => tell(fd, result), (fd, _folder, result) => tell(fd, result));
    for (const name of ['fd_advise', 'fd_sync', 'fd_datasync']) onFile(name, () => ERRNO_SUCCESS, ERRNO_SUCCESS);
    // The flags of a descriptor: only O_APPEND, on a file open for writing, as Linux's fcntl(F_SETFL) sets it.
    const setFlags = (fd: number, flags: number | bigint | undefined, file: boolean) => {
      const wanted = Number(flags);
      if (wanted === 0 || (file && wanted === FDFLAGS_APPEND && writes(fd))) {
        const slot = slotOf(fd);
        this.modes[slot] = (this.modes[slot]! & ~MODE_APPEND) | (wanted === FDFLAGS_APPEND ? MODE_APPEND : 0);
        return ERRNO_SUCCESS;
      }
      return ERRNO_NOTCAPABLE;
    };
    onFile('fd_fdstat_set_flags', (fd, _file, flags) => setFlags(fd, flags, true), (fd, _folder, flags) => setFlags(fd, flags, false));
    onFile('fd_readdir', () => ERRNO_NOTDIR, (fd, folder, buffer, length, cookie, used) => {
      // The folder's entries, as the OS lists them: `.` and `..`, then the names in the folder, in byte order. A listing
      // that goes on (a cookie after the first) reads the entries that the folder had when it began (POSIX leaves a
      // name made or removed since then open), so a run that removes what it lists sees each entry once.
      const slot = slotOf(fd);
      const next0 = Number(BigInt.asUintN(32, BigInt(cookie!)));
      let entries = next0 === 0 ? undefined : this.listings.get(slot);
      if (entries === undefined) {
        const encoder = new TextEncoder();
        const children = [...folder.children.values()].map((child) => ({
          bytes: encoder.encode(child.name),
          type: child.kind === 'dir' ? FILETYPE_DIRECTORY : FILETYPE_REGULAR_FILE,
          ino: child.ino,
        }));
        children.sort((a, b) => compareBytes(a.bytes, b.bytes));
        entries = [
          { bytes: encoder.encode('.'), type: FILETYPE_DIRECTORY, ino: folder.ino },
          { bytes: encoder.encode('..'), type: FILETYPE_DIRECTORY, ino: (folder.parent ?? folder).ino },
          ...children,
        ];
        this.listings.set(slot, entries);
      }
      const room = Number(length);
      const target = bytes(Number(buffer), room);
      let written = 0;
      for (let next = next0; next < entries.length && written < room; next++) {
        const entry = entries[next]!;
        const record = new Uint8Array(24 + entry.bytes.length);
        const header = new DataView(record.buffer);
        header.setBigUint64(0, BigInt(next + 1), true); // d_next: the cookie of the entry after this
        header.setBigUint64(8, entry.ino, true); // d_ino
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
      this.modes[slot] = 0;
      this.listings.delete(slot);
      const node = state > 0 ? this.fileNodes[state - 1] : undefined;
      if (node !== undefined) dropIfGone(node);
      return ERRNO_SUCCESS;
    };
    const merged: WasiImports = { ...fallback };
    // Only the functions that the module imports from the shim are replaced; the others stay out.
    for (const [name, call] of Object.entries(own)) if (name in fallback) merged[name] = call;
    return merged as WebAssembly.ModuleImports;
  }
}
