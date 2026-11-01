"use strict";
// The table of the files that a LOSAT WASI command opens, kept by the Node host itself in front of
// uvwasi: the threaded host (wasi_thread_host.js, every thread of a command) and the serial runner
// (run_losat_wasi.js).
const fs = require("fs");
// node:tty is loaded at its first use (an open that the table takes).
const tty = { isatty: descriptor => require("tty").isatty(descriptor) };

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbatlas.cpp:114-127
// ```c++
// CMemoryFile* CSeqDBAtlas::GetMemoryFile(const string& fileName)
// {
//     std::lock_guard<std::mutex> guard(m_FileMemMapMutex);
//     auto it = m_FileMemMap.find(fileName);
//     ...
//     CAtlasMappedFile* file(new CAtlasMappedFile(fileName));
//     m_FileMemMap[fileName].reset(file);
// ```
// NCBI opens a database file once for the process, and every thread reads it there. LOSAT does
// the same with one descriptor (LOSAT/src/seqdb/file.rs: a seek and a read under a lock), so a
// descriptor must name the same open file on every thread, as in the native process's table.
// Node's WASI keeps its table in the instance (uvwasi), and each thread worker has its own, so a
// descriptor that one thread opened was EBADF on the others (-db with 2 or more threads,
// 2026-10-26). The host therefore opens itself (Node's fs, whose descriptors belong to the
// process) each file that path_open opens from the preopened "/", enters it in a table that every
// thread shares (the OS descriptor, the file position, the rights and the kind) and gives the
// guest SHARED_FD_BASE + its slot; each thread then reads, writes, seeks, stats and closes such a
// descriptor itself, through the table.
// (One I/O worker that ran the file functions of every thread was tried first: each call waited
// until that worker was scheduled on the 16 busy CPUs, and BLASTN at 16 threads took 22.6 s
// instead of 9.3 s on ws, 2026-10-27.)
//
// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:225-228
// ```c++
//         CBlastFormat formatter(opt, *db_adapter,
//                                fmt_args->GetFormattedOutputChoice(),
//                                query_opts->GetParseDeflines(),
//                                m_CmdLineArgs->GetOutputStream(),
// ```
// The output (-out) is one stream of the process too. LOSAT creates the file at the report's
// first write, and a write from a thread other than the one that opened -out must work: a file
// opened by one thread and written by another gave EBADF, so opens for writing enter the table as
// well as opens for reading.
//
// NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:727-729
// ```c++
//             fstrm->open(AsString().c_str(),IOS_BASE::in | mode);
//             if ( !fstrm->is_open() ) {
//                 delete fstrm;
// ```
// NCBI opens a file argument by its name through the OS, so -query /dev/stdin and -out /dev/stdout
// name the process's standard input and output whatever they are. uvwasi resolves symbolic links
// itself and cannot follow /dev/stdout -> /proc/self/fd/1 -> "pipe:[N]": with a pipe there it gave
// ENOENT ("File is not accessible", at 1 thread too and in the serial runner, which uses this
// table too). The host opens these by name as the OS does.
//
// The table takes an open from the preopened "/" of a regular file (a new one for writing, with
// create), a character device or a FIFO, with no open flag but create, exclusive and truncate (for
// writing) and no descriptor flag but append (for writing); stat before the open decides whether
// the table takes it (a FIFO's open waits for its other end, as native), fstat after the open
// gives its kind. A regular file is read at the table's position, and written at the OS's
// position with one writev(), which moves as the native offset moves (on macOS /dev/stdout is a
// duplicate of descriptor 1 and shares its offset), until the guest seeks or reads it; from then on
// it is written at the table's position, the bytes of one write reserved before they are written
// as one write() moves the native offset (Node has no lseek). (A seek by one thread between another
// thread's write at the OS's position and its count could leave the position past that write;
// LOSAT does not seek a file that another thread writes.) An append leaves the position at the end
// of the file. A device or FIFO is read and written at the OS's position, and a full one
// (EAGAIN) is written again after 1 ms as blockingStdioWrite writes descriptors 1 and 2. A
// character device that is not a terminal seeks and tells 0, preads and pwrites, as Linux's
// /dev/null and /dev/zero do (LOSAT tells and seeks its FASTA inputs, LOSAT/src/blastinput/
// input_files.rs, so -query /dev/null finds it empty); a terminal or a FIFO is ESPIPE there.
// Every other open (from another directory descriptor, with other flags, of a directory, a socket
// or a block device, a symbolic link not followed, one that fails, with uvwasi's errno) stays with
// the calling thread's uvwasi.
const SHARED_FD_BASE = 0x40000000, SHARED_FDS = 1024;
const WASI_EBADF = 8, WASI_EINVAL = 28, WASI_EIO = 29, WASI_ENOTSUP = 58, WASI_ESPIPE = 70, WASI_ENOTCAPABLE = 76;
const WASI_ERRNO = { EACCES: 2, EAGAIN: 6, EBADF: WASI_EBADF, EDQUOT: 19, EEXIST: 20, EFBIG: 22, EINTR: 27, EINVAL: WASI_EINVAL,
  EIO: WASI_EIO, EISDIR: 31, ELOOP: 32, ENOENT: 44, ENOSPC: 51, ENOTDIR: 54, EOVERFLOW: 61, EPERM: 63, EPIPE: 64,
  EROFS: 69, ESPIPE: WASI_ESPIPE };
const RIGHT_FD_READ = 1n << 1n, RIGHT_FD_WRITE = 1n << 6n, RIGHT_FD_ALLOCATE = 1n << 8n, RIGHT_FD_FILESTAT_SET_SIZE = 1n << 22n;
const LOOKUP_SYMLINK_FOLLOW = 1, OFLAGS_CREAT = 1, OFLAGS_EXCL = 4, OFLAGS_TRUNC = 8, FDFLAGS_APPEND = 1, PREOPEN_ROOT_FD = 3;
// The filetypes of the table's files, as uvwasi gives them (a FIFO is a socket stream there).
const FILETYPE_CHARACTER_DEVICE = 2, FILETYPE_REGULAR_FILE = 4, FILETYPE_SOCKET_STREAM = 6;
// The table: per slot the OS descriptor (-1 free, -2 being entered), the position, the rights and
// the kind (the filetype, with KIND_APPEND for a file opened to append, KIND_OWN_POSITION once
// the guest seeked or read a regular file and KIND_TERMINAL for a terminal).
const KIND_APPEND = 0x100, KIND_OWN_POSITION = 0x200, KIND_TERMINAL = 0x400;
const STREAM_WAIT = new Int32Array(new SharedArrayBuffer(4));
const sharedFileViews = table => ({ os: new Int32Array(table, 0, SHARED_FDS),
  position: new BigInt64Array(table, 4 * SHARED_FDS, SHARED_FDS),
  rights: new BigUint64Array(table, 12 * SHARED_FDS, 2 * SHARED_FDS),
  kind: new Int32Array(table, 28 * SHARED_FDS, SHARED_FDS) });
function sharedFileTable() {
  const table = new SharedArrayBuffer(32 * SHARED_FDS);
  sharedFileViews(table).os.fill(-1);
  return table;
}

// The file functions of one thread with the shared table in front of its own uvwasi functions
// (`local`, which take the guest's signed i32 values themselves).
function sharedFileImports(table, memory, local) {
  const { os, position, rights, kind } = sharedFileViews(table);
  const slotOf = fd => (fd >>> 0) - SHARED_FD_BASE;
  const shared = fd => slotOf(fd) >= 0 && slotOf(fd) < SHARED_FDS;
  // A regular file has the table's position (the comment above SHARED_FD_BASE).
  const regular = slot => (Atomics.load(kind, slot) & 0xff) === FILETYPE_REGULAR_FILE;
  const ownPosition = slot => {
    if (!(Atomics.load(kind, slot) & KIND_OWN_POSITION)) Atomics.or(kind, slot, KIND_OWN_POSITION);
  };
  // A character device that is not a terminal: its offset is 0 (the comment above SHARED_FD_BASE).
  const atZero = slot => (Atomics.load(kind, slot) & (0xff | KIND_TERMINAL)) === FILETYPE_CHARACTER_DEVICE;
  // The guest's memory as this isolate sees it; memory.grow(0) updates the view when the guest
  // grew the memory on another thread (see the start of a thread below).
  const buffer = end => {
    if (end > memory.buffer.byteLength) memory.grow(0);
    return memory.buffer;
  };
  const view = (pointer, length) => new DataView(buffer(pointer + length), pointer, length);
  const lengthOf = (iovs, count) => {
    let total = 0n;
    const vectors = view(iovs, 8 * count);
    for (let i = 0; i < count; i++) total += BigInt(vectors.getUint32(8 * i + 4, true));
    return total;
  };
  // fs.readSync or fs.writeSync over the guest's vectors at `at` (null: the OS's position); an
  // error after some bytes ends the call with their count, as read() and write() do. writeSync
  // takes a Number (it writes at the OS's position when given a BigInt); readSync a BigInt.
  const transfer = (move, descriptor, iovs, count, at) => {
    let total = 0n;
    const vectors = view(iovs, 8 * count);
    for (let i = 0; i < count; i++) {
      const pointer = vectors.getUint32(8 * i, true), length = vectors.getUint32(8 * i + 4, true);
      if (length === 0) continue;
      let done;
      try {
        done = move(descriptor, new Uint8Array(buffer(pointer + length), pointer, length), 0, length,
          at === null ? null : move === fs.writeSync ? Number(at + total) : at + total);
      } catch (err) {
        if (total > 0n) break;
        throw err;
      }
      total += BigInt(done);
      if (done < length) break;
    }
    return total;
  };
  // fs.writevSync of the guest's vectors at the OS's position, as one writev(); a full stream
  // (EAGAIN) is written again after 1 ms.
  const writeAtOs = (descriptor, iovs, count) => {
    const vectors = view(iovs, 8 * count), buffers = [];
    for (let i = 0; i < count; i++) {
      const pointer = vectors.getUint32(8 * i, true), length = vectors.getUint32(8 * i + 4, true);
      if (length > 0) buffers.push(new Uint8Array(buffer(pointer + length), pointer, length));
    }
    if (buffers.length === 0) return 0n;
    for (;;) {
      try {
        return BigInt(fs.writevSync(descriptor, buffers, null));
      } catch (err) {
        if (!err || err.code !== "EAGAIN") throw err;
        Atomics.wait(STREAM_WAIT, 0, 0, 1);
      }
    }
  };
  const imports = {};
  // Any function given a shared descriptor where that is not one of the files: no right to it.
  for (const [name, call] of Object.entries(local)) {
    if (!/^(fd|path|sock)_/.test(name)) continue;
    const fds = { path_symlink: [2], path_rename: [0, 3], path_link: [0, 4], fd_renumber: [0, 1] }[name] || [0];
    imports[name] = (...args) => (fds.some(i => shared(args[i])) ? WASI_ENOTCAPABLE : call(...args));
  }
  const own = (name, op) => {
    const call = imports[name];
    imports[name] = (fd, ...args) => {
      if (!shared(fd)) return call(fd, ...args);
      const slot = slotOf(fd), descriptor = Atomics.load(os, slot);
      if (descriptor < 0) return WASI_EBADF;
      try {
        return op(slot, descriptor, ...args.map(arg => (typeof arg === "number" ? arg >>> 0 : arg))) ?? 0;
      } catch (err) {
        return WASI_ERRNO[err && err.code] ?? WASI_EIO;
      }
    };
  };
  const allowed = (slot, right) => (rights[2 * slot] & right) !== 0n;
  own("fd_read", (slot, descriptor, iovs, count, nread) => {
    if (!allowed(slot, RIGHT_FD_READ)) return WASI_ENOTCAPABLE;
    if (!regular(slot)) {
      view(nread, 4).setUint32(0, Number(transfer(fs.readSync, descriptor, iovs, count, null)), true);
      return;
    }
    ownPosition(slot);
    const got = transfer(fs.readSync, descriptor, iovs, count, Atomics.load(position, slot));
    Atomics.add(position, slot, got);
    view(nread, 4).setUint32(0, Number(got), true);
  });
  own("fd_pread", (slot, descriptor, iovs, count, offset, nread) => {
    if (!allowed(slot, RIGHT_FD_READ)) return WASI_ENOTCAPABLE;
    if (!regular(slot) && !atZero(slot)) return WASI_ESPIPE;
    view(nread, 4).setUint32(0, Number(transfer(fs.readSync, descriptor, iovs, count, BigInt.asUintN(64, offset))), true);
  });
  own("fd_write", (slot, descriptor, iovs, count, nwritten) => {
    if (!allowed(slot, RIGHT_FD_WRITE)) return WASI_ENOTCAPABLE;
    let wrote;
    const type = Atomics.load(kind, slot);
    if (type === (FILETYPE_REGULAR_FILE | KIND_OWN_POSITION)) {
      const length = lengthOf(iovs, count), at = Atomics.add(position, slot, length);
      try {
        wrote = transfer(fs.writeSync, descriptor, iovs, count, at);
      } finally {
        // The bytes not written given back, unless another write came after them.
        if (wrote !== length) Atomics.compareExchange(position, slot, at + length, at + (wrote ?? 0n));
      }
    } else {
      wrote = writeAtOs(descriptor, iovs, count);
      // An append leaves the position at the end of the file, as it leaves the native offset;
      // another write moves it by the bytes written, as the OS's.
      if (type & KIND_APPEND) Atomics.store(position, slot, fs.fstatSync(descriptor, { bigint: true }).size);
      else if ((type & 0xff) === FILETYPE_REGULAR_FILE) Atomics.add(position, slot, wrote);
    }
    view(nwritten, 4).setUint32(0, Number(wrote), true);
  });
  own("fd_pwrite", (slot, descriptor, iovs, count, offset, nwritten) => {
    if (!allowed(slot, RIGHT_FD_WRITE)) return WASI_ENOTCAPABLE;
    if (!regular(slot) && !atZero(slot)) return WASI_ESPIPE;
    view(nwritten, 4).setUint32(0, Number(transfer(fs.writeSync, descriptor, iovs, count, BigInt.asUintN(64, offset))), true);
  });
  own("fd_seek", (slot, descriptor, offset, whence, result) => {
    if (atZero(slot) && whence <= 2) {
      view(result, 8).setBigUint64(0, 0n, true);
      return;
    }
    if (!regular(slot)) return WASI_ESPIPE;
    ownPosition(slot);
    const from = whence === 0 ? 0n : whence === 1 ? Atomics.load(position, slot)
      : whence === 2 ? fs.fstatSync(descriptor, { bigint: true }).size : undefined;
    if (from === undefined || from + offset < 0n) return WASI_EINVAL;
    Atomics.store(position, slot, from + offset);
    view(result, 8).setBigUint64(0, from + offset, true);
  });
  own("fd_tell", (slot, descriptor, result) => {
    if (!regular(slot) && !atZero(slot)) return WASI_ESPIPE;
    view(result, 8).setBigUint64(0, atZero(slot) ? 0n : Atomics.load(position, slot), true);
  });
  own("fd_close", (slot, descriptor) => {
    if (Atomics.compareExchange(os, slot, descriptor, -1) !== descriptor) return WASI_EBADF;
    fs.closeSync(descriptor);
  });
  own("fd_filestat_get", (slot, descriptor, result) => {
    const stat = fs.fstatSync(descriptor, { bigint: true }), out = view(result, 64);
    for (const [at, value] of [[0, stat.dev], [8, stat.ino], [16, 0n], [24, stat.nlink], [32, stat.size],
      [40, stat.atimeNs], [48, stat.mtimeNs], [56, stat.ctimeNs]]) out.setBigUint64(at, value, true);
    out.setUint8(16, kind[slot] & 0xff);
  });
  own("fd_fdstat_get", (slot, descriptor, result) => {
    const out = view(result, 24);
    out.setBigUint64(0, 0n, true);
    out.setUint8(0, kind[slot] & 0xff);
    out.setUint16(2, kind[slot] & KIND_APPEND ? FDFLAGS_APPEND : 0, true);
    out.setBigUint64(8, rights[2 * slot], true);
    out.setBigUint64(16, rights[2 * slot + 1], true);
  });
  own("fd_fdstat_set_rights", (slot, descriptor, base, inheriting) => {
    [base, inheriting] = [BigInt.asUintN(64, base), BigInt.asUintN(64, inheriting)];
    if ((base & ~rights[2 * slot]) || (inheriting & ~rights[2 * slot + 1])) return WASI_ENOTCAPABLE;
    [rights[2 * slot], rights[2 * slot + 1]] = [base, inheriting];
  });
  own("fd_fdstat_set_flags", (slot, descriptor, flags) => (flags === (kind[slot] & KIND_APPEND ? FDFLAGS_APPEND : 0) ? 0 : WASI_ENOTSUP));
  own("fd_filestat_set_size", (slot, descriptor, size) => {
    if (!allowed(slot, RIGHT_FD_FILESTAT_SET_SIZE)) return WASI_ENOTCAPABLE;
    fs.ftruncateSync(descriptor, Number(BigInt.asUintN(64, size)));
  });
  // Grows the file to offset + length, as uvwasi does (Node has no posix_fallocate).
  own("fd_allocate", (slot, descriptor, offset, length) => {
    if (!allowed(slot, RIGHT_FD_ALLOCATE)) return WASI_ENOTCAPABLE;
    const end = BigInt.asUintN(64, offset) + BigInt.asUintN(64, length);
    if (fs.fstatSync(descriptor, { bigint: true }).size < end) fs.ftruncateSync(descriptor, Number(end));
  });
  own("fd_advise", () => 0);
  own("fd_sync", (slot, descriptor) => { fs.fsyncSync(descriptor); });
  own("fd_datasync", (slot, descriptor) => { fs.fdatasyncSync(descriptor); });
  own("fd_prestat_get", () => WASI_EBADF);
  own("fd_prestat_dir_name", () => WASI_EBADF);
  // A free slot marked as being entered, or -1 when the table is full.
  const reserve = () => {
    for (let slot = 0; slot < SHARED_FDS; slot++) {
      if (Atomics.load(os, slot) === -1 && Atomics.compareExchange(os, slot, -1, -2) === -1) return slot;
    }
    return -1;
  };
  const pathOpen = imports.path_open;
  imports.path_open = (...args) => {
    const [dirfd, dirflags, path, pathLength, oflags] = args.slice(0, 5).map(arg => arg >>> 0);
    const [base, inheriting, fdflags, opened] = [BigInt.asUintN(64, args[5]), BigInt.asUintN(64, args[6]), args[7] >>> 0, args[8] >>> 0];
    const write = (base & RIGHT_FD_WRITE) !== 0n;
    if (dirfd !== PREOPEN_ROOT_FD || (write
      ? (oflags & ~(OFLAGS_CREAT | OFLAGS_EXCL | OFLAGS_TRUNC)) || (fdflags & ~FDFLAGS_APPEND)
      : oflags !== 0 || fdflags !== 0)) return pathOpen(...args);
    let slot = -1, descriptor = -1;
    try {
      // The guest's bytes as they are (the path may create a file).
      const host = Buffer.concat([Buffer.from("/"), new Uint8Array(buffer(path + pathLength), path, pathLength)]);
      const follow = (dirflags & LOOKUP_SYMLINK_FOLLOW) !== 0;
      const stat = (follow ? fs.statSync : fs.lstatSync)(host, { throwIfNoEntry: false });
      const kept = stat => stat.isFile() || stat.isCharacterDevice() || stat.isFIFO();
      if (stat ? !kept(stat) : !(write && (oflags & OFLAGS_CREAT))) return pathOpen(...args);
      // A full table: the thread's own uvwasi opens it, as before (reserved before the open,
      // which may create the file).
      slot = reserve();
      if (slot < 0) return pathOpen(...args);
      const { O_RDONLY, O_WRONLY, O_RDWR, O_CREAT, O_EXCL, O_TRUNC, O_APPEND, O_NOFOLLOW } = fs.constants;
      descriptor = fs.openSync(host, (write ? (base & RIGHT_FD_READ ? O_RDWR : O_WRONLY) : O_RDONLY)
        | (oflags & OFLAGS_CREAT ? O_CREAT : 0) | (oflags & OFLAGS_EXCL ? O_EXCL : 0) | (oflags & OFLAGS_TRUNC ? O_TRUNC : 0)
        | (fdflags & FDFLAGS_APPEND ? O_APPEND : 0) | (follow ? 0 : O_NOFOLLOW));
      const now = fs.fstatSync(descriptor);
      if (!kept(now)) throw new Error("not a file of the table");
      Atomics.store(position, slot, 0n);
      [rights[2 * slot], rights[2 * slot + 1]] = [base, inheriting];
      Atomics.store(kind, slot, (now.isFile() ? FILETYPE_REGULAR_FILE : now.isFIFO() ? FILETYPE_SOCKET_STREAM
        : FILETYPE_CHARACTER_DEVICE) | (fdflags & FDFLAGS_APPEND ? KIND_APPEND : 0)
        | (now.isCharacterDevice() && tty.isatty(descriptor) ? KIND_TERMINAL : 0));
      Atomics.store(os, slot, descriptor);
      view(opened, 4).setUint32(0, SHARED_FD_BASE + slot, true);
      return 0;
    } catch {
      if (descriptor >= 0) fs.closeSync(descriptor);
      if (slot >= 0) Atomics.store(os, slot, -1);
      return pathOpen(...args);
    }
  };
  return imports;
}

module.exports = { sharedFileTable, sharedFileImports, sharedFileViews };
