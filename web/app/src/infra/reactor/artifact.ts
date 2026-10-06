// Reads the kind and the memory limits of an engine module from its bytes, before it is
// compiled (docs/web/abi_v2.md §2-§4). The JavaScript API does not expose the limits of
// an imported memory, and the threaded host must create the shared memory with exactly
// the maximum that the module declares (plan TD-7).

export type ReactorKind = 'serial-reactor' | 'threaded-reactor';

export interface MemoryLimits {
  /** In 64 KiB pages. */
  readonly initial: number;
  readonly maximum?: number;
  readonly shared: boolean;
}

export interface ReactorArtifact {
  readonly kind: ReactorKind;
  readonly memory: MemoryLimits;
}

/** The exports that both reactors must have (ABI v2 §4). */
export const ABI_EXPORTS: readonly string[] = Object.freeze([
  '_initialize',
  'memory',
  'losat_web2_abi_version',
  'losat_web2_alloc',
  'losat_web2_dealloc',
  'losat_web2_last_error_ptr',
  'losat_web2_last_error_len',
  'losat_web2_describe',
  'losat_web2_validate',
  'losat_web2_register',
  'losat_web2_release',
  'losat_web2_scan_begin',
  'losat_web2_scan_chunk',
  'losat_web2_scan_end',
  'losat_web2_run',
]);

class Reader {
  position = 0;
  constructor(private readonly bytes: Uint8Array) {}

  get done(): boolean {
    return this.position >= this.bytes.length;
  }

  byte(): number {
    if (this.position >= this.bytes.length) throw new Error('the engine module is truncated');
    return this.bytes[this.position++]!;
  }

  u32(): number {
    let value = 0;
    for (let shift = 0; shift < 35; shift += 7) {
      const byte = this.byte();
      value += (byte & 0x7f) * 2 ** shift;
      if ((byte & 0x80) === 0) return value;
    }
    throw new Error('the engine module has an invalid integer');
  }

  name(): string {
    const length = this.u32();
    const start = this.position;
    this.position += length;
    if (this.position > this.bytes.length) throw new Error('the engine module is truncated');
    return new TextDecoder().decode(this.bytes.subarray(start, this.position));
  }

  limits(): MemoryLimits {
    const flags = this.u32();
    if (flags & ~0x3) throw new Error('the engine module uses memory limits that LOSAT Web does not support');
    const initial = this.u32();
    const maximum = flags & 0x1 ? this.u32() : undefined;
    return { initial, ...(maximum === undefined ? {} : { maximum }), shared: (flags & 0x2) !== 0 };
  }

  skip(length: number): void {
    this.position += length;
  }
}

/** Checks that `bytes` is an ABI v2 reactor and returns its kind and memory limits. */
export function inspectReactor(bytes: Uint8Array): ReactorArtifact {
  const header = [0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00];
  if (bytes.length < 8 || header.some((value, i) => bytes[i] !== value)) {
    throw new Error('the engine file is not a WebAssembly module');
  }
  const reader = new Reader(bytes);
  reader.skip(8);
  const functionImports = new Set<string>();
  const exports = new Set<string>();
  let imported: MemoryLimits | undefined;
  let own: MemoryLimits | undefined;
  while (!reader.done) {
    const id = reader.byte();
    const size = reader.u32();
    const end = reader.position + size;
    if (id === 2) {
      const count = reader.u32();
      for (let i = 0; i < count; i++) {
        const module = reader.name();
        const field = reader.name();
        const kind = reader.byte();
        if (kind === 0) {
          reader.u32();
          functionImports.add(`${module}.${field}`);
        } else if (kind === 1) {
          reader.byte();
          reader.limits();
        } else if (kind === 2) {
          if (module !== 'env' || field !== 'memory') throw new Error(`the engine module imports memory ${module}.${field}`);
          imported = reader.limits();
        } else if (kind === 3) {
          reader.byte();
          reader.byte();
        } else {
          throw new Error(`the engine module has an import of kind ${kind}`);
        }
      }
    } else if (id === 5) {
      const count = reader.u32();
      if (count !== 1) throw new Error('the engine module must define one memory');
      own = reader.limits();
    } else if (id === 7) {
      const count = reader.u32();
      for (let i = 0; i < count; i++) {
        exports.add(reader.name());
        reader.byte();
        reader.u32();
      }
    }
    reader.position = end;
  }
  const missing = ABI_EXPORTS.filter((name) => !exports.has(name));
  if (missing.length > 0) throw new Error(`the engine module lacks the exports ${missing.join(', ')}`);
  const threaded = functionImports.has('wasi.thread-spawn');
  if (threaded) {
    if (imported === undefined || !imported.shared || imported.maximum === undefined || own !== undefined) {
      throw new Error('the threaded engine module must import one shared memory with a maximum');
    }
    if (!exports.has('wasi_thread_start')) throw new Error('the threaded engine module lacks wasi_thread_start');
    return { kind: 'threaded-reactor', memory: imported };
  }
  if (own === undefined || own.shared || imported !== undefined || exports.has('wasi_thread_start')) {
    throw new Error('the serial engine module must define one memory that is not shared');
  }
  return { kind: 'serial-reactor', memory: own };
}
