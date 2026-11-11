// Inflation port (DW-45, line GZ; swebok-improve SW8 D2): a source file named `.gz`, `.bz2` or `.zst`
// is inflated by the engine's own decoders, the ones the command line reads such a file with (ABI v2
// `inflate_begin` / `inflate_chunk` / `inflate_end`, docs/web/abi_v2.md §4, §5, §9). The engine's
// name rule decides the compression and the name of the inflated file; the application holds no
// copy of it.
//
// The real implementation is the adapter's serial reactor inside the Data worker; a build without
// the engine reads every file as it is.

/** The formats that the engine inflates, as `inflate_begin` names them. */
export type Compression = 'gzip' | 'bzip2' | 'zstd';

/** Receives inflated bytes; `bytes` is valid only during the call. */
export type InflatedSink = (bytes: Uint8Array) => void;

/** The inflation of one file, from `open` to `end` (or `cancel`). */
export interface Inflation {
  /** The format that the file's name announces. */
  readonly compression: Compression;
  /** The name without its suffix (`input` when nothing is left), for the inflated file. */
  readonly plainName: string;
  /**
   * Inflates the next bytes of the file; the inflated bytes go to the sink during the call.
   * Rejects with the engine's sentence for data that is not the named format, and then ends.
   */
  push(bytes: Uint8Array): Promise<void>;
  /** The end of the file: the rest of the inflated bytes, or the engine's sentence for data cut short or an empty file. */
  end(): Promise<void>;
  /** Frees an inflation that stops before its end; does nothing once it has ended or failed. */
  cancel(): Promise<void>;
}

export interface FileInflater {
  /** The inflation of the file named `name`, whose bytes go to `sink`, or undefined when the name is plain. */
  open(name: string, sink: InflatedSink): Promise<Inflation | undefined>;
}
