// The inflation of the adapter's serial reactor (ABI v2 `inflate_begin` / `inflate_chunk` /
// `inflate_end`, docs/web/abi_v2.md §4, §5, §9): the FileInflater of the Data worker. The
// engine names the compression and the plain name, inflates with the command line's decoders,
// and a failure is the command line's sentence.
//
// Each inflation runs on an instance of its own, opened for that file and dropped when the
// inflation ends (end, cancel, error or trap): what the decoders grow in the linear memory (a
// zstd window of up to 128 MiB, with the buffers of the file's slices) cannot be given back
// while the instance lives, so it goes with the instance and not into the session's shared one.
import type { Compression, FileInflater, InflatedSink, Inflation } from '../../ports/inflate';
import type { BegunInflation, ReactorAbi } from './abi';

export class ReactorInflater implements FileInflater {
  /**
   * `shared` is the session's instance (the scanner's and the checker's too); it only reads the
   * engine's name rule, so that a plain name opens no instance. `instantiate` opens a new
   * instance of the serial module; each inflation takes one.
   */
  constructor(
    private readonly shared: () => Promise<ReactorAbi>,
    private readonly instantiate: () => Promise<ReactorAbi>,
  ) {}

  async open(name: string, sink: InflatedSink): Promise<Inflation | undefined> {
    if (!(await this.named(name))) return undefined;
    const abi = await this.instantiate();
    const begun = abi.inflateBegin(name);
    // The engine named the file a moment ago: the same name gives an inflation again.
    return begun === undefined ? undefined : new ReactorInflation(abi, begun, sink);
  }

  /** Whether the engine's rule names the file as compressed: `inflate_begin` returns a handle. */
  private async named(name: string): Promise<boolean> {
    const abi = await this.shared();
    const begun = abi.inflateBegin(name);
    if (begun === undefined) return false;
    try {
      // Nothing was pushed: the end of an empty file fails with its sentence, which is not asked here.
      abi.inflateEnd(begun.handle, () => undefined);
    } catch {
      // The inflation is dropped whatever its outcome (§4).
    }
    return true;
  }
}

/**
 * One inflation, alone on the instance that began it. The instance is let go when the inflation
 * ends, fails or is cancelled (a trap ends it with the instance), so that the memory it grew is
 * collected with it.
 */
class ReactorInflation implements Inflation {
  readonly compression: Compression;
  readonly plainName: string;
  private readonly handle: number;
  private abi: ReactorAbi | undefined;

  constructor(
    abi: ReactorAbi,
    begun: BegunInflation,
    private readonly sink: InflatedSink,
  ) {
    this.abi = abi;
    this.handle = begun.handle;
    this.compression = begun.compression;
    this.plainName = begun.plainName;
  }

  async push(bytes: Uint8Array): Promise<void> {
    if (this.abi === undefined) throw new Error('the inflation has ended');
    try {
      this.abi.inflateChunk(this.handle, bytes, this.sink);
    } catch (error) {
      // The engine drops an inflation that fails.
      this.abi = undefined;
      throw error;
    }
  }

  async end(): Promise<void> {
    const abi = this.abi;
    if (abi === undefined) throw new Error('the inflation has ended');
    // `inflate_end` drops the inflation whatever its outcome.
    this.abi = undefined;
    abi.inflateEnd(this.handle, this.sink);
  }

  async cancel(): Promise<void> {
    const abi = this.abi;
    if (abi === undefined) return;
    this.abi = undefined;
    if (abi.stoppedBy !== undefined) return;
    try {
      abi.inflateEnd(this.handle, () => undefined);
    } catch {
      // A file read only in part fails to end; the reason it stopped is the one to report.
    }
  }
}
