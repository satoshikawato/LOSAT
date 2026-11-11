// The inflation of the adapter's serial reactor (ABI v2 `inflate_begin` / `inflate_chunk` /
// `inflate_end`, docs/web/abi_v2.md §4, §5, §9): the FileInflater of the Data worker. The
// engine names the compression and the plain name, inflates with the command line's decoders,
// and a failure is the command line's sentence.
import type { Compression, FileInflater, InflatedSink, Inflation } from '../../ports/inflate';
import type { BegunInflation, ReactorAbi } from './abi';

export class ReactorInflater implements FileInflater {
  constructor(private readonly reactor: () => Promise<ReactorAbi>) {}

  async open(name: string, sink: InflatedSink): Promise<Inflation | undefined> {
    const abi = await this.reactor();
    const begun = abi.inflateBegin(name);
    return begun === undefined ? undefined : new ReactorInflation(abi, begun, sink);
  }
}

/** One inflation on the instance that began it (a trap ends it with the instance). */
class ReactorInflation implements Inflation {
  readonly compression: Compression;
  readonly plainName: string;
  private readonly handle: number;
  private open = true;

  constructor(
    private readonly abi: ReactorAbi,
    begun: BegunInflation,
    private readonly sink: InflatedSink,
  ) {
    this.handle = begun.handle;
    this.compression = begun.compression;
    this.plainName = begun.plainName;
  }

  async push(bytes: Uint8Array): Promise<void> {
    if (!this.open) throw new Error('the inflation has ended');
    try {
      this.abi.inflateChunk(this.handle, bytes, this.sink);
    } catch (error) {
      // The engine drops an inflation that fails.
      this.open = false;
      throw error;
    }
  }

  async end(): Promise<void> {
    if (!this.open) throw new Error('the inflation has ended');
    // `inflate_end` drops the inflation whatever its outcome.
    this.open = false;
    this.abi.inflateEnd(this.handle, this.sink);
  }

  async cancel(): Promise<void> {
    if (!this.open) return;
    this.open = false;
    if (this.abi.stoppedBy !== undefined) return;
    try {
      this.abi.inflateEnd(this.handle, () => undefined);
    } catch {
      // A file read only in part fails to end; the reason it stopped is the one to report.
    }
  }
}
