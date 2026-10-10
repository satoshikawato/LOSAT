// Compression port (plan §5.8): session files are gzip (design §12.2), written and read with the
// browser's CompressionStream and DecompressionStream, without a dependency. Both directions
// stream: neither holds a whole file.
import type { ExportSink } from './download';

export interface Compression {
  /**
   * A sink that gzips the blocks written to it and hands each compressed block, in order, to
   * `out`. `close` ends the gzip stream (its last blocks reach `out` before it resolves); `abort`
   * drops it. A failure of `out` fails the next `write` or `close`.
   */
  gzip(out: (bytes: Uint8Array) => Promise<void>): ExportSink;
  /**
   * The decompressed bytes of a gzip file, in chunks, in order. The iteration throws when the
   * gzip data is damaged or cut short (often only at its end, where the checksum is); ending the
   * iteration early stops reading the file.
   */
  gunzip(file: Blob): AsyncIterable<Uint8Array>;
}
