// A test double of the engine's inflation (src/ports/inflate.ts) over Node's zlib, for the tests of
// the flows around it (the draft, the run coordinator). It names the compression by the suffix as
// the engine's rule does (`.gz`, `.bz2`, `.zst`, in any case) and inflates the whole file at its
// end. It is not the engine's decoder: tests/unit/engine-runtime.test.ts checks the reactor's bytes
// and sentences, and web/adapter/tests/inflate.rs the adapter's.
import { gunzipSync, zstdDecompressSync } from 'node:zlib';
import type { Compression, FileInflater, InflatedSink, Inflation } from '../../../src/ports/inflate';

const SUFFIXES: readonly (readonly [string, Compression])[] = [
  ['.gz', 'gzip'],
  ['.bz2', 'bzip2'],
  ['.zst', 'zstd'],
];

export class NodeInflater implements FileInflater {
  async open(name: string, sink: InflatedSink): Promise<Inflation | undefined> {
    const found = SUFFIXES.find(([suffix]) => name.toLowerCase().endsWith(suffix));
    if (found === undefined) return undefined;
    const [suffix, compression] = found;
    const parts: Uint8Array[] = [];
    return {
      compression,
      plainName: name.slice(0, name.length - suffix.length) || 'input',
      push: async (bytes) => {
        parts.push(bytes.slice());
      },
      end: async () => {
        if (compression === 'bzip2') throw new Error('NodeInflater: no bzip2 in Node');
        const file = Buffer.concat(parts);
        sink(compression === 'gzip' ? gunzipSync(file) : zstdDecompressSync(file));
      },
      cancel: async () => undefined,
    };
  }
}
