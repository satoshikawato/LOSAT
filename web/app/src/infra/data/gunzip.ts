// A source file named `.gz` is inflated when it is opened (DW-45, line GZ). The rule is the one of
// NCBI's makeblastdb: the name ends in `.gz`, in any case. The inflated bytes are kept as a new File
// under the name without `.gz`, so the engine never sees a plain text under a `.gz` name.
// fflate's streaming Gunzip does the work and reads every member of a concatenated file (bgzip):
// `DecompressionStream('gzip')` rejects concatenated members in Chromium ("Failed to fetch"),
// Firefox ("The operation was aborted") and WebKit ("Extra bytes past the end"), checked 2026-11-03.
import { Gunzip } from 'fflate';

/** Bytes read from the file at once, and inflated bytes gathered into one Blob part. */
const SLICE_BYTES = 8 * 1024 * 1024;
const PART_BYTES = 8 * 1024 * 1024;
const GZ_SUFFIX = /\.gz$/i;

export const isGzipName = (name: string): boolean => GZ_SUFFIX.test(name);

/** The name without `.gz`; `input` when nothing is left. */
export function plainName(name: string): string {
  const plain = name.replace(GZ_SUFFIX, '');
  return plain === '' ? 'input' : plain;
}

export async function gunzipFile(file: File): Promise<File> {
  const parts: Blob[] = [];
  let pending: Uint8Array[] = [];
  let pendingBytes = 0;
  const flush = () => {
    if (pendingBytes === 0) return;
    parts.push(new Blob(pending as BlobPart[]));
    pending = [];
    pendingBytes = 0;
  };
  const gunzip = new Gunzip((chunk) => {
    // fflate may reuse a buffer later, so the chunk is copied.
    pending.push(chunk.slice());
    pendingBytes += chunk.length;
    if (pendingBytes >= PART_BYTES) flush();
  });
  try {
    if (file.size === 0) throw new Error('the file is empty');
    for (let start = 0; start < file.size; start += SLICE_BYTES) {
      const end = Math.min(start + SLICE_BYTES, file.size);
      gunzip.push(new Uint8Array(await file.slice(start, end).arrayBuffer()), end === file.size);
    }
  } catch (error) {
    const message = error instanceof Error ? error.message : String(error);
    throw new Error(`${file.name} is named .gz but is not gzip data (${message})`);
  }
  flush();
  return new File(parts, plainName(file.name), { type: 'text/plain' });
}
