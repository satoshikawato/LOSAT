// A Downloader for tests: it keeps each saved file in memory, with the number of blocks it was
// written in, and saves nothing for a file that was aborted (ports/download.ts).
import type { Downloader, ExportSink } from '../../../src/ports/download';

export interface SavedFile {
  readonly name: string;
  readonly mime: string;
  readonly bytes: Uint8Array;
  /** The writes that made the file. */
  readonly blocks: number;
}

export function memoryDownloader(onSave: (file: SavedFile) => void): Downloader {
  return {
    open(name, mime): ExportSink {
      const parts: Uint8Array[] = [];
      let state: 'open' | 'closed' | 'aborted' = 'open';
      return {
        async write(bytes) {
          if (state !== 'open') throw new Error(`the file ${name} is no longer open`);
          parts.push(bytes.slice());
        },
        async close() {
          if (state !== 'open') throw new Error(`the file ${name} is no longer open`);
          state = 'closed';
          const bytes = new Uint8Array(parts.reduce((sum, part) => sum + part.length, 0));
          let offset = 0;
          for (const part of parts) {
            bytes.set(part, offset);
            offset += part.length;
          }
          onSave({ name, mime, bytes, blocks: parts.length });
        },
        abort() {
          if (state === 'open') state = 'aborted';
        },
      };
    },
  };
}
