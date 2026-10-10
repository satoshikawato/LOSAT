// Saving a file happens only on an explicit user action (design document §3.3). Every export
// is written as blocks, in order, into a sink (design §12.1, the Writer contract): the
// formatters, the extraction and the session file never build a whole file in memory
// (application/export-writer.ts collects text into blocks).

/** Where one export is written, block by block, in order. */
export interface ExportSink {
  /** Appends a block. The caller may reuse `bytes` once the promise settles. */
  write(bytes: Uint8Array): Promise<void>;
  /** Ends the file and hands it to the browser's download. */
  close(): Promise<void>;
  /** Drops what was written: nothing is saved. Safe after a failure, after `close`, and more than once. */
  abort(): void;
}

export interface Downloader {
  /** Starts a file; nothing is saved before `close`. */
  open(fileName: string, mimeType: string): ExportSink;
}
