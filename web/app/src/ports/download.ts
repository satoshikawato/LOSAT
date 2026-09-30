// Saving a file happens only on an explicit user action (design document §3.3).
export interface Downloader {
  save(fileName: string, bytes: Uint8Array, mimeType: string): void;
}
