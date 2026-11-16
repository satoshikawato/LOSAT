// Saving a file happens only on an explicit user action (design document §3.3).
export interface Downloader {
  save(fileName: string, bytes: Uint8Array, mimeType: string): void;
}

/**
 * Saving a file that this browser keeps (the saved copy of an Example set, line WF-3), on an
 * explicit user action: into a folder that the user picks where the browser can (the File
 * System Access API: Chromium, Edge), else as the browser's download. Nothing is sent.
 */
export interface FileSaver {
  /** Whether the browser writes into a folder that the user picks ("Save to folder…"). */
  readonly folders: boolean;
  /**
   * Asks for a folder, then writes the file there as `fileName`, replacing a file of that name.
   * `file` is called once the folder is picked (the picker opens on the user's click). Resolves
   * the folder's name, or undefined when the user closed the picker.
   */
  saveToFolder(fileName: string, file: () => Promise<Blob>): Promise<string | undefined>;
  /** The browser's download of the file. */
  download(fileName: string, file: Blob): void;
}
