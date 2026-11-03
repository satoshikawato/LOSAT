// Browser implementations of small platform services used by the composition root.
import type { Downloader, FileSaver } from '../../ports/download';

export const browserDownloader: Downloader = {
  save(fileName, bytes, mimeType) {
    downloadBlob(fileName, new Blob([bytes as BlobPart], { type: mimeType }));
  },
};

/** The File System Access API's folder picker, which only Chromium and Edge have (not on Android). */
type PickFolder = (options: { readonly id?: string; readonly mode?: 'read' | 'readwrite' }) => Promise<FileSystemDirectoryHandle>;

function folderPicker(): PickFolder | undefined {
  const pick = (globalThis as { showDirectoryPicker?: PickFolder }).showDirectoryPicker;
  return typeof pick === 'function' ? pick : undefined;
}

export const browserFileSaver: FileSaver = {
  get folders() {
    return folderPicker() !== undefined;
  },
  async saveToFolder(fileName, file) {
    const pick = folderPicker();
    if (pick === undefined) throw new Error('this browser cannot write into a folder');
    let folder: FileSystemDirectoryHandle;
    try {
      // The id makes the picker open at the folder chosen last time.
      folder = await pick.call(globalThis, { id: 'losat-examples', mode: 'readwrite' });
    } catch (error) {
      if (error instanceof DOMException && error.name === 'AbortError') return undefined;
      throw error;
    }
    const data = await file();
    const handle = await folder.getFileHandle(fileName, { create: true });
    // The writable replaces the file when it closes; a failed copy aborts it and leaves the old file.
    await data.stream().pipeTo(await handle.createWritable());
    return folder.name;
  },
  download: downloadBlob,
};

function downloadBlob(fileName: string, blob: Blob): void {
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = fileName;
  link.click();
  // Revoking immediately can cancel the download in some browsers.
  setTimeout(() => URL.revokeObjectURL(url), 60_000);
}

export async function sha256Hex(bytes: Uint8Array): Promise<string> {
  const digest = await crypto.subtle.digest('SHA-256', bytes as BufferSource);
  return Array.from(new Uint8Array(digest), (b) => b.toString(16).padStart(2, '0')).join('');
}
