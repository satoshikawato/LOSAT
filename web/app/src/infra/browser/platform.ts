// Browser implementations of small platform services used by the composition root.
import type { Downloader } from '../../ports/download';

export const browserDownloader: Downloader = {
  save(fileName, bytes, mimeType) {
    const url = URL.createObjectURL(new Blob([bytes as BlobPart], { type: mimeType }));
    const link = document.createElement('a');
    link.href = url;
    link.download = fileName;
    link.click();
    // Revoking immediately can cancel the download in some browsers.
    setTimeout(() => URL.revokeObjectURL(url), 60_000);
  },
};

export async function sha256Hex(bytes: Uint8Array): Promise<string> {
  const digest = await crypto.subtle.digest('SHA-256', bytes as BufferSource);
  return Array.from(new Uint8Array(digest), (b) => b.toString(16).padStart(2, '0')).join('');
}
