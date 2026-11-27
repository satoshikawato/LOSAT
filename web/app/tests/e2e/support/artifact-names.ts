// The names of the files a test leaves in its output directory. CI uploads web/app/test-results/ when a
// job fails (web.yml, actions/upload-artifact), and the upload refuses a whole artifact if one file's
// path has a character NTFS does not allow. The random walk once wrote its native run's inputs under
// their own names (`<svg onload="…">.fa`), and every failed shard that ran it then kept no trace,
// error context or screenshot of any test (line WQA's D24, 2026-11-27). Every test therefore fails
// when it leaves such a name (support/test.ts), and evidence that must keep such a name is written
// under uploadableName() with the real name beside it.
import { readdirSync } from 'node:fs';
import { join, relative } from 'node:path';

/** The characters actions/upload-artifact refuses in a file's path (" : < > | * ? CR LF). */
const REFUSED = /[":<>|*?\r\n]/;

/** `name` with each character the upload refuses, and `%`, as `%XX` (its UTF-8 bytes): a name the upload keeps. */
export function uploadableName(name: string): string {
  return name.replace(/[":<>|*?\r\n%]/g, (character) => `%${character.charCodeAt(0).toString(16).toUpperCase().padStart(2, '0')}`);
}

/** The paths, relative to `directory`, of the files and folders under it whose names the upload refuses. */
export function unuploadableNames(directory: string): string[] {
  const found: string[] = [];
  const visit = (folder: string): void => {
    let entries;
    try {
      entries = readdirSync(folder, { withFileTypes: true });
    } catch (error) {
      if ((error as NodeJS.ErrnoException).code === 'ENOENT') return;
      throw error;
    }
    for (const entry of entries) {
      const path = join(folder, entry.name);
      if (REFUSED.test(entry.name)) found.push(relative(directory, path));
      if (entry.isDirectory()) visit(path);
    }
  };
  visit(directory);
  return found;
}
