// The names of the files a test leaves in its output directory. CI uploads web/app/test-results/ when a
// job fails (web.yml, actions/upload-artifact), and the upload refuses the whole artifact if one file's
// path has one of nine characters. The random walk once wrote its native run's inputs under their own
// names (`<svg onload="…">.fa`), and every failed shard that ran it then kept no trace, error context or
// screenshot of any test (line WQA's D24, 2026-11-27). Evidence that must keep such a name is written
// under uploadableName() with the real name beside it, and every test fails when it leaves such a name
// (support/test.ts), after renaming it so that the upload still keeps everything.
import { existsSync, readdirSync, renameSync, writeFileSync } from 'node:fs';
import { basename, dirname, join, relative } from 'node:path';

/** The characters actions/upload-artifact refuses in a file's path: " : < > | * ? CR LF. */
const REFUSED = /[":<>|*?\r\n]/;

/** `name` with each character the upload refuses, and `%`, as `%XX` (its UTF-8 bytes): a name the upload keeps. */
export function uploadableName(name: string): string {
  return name.replace(/[":<>|*?\r\n%]/g, (character) => `%${character.charCodeAt(0).toString(16).toUpperCase().padStart(2, '0')}`);
}

/** The files and folders under `directory` whose names the upload refuses, deepest first, as absolute paths. */
function refusedEntries(directory: string): string[] {
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
      if (entry.isDirectory()) visit(path);
      if (REFUSED.test(entry.name)) found.push(path);
    }
  };
  visit(directory);
  return found;
}

/** The paths, relative to `directory`, of the files and folders under it whose names the upload refuses. */
export function unuploadableNames(directory: string): string[] {
  return refusedEntries(directory).map((path) => relative(directory, path));
}

/** The file, in a test's output directory, that maps renamed entries back to the names the test gave them. */
export const RENAMED_NAMES_FILE = 'renamed-names.json';

/**
 * Renames every file and folder under `directory` whose name the upload refuses to its uploadableName()
 * (with a number before the extension if that name is taken), and writes RENAMED_NAMES_FILE there: for
 * each, its path as the test left it (relative to `directory`) and the name it got. Returns those paths.
 */
export function renameUnuploadable(directory: string): string[] {
  const entries = refusedEntries(directory);
  if (entries.length === 0) return [];
  const renamed: Array<{ path: string; renamed: string }> = [];
  // Deepest first: a folder is renamed after the entries in it, whose paths then still hold.
  for (const path of entries) {
    const folder = dirname(path);
    const name = basename(path);
    let kept = uploadableName(name);
    for (let n = 2; existsSync(join(folder, kept)); n++) kept = uploadableName(name).replace(/(\.[^.]*)?$/, (extension) => `-${n}${extension}`);
    renameSync(path, join(folder, kept));
    renamed.push({ path: relative(directory, path), renamed: kept });
  }
  writeFileSync(join(directory, RENAMED_NAMES_FILE), `${JSON.stringify(renamed, null, 1)}\n`);
  return renamed.map((entry) => entry.path);
}
