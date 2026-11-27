// The native command line's makeblastdb runs that the e2e compares the browser's runs with (line NP WW4b1, line WW
// WW4b2): a run in a folder of its inputs with the browser's zone (`TZ`), its paths written as the engine's `/work`; what
// may differ between two runs of the same input (the seconds of the times that makeblastdb prints and stores as
// MM/DD/YYYY HH:MM:SS, and the elapsed time) masked; and the browser's run equal to the native run just before it or
// just after it, so that the minute of the index file's date and of the metadata file's time stays compared (a minute
// may pass between two runs, never two).
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { mkdirSync, mkdtempSync, readdirSync, readFileSync, rmSync, statSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join } from 'node:path';
import { NATIVE } from './native';

/** What a run wrote: its exit status, its streams and the files of its folder that it made or changed, by name. */
export interface Written {
  readonly status: number;
  readonly stdout: Buffer;
  readonly stderr: Buffer;
  readonly files: ReadonlyMap<string, Buffer>;
}

/** The files under `folder`, by their paths in it, in the byte order of the names. */
export function treeOf(folder: string, prefix = ''): Map<string, Buffer> {
  const files = new Map<string, Buffer>();
  for (const name of readdirSync(join(folder, prefix)).sort()) {
    const path = prefix === '' ? name : `${prefix}/${name}`;
    if (statSync(join(folder, path)).isDirectory()) for (const [inner, bytes] of treeOf(folder, path)) files.set(inner, bytes);
    else files.set(path, readFileSync(join(folder, path)));
  }
  return files;
}

/** A new folder holding `files` (each by its path in the folder). */
export function folderOf(files: ReadonlyMap<string, Buffer>): string {
  const folder = mkdtempSync(join(tmpdir(), 'losat-makeblastdb-'));
  for (const [name, bytes] of files) {
    mkdirSync(dirname(join(folder, name)), { recursive: true });
    writeFileSync(join(folder, name), bytes);
  }
  return folder;
}

/** The native command line's run of `argv` in a folder of `inputs` with `TZ=zone`, its paths written as `/work`. */
export function nativeMakeblastdb(argv: readonly string[], inputs: ReadonlyMap<string, Buffer>, zone: string): Written {
  const folder = folderOf(inputs);
  try {
    const run = spawnSync(NATIVE!, argv, { cwd: folder, env: { ...process.env, TZ: zone }, maxBuffer: 1 << 28 });
    if (run.error) throw run.error;
    const place = (bytes: Buffer) => Buffer.from(bytes.toString('latin1').split(folder).join('/work'), 'latin1');
    const files = new Map<string, Buffer>();
    for (const [name, bytes] of treeOf(folder)) if (!inputs.get(name)?.equals(bytes)) files.set(name, place(bytes));
    return { status: run.status ?? -1, stdout: place(run.stdout), stderr: place(run.stderr), files };
  } finally {
    rmSync(folder, { recursive: true, force: true });
  }
}

/** The bytes with the seconds of the times that makeblastdb prints and stores (MM/DD/YYYY HH:MM:SS) and its elapsed time masked. */
export function withoutSeconds(bytes: Buffer): string {
  return bytes
    .toString('latin1')
    .replace(/(\d{2}\/\d{2}\/\d{4} \d{2}:\d{2}):\d{2}/g, '$1:##')
    .replace(/( sequences? in )[0-9.e+-]+( seconds)/g, '$1#$2');
}

/** Each file as its name and the SHA-256 of its bytes `withoutSeconds`, in the order of the names. */
export function maskedFiles(files: ReadonlyMap<string, Buffer>): string[] {
  return [...files].map(([name, bytes]) => `${name} ${createHash('sha256').update(withoutSeconds(bytes), 'latin1').digest('hex')}`).sort();
}

/** A run as its exit status, its streams `withoutSeconds` and `maskedFiles`. */
export function maskedOf(written: Written): string[] {
  return [`status ${written.status}`, `stdout ${withoutSeconds(written.stdout)}`, `stderr ${withoutSeconds(written.stderr)}`, ...maskedFiles(written.files)];
}

/**
 * The native run (masked) that the browser's run (masked) is: the one of `natives` (the runs just before and just after
 * the browser's) equal to it, or the first, to show the difference.
 */
export function asNative(browser: readonly string[], natives: readonly (readonly string[])[]): readonly string[] {
  return natives.find((native) => JSON.stringify(native) === JSON.stringify(browser)) ?? natives[0]!;
}
