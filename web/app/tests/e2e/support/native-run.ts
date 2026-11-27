// What the real-engine tests of the tool tabs share (makeblastdb-tool.spec.ts, masking-tool.spec.ts): the files that a
// native `losat` command line made in a folder of the same input files, a tar archive read back, and the bytes of
// makeblastdb's files with the seconds of their dates masked (the one tolerance: the engine writes the dates of a
// database in the browser's time zone, the native run uses the same zone, and the seconds differ between two runs).
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { mkdirSync, mkdtempSync, readdirSync, readFileSync, rmSync, statSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { dirname, join } from 'node:path';
import type { Page } from '@playwright/test';
import { NATIVE } from './native';

/** The browser's time zone of the real-engine tests (`test.use({ timezoneId })`); the native runs use it as TZ. */
export const ZONE = 'Asia/Tokyo';

/** The names and bytes of the files in a tar archive (ustar, with pax `path` records). */
export function untar(bytes: Buffer): Map<string, Buffer> {
  const files = new Map<string, Buffer>();
  let at = 0;
  let paxPath: string | undefined;
  while (at + 512 <= bytes.length) {
    const header = bytes.subarray(at, at + 512);
    if (header.every((byte) => byte === 0)) break;
    const field = (from: number, length: number) => header.subarray(from, from + length).toString('utf8').replace(/\0.*$/s, '');
    const size = Number.parseInt(field(124, 12).trim(), 8);
    const type = String.fromCharCode(header[156]!);
    const body = bytes.subarray(at + 512, at + 512 + size);
    if (type === 'x') {
      paxPath = /\d+ path=(.*)\n/.exec(body.toString('utf8'))?.[1];
    } else {
      const prefix = field(345, 155);
      const name = paxPath ?? (prefix === '' ? field(0, 100) : `${prefix}/${field(0, 100)}`);
      files.set(name, Buffer.from(body));
      paxPath = undefined;
    }
    at += 512 + Math.ceil(size / 512) * 512;
  }
  return files;
}

export async function downloadAll(page: Page): Promise<{ readonly name: string; readonly files: Map<string, Buffer> }> {
  const [download] = await Promise.all([page.waitForEvent('download'), page.getByTestId('tool-made-files-download-all').click()]);
  return { name: download.suggestedFilename(), files: untar(readFileSync((await download.path())!)) };
}


/** The files under `folder`, by their paths in it. */
export function treeOf(folder: string, prefix = ''): Map<string, Buffer> {
  const files = new Map<string, Buffer>();
  for (const name of readdirSync(join(folder, prefix)).sort()) {
    const path = prefix === '' ? name : `${prefix}/${name}`;
    if (statSync(join(folder, path)).isDirectory()) for (const [inner, bytes] of treeOf(folder, path)) files.set(inner, bytes);
    else files.set(path, readFileSync(join(folder, path)));
  }
  return files;
}

export function folderOf(files: ReadonlyMap<string, Buffer>): string {
  const folder = mkdtempSync(join(process.env['TMPDIR'] ?? tmpdir(), 'losat-mk-tab-'));
  for (const [name, bytes] of files) {
    mkdirSync(dirname(join(folder, name)), { recursive: true });
    writeFileSync(join(folder, name), bytes);
  }
  return folder;
}

/** The bytes with the seconds of makeblastdb's times (MM/DD/YYYY HH:MM:SS) and its elapsed time masked. */
export function withoutSeconds(bytes: Buffer): string {
  return bytes
    .toString('latin1')
    .replace(/(\d{2}\/\d{2}\/\d{4} \d{2}:\d{2}):\d{2}/g, '$1:##')
    .replace(/(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}):\d{2}/g, '$1:##')
    .replace(/( sequences? in )[0-9.e+-]+( seconds)/g, '$1#$2');
}

export const masked = (files: ReadonlyMap<string, Buffer>) =>
  [...files].map(([name, bytes]) => `${name} ${createHash('sha256').update(withoutSeconds(bytes), 'latin1').digest('hex')}`).sort();

/**
 * The standard output and the files that the native `losat` command line `argv` (the program first) made, in a folder of
 * `inputs` (in the browser's time zone), and its exit status and standard error. With `place`, the folder's name in the
 * made files is the engine's `/work`.
 */
export function nativeRun(
  argv: readonly string[],
  inputs: ReadonlyMap<string, Buffer>,
  place = false,
): { readonly stdout: Buffer; readonly stderr: Buffer; readonly status: number | null; readonly made: Map<string, Buffer> } {
  const folder = folderOf(inputs);
  try {
    const run = spawnSync(NATIVE!, argv, { cwd: folder, env: { ...process.env, TZ: ZONE }, maxBuffer: 1 << 28 });
    const inWork = (bytes: Buffer) => (place ? Buffer.from(bytes.toString('latin1').split(folder).join('/work'), 'latin1') : bytes);
    const made = new Map<string, Buffer>();
    for (const [name, bytes] of treeOf(folder)) if (!inputs.has(name)) made.set(name, inWork(bytes));
    return { stdout: run.stdout, stderr: run.stderr, status: run.status, made };
  } finally {
    rmSync(folder, { recursive: true, force: true });
  }
}

/** The files that the native makeblastdb made from `inputs` (in a folder of them, in the browser's zone), by name. */
export function nativeMakeblastdb(argv: readonly string[], inputs: ReadonlyMap<string, Buffer>): Map<string, Buffer> {
  const run = nativeRun(argv, inputs, true);
  if (run.status !== 0) throw new Error(`native ${argv.join(' ')} exited ${run.status}: ${run.stderr.toString()}`);
  return run.made;
}
