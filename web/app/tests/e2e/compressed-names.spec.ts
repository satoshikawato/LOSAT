// Line WQA's D21 (2026-11-23, split from D20's S3-3): the command shown for a compressed query or subject names the file as the
// user gave it, and run from a folder with the user's files the native CLI writes the bytes that the application stored
// (the engine build, LOSAT_WEB_REACTORS, and LOSAT_WEB_NATIVE).
import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { gzipSync, zstdCompressSync } from 'node:zlib';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { chooseRun, fasta, openFiles, openSearch, program, storedOutputs, submit, waitStatus } from './support/search';

test.setTimeout(240_000);
test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine inflates no file: the engine build runs these tests');

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

/** An input of a run: its bytes under the name that the user's file has. */
interface UserFile {
  readonly name: string;
  readonly bytes: Buffer;
}

/** The native CLI's outfmt 6 for the command that the application showed, run from a folder with the user's files. */
function nativeShown(command: string, files: readonly UserFile[]): Buffer {
  mkdirSync(test.info().outputPath(), { recursive: true });
  const directory = mkdtempSync(join(test.info().outputPath(), 'user-'));
  for (const file of files) writeFileSync(join(directory, file.name), file.bytes);
  const words = command.split(' ');
  expect(words[0], command).toBe('losat');
  const run = spawnSync(NATIVE!, [...words.slice(1), '-num_threads', '1'], { cwd: directory, maxBuffer: 1 << 28 });
  expect(run.status, `${command}: ${run.stderr.toString('utf8')}`).toBe(0);
  return run.stdout;
}

async function nameOfCommand(page: Page, run: number): Promise<string> {
  await chooseRun(page, run);
  await page.getByTestId('results-view-outputs').click();
  await page.getByTestId('format-6').click();
  await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', `${run}:6`);
  return ((await page.getByTestId('result-command').textContent()) ?? '').trim();
}

test('the command shown for a compressed query or subject names the file as the user gave it, and runs on the CLI (S3-3, D21)', async ({ page }) => {
  const query = fasta('LC738884.fasta');
  const bz2 = spawnSync('bzip2', ['-9c'], { input: query, maxBuffer: 1 << 26 });
  expect(bz2.status, String(bz2.error ?? bz2.stderr)).toBe(0);
  const CASES: ReadonlyArray<{ readonly query: UserFile; readonly subject: UserFile }> = [
    { query: { name: 'q.fa.gz', bytes: gzipSync(query) }, subject: { name: 's.fa.zst', bytes: zstdCompressSync(query) } },
    { query: { name: 'q.fa.bz2', bytes: bz2.stdout }, subject: { name: 's.fa', bytes: query } },
  ];
  await program(page, 'blastn');
  for (const [index, c] of CASES.entries()) {
    const run = index + 1;
    if (run > 1) await page.getByTestId('tab-search').click();
    // The previous run's files go, the next run's come.
    if (run > 1) {
      await page.getByTestId('query-source-0-remove').click();
      await page.getByTestId('subject-source-0-remove').click();
    }
    await openFiles(page, 'query', [{ name: c.query.name, text: c.query.bytes }]);
    await openFiles(page, 'subject', [{ name: c.subject.name, text: c.subject.bytes }]);
    await submit(page);
    await waitStatus(page, run, 'completed');
    const expected = `losat blastn -query ${c.query.name} -subject ${c.subject.name} -outfmt 6`;
    // Outputs' command line, Run details' argument and command lines, and the printed argv of the run all name the user's files.
    expect(await nameOfCommand(page, run)).toBe(expected);
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-argv')).toHaveText(`blastn -query ${c.query.name} -subject ${c.subject.name}`);
    await expect(page.getByTestId('run-command-6')).toHaveText(expected);
    await expect(page.getByTestId('run-command-0')).toHaveText(expected.replace('-outfmt 6', '-outfmt 0'));
    await expect(page.getByTestId('run-command-7')).toHaveText(expected.replace('-outfmt 6', '-outfmt 7'));
    if (NATIVE === undefined) test.info().annotations.push({ type: 'skip', description: `The CLI's bytes: ${NO_NATIVE_REASON}` });
    else {
      const stored = (await storedOutputs(page, run)).outputs[6];
      expect(stored.length).toBeGreaterThan(0);
      expect(nativeShown(expected, [c.query, c.subject]).equals(stored), `${expected}: the native CLI's standard output`).toBe(true);
    }
  }
});
