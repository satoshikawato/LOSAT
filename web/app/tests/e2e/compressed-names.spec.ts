// Line WQA's D21 (2026-11-23, S3-3 of the third sweep and F1/F2 of the review of D20): a run of a compressed query or subject
// names the file as the user gave it where the command line writes a file name: the command shown, the run picker, outfmt 0 and 7's
// `Input:` line of a subject set, and SAM's command line. Run on the shown command in a folder with the user's files, the native CLI
// writes the bytes that the application stored, for outfmt 0, 6, 7 and 17 (the engine build, LOSAT_WEB_REACTORS and LOSAT_WEB_NATIVE).
import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { gzipSync, zstdCompressSync } from 'node:zlib';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { chooseRun, fasta, openFiles, openParameters, openSearch, pickRun, program, storedOutputs, submit, waitStatus } from './support/search';

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

/** The native CLI's standard output for the command that the application showed, run from a folder with the user's files. */
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

/** The command that Outputs shows for run `run`'s format `format` (a number, or 'custom' for the run's own -outfmt). */
async function shownCommand(page: Page, run: number, format: 6 | 'custom'): Promise<string> {
  await chooseRun(page, run);
  await page.getByTestId('results-view-outputs').click();
  await page.getByTestId(`format-${format}`).click();
  await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', `${run}:${format === 'custom' ? 4 : format}`);
  return ((await page.getByTestId('result-command').textContent()) ?? '').trim();
}

test('a compressed query or subject is named as the user gave it: the command, the run picker, outfmt 0, 6, 7 and 17 equal the CLI run on it (S3-3, F1, F2; D21)', async ({ page }) => {
  const query = fasta('LC738884.fasta');
  const bz2 = spawnSync('bzip2', ['-9c'], { input: query, maxBuffer: 1 << 26 });
  expect(bz2.status, String(bz2.error ?? bz2.stderr)).toBe(0);
  const CASES: ReadonlyArray<{ readonly query: UserFile; readonly subject: UserFile }> = [
    { query: { name: 'q.fa.gz', bytes: gzipSync(query) }, subject: { name: 's.fa.zst', bytes: zstdCompressSync(query) } },
    { query: { name: 'q.fa.bz2', bytes: bz2.stdout }, subject: { name: 's.fa', bytes: query } },
    { query: { name: 'Q.FA.ZST', bytes: zstdCompressSync(query) }, subject: { name: 's.fa.bz2', bytes: bz2.stdout } },
  ];
  await program(page, 'blastn');
  let run = 0;
  for (const c of CASES) {
    // The previous case's files go, this case's come.
    if (run > 0) {
      await page.getByTestId('tab-search').click();
      await page.getByTestId('query-source-0-remove').click();
      await page.getByTestId('subject-source-0-remove').click();
    }
    await openFiles(page, 'query', [{ name: c.query.name, text: c.query.bytes }]);
    await openFiles(page, 'subject', [{ name: c.subject.name, text: c.subject.bytes }]);
    const files = [c.query, c.subject];
    const shown = `losat blastn -query ${c.query.name} -subject ${c.subject.name}`;
    // A run of the form's defaults: outfmt 0, 6 and 7.
    await submit(page);
    run += 1;
    await waitStatus(page, run, 'completed');
    // The run picker, Outputs' command line, Run details' argument and command lines all name the user's files.
    await expect(page.getByTestId(`run-${run}-option`)).toHaveText(new RegExp(`^Run ${run} · BLASTN · ${c.query.name} vs ${c.subject.name}$`));
    expect(await shownCommand(page, run, 6)).toBe(`${shown} -outfmt 6`);
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-argv')).toHaveText(shown.replace('losat ', ''));
    for (const format of [0, 6, 7]) await expect(page.getByTestId(`run-command-${format}`)).toHaveText(`${shown} -outfmt ${format}`);
    if (NATIVE === undefined) test.info().annotations.push({ type: 'skip', description: `The CLI's bytes: ${NO_NATIVE_REASON}` });
    else {
      const stored = (await storedOutputs(page, run)).outputs;
      for (const format of [0, 6, 7] as const) {
        expect(stored[format].length).toBeGreaterThan(0);
        expect(nativeShown(`${shown} -outfmt ${format}`, files).equals(stored[format]), `${shown} -outfmt ${format}: the native CLI's standard output`).toBe(true);
      }
    }
    // A run of its own -outfmt 17 (SAM, whose command line names both files).
    await page.getByTestId('tab-search').click();
    await openParameters(page);
    await page.getByTestId('param-extra').fill('-outfmt 17');
    await submit(page);
    run += 1;
    await waitStatus(page, run, 'completed');
    expect(await shownCommand(page, run, 'custom')).toBe(`${shown} -outfmt 17`);
    if (NATIVE !== undefined) {
      const text = (await page.getByTestId('result-output').textContent()) ?? '';
      expect(text).toContain('@PG');
      expect(text, `${shown} -outfmt 17: the native CLI's standard output`).toBe(nativeShown(`${shown} -outfmt 17`, files).toString('utf8'));
    }
    await page.getByTestId('tab-search').click();
    await page.getByTestId('param-extra').fill('');
  }
  await pickRun(page, 1);
});
