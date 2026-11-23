// Line WQA's D14 (2026-11-20): the words with which the command line prints its help (`-h`, `-help`,
// `-help-full`, `-xmlhelp`) or its version (`-version` and its full forms), typed in "Additional
// options" or as a value that the command line reads as one, run as the native CLI runs them: the run
// needs no input, prints the command line's text and searches nothing. The results say so, the Outputs
// show the text under "Help" or "Version" and save it. With LOSAT_WEB_NATIVE, the saved text is the
// bytes that the native CLI of the same commit writes on its standard output for the command shown
// (web/adapter/tests/v_abi.js compares every program and word through the reactors).
import { spawnSync } from 'node:child_process';
import { mkdirSync, readFileSync } from 'node:fs';
import { expect, type Download, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { chooseRun, openParameters, openSearch, program, submit, waitStatus } from './support/search';

test.setTimeout(180_000);
test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine prints no help or version: the engine build runs these tests');

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

/** The bytes of the file that `click` downloads. */
async function downloaded(page: Page, click: () => Promise<void>): Promise<{ readonly name: string; readonly bytes: Buffer }> {
  const [download]: [Download, void] = await Promise.all([page.waitForEvent('download'), click()]);
  return { name: download.suggestedFilename(), bytes: readFileSync((await download.path())!) };
}

/**
 * What the native CLI writes on its standard output for the command shown (`losat` first, words without
 * spaces), which ends with exit status 0 and writes nothing on its standard error: undefined without
 * LOSAT_WEB_NATIVE.
 */
function nativeText(command: string): Buffer | undefined {
  if (NATIVE === undefined) return undefined;
  const directory = test.info().outputPath();
  mkdirSync(directory, { recursive: true });
  const run = spawnSync(NATIVE, command.split(' ').slice(1), { cwd: directory, maxBuffer: 1 << 26 });
  expect(run.status, command).toBe(0);
  expect(run.stderr.toString('utf8'), command).toBe('');
  return run.stdout;
}

const CASES = [
  { id: 'blastn', extra: '-help', kind: 'help', label: 'Help', starts: 'USAGE\n', ext: 'txt' },
  { id: 'blastp', extra: '-xmlhelp', kind: 'help', label: 'Help', starts: '<?xml version="1.0"?>', ext: 'xml' },
  { id: 'blastx', extra: '-version-full-json', kind: 'version', label: 'Version', starts: '{\n', ext: 'json' },
  // A toolkit word as the value of a free-text field (the Low complexity regions filter's arguments): the
  // command line prints its version for `-version` anywhere in its words.
  { id: 'tblastx', seg: '-version', kind: 'version', label: 'Version', starts: 'tblastx: ', ext: 'txt' },
] as const;

test("help and version words print the command line's text without inputs, and search nothing", async ({ page }, testInfo) => {
  if (NATIVE === undefined) testInfo.annotations.push({ type: 'skip', description: `The bytes of the text: ${NO_NATIVE_REASON}` });
  for (const [index, c] of CASES.entries()) {
    const run = index + 1;
    const notice = `This run printed the command line's ${c.kind} and ran no search. The text is in Outputs.`;
    if (run > 1) await page.getByTestId('tab-search').click();
    await program(page, c.id);
    await openParameters(page);
    await page.getByTestId('param-extra').fill('extra' in c ? c.extra : '');
    if ('seg' in c) await page.getByTestId('param-seg-value').fill(c.seg);
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
    // No query and no subject: the command line reads none for these words.
    await submit(page);
    await waitStatus(page, run, 'completed');
    // The run picker names the mode alone: the run read no file (line WQA's D21, the review of D20's F2, 2026-11-23).
    await expect(page.getByTestId(`run-${run}-option`)).toHaveText(new RegExp(`^Run ${run} · [A-Z]+$`));
    await chooseRun(page, run);
    await expect(page.getByTestId('results-notice')).toHaveAttribute('data-kind', `printed-${c.kind}`);
    await expect(page.getByTestId('results-notice')).toHaveText(notice);
    await expect(page.getByTestId('section-graphic')).toHaveCount(0);
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-printed')).toHaveText(`the command line's ${c.kind}: no search`);

    // The Outputs show the text in place of the formats, with the command that prints it.
    await page.getByTestId('results-view-outputs').click();
    const tab = page.getByTestId('format-printed');
    await expect(tab).toHaveText(c.label);
    await expect(tab).toHaveAttribute('aria-selected', 'true');
    await expect(page.getByTestId('format-6')).toHaveCount(0);
    await expect(page.getByTestId('printed-note')).toHaveText('The command above prints this text and runs no search.');
    const command = (await page.getByTestId('result-command').textContent()) ?? '';
    expect(command.startsWith(`losat ${c.id} `), command).toBe(true);
    const saved = await downloaded(page, () => page.getByTestId('export-printed').click());
    expect(saved.name).toMatch(new RegExp(`^losat-run${run}-${c.id}(-[a-z-]+)?-\\d{8}T\\d{6}Z\\.${c.kind}\\.${c.ext}$`));
    expect(saved.bytes.toString('utf8').startsWith(c.starts), saved.bytes.toString('utf8').slice(0, 80)).toBe(true);
    await expect(page.getByTestId('printed-output')).toContainText(c.starts.trim());
    // "Download All" lists the text alone, and saves the same bytes under the same name.
    await page.getByTestId('results-download-all').click();
    const menu = page.getByTestId('results-download-menu');
    await expect(menu.getByTestId('results-download-printed')).toHaveText(c.label);
    await expect(menu.getByTestId('results-download-6')).toHaveCount(0);
    const again = await downloaded(page, () => menu.getByTestId('results-download-printed').click());
    expect(again).toEqual(saved);
    const native = nativeText(command);
    if (native !== undefined) expect(saved.bytes.equals(native), `${command}: the native CLI's standard output`).toBe(true);
  }
});
