// A tool form without its input (line NP WW4f, review R0's F1, 2026-11-28). Since WW4f the engine's check of a masking
// tool's argv opens `-in` as the run does (and reads its start); the form's check names an input that the form does not
// have yet by a stand-in, an empty `input.fa` read as FASTA. So the Masking tab (dustmasker, segmasker,
// convert2blastmask) and WindowMasker's Mask stage show no refusal when they open, with a FASTA input or a BLAST database
// not named yet (`-infmt blastdb`), as before WW4f; Run LOSAT asks for the input; and the check still reports the form's
// other errors: with the real engine, an `-out` in a folder that the run does not have, which segmasker,
// convert2blastmask and WindowMasker's masking open after `-in`.
//
// LOSAT_E2E_SCREENSHOT_DIR=<dir> saves the Masking tab of each tool and WindowMasker's Mask stage once their check
// answered (<browser>-empty-<tool>.png).
import { mkdirSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { described, openSearch } from './support/search';

const SCREENSHOT_DIR = process.env['LOSAT_E2E_SCREENSHOT_DIR'];
const FASTA_ASKED = 'Enter sequences or open FASTA files.';
const NAME_ASKED = 'Enter an input name or open files.';
const OUT = 'nodir/m.out';

test.setTimeout(120_000);

async function shot(page: Page, name: string): Promise<void> {
  if (!SCREENSHOT_DIR) return;
  mkdirSync(SCREENSHOT_DIR, { recursive: true });
  await page.screenshot({ path: join(SCREENSHOT_DIR, `${test.info().project.name}-empty-${name}.png`), fullPage: true, animations: 'disabled' });
}

/** The check of the tab (`ms`, `wm`) has answered (`ok` or `invalid`). */
async function answered(page: Page, prefix: 'ms' | 'wm'): Promise<void> {
  await expect(page.getByTestId(`${prefix}-argv-validation`)).toHaveAttribute('data-state', /^(ok|invalid)$/, { timeout: 30_000 });
}

/** The check of the tab (`ms`, `wm`) has answered `state`; for `ok` no refusal is shown (else the failure names it). */
async function checkIs(page: Page, prefix: 'ms' | 'wm', state: 'ok' | 'invalid'): Promise<void> {
  const validation = page.getByTestId(`${prefix}-argv-validation`);
  await answered(page, prefix);
  if (state === 'ok') await expect(validation).toHaveText('');
  await expect(validation).toHaveAttribute('data-state', state);
}

/** Run LOSAT says `message` and queues no run. */
async function runAsks(page: Page, prefix: 'ms' | 'wm', message: string): Promise<void> {
  await page.getByTestId(`${prefix}-run-losat`).click();
  await expect(page.getByTestId(`${prefix}-message`)).toHaveText(message);
  await expect(page.getByTestId('run-1-option')).toHaveCount(0);
}

/** With the real engine, an `-out` in a missing folder is the check's refusal, which names it; the form is then clean again. */
async function outRefused(page: Page, prefix: 'ms' | 'wm'): Promise<void> {
  if (!BUILD_HAS_ENGINE) return;
  const out = page.getByTestId(`${prefix}-out`);
  const before = await out.inputValue();
  await out.fill(OUT);
  await checkIs(page, prefix, 'invalid');
  await expect(page.getByTestId(`${prefix}-argv-validation`)).toContainText(OUT);
  await out.fill(before);
  await checkIs(page, prefix, 'ok');
}

const MASKING = [
  { tool: 'dustmasker', database: true },
  { tool: 'segmasker', database: true },
  { tool: 'convert2blastmask', database: false },
] as const;

for (const { tool, database } of MASKING) {
  test(`${tool}: the empty Masking form shows no refusal, Run LOSAT asks for the input, and the check reports -out`, async ({ page }) => {
    await openSearch(page);
    await page.getByTestId('tab-masking').click();
    await described(page, 'ms-panel', 'ready');
    if (tool !== 'dustmasker') {
      await page.getByTestId(`ms-tool-${tool}`).check();
      await described(page, 'ms-panel', 'ready');
    }
    if (tool === 'convert2blastmask') {
      // Its mandatory options are the check's refusal until typed (WW4c), as the command line's.
      await checkIs(page, 'ms', 'invalid');
      await expect(page.getByTestId('ms-argv-validation')).toContainText('masking_algorithm');
      await page.getByTestId('ms-param-masking_algorithm').fill('dust');
      await page.getByTestId('ms-param-masking_options').fill('window=64');
    }
    await answered(page, 'ms');
    await shot(page, tool);
    await checkIs(page, 'ms', 'ok');
    await runAsks(page, 'ms', FASTA_ASKED);
    await outRefused(page, 'ms');
    if (!database) return;

    // A BLAST database not named yet: no refusal either, and Run LOSAT asks for its name.
    // dustmasker's -infmt is a text field, segmasker's a list (their -xmlhelp).
    const infmt = page.getByTestId('ms-param-infmt');
    if ((await infmt.evaluate((element) => element.tagName)) === 'SELECT') await infmt.selectOption('blastdb');
    else await infmt.fill('blastdb');
    await expect(page.getByTestId('ms-in-name')).toBeVisible();
    await checkIs(page, 'ms', 'ok');
    await runAsks(page, 'ms', NAME_ASKED);
    await outRefused(page, 'ms');
  });
}

test("WindowMasker's Mask stage: the empty form shows no refusal, Run LOSAT asks for the counts and the input, and the check reports -out", async ({
  page,
}) => {
  await openSearch(page);
  await page.getByTestId('tab-windowmasker').click();
  await described(page, 'wm-panel', 'ready');
  await page.getByTestId('wm-stage-mask').check();
  await answered(page, 'wm');
  await shot(page, 'windowmasker-mask');
  await checkIs(page, 'wm', 'ok');
  await runAsks(page, 'wm', 'Choose a counts file for -ustat.');
  await page.getByTestId('wm-ustat-file').setInputFiles({ name: 'c.counts', mimeType: 'text/plain', buffer: Buffer.from('') });
  await expect(page.getByTestId('wm-ustat-chosen')).toContainText('c.counts');
  await checkIs(page, 'wm', 'ok');
  await runAsks(page, 'wm', FASTA_ASKED);
  await outRefused(page, 'wm');

  // A BLAST database not named yet.
  await page.getByTestId('wm-param-infmt').selectOption('blastdb');
  await expect(page.getByTestId('wm-in-name')).toBeVisible();
  await checkIs(page, 'wm', 'ok');
  await runAsks(page, 'wm', NAME_ASKED);
  await outRefused(page, 'wm');
});
