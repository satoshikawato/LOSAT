// Screen records of the search screen for the visual review (S12; session README
// "画面レビュー"). Not a test: it runs only with LOSAT_WEB_SCREENS=<directory>, in the
// browsers chosen with --project, and writes full-page PNGs of the same states at a desktop
// size (1280 x 900) and a phone size (390 x 844).
import { mkdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, test, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';

const SCREENS = process.env['LOSAT_WEB_SCREENS'] || undefined;
test.skip(SCREENS === undefined, 'screen records are taken only with LOSAT_WEB_SCREENS');
test.skip(!BUILD_HAS_ENGINE, "the screens show the engine's messages: build with LOSAT_WEB_REACTORS");
test.setTimeout(600_000);

const SIZES = [
  { name: 'desktop', width: 1280, height: 900 },
  { name: 'phone', width: 390, height: 844 },
] as const;

async function shoot(page: Page, browser: string, size: string, name: string): Promise<void> {
  const directory = join(SCREENS!, browser);
  mkdirSync(directory, { recursive: true });
  // Let the elapsed time and the engine's checks settle.
  await page.waitForTimeout(300);
  await page.screenshot({ path: join(directory, `${size}-${name}.png`), fullPage: true });
}

async function ready(page: Page, role: string, index = 0): Promise<void> {
  await expect(page.getByTestId(`${role}-source-${index}`)).toHaveAttribute('data-status', /ready|failed/, { timeout: 60_000 });
  const check = page.getByTestId(`${role}-source-${index}-check`);
  if ((await check.count()) > 0) await expect(check).not.toHaveAttribute('data-check', 'pending', { timeout: 60_000 });
}

for (const size of SIZES) {
  test(`search screen at the ${size.name} size`, async ({ page, browserName }) => {
    await page.setViewportSize({ width: size.width, height: size.height });
    await page.goto('/');
    await expect(page.getByTestId('storage-status')).toBeVisible();
    await shoot(page, browserName, size.name, '01-empty');

    // Inputs: a query with a duplicate ID, a refused record and a protein-like record; one subject record with a region.
    await page.getByTestId('query-input').fill(
      '>contig_1 assembled contig\nACGTACGTTTGACCATGGCATGCATGCATTTAGGCCAAGTACGATCGATCG\n' +
        '>contig_1 second copy\nACGTTGCAACGTTGCAACGTTGCAAGGT\n' +
        '>gene_x\nMKLVVLAAGGHHKLMKLVVLAAGG\n',
    );
    await page.getByTestId('subject-files').setInputFiles({
      name: 'LC738884.fasta',
      mimeType: 'text/plain',
      buffer: readFileSync(join(REPOSITORY, 'LOSAT/tests/fasta/LC738884.fasta')),
    });
    await ready(page, 'query');
    await ready(page, 'subject');
    await page.getByTestId('subject-region-start').fill('1001');
    await page.getByTestId('subject-region-stop').fill('25000');
    await page.getByTestId('param-evalue').fill('1e-5');
    await page.getByTestId('param-task').selectOption('blastn');
    await expect(page.getByTestId('argv-validation')).not.toHaveAttribute('data-state', 'checking');
    await shoot(page, browserName, size.name, '02-inputs-refused-region');

    // Exclude the refused record and queue a group of separate searches behind a running one.
    await page.getByTestId('query-source-0-exclude-refused').click();
    await ready(page, 'query');
    await page.getByTestId('add-to-queue').click();
    await expect(page.getByTestId('run-1')).toBeVisible();
    await page.getByTestId('subject-source-0-remove').click();
    await page.getByTestId('subject-files').setInputFiles([
      { name: 'a.fa', mimeType: 'text/plain', buffer: Buffer.from('>a\nACGTACGTTTGACCATGGCATGCATGCATTTAGG\n') },
      { name: 'b.fa', mimeType: 'text/plain', buffer: Buffer.from('>b\nGGCATGCATGCATTTAGGCCAAGTACGATCGATCG\n') },
    ]);
    await ready(page, 'subject', 0);
    await ready(page, 'subject', 1);
    await page.getByTestId('subject-mode-separate').check();
    await page.getByTestId('add-to-queue').click();
    await expect(page.getByTestId('run-3')).toBeVisible();
    await expect(page.getByTestId('run-1-status')).toHaveText(/completed|failed/, { timeout: 120_000 });
    await expect(page.getByTestId('run-3-status')).toHaveText(/completed|failed/, { timeout: 120_000 });
    await page.getByTestId('run-1').locator('summary', { hasText: 'Details' }).click();
    await shoot(page, browserName, size.name, '03-separate-queue-details');

    // Parameters of the other programs.
    for (const program of ['blastp', 'tblastn', 'tblastx'] as const) {
      await page.getByTestId(`program-${program}`).check();
      await expect(page.getByTestId('parameter-form')).toBeVisible();
      if (program === 'tblastx') await page.getByTestId('param-db_gencode').selectOption('11');
      await expect(page.getByTestId('argv-validation')).not.toHaveAttribute('data-state', 'checking');
      await shoot(page, browserName, size.name, `04-${program}`);
    }
    await page.getByTestId('program-blastx').check();
    await shoot(page, browserName, size.name, '05-blastx-unavailable');

    if (BUILD_HAS_ENGINE) {
      // A search in progress, and the notice after the page was hidden.
      await page.getByTestId('program-blastp').check();
      for (const role of ['query', 'subject'] as const) {
        const sources = page.locator(`[data-testid^="${role}-source-"][data-status]`);
        while ((await sources.count()) > 0) await page.getByTestId(`${role}-source-0-remove`).click();
        await page.getByTestId(`${role}-files`).setInputFiles({
          name: 'NZ_CP006932.faa',
          mimeType: 'text/plain',
          buffer: readFileSync(join(REPOSITORY, 'LOSAT/tests/fasta/NZ_CP006932.faa')),
        });
        await ready(page, role);
      }
      await page.getByTestId('add-to-queue').click();
      const run = page.locator('[data-testid^="run-"][data-status="running"]');
      await expect(run).toHaveCount(1, { timeout: 60_000 });
      await page.evaluate(() => {
        Object.defineProperty(document, 'visibilityState', { get: () => 'hidden', configurable: true });
        document.dispatchEvent(new Event('visibilitychange'));
      });
      await page.waitForTimeout(2100);
      await page.evaluate(() => {
        Object.defineProperty(document, 'visibilityState', { get: () => 'visible', configurable: true });
        document.dispatchEvent(new Event('visibilitychange'));
      });
      await expect(page.getByTestId('resume-data')).not.toHaveText('Checking the stored results…');
      await shoot(page, browserName, size.name, '06-running-resume');
      await run.getByRole('button', { name: 'Cancel' }).first().click();
    }
  });
}
