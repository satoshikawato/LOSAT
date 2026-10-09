// Screen records of the search screen for the visual review (S12; session README
// "画面レビュー"). Not a test: it runs only with LOSAT_WEB_SCREENS=<directory>, in the
// browsers chosen with --project, and writes full-page PNGs of the same states at a desktop
// size (1280 x 900) and a phone size (390 x 844). S13 added the states of the results screen
// (files 07 to 18, which sort after the search screen's 01 to 06); 18 is the window, not the
// whole page, right after "Open results" (S13 screen review M3).
import { mkdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, test, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { fasta, openFiles, program, submit, waitStatus } from './support/search';

const SCREENS = process.env['LOSAT_WEB_SCREENS'] || undefined;
test.skip(SCREENS === undefined, 'screen records are taken only with LOSAT_WEB_SCREENS');
test.skip(!BUILD_HAS_ENGINE, "the screens show the engine's messages: build with LOSAT_WEB_REACTORS");
test.setTimeout(600_000);

const SIZES = [
  { name: 'desktop', width: 1280, height: 900 },
  { name: 'phone', width: 390, height: 844 },
] as const;

async function shoot(page: Page, browser: string, size: string, name: string, fullPage = true): Promise<void> {
  const directory = join(SCREENS!, browser);
  mkdirSync(directory, { recursive: true });
  // Let the elapsed time and the engine's checks settle.
  await page.waitForTimeout(300);
  await page.screenshot({ path: join(directory, `${size}-${name}.png`), fullPage });
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
    // After the engine has checked the default options (S13 screen review L6).
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', /^(ok|invalid)$/, { timeout: 60_000 });
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

    // A run that fails: TBLASTX refuses at run time a subject title with an HTML character
    // reference, which NCBI decodes in outfmt 0 (docs/web/abi_v2.md §4).
    let sequence = '';
    let state = 7;
    for (let i = 0; i < 300; i++) {
      state = (Math.imul(state, 1103515245) + 12345) >>> 0;
      sequence += 'ACGT'[(state >>> 16) % 4];
    }
    await page.getByTestId('program-tblastx').check();
    for (const role of ['query', 'subject'] as const) {
      const sources = page.locator(`[data-testid^="${role}-source-"][data-status]`);
      while ((await sources.count()) > 0) await page.getByTestId(`${role}-source-0-remove`).click();
    }
    await page.getByTestId('query-input').fill(`>q1\n${sequence}\n`);
    await page.getByTestId('subject-input').fill(`>s1 alpha &amp; beta\n${sequence}\n`);
    await ready(page, 'query');
    await ready(page, 'subject');
    await page.getByTestId('add-to-queue').click();
    await expect(page.locator('li[data-testid^="run-"][data-status="failed"]')).toHaveCount(1, { timeout: 60_000 });

    // A search in progress, one cancelled before it started, one waiting, and the notice
    // after the page was hidden.
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
    const running = page.locator('li[data-testid^="run-"][data-status="running"]');
    const queued = page.locator('li[data-testid^="run-"][data-status="queued"]');
    await page.getByTestId('add-to-queue').click();
    await expect(running).toHaveCount(1, { timeout: 60_000 });
    await page.getByTestId('add-to-queue').click();
    await expect(queued).toHaveCount(1);
    await queued.getByRole('button', { name: 'Cancel' }).click();
    await page.getByTestId('param-evalue').fill('1e-10');
    await page.getByTestId('add-to-queue').click();
    await expect(queued).toHaveCount(1);
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
    await shoot(page, browserName, size.name, '06-running-queued-failed-resume');
    await queued.getByRole('button', { name: 'Cancel' }).click();
    await running.getByRole('button', { name: 'Cancel' }).first().click();
  });
}

// --- the results screen (S13) --------------------------------------------------------------

/** Removes every input source of a role (the search form keeps the sources of the last search). */
async function clearInputs(page: Page): Promise<void> {
  for (const role of ['query', 'subject'] as const) {
    const sources = page.locator(`[data-testid^="${role}-source-"][data-status]`);
    while ((await sources.count()) > 0) await page.getByTestId(`${role}-source-0-remove`).click();
  }
}

/** Searches a pair of FASTA files of LOSAT/tests/fasta with the program's defaults, and waits until it completes. */
async function search(page: Page, id: string, number: number, query: string, subject: string): Promise<void> {
  await clearInputs(page);
  await program(page, id);
  await openFiles(page, 'query', [{ name: query.split('/').pop()!, text: fasta(query) }]);
  await openFiles(page, 'subject', [{ name: subject.split('/').pop()!, text: fasta(subject) }]);
  await submit(page);
  await waitStatus(page, number, 'completed');
}

/** Opens a completed run from the queue and waits until its first query's first HSP is read. */
async function openResults(page: Page, number: number): Promise<void> {
  await page.getByTestId(`run-${number}-open`).click();
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', String(number), { timeout: 60_000 });
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-state', 'ready', { timeout: 60_000 });
}

for (const size of SIZES) {
  test(`results screen at the ${size.name} size`, async ({ page, browserName }) => {
    await page.setViewportSize({ width: size.width, height: size.height });
    await page.goto('/');
    await expect(page.getByTestId('storage-status')).toBeVisible();
    const hspRows = page.getByTestId('hsp-list').locator('[data-testid^="hsp-row-"]');

    // Run 1: three queries (one without hits, one on the minus strand) against six subjects;
    // run 2: one query with 260 subjects, of which outfmt 0 shows 250 alignments; run 3: TBLASTX.
    await search(page, 'blastn', 1, 'outfmt0/multi_query.fasta', 'outfmt0/multi_subject.fasta');
    await clearInputs(page);
    await page.getByTestId('param-task').selectOption('blastn');
    await openFiles(page, 'query', [{ name: 'many_query.fasta', text: fasta('outfmt0/many_query.fasta') }]);
    await openFiles(page, 'subject', [{ name: 'many_subject.fasta', text: fasta('outfmt0/many_subject.fasta') }]);
    await submit(page);
    await waitStatus(page, 2, 'completed');
    await search(page, 'tblastx', 3, 'outfmt0/tblastx_ambig_query.fasta', 'outfmt0/tblastx_ambig_subject.fasta');

    // The hits view: an HSP selected, with its alignment from outfmt 0.
    await openResults(page, 1);
    await hspRows.first().click();
    await expect(hspRows.first()).toHaveAttribute('aria-pressed', 'true');
    await expect(page.getByTestId('detail-section')).toBeVisible();
    await shoot(page, browserName, size.name, '07-results-hits-alignment');

    // The dot plot of a pair with an HSP on each strand (subject msD), the second HSP selected.
    await page.getByTestId('subject-list').locator('[data-testid^="subject-row-"]', { hasText: 'msD' }).click();
    await expect(hspRows).toHaveCount(2);
    await hspRows.nth(1).click();
    await page.getByTestId('pane-dotplot').click();
    const canvas = page.getByTestId('dotplot-canvas');
    await expect(canvas).toHaveAttribute('data-segments', '2');
    await expect(canvas).toHaveAttribute('data-selected', /^0:\d+$/);
    await shoot(page, browserName, size.name, '08-results-dotplot');
    await page.getByTestId('pane-alignment').click();

    // An HSP that outfmt 0 does not show: the last subject of run 2, in the engine's order.
    await openResults(page, 2);
    await page.getByTestId('subject-sort-order').click();
    await expect(page.getByTestId('subject-sort-order').locator('..')).toHaveAttribute('aria-sort', 'descending');
    await page.getByTestId('subject-list').evaluate((element) => (element.scrollTop = 0));
    const last = page.getByTestId('subject-list').locator('[data-testid^="subject-row-"]').first();
    await expect(last).toHaveAttribute('data-order', '260');
    await last.click();
    await expect(last).toHaveAttribute('aria-pressed', 'true');
    await expect(page.getByTestId('detail-not-in-outfmt0')).toBeVisible();
    await shoot(page, browserName, size.name, '09-results-hsp-not-in-outfmt0');

    // A query without hits.
    await openResults(page, 1);
    await page.getByTestId('query-row-1').click();
    await expect(page.getByTestId('query-row-1')).toHaveAttribute('aria-pressed', 'true');
    await expect(page.locator('[data-testid="results-notice"][data-kind="no-hits"]')).toBeVisible();
    await shoot(page, browserName, size.name, '10-results-query-without-hits');

    // HSPs hidden by the view filters: the notice that tells them from a query without hits.
    await page.getByTestId('query-row-0').click();
    await expect(page.getByTestId('query-row-0')).toHaveAttribute('aria-pressed', 'true');
    await page.getByTestId('filter-subject').fill('no-such-subject');
    await page.getByTestId('filter-apply').click();
    await expect(page.locator('[data-testid="results-notice"][data-kind="filtered-out"]')).toBeVisible();
    await shoot(page, browserName, size.name, '11-results-filtered-out');
    await page.getByTestId('results-notice-clear').click();
    await expect(page.getByTestId('hsp-table')).toBeVisible();

    // The run's details with the verification badge, then its outputs.
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('verification-details')).toBeVisible();
    await shoot(page, browserName, size.name, '12-results-run-details');
    await page.getByTestId('results-view-outputs').click();
    await page.getByTestId('format-0').click();
    await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', '1:0');
    await shoot(page, browserName, size.name, '13-results-outputs');

    // Frames of a translated search.
    await openResults(page, 3);
    await page.getByTestId('results-view-hits').click();
    await expect(page.locator('[data-field="frames"]').first()).toBeVisible();
    await shoot(page, browserName, size.name, '14-results-translated-frames');

    // The queue: finished runs with "Open results", a search in progress, one cancelled before
    // it started, and one waiting (queued after the cancel). BLASTP of a bacterial proteome
    // against itself keeps the engine busy long enough.
    await page.getByTestId('tab-search').click();
    await clearInputs(page);
    await program(page, 'blastp');
    for (const role of ['query', 'subject'] as const) {
      await openFiles(page, role, [{ name: 'NZ_CP006932.faa', text: fasta('NZ_CP006932.faa') }]);
    }
    const running = page.locator('li[data-testid^="run-"][data-status="running"]');
    const queued = page.locator('li[data-testid^="run-"][data-status="queued"]');
    await submit(page);
    await expect(running).toHaveCount(1, { timeout: 60_000 });
    await submit(page);
    await expect(queued).toHaveCount(1);
    await queued.getByRole('button', { name: 'Cancel' }).click();
    await expect(page.getByTestId('run-5-status')).toHaveText('cancelled');
    await submit(page);
    await expect(queued).toHaveCount(1);
    // The queue shows beside either tab; the results tab keeps the page short, so that the queue is most of it.
    await page.getByTestId('tab-results').click();
    await shoot(page, browserName, size.name, '15-results-queue-open-waiting-cancelled');

    // The results screen for the cancelled run.
    const select = page.getByTestId('results-run');
    const label = await select.locator('option', { hasText: /^Run 5 · / }).textContent();
    await select.selectOption({ label: label!.trim() });
    await expect(page.getByTestId('results-status')).toHaveAttribute('data-run-status', 'cancelled');
    await shoot(page, browserName, size.name, '16-results-cancelled-run');

    // The dot plot of a translated search (TBLASTX, run 3): axes in nt, the frames of the selected HSP.
    const run3 = await select.locator('option', { hasText: /^Run 3 · / }).textContent();
    await select.selectOption({ label: run3!.trim() });
    await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '3', { timeout: 60_000 });
    await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-state', 'ready', { timeout: 60_000 });
    await page.getByTestId('pane-dotplot').click();
    await expect(page.getByTestId('dotplot-canvas')).toHaveAttribute('data-segments', /^[1-9]/);
    await expect(page.getByTestId('dotplot-selected')).toContainText('frames');
    await shoot(page, browserName, size.name, '17-results-translated-dotplot');
    await page.getByTestId('pane-alignment').click();

    // "Open results" from the search tab: the window right after the tap shows the results' heading.
    await page.getByTestId('tab-search').click();
    await page.getByTestId('run-1-open').click();
    await expect(page.getByTestId('results-heading')).toBeFocused();
    await expect(page.getByTestId('results-heading')).toBeInViewport({ ratio: 1 });
    await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-state', 'ready', { timeout: 60_000 });
    await shoot(page, browserName, size.name, '18-results-open-from-queue-window', false);
    await queued.getByRole('button', { name: 'Cancel' }).click();
    await running.getByRole('button', { name: 'Cancel' }).first().click();
  });
}
