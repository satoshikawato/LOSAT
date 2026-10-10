// Screen records of the search screen for the visual review (S12; session README
// "画面レビュー"). Not a test: it runs only with LOSAT_WEB_SCREENS=<directory>, in the
// browsers chosen with --project, and writes full-page PNGs of the same states at a desktop
// size (1280 x 900) and a phone size (390 x 844). S13 added the states of the results screen
// (files 07 to 18, which sort after the search screen's 01 to 06); 18 is the window, not the
// whole page, right after "Open results" (S13 screen review M3). W4b added 19 to 23: the
// Descriptions, the Graphic Summary, the Alignments with two Ranges, the dot plot's popup and a
// TBLASTN dot plot. S14 added 24 to 28: the Descriptions' marks, the Alignments with a Range "In
// candidates", the popup of a TBLASTN HSP (a W4b state not recorded before), the candidate tray
// with candidates of two runs, notes and Origins, and the tray after an extraction cut at a
// record's end. S15 (W6) added 29 to 40 in a test of their own, after which the states 01 to 28
// change only where W6 changed the screen (the tray's Subject ID on phones, the notice on a source
// read again, the reproduction panel in Run details, the Outputs tab's "LOSAT Web formats", the
// dot plot's "Download SVG"): the Outputs tab after a CSV export, Run details with the
// reproduction panel of a run with a subject record left out and of a TBLASTN run with genetic
// code 4, the search form after "Load settings..." with an item not applied, the session panel
// before saving, the queue and the results header after opening a session file, Run details of a
// loaded run without its originals, the tray's extract form that needs the original FASTA, a
// refused session file, Run details after a refused and an accepted re-attachment, the search form
// after Edit Search, and two W5 states not recorded before: the empty tray and a flank error.
//
// Since the Owner's instruction of 2026-10-10 the results' first tab is NCBI's classic one page
// (the Graphic Summary, the Descriptions directly under it, the selected subject's Alignments), so
// the states 07, 09 to 11, 14, 18 to 21, 24 and 25 show that page whole (their names say what the
// state is about); the states of the other tabs (08, 12, 13, 15 to 17, 22, 23, 26) are as before.
import { mkdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, test, type Locator, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { fasta, openFiles, openParameters, program, showRecords, submit, task, waitStatus } from './support/search';

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

    // Inputs: a query with a duplicate ID and a protein-like record, and a query file that the
    // index scan refuses (a gap line: since S14 the scan reads as the engine's reader does and
    // refuses what the engine refuses, before any check); one subject record with a region.
    await page.getByTestId('query-input').fill(
      '>contig_1 assembled contig\nACGTACGTTTGACCATGGCATGCATGCATTTAGGCCAAGTACGATCGATCG\n' +
        '>contig_1 second copy\nACGTTGCAACGTTGCAACGTTGCAAGGT\n' +
        '>gene_x\nMKLVVLAAGGHHKLMKLVVLAAGG\n',
    );
    await page.getByTestId('query-files').setInputFiles({
      name: 'contig_2.fa',
      mimeType: 'text/plain',
      buffer: Buffer.from('>contig_2\nACGTACGTTTGACCATGGCA\n>?10\nTGCATTTAGG\n'),
    });
    await page.getByTestId('subject-files').setInputFiles({
      name: 'LC738884.fasta',
      mimeType: 'text/plain',
      buffer: readFileSync(join(REPOSITORY, 'LOSAT/tests/fasta/LC738884.fasta')),
    });
    await ready(page, 'query');
    await ready(page, 'query', 1);
    await ready(page, 'subject');
    await page.getByTestId('subject-region-start').fill('1001');
    await page.getByTestId('subject-region-stop').fill('25000');
    // Algorithm parameters open from here on (W4b), with a changed value.
    await openParameters(page);
    await page.getByTestId('param-evalue').fill('1e-5');
    await task(page, 'blastn');
    await page.getByTestId('job-title').fill('Contigs against LC738884');
    await expect(page.getByTestId('argv-validation')).not.toHaveAttribute('data-state', 'checking');
    // The engine's reader now refuses at the file (above), not at a record of an indexed input, so
    // the refused-record mark and "Exclude record" are covered by the FakeEngine E2E only
    // (search.spec.ts). The user's own exclusion stays on screen: the protein-like record, struck
    // through in the list of records, left out of the following searches until a program of the
    // other reader kind reads the source again (which clears it and says so on the source).
    await showRecords(page, 'query');
    await page.getByTestId('query-source-0-record-2').uncheck();
    await expect(page.getByTestId('query-source-0-summary')).toContainText('(2 included)');
    await shoot(page, browserName, size.name, '02-inputs-refused-region');

    // Remove the file that cannot be read and queue a group of separate searches behind a running one.
    await page.getByTestId('query-source-1-remove').click();
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

    // A run that fails: a TBLASTX subject record without residues. `register` and `validate`
    // accept it, and the run ends with NCBI's `BLAST engine error: The average subject length
    // is too short`, as the CLI does.
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
    await page.getByTestId('subject-input').fill('>s1 alpha\n');
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

/**
 * Presses a table's row with the mouse at its part in view, as a finger would (W4b screen review
 * L5): Playwright's click() scrolls the whole row's button into view first, which scrolled a
 * phone's table sideways in the records, though the app keeps it where it is (results.spec.ts).
 */
async function pressRow(page: Page, row: Locator): Promise<void> {
  await row.scrollIntoViewIfNeeded();
  await row.evaluate((element) => {
    const scroll = element.closest<HTMLElement>('.table-scroll');
    if (scroll !== null) scroll.scrollLeft = 0;
  });
  const box = (await row.boundingBox())!;
  await page.mouse.click(box.x + 24, box.y + box.height / 2);
}

/** Opens a completed run from the queue and waits until its first query's first HSP is read. */
async function openResults(page: Page, number: number): Promise<void> {
  await page.getByTestId(`run-${number}-open`).click();
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', String(number), { timeout: 60_000 });
  await page.getByTestId('results-view-hits').click();
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-state', 'ready', { timeout: 60_000 });
}

for (const size of SIZES) {
  test(`results screen at the ${size.name} size`, async ({ page, browserName }) => {
    await page.setViewportSize({ width: size.width, height: size.height });
    await page.goto('/');
    await expect(page.getByTestId('storage-status')).toBeVisible();
    const hspRows = page.getByTestId('hsp-list').locator('[data-testid^="hsp-row-"]');
    const subjectRow = (text: string) => page.getByTestId('subject-list').locator('[data-testid^="subject-row-"]', { hasText: text });

    // Run 1: three queries (one without hits, one on the minus strand) against six subjects;
    // run 2: one query with 260 subjects, of which outfmt 0 shows 250 alignments; run 3: TBLASTX;
    // run 4: TBLASTN (W4b: amino acids against nucleotides on the dot plot).
    await search(page, 'blastn', 1, 'outfmt0/multi_query.fasta', 'outfmt0/multi_subject.fasta');
    await clearInputs(page);
    await task(page, 'blastn');
    await openFiles(page, 'query', [{ name: 'many_query.fasta', text: fasta('outfmt0/many_query.fasta') }]);
    await openFiles(page, 'subject', [{ name: 'many_subject.fasta', text: fasta('outfmt0/many_subject.fasta') }]);
    await page.getByTestId('job-title').fill('Many subjects');
    await submit(page);
    await waitStatus(page, 2, 'completed');
    await page.getByTestId('job-title').fill('');
    await search(page, 'tblastx', 3, 'outfmt0/tblastx_ambig_query.fasta', 'outfmt0/tblastx_ambig_subject.fasta');
    await search(page, 'tblastn', 4, 'outfmt0/e2e_protein_query.faa', 'outfmt0/e2e_amb_subject.fna');

    // The Alignments: the first subject's block, an HSP selected with its section of outfmt 0.
    await openResults(page, 1);
    await pressRow(page, hspRows.first());
    await expect(hspRows.first()).toHaveAttribute('aria-pressed', 'true');
    await expect(page.getByTestId('detail-section')).toBeVisible();
    await shoot(page, browserName, size.name, '07-results-hits-alignment');

    // The dot plot of a pair with an HSP on each strand (subject msD), the second HSP selected.
    await page.getByTestId('results-view-hits').click();
    await subjectRow('msD').click();
    await page.getByTestId('pane-dotplot').click();
    await expect(hspRows).toHaveCount(2);
    await pressRow(page, hspRows.nth(1));
    const canvas = page.getByTestId('dotplot-canvas');
    await expect(canvas).toHaveAttribute('data-segments', '2');
    await expect(canvas).toHaveAttribute('data-selected', /^0:\d+$/);
    await shoot(page, browserName, size.name, '08-results-dotplot');

    // An HSP that outfmt 0 does not show: the last subject of run 2, in the engine's order.
    await openResults(page, 2);
    await page.getByTestId('results-view-hits').click();
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
    await page.getByTestId('results-view-hits').click();
    await page.getByTestId('filter-subject').fill('no-such-subject');
    await page.getByTestId('filter-apply').click();
    await expect(page.locator('[data-testid="results-notice"][data-kind="filtered-out"]')).toBeVisible();
    await shoot(page, browserName, size.name, '11-results-filtered-out');
    await page.getByTestId('results-notice-clear').click();
    await expect(page.getByTestId('subject-table')).toBeVisible();

    // The run's details with the verification badge, then its outputs.
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('verification-details')).toBeVisible();
    await shoot(page, browserName, size.name, '12-results-run-details');
    await page.getByTestId('results-view-outputs').click();
    await page.getByTestId('format-0').click();
    await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', '1:0');
    await shoot(page, browserName, size.name, '13-results-outputs');

    // Frames of a translated search, in the HSP table of the Alignments.
    await openResults(page, 3);
    await expect(page.locator('[data-field="frames"]').first()).toBeVisible();
    await shoot(page, browserName, size.name, '14-results-translated-frames');

    // W4b states. The Descriptions of run 1, under the header block, "Filter Results" and "Results for".
    await openResults(page, 1);
    await page.getByTestId('results-view-hits').click();
    await expect(page.getByTestId('subject-list')).toBeVisible();
    await shoot(page, browserName, size.name, '19-results-descriptions');

    // The Graphic Summary of the query with 260 subjects (run 2), the popover of an HSP.
    await openResults(page, 2);
    const graphic = page.getByTestId('graphic-canvas');
    await expect(graphic).toHaveAttribute('data-rows', '100');
    const bars = JSON.parse((await graphic.getAttribute('data-targets'))!) as { x: number; y: number }[];
    await graphic.hover({ position: { x: bars[2]!.x, y: bars[2]!.y } });
    await expect(page.getByTestId('graphic-popover')).toBeVisible();
    await shoot(page, browserName, size.name, '20-results-graphic-summary');

    // The Alignments of a subject with two Ranges (run 1, subject msD), the second selected.
    await openResults(page, 1);
    await page.getByTestId('results-view-hits').click();
    await subjectRow('msD').click();
    await expect(page.locator('[data-testid^="range-0-"]')).toHaveCount(2);
    await page.locator('[data-testid^="range-0-"]').first().getByTestId('range-next').click();
    await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-state', 'ready');
    await expect(page.getByTestId('range-section')).toHaveAttribute('data-state', 'ready');
    await shoot(page, browserName, size.name, '21-results-alignments-two-ranges');

    // The dot plot's popup of the selected HSP (the pair of state 08).
    await page.getByTestId('pane-dotplot').click();
    await canvas.focus();
    await page.keyboard.press('Enter');
    await expect(page.getByTestId('dotplot-popup')).toBeVisible();
    await shoot(page, browserName, size.name, '22-results-dotplot-popup');
    await page.keyboard.press('Escape');

    // The dot plot of TBLASTN (run 4): amino acids along the query, nucleotides along the subject.
    await openResults(page, 4);
    await page.getByTestId('pane-dotplot').click();
    await expect(canvas).toHaveAttribute('data-segments', /^[1-9]/);
    await expect(page.getByTestId('dotplot-selected')).toContainText('subject frame');
    await shoot(page, browserName, size.name, '23-results-tblastn-dotplot');

    // S14. The popup of the TBLASTN HSP, which adds it to the candidates.
    const confirmation = page.getByTestId('candidates-added');
    /** The confirmation of an addition goes after a few seconds: the records show the screen without it. */
    const settledConfirmation = () => expect(confirmation).toHaveCount(0, { timeout: 10_000 });
    await canvas.focus();
    await page.keyboard.press('Enter');
    await expect(page.getByTestId('dotplot-popup')).toBeVisible();
    await shoot(page, browserName, size.name, '26-results-tblastn-dotplot-popup');
    await page.getByTestId('dotplot-popup-add').click();
    await expect(page.getByTestId('dotplot-popup-add')).toHaveText('In candidates');
    await page.keyboard.press('Escape');

    // The Descriptions of run 1 with two rows marked (not msD, whose Ranges state 25 adds one by one).
    await openResults(page, 1);
    await page.getByTestId('results-view-hits').click();
    const marks = page.getByTestId('subject-list').locator('.marked-line').filter({ hasNotText: 'msD' }).locator('[data-testid^="subject-mark-"]');
    await marks.nth(0).check();
    await marks.nth(1).check();
    await expect(page.getByTestId('descriptions-selected')).toHaveText('2 sequences selected');
    await settledConfirmation();
    await shoot(page, browserName, size.name, '24-results-descriptions-marks');
    await page.getByTestId('descriptions-add-candidates').click();
    await expect(confirmation).toContainText('added to Candidates');

    // The Alignments of msD (two Ranges): the first added, "In candidates"; the second not.
    await subjectRow('msD').click();
    const ranges = page.locator('[data-testid^="range-0-"]');
    await expect(ranges).toHaveCount(2);
    await ranges.first().locator('[data-testid^="range-add-"]').click();
    await expect(ranges.first().locator('[data-testid^="range-add-"]')).toHaveText('In candidates');
    await expect(ranges.nth(1).locator('[data-testid^="range-add-"]')).toHaveText('Add to candidates');
    await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-state', 'ready');
    // The second Range's section is read when it comes into view.
    await ranges.nth(1).scrollIntoViewIfNeeded();
    await expect(page.getByTestId('range-section')).toHaveAttribute('data-state', 'ready');
    await settledConfirmation();
    await shoot(page, browserName, size.name, '25-results-alignments-in-candidates');

    // The tray: candidates of runs 1 and 4, two notes, and the Origins.
    await page.getByTestId('tab-candidates').click();
    await expect(page.getByTestId('candidate-list')).toBeVisible();
    await page.getByTestId('candidate-note-1').fill('TBLASTN hit to check');
    await page.getByTestId('candidate-note-2').fill('compare with msA');
    await expect(page.getByTestId('candidate-origins').locator('li[data-testid^="candidate-origin-"]')).toHaveCount(2);
    await shoot(page, browserName, size.name, '27-candidates-two-runs-notes-origins');

    // An extraction with flanks that the records' ends cut: the summary lists the requested and written ranges.
    await page.getByTestId('extract-region-flanked').check();
    await page.getByTestId('extract-flank-left').fill('5000');
    await page.getByTestId('extract-flank-right').fill('5000');
    const download = page.waitForEvent('download');
    await page.getByTestId('extract-download').click();
    await download;
    await expect(page.getByTestId('extract-clipped').first()).toBeVisible();
    await page.getByTestId('extract-summary').scrollIntoViewIfNeeded();
    await shoot(page, browserName, size.name, '28-candidates-extraction-clipped');

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
    await expect(page.getByTestId('run-6-status')).toHaveText('cancelled');
    await submit(page);
    await expect(queued).toHaveCount(1);
    // The queue shows beside either tab, at the same place (W4b). State 15 shows run 4's Dot Plot
    // as S13b recorded it (24 to 28 left run 1's Alignments open).
    await openResults(page, 4);
    await page.getByTestId('pane-dotplot').click();
    await expect(canvas).toHaveAttribute('data-segments', /^[1-9]/);
    await page.getByTestId('tab-results').click();
    await shoot(page, browserName, size.name, '15-results-queue-open-waiting-cancelled');

    // The results screen for the cancelled run.
    const select = page.getByTestId('results-run');
    const label = await select.locator('option', { hasText: /^Run 6 · / }).textContent();
    await select.selectOption({ label: label!.trim() });
    await expect(page.getByTestId('results-status')).toHaveAttribute('data-run-status', 'cancelled');
    await shoot(page, browserName, size.name, '16-results-cancelled-run');

    // The dot plot of a translated search (TBLASTX, run 3): axes in nt, the frames of the selected HSP.
    const run3 = await select.locator('option', { hasText: /^Run 3 · / }).textContent();
    await select.selectOption({ label: run3!.trim() });
    await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '3', { timeout: 60_000 });
    await page.getByTestId('pane-dotplot').click();
    await expect(page.getByTestId('dotplot-canvas')).toHaveAttribute('data-segments', /^[1-9]/);
    await expect(page.getByTestId('dotplot-selected')).toContainText('frames');
    await shoot(page, browserName, size.name, '17-results-translated-dotplot');
    await page.getByTestId('results-view-hits').click();

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

// --- outputs, reproduction and session files (S15) -------------------------------------------

/** The bytes of the file that clicking a button saves (the browser's download). */
async function saveDownload(page: Page, testid: string): Promise<Buffer> {
  const download = page.waitForEvent('download');
  await page.getByTestId(testid).click();
  return readFileSync((await (await download).path())!);
}

/** A FASTA text with one residue changed: the `at`-th letter (0-based) of the record after `title`'s line. */
function oneResidueChanged(text: string, title: string, at: number): string {
  const start = text.indexOf('\n', text.indexOf(`>${title}`)) + 1;
  let position = start;
  for (let letters = 0; ; position++) {
    if (text[position] === '\n') continue;
    if (letters++ === at) break;
  }
  const letter = text[position]!;
  return text.slice(0, position) + (letter.toUpperCase() === 'A' ? 'C' : 'A') + text.slice(position + 1);
}

for (const size of SIZES) {
  test(`outputs, reproduction and session files at the ${size.name} size`, async ({ page, browserName }) => {
    await page.setViewportSize({ width: size.width, height: size.height });
    await page.goto('/');
    await expect(page.getByTestId('storage-status')).toBeVisible();
    const multiQuery = fasta('outfmt0/multi_query.fasta');
    const multiSubject = fasta('outfmt0/multi_subject.fasta');

    // Run 1: BLASTN of three queries against six subjects with the fifth (msE) left out; run 2:
    // TBLASTN with genetic code 4 for the subject.
    await program(page, 'blastn');
    await openFiles(page, 'query', [{ name: 'multi_query.fasta', text: multiQuery }]);
    await openFiles(page, 'subject', [{ name: 'multi_subject.fasta', text: multiSubject }]);
    await showRecords(page, 'subject');
    await page.getByTestId('subject-source-0-record-4').uncheck();
    await expect(page.getByTestId('subject-source-0-summary')).toContainText('(5 included)');
    await page.getByTestId('job-title').fill('Five of six subjects');
    await submit(page);
    await waitStatus(page, 1, 'completed');
    await page.getByTestId('job-title').fill('');
    await clearInputs(page);
    await program(page, 'tblastn');
    await openFiles(page, 'query', [{ name: 'e2e_protein_query.faa', text: fasta('outfmt0/e2e_protein_query.faa') }]);
    await openFiles(page, 'subject', [{ name: 'e2e_amb_subject.fna', text: fasta('outfmt0/e2e_amb_subject.fna') }]);
    await openParameters(page);
    await page.getByTestId('param-db_gencode').selectOption('4');
    await submit(page);
    await waitStatus(page, 2, 'completed');

    // The Outputs tab after a CSV export of the whole run: the summary under "LOSAT Web formats".
    await openResults(page, 1);
    await page.getByTestId('results-view-outputs').click();
    await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', '1:6');
    await saveDownload(page, 'export-csv');
    await expect(page.getByTestId('export-summary')).toHaveAttribute('data-format', 'csv');
    await shoot(page, browserName, size.name, '29-results-outputs-own-formats');

    // Run details: how to reproduce run 1 (a subject record left out) and run 2 (genetic code 4: the approved exception).
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-reproduce')).toBeVisible();
    await expect(page.getByTestId('run-input-file-subject')).toContainText('multi_subject.fasta');
    await shoot(page, browserName, size.name, '30-results-run-details-reproduce');
    await openResults(page, 2);
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-ncbi-exception')).toBeVisible();
    await shoot(page, browserName, size.name, '30b-results-run-details-reproduce-gencode');

    // "Load settings...": a BLASTN file with a query region, which a query of three records cannot take.
    await page.getByTestId('tab-search').click();
    await clearInputs(page);
    await program(page, 'blastn');
    await openFiles(page, 'query', [{ name: 'multi_query.fasta', text: multiQuery }]);
    await openFiles(page, 'subject', [{ name: 'multi_subject.fasta', text: multiSubject }]);
    const settings = JSON.parse((await saveDownload(page, 'settings-save')).toString('utf8')) as { options: string[] };
    settings.options.push('-query_loc', '3-40');
    await page.getByTestId('settings-load').setInputFiles({
      name: 'losat-settings-blastn.json',
      mimeType: 'application/json',
      buffer: Buffer.from(JSON.stringify(settings)),
    });
    await expect(page.getByTestId('settings-message')).toContainText('Not applied: -query_loc 3-40');
    await expect(page.getByTestId('settings-message')).toHaveAttribute('data-partial', 'true');
    await shoot(page, browserName, size.name, '31-search-settings-loaded-not-applied');

    // The session panel before saving, with candidates of run 1 and a note.
    await openResults(page, 1);
    await page.getByTestId('results-view-hits').click();
    await page.getByTestId('descriptions-select-all').check();
    await page.getByTestId('descriptions-add-candidates').click();
    await expect(page.getByTestId('candidates-added')).toContainText('added to Candidates');
    await page.getByTestId('tab-candidates').click();
    await page.getByTestId('candidate-note-1').fill('check the minus-strand copy');
    await expect(page.getByTestId('candidates-added')).toHaveCount(0, { timeout: 10_000 });
    await page.getByTestId('session-panel').scrollIntoViewIfNeeded();
    await shoot(page, browserName, size.name, '32-session-panel-before-saving');
    const session = await saveDownload(page, 'session-save');
    await expect(page.getByTestId('session-message')).toContainText('Saved ');

    // A new working session: the empty tray (a W5 state), then a refused session file.
    await page.reload();
    await expect(page.getByTestId('storage-status')).toBeVisible();
    await page.getByTestId('tab-candidates').click();
    await expect(page.getByTestId('candidates-empty')).toBeVisible();
    await shoot(page, browserName, size.name, '39-candidates-empty');
    const openSession = async (name: string, buffer: Buffer) => {
      await page.getByTestId('session-file').setInputFiles({ name, mimeType: 'application/gzip', buffer });
      await expect(page.getByTestId('session-message')).toBeVisible({ timeout: 60_000 });
    };
    await openSession('cut.losat-session.gz', session.subarray(0, Math.floor(session.length / 2)));
    await expect(page.getByTestId('session-message')).toHaveAttribute('data-error', 'true');
    await page.getByTestId('session-panel').scrollIntoViewIfNeeded();
    await shoot(page, browserName, size.name, '36-session-refused');

    // The session opened: the queue says where the runs come from, and so does the results header.
    await openSession('w6-states.losat-session.gz', session);
    await expect(page.getByTestId('session-message')).toHaveAttribute('data-error', 'false', { timeout: 60_000 });
    await expect(page.getByTestId('run-2-origin')).toBeVisible();
    await openResults(page, 1);
    await expect(page.getByTestId('results-origin')).toBeVisible();
    await shoot(page, browserName, size.name, '33-session-opened-queue-results-origin');

    // Run details of a loaded run, its originals not attached.
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-original-subject')).toHaveAttribute('data-attached', 'false');
    await shoot(page, browserName, size.name, '34-results-run-details-loaded');

    // The tray's extract form: the sequences need the original FASTA; then a flank that is not a number (a W5 state).
    await page.getByTestId('tab-candidates').click();
    await expect(page.getByTestId('extract-originals')).toBeVisible();
    await page.getByTestId('extract-form').scrollIntoViewIfNeeded();
    await shoot(page, browserName, size.name, '35-candidates-extract-needs-original');
    await page.getByTestId('extract-region-flanked').check();
    await page.getByTestId('extract-flank-left').fill('-5');
    await expect(page.getByTestId('extract-flank-error')).toBeVisible();
    await shoot(page, browserName, size.name, '40-candidates-flank-error');
    await page.getByTestId('extract-flank-left').fill('0');

    // Run details after re-attachment: the subject refused with one residue changed, then accepted;
    // the query refused (its message stays).
    await openResults(page, 1);
    await page.getByTestId('results-view-details').click();
    const attach = (role: string, name: string, text: string | Buffer) =>
      page.getByTestId(`run-attach-${role}-files`).setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(text) });
    await attach('subject', 'multi_subject.fasta', oneResidueChanged(multiSubject.toString('utf8'), 'msA', 100));
    await expect(page.getByTestId('run-attach-subject-message')).toBeVisible();
    await attach('subject', 'multi_subject.fasta', multiSubject);
    await expect(page.getByTestId('run-original-subject')).toHaveAttribute('data-attached', 'true');
    await attach('query', 'multi_query.fasta', oneResidueChanged(multiQuery.toString('utf8'), 'mq1', 10));
    await expect(page.getByTestId('run-attach-query-message')).toBeVisible();
    await shoot(page, browserName, size.name, '37-results-run-details-reattached');

    // Edit Search of run 2: its settings in the search form, and the message that says so.
    await openResults(page, 2);
    await page.getByTestId('edit-search').click();
    await expect(page.getByTestId('tab-search')).toHaveAttribute('aria-pressed', 'true');
    await expect(page.getByTestId('settings-message')).toContainText('The search form has the settings of Run 2.');
    await shoot(page, browserName, size.name, '38-search-after-edit-search');
  });
}
