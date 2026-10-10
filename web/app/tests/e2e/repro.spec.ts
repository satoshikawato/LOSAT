// Settings files, "Edit Search" and the reproduction panel of a run (S15 items 3 and 5;
// docs/web/ncbi_ui_mapping.md "Edit Search"), through the real application: the search form's
// settings saved, changed and loaded back, broken files refused; "Edit Search" puts a run's
// settings and Job Title in the form and shows the Search tab without searching; Run details
// gives the LOSAT and NCBI commands from the run's argv and saves the run's input FASTA, the
// bytes that the engine searched, under the argv's names.
import { readFile } from 'node:fs/promises';
import { expect, test, type Page } from '@playwright/test';
import { chooseRun, openFiles, openParameters, paste, program, submit, waitStatus } from './support/search';

test.setTimeout(120_000);

test.beforeEach(async ({ page }) => {
  await page.goto('/');
});

const QUERY = '>q1 query\nACGTACGTTGCAACGTGGCCAATTGGCCAACGTACGTTGCAACGTGGCCAATT\n';
const SUBJECT_A = '>s1 first\nACGTACGTTGCAACGTGGCCAATTGGCC\n>s2\nTTGGCCAACGTACGTTGCAACG\n';
const SUBJECT_B = '>s3 third\nGGCCAATTGGCCAACGTACGTTGCAACGTGG\n';

/** The name and the bytes of the file that clicking a button saves (the browser's download). */
async function saved(page: Page, testid: string): Promise<{ name: string; bytes: Buffer }> {
  const download = page.waitForEvent('download');
  await page.getByTestId(testid).click();
  const file = await download;
  return { name: file.suggestedFilename(), bytes: await readFile((await file.path())!) };
}

/** The page is no wider than a phone's window (as results.spec.ts checks the results screen). */
async function expectNoSideScroll(page: Page, state: string): Promise<void> {
  await page.setViewportSize({ width: 390, height: 844 });
  const { scrollWidth, innerWidth } = await page.evaluate(() => ({ scrollWidth: document.documentElement.scrollWidth, innerWidth: window.innerWidth }));
  expect(scrollWidth, `${state}: the page is wider than the window`).toBeLessThanOrEqual(innerWidth);
  await page.setViewportSize({ width: 1280, height: 900 });
}

async function validated(page: Page): Promise<void> {
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
}

test('settings: saved from the form, the form changed, loaded back; broken and newer files are refused', async ({ page }) => {
  await paste(page, 'query', QUERY);
  await paste(page, 'subject', SUBJECT_B);
  await page.getByTestId('param-task-blastn').check();
  await openParameters(page);
  await page.getByTestId('param-evalue').fill('1e-5');
  await page.getByTestId('param-penalty').fill('-3');
  await page.getByTestId('param-lcase_masking').check();
  await page.getByTestId('query-region-start').fill('3');
  await page.getByTestId('query-region-stop').fill('40');
  await page.getByTestId('threads').selectOption('2');
  await page.getByTestId('job-title').fill('Not in the file');
  await validated(page);

  const file = await saved(page, 'settings-save');
  expect(file.name).toBe('losat-settings-blastn.json');
  const settings = JSON.parse(file.bytes.toString('utf8'));
  expect(settings).toMatchObject({
    format: 'LOSAT Web search settings',
    schema: 1,
    program: 'blastn',
    options: ['-task', 'blastn', '-evalue', '1e-5', '-penalty', '-3', '-lcase_masking', '-query_loc', '3-40'],
    threads: 2,
  });
  // Only the search conditions: no input, name or title.
  expect(file.bytes.toString('utf8')).not.toMatch(/query\.fa|subject\.fa|Not in the file|ACGT/);
  await expect(page.getByTestId('settings-message')).toHaveText('Saved losat-settings-blastn.json: the program, the options and the threads.');
  await expectNoSideScroll(page, 'the search form with its settings line');

  // Change the form, then load the file back: the form has its conditions again; inputs and title stay.
  await page.getByTestId('param-task-megablast').check();
  await page.getByTestId('param-evalue').fill('5');
  await page.getByTestId('param-penalty').fill('');
  await page.getByTestId('param-word_size').fill('20');
  await page.getByTestId('param-lcase_masking').uncheck();
  await page.getByTestId('query-region-start').fill('');
  await page.getByTestId('threads').selectOption('auto');
  await page.getByTestId('settings-load').setInputFiles({ name: 'mine.json', mimeType: 'application/json', buffer: file.bytes });
  await expect(page.getByTestId('settings-message')).toHaveText(
    'Loaded mine.json: BLASTN, 9 words of options, threads 2. The inputs and the Job Title did not change.',
  );
  await expect(page.getByTestId('settings-message')).toHaveAttribute('data-kind', 'info');
  await expect(page.getByTestId('param-task-blastn')).toBeChecked();
  await expect(page.getByTestId('param-evalue')).toHaveValue('1e-5');
  await expect(page.getByTestId('param-penalty')).toHaveValue('-3');
  await expect(page.getByTestId('param-word_size')).toHaveValue('');
  await expect(page.getByTestId('param-lcase_masking')).toBeChecked();
  await expect(page.getByTestId('query-region-start')).toHaveValue('3');
  await expect(page.getByTestId('query-region-stop')).toHaveValue('40');
  await expect(page.getByTestId('threads')).toHaveValue('2');
  await expect(page.getByTestId('job-title')).toHaveValue('Not in the file');
  await expect(page.getByTestId('query-input')).toHaveValue(QUERY);
  await validated(page);

  // A broken file and a file of a newer schema are refused; the form stays as it is.
  await page.getByTestId('param-evalue').fill('2');
  await page.getByTestId('settings-load').setInputFiles({ name: 'broken.json', mimeType: 'application/json', buffer: Buffer.from('{"format": ') });
  await expect(page.getByTestId('settings-message')).toHaveAttribute('data-kind', 'error');
  await expect(page.getByTestId('settings-message')).toContainText('broken.json was not loaded: The file is not JSON (');
  const newer = Buffer.from(JSON.stringify({ ...settings, schema: 2 }));
  await page.getByTestId('settings-load').setInputFiles({ name: 'newer.json', mimeType: 'application/json', buffer: newer });
  await expect(page.getByTestId('settings-message')).toHaveText(
    'newer.json was not loaded: The file has schema 2: a newer LOSAT Web saved it. This one reads schema 1.',
  );
  const reserved = Buffer.from(JSON.stringify({ ...settings, options: ['-out', 'x.txt'] }));
  await page.getByTestId('settings-load').setInputFiles({ name: 'reserved.json', mimeType: 'application/json', buffer: reserved });
  await expect(page.getByTestId('settings-message')).toContainText('reserved.json was not loaded: "options" word 1 is -out, which LOSAT Web sets itself');
  await expect(page.getByTestId('param-evalue')).toHaveValue('2');
  await expect(page.getByTestId('param-task-blastn')).toBeChecked();
});

test("Run details: the commands from the run's argv, the run's input FASTA and settings; Edit Search fills the form without searching", async ({
  page,
}) => {
  await paste(page, 'query', QUERY);
  await openFiles(page, 'subject', [
    { name: 'a.fa', text: SUBJECT_A },
    { name: 'b.fa', text: SUBJECT_B },
  ]);
  await openParameters(page);
  await page.getByTestId('param-evalue').fill('1e-5');
  await page.getByTestId('job-title').fill('Repro run');
  await validated(page);
  await submit(page);
  await waitStatus(page, 1, 'completed');

  // The form changes after the run: the panel still shows the run's own commands.
  await page.getByTestId('param-evalue').fill('7');
  await chooseRun(page, 1);
  await page.getByTestId('results-view-details').click();
  for (const format of [0, 6, 7]) {
    await expect(page.getByTestId(`run-command-${format}`)).toHaveText(
      `LOSAT blastn -query query.fa -subject combined_subject.fa -evalue 1e-5 -outfmt ${format}`,
    );
    await expect(page.getByTestId(`run-ncbi-command-${format}`)).toHaveText(
      `blastn -query query.fa -subject combined_subject.fa -evalue 1e-5 -outfmt ${format}`,
    );
  }
  await expect(page.getByTestId('run-ncbi-unavailable')).toHaveCount(0);
  await expect(page.getByTestId('run-reproduce-notes')).toContainText('Put the files query.fa and combined_subject.fa in one folder');
  await expect(page.getByTestId('run-reproduce-notes')).toContainText('do not set -num_threads');
  await expect(page.getByTestId('run-input-file-query')).toContainText('query.fa is the pasted query text, as the run searched it (1 record).');
  await expect(page.getByTestId('run-input-file-subject')).toContainText('combined_subject.fa joins the 2 subject inputs (a.fa, b.fa) in the order chosen');
  await expectNoSideScroll(page, 'the reproduction panel');

  // The input FASTA: the bytes that the engine searched, named as in the argv.
  const query = await saved(page, 'run-input-save-query');
  expect(query.name).toBe('query.fa');
  expect(query.bytes.toString('utf8')).toBe(QUERY);
  const subject = await saved(page, 'run-input-save-subject');
  expect(subject.name).toBe('combined_subject.fa');
  expect(subject.bytes.toString('utf8')).toBe(SUBJECT_A + SUBJECT_B);
  await expect(page.getByTestId('run-files-message')).toHaveText('Saved combined_subject.fa: the subject that Run 1 searched.');
  const settings = await saved(page, 'run-settings-save');
  expect(settings.name).toBe('losat-settings-blastn.json');
  expect(JSON.parse(settings.bytes.toString('utf8'))).toMatchObject({ program: 'blastn', options: ['-evalue', '1e-5'], threads: 'auto' });

  // Edit Search: the run's settings and Job Title in the form, the Search tab shown, the inputs kept, nothing queued.
  await page.getByTestId('tab-search').click();
  await page.getByTestId('job-title').fill('Other title');
  await page.getByTestId('threads').selectOption('3');
  await page.getByTestId('tab-results').click();
  await page.getByTestId('edit-search').click();
  await expect(page.getByTestId('tab-search')).toHaveAttribute('aria-pressed', 'true');
  await expect(page.getByTestId('settings-message')).toHaveText("The search form has the settings of Run 1. The inputs are the form's own.");
  await expect(page.getByTestId('settings-message')).toBeInViewport();
  await expect(page.getByTestId('param-evalue')).toHaveValue('1e-5');
  await expect(page.getByTestId('job-title')).toHaveValue('Repro run');
  await expect(page.getByTestId('threads')).toHaveValue('auto');
  await expect(page.getByTestId('query-input')).toHaveValue(QUERY);
  await expect(page.locator('[data-testid^="subject-source-"][data-status]')).toHaveCount(2);
  await expect(page.getByTestId('run-2')).toHaveCount(0);
});

test('a run that NCBI BLAST+ cannot run as LOSAT did gets the reason instead of an NCBI command', async ({ page }) => {
  await program(page, 'tblastn');
  await paste(page, 'query', '>p1\nMKVLAAGIVGLLLAHHKKEEDDPPWWRRSS\n');
  await paste(page, 'subject', SUBJECT_B);
  await page.getByTestId('param-db_gencode').selectOption('32');
  await validated(page);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-command-6')).toHaveText('LOSAT tblastn -query query.fa -subject subject.fa -db_gencode 32 -outfmt 6');
  await expect(page.getByTestId('run-ncbi-command-6')).toHaveCount(0);
  await expect(page.getByTestId('run-ncbi-unavailable')).toHaveText(
    'NCBI BLAST+ 2.17.0 does not accept -db_gencode 32: its command line takes the genetic codes 1-6, 9-16, 21-31 and 33. ' +
      'LOSAT searched with it (PD-TLOSAN-LOCAL-GENCODE-32), so there is no NCBI command to compare with.',
  );
  await expect(page.getByTestId('run-ncbi-exception')).toContainText('Approved exception (PD-TLOSAN-LOCAL-GENCODE-32)');
  await expectNoSideScroll(page, 'the reproduction panel without an NCBI command');
});
