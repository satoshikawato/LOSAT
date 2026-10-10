// Session files through the real application (S15 items 4 and 6, design §12.2, REQ-23): runs,
// candidates and notes are saved, the page is loaded again (a new working session), and the file
// is opened: the runs, their results and the tray come back and no search runs; the extraction of
// sequences waits for the original FASTA, which is attached only when it matches; a damaged file
// is refused and changes nothing; names, IDs, titles and notes that look like HTML, scripts or
// URLs are shown as text. The FakeEngine build: its HSPs are known here.
import { readFile } from 'node:fs/promises';
import { expect, test, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { chooseRun, openFiles, paste, program, showOutput, submit, waitStatus } from './support/search';

test.setTimeout(180_000);
test.skip(BUILD_HAS_ENGINE, "the expected HSPs and files follow the FakeEngine's");

test.beforeEach(async ({ page }) => {
  await page.goto('/');
});

/** Pseudo-random nucleotides: the same letters in every run. */
function dna(seed: number, length: number): string {
  let state = seed;
  let out = '';
  for (let i = 0; i < length; i++) {
    state = (Math.imul(state, 1103515245) + 12345) >>> 0;
    out += 'ACGT'[(state >>> 16) % 4];
  }
  return out;
}

const record = (title: string, letters: string) => `>${title}\n${letters.match(/.{1,30}/g)!.join('\n')}\n`;
const S1 = dna(51, 40);
const S2 = dna(52, 36);
const SUBJECTS = record('s1 first subject', S1) + record('s2', S2);
const QUERIES = record('q1 first query', dna(61, 30)) + record('q2', dna(62, 24));
const PROTEIN = record('p1', 'MKVLAAGIVGLLLAHHKKEE');

interface Download {
  readonly name: string;
  readonly bytes: Buffer;
}

/** The file that clicking a button saves (the browser's download). */
async function downloaded(page: Page, testid: string): Promise<Download> {
  const download = page.waitForEvent('download');
  await page.getByTestId(testid).click();
  const file = await download;
  return { name: file.suggestedFilename(), bytes: await readFile((await file.path())!) };
}

async function search(page: Page, number: number, title: string): Promise<void> {
  await page.getByTestId('job-title').fill(title);
  await submit(page);
  await waitStatus(page, number, 'completed');
}

/** Run 1: BLASTN of the pasted queries against subjects.fa (opened as a file). */
async function blastnRun(page: Page, subjects = SUBJECTS, title = 'First run'): Promise<void> {
  await program(page, 'blastn');
  await paste(page, 'query', QUERIES);
  await openFiles(page, 'subject', [{ name: 'subjects.fa', text: subjects }]);
  await search(page, 1, title);
}

/** Adds every HSP of the selected query of the shown run from the Descriptions' marks. */
async function addAll(page: Page, expected: string): Promise<void> {
  await page.getByTestId('results-view-hits').click();
  await page.getByTestId('descriptions-select-all').check();
  await page.getByTestId('descriptions-add-candidates').click();
  await expect(page.getByTestId('candidates-added')).toHaveText(expected);
}

async function openSession(page: Page, file: Download, name = file.name): Promise<void> {
  await page.getByTestId('session-file').setInputFiles({ name, mimeType: 'application/gzip', buffer: file.bytes });
  await expect(page.getByTestId('session-busy')).toHaveCount(0, { timeout: 30_000 });
}

const message = (page: Page) => page.getByTestId('session-message');
const queueRuns = (page: Page) => page.getByTestId('queue').locator('li[data-status]');

test('saves the runs, results, candidates and notes, and opens them in a new page without searching', async ({ page }) => {
  await blastnRun(page);
  await page.getByTestId('tab-results').click();
  await addAll(page, '3 HSPs added to Candidates.');
  await page.getByTestId('tab-search').click();
  await program(page, 'tblastn');
  await paste(page, 'query', PROTEIN);
  await search(page, 2, 'Second run');
  await chooseRun(page, 2);
  await addAll(page, '3 HSPs added to Candidates.');
  await page.getByTestId('tab-candidates').click();
  await page.getByTestId('candidate-note-2').fill('check the second exon');
  const aligned = await downloaded(page, 'extract-aligned');
  await showOutput(page, 1, 6);
  const out6 = await page.getByTestId('result-output').textContent();
  await showOutput(page, 2, 0);
  const out0 = await page.getByTestId('result-output').textContent();

  await expect(page.getByTestId('session-include-candidates')).toBeChecked();
  await expect(page.getByTestId('session-save-note')).toContainText('2 completed runs will be saved.');
  const session = await downloaded(page, 'session-save');
  expect(session.name).toMatch(/^losat-session-\d{8}-\d{6}\.losat-session\.gz$/);
  expect([...session.bytes.subarray(0, 2)]).toEqual([0x1f, 0x8b]);
  await expect(message(page)).toHaveText(`Saved ${session.name}: 2 runs and 6 candidates with their notes.`);

  // A new working session: nothing is there until the file is opened.
  await page.reload();
  await expect(queueRuns(page)).toHaveCount(0);
  await openSession(page, session);
  await expect(message(page)).toHaveText(
    `Opened ${session.name}: 2 runs (runs 1 to 2 here) and 6 candidates with their notes. Nothing was searched again.`,
  );
  await expect(message(page)).toHaveAttribute('data-error', 'false');
  // The queue has the two loaded runs and nothing else: no search was queued or run.
  await expect(queueRuns(page)).toHaveCount(2);
  for (const n of [1, 2]) {
    await expect(page.getByTestId(`run-${n}-status`)).toHaveText('completed');
    await expect(page.getByTestId(`run-${n}-origin`)).toHaveText(`From ${session.name}, run ${n} there (not searched again)`);
  }
  await expect(page.getByTestId('run-2-title')).toHaveText('Second run');
  await expect(page.getByTestId('queue').locator('li[data-status="queued"], li[data-status="running"], li[data-status="preparing"]')).toHaveCount(0);

  // The results are those saved, byte for byte.
  await showOutput(page, 1, 6);
  await expect(page.getByTestId('result-output')).toHaveText(out6!);
  await expect(page.getByTestId('results-origin')).toHaveText(`from ${session.name}, run 1 there (not searched again)`);
  await showOutput(page, 2, 0);
  await expect(page.getByTestId('result-output')).toHaveText(out0!);

  // The tray and its notes; the sequences wait for the originals, the aligned rows do not.
  await expect(page.getByTestId('tab-candidates-count')).toHaveText('6');
  await page.getByTestId('tab-candidates').click();
  await expect(page.getByTestId('candidate-note-2')).toHaveValue('check the second exon');
  await expect(page.getByTestId('candidate-1').locator('[data-field="run"]')).toHaveText('Run 1 First run');
  await expect(page.getByTestId('extract-originals')).toHaveText(
    'Runs 1 and 2 were loaded from a session file; choose their original subject FASTA in Run details to extract sequences.',
  );
  await expect(page.getByTestId('extract-download')).toBeDisabled();
  const again = await downloaded(page, 'extract-aligned');
  expect(again.name).toBe('losat-candidates-aligned.fa');
  expect(again.bytes.equals(aligned.bytes)).toBe(true);
  await expect(queueRuns(page)).toHaveCount(2);
});

test('a damaged or foreign file is refused with its reason and changes nothing', async ({ page }) => {
  await blastnRun(page);
  await page.getByTestId('tab-results').click();
  await addAll(page, '3 HSPs added to Candidates.');
  const session = await downloaded(page, 'session-save');

  // Cut in the middle: the gzip data ends early.
  await openSession(page, { name: 'cut.losat-session.gz', bytes: session.bytes.subarray(0, Math.floor(session.bytes.length / 2)) });
  await expect(message(page)).toHaveAttribute('data-error', 'true');
  await expect(message(page)).toContainText('cut.losat-session.gz was not opened, and nothing was loaded.');
  await expect(message(page)).toContainText(/The session file is (damaged|incomplete)/);
  // A byte changed inside the compressed data.
  const flipped = Buffer.from(session.bytes);
  flipped[Math.floor(flipped.length / 2)]! ^= 0x40;
  await openSession(page, { name: 'flipped.gz', bytes: flipped });
  await expect(message(page)).toContainText('flipped.gz was not opened, and nothing was loaded.');
  // A FASTA file is not a session file.
  await openSession(page, { name: 'subjects.fa', bytes: Buffer.from(SUBJECTS) });
  await expect(message(page)).toHaveText(
    'subjects.fa was not opened, and nothing was loaded. This is not a LOSAT Web session file: it is not gzip data (a session file is saved as a .losat-session.gz file).',
  );
  // Nothing changed: the one run of this page, its three candidates.
  await expect(queueRuns(page)).toHaveCount(1);
  await expect(page.getByTestId('run-1-origin')).toHaveCount(0);
  await expect(page.getByTestId('tab-candidates-count')).toHaveText('3');
  // The intact file opens after the refusals.
  await openSession(page, session);
  await expect(message(page)).toHaveAttribute('data-error', 'false');
  await expect(queueRuns(page)).toHaveCount(2);
  await expect(page.getByTestId('run-2-origin')).toHaveText(`From ${session.name}, run 1 there (not searched again)`);
});

test('names, IDs, titles and notes that look like HTML, scripts or URLs are shown as text, and run nothing', async ({ page }) => {
  const dialogs: string[] = [];
  page.on('dialog', (dialog) => {
    dialogs.push(dialog.message());
    void dialog.dismiss();
  });
  const id = 'q<img/src=x/onerror=window.__pwned=1>';
  const title = '<script>window.__pwned=3</script> javascript:window.__pwned=4';
  const fileName = '<img src=x onerror=window.__pwned=2>.fa';
  const note = '<img src=x onerror="window.__pwned=5"> javascript:alert(5)';
  await program(page, 'blastn');
  await paste(page, 'query', record(`${id} description`, dna(61, 30)));
  await openFiles(page, 'subject', [{ name: fileName, text: SUBJECTS }]);
  await search(page, 1, title);
  await page.getByTestId('tab-results').click();
  await addAll(page, '3 HSPs added to Candidates.');
  await page.getByTestId('tab-candidates').click();
  await page.getByTestId('candidate-note-1').fill(note);
  const session = await downloaded(page, 'session-save');

  await page.reload();
  const sessionName = '<img src=x onerror=window.__pwned=6>.gz';
  await openSession(page, session, sessionName);
  await expect(message(page)).toContainText(`Opened ${sessionName}: 1 run`);
  await expect(page.getByTestId('run-1-title')).toHaveText(title);
  await expect(page.getByTestId('run-1-origin')).toHaveText(`From ${sessionName}, run 1 there (not searched again)`);
  await expect(page.getByTestId('queue')).toContainText(`query.fa vs ${fileName}`);
  await page.getByTestId('run-1-open').click();
  await expect(page.getByTestId('results-job-title')).toHaveText(title);
  await expect(page.getByTestId('results-query-id')).toHaveText(id);
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-input-subject')).toContainText(`${fileName}: 2 records`);
  await expect(page.getByTestId('run-origin-file')).toContainText(`From ${sessionName}, run 1 there`);
  await page.getByTestId('tab-candidates').click();
  await expect(page.getByTestId('candidate-note-1')).toHaveValue(note);
  await expect(page.getByTestId('candidate-1').locator('[data-field="query"]')).toHaveText(id);

  expect(await page.evaluate(() => (window as unknown as { __pwned?: unknown }).__pwned)).toBeUndefined();
  await expect(page.locator('img[src="x"], script:not([src])')).toHaveCount(0);
  expect(dialogs).toEqual([]);
});

test('Run details attaches the original FASTA only when it matches; then sequences are extracted as before saving', async ({ page }) => {
  await blastnRun(page);
  await page.getByTestId('tab-results').click();
  await addAll(page, '3 HSPs added to Candidates.');
  await page.getByTestId('tab-candidates').click();
  const before = await downloaded(page, 'extract-download');
  expect(before.name).toBe('losat-candidates.fa');
  const session = await downloaded(page, 'session-save');

  await page.reload();
  await openSession(page, session);
  await page.getByTestId('run-1-open').click();
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-origin')).toBeVisible();
  await expect(page.getByTestId('run-input-subject')).toContainText('subjects.fa: 2 records');
  const subject = page.getByTestId('run-original-subject');
  await expect(subject).toHaveAttribute('data-attached', 'false');
  await expect(page.getByTestId('run-original-subject-missing')).toContainText('the subject sequences of this run\'s candidates (hit regions, flanks, complete sequences) cannot be extracted');

  // A file of the same name with one residue changed: refused, and the record is named.
  const changed = SUBJECTS.replace(S2.slice(0, 30), `${S2[0] === 'A' ? 'C' : 'A'}${S2.slice(1, 30)}`);
  expect(changed).not.toBe(SUBJECTS);
  await page.getByTestId('run-attach-subject-files').setInputFiles({ name: 'subjects.fa', mimeType: 'text/plain', buffer: Buffer.from(changed) });
  await expect(page.getByTestId('run-attach-subject-message')).toHaveText(
    'The subject FASTA was not attached to run 1: record 2 ("s2") differs from the saved run\'s record: its bytes have another SHA-256.',
  );
  await expect(subject).toHaveAttribute('data-attached', 'false');
  // One record more.
  await page
    .getByTestId('run-attach-subject-files')
    .setInputFiles({ name: 'subjects.fa', mimeType: 'text/plain', buffer: Buffer.from(`${SUBJECTS}>s3\nACGT\n`) });
  await expect(page.getByTestId('run-attach-subject-message')).toContainText('"subjects.fa" has 3 records, but file 1 of the saved input ("subjects.fa") had 2');
  await page.getByTestId('tab-candidates').click();
  await expect(page.getByTestId('extract-originals')).toHaveText(
    'Run 1 was loaded from a session file; choose its original subject FASTA in Run details to extract sequences.',
  );

  // The right file, even under another name.
  await page.getByTestId('tab-results').click();
  await page.getByTestId('results-view-details').click();
  await page.getByTestId('run-attach-subject-files').setInputFiles({ name: 'renamed.fa', mimeType: 'text/plain', buffer: Buffer.from(SUBJECTS) });
  await expect(subject).toHaveAttribute('data-attached', 'true');
  await expect(page.getByTestId('run-original-subject-attached')).toContainText('Attached: renamed.fa.');
  await expect(page.getByTestId('run-attach-subject-message')).toHaveCount(0);
  await expect(page.getByTestId('run-original-query')).toHaveAttribute('data-attached', 'false');

  await page.getByTestId('tab-candidates').click();
  await expect(page.getByTestId('extract-originals')).toHaveCount(0);
  const after = await downloaded(page, 'extract-download');
  expect(after.bytes.toString()).toBe(before.bytes.toString());
  // The query was pasted and is not attached: its sequences still wait for it.
  await page.getByTestId('extract-role-query').check();
  await expect(page.getByTestId('extract-originals')).toHaveText(
    'Run 1 was loaded from a session file; choose its original query FASTA in Run details to extract sequences.',
  );
  await expect(page.getByTestId('extract-download')).toBeDisabled();
});
