// Session files through the real application (S15 items 4 and 6, design §12.2, REQ-23): runs,
// candidates and notes are saved, the page is loaded again (a new working session), and the file
// is opened: the runs, their results and the tray come back and no search runs; the extraction of
// sequences waits for the original FASTA, which is attached only when it matches; a damaged file
// is refused and changes nothing; names, IDs, titles and notes that look like HTML, scripts or
// URLs are shown as text. The FakeEngine build: its HSPs are known here.
//
// The engine build (plan §7, the S15 row's completion conditions) saves real searches - a BLASTN
// run of two subject files joined and a TBLASTN run with a subject record left out - and opens
// them in a new page: no search runs (no run passes through the queue's other states, and no
// engine worker starts); the outputs 0, 6 and 7 exported from the loaded runs are the bytes
// exported before saving and the native CLI's for the run's argv and input FASTA; what needs the
// original FASTA says so, and the rest works; the same files in the search form attach nothing,
// a file with one residue changed is refused, and the right files attach, after which the
// extraction gives the bytes it gave before saving; damaged files are refused and change nothing.
import { createHash } from 'node:crypto';
import { mkdirSync, writeFileSync } from 'node:fs';
import { readFile } from 'node:fs/promises';
import { join } from 'node:path';
import { expect, test, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NATIVE, nativeExpectation, NO_NATIVE_REASON } from './support/native';
import { chooseRun, fasta, openFiles, paste, program, showOutput, showRecords, submit, waitStatus } from './support/search';

test.setTimeout(180_000);

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
async function addAll(page: Page, expected: string | RegExp): Promise<void> {
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

// --- FakeEngine build -----------------------------------------------------------------------------

test.describe('FakeEngine build', () => {
  test.skip(BUILD_HAS_ENGINE, "the expected HSPs and files follow the FakeEngine's");

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
});

// --- engine build --------------------------------------------------------------------------------

test.describe('engine build', () => {
  test.skip(!BUILD_HAS_ENGINE, 'real searches need the engine (LOSAT_WEB_REACTORS)');
  // Two searches, the native CLI, and the session file's round trip.
  test.setTimeout(300_000);

  type Format = 0 | 6 | 7;
  const FORMATS: readonly Format[] = [0, 6, 7];
  const reverseComplement = (letters: string) =>
    [...letters]
      .reverse()
      .map((letter) => ({ A: 'T', C: 'G', G: 'C', T: 'A' })[letter] ?? letter)
      .join('');
  const sha256 = (bytes: Uint8Array) => createHash('sha256').update(bytes).digest('hex');

  // BLASTN: r1 lies forward in sA (part1.fa) and reverse-complemented in sB (part2.fa), r2 in sB,
  // and sC is unrelated; the two subject files are joined (combined_subject.fa).
  const R1 = dna(71, 150);
  const R2 = dna(72, 120);
  const READS = record('r1 first read', R1) + record('r2', R2);
  const SA = dna(73, 80) + R1 + dna(74, 60);
  const PART1 = record('sA first part', SA);
  /** part1.fa with one residue of r1's copy changed. */
  const PART1_CHANGED = record('sA first part', SA.slice(0, 99) + (SA[99] === 'A' ? 'C' : 'A') + SA.slice(100));
  const PART2 = record('sB second part', dna(75, 50) + reverseComplement(R1) + dna(76, 40) + R2 + dna(77, 30)) + record('sC unrelated', dna(78, 200));
  // TBLASTN: eight proteins against five windows of a genome, the fifth left out.
  const PROTEINS = fasta('outfmt0/e2e_protein_query.faa');
  const WINDOWS = fasta('outfmt0/e2e_amb_subject.fna');
  const NOTE = 'the minus-strand copy in part2.fa';

  async function clearInputs(page: Page): Promise<void> {
    for (const role of ['query', 'subject'] as const) {
      const sources = page.locator(`[data-testid^="${role}-source-"][data-status]`);
      while ((await sources.count()) > 0) await page.getByTestId(`${role}-source-0-remove`).click();
    }
  }

  async function details(page: Page, run: number): Promise<void> {
    await chooseRun(page, run);
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-details')).toBeVisible();
  }

  /** The file that "Export outfmt N" of the Outputs view saves, for each format of a run. */
  async function exportedOutputs(page: Page, run: number): Promise<Record<Format, Buffer>> {
    const files: Partial<Record<Format, Buffer>> = {};
    for (const format of FORMATS) {
      await showOutput(page, run, format);
      const file = await downloaded(page, 'export-output');
      expect(file.name).toMatch(new RegExp(`^losat-run${run}-t?blastn\\.outfmt${format}\\.txt$`));
      files[format] = file.bytes;
    }
    return files as Record<Format, Buffer>;
  }

  interface RunInputs {
    readonly argv: readonly string[];
    readonly query: Download;
    readonly subject: Download;
    /** What the reproduction panel says of each input's origin. */
    readonly relations: readonly string[];
  }

  const relations = (page: Page) =>
    Promise.all((['query', 'subject'] as const).map(async (role) => (await page.getByTestId(`run-input-file-${role}`).textContent())!.trim()));

  /** A run's argv, and its input FASTA as the reproduction panel saves them: the bytes searched, under the argv's names. */
  async function runInputs(page: Page, run: number): Promise<RunInputs> {
    await details(page, run);
    const argv = (await page.getByTestId('run-argv').textContent())!.split(' ');
    const query = await downloaded(page, 'run-input-save-query');
    const subject = await downloaded(page, 'run-input-save-subject');
    expect(argv.slice(1, 5)).toEqual(['-query', query.name, '-subject', subject.name]);
    return { argv, query, subject, relations: await relations(page) };
  }

  /** The native CLI's outputs (SHA-256 by format) of the argv, run in a folder that holds the run's input FASTA. */
  function nativeOutputs(run: number, inputs: RunInputs): Record<Format, string> {
    const folder = test.info().outputPath(`native-run${run}`);
    mkdirSync(folder, { recursive: true });
    for (const file of [inputs.query, inputs.subject]) writeFileSync(join(folder, file.name), file.bytes);
    return nativeExpectation(inputs.argv, folder).sha256;
  }

  /** CSV, JSON and the report of a whole run, as the Outputs tab's "LOSAT Web formats" save them. */
  async function ownFormats(page: Page, run: number): Promise<{ hsps: number; csv: Download; json: Download; report: Download }> {
    await showOutput(page, run, 6);
    const hsps = Number(await page.getByTestId('export-scope-all').getAttribute('data-count'));
    expect(hsps).toBeGreaterThan(0);
    await expect(page.getByTestId('export-scope-all')).toBeChecked();
    const files: Partial<Record<'csv' | 'json' | 'report', Download>> = {};
    for (const format of ['csv', 'json', 'report'] as const) {
      files[format] = await downloaded(page, `export-${format}`);
      const summary = page.getByTestId('export-summary');
      await expect(summary).toHaveAttribute('data-format', format);
      await expect(summary).toHaveAttribute('data-hsps', String(hsps));
      await expect(summary).toContainText(`Saved ${files[format]!.name}: `);
    }
    return { hsps, csv: files.csv!, json: files.json!, report: files.report! };
  }

  /** Records every status that the queue shows from now on: a MutationObserver sees each change. */
  async function watchQueue(page: Page): Promise<void> {
    await page.evaluate(() => {
      const seen: string[] = [];
      (window as unknown as { __statuses: string[] }).__statuses = seen;
      const look = () => {
        for (const element of document.querySelectorAll('[data-testid="queue"] [data-status]')) {
          const status = element.getAttribute('data-status')!;
          if (!seen.includes(status)) seen.push(status);
        }
      };
      new MutationObserver(look).observe(document.body, { subtree: true, childList: true, attributes: true, attributeFilter: ['data-status'] });
      look();
    });
  }
  const statusesSeen = (page: Page) => page.evaluate(() => (window as unknown as { __statuses: string[] }).__statuses);

  test("real searches saved and opened in a new page: nothing is searched, the outputs are those before saving and the native CLI's, and the originals attach only when chosen and matching", async ({
    page,
  }) => {
    // The engine's workers (the engine and its threads) start with a search only.
    const engineWorkers: string[] = [];
    page.on('worker', (worker) => {
      if (/engine-worker|thread-worker/.test(worker.url())) engineWorkers.push(worker.url());
    });

    // Run 1: BLASTN of a query file against two subject files joined.
    await program(page, 'blastn');
    await openFiles(page, 'query', [{ name: 'reads.fa', text: READS }]);
    await openFiles(page, 'subject', [
      { name: 'part1.fa', text: PART1 },
      { name: 'part2.fa', text: PART2 },
    ]);
    await search(page, 1, 'Joined subjects');
    // Run 2: TBLASTN against the windows, the fifth left out.
    await clearInputs(page);
    await program(page, 'tblastn');
    await openFiles(page, 'query', [{ name: 'e2e_protein_query.faa', text: PROTEINS }]);
    await openFiles(page, 'subject', [{ name: 'e2e_amb_subject.fna', text: WINDOWS }]);
    await showRecords(page, 'subject');
    await page.getByTestId('subject-source-0-record-4').uncheck();
    await expect(page.getByTestId('subject-source-0-summary')).toContainText('(4 included)');
    await search(page, 2, 'Fifth window left out');
    expect(engineWorkers.length, 'the searches started the engine worker').toBeGreaterThan(0);

    // Candidates of both runs (the first query with hits of each), and a note.
    await chooseRun(page, 1);
    await addAll(page, /^[2-9] HSPs added to Candidates\.$/);
    await chooseRun(page, 2);
    await addAll(page, /^\d+ HSPs? added to Candidates\.$/);
    const candidates = (await page.getByTestId('tab-candidates-count').textContent())!;
    await page.getByTestId('tab-candidates').click();
    await page.getByTestId('candidate-note-1').fill(NOTE);
    const sequences = await downloaded(page, 'extract-download');
    expect(sequences.name).toBe('losat-candidates.fa');
    const aligned = await downloaded(page, 'extract-aligned');

    // Before saving: the outputs as exported, the runs' input FASTA, the native CLI's outputs, and LOSAT Web's own files.
    const before: Record<number, Record<Format, Buffer>> = {};
    const inputs: Record<number, RunInputs> = {};
    for (const run of [1, 2]) {
      before[run] = await exportedOutputs(page, run);
      inputs[run] = await runInputs(page, run);
      if (NATIVE === undefined) continue;
      const native = nativeOutputs(run, inputs[run]!);
      for (const format of FORMATS) expect(sha256(before[run]![format]), `run ${run} outfmt ${format} = the native CLI's`).toBe(native[format]);
    }
    if (NATIVE === undefined) test.info().annotations.push({ type: 'not compared with the native CLI', description: NO_NATIVE_REASON });
    expect(inputs[1]!.relations[1]).toContain('combined_subject.fa joins the 2 subject inputs (part1.fa, part2.fa) in the order chosen: 3 records.');
    expect(inputs[2]!.relations[1]).toContain(
      'e2e_amb_subject.fna has the 4 records that the run searched; the records left out of the chosen subject are not in it.',
    );
    const own = await ownFormats(page, 1);

    await expect(page.getByTestId('session-include-candidates')).toBeChecked();
    const session = await downloaded(page, 'session-save');
    await expect(message(page)).toHaveText(`Saved ${session.name}: 2 runs and ${candidates} candidates with their notes.`);

    // A new working session; from here on no run may pass through the queue, and no engine worker may start.
    await page.reload();
    await expect(queueRuns(page)).toHaveCount(0);
    const workersBefore = engineWorkers.length;
    await watchQueue(page);
    await openSession(page, session);
    await expect(message(page)).toHaveText(
      `Opened ${session.name}: 2 runs (runs 1 to 2 here) and ${candidates} candidates with their notes. Nothing was searched again.`,
    );
    await expect(queueRuns(page)).toHaveCount(2);
    for (const n of [1, 2]) {
      await expect(page.getByTestId(`run-${n}-status`)).toHaveText('completed');
      await expect(page.getByTestId(`run-${n}-origin`)).toHaveText(`From ${session.name}, run ${n} there (not searched again)`);
    }

    // The outputs of the loaded runs: the bytes exported before saving (so the native CLI's too);
    // the reproduction panel says where the inputs came from as it did before saving.
    for (const run of [1, 2]) {
      const after = await exportedOutputs(page, run);
      for (const format of FORMATS) expect(after[format].equals(before[run]![format]), `run ${run} outfmt ${format}: the bytes before saving`).toBe(true);
      await details(page, run);
      expect(await relations(page)).toEqual(inputs[run]!.relations);
    }
    // LOSAT Web's own files need no original FASTA: the same HSPs.
    const ownAfter = await ownFormats(page, 1);
    expect(ownAfter.hsps).toBe(own.hsps);
    expect(ownAfter.csv.bytes.equals(own.csv.bytes)).toBe(true);
    expect(JSON.parse(ownAfter.json.bytes.toString('utf8')).hsps).toEqual(JSON.parse(own.json.bytes.toString('utf8')).hsps);
    const reportText = (report: Download) => report.bytes.toString('utf8').replace(/\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(\.\d+)?(Z|[+-]\d\d:\d\d)?/g, '<time>');
    expect(reportText(ownAfter.report)).toBe(reportText(own.report));

    // What needs the original FASTA says so: the sequences of the tray and the run's input; the aligned rows do not need it.
    await page.getByTestId('tab-candidates').click();
    await expect(page.getByTestId('tab-candidates-count')).toHaveText(candidates);
    await expect(page.getByTestId('candidate-note-1')).toHaveValue(NOTE);
    const requirement = 'Runs 1 and 2 were loaded from a session file; choose their original subject FASTA in Run details to extract sequences.';
    await expect(page.getByTestId('extract-originals')).toHaveText(requirement);
    await expect(page.getByTestId('extract-download')).toBeDisabled();
    const alignedAfter = await downloaded(page, 'extract-aligned');
    expect(alignedAfter.bytes.equals(aligned.bytes)).toBe(true);
    await details(page, 1);
    await page.getByTestId('run-input-save-subject').click();
    await expect(page.getByTestId('run-files-message')).toHaveText(
      'combined_subject.fa was not saved: Run 1 was loaded from a session file; choose its original subject FASTA in Run details to save the input it searched.',
    );
    await expect(page.getByTestId('run-files-message')).toHaveAttribute('data-kind', 'error');

    // The same files in the search form attach nothing.
    await page.getByTestId('tab-search').click();
    await program(page, 'blastn');
    await openFiles(page, 'query', [{ name: 'reads.fa', text: READS }]);
    await openFiles(page, 'subject', [
      { name: 'part1.fa', text: PART1 },
      { name: 'part2.fa', text: PART2 },
    ]);
    await details(page, 1);
    for (const role of ['query', 'subject']) await expect(page.getByTestId(`run-original-${role}`)).toHaveAttribute('data-attached', 'false');
    await page.getByTestId('tab-candidates').click();
    await expect(page.getByTestId('extract-originals')).toHaveText(requirement);

    // A file with one residue changed is refused, and named; the right files attach.
    await details(page, 1);
    const attach = (testid: string, files: ReadonlyArray<{ name: string; text: string | Buffer }>) =>
      page.getByTestId(testid).setInputFiles(files.map((file) => ({ name: file.name, mimeType: 'text/plain', buffer: Buffer.from(file.text) })));
    await attach('run-attach-subject-files', [
      { name: 'part1.fa', text: PART1_CHANGED },
      { name: 'part2.fa', text: PART2 },
    ]);
    await expect(page.getByTestId('run-attach-subject-message')).toHaveText(
      'The subject FASTA was not attached to run 1: record 1 ("sA") differs from the saved run\'s record: its bytes have another SHA-256.',
    );
    await expect(page.getByTestId('run-original-subject')).toHaveAttribute('data-attached', 'false');
    await attach('run-attach-subject-files', [
      { name: 'part1.fa', text: PART1 },
      { name: 'part2.fa', text: PART2 },
    ]);
    await expect(page.getByTestId('run-original-subject')).toHaveAttribute('data-attached', 'true');
    await expect(page.getByTestId('run-original-subject-attached')).toContainText('Attached: part1.fa, part2.fa.');
    await expect(page.getByTestId('run-attach-subject-message')).toHaveCount(0);
    // The run's input FASTA is the one saved before: the attached files rebuilt with the run's records.
    const subjectAgain = await downloaded(page, 'run-input-save-subject');
    expect(subjectAgain.name).toBe('combined_subject.fa');
    expect(subjectAgain.bytes.equals(inputs[1]!.subject.bytes)).toBe(true);
    // Run 2's file, whose fifth record the session file says to leave out.
    await details(page, 2);
    await attach('run-attach-subject-files', [{ name: 'e2e_amb_subject.fna', text: WINDOWS }]);
    await expect(page.getByTestId('run-original-subject')).toHaveAttribute('data-attached', 'true');

    // The extraction gives the bytes that it gave before saving.
    await page.getByTestId('tab-candidates').click();
    await expect(page.getByTestId('extract-originals')).toHaveCount(0);
    const sequencesAfter = await downloaded(page, 'extract-download');
    expect(sequencesAfter.bytes.toString('utf8')).toBe(sequences.bytes.toString('utf8'));

    // Damaged session files are refused and change nothing.
    await openSession(page, { name: 'cut.losat-session.gz', bytes: session.bytes.subarray(0, Math.floor(session.bytes.length / 2)) });
    await expect(message(page)).toHaveAttribute('data-error', 'true');
    await expect(message(page)).toContainText('cut.losat-session.gz was not opened, and nothing was loaded.');
    const flipped = Buffer.from(session.bytes);
    flipped[Math.floor(flipped.length / 2)]! ^= 0x40;
    await openSession(page, { name: 'flipped.losat-session.gz', bytes: flipped });
    await expect(message(page)).toHaveAttribute('data-error', 'true');
    await expect(message(page)).toContainText('flipped.losat-session.gz was not opened, and nothing was loaded.');
    await expect(queueRuns(page)).toHaveCount(2);
    await expect(page.getByTestId('tab-candidates-count')).toHaveText(candidates);
    await details(page, 1);
    await expect(page.getByTestId('run-original-subject')).toHaveAttribute('data-attached', 'true');

    // Nothing was searched from the reload on.
    expect(await statusesSeen(page)).toEqual(['completed']);
    expect(engineWorkers.slice(workersBefore)).toEqual([]);
  });
});
