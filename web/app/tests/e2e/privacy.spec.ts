// V-PRIV (plan §5.10, §6.2; design §13.1, §13.3): research data never goes over the network. A
// unique canary of each test is put into everything that is research data - the file names, the
// deflines, the sequences (a nucleotide 40-mer, a protein motif), the job title, the candidates'
// notes, the names of the settings and session files the user opens, and inputs that are refused -
// and each path of the application is taken: a search and its results (Graphic Summary,
// Descriptions, Alignments, Dot Plot), every export, a reload, the saved session opened and its
// originals attached again; a cancelled search, a refused input, a refused option and a failed run.
// Every request is captured at the server and in the browser (support/privacy.ts), and none may go
// to another origin, send a body, ask for a path that is not a file of the site, carry a query
// string (the application uses none; no browser added one by itself) or carry a canary in its URL,
// a header or its body (as written, percent-encoded or in base64); no CSP violation may occur, and
// no page title or address may hold a canary, a query or a fragment. Downloads are made in the page:
// their contents hold the canaries, which shows that the data was there and stayed.
//
// The FakeEngine build and the engine build (CI's browser-engine job) run the same flows; the engine
// build searches for real, and only it can refuse an option or fail a run.
//
// The browsers differ in what they show Playwright (S16), which is why the server logs too: Chromium
// reports every request (the nested thread workers' scripts and the CSP reports included); Firefox
// leaves out the thread workers' scripts, its own /favicon.ico and the CSP reports; WebKit reports
// the page's blob: reads as requests, and sends no CSP report for a violation in a worker (the CSP
// still blocks it: the capture's own test shows that nothing reaches the other origin).
//
// Requests to another origin never reach the capture's server, so for them the proof is the CSP
// (`connect-src 'self'` on every response of the site, workers' scripts included). The capture's own
// test plants an attempt in the page, in a dedicated worker, in a worker that a worker starts and in
// a service worker (code review 2, L3): in all three browsers each fetch fails and the other origin,
// which logs every request, receives nothing. Chromium and Firefox also send a CSP report from each
// worker; WebKit sends none from any worker, so there the evidence is the failed fetch with the
// empty log: the request never left (one that CORS refused would still reach the other origin).
//
// TD-11 (chromium project only: it reads files): no file that the site serves holds a path of the
// build machine (the checkout, CARGO_HOME, RUSTUP_HOME, HOME) or the user name.
//
// The delivery paths (S16): an update while the tab holds research data (the notice, the download
// and check, the switch; the service worker's caches hold only the releases' files and the shell),
// "Try an example" and the self-test (they fetch release files only), and measurement (nothing
// under /api/ before consent or after "Don't allow"; after "Allow" only the fixed POST /api/visit).
import { execSync } from 'node:child_process';
import { mkdtempSync, realpathSync, rmSync } from 'node:fs';
import { readFile } from 'node:fs/promises';
import { homedir, tmpdir, userInfo } from 'node:os';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { expect, type Page } from '@playwright/test';
import { readSiteIndex } from '../../build/site-files';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { DIST, PrivacyCapture, siteFiles, startCaptureSite, test } from './support/privacy';
import {
  canaryForms,
  filesUnder,
  findCanaries,
  machinePatterns,
  makeCanary,
  scanFiles,
  userNameInPathsOnly,
  type Canary,
  type Exchange,
  type MachineContext,
} from './support/privacy-scan';
import { fasta, openFiles, openParameters, paste, program, settled, showOutput, submit, waitStatus } from './support/search';
import { assemble, buildSite, changeServiceWorker, expectedCaches } from './support/sites';

test.setTimeout(BUILD_HAS_ENGINE ? 300_000 : 180_000);

/** Pseudo-random letters (a linear congruential generator): the same in every run. */
function letters(alphabet: string, seed: number, length: number): string {
  let state = seed;
  let out = '';
  for (let i = 0; i < length; i++) {
    state = (Math.imul(state, 1103515245) + 12345) >>> 0;
    out += alphabet[(state >>> 16) % alphabet.length];
  }
  return out;
}
const dna = (seed: number, length: number) => letters('ACGT', seed, length);
const protein = (seed: number, length: number) => letters('ACDEFGHIKLMNPQRSTVWY', seed, length);
/** A FASTA record with 80 letters a line, so that the canary's 40-mer is not split (it starts a line). */
const record = (title: string, sequence: string) => `>${title}\n${sequence.match(/.{1,80}/g)!.join('\n')}\n`;

/** The name and bytes of the file that clicking a button saves (the browser's download, made in the page). */
async function saved(page: Page, testid: string): Promise<{ name: string; bytes: Buffer }> {
  const download = page.waitForEvent('download');
  await page.getByTestId(testid).click();
  const file = await download;
  return { name: file.suggestedFilename(), bytes: await readFile((await file.path())!) };
}

async function clearInputs(page: Page): Promise<void> {
  for (const role of ['query', 'subject'] as const) {
    const sources = page.locator(`[data-testid^="${role}-source-"][data-status]`);
    while ((await sources.count()) > 0) await page.getByTestId(`${role}-source-0-remove`).click();
  }
}

/**
 * Waits until the application's start-up has finished: the storage cleanup, and the service
 * worker's storing of the release ("Offline ready"), whose progress line in the header would move
 * the page under a click (WebKit then fires no click).
 */
async function started(page: Page): Promise<void> {
  await expect(page.getByTestId('storage-status')).toHaveAttribute('data-cleanup', 'done');
  await expect(page.getByTestId('offline-status')).toHaveAttribute('data-state', 'ready', { timeout: 60_000 });
}

/** Opens the application from the capture's server and waits until its start-up has finished. */
async function openApp(page: Page, capture: PrivacyCapture): Promise<void> {
  await page.goto(`${capture.origin}/`);
  await started(page);
}

/**
 * The capture saw the site's own requests: the boot page, the release's script and the Data
 * worker's (and the engine's), and the service worker.
 */
function expectSiteRequests(capture: PrivacyCapture): void {
  const paths = capture.serverPaths();
  for (const path of ['/', '/boot.js', '/sw.js']) expect(paths).toContain(path);
  const release = (pattern: string) => expect(paths.some((path) => new RegExp(`^/r/[0-9a-z.-]+/${pattern}$`).test(path)), pattern).toBe(true);
  release('assets/index-[\\w-]+\\.js');
  release('assets/data-worker-[\\w-]+\\.js');
  if (BUILD_HAS_ENGINE) {
    release('assets/engine-worker-[\\w-]+\\.js');
    release('engine/losat-web-(serial|threads)-[0-9a-f]+\\.wasm');
  }
}

/** The research data of the main flow, with the canary in every part. */
function searchInputs(canary: Canary) {
  const t = canary.token;
  const shared = dna(81, 40);
  return {
    query: {
      name: `${t}-query.fa`,
      text: record(`${t}_q1 ${t} first query`, canary.nucleotides + shared) + record(`${t}_q2 second query`, dna(82, 50)),
    },
    subject: {
      name: `${t}-subject.fa`,
      // The 40-mer starts the subject's hit region: in the engine build's HSP, and in the FakeEngine's (the first half).
      text: record(`${t}_s1 ${t} first subject`, canary.nucleotides + shared + dna(83, 30)) + record(`${t}_s2`, dna(84, 60)),
    },
    subjectId: `${t}_s1`,
    title: `${t} job title`,
    note: `${t} candidate note`,
  };
}

test('a search, its results and every export, a reload, the saved session and its originals send no research data', async ({
  page,
  capture,
}) => {
  const canary = makeCanary();
  const t = canary.token;
  const inputs = searchInputs(canary);
  await openApp(page, capture);

  // --- a search ---
  await program(page, 'blastn');
  await openFiles(page, 'query', [inputs.query]);
  await openFiles(page, 'subject', [inputs.subject]);
  await page.getByTestId('job-title').fill(inputs.title);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await capture.checkpoint(page);

  // --- its results: the Graphic Summary, the Descriptions and the Alignments on one page, the Dot Plot ---
  await page.getByTestId('tab-results').click();
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
  await expect(page.getByTestId('graphic-canvas')).toHaveAttribute('data-drawn', /^[1-9]/);
  await expect(page.getByTestId('subject-list')).toContainText(`${t}_s1`);
  await page.getByTestId('results-alignments').scrollIntoViewIfNeeded();
  await expect(page.getByTestId('hsp-detail')).toHaveAttribute('data-state', 'ready');
  await page.getByTestId('pane-dotplot').click();
  await expect(page.getByTestId('dotplot-canvas')).toHaveAttribute('data-segments', /^[1-9]/);
  await capture.checkpoint(page);

  // --- every export: each is made in the page and saved by the browser ---
  const svg = await saved(page, 'dotplot-svg');
  expect(svg.name).toMatch(/^losat-run1-dotplot-q1-s1\.svg$/);
  const outputs: Buffer[] = [];
  for (const format of [0, 6, 7] as const) {
    await showOutput(page, 1, format, false);
    outputs.push((await saved(page, 'export-output')).bytes);
  }
  expect(outputs[1]!.toString()).toContain(`${t}_q1`);
  const own: Record<string, Buffer> = {};
  for (const format of ['csv', 'json', 'report'] as const) {
    own[format] = (await saved(page, `export-${format}`)).bytes;
    await expect(page.getByTestId('export-summary')).toHaveAttribute('data-format', format);
  }
  expect(own['csv']!.toString()).toContain(`${t}_s1`);
  expect(own['report']!.toString()).toContain(inputs.title);

  // The run's input FASTA and its settings (Run details).
  await page.getByTestId('results-view-details').click();
  const runQuery = await saved(page, 'run-input-save-query');
  expect(runQuery.bytes.toString()).toContain(canary.nucleotides);
  await saved(page, 'run-input-save-subject');
  await saved(page, 'run-settings-save');

  // Candidates with a note: the sequences and the aligned rows.
  await page.getByTestId('results-view-hits').click();
  await page.getByTestId('descriptions-select-all').check();
  await page.getByTestId('descriptions-add-candidates').click();
  await expect(page.getByTestId('candidates-added')).toHaveText(/^\d+ HSPs? added to Candidates\.$/);
  await page.getByTestId('tab-candidates').click();
  await page.getByTestId('candidate-note-1').fill(inputs.note);
  const extracted = await saved(page, 'extract-download');
  expect(extracted.bytes.toString()).toContain(canary.nucleotides);
  await saved(page, 'extract-aligned');
  await capture.checkpoint(page);

  // The search form's settings, saved and opened again under a name with the canary; a broken one refused.
  await page.getByTestId('tab-search').click();
  const settings = await saved(page, 'settings-save');
  await page.getByTestId('settings-load').setInputFiles({ name: `${t}-settings.json`, mimeType: 'application/json', buffer: settings.bytes });
  await expect(page.getByTestId('settings-message')).toContainText(`Loaded ${t}-settings.json`);
  await page
    .getByTestId('settings-load')
    .setInputFiles({ name: `${t}-broken.json`, mimeType: 'application/json', buffer: Buffer.from(`{"format": "${t}`) });
  await expect(page.getByTestId('settings-message')).toHaveAttribute('data-kind', 'error');
  await expect(page.getByTestId('settings-message')).toContainText(`${t}-broken.json was not loaded`);

  // The session file: the runs, the candidates and their notes.
  const session = await saved(page, 'session-save');
  expect(gunzipSync(session.bytes).toString('latin1')).toContain(inputs.note);
  await capture.checkpoint(page);

  // --- a reload: a new working session; the session file opened under a name with the canary ---
  await page.reload();
  await started(page);
  await expect(page.getByTestId('queue').locator('li[data-status]')).toHaveCount(0);
  const sessionName = `${t}.losat-session.gz`;
  await page.getByTestId('session-file').setInputFiles({ name: sessionName, mimeType: 'application/gzip', buffer: session.bytes });
  await expect(page.getByTestId('session-busy')).toHaveCount(0, { timeout: 30_000 });
  await expect(page.getByTestId('session-message')).toContainText(`Opened ${sessionName}: 1 run`);
  await expect(page.getByTestId('session-message')).toHaveAttribute('data-error', 'false');

  // The originals attached again, under other names with the canary; then the sequences are extracted as before.
  await page.getByTestId('run-1-open').click();
  await page.getByTestId('results-view-details').click();
  for (const role of ['query', 'subject'] as const) {
    await page
      .getByTestId(`run-attach-${role}-files`)
      .setInputFiles({ name: `${t}-${role}-original.fa`, mimeType: 'text/plain', buffer: Buffer.from(inputs[role].text) });
    await expect(page.getByTestId(`run-original-${role}`)).toHaveAttribute('data-attached', 'true');
  }
  await page.getByTestId('tab-candidates').click();
  await expect(page.getByTestId('candidate-note-1')).toHaveValue(inputs.note);
  await expect(page.getByTestId('extract-originals')).toHaveCount(0);
  expect((await saved(page, 'extract-download')).bytes.equals(extracted.bytes)).toBe(true);
  await capture.checkpoint(page);

  expect(await capture.findings(canary)).toEqual([]);
  expectSiteRequests(capture);
});

test('a cancelled search, a refused input, a refused option and a failed run send no research data', async ({ page, capture }) => {
  const canary = makeCanary();
  const t = canary.token;
  const motif = record(`${t}_p1 ${t} protein`, protein(91, 30) + canary.protein + protein(92, 30));
  await openApp(page, capture);

  // --- a search cancelled while it runs ---
  await program(page, 'blastp');
  await page.getByTestId('job-title').fill(`${t} cancelled`);
  if (BUILD_HAS_ENGINE) {
    // A proteome against itself (about 500,000 residues, as search.spec.ts's slow search), the motif's record first.
    const proteome = { name: `${t}-proteome.faa`, text: motif + fasta('NZ_CP006932.faa').toString('latin1') };
    await openFiles(page, 'query', [proteome]);
    await openFiles(page, 'subject', [proteome]);
    await submit(page);
    await expect(page.getByTestId('run-1-status')).toHaveText('running', { timeout: 60_000 });
    await page.getByTestId('run-1-cancel').click();
  } else {
    // The FakeEngine's run takes three phases of 50 ms: Cancel is pressed in the page as soon as it runs.
    await openFiles(page, 'query', [{ name: `${t}-motif.faa`, text: motif }]);
    await openFiles(page, 'subject', [{ name: `${t}-motif-subject.faa`, text: motif }]);
    await page.evaluate(() => {
      const observer = new MutationObserver(() => {
        if (document.querySelector('[data-testid="run-1-status"]')?.textContent !== 'running') return;
        observer.disconnect();
        document.querySelector<HTMLButtonElement>('[data-testid="run-1-cancel"]')?.click();
      });
      observer.observe(document.body, { subtree: true, childList: true, characterData: true });
    });
    await submit(page);
  }
  await waitStatus(page, 1, 'cancelled');
  // The next search runs (in the engine build, on the engine started again after the cancel).
  await clearInputs(page);
  await openFiles(page, 'query', [{ name: `${t}-motif.faa`, text: motif }]);
  await openFiles(page, 'subject', [{ name: `${t}-motif-subject.faa`, text: record(`${t}_p2`, protein(93, 40) + canary.protein) }]);
  await page.getByTestId('job-title').fill(`${t} after the cancel`);
  await submit(page);
  await waitStatus(page, 2, 'completed');
  await capture.checkpoint(page);

  // --- an input that is refused: the canary in its name and its content ---
  await program(page, 'blastn');
  await clearInputs(page);
  const refused = BUILD_HAS_ENGINE ? `>${t}_bad ${t}\n${canary.nucleotides}\n>?10\nACGT\n` : `>${t}_bad ${t}\n${canary.nucleotides}!ACGT\n`;
  await openFiles(page, 'query', [{ name: `${t}-refused.fa`, text: refused }]);
  if (BUILD_HAS_ENGINE) await expect(page.getByTestId('query-source-0-error')).toContainText("is a gap line ('>?')");
  else await expect(page.getByTestId('query-source-0-check')).toHaveAttribute('data-check', 'refused');
  await paste(page, 'subject', record(`${t}_s1`, canary.nucleotides));
  await submit(page);
  await expect(page.getByTestId('search-message')).toContainText(`${t}-refused.fa`);
  await expect(page.getByTestId('run-3')).toHaveCount(0);
  await capture.checkpoint(page);

  if (BUILD_HAS_ENGINE) {
    // --- an option that the engine refuses: the canary as its value ---
    await program(page, 'tblastx');
    await clearInputs(page);
    await openParameters(page);
    await page.getByTestId('param-seg').fill(t);
    await expect(page.getByTestId('argv-validation')).toContainText('The engine refuses these options');
    await paste(page, 'query', record(`${t}_q1 ${t}`, dna(7, 130) + canary.nucleotides + dna(8, 130)));
    await paste(page, 'subject', record(`${t}_s1 ${t}`, dna(9, 300)));
    await submit(page);
    await expect(page.getByTestId('search-message')).toBeVisible();
    await expect(page.getByTestId('run-3')).toHaveCount(0);
    await page.getByTestId('param-seg').fill('');
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');

    // --- a run that fails: a TBLASTX subject record without residues (results.spec.ts's failing run) ---
    await clearInputs(page);
    await paste(page, 'query', record(`${t}_q1 ${t}`, dna(7, 130) + canary.nucleotides + dna(8, 130)));
    await paste(page, 'subject', `>${t}_empty ${t}\n`);
    await settled(page, 'subject');
    await page.getByTestId('job-title').fill(`${t} failed`);
    await submit(page);
    await waitStatus(page, 3, 'failed');
    await capture.checkpoint(page);
  }

  expect(await capture.findings(canary)).toEqual([]);
  expectSiteRequests(capture);
});

// --- delivery: update, onboarding, measurement ----------------------------------------------------

/** The work of a tab: a BLASTN search of the canary's inputs as run `run`, its HSPs in Candidates with a note. */
async function researchWork(page: Page, inputs: ReturnType<typeof searchInputs>, run: number): Promise<void> {
  await page.getByTestId('tab-search').click();
  await program(page, 'blastn');
  await clearInputs(page);
  await openFiles(page, 'query', [inputs.query]);
  await openFiles(page, 'subject', [inputs.subject]);
  await page.getByTestId('job-title').fill(inputs.title);
  await submit(page);
  await waitStatus(page, run, 'completed');
  await page.getByTestId(`run-${run}-open`).click();
  await page.getByTestId('results-view-hits').click();
  await expect(page.getByTestId('subject-list')).toContainText(inputs.subjectId);
  await page.getByTestId('descriptions-select-all').check();
  await page.getByTestId('descriptions-add-candidates').click();
  await expect(page.getByTestId('candidates-added')).toHaveText(/^\d+ HSPs? added to Candidates\.$/);
  await page.getByTestId('tab-candidates').click();
  await page.getByTestId('candidate-note-1').fill(inputs.note);
}

/** Every cache of the origin, read from the page: by name, each entry's URL (path and query) and the SHA-256 of its body. */
function storedFiles(page: Page): Promise<Record<string, Record<string, string>>> {
  return page.evaluate(async () => {
    const hex = (bytes: ArrayBuffer) => Array.from(new Uint8Array(bytes), (byte) => byte.toString(16).padStart(2, '0')).join('');
    const stored: Record<string, Record<string, string>> = {};
    for (const name of await caches.keys()) {
      const cache = await caches.open(name);
      const files: Record<string, string> = {};
      for (const request of await cache.keys()) {
        const url = new URL(request.url);
        const response = await cache.match(request);
        files[url.pathname + url.search] = hex(await crypto.subtle.digest('SHA-256', await response!.arrayBuffer()));
      }
      stored[name] = files;
    }
    return stored;
  });
}

/**
 * What the caches hold that is not a file of a release (its manifest's files and the manifest) or of
 * the shell, byte for byte (support/sites.ts `expectedCaches`): each such entry, and a canary in a
 * cache's name or an entry's URL. Bytes equal to a release's or the shell's hold no canary.
 */
function cacheFindings(stored: Record<string, Record<string, string>>, releases: ReadonlyArray<readonly [site: string, release: string]>, canary: Canary): string[] {
  const allowed = new Map<string, Array<Record<string, string>>>();
  for (const [site, release] of releases) {
    for (const [name, files] of Object.entries(expectedCaches(site, release))) allowed.set(name, [...(allowed.get(name) ?? []), files]);
  }
  const forms = canaryForms(canary);
  const found: string[] = [];
  for (const [name, files] of Object.entries(stored)) {
    const expected = allowed.get(name);
    if (expected === undefined || findCanaries(name, forms).length > 0) {
      found.push(`a cache of neither a release nor the shell: ${name}`);
      continue;
    }
    for (const [url, sha256] of Object.entries(files)) {
      if (findCanaries(url, forms).length > 0) found.push(`a canary in the URL of a cached entry: ${name} ${url}`);
      else if (!expected.some((files) => files[url] === sha256)) found.push(`${name} holds ${url}, which is not that file of the release or the shell`);
    }
  }
  return found;
}

test('an update while the tab holds research data: notified, downloaded and checked, switched; the caches hold only release files and the shell', async ({
  page,
  capture,
}) => {
  test.setTimeout(BUILD_HAS_ENGINE ? 420_000 : 240_000);
  const canary = makeCanary();
  const t = canary.token;
  const inputs = searchInputs(canary);
  // Release A is the build under test (dist/); B is built like it, and its deployment keeps A.
  const a = readSiteIndex(DIST).latest;
  const b = 'vpriv-b';
  const work = mkdtempSync(join(tmpdir(), 'losat-vpriv-update-'));
  try {
    const ab = join(work, 'ab');
    buildSite(ab, b);
    assemble(ab, [DIST]);
    changeServiceWorker(ab, 'the V-PRIV deployment of vpriv-b');
    const releases = [[DIST, a], [ab, a], [ab, b]] as const;

    await openApp(page, capture);
    await expect(page.getByTestId('release-id')).toHaveText(a);
    expect(cacheFindings(await storedFiles(page), releases, canary)).toEqual([]);

    // B is deployed; the tab is notified, and its work holds the canaries while it downloads and switches.
    capture.server.deploy(ab);
    await page.reload();
    await started(page);
    await expect(page.getByTestId('release-id')).toHaveText(a);
    await expect(page.getByTestId('update-notice')).toHaveAttribute('data-state', 'available', { timeout: 30_000 });
    await researchWork(page, inputs, 1);
    const session = await saved(page, 'session-save');
    await capture.checkpoint(page);
    await page.getByTestId('update-review').click();
    await expect(page.getByTestId('update-impact-saved')).toContainText('1 completed run');
    await page.getByTestId('update-download').click();
    await expect(page.getByTestId('update-downloaded')).toBeVisible({ timeout: 60_000 });
    expect(cacheFindings(await storedFiles(page), releases, canary)).toEqual([]);
    expect(Object.keys(await storedFiles(page))).toContain(`losat-web-release-${b}`);
    await capture.checkpoint(page);
    await page.getByTestId('update-switch').click();
    await expect(page.getByTestId('release-id')).toHaveText(b, { timeout: 60_000 });
    await started(page);

    // The work saved before the switch, opened in B under a name with the canary.
    const sessionName = `${t}-before-the-update.losat-session.gz`;
    await page.getByTestId('session-file').setInputFiles({ name: sessionName, mimeType: 'application/gzip', buffer: session.bytes });
    await expect(page.getByTestId('session-busy')).toHaveCount(0, { timeout: 30_000 });
    await expect(page.getByTestId('session-message')).toContainText(`Opened ${sessionName}: 1 run`);
    expect(cacheFindings(await storedFiles(page), releases, canary)).toEqual([]);
    await capture.checkpoint(page);

    expect(await capture.findings(canary)).toEqual([]);
    expect(capture.serverPaths().some((path) => path.startsWith(`/r/${b}/`))).toBe(true);
  } finally {
    rmSync(work, { recursive: true, force: true });
  }
});

/** Chooses a program and puts its example in the form, answering the question if inputs would be replaced. */
async function loadExample(page: Page, id: string): Promise<void> {
  await program(page, id);
  await page.getByTestId('try-example').click();
  const loaded = page.getByTestId('example-loaded').filter({ hasText: `example-${id}-query.fasta` });
  const confirm = page.getByTestId('try-example-confirm');
  await expect(loaded.or(confirm)).toBeVisible();
  if (await confirm.isVisible()) await confirm.click();
  await expect(loaded).toBeVisible();
  await settled(page, 'query');
  await settled(page, 'subject');
}

test('"Try an example" and the self-test fetch only files of the release while the form holds research data', async ({ page, capture }) => {
  const canary = makeCanary();
  const t = canary.token;
  await openApp(page, capture);
  const release = (await page.getByTestId('release-id').textContent())!.trim();
  await program(page, 'blastn');
  await paste(page, 'query', record(`${t}_q1 ${t}`, canary.nucleotides));
  await page.getByTestId('job-title').fill(`${t} before the examples`);

  const mark = await capture.mark();
  for (const id of ['blastn', 'blastp', 'tblastn', 'tblastx']) await loadExample(page, id);
  if (BUILD_HAS_ENGINE) {
    await page.getByTestId('self-test').locator('summary').click();
    await page.getByTestId('self-test-run').click();
    await expect(page.getByTestId('self-test-summary')).toHaveText(/^All \d+ cases passed/, { timeout: 180_000 });
  }
  await capture.checkpoint(page);
  const fetched = (await capture.since(mark))
    .filter((exchange) => !/^(blob|data):/.test(exchange.url))
    .map((exchange) => `${exchange.method} ${new URL(exchange.url).pathname}${new URL(exchange.url).search}`);
  expect(fetched.some((line) => line.startsWith(`GET /r/${release}/assets/example-blastn-query-`))).toBe(true);
  // Besides the release's files, only the service worker's script, which the browser checks for an
  // update by itself (Firefox does during the self-test), may be fetched.
  expect(fetched.filter((line) => !line.startsWith(`GET /r/${release}/`) && line !== 'GET /sw.js')).toEqual([]);
  expect(await capture.findings(canary)).toEqual([]);
});

const VISIT_ENDPOINT = '/api/visit';
const CLIENT_ID = /^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$/;

/** Whether a request is the one event of measurement: POST /api/visit with exactly {"v":1,"release":..,"client":..}. */
function isVisit(exchange: Exchange, release: string, client: string): boolean {
  const url = new URL(exchange.url);
  return (
    exchange.method === 'POST' &&
    url.pathname === VISIT_ENDPOINT &&
    url.search === '' &&
    Buffer.from(exchange.body).toString('utf8') === JSON.stringify({ v: 1, release, client })
  );
}

test.describe('measurement', () => {
  // A site whose releases.json names the endpoint, as one built with LOSAT_WEB_ANALYTICS_ENDPOINT=/api/visit.
  test.use({ captureSite: { analyticsEndpoint: VISIT_ENDPOINT } });

  test('nothing under /api/ before consent or after "Don\'t allow"; after "Allow" only the fixed POST /api/visit, while the work holds research data', async ({
    page,
    capture,
  }) => {
    const canary = makeCanary();
    const inputs = searchInputs(canary);
    const release = readSiteIndex(DIST).latest;
    const api = () => capture.server.requests.filter((exchange) => new URL(exchange.url).pathname.startsWith('/api/'));

    // Before consent: the question, and the work of the tab.
    await openApp(page, capture);
    await expect(page.getByTestId('analytics-consent')).toBeVisible();
    await researchWork(page, inputs, 1);
    await saved(page, 'session-save');
    await capture.checkpoint(page);
    expect(api()).toEqual([]);

    // "Don't allow", a reload, and the work again.
    await page.getByTestId('analytics-deny').click();
    await expect(page.getByTestId('analytics-setting')).toHaveText('Usage statistics: off');
    await page.reload();
    await started(page);
    await expect(page.getByTestId('analytics-consent')).toHaveCount(0);
    await researchWork(page, inputs, 1);
    await capture.checkpoint(page);
    expect(api()).toEqual([]);

    // "Allow": one event for this page load; the work of the tab; a reload, one event more.
    await page.getByTestId('analytics-change').click();
    await page.getByTestId('analytics-allow').click();
    await expect(page.getByTestId('analytics-setting')).toHaveText('Usage statistics: on');
    await expect.poll(() => api().length).toBe(1);
    await researchWork(page, inputs, 2);
    await saved(page, 'session-save');
    await page.reload();
    await started(page);
    await expect.poll(() => api().length).toBe(2);
    await researchWork(page, inputs, 1);
    await capture.checkpoint(page);

    const visits = api();
    expect(visits).toHaveLength(2);
    const client = (JSON.parse(Buffer.from(visits[0]!.body).toString('utf8')) as { client: string }).client;
    expect(client).toMatch(CLIENT_ID);
    for (const visit of visits) {
      expect(isVisit(visit, release, client)).toBe(true);
      expect(visit.headers['content-type']).toBe('application/json');
      expect(visit.headers['cookie']).toBeUndefined();
      expect(visit.headers['referer']).toBeUndefined();
    }
    expect(await capture.findings(canary, { allowBody: (exchange) => isVisit(exchange, release, client) })).toEqual([]);
  });
});

/**
 * The planted sources of the capture's own test that try a request to another origin from outside
 * the page (code review 2, L3): a dedicated worker, a worker that a worker starts, and a service
 * worker (at its own scope, in its install event). Each says whether its fetch went through.
 */
const PLANTED_SOURCES = ['worker', 'nested-worker', 'service-worker'] as const;
type PlantedSource = (typeof PLANTED_SOURCES)[number];
const PLANTED_FILES = new Map([
  ['__privacy/planted-worker.js', "fetch(new URL(location.href).searchParams.get('to')).then(() => postMessage('sent'), () => postMessage('failed'));"],
  [
    '__privacy/planted-nest.js',
    "const nested = new Worker('/__privacy/planted-worker.js' + location.search); nested.onmessage = (event) => postMessage(event.data); nested.onerror = () => postMessage('the nested worker did not start');",
  ],
  [
    '__privacy/planted-sw.js',
    `self.addEventListener('install', (event) => event.waitUntil(
       fetch(new URL(location.href).searchParams.get('to')).then(() => 'sent', () => 'failed').then(async (result) => {
         for (const client of await self.clients.matchAll({ includeUncontrolled: true, type: 'window' })) client.postMessage({ planted: result });
       })));`,
  ],
]);

test('the capture finds what it looks for: a canary in a URL or a header (also in base64), a query string, a path that is not the site\'s, CSP violations in a page and in workers; nothing reaches another origin from a page, a worker, a nested worker or a service worker', async ({
  context,
  browserName,
}, testInfo) => {
  // Another origin that logs what reaches it, and the planted workers that only this test serves:
  // the page and each worker try a request to that origin, which the CSP blocks.
  const foreign = await startCaptureSite(DIST);
  const server = await startCaptureSite(DIST, { extra: PLANTED_FILES });
  const capture = new PrivacyCapture(server);
  try {
    await capture.attach(context);
    const page = await context.newPage();
    await openApp(page, capture);
    const canary = makeCanary();
    const outcomes = await page.evaluate(
      async ({ token, to }) => {
        await fetch(`/__privacy/planted?note=${token}`);
        // The canary in base64, in a query and in a header (code review 2, L2).
        const encoded = btoa(`>${token} a record`);
        await fetch(`/__privacy/planted-b64?d=${encoded}`, { headers: { 'x-d': encoded } });
        await fetch(`${to}/${token}`).catch(() => undefined);
        const target = (source: string) => `?to=${encodeURIComponent(`${to}/from-a-${source}`)}`;
        const answer = (worker: Worker) =>
          new Promise<string>((resolve) => {
            worker.onmessage = (event: MessageEvent<unknown>) => resolve(String(event.data));
            worker.onerror = () => resolve('the worker did not start');
          });
        const worker = await answer(new Worker(`/__privacy/planted-worker.js${target('worker')}`));
        const nested = await answer(new Worker(`/__privacy/planted-nest.js${target('nested-worker')}`));
        const fromServiceWorker = new Promise<string>((resolve) => {
          navigator.serviceWorker.addEventListener('message', (event: MessageEvent<{ planted?: string } | null>) => {
            if (typeof event.data?.planted === 'string') resolve(event.data.planted);
          });
          setTimeout(() => resolve('no answer within 20 s'), 20_000);
        });
        navigator.serviceWorker.startMessages();
        await navigator.serviceWorker.register(`/__privacy/planted-sw.js${target('service-worker')}`, { scope: '/__privacy/' });
        return { worker, 'nested-worker': nested, 'service-worker': await fromServiceWorker };
      },
      { token: canary.token, to: foreign.origin },
    );
    // Each planted fetch failed: the CSP blocked it (the other origin would have answered).
    expect(outcomes).toEqual({ worker: 'failed', 'nested-worker': 'failed', 'service-worker': 'failed' });

    const reports = async () => (await capture.findings(canary)).filter((line) => line.startsWith('a CSP violation report'));
    const reported = async (source: PlantedSource) => (await reports()).some((line) => line.includes(` blocked ${foreign.origin}/from-a-${source} in `));
    // Browsers send CSP reports in their own time. Chromium and Firefox send one from each planted
    // worker; WebKit sends none from any worker (S16: WebKit 26.x), where the evidence is the
    // worker's failed fetch and the other origin's empty log below.
    if (browserName !== 'webkit') await expect.poll(async () => (await Promise.all(PLANTED_SOURCES.map(reported))).every(Boolean), { timeout: 15_000 }).toBe(true);
    const evidence = Object.fromEntries(await Promise.all(PLANTED_SOURCES.map(async (source) => [source, { fetch: outcomes[source], cspReport: await reported(source) }] as const)));
    testInfo.annotations.push({ type: 'planted workers', description: JSON.stringify(evidence) });
    await expect.poll(async () => (await reports()).length, { timeout: 15_000 }).toBeGreaterThanOrEqual(1);
    testInfo.annotations.push({ type: 'CSP reports', description: (await reports()).join('\n') });
    expect(await capture.findings(canary)).toEqual(
      expect.arrayContaining([
        expect.stringMatching(/^canary \(token\) in the URL: server GET http:\/\/127\.0\.0\.1:\d+\/__privacy\/planted\?note=LOSATPRIV/),
        expect.stringMatching(/^canary \(token\) in the URL: browser GET http:\/\/127\.0\.0\.1:\d+\/__privacy\/planted\?note=LOSATPRIV/),
        expect.stringMatching(/^a request for a path that is not a file of the site: server GET .*\/__privacy\/planted\?note=/),
        expect.stringMatching(/^canary \(token, base64\) in the URL: server GET .*\/__privacy\/planted-b64\?d=/),
        expect.stringMatching(/^canary \(token, base64\) in a header: server GET .*\/__privacy\/planted-b64\?d=/),
        expect.stringMatching(/^a query string on a request to the site: server GET .*\/__privacy\/planted-b64\?d=/),
        expect.stringMatching(new RegExp(`^a CSP violation in a page: connect-src blocked ${foreign.origin}/LOSATPRIV`)),
        expect.stringMatching(new RegExp(`^a CSP violation report: connect-src blocked ${foreign.origin}/LOSATPRIV\\w+ in ${server.origin}/$`)),
      ]),
    );
    // The CSP held: nothing reached the other origin, from the page or from any worker.
    expect(foreign.requests).toEqual([]);
  } finally {
    await capture.report(testInfo);
    await server.close();
    await foreign.close();
  }
});

// --- TD-11 ---------------------------------------------------------------------------------------------

const APP = fileURLToPath(new URL('../../', import.meta.url));

/** The build machine as this process sees it: in CI, /home/runner/work/LOSAT/LOSAT and the user `runner`. */
function machineContext(): MachineContext {
  const git = (args: string) => {
    try {
      return execSync(`git ${args}`, { cwd: APP, encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] }).trim();
    } catch {
      return '';
    }
  };
  const real = (path: string) => {
    try {
      return realpathSync(path);
    } catch {
      return path;
    }
  };
  const top = git('rev-parse --show-toplevel');
  const common = git('rev-parse --git-common-dir');
  const checkouts = [REPOSITORY, real(REPOSITORY)];
  if (top !== '') checkouts.push(top, real(top));
  // A worktree's main checkout: the folder of the shared .git.
  if (common !== '') checkouts.push(dirname(real(resolve(APP, common))));
  const home = process.env['HOME'] || homedir();
  return {
    user: userInfo().username,
    home,
    cargoHome: process.env['CARGO_HOME'] || join(home, '.cargo'),
    rustupHome: process.env['RUSTUP_HOME'] || join(home, '.rustup'),
    checkouts: [...new Set(checkouts.map((path) => path.replace(/[/\\]+$/, '')))],
  };
}

test('TD-11: no file that the site serves holds a path of the build machine or the user name', async ({ browserName }, testInfo) => {
  test.skip(browserName !== 'chromium', 'reads the built files: once, not in every browser');
  const machine = machineContext();
  // The files that the site lists (the shell, the other site-level files, each release's manifest
  // and its files) are every file of the build but the `_headers` that Pages reads.
  const files = siteFiles(DIST).map(({ file }) => file);
  expect([...files].sort()).toEqual(filesUnder(DIST).filter((file) => file !== '_headers'));
  for (const file of ['index.html', 'boot.js', 'sw.js', 'releases.json']) expect(files).toContain(file);
  expect(files.filter((file) => file.endsWith('.js')).length).toBeGreaterThan(4);
  if (BUILD_HAS_ENGINE) expect(files.filter((file) => file.endsWith('.wasm'))).toHaveLength(2);
  const found = scanFiles(DIST, files, machinePatterns(machine));
  // A common user name (CI's `runner`) counts only in a path; how often it occurs as a bare word is recorded.
  if (userNameInPathsOnly(machine.user)) {
    const word = new RegExp(`(?<![A-Za-z0-9_])${machine.user.replace(/[.*+?^${}()|[\]\\]/g, '\\$&')}(?![A-Za-z0-9_])`, 'g');
    const bare = scanFiles(DIST, files, [{ label: 'bare word', pattern: word }]);
    testInfo.annotations.push({ type: 'TD-11', description: `user name "${machine.user}" as a bare word: ${bare.length} (only paths count)` });
    console.log(`TD-11: the user name "${machine.user}" occurs ${bare.length} times as a bare word in ${files.length} files`);
  }
  testInfo.annotations.push({ type: 'TD-11', description: `${files.length} files; checked ${JSON.stringify(machine)}` });
  expect(found).toEqual([]);
});
