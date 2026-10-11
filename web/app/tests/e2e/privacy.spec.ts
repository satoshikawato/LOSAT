// V-PRIV (plan §5.10, §6.2; design §13.1, §13.3): research data never goes over the network. A
// unique canary of each test is put into everything that is research data - the file names, the
// deflines, the sequences (a nucleotide 40-mer, a protein motif), the job title, the candidates'
// notes, the names of the settings and session files the user opens, and inputs that are refused -
// and each path of the application is taken: a search and its results (Graphic Summary,
// Descriptions, Alignments, Dot Plot), every export, a reload, the saved session opened and its
// originals attached again; a cancelled search, a refused input, a refused option and a failed run.
// Every request is captured at the server and in the browser (support/privacy.ts), and none may go
// to another origin, send a body, ask for a path that is not a file of the site or carry a canary;
// no CSP violation may occur, and no page title or address may hold a canary. Downloads are made in
// the page: their contents hold the canaries, which shows that the data was there and stayed.
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
// TD-11 (chromium project only: it reads files): no file that the site serves holds a path of the
// build machine (the checkout, CARGO_HOME, RUSTUP_HOME, HOME) or the user name.
//
// To come with the work packages that add them (S16): the update path (service worker, explicit
// update; the caches hold only release files) and the measurement path (consent, /api/visit).
import { execSync } from 'node:child_process';
import { realpathSync } from 'node:fs';
import { readFile } from 'node:fs/promises';
import { homedir, userInfo } from 'node:os';
import { dirname, join, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { gunzipSync } from 'node:zlib';
import { expect, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { DIST, PrivacyCapture, startSiteServer, test } from './support/privacy';
import {
  filesUnder,
  machinePatterns,
  makeCanary,
  scanFiles,
  userNameInPathsOnly,
  type Canary,
  type MachineContext,
} from './support/privacy-scan';
import { fasta, openFiles, openParameters, paste, program, settled, showOutput, submit, waitStatus } from './support/search';

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

/** Opens the application from the capture's server and waits until its start-up has finished. */
async function openApp(page: Page, capture: PrivacyCapture): Promise<void> {
  await page.goto(`${capture.origin}/`);
  await expect(page.getByTestId('storage-status')).toHaveAttribute('data-cleanup', 'done');
}

/** The capture saw the site's own requests: the page, its script and the Data worker's (and the engine's). */
function expectSiteRequests(capture: PrivacyCapture): void {
  const paths = capture.serverPaths();
  expect(paths).toContain('/');
  expect(paths.some((path) => /^\/assets\/index-[\w-]+\.js$/.test(path))).toBe(true);
  expect(paths.some((path) => /^\/assets\/data-worker-[\w-]+\.js$/.test(path))).toBe(true);
  if (BUILD_HAS_ENGINE) {
    expect(paths.some((path) => /^\/assets\/engine-worker-[\w-]+\.js$/.test(path))).toBe(true);
    expect(paths.some((path) => /^\/engine\/losat-web-(serial|threads)-[0-9a-f]+\.wasm$/.test(path))).toBe(true);
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
  await expect(page.getByTestId('storage-status')).toHaveAttribute('data-cleanup', 'done');
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

test('the capture finds what it looks for: a canary in a URL, a path that is not the site\'s, CSP violations in a page and in a worker', async ({
  context,
  browserName,
}, testInfo) => {
  // Another origin that logs what reaches it, and a worker that only this test serves: the page and
  // the worker each try a request to that origin, which the CSP blocks.
  const foreign = await startSiteServer(DIST);
  const worker = "fetch(new URL(location.href).searchParams.get('to')).then(() => postMessage('sent'), () => postMessage('failed'));";
  const server = await startSiteServer(DIST, { extra: new Map([['/__privacy/planted-worker.js', { type: 'text/javascript', body: worker }]]) });
  const capture = new PrivacyCapture(server);
  try {
    await capture.attach(context);
    const page = await context.newPage();
    await openApp(page, capture);
    const canary = makeCanary();
    await page.evaluate(
      async ({ token, to }) => {
        await fetch(`/__privacy/planted?note=${token}`);
        await fetch(`${to}/${token}`).catch(() => undefined);
        const planted = new Worker(`/__privacy/planted-worker.js?to=${encodeURIComponent(`${to}/from-a-worker`)}`);
        await new Promise((resolve) => (planted.onmessage = resolve));
      },
      { token: canary.token, to: foreign.origin },
    );
    const reports = async () => (await capture.findings(canary)).filter((line) => line.startsWith('a CSP violation report'));
    // Browsers send CSP reports in their own time. WebKit sends none for a worker (S16: seen in WebKit 26.x).
    const workerReport = expect.stringMatching(/^a CSP violation report: connect-src blocked .* in .*\/__privacy\/planted-worker\.js/);
    if (browserName !== 'webkit') await expect.poll(reports, { timeout: 15_000 }).toEqual(expect.arrayContaining([workerReport]));
    await expect.poll(async () => (await reports()).length, { timeout: 15_000 }).toBeGreaterThanOrEqual(1);
    testInfo.annotations.push({ type: 'CSP reports', description: (await reports()).join('\n') });
    expect(await capture.findings(canary)).toEqual(
      expect.arrayContaining([
        expect.stringMatching(/^canary \(token\) in the URL: server GET http:\/\/127\.0\.0\.1:\d+\/__privacy\/planted\?note=LOSATPRIV/),
        expect.stringMatching(/^canary \(token\) in the URL: browser GET http:\/\/127\.0\.0\.1:\d+\/__privacy\/planted\?note=LOSATPRIV/),
        expect.stringMatching(/^a request for a path that is not a file of the site: server GET .*\/__privacy\/planted\?note=/),
        expect.stringMatching(new RegExp(`^a CSP violation in a page: connect-src blocked ${foreign.origin}/LOSATPRIV`)),
        expect.stringMatching(new RegExp(`^a CSP violation report: connect-src blocked ${foreign.origin}/LOSATPRIV\\w+ in ${server.origin}/$`)),
      ]),
    );
    // The CSP held: nothing reached the other origin.
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
  const files = filesUnder(DIST);
  expect(files).toContain('index.html');
  expect(files.some((file) => file.endsWith('.js'))).toBe(true);
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
