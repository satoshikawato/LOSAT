// The makeblastdb main tab (line WW, WW4b2): its form from `-xmlhelp` (the options in its groups, the command beside Run
// LOSAT, a missing -dbtype marked in its field), the layout at 1280 and 390 px (no horizontal scroll, no "?" alone on a
// line), a run's made files (each downloads, or all as one tar), "Use as Database (-db)" and "Use as Repeats database
// (-filtering_db)", which open Search with the database chosen, and a WindowMasker run's mask data handed to -mask_data
// ("Use as -mask_data in makeblastdb"), all without a request outside the page's origin.
//
// With the real engine (LOSAT_WEB_REACTORS) and the native binary (LOSAT_WEB_NATIVE): a nucleotide version 5 database
// with -parse_seqids made from pasted FASTA, and a protein database with -mask_data (seg) made from opened files, are the
// files of the native command line's run just before or just after the browser's byte for byte but the seconds of their
// times (run in a folder of the same files, in the browser's time zone; support/makeblastdb-native.ts), and the searches
// that "Use as Database (-db)" starts (blastn; blastp with -db_soft_mask 21) on the database in its run's folder
// (run1/db, run3/pdb) write the native command line's outfmt 0 and 6 on the files that the browser made, in that folder;
// and a WindowMasker run's mask data handed over ("Use as -mask_data in makeblastdb") makes the database that the native
// windowmasker and makeblastdb make on the same files.
//
// LOSAT_E2E_SCREENSHOT_DIR=<dir> saves the screens (empty tab, filled form, finished run, Search after "Use as Database
// (-db)") at 1280 and 390 px, in the light and the dark colour scheme.
import { spawnSync } from 'node:child_process';
import { mkdirSync, readFileSync, rmSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE, ORIGIN } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { asNative, folderOf, maskedFiles, nativeMakeblastdb, type Written } from './support/makeblastdb-native';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { openSearch, paste, pickRun, result, waitStatus } from './support/search';
import { test } from './support/test';

test.setTimeout(300_000);
// The engine writes the dates of a database in the browser's zone; the native run below uses the same.
const ZONE = 'Asia/Tokyo';
test.use({ timezoneId: ZONE });

const LOSAT = join(REPOSITORY, 'LOSAT');
const SCREENSHOT_DIR = process.env['LOSAT_E2E_SCREENSHOT_DIR'];
const DESKTOP = { width: 1280, height: 900 };
const MOBILE = { width: 390, height: 844 };
const NUCL_FASTA = readFileSync(join(LOSAT, 'tests/dropin_corpus/ns.fna'), 'utf8');

function trackNet(page: Page): string[] {
  const external: string[] = [];
  page.on('request', (req) => {
    const url = req.url();
    if (!url.startsWith(ORIGIN) && !url.startsWith('blob:') && !url.startsWith('data:')) external.push(url);
  });
  return external;
}

async function openMakeblastdb(page: Page): Promise<void> {
  await openSearch(page);
  await page.getByTestId('tab-makeblastdb').click();
  await expect(page.getByTestId('mk-panel')).toBeVisible();
  await expect(page.getByTestId('mk-panel')).toHaveAttribute('data-describe', 'ready');
}

async function openOptions(page: Page): Promise<void> {
  const toggle = page.getByTestId('mk-options-toggle');
  if ((await toggle.getAttribute('aria-expanded')) !== 'true') await toggle.click();
  await expect(page.getByTestId('mk-parameter-form')).toBeVisible();
}

async function settled(page: Page): Promise<void> {
  await expect(page.getByTestId('mk-argv-validation')).not.toHaveAttribute('data-state', 'checking');
}

/** Saves the page at 1280 and 390 px in the light and the dark colour scheme, back at 1280 px light after. */
async function shots(page: Page, name: string): Promise<void> {
  if (!SCREENSHOT_DIR) return;
  mkdirSync(SCREENSHOT_DIR, { recursive: true });
  const project = test.info().project.name;
  for (const scheme of ['light', 'dark'] as const) {
    await page.emulateMedia({ colorScheme: scheme });
    for (const viewport of [DESKTOP, MOBILE]) {
      await page.setViewportSize(viewport);
      await page.screenshot({ path: join(SCREENSHOT_DIR, `${project}-${name}-${viewport.width}-${scheme}.png`), fullPage: true });
    }
  }
  await page.emulateMedia({ colorScheme: 'light' });
  await page.setViewportSize(DESKTOP);
}

/** No horizontal scroll, and every visible "?" of the panel follows the text before it on its line. */
async function expectLayout(page: Page, panel: string, viewport: { readonly width: number; readonly height: number }): Promise<void> {
  await page.setViewportSize(viewport);
  const report = await page.evaluate((testid) => {
    const overflow = document.documentElement.scrollWidth - document.documentElement.clientWidth;
    const visible = (el: Element): boolean => {
      if (el.classList.contains('visually-hidden')) return false;
      const r = el.getBoundingClientRect();
      return r.width > 0.5 && r.height > 0.5;
    };
    const lastRect = (start: Node, endBefore?: Node): DOMRect | undefined => {
      const range = document.createRange();
      if (endBefore !== undefined) {
        range.setStartBefore(start);
        range.setEndBefore(endBefore);
      } else {
        range.selectNodeContents(start);
      }
      return [...range.getClientRects()].filter((r) => r.width > 0.5 && r.height > 0.5).at(-1);
    };
    const previous = (node: Element): Element | null => {
      let cur: Element | null = node.previousElementSibling;
      while (cur !== null && !visible(cur)) cur = cur.previousElementSibling;
      if (cur !== null) return cur;
      const parent = node.parentElement;
      return parent !== null && parent.classList.contains('control-tip') ? previous(parent) : null;
    };
    const alone: string[] = [];
    let tips = 0;
    for (const btn of document.querySelectorAll<HTMLElement>(`[data-testid="${testid}"] .help-toggle`)) {
      const br = btn.getBoundingClientRect();
      if (br.width < 0.5 || br.height < 0.5) continue;
      tips += 1;
      const tip = btn.closest('.help-tip') ?? btn;
      const wrapper = tip.closest('.hint-line, .param-hint-line');
      let before: DOMRect | undefined;
      if (wrapper !== null && wrapper.firstChild !== null && wrapper.firstChild !== tip) {
        before = lastRect(wrapper.firstChild, tip);
      } else {
        const prev = previous(tip);
        if (prev !== null) before = prev.matches('.param-hint, .hint') ? lastRect(prev) : prev.getBoundingClientRect();
      }
      const sameLine = before !== undefined && br.top < before.bottom - 1 && br.bottom > before.top + 1;
      const follows = before !== undefined && br.left >= before.right - 2;
      if (!sameLine || !follows) alone.push(btn.getAttribute('aria-label') || 'help-tip');
    }
    const cut: string[] = [];
    for (const el of document.querySelectorAll<HTMLElement>(`[data-testid="${testid}"] button, [data-testid="${testid}"] label`)) {
      if (visible(el) && el.scrollWidth > el.clientWidth + 1 && getComputedStyle(el).overflow === 'hidden') cut.push(el.textContent ?? '');
    }
    return { overflow, alone, tips, cut };
  }, panel);
  expect(report.overflow, `horizontal scroll at ${viewport.width}px`).toBeLessThanOrEqual(1);
  expect(report.alone, `"?" alone on a line at ${viewport.width}px`).toEqual([]);
  expect(report.cut, `cut-off labels at ${viewport.width}px`).toEqual([]);
}

/** The names and bytes of the files in a tar archive (ustar, with pax `path` records). */
function untar(bytes: Buffer): Map<string, Buffer> {
  const files = new Map<string, Buffer>();
  let at = 0;
  let paxPath: string | undefined;
  while (at + 512 <= bytes.length) {
    const header = bytes.subarray(at, at + 512);
    if (header.every((byte) => byte === 0)) break;
    const field = (from: number, length: number) => header.subarray(from, from + length).toString('utf8').replace(/\0.*$/s, '');
    const size = Number.parseInt(field(124, 12).trim(), 8);
    const type = String.fromCharCode(header[156]!);
    const body = bytes.subarray(at + 512, at + 512 + size);
    if (type === 'x') {
      paxPath = /\d+ path=(.*)\n/.exec(body.toString('utf8'))?.[1];
    } else {
      const prefix = field(345, 155);
      const name = paxPath ?? (prefix === '' ? field(0, 100) : `${prefix}/${field(0, 100)}`);
      files.set(name, Buffer.from(body));
      paxPath = undefined;
    }
    at += 512 + Math.ceil(size / 512) * 512;
  }
  return files;
}

async function downloadAll(page: Page): Promise<{ readonly name: string; readonly files: Map<string, Buffer> }> {
  const [download] = await Promise.all([page.waitForEvent('download'), page.getByTestId('tool-made-files-download-all').click()]);
  return { name: download.suggestedFilename(), files: untar(readFileSync((await download.path())!)) };
}

/** Runs the makeblastdb form as run `run` and waits until the run's view shows it finished. */
async function runMakeblastdb(page: Page, run: number): Promise<void> {
  await settled(page);
  await expect(page.getByTestId('mk-argv-validation')).toHaveAttribute('data-state', 'ok');
  await page.getByTestId('mk-run-losat').click();
  await waitStatus(page, run, 'completed');
  await expect(page.getByTestId('tool-run-view')).toHaveAttribute('data-run', String(run));
  await expect(page.getByTestId('tool-made-files')).toBeVisible();
}

test('makeblastdb tab: the options in the groups of -xmlhelp, the command, a missing -dbtype in its field, the layout', async ({ page }) => {
  const external = trackNet(page);
  await openMakeblastdb(page);
  await expect(page.getByTestId('mk-command')).toHaveText(/^losat makeblastdb/);
  await shots(page, 'mk-empty');

  await openOptions(page);
  const legends = await page.getByTestId('mk-parameter-form').locator('.section-legend').allTextContents();
  expect(legends.map((text) => text.trim())).toEqual(
    expect.arrayContaining(['Configuration options', 'Sequence masking options', 'Output options', 'Taxonomy options']),
  );
  for (const flag of ['title', 'parse_seqids', 'hash_index', 'mask_data', 'mask_id', 'gi_mask', 'blastdb_version', 'max_file_sz', 'taxid', 'taxid_map']) {
    await expect(page.getByTestId(`mk-param-${flag}-field`), flag).toBeVisible();
  }
  await expectLayout(page, 'mk-panel', DESKTOP);
  await expectLayout(page, 'mk-panel', MOBILE);
  await page.setViewportSize(DESKTOP);

  // Sequences without -dbtype: the CLI's refusal names -dbtype, whose radio buttons are marked.
  await page.getByTestId('mk-fasta-input').fill(NUCL_FASTA);
  await expect(page.getByTestId('mk-argv-validation')).toHaveAttribute('data-state', 'invalid', { timeout: 30_000 });
  await expect(page.getByTestId('mk-argv-validation')).toContainText('dbtype');
  await expect(page.getByTestId('mk-param-dbtype')).toHaveAttribute('aria-invalid', 'true');
  await page.getByTestId('mk-param-dbtype-nucl').check();
  await expect(page.getByTestId('mk-argv-validation')).toHaveAttribute('data-state', 'ok', { timeout: 30_000 });
  await expect(page.getByTestId('mk-param-dbtype')).not.toHaveAttribute('aria-invalid', 'true');

  await page.getByTestId('mk-param-out').fill('mydb');
  await page.getByTestId('mk-param-parse_seqids').check();
  await page.getByTestId('mk-param-title').fill('my genome');
  await expect(page.getByTestId('mk-command')).toHaveText("losat makeblastdb -in input.fa -dbtype nucl -title 'my genome' -parse_seqids -out mydb");
  await expect(page.getByTestId('mk-options-changed')).toHaveText('(2 changed)');
  await settled(page);
  await shots(page, 'mk-filled');
  await expectLayout(page, 'mk-panel', MOBILE);
  await page.setViewportSize(DESKTOP);

  // Restore default options clears the options of the groups and keeps the input, -dbtype and -out.
  await page.getByTestId('mk-restore-defaults').click();
  await expect(page.getByTestId('mk-command')).toHaveText('losat makeblastdb -in input.fa -dbtype nucl -out mydb');

  expect(external, 'requests outside the origin').toEqual([]);
});

test('makeblastdb run: the files it made, Download all, Use as Database (-db), Use as Repeats database (-filtering_db)', async ({ page }) => {
  const external = trackNet(page);
  await openMakeblastdb(page);
  await page.getByTestId('mk-fasta-input').fill(NUCL_FASTA);
  await page.getByTestId('mk-param-dbtype-nucl').check();
  await page.getByTestId('mk-param-out').fill('db/mydb');
  await runMakeblastdb(page, 1);

  const names = (await page.getByTestId('tool-made-files-list').locator('.source-name').allTextContents()).map((text) => text.trim());
  expect(names).toEqual(expect.arrayContaining(['db/mydb.nhr', 'db/mydb.nin', 'db/mydb.nsq']));
  const index = names.indexOf('db/mydb.nin');
  const [one] = await Promise.all([page.waitForEvent('download'), page.getByTestId(`tool-made-file-download-${index}`).click()]);
  expect(one.suggestedFilename()).toBe('mydb.nin');
  const all = await downloadAll(page);
  expect(all.name).toMatch(/^losat-run1-makeblastdb-.*db_mydb\.tar$/);
  expect([...all.files.keys()].sort()).toEqual([...names].sort());
  await expectLayout(page, 'tool-run-view', MOBILE);
  await page.setViewportSize(DESKTOP);
  await shots(page, 'mk-finished');

  await page.getByTestId('tool-use-as-database').click();
  await expect(page.getByTestId('search-panel')).toBeVisible();
  await expect(page.getByTestId('program-blastn')).toBeChecked();
  await expect(page.getByTestId('subject-source-database')).toBeChecked();
  await expect(page.getByTestId('subject-db-selected')).toContainText('run1/db/mydb');
  await shots(page, 'mk-search-after-use-as-database');

  await pickRun(page, 1);
  await page.getByTestId('tool-use-as-repeats-database').click();
  await expect(page.getByTestId('search-panel')).toBeVisible();
  await expect(page.getByTestId('param-filtering_db-selected')).toContainText('run1/db/mydb');

  // A protein database: Use as Database (-db) chooses blastp, and there is no Repeats database button.
  await page.getByTestId('tab-makeblastdb').click();
  await page.getByTestId('mk-fasta-input').fill('>p1\nMKVLAAGIVGLLLAGCSSHE\n');
  await page.getByTestId('mk-param-dbtype-prot').check();
  await page.getByTestId('mk-param-out').fill('pdb');
  await runMakeblastdb(page, 2);
  await expect(page.getByTestId('tool-use-as-repeats-database')).toHaveCount(0);
  await page.getByTestId('tool-use-as-database').click();
  await expect(page.getByTestId('program-blastp')).toBeChecked();
  await expect(page.getByTestId('subject-db-selected')).toContainText('run2/pdb');

  expect(external, 'requests outside the origin').toEqual([]);
});

test("WindowMasker's mask data goes to makeblastdb's -mask_data (Use as -mask_data in makeblastdb)", async ({ page }) => {
  test.skip(BUILD_HAS_ENGINE, "the fake engine masks without counts; the real engine's hand-over is the test \"real engine: WindowMasker's mask data …\"");
  await openSearch(page);
  await page.getByTestId('tab-windowmasker').click();
  await page.getByTestId('wm-stage-mask').check();
  await page.getByTestId('wm-ustat-file').setInputFiles({ name: 'counts.obinary', mimeType: 'application/octet-stream', buffer: Buffer.from('counts') });
  await page.getByTestId('wm-fasta-input').fill('>a\nACGTACGTACGT\n');
  await page.getByTestId('wm-options-toggle').click();
  await page.getByTestId('wm-param-outfmt').selectOption('maskinfo_asn1_bin');
  await expect(page.getByTestId('wm-argv-validation')).not.toHaveAttribute('data-state', 'checking');
  await page.getByTestId('wm-run-losat').click();
  await waitStatus(page, 1, 'completed');
  await page.getByTestId('tool-use-as-mask-data').click();

  await expect(page.getByTestId('mk-panel')).toBeVisible();
  await expect(page.getByTestId('mk-parameter-form')).toBeVisible();
  await expect(page.getByTestId('mk-param-mask_data')).toHaveValue('windowmasker-run1.asnb');
  await expect(page.getByTestId('mk-run-files')).toContainText('windowmasker-run1.asnb');
  await expect(page.getByTestId('mk-options-changed')).toHaveText('(1 changed)');

  // A later run's mask data under the name of an earlier run's takes its place (WW4b2 R1 finding 4).
  await page.getByTestId('tab-windowmasker').click();
  await page.getByTestId('wm-out').fill('wm.asnb');
  for (const run of [2, 3]) {
    await page.getByTestId('tab-windowmasker').click();
    await expect(page.getByTestId('wm-argv-validation')).not.toHaveAttribute('data-state', 'checking');
    await page.getByTestId('wm-run-losat').click();
    await waitStatus(page, run, 'completed');
    await page.getByTestId('tool-use-as-mask-data').click();
    await expect(page.getByTestId('mk-param-mask_data')).toHaveValue('windowmasker-run1.asnb,wm.asnb');
    await expect(page.getByTestId('mk-run-files-list').locator('.source-name')).toHaveText(['windowmasker-run1.asnb', 'wm.asnb']);
    await expect(page.getByTestId('mk-run-files-error')).toHaveCount(0);
    await expect(page.getByTestId('mk-message')).toHaveCount(0);
  }

  // A hand-over that the run directory refuses (a file that the user opened holds the name) says why at the top of the
  // tab, where it opens, after a run from the bottom Run LOSAT too.
  await page.getByTestId('mk-run-file-remove-1').click();
  await page.getByTestId('mk-run-files-input').setInputFiles({ name: 'wm.asnb', mimeType: 'application/octet-stream', buffer: Buffer.from('mine') });
  await page.getByTestId('mk-fasta-input').fill(NUCL_FASTA);
  await page.getByTestId('mk-param-dbtype-nucl').check();
  await settled(page);
  await page.getByTestId('mk-run-losat-bottom').click();
  await waitStatus(page, 4, 'completed');
  await page.getByTestId('tab-windowmasker').click();
  await page.getByTestId('wm-run-losat').click();
  await waitStatus(page, 5, 'completed');
  await page.getByTestId('tool-use-as-mask-data').click();
  await expect(page.getByTestId('mk-message')).toHaveText('A file named "wm.asnb" is already in the run directory.');
  await expect(page.getByTestId('mk-message-bottom')).toHaveCount(0);
  await expect(page.getByTestId('mk-param-mask_data')).toHaveValue('windowmasker-run1.asnb');
});

// --- the real engine, against the native command line ----------------------------------------------------------------

/** The native makeblastdb's run of `argv` in a folder of `inputs`, in the browser's zone; a run that fails stops the test. */
function nativeMade(argv: readonly string[], inputs: ReadonlyMap<string, Buffer>): Written {
  const run = nativeMakeblastdb(argv, inputs, ZONE);
  if (run.status !== 0) throw new Error(`native ${argv.join(' ')} exited ${run.status}: ${run.stderr.toString()}`);
  return run;
}

/** `files` under `folder`: where "Use as Database (-db)" puts a run's database in Search (run<N>/, one folder per run). */
const inFolder = (folder: string, files: ReadonlyMap<string, Buffer>) => [...files].map(([name, bytes]) => [`${folder}/${name}`, bytes] as const);

/** The native search's outfmt `format` in a folder of `files` (the database that the browser made, and the query). */
function nativeSearch(argv: readonly string[], files: ReadonlyMap<string, Buffer>, format: 0 | 6): string {
  const folder = folderOf(files);
  try {
    const run = spawnSync(NATIVE!, [...argv, '-outfmt', String(format), '-num_threads', '1'], { cwd: folder, maxBuffer: 1 << 28 });
    if (run.status !== 0) throw new Error(`native ${argv.join(' ')} exited ${run.status}: ${run.stderr.toString()}`);
    return run.stdout.toString('utf8');
  } finally {
    rmSync(folder, { recursive: true, force: true });
  }
}

test('real engine: the databases of the tab are the native command line\'s, and their searches are too', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, 'requires built Wasm reactors');
  test.skip(NATIVE === undefined, NO_NATIVE_REASON ?? 'requires the native losat binary');
  const external = trackNet(page);
  await openMakeblastdb(page);

  // 1. A nucleotide version 5 database with -parse_seqids from pasted FASTA; blastn -db on it.
  await page.getByTestId('mk-fasta-input').fill(NUCL_FASTA);
  await page.getByTestId('mk-param-dbtype-nucl').check();
  await page.getByTestId('mk-param-out').fill('db');
  await openOptions(page);
  await page.getByTestId('mk-param-parse_seqids').check();
  await expect(page.getByTestId('mk-command')).toHaveText('losat makeblastdb -in input.fa -dbtype nucl -parse_seqids -out db');
  const nuclArgv = ['makeblastdb', '-in', 'input.fa', '-dbtype', 'nucl', '-parse_seqids', '-out', 'db'];
  const nuclInputs = new Map([['input.fa', Buffer.from(NUCL_FASTA)]]);
  const nuclBefore = nativeMade(nuclArgv, nuclInputs);
  await runMakeblastdb(page, 1);
  const nucl = await downloadAll(page);
  const nuclAfter = nativeMade(nuclArgv, nuclInputs);
  expect(maskedFiles(nucl.files)).toEqual(asNative(maskedFiles(nucl.files), [maskedFiles(nuclBefore.files), maskedFiles(nuclAfter.files)]));

  const query = readFileSync(join(LOSAT, 'tests/dropin_corpus/nq.fna'));
  await page.getByTestId('tool-use-as-database').click();
  await expect(page.getByTestId('subject-db-selected')).toContainText('run1/db');
  await paste(page, 'query', query.toString('utf8'));
  await expect(page.getByTestId('argv-validation')).not.toHaveAttribute('data-state', 'checking');
  await page.getByTestId('run-losat').click();
  await waitStatus(page, 2, 'completed');
  const nuclFiles = new Map([...inFolder('run1', nucl.files), ['q.fna', query]]);
  for (const format of [0, 6] as const) {
    const shown = await result(page, 2, format);
    expect(shown.output, `blastn -db outfmt ${format}`).toBe(nativeSearch(['blastn', '-db', 'run1/db', '-query', 'q.fna'], nuclFiles, format));
  }

  // 2. A protein database with -mask_data (seg) from opened files; blastp -db_soft_mask 21 on it.
  const protein = readFileSync(join(LOSAT, 'tests/seqdb/mask_prot.faa'));
  const seg = readFileSync(join(LOSAT, 'tests/seqdb/mask_seg.asnb'));
  await page.getByTestId('tab-makeblastdb').click();
  await page.getByTestId('mk-fasta-clear').click();
  await page.getByTestId('mk-input-files').setInputFiles({ name: 'prot.faa', mimeType: 'text/plain', buffer: protein });
  await page.getByTestId('mk-param-dbtype-prot').check();
  await page.getByTestId('mk-param-out').fill('pdb');
  await page.getByTestId('mk-param-parse_seqids').uncheck();
  await page.getByTestId('mk-param-mask_data-file').setInputFiles({ name: 'seg.asnb', mimeType: 'application/octet-stream', buffer: seg });
  await expect(page.getByTestId('mk-command')).toHaveText('losat makeblastdb -in prot.faa -dbtype prot -mask_data seg.asnb -out pdb');
  await shots(page, 'mk-mask-data-form');
  const protArgv = ['makeblastdb', '-in', 'prot.faa', '-dbtype', 'prot', '-mask_data', 'seg.asnb', '-out', 'pdb'];
  const protInputs = new Map([
    ['prot.faa', protein],
    ['seg.asnb', seg],
  ]);
  const protBefore = nativeMade(protArgv, protInputs);
  await runMakeblastdb(page, 3);
  await expect(page.getByTestId('tool-output-text')).toContainText('Adding sequences from FASTA');
  await shots(page, 'mk-mask-data-finished');
  const prot = await downloadAll(page);
  const protAfter = nativeMade(protArgv, protInputs);
  expect(maskedFiles(prot.files)).toEqual(asNative(maskedFiles(prot.files), [maskedFiles(protBefore.files), maskedFiles(protAfter.files)]));

  await page.getByTestId('tool-use-as-database').click();
  await expect(page.getByTestId('program-blastp')).toBeChecked();
  await expect(page.getByTestId('subject-db-selected')).toContainText('run3/pdb');
  await paste(page, 'query', protein.toString('utf8'));
  await page.getByTestId('subject-db-filters-toggle').click();
  await page.getByTestId('db-mask-mode').selectOption('soft');
  await page.getByTestId('db-mask-algorithm').fill('21');
  await expect(page.getByTestId('argv-validation')).not.toHaveAttribute('data-state', 'checking');
  await page.getByTestId('run-losat').click();
  await waitStatus(page, 4, 'completed');
  const protFiles = new Map([...inFolder('run3', prot.files), ['q.faa', protein]]);
  for (const format of [0, 6] as const) {
    const shown = await result(page, 4, format);
    expect(shown.command, 'blastp command').toContain('-db_soft_mask 21');
    expect(shown.output, `blastp -db -db_soft_mask 21 outfmt ${format}`).toBe(
      nativeSearch(['blastp', '-db', 'run3/pdb', '-query', 'q.faa', '-db_soft_mask', '21'], protFiles, format),
    );
  }

  expect(external, 'requests outside the origin').toEqual([]);
});

test("real engine: WindowMasker's mask data handed to makeblastdb's -mask_data makes the native command line's database", async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, 'requires built Wasm reactors');
  test.skip(NATIVE === undefined, NO_NATIVE_REASON ?? 'requires the native losat binary');
  const external = trackNet(page);
  const fasta = readFileSync(join(LOSAT, 'tests/windowmasker/in/multi.fa'));
  const counts = readFileSync(join(LOSAT, 'tests/winmask/counts/dm2L_u9.binary'));

  // WindowMasker's Mask with -outfmt maskinfo_asn1_bin (to its standard output), then "Use as -mask_data in makeblastdb".
  await openSearch(page);
  await page.getByTestId('tab-windowmasker').click();
  await page.getByTestId('wm-stage-mask').check();
  await page.getByTestId('wm-ustat-file').setInputFiles({ name: 'dm2L_u9.binary', mimeType: 'application/octet-stream', buffer: counts });
  await page.getByTestId('wm-fasta-files').setInputFiles({ name: 'multi.fa', mimeType: 'application/octet-stream', buffer: fasta });
  await expect(page.locator('[data-testid="wm-fasta-list"] > li').first()).toHaveAttribute('data-status', 'ready', { timeout: 60_000 });
  await page.getByTestId('wm-options-toggle').click();
  await page.getByTestId('wm-param-outfmt').selectOption('maskinfo_asn1_bin');
  await expect(page.getByTestId('wm-command')).toHaveText('losat windowmasker -ustat dm2L_u9.binary -in multi.fa -outfmt maskinfo_asn1_bin');
  await expect(page.getByTestId('wm-argv-validation')).not.toHaveAttribute('data-state', 'checking');
  await page.getByTestId('wm-run-losat').click();
  await waitStatus(page, 1, 'completed');
  await page.getByTestId('tool-use-as-mask-data').click();
  await expect(page.getByTestId('mk-param-mask_data')).toHaveValue('windowmasker-run1.asnb');

  // makeblastdb on the same FASTA with that mask data.
  await page.getByTestId('mk-input-files').setInputFiles({ name: 'multi.fa', mimeType: 'text/plain', buffer: fasta });
  await page.getByTestId('mk-param-dbtype-nucl').check();
  await page.getByTestId('mk-param-out').fill('wmdb');
  await expect(page.getByTestId('mk-command')).toHaveText('losat makeblastdb -in multi.fa -dbtype nucl -mask_data windowmasker-run1.asnb -out wmdb');

  // The native command line's windowmasker and makeblastdb on the same files.
  const maskFolder = folderOf(new Map([['dm2L_u9.binary', counts], ['multi.fa', fasta]]));
  let mask: Buffer;
  try {
    const run = spawnSync(NATIVE!, ['windowmasker', '-ustat', 'dm2L_u9.binary', '-in', 'multi.fa', '-outfmt', 'maskinfo_asn1_bin'], { cwd: maskFolder, maxBuffer: 1 << 28 });
    if (run.status !== 0) throw new Error(`native windowmasker exited ${run.status}: ${run.stderr.toString()}`);
    mask = run.stdout;
  } finally {
    rmSync(maskFolder, { recursive: true, force: true });
  }
  const argv = ['makeblastdb', '-in', 'multi.fa', '-dbtype', 'nucl', '-mask_data', 'windowmasker-run1.asnb', '-out', 'wmdb'];
  const inputs = new Map([
    ['multi.fa', fasta],
    ['windowmasker-run1.asnb', mask],
  ]);
  const before = nativeMade(argv, inputs);
  await runMakeblastdb(page, 2);
  const made = await downloadAll(page);
  const after = nativeMade(argv, inputs);
  // The mask's column (.naa, .nab, .nac) is among the files: the mask data was read.
  expect([...made.files.keys()]).toEqual(expect.arrayContaining(['wmdb.naa', 'wmdb.nab', 'wmdb.nac']));
  expect(maskedFiles(made.files)).toEqual(asNative(maskedFiles(made.files), [maskedFiles(before.files), maskedFiles(after.files)]));

  expect(external, 'requests outside the origin').toEqual([]);
});
