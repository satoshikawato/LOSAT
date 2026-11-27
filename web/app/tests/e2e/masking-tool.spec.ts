// End-to-end tests for the Masking main tab (WW4c): the tool switcher (DustMasker, SegMasker,
// convert2blastmask) with each tool's -xmlhelp sentence and options, the FASTA input shared by the
// tools (paste, files, Separate runs), a run of each tool with its outputs and their download,
// convert2blastmask's required options refused in place until typed, the layout at 1280 and 390 px
// (no horizontal overflow, no orphaned `?` help tips), and zero external network requests.
//
// With the real engine (LOSAT_WEB_REACTORS) and the native binary (LOSAT_WEB_NATIVE): each tool runs once through the
// form with a maskinfo `-outfmt` (a different one each; two write to the standard output, one to a -out file), and
// its download is the native command line's bytes (the argv shown beside Run LOSAT, in a folder of the same files);
// "Use as -mask_data in makeblastdb" then hands that file to the makeblastdb tab, whose database files are the native
// makeblastdb's but the seconds of their dates (support/native-run.ts). Without the engine the same steps run on the
// FakeEngine's outputs and nothing is compared.
//
// LOSAT_E2E_SCREENSHOT_DIR=<dir> saves the Masking tab empty, filled and finished at 1280 and 390 px in the light and
// the dark colour scheme (<browser>-ms-<state>-<width>-<scheme>.png), and in Chromium each tool's finished run.
import { mkdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { gzipSync } from 'node:zlib';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE, ORIGIN } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { NATIVE } from './support/native';
import { downloadAll, masked, nativeMakeblastdb, nativeRun, ZONE } from './support/native-run';
import { openSearch, waitStatus } from './support/search';

const SCREENSHOT_DIR = process.env['LOSAT_E2E_SCREENSHOT_DIR'];
const FAKE_MARKER = 'FAKE ENGINE OUTPUT';
/** The name of the .gz input in the list: inflated by the engine build's Data worker only. */
const B_NAME = BUILD_HAS_ENGINE ? 'b.fa' : 'b.fa.gz';
const DESKTOP_VIEWPORT = { width: 1280, height: 900 };
const MOBILE_VIEWPORT = { width: 390, height: 844 };
const QUERY = '>q1 masking test\nACGTACGTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACGTTGCAGGCATTGCAACGT\n';

test.setTimeout(300_000);

function trackNet(page: Page): string[] {
  const external: string[] = [];
  page.on('request', (req) => {
    const url = req.url();
    if (!url.startsWith(ORIGIN) && !url.startsWith('blob:') && !url.startsWith('data:')) external.push(url);
  });
  return external;
}

async function openMasking(page: Page): Promise<void> {
  await page.getByTestId('tab-masking').click();
  await expect(page.getByTestId('tab-masking')).toHaveAttribute('aria-selected', 'true');
  await expect(page.getByTestId('ms-panel')).toHaveAttribute('data-describe', 'ready', { timeout: 60_000 });
}

async function chooseTool(page: Page, tool: 'dustmasker' | 'segmasker' | 'convert2blastmask'): Promise<void> {
  await page.getByTestId(`ms-tool-${tool}`).check();
  await expect(page.getByTestId('ms-panel')).toHaveAttribute('data-tool', tool);
  await expect(page.getByTestId('ms-panel')).toHaveAttribute('data-describe', 'ready', { timeout: 60_000 });
}

async function validIs(page: Page, state: 'ok' | 'invalid'): Promise<void> {
  await expect(page.getByTestId('ms-argv-validation')).toHaveAttribute('data-state', state, { timeout: 30_000 });
}

async function openOptions(page: Page): Promise<void> {
  const toggle = page.getByTestId('ms-options-toggle');
  if ((await toggle.getAttribute('aria-expanded')) !== 'true') await toggle.click();
  await expect(page.getByTestId('ms-parameter-form')).toBeVisible();
}

/** Runs the form as run `run` of the session and waits until it completes on the Results tab. */
async function runAndWait(page: Page, run: number): Promise<void> {
  await validIs(page, 'ok');
  await page.getByTestId('ms-run-losat').click();
  await waitStatus(page, run, 'completed');
  await expect(page.getByTestId('tab-results')).toHaveAttribute('aria-selected', 'true');
  await expect(page.getByTestId('tool-run-view')).toHaveAttribute('data-run', String(run));
}

/** The text of the output shown, and the bytes of its download. */
async function shownOutput(page: Page): Promise<{ readonly text: string; readonly bytes: Buffer }> {
  await expect(page.getByTestId('tool-output-text')).toBeVisible();
  const text = (await page.getByTestId('tool-output-text').textContent()) ?? '';
  const [download] = await Promise.all([page.waitForEvent('download'), page.getByTestId('tool-output-download').click()]);
  const path = await download.path();
  if (path === null) throw new Error('download path was null');
  return { text, bytes: readFileSync(path) };
}

async function maybeShot(page: Page, name: string): Promise<void> {
  if (!SCREENSHOT_DIR) return;
  mkdirSync(SCREENSHOT_DIR, { recursive: true });
  // Transitions (a tab that is still taking its colour) are finished before the screen is taken.
  await page.screenshot({ path: join(SCREENSHOT_DIR, `${test.info().project.name}-${name}`), fullPage: true, animations: 'disabled' });
}

/** Saves the page at 1280 and 390 px in the light and the dark colour scheme, back at 1280 px light after. */
async function shots(page: Page, state: string): Promise<void> {
  if (!SCREENSHOT_DIR) return;
  for (const scheme of ['light', 'dark'] as const) {
    await page.emulateMedia({ colorScheme: scheme });
    for (const viewport of [DESKTOP_VIEWPORT, MOBILE_VIEWPORT]) {
      await page.setViewportSize(viewport);
      await maybeShot(page, `ms-${state}-${viewport.width}-${scheme}.png`);
    }
  }
  await page.emulateMedia({ colorScheme: 'light' });
  await page.setViewportSize(DESKTOP_VIEWPORT);
}

/** No horizontal overflow, and every visible `?` of the Masking panel follows its content on the same line. */
async function assertLayout(page: Page, viewport: { readonly width: number; readonly height: number }): Promise<void> {
  await page.setViewportSize(viewport);
  const report = await page.evaluate(() => {
    const docOverflow = document.documentElement.scrollWidth - document.documentElement.clientWidth;
    const panel = document.querySelector('[data-testid="ms-panel"]')!;
    const panelRect = panel.getBoundingClientRect();
    const outside: string[] = [];
    for (const el of panel.querySelectorAll<HTMLElement>('.program-tab, input, select, textarea, button')) {
      const r = el.getBoundingClientRect();
      if (r.width < 0.5 || el.classList.contains('visually-hidden')) continue;
      if (r.left < panelRect.left - 1 || r.right > panelRect.right + 1) outside.push(el.dataset['testid'] ?? el.textContent?.trim() ?? el.tagName);
    }
    const orphanedTips: string[] = [];
    let visibleCount = 0;
    for (const btn of panel.querySelectorAll<HTMLElement>('.help-toggle')) {
      const br = btn.getBoundingClientRect();
      if (br.width < 0.5 || br.height < 0.5) continue;
      visibleCount += 1;
      const tip = btn.closest('.help-tip') ?? btn;
      let prev = tip.previousElementSibling;
      while (prev !== null && prev.getBoundingClientRect().width < 0.5) prev = prev.previousElementSibling;
      const wrapper = tip.closest('.hint-line, .param-hint-line');
      let prevRect: DOMRect | undefined;
      if (wrapper !== null && wrapper.firstChild !== null && wrapper.firstChild !== tip) {
        const range = document.createRange();
        range.setStartBefore(wrapper.firstChild);
        range.setEndBefore(tip);
        prevRect = [...range.getClientRects()].filter((r) => r.width > 0.5 && r.height > 0.5).at(-1);
      } else if (prev !== null) {
        prevRect = prev.getBoundingClientRect();
      }
      const overlapsY = prevRect !== undefined && br.top < prevRect.bottom - 1 && br.bottom > prevRect.top + 1;
      const followsX = prevRect !== undefined && br.left >= prevRect.right - 2;
      if (!overlapsY || !followsX) orphanedTips.push(btn.getAttribute('aria-label') || 'help-tip');
    }
    return { docOverflow, outside, orphanedTips, visibleCount };
  });
  expect(report.docOverflow, `horizontal overflow at ${viewport.width}px`).toBeLessThanOrEqual(1);
  expect(report.outside, `controls outside the panel at ${viewport.width}px`).toEqual([]);
  expect(report.visibleCount, `visible help tips at ${viewport.width}px`).toBeGreaterThan(0);
  expect(report.orphanedTips, `orphaned ? help tips at ${viewport.width}px`).toEqual([]);
  // Each tool's name stays inside its tab ("convert2blastmask" passed its tab's border at 390 px).
  const clipped = await page.getByTestId('ms-tool-tabs').evaluate((tabs) =>
    [...tabs.querySelectorAll<HTMLElement>('.program-tab')]
      .filter((tab) => tab.scrollWidth > tab.clientWidth + 1)
      .map((tab) => tab.textContent?.trim() ?? ''),
  );
  expect(clipped, `tool names wider than their tabs at ${viewport.width}px`).toEqual([]);
}

test('Masking tab: tool switcher, each tool run on the shared input, outputs and download, required options refused in place', async ({
  page,
}) => {
  const external = trackNet(page);
  await openSearch(page);
  await openMasking(page);

  // DustMasker is chosen first, with NCBI's sentence and its options from -xmlhelp.
  await expect(page.getByTestId('ms-tool-dustmasker')).toBeChecked();
  await expect(page.getByTestId('ms-tool-summary')).toHaveText('Low complexity region masker based on Symmetric DUST algorithm');
  await openOptions(page);
  await expect(page.getByTestId('ms-options-toggle')).toContainText('DustMasker options');
  for (const name of ['window', 'level', 'linker', 'parse_seqids', 'hard_masking']) {
    await expect(page.getByTestId(`ms-param-${name}-field`)).toBeVisible();
  }
  await expect(page.getByTestId('ms-required-block')).toHaveCount(0);
  await expect(page.getByTestId('ms-param-outfmt')).toBeVisible();

  // Run 1: DustMasker on the pasted FASTA, to the standard output.
  await page.getByTestId('ms-fasta-input').fill(QUERY);
  await expect(page.getByTestId('ms-command')).toHaveText('losat dustmasker -in input.fa');
  await page.getByTestId('ms-param-window').fill('32');
  await expect(page.getByTestId('ms-command')).toHaveText('losat dustmasker -in input.fa -window 32');
  await expect(page.getByTestId('ms-options-changed')).toHaveText('(1 changed)');
  await runAndWait(page, 1);
  await expect(page.getByTestId('results-program')).toHaveText('DustMasker');
  await expect(page.getByTestId('result-command')).toHaveText('losat dustmasker -in input.fa -window 32');
  const dust = await shownOutput(page);
  expect(dust.bytes.length).toBeGreaterThan(0);
  if (!BUILD_HAS_ENGINE) expect(dust.text).toContain(FAKE_MARKER);
  // DustMasker's -outfmt default (interval) is not masking information: no hand-over to makeblastdb.
  await expect(page.getByTestId('tool-use-as-mask-data')).toHaveCount(0);

  // Run 2: SegMasker on the same input (shared), with an -out file of its own.
  await openMasking(page);
  await chooseTool(page, 'segmasker');
  await expect(page.getByTestId('ms-tool-summary')).toHaveText('Low complexity region masker based on the SEG algorithm');
  await expect(page.getByTestId('ms-fasta-input')).toHaveValue(QUERY);
  await expect(page.getByTestId('ms-options-toggle')).toContainText('SegMasker options');
  for (const name of ['window', 'locut', 'hicut', 'parse_seqids']) await expect(page.getByTestId(`ms-param-${name}-field`)).toBeVisible();
  // SegMasker's own -window default from its -xmlhelp, not DustMasker's value; the command has no -window word.
  await expect(page.getByTestId('ms-param-window')).toHaveValue('12');
  await expect(page.getByTestId('ms-param-window-field')).not.toHaveClass(/changed/);
  await expect(page.getByTestId('ms-out')).toHaveValue('');
  await page.getByTestId('ms-out').fill('seg.txt');
  await expect(page.getByTestId('ms-command')).toHaveText('losat segmasker -in input.fa -out seg.txt');
  await runAndWait(page, 2);
  await expect(page.getByTestId('results-program')).toHaveText('SegMasker');
  await expect(page.getByTestId('tool-tab-file-0')).toHaveText('-out seg.txt');
  await page.getByTestId('tool-tab-file-0').click();
  const seg = await shownOutput(page);
  expect(seg.bytes.length).toBeGreaterThan(0);

  // convert2blastmask: its required options are refused in place until typed, as the command line refuses them.
  await openMasking(page);
  await chooseTool(page, 'convert2blastmask');
  await expect(page.getByTestId('ms-required-block')).toBeVisible();
  await expect(page.getByTestId('ms-param-masking_algorithm')).toHaveValue('');
  await expect(page.getByTestId('ms-param-masking_options')).toHaveValue('');
  await validIs(page, 'invalid');
  await expect(page.getByTestId('ms-argv-validation')).toContainText('masking_algorithm');
  await expect(page.getByTestId('ms-param-masking_algorithm')).toHaveAttribute('aria-invalid', 'true');
  await page.getByTestId('ms-run-losat').click();
  await expect(page.getByTestId('ms-message')).toContainText('masking_algorithm');
  await expect(page.getByTestId('run-3-option')).toHaveCount(0);
  await page.getByTestId('ms-param-masking_algorithm').fill('dust');
  await page.getByTestId('ms-param-masking_options').fill('window=64');
  await expect(page.getByTestId('ms-command')).toHaveText(
    'losat convert2blastmask -in input.fa -masking_algorithm dust -masking_options window=64',
  );
  await runAndWait(page, 3);
  await expect(page.getByTestId('results-program')).toHaveText('convert2blastmask');
  await expect(page.getByTestId('tool-tab-stdout')).toBeVisible();
  await page.getByTestId('tool-tab-stdout').click();
  const c2b = await shownOutput(page);
  expect(c2b.bytes.length).toBeGreaterThan(0);
  // convert2blastmask's -outfmt default (maskinfo_asn1_text, from its -xmlhelp) is masking information: its standard
  // output goes to makeblastdb's -mask_data, named for the run.
  await expect(page.getByTestId('tool-use-as-mask-data')).toBeVisible();
  await page.getByTestId('tool-use-as-mask-data').click();
  await expect(page.getByTestId('mk-panel')).toBeVisible();
  await expect(page.getByTestId('mk-param-mask_data')).toHaveValue('convert2blastmask-run3.asnt');
  await expect(page.getByTestId('mk-run-files')).toContainText('convert2blastmask-run3.asnt');
  await page.getByTestId('tab-results').click();
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-snapshot')).toContainText('convert2blastmask');

  // Back on DustMasker, its options and -out are as they were.
  await openMasking(page);
  await chooseTool(page, 'dustmasker');
  await expect(page.getByTestId('ms-param-window')).toHaveValue('32');
  await expect(page.getByTestId('ms-out')).toHaveValue('');

  expect(external, 'zero external network requests').toEqual([]);
});

test('Masking tab: FASTA files (a .gz read) shared by the tools, Separate runs, and the layout at 1280 and 390 px', async ({ page }) => {
  await openSearch(page);
  await openMasking(page);

  await page.getByTestId('ms-fasta-files').setInputFiles([
    { name: 'a.fa', mimeType: 'text/plain', buffer: Buffer.from('>a\nACGTACGTACGT\n') },
    { name: 'b.fa.gz', mimeType: 'application/gzip', buffer: gzipSync('>b\nTTTTTTTTTTTTTTTTTTTT\n') },
  ]);
  await expect(page.getByTestId('ms-fasta-name-0')).toHaveText('a.fa');
  // The engine build's Data worker inflates the .gz into b.fa; the FakeEngine build's inflates nothing (compressed-names.spec.ts).
  await expect(page.getByTestId('ms-fasta-list').locator('[data-status="ready"]')).toHaveCount(2, { timeout: 30_000 });
  await expect(page.getByTestId('ms-fasta-name-1')).toHaveText(B_NAME);
  if (BUILD_HAS_ENGINE) await expect(page.getByTestId('ms-fasta-inflated-1')).toHaveText('from b.fa.gz');
  await expect(page.getByTestId('ms-command')).toHaveText('losat dustmasker -in combined_input.fa');
  await page.getByTestId('ms-fasta-mode-separate').check();
  await expect(page.getByTestId('ms-run-losat')).toHaveText('Run LOSAT (2 runs)');
  await expect(page.getByTestId('ms-command')).toHaveText('losat dustmasker -in a.fa');

  // The input stays when the tool changes.
  await chooseTool(page, 'segmasker');
  await expect(page.getByTestId('ms-fasta-name-1')).toHaveText(B_NAME);
  await expect(page.getByTestId('ms-run-losat')).toHaveText('Run LOSAT (2 runs)');
  await validIs(page, 'ok');
  await page.getByTestId('ms-run-losat').click();
  await waitStatus(page, 1, 'completed');
  await waitStatus(page, 2, 'completed');
  await expect(page.getByTestId('result-command')).toHaveText('losat segmasker -in a.fa');

  // The layout of each tool at 1280 and 390 px, its options open.
  await openMasking(page);
  for (const tool of ['dustmasker', 'segmasker', 'convert2blastmask'] as const) {
    await chooseTool(page, tool);
    await openOptions(page);
    for (const [viewport, width] of [
      [DESKTOP_VIEWPORT, 1280],
      [MOBILE_VIEWPORT, 390],
    ] as const) {
      await assertLayout(page, viewport);
      for (const scheme of ['light', 'dark'] as const) {
        await page.emulateMedia({ colorScheme: scheme });
        await maybeShot(page, `ms-${tool}-${width}-${scheme}.png`);
      }
      await page.emulateMedia({ colorScheme: 'light' });
    }
    await page.setViewportSize(DESKTOP_VIEWPORT);
  }
});

// --- each tool through the form, its download against the native command line, then makeblastdb's -mask_data ---------

const LOSAT = join(REPOSITORY, 'LOSAT');

interface MaskingCase {
  readonly tool: 'dustmasker' | 'segmasker' | 'convert2blastmask';
  /** The FASTA of LOSAT/tests/<tool>/in/, pasted or opened as a file. */
  readonly input: string;
  readonly open: boolean;
  /** The `-outfmt` chosen in the form; undefined leaves the tool's default (convert2blastmask: maskinfo_asn1_text). */
  readonly outfmt: string | undefined;
  /** The extension of the mask data that the run names for makeblastdb (`maskDataName`). */
  readonly extension: string;
  /** A `-out` file name; undefined writes to the standard output. */
  readonly out: string | undefined;
  readonly dbtype: 'nucl' | 'prot';
  /** makeblastdb reads mask data in ASN.1 only: its run with a maskinfo_xml file is refused, as the command line refuses it. */
  readonly refusedByMakeblastdb?: true;
  /** The command beside Run LOSAT: the words in the tool's usage order. */
  readonly command: string;
}

const MASKING_CASES: readonly MaskingCase[] = [
  { tool: 'dustmasker', input: 'basic.fa', open: false, outfmt: 'maskinfo_asn1_bin', extension: 'asnb', out: 'dust.asnb', dbtype: 'nucl',
    command: 'losat dustmasker -in input.fa -out dust.asnb -outfmt maskinfo_asn1_bin -parse_seqids' },
  { tool: 'segmasker', input: 'basic.fa', open: true, outfmt: 'maskinfo_xml', extension: 'xml', out: undefined, dbtype: 'prot', refusedByMakeblastdb: true,
    command: 'losat segmasker -in basic.fa -parse_seqids -outfmt maskinfo_xml' },
  { tool: 'convert2blastmask', input: 'lc.fa', open: false, outfmt: undefined, extension: 'asnt', out: undefined, dbtype: 'nucl',
    command: 'losat convert2blastmask -in input.fa -parse_seqids -masking_algorithm dust -masking_options window=64' },
];

test.describe('Masking tab: each tool run, its output, and its mask data in makeblastdb', () => {
  // The engine writes the dates of a database in the browser's zone; the native runs use the same.
  test.use({ timezoneId: ZONE });

  for (const [index, c] of MASKING_CASES.entries()) {
    test(`${c.tool}: the download is the native command line's bytes, and makeblastdb's database made with it too`, async ({ page }) => {
      const real = BUILD_HAS_ENGINE && NATIVE !== undefined;
      const external = trackNet(page);
      const fasta = readFileSync(join(LOSAT, 'tests', c.tool, 'in', c.input));
      const first = index === 0;
      await openSearch(page);
      await openMasking(page);
      if (first) await shots(page, 'empty');
      await chooseTool(page, c.tool);

      // The form: the FASTA, convert2blastmask's required options, -outfmt, -parse_seqids, -out.
      if (c.open) {
        await page.getByTestId('ms-fasta-files').setInputFiles({ name: c.input, mimeType: 'text/plain', buffer: fasta });
        await expect(page.getByTestId('ms-fasta-list').locator('[data-status="ready"]')).toHaveCount(1, { timeout: 30_000 });
      } else {
        await page.getByTestId('ms-fasta-input').fill(fasta.toString('utf8'));
      }
      if (c.tool === 'convert2blastmask') {
        await page.getByTestId('ms-param-masking_algorithm').fill('dust');
        await page.getByTestId('ms-param-masking_options').fill('window=64');
      }
      await openOptions(page);
      await page.getByTestId('ms-param-parse_seqids').check();
      if (c.outfmt === undefined) {
        await expect(page.getByTestId('ms-param-outfmt')).toHaveValue('maskinfo_asn1_text');
      } else {
        await page.getByTestId('ms-param-outfmt').selectOption(c.outfmt);
      }
      if (c.out !== undefined) await page.getByTestId('ms-out').fill(c.out);
      await validIs(page, 'ok');
      const inName = c.open ? c.input : 'input.fa';
      await expect(page.getByTestId('ms-command')).toHaveText(c.command);
      await assertLayout(page, MOBILE_VIEWPORT);
      await assertLayout(page, DESKTOP_VIEWPORT);
      if (first) await shots(page, 'filled');

      // Run: the output that holds the mask data (the -out file, else the standard output), and its download.
      await runAndWait(page, 1);
      await expect(page.getByTestId('results-program')).toHaveText(
        { dustmasker: 'DustMasker', segmasker: 'SegMasker', convert2blastmask: 'convert2blastmask' }[c.tool],
      );
      await page.getByTestId(c.out === undefined ? 'tool-tab-stdout' : 'tool-tab-file-0').click();
      const argv = ((await page.getByTestId('result-command').textContent()) ?? '').trim().split(/\s+/);
      expect(argv.slice(0, 2), 'the command beside the run').toEqual(['losat', c.tool]);
      const made = await downloadActive(page);
      expect(made.length).toBeGreaterThan(0);
      if (test.info().project.name === 'chromium') await maybeShot(page, `ms-finished-${c.tool}-1280-light.png`);
      if (first) await shots(page, 'finished');
      await page.setViewportSize(MOBILE_VIEWPORT);
      expect(await page.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth), 'run view overflow at 390px').toBeLessThanOrEqual(1);
      await page.setViewportSize(DESKTOP_VIEWPORT);
      if (real) {
        const native = nativeRun(argv.slice(1), new Map([[inName, fasta]]));
        expect(native.status, `native ${argv.slice(1).join(' ')}: ${native.stderr.toString()}`).toBe(0);
        expect(made.equals(c.out === undefined ? native.stdout : native.made.get(c.out)!), `${c.tool} ${argv.slice(1).join(' ')}: the downloaded bytes`).toBe(true);
        if (c.out !== undefined) expect([...native.made.keys()], 'the files the native run made').toEqual([c.out]);
      }

      // Use as -mask_data in makeblastdb: the file goes to -mask_data, the same FASTA is -in, and the database is built.
      const maskName = c.out ?? `${c.tool}-run1.${c.extension}`;
      await page.getByTestId('tool-use-as-mask-data').click();
      await expect(page.getByTestId('mk-panel')).toBeVisible();
      await expect(page.getByTestId('mk-param-mask_data')).toHaveValue(maskName);
      await expect(page.getByTestId('mk-run-files')).toContainText(maskName);
      await page.getByTestId('mk-fasta-input').fill(fasta.toString('utf8'));
      await page.getByTestId(`mk-param-dbtype-${c.dbtype}`).check();
      await page.getByTestId('mk-param-out').fill('db');
      await page.getByTestId('mk-param-parse_seqids').check();
      await expect(page.getByTestId('mk-command')).toHaveText(
        `losat makeblastdb -in input.fa -dbtype ${c.dbtype} -parse_seqids -mask_data ${maskName} -out db`,
      );
      await expect(page.getByTestId('mk-argv-validation')).toHaveAttribute('data-state', 'ok', { timeout: 30_000 });
      const makeblastdbArgv = ['makeblastdb', '-in', 'input.fa', '-dbtype', c.dbtype, '-parse_seqids', '-mask_data', maskName, '-out', 'db'];
      const makeblastdbInputs = new Map([
        ['input.fa', fasta],
        [maskName, made],
      ]);
      await page.getByTestId('mk-run-losat').click();
      if (c.refusedByMakeblastdb && BUILD_HAS_ENGINE) {
        // The command line refuses the file (exit 1, its message on the standard error, no database): so does the run.
        await waitStatus(page, 2, 'failed');
        await page.getByTestId('tool-tab-stderr').click();
        const shown = (await page.getByTestId('tool-output-text').textContent()) ?? '';
        expect(shown).toContain(`Unknown encoding for mask data in '${maskName}'`);
        if (real) {
          const native = nativeRun(makeblastdbArgv, makeblastdbInputs, true);
          expect(native.status, 'native makeblastdb refuses the file').toBe(1);
          expect(shown, 'the run\'s standard error is the native one').toBe(native.stderr.toString());
          expect([...native.made.keys()], 'no database files').toEqual([]);
        }
        expect(external, 'requests outside the origin').toEqual([]);
        return;
      }
      await waitStatus(page, 2, 'completed');
      await expect(page.getByTestId('tool-run-view')).toHaveAttribute('data-run', '2');
      await expect(page.getByTestId('tool-made-files')).toBeVisible();
      const database = await downloadAll(page);
      const names = [...database.files.keys()];
      expect(names.some((name) => /^db\.[np]in$/.test(name)), `database files: ${names.join(' ')}`).toBe(true);
      if (real) {
        const native = nativeMakeblastdb(makeblastdbArgv, makeblastdbInputs);
        expect(masked(database.files), `${c.tool}: makeblastdb -mask_data ${maskName}`).toEqual(masked(native));
      }

      expect(external, 'requests outside the origin').toEqual([]);
    });
  }
});

/** The bytes that the active output's download saves. */
async function downloadActive(page: Page): Promise<Buffer> {
  const [download] = await Promise.all([page.waitForEvent('download'), page.getByTestId('tool-output-download').click()]);
  const path = await download.path();
  if (path === null) throw new Error('download path was null');
  return readFileSync(path);
}
