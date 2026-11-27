// End-to-end tests for the Masking main tab (WW4c): the tool switcher (DustMasker, SegMasker,
// convert2blastmask) with each tool's -xmlhelp sentence and options, the FASTA input shared by the
// tools (paste, files, Separate runs), a run of each tool with its outputs and their download,
// convert2blastmask's required options refused in place until typed, the layout at 1280 and 390 px
// (no horizontal overflow, no orphaned `?` help tips), and zero external network requests.
import { mkdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { gzipSync } from 'node:zlib';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE, ORIGIN } from './support/browser';
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
  await page.screenshot({ path: join(SCREENSHOT_DIR, name), fullPage: true });
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
