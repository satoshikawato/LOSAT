// End-to-end tests for the Format Archive (`blast_formatter`) main tab (RG-F3p7):
// -xmlhelp loading, archive file / paste / session -outfmt 11 shortcut ("Format archive…"),
// local -db archive with run directory folder, byte-for-byte parity against native
// `losat blast_formatter` for every -outfmt choice the tab offers (each exits 0 on the engine) and -out,
// validateTool refusal of -rid, responsive layout at 1280 px & 390 px (no orphaned `?` help tips, no
// horizontal overflow, also of the Outputs with "Format archive…"), and zero external network requests.
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { mkdirSync, readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE, ORIGIN } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { described, openSearch, paste, program, submit, waitStatus } from './support/search';

const BF_DIR = join(REPOSITORY, 'LOSAT/tests/blast_formatter');
const SCREENSHOT_DIR = process.env['LOSAT_E2E_SCREENSHOT_DIR'];

test.setTimeout(300_000);

const sha256 = (bytes: Uint8Array | string): string => createHash('sha256').update(bytes).digest('hex');

function runNativeBf(
  args: readonly string[],
  cwd: string,
  stdin?: Buffer,
): { readonly status: number; readonly stdout: Buffer; readonly stderr: Buffer } {
  if (NATIVE === undefined) throw new Error(NO_NATIVE_REASON);
  const result = spawnSync(NATIVE, ['blast_formatter', ...args], {
    cwd,
    input: stdin,
    maxBuffer: 1 << 28,
  });
  if (result.error) throw result.error;
  return { status: result.status ?? -1, stdout: result.stdout, stderr: result.stderr };
}

function trackNet(page: Page): string[] {
  const external: string[] = [];
  page.on('request', (req) => {
    const url = req.url();
    if (!url.startsWith(ORIGIN) && !url.startsWith('blob:') && !url.startsWith('data:')) {
      external.push(url);
    }
  });
  return external;
}

const DESKTOP_VIEWPORT = { width: 1280, height: 900 };
const MOBILE_VIEWPORT = { width: 390, height: 844 };

async function clickBfTab(page: Page): Promise<void> {
  await page.getByTestId('tab-blast-formatter').click();
}

async function expectBfPanelVisible(page: Page): Promise<void> {
  await expect(page.getByTestId('bf-panel')).toBeVisible();
}

async function openFormatArchive(page: Page): Promise<void> {
  await openSearch(page);
  await clickBfTab(page);
  await expectBfPanelVisible(page);
  await described(page, 'bf-panel', 'ready');
}

/** Runs the Format Archive form as run `run` of the session, and waits until it exits 0. */
async function runBfAndWait(page: Page, run: number): Promise<void> {
  await page.getByTestId('bf-run-losat').click();
  await waitStatus(page, run, 'completed');
  await expect(page.getByTestId('tool-run-view')).toHaveAttribute('data-run', String(run));
}

/** The width by which the page is wider than its window (0 when it fits). */
async function pageOverflow(page: Page): Promise<number> {
  return page.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth);
}

async function expectTextOutputVisible(page: Page): Promise<void> {
  await expect(page.getByTestId('tool-output-text')).toBeVisible();
}

async function maybeShot(page: Page, name: string, viewport?: { readonly width: number; readonly height: number }): Promise<void> {
  if (!SCREENSHOT_DIR) return;
  if (viewport !== undefined) await page.setViewportSize(viewport);
  mkdirSync(SCREENSHOT_DIR, { recursive: true });
  await page.screenshot({ path: join(SCREENSHOT_DIR, name), fullPage: true });
}

async function resetDesktopViewport(page: Page): Promise<void> {
  await page.setViewportSize(DESKTOP_VIEWPORT);
}

function expectNoExternalRequests(external: readonly string[]): void {
  expect(external, 'zero external network requests').toEqual([]);
}

/** Reads the active tool output tab's full downloaded bytes via Playwright download interception. */
async function downloadActiveToolOutput(page: Page): Promise<Buffer> {
  const [download] = await Promise.all([
    page.waitForEvent('download'),
    page.getByTestId('tool-output-download').click(),
  ]);
  const path = await download.path();
  if (path === null) throw new Error('download path was null');
  return readFileSync(path);
}

/** Asserts every visible `?` (.help-toggle) shares a line with its preceding content and no horizontal overflow occurs. */
async function assertLayoutAndHelpTips(page: Page, viewport: { width: number; height: number }): Promise<void> {
  await page.setViewportSize(viewport);
  const report = await page.evaluate(() => {
    const docOverflow = document.documentElement.scrollWidth - document.documentElement.clientWidth;
    const isVisible = (el: Element): boolean => {
      if (el.classList.contains('visually-hidden')) return false;
      const r = el.getBoundingClientRect();
      return r.width > 0.5 && r.height > 0.5;
    };
    const lastRangeRect = (start: Node, endBefore?: Node): DOMRect | undefined => {
      const range = document.createRange();
      if (endBefore !== undefined) {
        range.setStartBefore(start);
        range.setEndBefore(endBefore);
      } else {
        range.selectNodeContents(start);
      }
      return [...range.getClientRects()].filter((r) => r.width > 0.5 && r.height > 0.5).at(-1);
    };
    const findPrev = (node: Element): Element | null => {
      let cur: Element | null = node.previousElementSibling;
      while (cur !== null && !isVisible(cur)) cur = cur.previousElementSibling;
      if (cur !== null) return cur;
      const parent = node.parentElement;
      if (parent !== null && parent.classList.contains('control-tip')) return findPrev(parent);
      return null;
    };
    const orphanedTips: string[] = [];
    let visibleCount = 0;
    for (const btn of document.querySelectorAll<HTMLElement>('[data-testid="bf-panel"] .help-toggle')) {
      const br = btn.getBoundingClientRect();
      if (br.width < 0.5 || br.height < 0.5) continue;
      visibleCount += 1;
      const tip = btn.closest('.help-tip') ?? btn;
      const wrapper = tip.closest('.hint-line, .param-hint-line');
      let prevRect: DOMRect | undefined;
      if (wrapper !== null && wrapper.firstChild !== null && wrapper.firstChild !== tip) {
        prevRect = lastRangeRect(wrapper.firstChild, tip);
      } else {
        const prev = findPrev(tip);
        if (prev !== null) {
          prevRect = prev.matches('.param-hint, .hint') ? lastRangeRect(prev) : prev.getBoundingClientRect();
        }
      }
      const overlapsY = prevRect !== undefined && br.top < prevRect.bottom - 1 && br.bottom > prevRect.top + 1;
      const followsX = prevRect !== undefined && br.left >= prevRect.right - 2;
      if (!overlapsY || !followsX) {
        orphanedTips.push(btn.getAttribute('aria-label') || 'help-tip');
      }
    }
    return { docOverflow, orphanedTips, visibleCount };
  });
  expect(report.docOverflow, `horizontal overflow at ${viewport.width}px`).toBeLessThanOrEqual(1);
  expect(report.visibleCount, `visible help tips at ${viewport.width}px`).toBeGreaterThan(0);
  expect(report.orphanedTips, `orphaned ? help tips at ${viewport.width}px`).toEqual([]);
}

test('Format Archive tab layout, -xmlhelp options, responsive viewports and -rid refusal', async ({ page }) => {
  const external = trackNet(page);
  await openFormatArchive(page);

  // Expand Formatting options so all parameter help tips are visible and checked.
  await page.getByTestId('bf-options-toggle').click();
  await expect(page.getByTestId('bf-parameter-form')).toBeVisible();

  // No search of this session wrote -outfmt 11 yet: the session source shows (0) and cannot be chosen.
  await expect(page.getByTestId('bf-source-session')).toBeDisabled();
  await expect(page.getByTestId('bf-session-count')).toHaveText('(0)');

  await assertLayoutAndHelpTips(page, DESKTOP_VIEWPORT);
  await maybeShot(page, 'bf-tab-1280.png');

  await assertLayoutAndHelpTips(page, MOBILE_VIEWPORT);
  await maybeShot(page, 'bf-tab-390.png');

  // -rid in Additional options is refused in-place by validateTool.
  await resetDesktopViewport(page);
  await page.getByTestId('bf-param-extra').fill('-rid ABCDEFGH01R');
  await expect(page.getByTestId('bf-argv-validation')).toHaveAttribute('data-state', 'invalid');
  await expect(page.getByTestId('bf-argv-validation')).toContainText("BLAST_FORMATTER");

  expectNoExternalRequests(external);
});

test('real engine: Format Archive matches native losat blast_formatter byte-for-byte (-subject, -db, outfmt 8/9/12, and Search -outfmt 11 shortcut)', async ({
  page,
}) => {
  test.skip(!BUILD_HAS_ENGINE, 'requires built Wasm reactors');
  test.skip(NATIVE === undefined, NO_NATIVE_REASON ?? 'requires native losat binary');
  const external = trackNet(page);
  await openFormatArchive(page);

  // 1. -subject archive (bn.asn) in every -outfmt that the tab offers (each exits 0), via the file picker.
  const bnPath = join(BF_DIR, 'in/bn.asn');
  await page.getByTestId('bf-archive-file-input').setInputFiles(bnPath);
  await expect(page.getByTestId('bf-archive-chosen')).toContainText('bn.asn');

  const offered = await page
    .getByTestId('bf-outfmt-select')
    .locator('option')
    .evaluateAll((options) => options.map((option) => (option as HTMLOptionElement).value).filter((value) => value !== 'other'));
  expect(offered.length, 'the -outfmt choices').toBeGreaterThan(0);
  let run = 0;
  for (const outfmt of offered) {
    await clickBfTab(page);
    await page.getByTestId('bf-outfmt-select').selectOption(outfmt);
    await runBfAndWait(page, ++run);
    await expect(page.getByTestId('tool-output-download')).toBeVisible();

    if (outfmt === '9') {
      await expect(page.getByTestId('tool-output-binary')).toBeVisible();
    } else {
      await expectTextOutputVisible(page);
    }
    const webBytes = await downloadActiveToolOutput(page);
    const native = runNativeBf(['-archive', 'in/bn.asn', '-outfmt', outfmt], BF_DIR);
    expect(native.status).toBe(0);
    expect(sha256(webBytes), `bn.asn outfmt ${outfmt}`).toBe(sha256(native.stdout));
  }

  // Capture results view screenshot at 1280 px and 390 px if requested.
  if (SCREENSHOT_DIR) {
    await maybeShot(page, 'bf-result-1280.png', DESKTOP_VIEWPORT);
    await maybeShot(page, 'bf-result-390.png', MOBILE_VIEWPORT);
    await resetDesktopViewport(page);
  }

  // Verify Support badge on Run details tab for completed blast_formatter run.
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('verification-badge')).toContainText('Supported');

  // 2. Local -db archive (db.asn + in/db/s3.*) with -out report.json.
  await clickBfTab(page);
  await page.getByTestId('bf-archive-file-input').setInputFiles(join(BF_DIR, 'in/db.asn'));
  const dbDir = join(BF_DIR, 'in/db');
  const dbFiles = readdirSync(dbDir)
    .sort()
    .map((name) => ({
      name: `in/db/${name}`,
      mimeType: 'application/octet-stream',
      buffer: readFileSync(join(dbDir, name)),
    }));
  // Inject relative-path files into the run-files input via DataTransfer so webkitRelativePath is preserved.
  await page.evaluate(
    (items) => {
      const input = document.querySelector<HTMLInputElement>('[data-testid="bf-run-files-input"]')!;
      const dt = new DataTransfer();
      for (const item of items) {
        const file = new File([new Uint8Array(item.bytes)], item.base);
        Object.defineProperty(file, 'webkitRelativePath', { value: item.rel });
        dt.items.add(file);
      }
      input.files = dt.files;
      input.dispatchEvent(new Event('change', { bubbles: true }));
    },
    dbFiles.map((f) => ({ rel: f.name, base: f.name.split('/').at(-1)!, bytes: [...f.buffer] })),
  );
  await expect(page.getByTestId('bf-run-files-list')).toContainText('in/db/s3.nin');
  await page.getByTestId('bf-outfmt-select').selectOption('12');
  await page.getByTestId('bf-out').fill('report.json');
  await runBfAndWait(page, ++run);
  await expect(page.getByTestId('tool-tab-file-0')).toContainText('-out report.json');
  const dbWebBytes = await downloadActiveToolOutput(page);
  const dbNative = runNativeBf(['-archive', 'in/db.asn', '-outfmt', '12'], BF_DIR);
  expect(dbNative.status).toBe(0);
  expect(sha256(dbWebBytes), 'db.asn outfmt 12').toBe(sha256(dbNative.stdout));

  // 3. End-to-end Search (-outfmt 11) -> "Format archive…" button -> Format Archive -> byte parity!
  await page.getByTestId('tab-search').click();
  await program(page, 'blastn');
  await paste(page, 'query', '>q1\nACGTTGCAAGGCTTAACCGGTTAACCGGATCGATCGATCG\n');
  await paste(page, 'subject', '>s1\nACGTTGCAAGGCTTAACCGGTTAACCGGATCGATCGATCG\n');
  await page.getByTestId('algorithm-parameters-toggle').click();
  await page.getByTestId('param-extra').fill('-outfmt 11');
  await submit(page);
  await waitStatus(page, ++run, 'completed');
  await page.getByTestId('results-view-outputs').click();
  await page.getByTestId('format-custom').click();
  await expect(page.getByTestId('format-archive-action')).toBeVisible();
  const archiveText = await page.getByTestId('result-output').textContent();
  expect(archiveText).toContain('Blast4-archive');
  // "Format archive…" beside Export fits a phone's width.
  await maybeShot(page, 'bf-outputs-1280.png', DESKTOP_VIEWPORT);
  await page.setViewportSize(MOBILE_VIEWPORT);
  expect(await pageOverflow(page), 'Outputs at 390px').toBeLessThanOrEqual(1);
  await maybeShot(page, 'bf-outputs-390.png');
  await resetDesktopViewport(page);

  // Click "Format archive…" -> opens Format Archive tab with session archive preselected.
  await page.getByTestId('format-archive-action').click();
  await expectBfPanelVisible(page);
  await expect(page.getByTestId('bf-source-session')).toBeChecked();
  // Clear the -out of step 2 so that the session archive is formatted to the standard output in outfmt 8.
  await page.getByTestId('bf-out').fill('');
  await page.getByTestId('bf-outfmt-select').selectOption('8');
  await runBfAndWait(page, ++run);
  await expectTextOutputVisible(page);
  const roundtripBytes = await downloadActiveToolOutput(page);
  const roundtripNative = runNativeBf(['-archive', '-', '-outfmt', '8'], BF_DIR, Buffer.from(archiveText ?? '', 'utf8'));
  expect(roundtripNative.status).toBe(0);
  expect(sha256(roundtripBytes), 'Search -outfmt 11 -> Format archive… -> outfmt 8').toBe(sha256(roundtripNative.stdout));

  expectNoExternalRequests(external);
});
