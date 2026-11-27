// End-to-end tests for the Format Archive (`blast_formatter`) main tab (RG-F3p7):
// -xmlhelp loading, archive file / paste / session -outfmt 11 shortcut ("Format archive…"),
// local -db archive with run directory folder, byte-for-byte parity against native
// `losat blast_formatter` across -outfmt 8, 9 (binary), 12 and -out, validateTool refusal of -rid,
// responsive layout at 1280 px & 390 px (no orphaned `?` help tips, no horizontal overflow),
// and zero external network requests.
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { mkdirSync, readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE, ORIGIN } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { openSearch, paste, program, submit, waitStatus } from './support/search';

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

async function openFormatArchive(page: Page): Promise<void> {
  await openSearch(page);
  await page.getByTestId('tab-blast-formatter').click();
  await expect(page.getByTestId('bf-panel')).toBeVisible();
  await expect(page.getByTestId('bf-panel')).toHaveAttribute('data-describe', 'ready');
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

  await assertLayoutAndHelpTips(page, { width: 1280, height: 900 });
  if (SCREENSHOT_DIR) {
    mkdirSync(SCREENSHOT_DIR, { recursive: true });
    await page.screenshot({ path: join(SCREENSHOT_DIR, 'bf-tab-1280.png'), fullPage: true });
  }

  await assertLayoutAndHelpTips(page, { width: 390, height: 844 });
  if (SCREENSHOT_DIR) {
    await page.screenshot({ path: join(SCREENSHOT_DIR, 'bf-tab-390.png'), fullPage: true });
  }

  // -rid in Additional options is refused in-place by validateTool.
  await page.setViewportSize({ width: 1280, height: 900 });
  await page.getByTestId('bf-param-extra').fill('-rid ABCDEFGH01R');
  await expect(page.getByTestId('bf-argv-validation')).toHaveAttribute('data-state', 'invalid');
  await expect(page.getByTestId('bf-argv-validation')).toContainText("BLAST_FORMATTER");

  expect(external, 'zero external network requests').toEqual([]);
});

test('real engine: Format Archive matches native losat blast_formatter byte-for-byte (-subject, -db, outfmt 8/9/12, and Search -outfmt 11 shortcut)', async ({
  page,
}) => {
  test.skip(!BUILD_HAS_ENGINE, 'requires built Wasm reactors');
  test.skip(NATIVE === undefined, NO_NATIVE_REASON ?? 'requires native losat binary');
  const external = trackNet(page);
  await openFormatArchive(page);

  // 1. -subject archive (bn.asn) in outfmt 8, 9 (binary), 12 via file picker.
  const bnPath = join(BF_DIR, 'in/bn.asn');
  await page.getByTestId('bf-archive-file-input').setInputFiles(bnPath);
  await expect(page.getByTestId('bf-archive-chosen')).toContainText('bn.asn');

  for (const outfmt of ['8', '9', '12'] as const) {
    await page.getByTestId('tab-blast-formatter').click();
    await page.getByTestId('bf-outfmt-select').selectOption(outfmt);
    await page.getByTestId('bf-run-losat').click();
    await expect(page.getByTestId('tool-run-view')).toBeVisible();
    await expect(page.getByTestId('tool-output-download')).toBeVisible();

    if (outfmt === '9') {
      await expect(page.getByTestId('tool-output-binary')).toBeVisible();
    } else {
      await expect(page.getByTestId('tool-output-text')).toBeVisible();
    }
    const webBytes = await downloadActiveToolOutput(page);
    const native = runNativeBf(['-archive', 'in/bn.asn', '-outfmt', outfmt], BF_DIR);
    expect(native.status).toBe(0);
    expect(sha256(webBytes), `bn.asn outfmt ${outfmt}`).toBe(sha256(native.stdout));
  }

  // Capture results view screenshot at 1280 px and 390 px if requested.
  if (SCREENSHOT_DIR) {
    mkdirSync(SCREENSHOT_DIR, { recursive: true });
    await page.setViewportSize({ width: 1280, height: 900 });
    await page.screenshot({ path: join(SCREENSHOT_DIR, 'bf-result-1280.png'), fullPage: true });
    await page.setViewportSize({ width: 390, height: 844 });
    await page.screenshot({ path: join(SCREENSHOT_DIR, 'bf-result-390.png'), fullPage: true });
    await page.setViewportSize({ width: 1280, height: 900 });
  }

  // Verify Support badge on Run details tab for completed blast_formatter run.
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('verification-badge')).toContainText('Supported');

  // 2. Local -db archive (db.asn + in/db/s3.*) with -out report.json.
  await page.getByTestId('tab-blast-formatter').click();
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
  await page.getByTestId('bf-run-losat').click();
  await expect(page.getByTestId('tool-run-view')).toBeVisible();
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
  await waitStatus(page, 5, 'completed');
  await page.getByTestId('results-view-outputs').click();
  await page.getByTestId('format-custom').click();
  await expect(page.getByTestId('format-archive-action')).toBeVisible();
  const archiveText = await page.getByTestId('result-output').textContent();
  expect(archiveText).toContain('Blast4-archive');

  // Click "Format archive…" -> opens Format Archive tab with session archive preselected.
  await page.getByTestId('format-archive-action').click();
  await expect(page.getByTestId('bf-panel')).toBeVisible();
  await expect(page.getByTestId('bf-source-session')).toBeChecked();
  // Clear previous -out and run files from step 2 so we format the session archive to stdout in outfmt 8.
  await page.getByTestId('bf-out').fill('');
  await page.getByTestId('bf-outfmt-select').selectOption('8');
  await page.getByTestId('bf-run-losat').click();
  await expect(page.getByTestId('tool-run-view')).toBeVisible();
  await expect(page.getByTestId('tool-output-text')).toBeVisible();
  const roundtripBytes = await downloadActiveToolOutput(page);
  const roundtripNative = runNativeBf(['-archive', '-', '-outfmt', '8'], BF_DIR, Buffer.from(archiveText ?? '', 'utf8'));
  expect(roundtripNative.status).toBe(0);
  expect(sha256(roundtripBytes), 'Search -outfmt 11 -> Format archive… -> outfmt 8').toBe(sha256(roundtripNative.stdout));

  expect(external, 'zero external network requests').toEqual([]);
});
