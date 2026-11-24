// End-to-end tests for saving browser-made WindowMasker counts under a taxonomy ID,
// ReplaceCountsDialog, CountsLine "Made in this browser", cross-tab BroadcastChannel sync and
// Web Lock busy refusal, -ustat/-in saved counts pickers, Supported / Not yet supported badges
// on tool runs, Download output/warnings button labels, inflated .gz session-file ref-counting
// when removed while queued or running, and the 64 MiB save page JS heap limit in Chromium.
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { closeSync, copyFileSync, mkdirSync, openSync, readFileSync, writeFileSync, writeSync } from 'node:fs';
import { join } from 'node:path';
import { gzipSync } from 'node:zlib';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE, ORIGIN } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, nativeExpectation, NO_NATIVE_REASON } from './support/native';
import { reload } from './support/navigate';
import { E2E_ORIGIN } from './support/origin';
import {
  chooseRun,
  described,
  openFiles,
  openParameters,
  openSearch,
  pickRun,
  program,
  storedOutputs,
  submit,
  task,
  waitStatus,
} from './support/search';

const COUNTS_FIXTURE = 'LOSAT/tests/winmask/counts/mg1655_u9.obinary';
const BLASTN_REG_DIR = join(REPOSITORY, 'LOSAT/tests/fixtures/blastn_regression/inputs');
const FLAG = 'window_masker_db';
const TAXID = 562;
const COUNTS_URL = `https://ftp.ncbi.nlm.nih.gov/blast/windowmasker_files/${TAXID}/wmasker.obinary`;
const LAST_MODIFIED = 'Thu, 01 Oct 2026 08:00:00 GMT';
const MIB = 1024 * 1024;

test.setTimeout(600_000);

const sha256 = (bytes: Uint8Array | string): string => createHash('sha256').update(bytes).digest('hex');

interface StubSite {
  gate: Promise<void> | undefined;
  readonly away: string[];
}

/** Stubs NCBI FTP for taxonomy ID 562 using the repository fixture (never contacts NCBI). */
async function stubCountsSite(page: Page): Promise<StubSite> {
  const fixtureBytes = readFileSync(join(REPOSITORY, COUNTS_FIXTURE));
  const site: StubSite = { gate: undefined, away: [] };
  page.on('request', (request) => {
    const url = request.url();
    if (/^https?:/.test(url) && !url.startsWith(E2E_ORIGIN)) site.away.push(url);
  });
  await page.route(/^https:\/\/ftp\.ncbi\.nlm\.nih\.gov\//, async (route) => {
    const cors = {
      'Access-Control-Allow-Origin': '*',
      'Access-Control-Expose-Headers': 'Content-Range, Content-Length, Last-Modified',
    };
    if (route.request().url() !== COUNTS_URL) {
      await route.fulfill({ status: 404, headers: cors, body: 'Not Found' });
      return;
    }
    const gate = site.gate;
    site.gate = undefined;
    await gate;
    await route.fulfill({
      status: 200,
      headers: {
        ...cors,
        'Content-Type': 'application/octet-stream',
        'Content-Length': String(fixtureBytes.length),
        'Last-Modified': LAST_MODIFIED,
      },
      body: fixtureBytes,
    });
  });
  return site;
}

function trackExternalNet(page: Page): string[] {
  const external: string[] = [];
  page.on('request', (req) => {
    const url = req.url();
    if (!url.startsWith(ORIGIN) && !url.startsWith('blob:') && !url.startsWith('data:')) external.push(url);
  });
  return external;
}

async function openWindowMasker(page: Page): Promise<void> {
  await openSearch(page);
  await page.getByTestId('tab-windowmasker').click();
  await expect(page.getByTestId('wm-panel')).toHaveAttribute('data-describe', 'ready', { timeout: 60_000 });
}

async function showWmTab(page: Page): Promise<void> {
  await page.getByTestId('tab-windowmasker').click();
  await expect(page.getByTestId('wm-panel')).toHaveAttribute('data-describe', 'ready', { timeout: 60_000 });
}

async function openWmOptions(page: Page): Promise<void> {
  const toggle = page.getByTestId('wm-options-toggle');
  if ((await toggle.getAttribute('aria-expanded')) !== 'true') await toggle.click();
  await expect(page.getByTestId('wm-parameter-form')).toBeVisible();
}

async function submitWm(page: Page): Promise<void> {
  await expect(page.getByTestId('wm-argv-validation')).not.toHaveAttribute('data-state', 'checking');
  await page.getByTestId('wm-run-losat').click();
}

async function clearWmFasta(page: Page): Promise<void> {
  const clearBtn = page.getByTestId('wm-fasta-clear');
  if ((await clearBtn.count()) > 0) await clearBtn.click();
  const list = page.locator('[data-testid="wm-fasta-list"] > li');
  while ((await list.count()) > 0) await page.getByTestId('wm-fasta-remove-0').click();
}

async function addWmFastaFile(page: Page, name: string, buffer: Buffer): Promise<void> {
  const before = await page.locator('[data-testid="wm-fasta-list"] > li').count();
  await page.getByTestId('wm-fasta-files').setInputFiles([{ name, mimeType: 'application/octet-stream', buffer }]);
  await expect(page.locator('[data-testid="wm-fasta-list"] > li').nth(before)).toHaveAttribute('data-status', 'ready', {
    timeout: 60_000,
  });
}

async function downloadToolActiveTab(page: Page): Promise<Buffer> {
  const [download] = await Promise.all([
    page.waitForEvent('download'),
    page.getByTestId('tool-output-download').click(),
  ]);
  return readFileSync((await download.path())!);
}

/** Deterministic nucleotide genome with repeated motifs so WindowMasker masks real intervals. */
function makeSyntheticGenome(targetBytes: number): Buffer {
  const alphabet = 'ACGT';
  const motif = 'ACGTTGCAAAGCTTAGGCTA'.repeat(4);
  const lines: string[] = ['>synth_chr1'];
  let written = lines[0]!.length + 1;
  let state = 0x13579bdf;
  let lineIndex = 0;
  while (written < targetBytes) {
    let line: string;
    if (lineIndex % 8 === 0) {
      line = motif;
    } else {
      line = '';
      for (let i = 0; i < 80; i++) {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        line += alphabet[(state >>> 30) & 3]!;
      }
    }
    lines.push(line);
    written += line.length + 1;
    lineIndex++;
  }
  return Buffer.from(`${lines.join('\n')}\n`);
}

function runNativeWindowMasker(
  args: readonly string[],
  cwd: string,
  outFile?: string,
): { readonly status: number; readonly stdout: Buffer; readonly stderr: Buffer; readonly file?: Buffer } {
  if (NATIVE === undefined) throw new Error(NO_NATIVE_REASON);
  const res = spawnSync(NATIVE, ['windowmasker', ...args], { cwd, maxBuffer: 1 << 30 });
  if (res.error) throw res.error;
  return {
    status: res.status ?? -1,
    stdout: res.stdout,
    stderr: res.stderr,
    ...(outFile === undefined ? {} : { file: readFileSync(join(cwd, outFile)) }),
  };
}

test('Save, Replace (over download and made copy, Escape cancels to Save), CountsLine after reload, saved picker, and Delete', async ({
  page,
}) => {
  const site = await stubCountsSite(page);
  await openSearch(page);
  const backend = await page.getByTestId('storage-status').getAttribute('data-backend');

  // 1) Download taxid 562 from stubbed FTP first so Replace over a download is tested.
  await program(page, 'blastn');
  await openParameters(page);
  await page.getByTestId(`param-${FLAG}`).check();
  await page.getByTestId(`param-${FLAG}-taxid`).fill(String(TAXID));
  await page.getByTestId('counts-download').click();
  await expect(page.getByTestId('counts-note')).toHaveText(/^(?:Saved in this browser|Kept in this tab only) \(2026-10-01, 568\.6[\s\u00a0]kB\)$/);
  expect(site.away).toEqual([COUNTS_URL]);

  // 2) Make counts in WindowMasker tab (run 1).
  await showWmTab(page);
  await page.getByTestId('wm-fasta-input').fill('>g1\nACGTTGCAAAGCTTAGGCTAACGTTGCAAAGCTTAGGCTAACGTTGCAAAGCTTAGGCTA\n');
  await submitWm(page);
  await waitStatus(page, 1, 'completed', 120_000);
  await pickRun(page, 1);
  await expect(page.getByTestId('tool-run-view')).toHaveAttribute('data-run', '1');

  // Check §4.6 Download warnings button label when stderr tab is present.
  if ((await page.getByTestId('tool-tab-stderr').count()) > 0) {
    await page.getByTestId('tool-tab-stderr').click();
    await expect(page.getByTestId('tool-output-download')).toContainText('Download warnings (');
  }

  // Open -out file tab and trigger Replace dialog over the downloaded copy.
  await page.getByTestId('tool-tab-file-0').click();
  await page.getByTestId('tool-save-taxid-input').fill(String(TAXID));
  await page.getByTestId('tool-save-taxid-button').click();
  const dialog = page.getByTestId('replace-counts-dialog');
  await expect(dialog).toBeVisible();
  await expect(page.getByTestId('replace-counts-message')).toHaveText(
    /^Counts for 562 are saved \(ftp\.ncbi\.nlm\.nih\.gov, 2026-10-01, 568\.6[\s\u00a0]kB\)\. Replace them with wmasker\.ascii of run 1\?$/,
  );

  // Escape cancels and returns focus to Save.
  await page.keyboard.press('Escape');
  await expect(dialog).not.toBeVisible();
  await expect(page.getByTestId('tool-save-taxid-button')).toBeFocused();

  // Save again and confirm Replace.
  await page.getByTestId('tool-save-taxid-button').click();
  await expect(dialog).toBeVisible();
  await page.getByTestId('replace-counts-confirm').click();
  await expect(dialog).not.toBeVisible();
  await expect(page.getByTestId('tool-save-status')).toHaveText(
    /^(?:Saved for taxonomy ID 562\.|Kept in this tab only for taxonomy ID 562\.)$/,
  );

  // 3) Make a second counts run (run 2) and Replace over the made copy.
  await showWmTab(page);
  await submitWm(page);
  await waitStatus(page, 2, 'completed', 120_000);
  await pickRun(page, 2);
  await expect(page.getByTestId('tool-run-view')).toHaveAttribute('data-run', '2');
  await page.getByTestId('tool-tab-file-0').click();
  await page.getByTestId('tool-save-taxid-input').fill(String(TAXID));
  await page.getByTestId('tool-save-taxid-button').click();
  await expect(dialog).toBeVisible();
  await expect(page.getByTestId('replace-counts-message')).toHaveText(
    /^Counts for 562 are saved \(made in this browser by run 1, \d{4}-\d{2}-\d{2}, .+\)\. Replace them with wmasker\.ascii of run 2\?$/,
  );
  await page.getByTestId('replace-counts-confirm').click();
  await expect(page.getByTestId('tool-save-status')).toHaveText(
    /^(?:Saved for taxonomy ID 562\.|Kept in this tab only for taxonomy ID 562\.)$/,
  );

  // 4) "Use in BLASTN" switches to -window_masker_taxid 562 and shows Made in this browser.
  await page.getByTestId('tool-use-in-blastn').click();
  await expect(page.getByTestId('tab-search')).toHaveAttribute('aria-selected', 'true');
  await expect(page.getByTestId(`param-${FLAG}`)).toBeChecked();
  await expect(page.getByTestId(`param-${FLAG}-by`)).toHaveValue('taxid');
  await expect(page.getByTestId(`param-${FLAG}-taxid`)).toHaveValue(String(TAXID));
  await expect(page.getByTestId('counts-note')).toHaveText(
    /^(?:Made in this browser|Made in this tab only) \(\d{4}-\d{2}-\d{2}, .+\)$/,
  );
  await expect(page.getByTestId('counts-download')).toHaveCount(0);
  await expect(page.getByTestId('counts-delete')).toBeVisible();

  // 5) After reload (OPFS), CountsLine says "Made in this browser" and -ustat / -in selects list it.
  if (backend === 'opfs') {
    await reload(page);
    await described(page);
    await program(page, 'blastn');
    await openParameters(page);
    await page.getByTestId(`param-${FLAG}`).check();
    await page.getByTestId(`param-${FLAG}-taxid`).fill(String(TAXID));
    await expect(page.getByTestId('counts-note')).toHaveText(/^Made in this browser \(\d{4}-\d{2}-\d{2}, .+\)$/);
    await expect(page.getByTestId('counts-download')).toHaveCount(0);
    expect(site.away).toEqual([COUNTS_URL]);

    await showWmTab(page);
    await page.getByTestId('wm-stage-mask').check();
    const ustatOption = page.locator('[data-testid="wm-ustat-select"] option[value="saved:562"]');
    await expect(ustatOption).toHaveCount(1);
    await expect(ustatOption).toContainText('windowmasker/562/wmasker.obinary (made in this browser,');

    await page.getByTestId('wm-stage-convert').check();
    const convertOption = page.locator('[data-testid="wm-convert-select"] option[value="saved:562"]');
    await expect(convertOption).toHaveCount(1);
    await expect(convertOption).toContainText('windowmasker/562/wmasker.obinary (made in this browser,');

    await page.getByTestId('tab-search').click();
  }

  // 6) Delete removes the made copy and brings Download back.
  await page.getByTestId('counts-delete').click();
  await expect(page.getByTestId('counts-note')).toHaveText('Downloads wmasker.obinary for 562 from ftp.ncbi.nlm.nih.gov');
  await expect(page.getByTestId('counts-download')).toBeVisible();
  expect(site.away).toEqual([COUNTS_URL]);
});

test('save is refused while another tab downloads the same ID, and cross-tab changes sync across pages', async ({
  context,
  page,
}) => {
  // Pause the streamed response body after the first 100 KB on page 1 so BrowserCountsDownloader
  // acquires the real OPFS Web Lock (`beginCounts`) and stays mid-download until released.
  await page.addInitScript((targetUrl: string) => {
    const origFetch = window.fetch.bind(window);
    const win = window as Window & { __releaseCountsBody?: () => void };
    window.fetch = async (input: RequestInfo | URL, init?: RequestInit) => {
      const res = await origFetch(input, init);
      const url = typeof input === 'string' ? input : input instanceof URL ? input.href : input.url;
      if (url !== targetUrl || res.body === null) return res;
      const bytes = new Uint8Array(await res.arrayBuffer());
      const gate = new Promise<void>((resolve) => {
        win.__releaseCountsBody = resolve;
      });
      const split = Math.min(100_000, bytes.byteLength);
      const stream = new ReadableStream<Uint8Array>({
        async start(controller) {
          controller.enqueue(bytes.subarray(0, split));
          await gate;
          if (split < bytes.byteLength) controller.enqueue(bytes.subarray(split));
          controller.close();
        },
      });
      return new Response(stream, { status: res.status, statusText: res.statusText, headers: res.headers });
    };
  }, COUNTS_URL);

  await stubCountsSite(page);
  await openSearch(page);
  const backend = await page.getByTestId('storage-status').getAttribute('data-backend');
  test.skip(backend !== 'opfs', 'cross-tab Web Locks and OPFS persistence require OPFS');

  // Page 1 starts downloading taxid 562, receives headers + first 100 KB, and holds the Web Lock open.
  await program(page, 'blastn');
  await openParameters(page);
  await page.getByTestId(`param-${FLAG}`).check();
  await page.getByTestId(`param-${FLAG}-taxid`).fill(String(TAXID));
  await page.getByTestId('counts-download').click();
  await expect(page.getByTestId('counts-note')).toHaveText(
    /^Downloading 100\.0[\s\u00a0]kB of 568\.6[\s\u00a0]kB from ftp\.ncbi\.nlm\.nih\.gov…$/,
  );

  const page2 = await context.newPage();
  try {
    await openWindowMasker(page2);
    await page2.getByTestId('wm-fasta-input').fill('>g1\nACGTTGCAAAGCTTAGGCTAACGTTGCAAAGCTTAGGCTA\n');
    await submitWm(page2);
    await waitStatus(page2, 1, 'completed', 120_000);
    await pickRun(page2, 1);
    await page2.getByTestId('tool-tab-file-0').click();
    await page2.getByTestId('tool-save-taxid-input').fill(String(TAXID));
    await page2.getByTestId('tool-save-taxid-button').click();

    // Refused in place while page 1 holds the download lock.
    await expect(page2.getByTestId('tool-save-error')).toHaveText(
      `Another tab is downloading the counts for ${TAXID}.`,
    );
    await expect(page2.getByTestId('replace-counts-dialog')).not.toBeVisible();

    // Release page 1's download body -> page 1 finishes.
    await page.evaluate(() => (window as Window & { __releaseCountsBody?: () => void }).__releaseCountsBody?.());
    await expect(page.getByTestId('counts-note')).toHaveText(/^Saved in this browser \(2026-10-01, 568\.6[\s\u00a0]kB\)$/);

    // Page 2 saves and replaces the download -> page 1 updates via BroadcastChannel "changed".
    await page2.getByTestId('tool-save-taxid-button').click();
    await expect(page2.getByTestId('replace-counts-dialog')).toBeVisible();
    await page2.getByTestId('replace-counts-confirm').click();
    await expect(page2.getByTestId('tool-save-status')).toHaveText('Saved for taxonomy ID 562.');
    await expect(page.getByTestId('counts-note')).toHaveText(/^Made in this browser \(\d{4}-\d{2}-\d{2}, .+\)$/, {
      timeout: 30_000,
    });

    // Page 1 deletes -> page 2's Mask saved-counts select drops saved:562.
    await showWmTab(page2);
    await page2.getByTestId('wm-stage-mask').check();
    await expect(page2.locator('[data-testid="wm-ustat-select"] option[value="saved:562"]')).toHaveCount(1);
    await page.getByTestId('counts-delete').click();
    await expect(page.getByTestId('counts-download')).toBeVisible();
    await expect(page2.locator('[data-testid="wm-ustat-select"] option[value="saved:562"]')).toHaveCount(0, {
      timeout: 30_000,
    });
  } finally {
    await page2.close();
  }
});

test.describe('real reactors and native CLI parity', () => {
  test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
  test.skip(NATIVE === undefined, NO_NATIVE_REASON);

  test('made counts (ascii, binary, obinary) saved for 562 match native BLASTN -window_masker_taxid before and after reload, and -ustat from saved copy matches native', async ({
    page,
  }, testInfo) => {
    const away = trackExternalNet(page);
    const nativeDir = testInfo.outputPath('native-taxid');
    const wmDir = join(nativeDir, 'wm', String(TAXID));
    mkdirSync(wmDir, { recursive: true });

    const genomeBytes = makeSyntheticGenome(40 * 1024);
    // Derive query & subject from the generated genome so hits exist and repeats get masked.
    const genomeText = genomeBytes.toString('utf8');
    const seqLines = genomeText.split('\n').slice(1).filter((l) => l.length > 0);
    const queryBytes = Buffer.from(`>q1\n${seqLines.slice(0, 12).join('\n')}\n`);
    const subjectBytes = Buffer.from(`>s1\n${seqLines.slice(0, 60).join('\n')}\n`);
    writeFileSync(join(nativeDir, 'genome.fa'), genomeBytes);
    writeFileSync(join(nativeDir, 'query.fa'), queryBytes);
    writeFileSync(join(nativeDir, 'subject.fa'), subjectBytes);

    await openWindowMasker(page);
    const backend = await page.getByTestId('storage-status').getAttribute('data-backend');
    await addWmFastaFile(page, 'genome.fa', genomeBytes);

    let runNo = 0;
    let lastCountsBytes!: Buffer;
    let lastBlastnExpect!: ReturnType<typeof nativeExpectation>;

    for (const fmt of ['ascii', 'binary', 'obinary'] as const) {
      await showWmTab(page);
      await page.getByTestId('wm-stage-mk_counts').check();
      await page.getByTestId('wm-param-sformat').selectOption(fmt);
      await expect(page.getByTestId('wm-out')).toHaveValue(`wmasker.${fmt}`);
      await submitWm(page);
      runNo++;
      const mkRun = runNo;
      await waitStatus(page, mkRun, 'completed', 180_000);
      await pickRun(page, mkRun);
      await page.getByTestId('tool-tab-file-0').click();
      const countsBytes = await downloadToolActiveTab(page);
      lastCountsBytes = countsBytes;

      // Verify counts bytes match native windowmasker -mk_counts -sformat <fmt>.
      const nativeMk = runNativeWindowMasker(
        ['-mk_counts', '-sformat', fmt, '-in', 'genome.fa', '-out', `wmasker.${fmt}`],
        nativeDir,
        `wmasker.${fmt}`,
      );
      expect(nativeMk.status).toBe(0);
      expect(sha256(countsBytes), `${fmt} counts sha256`).toBe(sha256(nativeMk.file!));

      // Save for taxonomy ID 562 (confirm Replace after first format).
      await page.getByTestId('tool-save-taxid-input').fill(String(TAXID));
      await page.getByTestId('tool-save-taxid-button').click();
      if (fmt !== 'ascii') {
        await expect(page.getByTestId('replace-counts-dialog')).toBeVisible();
        await page.getByTestId('replace-counts-confirm').click();
      }
      await expect(page.getByTestId('tool-save-status')).toHaveText(
        /^(?:Saved for taxonomy ID 562\.|Kept in this tab only for taxonomy ID 562\.)$/,
      );

      // Copy the same counts file into native WINDOW_MASKER_PATH/<taxid>/wmasker.obinary and run native BLASTN.
      writeFileSync(join(wmDir, 'wmasker.obinary'), countsBytes);
      const nativeBn = nativeExpectation(
        ['blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-task', 'blastn', '-window_masker_taxid', String(TAXID)],
        nativeDir,
        { WINDOW_MASKER_PATH: join(nativeDir, 'wm') },
      );
      lastBlastnExpect = nativeBn;

      // Use in BLASTN -> switches to -window_masker_taxid 562 -> run BLASTN -> compare outfmt 0, 6, 7 & stderr.
      await page.getByTestId('tool-use-in-blastn').click();
      await expect(page.getByTestId(`param-${FLAG}-by`)).toHaveValue('taxid');
      await expect(page.getByTestId(`param-${FLAG}-taxid`)).toHaveValue(String(TAXID));
      await expect(page.getByTestId('counts-note')).toHaveText(
        /^(?:Made in this browser|Made in this tab only) \(\d{4}-\d{2}-\d{2}, .+\)$/,
      );
      if (fmt === 'ascii') {
        await openFiles(page, 'query', [{ name: 'query.fa', text: queryBytes }]);
        await openFiles(page, 'subject', [{ name: 'subject.fa', text: subjectBytes }]);
        await task(page, 'blastn');
        await page.getByTestId('threads').selectOption('1');
      }
      await submit(page);
      runNo++;
      const bnRun = runNo;
      await waitStatus(page, bnRun, 'completed', 180_000);
      const { outputs, warnings } = await storedOutputs(page, bnRun);
      expect({ fmt, 0: sha256(outputs[0]), 6: sha256(outputs[6]), 7: sha256(outputs[7]) }).toEqual({
        fmt,
        0: nativeBn.sha256[0],
        6: nativeBn.sha256[6],
        7: nativeBn.sha256[7],
      });
      expect(sha256(warnings), `${fmt} BLASTN warnings`).toBe(nativeBn.stderrSha256);

      // Verify Run details show made origin without a URL line.
      await chooseRun(page, bnRun);
      await page.getByTestId('results-view-details').click();
      await expect(page.getByTestId('run-window-masker')).toContainText(
        `${TAXID}: wmasker.obinary made in this browser by run ${mkRun}`,
      );
      await expect(page.getByTestId('run-window-masker')).toContainText(`SHA-256 ${sha256(countsBytes)}`);
      await expect(page.getByTestId('run-window-masker-url')).toHaveCount(0);
    }

    // After reload (OPFS): CountsLine says "Made in this browser", BLASTN runs identically with 0 requests,
    // and WindowMasker Mask (-ustat from saved copy) equals native `losat windowmasker -ustat`.
    if (backend === 'opfs') {
      await reload(page);
      await described(page);
      await program(page, 'blastn');
      await openParameters(page);
      await page.getByTestId(`param-${FLAG}`).check();
      await page.getByTestId(`param-${FLAG}-taxid`).fill(String(TAXID));
      await expect(page.getByTestId('counts-note')).toHaveText(/^Made in this browser \(\d{4}-\d{2}-\d{2}, .+\)$/);
      await openFiles(page, 'query', [{ name: 'query.fa', text: queryBytes }]);
      await openFiles(page, 'subject', [{ name: 'subject.fa', text: subjectBytes }]);
      await task(page, 'blastn');
      await page.getByTestId('threads').selectOption('1');
      await submit(page);
      await waitStatus(page, 1, 'completed', 180_000);
      const afterReload = await storedOutputs(page, 1);
      expect({ 0: sha256(afterReload.outputs[0]), 6: sha256(afterReload.outputs[6]), 7: sha256(afterReload.outputs[7]) }).toEqual(
        lastBlastnExpect.sha256,
      );
      expect(sha256(afterReload.warnings)).toBe(lastBlastnExpect.stderrSha256);

      // Mask using saved taxonomy-ID copy in WindowMasker -ustat picker.
      await showWmTab(page);
      await page.getByTestId('wm-stage-mask').check();
      await page.getByTestId('wm-ustat-select').selectOption(`saved:${TAXID}`);
      await expect(page.getByTestId('wm-ustat-chosen')).toContainText(`windowmasker/${TAXID}/wmasker.obinary`);
      await addWmFastaFile(page, 'query.fa', queryBytes);
      await submitWm(page);
      await waitStatus(page, 2, 'completed', 180_000);
      await pickRun(page, 2);
      await page.getByTestId('results-view-outputs').click();
      await page.getByTestId('tool-tab-stdout').click();
      await expect(page.getByTestId('tool-output-download')).toContainText('Download output (');
      const maskStdout = await downloadToolActiveTab(page);

      const ustatDir = join(nativeDir, 'ustat-check');
      mkdirSync(join(ustatDir, 'windowmasker', String(TAXID)), { recursive: true });
      writeFileSync(join(ustatDir, 'windowmasker', String(TAXID), 'wmasker.obinary'), lastCountsBytes);
      writeFileSync(join(ustatDir, 'query.fa'), queryBytes);
      const nativeMask = runNativeWindowMasker(
        ['-ustat', `windowmasker/${TAXID}/wmasker.obinary`, '-in', 'query.fa'],
        ustatDir,
      );
      expect(nativeMask.status).toBe(0);
      expect(sha256(maskStdout), '-ustat from saved copy stdout').toBe(sha256(nativeMask.stdout));
    }

    expect(away, 'zero external network requests').toEqual([]);
  });

  test('Supported / Not yet supported badge on tool runs, and .gz input removed while queued or running completes identically to native', async ({
    page,
  }, testInfo) => {
    const away = trackExternalNet(page);
    await openWindowMasker(page);

    // 1) Supported badge for default -mk_counts (in cases.tsv).
    await page.getByTestId('wm-fasta-input').fill('>t\nACGTACGTACGTACGTACGT\n');
    await submitWm(page);
    await waitStatus(page, 1, 'completed', 120_000);
    await pickRun(page, 1);
    const badge = page.getByTestId('verification-badge');
    await expect(badge).toHaveText('Supported');
    await badge.click();
    await expect(page.getByTestId('verification-details')).toContainText('Supported');

    // 2) Not yet supported badge for an option set outside cases.tsv (-unit 11).
    await showWmTab(page);
    await openWmOptions(page);
    await page.getByTestId('wm-param-unit').fill('11');
    await submitWm(page);
    await waitStatus(page, 2, 'completed', 120_000);
    await pickRun(page, 2);
    await expect(page.getByTestId('verification-badge')).toHaveText('Not yet supported');

    // 3) §4.7 ①: .gz input removed from form while its run is queued behind a running run.
    const nativeDir = testInfo.outputPath('native-gz');
    mkdirSync(nativeDir, { recursive: true });
    const s60kPath = join(BLASTN_REG_DIR, 's60k.fa');
    copyFileSync(s60kPath, join(nativeDir, 's60k.fa'));
    const nativeS60k = runNativeWindowMasker(
      ['-mk_counts', '-in', 's60k.fa', '-out', 'wmasker.ascii'],
      nativeDir,
      'wmasker.ascii',
    );
    expect(nativeS60k.status).toBe(0);

    const longPath = testInfo.outputPath('long-blocker.fa');
    writeFileSync(longPath, makeSyntheticGenome(4 * MIB));
    await showWmTab(page);
    await clearWmFasta(page);
    await page.getByTestId('wm-restore-defaults').click();
    await openWmOptions(page);
    await page.getByTestId('wm-param-unit').fill('14');
    await page.getByTestId('wm-param-mem').fill('1');
    await page.getByTestId('wm-fasta-files').setInputFiles(longPath);
    await expect(page.locator('[data-testid="wm-fasta-list"] > li').first()).toHaveAttribute('data-status', 'ready', {
      timeout: 60_000,
    });
    await submitWm(page);
    await expect(page.getByTestId('run-3-option')).toHaveAttribute('data-status', 'running', { timeout: 60_000 });

    // Queue run 4 with s60k.fa.gz while run 3 occupies the Engine worker (submit directly without
    // waiting on background validateTool), then remove it from the form while queued.
    await showWmTab(page);
    await clearWmFasta(page);
    await page.getByTestId('wm-restore-defaults').click();
    const s60kGz = gzipSync(readFileSync(s60kPath));
    await addWmFastaFile(page, 's60k.fa.gz', s60kGz);
    await page.getByTestId('wm-run-losat').click();
    await expect(page.getByTestId('run-4-option')).toHaveAttribute('data-status', 'queued');

    await showWmTab(page);
    await page.getByTestId('wm-fasta-remove-0').click();
    await expect(page.locator('[data-testid="wm-fasta-list"] > li')).toHaveCount(0);

    // Cancel run 3; run 4 completes and matches native byte-for-byte.
    await pickRun(page, 3);
    await page.getByTestId('run-status-cancel').click();
    await waitStatus(page, 3, 'cancelled', 60_000);
    await waitStatus(page, 4, 'completed', 180_000);
    await pickRun(page, 4);
    await page.getByTestId('tool-tab-file-0').click();
    const queuedOut = await downloadToolActiveTab(page);
    expect(sha256(queuedOut), 'queued removed .gz output sha256').toBe(sha256(nativeS60k.file!));

    // 4) §4.7 ①: .gz input removed from form while its own run is running.
    // Multi-pass (-unit 12 -mem 1 on 512 KiB) runs ~3-5 s in Wasm so Remove executes mid-run while
    // windowmasker still makes additional passes over the OPFS session file.
    const medGenome = makeSyntheticGenome(512 * 1024);
    writeFileSync(join(nativeDir, 'med.fa'), medGenome);
    const nativeMed = runNativeWindowMasker(
      ['-mk_counts', '-unit', '12', '-mem', '1', '-in', 'med.fa', '-out', 'med.ascii'],
      nativeDir,
      'med.ascii',
    );
    expect(nativeMed.status).toBe(0);

    await showWmTab(page);
    await clearWmFasta(page);
    await openWmOptions(page);
    await page.getByTestId('wm-param-unit').fill('12');
    await page.getByTestId('wm-param-mem').fill('1');
    await page.getByTestId('wm-out').fill('med.ascii');
    await addWmFastaFile(page, 'med.fa.gz', gzipSync(medGenome));
    await submitWm(page);
    await expect(page.getByTestId('run-5-option')).toHaveAttribute('data-status', 'running', { timeout: 60_000 });
    // Remove the inflated .gz entry from the form while run 5 is actively running, then return to Results.
    await showWmTab(page);
    await page.getByTestId('wm-fasta-remove-0').click();
    await expect(page.locator('[data-testid="wm-fasta-list"] > li')).toHaveCount(0);
    await page.getByTestId('tab-results').click();

    await waitStatus(page, 5, 'completed', 240_000);
    await pickRun(page, 5);
    await page.getByTestId('tool-tab-file-0').click();
    const runningOut = await downloadToolActiveTab(page);
    expect(sha256(runningOut), 'running removed .gz output sha256').toBe(sha256(nativeMed.file!));

    expect(away, 'zero external network requests').toEqual([]);
  });

  test('saving a generated 64 MiB counts file grows page JS heap by < 32 MiB (Chromium)', async ({
    page,
    browserName,
  }, testInfo) => {
    test.skip(browserName !== 'chromium', 'performance.memory.usedJSHeapSize is Chromium-only');
    test.setTimeout(360_000);

    // Generate a valid 13-mer ASCII counts file (~28 MiB) that `-convert -sformat obinary` turns into ~69 MiB (> 64 MiB).
    const asciiPath = testInfo.outputPath('counts64.ascii');
    const fd = openSync(asciiPath, 'w');
    try {
      writeSync(fd, '13\n');
      const chunkLines: string[] = [];
      for (let i = 0; i < 2_900_000; i++) {
        chunkLines.push(`${i.toString(16)} ${25 + (i & 1)}\n`);
        if (chunkLines.length === 100_000) {
          writeSync(fd, chunkLines.join(''));
          chunkLines.length = 0;
        }
      }
      if (chunkLines.length > 0) writeSync(fd, chunkLines.join(''));
      writeSync(fd, '>t_low 10\n>t_extend 10\n>t_threshold 50\n>t_high 100\n');
    } finally {
      closeSync(fd);
    }

    await openWindowMasker(page);
    await page.getByTestId('wm-stage-convert').check();
    await page.getByTestId('wm-param-sformat').selectOption('obinary');
    await expect(page.getByTestId('wm-out')).toHaveValue('wmasker.obinary');
    await page.getByTestId('wm-convert-file').setInputFiles(asciiPath);
    await expect(page.getByTestId('wm-convert-chosen')).toContainText('counts64.ascii');
    await submitWm(page);
    await waitStatus(page, 1, 'completed', 240_000);
    await pickRun(page, 1);
    await page.getByTestId('tool-tab-file-0').click();
    await expect(page.getByTestId('tool-output-binary')).toContainText('MB');

    const readHeap = async (): Promise<number> => {
      const value = await page.evaluate(
        () => (performance as Performance & { memory?: { usedJSHeapSize: number } }).memory?.usedJSHeapSize ?? -1,
      );
      expect(value).toBeGreaterThan(0);
      return value;
    };

    const heapBefore = await readHeap();
    await page.evaluate(() => {
      const win = window as Window & { __savePeak?: number; __saveTimer?: number };
      const perf = performance as Performance & { memory?: { usedJSHeapSize: number } };
      win.__savePeak = perf.memory?.usedJSHeapSize ?? 0;
      win.__saveTimer = window.setInterval(() => {
        const cur = perf.memory?.usedJSHeapSize ?? 0;
        if (cur > (win.__savePeak ?? 0)) win.__savePeak = cur;
      }, 100);
    });

    await page.getByTestId('tool-save-taxid-input').fill(String(TAXID));
    await page.getByTestId('tool-save-taxid-button').click();
    await expect(page.getByTestId('tool-save-status')).toHaveText('Saved for taxonomy ID 562.', { timeout: 120_000 });

    const peakHeap = await page.evaluate(() => {
      const win = window as Window & { __savePeak?: number; __saveTimer?: number };
      if (win.__saveTimer !== undefined) window.clearInterval(win.__saveTimer);
      const perf = performance as Performance & { memory?: { usedJSHeapSize: number } };
      return Math.max(win.__savePeak ?? 0, perf.memory?.usedJSHeapSize ?? 0);
    });
    const heapAfter = await readHeap();
    const peakDelta = peakHeap - heapBefore;
    const endDelta = heapAfter - heapBefore;
    const toMiB = (b: number) => `${(b / MIB).toFixed(2)} MiB`;
    const binaryLabel = ((await page.getByTestId('tool-output-binary').textContent()) ?? '').trim();
    console.log(
      `[64 MiB save heap] output="${binaryLabel}" before=${toMiB(heapBefore)} peak=${toMiB(peakHeap)} after=${toMiB(heapAfter)} ` +
        `peakDelta=${toMiB(peakDelta)} endDelta=${toMiB(endDelta)}`,
    );

    expect(peakDelta, `peak JS heap growth ${toMiB(peakDelta)}`).toBeLessThan(32 * MIB);
    expect(endDelta, `end JS heap growth ${toMiB(endDelta)}`).toBeLessThan(32 * MIB);
  });
});
