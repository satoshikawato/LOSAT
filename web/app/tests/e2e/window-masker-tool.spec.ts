// End-to-end tests for the WindowMasker main tab: stage switching, -xmlhelp option filtering,
// BLASTN "Make counts…" shortcut, validateTool in-place refusals, M5 parity rows against native
// `losat windowmasker` and frozen NCBI bytes, "Use in BLASTN" / "Mask with these counts",
// compressed .gz inputs, cancel & shared queue ordering, -h/-version-full, -out/-logfile sharing,
// zero external network requests, and the 64 MiB page JS heap limit in Chromium.
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { copyFileSync, mkdirSync, readdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { gzipSync } from 'node:zlib';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE, ORIGIN } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, nativeExpectation, NO_NATIVE_REASON } from './support/native';
import {
  chooseRun,
  fasta,
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

const WM_DIR = join(REPOSITORY, 'LOSAT/tests/windowmasker');
const WINMASK_DIR = join(REPOSITORY, 'LOSAT/tests/winmask');
const BLASTN_REG_DIR = join(REPOSITORY, 'LOSAT/tests/fixtures/blastn_regression/inputs');
const MIB = 1024 * 1024;

test.setTimeout(600_000);

const sha256 = (bytes: Uint8Array | string): string => createHash('sha256').update(bytes).digest('hex');

interface ToolRunBytes {
  readonly status: 'completed' | 'failed';
  readonly exitStatus: number;
  readonly stdout: Buffer;
  readonly stderr: Buffer;
  readonly files: ReadonlyMap<string, { readonly label: string; readonly bytes: Buffer }>;
}

interface FrozenNcbiRow {
  readonly rc: number;
  readonly stdout: Buffer;
  readonly stderr: Buffer;
  readonly files: ReadonlyMap<string, Buffer>;
}

/** Reads frozen NCBI BLAST+ 2.17.0 windowmasker outputs from LOSAT/tests/windowmasker/ncbi/. */
function readFrozenNcbi(row: string): FrozenNcbiRow {
  const base = join(WM_DIR, 'ncbi', row);
  const rc = Number(readFileSync(`${base}.rc`, 'utf8').trim());
  const stderr = readFileSync(`${base}.err`);
  const rawOut = readFileSync(`${base}.out`);
  const marker = Buffer.from('\n@OUT@/');
  const firstAt = rawOut.indexOf(marker);
  if (firstAt < 0) {
    return { rc, stdout: rawOut, stderr, files: new Map() };
  }
  const stdout = rawOut.subarray(0, firstAt);
  const frozenFiles = new Map<string, Buffer>();
  let pos = firstAt;
  while (pos < rawOut.length) {
    if (!rawOut.subarray(pos, pos + marker.length).equals(marker)) {
      throw new Error(`unexpected @OUT@ framing in ${base}.out at offset ${pos}`);
    }
    const headerStart = pos + marker.length;
    const headerEnd = rawOut.indexOf(0x0a, headerStart);
    const header = rawOut.subarray(headerStart, headerEnd).toString('utf8');
    const space = header.lastIndexOf(' ');
    const name = header.slice(0, space);
    const length = Number(header.slice(space + 1));
    const bodyStart = headerEnd + 1;
    frozenFiles.set(name, Buffer.from(rawOut.subarray(bodyStart, bodyStart + length)));
    pos = bodyStart + length;
  }
  return { rc, stdout, stderr, files: frozenFiles };
}

/** Runs native `losat windowmasker` in a staged directory holding `inputs`. */
function runNativeTool(
  args: readonly string[],
  cwd: string,
  expectedFiles: readonly string[] = [],
): { readonly status: number; readonly stdout: Buffer; readonly stderr: Buffer; readonly files: ReadonlyMap<string, Buffer> } {
  if (NATIVE === undefined) throw new Error(NO_NATIVE_REASON);
  const result = spawnSync(NATIVE, ['windowmasker', ...args], { cwd, maxBuffer: 1 << 30 });
  if (result.error) throw result.error;
  const files = new Map<string, Buffer>();
  for (const name of expectedFiles) {
    files.set(name, readFileSync(join(cwd, name)));
  }
  return { status: result.status ?? -1, stdout: result.stdout, stderr: result.stderr, files };
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

function assertNoNet(external: readonly string[]): void {
  expect(external, 'zero external network requests').toEqual([]);
}

function ensureDir(dir: string): string {
  mkdirSync(dir, { recursive: true });
  return dir;
}

async function showWmTab(page: Page): Promise<void> {
  await page.getByTestId('tab-windowmasker').click();
}

async function expectWmTabSelected(page: Page): Promise<void> {
  await expect(page.getByTestId('tab-windowmasker')).toHaveAttribute('aria-selected', 'true');
}

async function expectWmReady(page: Page): Promise<void> {
  await expect(page.getByTestId('wm-panel')).toHaveAttribute('data-describe', 'ready', { timeout: 60_000 });
}

async function waitValidationSettled(page: Page): Promise<void> {
  await expect(page.getByTestId('wm-argv-validation')).not.toHaveAttribute('data-state', 'checking');
}

async function validIs(page: Page, state: 'ok' | 'invalid'): Promise<void> {
  await expect(page.getByTestId('wm-argv-validation')).toHaveAttribute('data-state', state, { timeout: 30_000 });
}

/** Opens the WindowMasker main tab and waits for `-xmlhelp` description to load. */
async function openWindowMasker(page: Page): Promise<void> {
  await openSearch(page);
  await showWmTab(page);
  await expectWmReady(page);
}

async function openWmOptions(page: Page): Promise<void> {
  const toggle = page.getByTestId('wm-options-toggle');
  if ((await toggle.getAttribute('aria-expanded')) !== 'true') await toggle.click();
  await expect(page.getByTestId('wm-parameter-form')).toBeVisible();
}

async function selectWmStage(page: Page, stage: 'mk_counts' | 'mask' | 'convert'): Promise<void> {
  await page.getByTestId(`wm-stage-${stage}`).check();
  await waitValidationSettled(page);
}

async function setSformat(page: Page, value: 'ascii' | 'binary' | 'oascii' | 'obinary'): Promise<void> {
  await page.getByTestId('wm-param-sformat').selectOption(value);
}

async function setOutName(page: Page, value: string): Promise<void> {
  await page.getByTestId('wm-out').fill(value);
}

async function setExtra(page: Page, value: string): Promise<void> {
  await page.getByTestId('wm-param-extra').fill(value);
}

async function restoreDefaults(page: Page): Promise<void> {
  await page.getByTestId('wm-restore-defaults').click();
}

async function noField(page: Page, name: string): Promise<void> {
  await expect(page.getByTestId(`wm-param-${name}-field`)).toHaveCount(0);
}

async function outIs(page: Page, value: string): Promise<void> {
  await expect(page.getByTestId('wm-out')).toHaveValue(value);
}

async function pickCounts(page: Page, kind: 'convert' | 'ustat', name: string, filePath: string): Promise<void> {
  await page.getByTestId(`wm-${kind}-file`).setInputFiles({
    name,
    mimeType: 'application/octet-stream',
    buffer: readFileSync(filePath),
  });
}

async function useDb(page: Page, folder: string, expectedName: string): Promise<void> {
  await page.getByTestId('wm-param-infmt').selectOption('blastdb');
  await page.getByTestId('wm-in-folder').setInputFiles(folder);
  await expect(page.getByTestId('wm-in-name')).toHaveValue(expectedName);
}

async function openOutputTab(page: Page, which: 'stdout' | 'file-0'): Promise<void> {
  await page.getByTestId(`tool-tab-${which}`).click();
}

async function expectBinaryNotice(page: Page): Promise<void> {
  await expect(page.getByTestId('tool-output-binary')).toContainText('Binary,');
}

async function clickMaskWithCounts(page: Page): Promise<void> {
  await page.getByTestId('tool-mask-with-counts').click();
}

async function clickCancelRun(page: Page): Promise<void> {
  await page.getByTestId('run-status-cancel').click();
}

async function setFastaFilePath(page: Page, filePath: string): Promise<void> {
  await page.getByTestId('wm-fasta-files').setInputFiles(filePath);
  await expect(page.locator('[data-testid="wm-fasta-list"] > li').first()).toHaveAttribute('data-status', 'ready', {
    timeout: 60_000,
  });
}

async function clearWmFasta(page: Page): Promise<void> {
  const clearBtn = page.getByTestId('wm-fasta-clear');
  if ((await clearBtn.count()) > 0) await clearBtn.click();
  const list = page.locator('[data-testid="wm-fasta-list"] > li');
  while ((await list.count()) > 0) await page.getByTestId('wm-fasta-remove-0').click();
}

async function clearWmRunFiles(page: Page): Promise<void> {
  await openWmOptions(page);
  const items = page.locator('[data-testid="wm-run-files-list"] > li');
  while ((await items.count()) > 0) await page.getByTestId('wm-run-file-remove-0').click();
}

async function addWmFastaFiles(page: Page, files: ReadonlyArray<{ name: string; buffer: Buffer }>): Promise<void> {
  const before = await page.locator('[data-testid="wm-fasta-list"] > li').count();
  await page.getByTestId('wm-fasta-files').setInputFiles(
    files.map((file) => ({ name: file.name, mimeType: 'application/octet-stream', buffer: file.buffer })),
  );
  for (let i = 0; i < files.length; i++) {
    await expect(page.locator('[data-testid="wm-fasta-list"] > li').nth(before + i)).toHaveAttribute('data-status', 'ready', {
      timeout: 60_000,
    });
  }
}

async function submitWm(page: Page): Promise<void> {
  await waitValidationSettled(page);
  await page.getByTestId('wm-run-losat').click();
}

/** Downloads every tab of a completed or failed tool run from Results. */
async function collectToolOutputs(page: Page, runNumber: number, expectedStatus: 'completed' | 'failed' = 'completed'): Promise<ToolRunBytes> {
  await waitStatus(page, runNumber, expectedStatus, 300_000);
  await pickRun(page, runNumber);
  await expect(page.getByTestId('tool-run-view')).toHaveAttribute('data-run', String(runNumber), { timeout: 60_000 });
  await page.getByTestId('results-view-outputs').click();

  let exitStatus = 0;
  if (expectedStatus === 'failed') {
    const text = ((await page.getByTestId('tool-exit-status').textContent()) ?? '').trim();
    const match = /Exit status (\d+)/.exec(text);
    exitStatus = match === null ? -1 : Number(match[1]);
  }

  const downloadActive = async (): Promise<Buffer> => {
    const [download] = await Promise.all([page.waitForEvent('download'), page.getByTestId('tool-output-download').click()]);
    return readFileSync((await download.path())!);
  };

  let stdout: Buffer = Buffer.alloc(0);
  if ((await page.getByTestId('tool-tab-stdout').count()) > 0) {
    await openOutputTab(page, 'stdout');
    stdout = await downloadActive();
  }

  const files = new Map<string, { readonly label: string; readonly bytes: Buffer }>();
  const fileTabs = page.locator('[data-testid^="tool-tab-file-"]');
  const fileCount = await fileTabs.count();
  for (let i = 0; i < fileCount; i++) {
    const tab = page.getByTestId(`tool-tab-file-${i}`);
    const label = ((await tab.textContent()) ?? '').trim();
    await tab.click();
    const bytes = await downloadActive();
    const fileName = label.replace(/^(?:-out,\s*-logfile|-out|-logfile)\s+/, '');
    files.set(fileName, { label, bytes });
  }

  let stderr: Buffer = Buffer.alloc(0);
  if ((await page.getByTestId('tool-tab-stderr').count()) > 0) {
    await page.getByTestId('tool-tab-stderr').click();
    stderr = await downloadActive();
  }

  return { status: expectedStatus, exitStatus, stdout, stderr, files };
}

/** Copies a BLAST database prefix (`wm.*`) from LOSAT/tests/windowmasker/db/ into `<dir>/db/`. */
function copyBlastDbFolder(targetParent: string, prefix: string): string {
  const dbDir = ensureDir(join(targetParent, 'db'));
  const sourceDir = join(WM_DIR, 'db');
  for (const entry of readdirSync(sourceDir)) {
    if (entry.startsWith(`${prefix}.`)) {
      copyFileSync(join(sourceDir, entry), join(dbDir, entry));
    }
  }
  return dbDir;
}

/** Writes a deterministic nucleotide FASTA file of at least `targetBytes` bytes. */
function writeSyntheticFasta(filePath: string, targetBytes: number): number {
  const alphabet = 'ACGT';
  const lineLen = 80;
  const motif = 'ACGTTGCAAAGCTTAGGCTA'.repeat(4);
  const chunks: Buffer[] = [];
  let written = 0;
  let seqIdx = 0;
  let state = 0x12345678;
  const recordTarget = 1024 * 1024;
  while (written < targetBytes) {
    seqIdx++;
    const header = Buffer.from(`>seq_${seqIdx}\n`);
    chunks.push(header);
    written += header.length;
    let recBytes = 0;
    const lines: string[] = [];
    while (recBytes < recordTarget && written + recBytes < targetBytes) {
      if ((recBytes / (lineLen + 1)) % 32 === 0) {
        lines.push(motif);
      } else {
        let line = '';
        for (let i = 0; i < lineLen; i++) {
          state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
          line += alphabet[(state >>> 30) & 3]!;
        }
        lines.push(line);
      }
      recBytes += lineLen + 1;
    }
    const body = Buffer.from(`${lines.join('\n')}\n`);
    chunks.push(body);
    written += body.length;
  }
  const full = Buffer.concat(chunks);
  writeFileSync(filePath, full);
  return full.length;
}

// --- UI behavior (runs in both FakeEngine and real engine builds) --------------------------

test('WindowMasker tab: stages, -xmlhelp option filtering, -out default, BLASTN shortcut, duplicate refusal, and Back history', async ({
  page,
}) => {
  await openSearch(page);
  await program(page, 'blastn');
  await openParameters(page);

  // BLASTN repeats row: "Make counts…" appears only when "Counts file" is selected.
  await page.getByTestId('param-window_masker_db').check();
  await expect(page.getByTestId('param-window_masker_db-make-counts')).toHaveCount(0);
  await page.getByTestId('param-window_masker_db-by').selectOption('file');
  await expect(page.getByTestId('param-window_masker_db-make-counts')).toBeVisible();
  await page.getByTestId('param-window_masker_db-make-counts').click();

  // Lands on WindowMasker tab, Make counts stage, with -xmlhelp loaded.
  await expectWmTabSelected(page);
  await expectWmReady(page);
  await expect(page.getByTestId('wm-stage-mk_counts')).toBeChecked();

  // Default -out follows -sformat until edited; '-' sends to stdout.
  await outIs(page, 'wmasker.ascii');
  await setSformat(page, 'obinary');
  await outIs(page, 'wmasker.obinary');
  await expect(page.getByTestId('wm-command')).toContainText('-out wmasker.obinary');
  await setOutName(page, 'custom.counts');
  await setSformat(page, 'binary');
  await outIs(page, 'custom.counts');
  await setOutName(page, '-');
  await expect(page.getByTestId('wm-command')).not.toContainText('-out ');

  // Stage option filtering driven by -xmlhelp exclusions.
  await openWmOptions(page);
  await expect(page.getByTestId('wm-param-checkdup-field')).toBeVisible();
  await expect(page.getByTestId('wm-param-unit-field')).toBeVisible();
  await noField(page, 'dust');
  await noField(page, 'window');

  // Duplicate file refusal in Files of the run.
  const dup = { name: 'ids.txt', mimeType: 'text/plain', buffer: Buffer.from('seq1\n') };
  const runFiles = page.getByTestId('wm-run-files-input');
  await runFiles.setInputFiles(dup);
  await expect(page.getByTestId('wm-run-file-name-0')).toHaveText('ids.txt');
  await runFiles.setInputFiles(dup);
  await expect(page.getByTestId('wm-run-files-error')).toContainText('already in the run directory');
  await page.getByTestId('wm-run-file-remove-0').click();

  // Mask stage options.
  await selectWmStage(page, 'mask');
  await expect(page.getByTestId('wm-param-dust-field')).toBeVisible();
  await expect(page.getByTestId('wm-param-window-field')).toBeVisible();
  await noField(page, 'checkdup');
  await noField(page, 'unit');

  // Convert counts stage options.
  await selectWmStage(page, 'convert');
  await expect(page.getByTestId('wm-param-meta-field')).toBeVisible();
  await expect(page.getByTestId('wm-param-smem-field')).toBeVisible();
  await noField(page, 'dust');
  await noField(page, 'checkdup');

  // Argument error refused in place by validateTool: -unit 0 and -dust T on Make counts.
  await selectWmStage(page, 'mk_counts');
  await page.getByTestId('wm-param-unit').fill('0');
  await validIs(page, 'invalid');
  await expect(page.getByTestId('wm-param-unit')).toHaveAttribute('aria-invalid', 'true');
  await page.getByTestId('wm-fasta-input').fill('>q\nACGTACGTACGT\n');
  await submitWm(page);
  await expect(page.getByTestId('run-1-option')).toHaveCount(0);

  await page.getByTestId('wm-param-unit').fill('');
  await setExtra(page, '-dust T');
  await validIs(page, 'invalid');
  await expect(page.getByTestId('wm-param-extra')).toHaveAttribute('aria-invalid', 'true');
  await setExtra(page, '');
  await validIs(page, 'ok');

  // Submitting opens Results; browser Back returns to WindowMasker.
  await setOutName(page, 'wmasker.ascii');
  await setSformat(page, 'ascii');
  await submitWm(page);
  await waitStatus(page, 1, 'completed');
  await expect(page.getByTestId('tab-results')).toHaveAttribute('aria-selected', 'true');
  await page.goBack();
  await expectWmTabSelected(page);
});

// --- Real engine + native CLI parity tests -------------------------------------------------

test.describe('WindowMasker with real reactors and native CLI', () => {
  test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
  test.skip(NATIVE === undefined, NO_NATIVE_REASON);

  test('M5 -mk_counts (4 sformats, -checkdup, -infmt blastdb via Open folder) and -convert match native CLI and frozen NCBI bytes', async ({
    page,
  }, testInfo) => {
    const net = trackNet(page);

    await openWindowMasker(page);
    await openWmOptions(page);

    const dm60kPath = join(WINMASK_DIR, 'fasta/dm2L_60k.fa');
    const dm60kBuf = readFileSync(dm60kPath);
    const nativeDir = ensureDir(testInfo.outputPath('native-mk'));
    copyFileSync(dm60kPath, join(nativeDir, 'dm2L_60k.fa'));

    let runNum = 0;

    // 1) -mk_counts across all 4 -sformat values (file output c.counts -> k_*_file rows).
    await addWmFastaFiles(page, [{ name: 'dm2L_60k.fa', buffer: dm60kBuf }]);
    await page.getByTestId('wm-param-unit').fill('5');
    await setOutName(page, 'c.counts');

    for (const [sformat, row] of [
      ['ascii', 'k_ascii_file'],
      ['binary', 'k_binary_file'],
      ['oascii', 'k_oascii_file'],
      ['obinary', 'k_obinary_file'],
    ] as const) {
      await showWmTab(page);
      await setSformat(page, sformat);
      await submitWm(page);
      runNum++;
      const actual = await collectToolOutputs(page, runNum);
      const native = runNativeTool(
        ['-mk_counts', '-in', 'dm2L_60k.fa', '-unit', '5', '-sformat', sformat, '-out', 'c.counts'],
        nativeDir,
        ['c.counts'],
      );
      const frozen = readFrozenNcbi(row);
      expect(actual.exitStatus, row).toBe(frozen.rc);
      expect(sha256(actual.stdout), `${row} stdout`).toBe(sha256(native.stdout));
      expect(sha256(actual.stdout), `${row} frozen stdout`).toBe(sha256(frozen.stdout));
      expect(sha256(actual.stderr), `${row} stderr`).toBe(sha256(native.stderr));
      expect(sha256(actual.stderr), `${row} frozen stderr`).toBe(sha256(frozen.stderr));
      const actualFile = actual.files.get('c.counts');
      expect(actualFile, `${row} c.counts present`).toBeDefined();
      expect(sha256(actualFile!.bytes), `${row} vs native`).toBe(sha256(native.files.get('c.counts')!));
      expect(sha256(actualFile!.bytes), `${row} vs NCBI`).toBe(sha256(frozen.files.get('c.counts')!));
      if (sformat === 'binary' || sformat === 'obinary') {
        await openOutputTab(page, 'file-0');
        await expectBinaryNotice(page);
      }
    }

    // 2) -mk_counts -checkdup T (cd_one row, stdout via -out -).
    await showWmTab(page);
    await clearWmFasta(page);
    const cdBuf = readFileSync(join(WM_DIR, 'in/cd.fa'));
    copyFileSync(join(WM_DIR, 'in/cd.fa'), join(nativeDir, 'cd.fa'));
    await addWmFastaFiles(page, [{ name: 'cd.fa', buffer: cdBuf }]);
    await page.getByTestId('wm-param-checkdup').selectOption('true');
    await page.getByTestId('wm-param-unit').fill('3');
    await setSformat(page, 'binary');
    await setOutName(page, '-');
    await submitWm(page);
    runNum++;
    const cdActual = await collectToolOutputs(page, runNum);
    const cdNative = runNativeTool(['-mk_counts', '-in', 'cd.fa', '-checkdup', 'T', '-unit', '3', '-sformat', 'binary'], nativeDir);
    const cdFrozen = readFrozenNcbi('cd_one');
    expect(sha256(cdActual.stdout), 'cd_one vs native').toBe(sha256(cdNative.stdout));
    expect(sha256(cdActual.stdout), 'cd_one vs NCBI').toBe(sha256(cdFrozen.stdout));
    expect(sha256(cdActual.stderr), 'cd_one stderr').toBe(sha256(cdFrozen.stderr));

    // 3) -mk_counts -infmt blastdb via Open folder… (k_blastdb row).
    await showWmTab(page);
    await restoreDefaults(page);
    const dbFolder = copyBlastDbFolder(testInfo.outputPath('blastdb-mk'), 'wm');
    copyBlastDbFolder(nativeDir, 'wm');
    await useDb(page, dbFolder, 'db/wm');
    await setSformat(page, 'binary');
    await setOutName(page, '-');
    await submitWm(page);
    runNum++;
    const dbActual = await collectToolOutputs(page, runNum);
    const dbNative = runNativeTool(['-mk_counts', '-infmt', 'blastdb', '-in', 'db/wm', '-sformat', 'binary'], nativeDir);
    const dbFrozen = readFrozenNcbi('k_blastdb');
    expect(sha256(dbActual.stdout), 'k_blastdb vs native').toBe(sha256(dbNative.stdout));
    expect(sha256(dbActual.stdout), 'k_blastdb vs NCBI').toBe(sha256(dbFrozen.stdout));
    expect(sha256(dbActual.stderr), 'k_blastdb stderr').toBe(sha256(dbFrozen.stderr));

    // 4) -convert: v_obinary_to_binary (file), v_ascii_to_obinary (stdout), v_meta_new (stdout + -meta).
    await showWmTab(page);
    await selectWmStage(page, 'convert');
    const u4ObPath = join(WM_DIR, 'counts/u4.obinary');
    const u4AscPath = join(WM_DIR, 'counts/u4.ascii');
    copyFileSync(u4ObPath, join(nativeDir, 'u4.obinary'));
    copyFileSync(u4AscPath, join(nativeDir, 'u4.ascii'));

    const pickObCounts = () => pickCounts(page, 'convert', 'u4.obinary', u4ObPath);
    await pickObCounts();
    await setSformat(page, 'binary');
    await setOutName(page, 'c.counts');
    await submitWm(page);
    runNum++;
    const convFileActual = await collectToolOutputs(page, runNum);
    const convFileNative = runNativeTool(['-convert', '-in', 'u4.obinary', '-sformat', 'binary', '-out', 'c.counts'], nativeDir, ['c.counts']);
    const convFileFrozen = readFrozenNcbi('v_obinary_to_binary');
    expect(sha256(convFileActual.files.get('c.counts')!.bytes), 'v_obinary_to_binary vs native').toBe(
      sha256(convFileNative.files.get('c.counts')!),
    );
    expect(sha256(convFileActual.files.get('c.counts')!.bytes), 'v_obinary_to_binary vs NCBI').toBe(
      sha256(convFileFrozen.files.get('c.counts')!),
    );

    await showWmTab(page);
    await pickCounts(page, 'convert', 'u4.ascii', u4AscPath);
    await setSformat(page, 'obinary');
    await setOutName(page, '-');
    await submitWm(page);
    runNum++;
    const convStdActual = await collectToolOutputs(page, runNum);
    const convStdFrozen = readFrozenNcbi('v_ascii_to_obinary');
    expect(sha256(convStdActual.stdout), 'v_ascii_to_obinary vs NCBI').toBe(sha256(convStdFrozen.stdout));

    await showWmTab(page);
    await pickObCounts();
    await setSformat(page, 'ascii');
    await page.getByTestId('wm-param-meta').fill('replaced note');
    await setOutName(page, '-');
    await submitWm(page);
    runNum++;
    const convMetaActual = await collectToolOutputs(page, runNum);
    const convMetaFrozen = readFrozenNcbi('v_meta_new');
    expect(sha256(convMetaActual.stdout), 'v_meta_new vs NCBI').toBe(sha256(convMetaFrozen.stdout));

    assertNoNet(net);
  });

  test('M5 -ustat across all 8 -outfmt formats, -dust, -ids, -infmt blastdb, and -out/-logfile sharing match native CLI and frozen NCBI bytes', async ({
    page,
  }, testInfo) => {
    await openWindowMasker(page);
    await selectWmStage(page, 'mask');
    await openWmOptions(page);

    const countsPath = join(WINMASK_DIR, 'counts/dm2L_u9.binary');
    const multiPath = join(WM_DIR, 'in/multi.fa');
    const idsPath = join(WM_DIR, 'in/ids_a.txt');
    const nativeDir = ensureDir(testInfo.outputPath('native-mask'));
    copyFileSync(countsPath, join(nativeDir, 'dm2L_u9.binary'));
    copyFileSync(multiPath, join(nativeDir, 'multi.fa'));
    copyFileSync(idsPath, join(nativeDir, 'ids_a.txt'));

    await addWmFastaFiles(page, [{ name: 'multi.fa', buffer: readFileSync(multiPath) }]);
    await pickCounts(page, 'ustat', 'dm2L_u9.binary', countsPath);

    let runNum = 0;

    // Every -outfmt of -ustat.
    const outfmts = [
      'interval',
      'fasta',
      'seqloc_asn1_bin',
      'seqloc_asn1_text',
      'seqloc_xml',
      'maskinfo_asn1_bin',
      'maskinfo_asn1_text',
      'maskinfo_xml',
    ] as const;

    for (const fmt of outfmts) {
      await showWmTab(page);
      await page.getByTestId('wm-param-outfmt').selectOption(fmt);
      await submitWm(page);
      runNum++;
      const outfmtActual = await collectToolOutputs(page, runNum);
      const native = runNativeTool(['-ustat', 'dm2L_u9.binary', '-in', 'multi.fa', '-outfmt', fmt], nativeDir);
      const frozen = readFrozenNcbi(`o_${fmt}`);
      expect(sha256(outfmtActual.stdout), `o_${fmt} vs native`).toBe(sha256(native.stdout));
      expect(sha256(outfmtActual.stdout), `o_${fmt} vs NCBI`).toBe(sha256(frozen.stdout));
      expect(sha256(outfmtActual.stderr), `o_${fmt} stderr`).toBe(sha256(frozen.stderr));
      if (fmt.endsWith('_bin')) {
        await openOutputTab(page, 'stdout');
        await expectBinaryNotice(page);
      }
    }

    // -dust T with maskinfo_asn1_text (d_multi_mi).
    await showWmTab(page);
    await page.getByTestId('wm-param-outfmt').selectOption('maskinfo_asn1_text');
    await page.getByTestId('wm-param-dust').selectOption('true');
    await submitWm(page);
    runNum++;
    const dustActual = await collectToolOutputs(page, runNum);
    const dustFrozen = readFrozenNcbi('d_multi_mi');
    expect(sha256(dustActual.stdout), 'd_multi_mi vs NCBI').toBe(sha256(dustFrozen.stdout));

    // -ids via file picker (i_ids_a).
    await showWmTab(page);
    await restoreDefaults(page);
    await page.getByTestId('wm-param-outfmt').selectOption('interval');
    await page.getByTestId('wm-param-ids-file').setInputFiles({
      name: 'ids_a.txt',
      mimeType: 'text/plain',
      buffer: readFileSync(idsPath),
    });
    await submitWm(page);
    runNum++;
    const idsActual = await collectToolOutputs(page, runNum);
    const idsFrozen = readFrozenNcbi('i_ids_a');
    expect(sha256(idsActual.stdout), 'i_ids_a vs NCBI').toBe(sha256(idsFrozen.stdout));

    // -out x.txt -logfile x.txt (single file opened twice, POSIX offset/truncate).
    await showWmTab(page);
    await restoreDefaults(page);
    await setOutName(page, 'x.txt');
    await setExtra(page, '-logfile x.txt');
    await submitWm(page);
    runNum++;
    const sharedActual = await collectToolOutputs(page, runNum);
    const sharedNative = runNativeTool(
      ['-ustat', 'dm2L_u9.binary', '-in', 'multi.fa', '-out', 'x.txt', '-logfile', 'x.txt'],
      nativeDir,
      ['x.txt'],
    );
    const sharedFile = sharedActual.files.get('x.txt');
    expect(sharedFile?.label).toBe('-out, -logfile x.txt');
    expect(sha256(sharedFile!.bytes), '-out/-logfile shared file').toBe(sha256(sharedNative.files.get('x.txt')!));
    expect(sha256(sharedActual.stderr), '-out/-logfile stderr').toBe(sha256(sharedNative.stderr));

    // -infmt blastdb with Open folder… (b_db).
    await showWmTab(page);
    await restoreDefaults(page);
    await clearWmRunFiles(page);
    await setOutName(page, '');
    const dbFolder = copyBlastDbFolder(testInfo.outputPath('blastdb-mask'), 'wm');
    copyBlastDbFolder(nativeDir, 'wm');
    await useDb(page, dbFolder, 'db/wm');
    await submitWm(page);
    runNum++;
    const bdbActual = await collectToolOutputs(page, runNum);
    const bdbFrozen = readFrozenNcbi('b_db');
    expect(sha256(bdbActual.stdout), 'b_db vs NCBI').toBe(sha256(bdbFrozen.stdout));
  });

  test('Use in BLASTN, Mask with these counts, .gz input, -h/-version-full, failure exit status, cancel, and shared queue', async ({
    page,
  }, testInfo) => {
    const net = trackNet(page);

    await openWindowMasker(page);
    await openWmOptions(page);

    const nativeDir = ensureDir(testInfo.outputPath('native-flow'));
    const s60kPath = join(BLASTN_REG_DIR, 's60k.fa');
    const queryPath = join(BLASTN_REG_DIR, 'wm_query.fa');
    const query = { name: 'wm_query.fa', buffer: readFileSync(queryPath) };
    copyFileSync(s60kPath, join(nativeDir, 's60k.fa'));
    copyFileSync(queryPath, join(nativeDir, 'wm_query.fa'));

    // 1) Make counts (.gz input!) -> byte-identical to native on uncompressed s60k.fa.
    const s60kGz = gzipSync(readFileSync(s60kPath));
    await addWmFastaFiles(page, [{ name: 's60k.fa.gz', buffer: s60kGz }]);
    await expect(page.getByTestId('wm-fasta-name-0')).toHaveText('s60k.fa');
    await expect(page.getByTestId('wm-fasta-inflated-0')).toHaveText('from s60k.fa.gz');
    await page.getByTestId('wm-param-unit').fill('9');
    await setSformat(page, 'obinary');
    await outIs(page, 'wmasker.obinary');
    await submitWm(page);
    const countsRun = await collectToolOutputs(page, 1);
    const nativeCounts = runNativeTool(
      ['-mk_counts', '-in', 's60k.fa', '-unit', '9', '-sformat', 'obinary', '-out', 'wmasker.obinary'],
      nativeDir,
      ['wmasker.obinary'],
    );
    expect(sha256(countsRun.files.get('wmasker.obinary')!.bytes), '.gz counts vs native').toBe(
      sha256(nativeCounts.files.get('wmasker.obinary')!),
    );

    // 2) "Use in BLASTN": outfmt 0, 6 and 7 equal native `losat blastn -window_masker_db wmasker.obinary`.
    await openOutputTab(page, 'file-0');
    await page.getByTestId('tool-use-in-blastn').click();
    await expect(page.getByTestId('tab-search')).toHaveAttribute('aria-selected', 'true');
    await expect(page.getByTestId('program-blastn')).toBeChecked();
    await expect(page.getByTestId('param-window_masker_db')).toBeChecked();
    await expect(page.getByTestId('param-window_masker_db-name')).toContainText('wmasker.obinary');

    await openFiles(page, 'query', [{ name: 'wm_query.fa', text: readFileSync(queryPath) }]);
    await openFiles(page, 'subject', [{ name: 's60k.fa', text: readFileSync(s60kPath) }]);
    await task(page, 'blastn');
    await submit(page);
    await waitStatus(page, 2, 'completed');
    const searchOut = await storedOutputs(page, 2);
    const nativeSearch = nativeExpectation(
      ['blastn', '-query', 'wm_query.fa', '-subject', 's60k.fa', '-task', 'blastn', '-window_masker_db', 'wmasker.obinary'],
      nativeDir,
    );
    expect({ 0: sha256(searchOut.outputs[0]), 6: sha256(searchOut.outputs[6]), 7: sha256(searchOut.outputs[7]) }).toEqual(
      nativeSearch.sha256,
    );
    expect(sha256(searchOut.warnings)).toBe(nativeSearch.stderrSha256);

    // 3) "Mask with these counts": preselects session counts on Mask stage and matches native.
    await pickRun(page, 1);
    await openOutputTab(page, 'file-0');
    await clickMaskWithCounts(page);
    await expectWmTabSelected(page);
    await expect(page.getByTestId('wm-stage-mask')).toBeChecked();
    await expect(page.getByTestId('wm-ustat-chosen')).toContainText('wmasker.obinary');
    await clearWmFasta(page);
    await addWmFastaFiles(page, [query]);
    await submitWm(page);
    const maskOut = await collectToolOutputs(page, 3);
    const nativeMask = runNativeTool(['-ustat', 'wmasker.obinary', '-in', 'wm_query.fa'], nativeDir);
    expect(sha256(maskOut.stdout), 'Mask with made counts vs native').toBe(sha256(nativeMask.stdout));

    // 4) -h and -version-full in Additional options (no inputs required).
    await showWmTab(page);
    await selectWmStage(page, 'mk_counts');
    await clearWmFasta(page);
    await restoreDefaults(page);
    await setExtra(page, '-h');
    await submitWm(page);
    const helpOut = await collectToolOutputs(page, 4);
    const nativeHelp = runNativeTool(['-mk_counts', '-out', 'wmasker.obinary', '-h'], nativeDir);
    expect(sha256(helpOut.stdout), '-h vs native').toBe(sha256(nativeHelp.stdout));

    await showWmTab(page);
    await setExtra(page, '-version-full');
    await submitWm(page);
    const verOut = await collectToolOutputs(page, 5);
    const nativeVer = runNativeTool(['-mk_counts', '-out', 'wmasker.obinary', '-version-full'], nativeDir);
    expect(sha256(verOut.stdout), '-version-full vs native').toBe(sha256(nativeVer.stdout));

    // 5) Failed tool run (bad counts file -> Exit status 3, Warnings matches native stderr).
    await showWmTab(page);
    await setExtra(page, '');
    await selectWmStage(page, 'mask');
    const garbagePath = join(WM_DIR, 'in/garbage.counts');
    copyFileSync(garbagePath, join(nativeDir, 'garbage.counts'));
    await pickCounts(page, 'ustat', 'garbage.counts', garbagePath);
    await addWmFastaFiles(page, [query]);
    await submitWm(page);
    const failOut = await collectToolOutputs(page, 6, 'failed');
    const nativeFail = runNativeTool(['-ustat', 'garbage.counts', '-in', 'wm_query.fa'], nativeDir);
    expect(failOut.exitStatus).toBe(nativeFail.status);
    expect(sha256(failOut.stderr), 'failed stderr vs native').toBe(sha256(nativeFail.stderr));

    // 6) Cancel during a long -mk_counts (generated multi-pass input).
    const bigPath = testInfo.outputPath('cancel-input.fa');
    writeSyntheticFasta(bigPath, 4 * MIB);
    await showWmTab(page);
    await selectWmStage(page, 'mk_counts');
    await clearWmFasta(page);
    await page.getByTestId('wm-param-unit').fill('14');
    await page.getByTestId('wm-param-mem').fill('1');
    await setFastaFilePath(page, bigPath);
    await submitWm(page);
    await expect(page.getByTestId('run-7-option')).toHaveAttribute('data-status', 'running', { timeout: 60_000 });
    await clickCancelRun(page);
    await waitStatus(page, 7, 'cancelled');

    // 7) Tool run queued behind a running BLASTP search, then completes once the search is cancelled.
    // Prepare the WindowMasker form first so switching and submitting takes <200ms while BLASTP is running.
    await showWmTab(page);
    await clearWmFasta(page);
    await restoreDefaults(page);
    await setSformat(page, 'ascii');
    await page.getByTestId('wm-fasta-input').fill('>tiny\nACGTACGTACGTACGTACGT\n');
    await waitValidationSettled(page);

    await page.getByTestId('tab-search').click();
    await program(page, 'blastp');
    await page.getByTestId('threads').selectOption('1');
    for (const role of ['query', 'subject'] as const) {
      const existing = page.locator(`[data-testid^="${role}-source-"][data-status]`);
      while ((await existing.count()) > 0) await page.getByTestId(`${role}-source-0-remove`).click();
      await openFiles(page, role, [{ name: 'NZ_CP006932.faa', text: fasta('NZ_CP006932.faa') }]);
    }
    await submit(page);
    await expect(page.getByTestId('run-8-option')).toHaveAttribute('data-status', 'running', { timeout: 60_000 });

    await showWmTab(page);
    await submitWm(page);
    await expect(page.getByTestId('run-9-option')).toHaveAttribute('data-status', 'queued');
    await expect(page.getByTestId('run-status')).toContainText('Waiting for run 8');

    // Cancel BLASTP run 8; queued WindowMasker run 9 proceeds and completes.
    await pickRun(page, 8);
    await clickCancelRun(page);
    await waitStatus(page, 8, 'cancelled');
    await waitStatus(page, 9, 'completed');
    await pickRun(page, 9);
    await expect(page.getByTestId('tool-run-view')).toHaveAttribute('data-run', '9');
    await chooseRun(page, 2);

    assertNoNet(net);
  });

  test('64 MiB FASTA through -mk_counts -sformat ascii and -ustat -outfmt fasta grows page JS heap by < 32 MiB', async ({
    page,
    browserName,
  }, testInfo) => {
    test.skip(browserName !== 'chromium', 'performance.memory.usedJSHeapSize is Chromium-only');

    const fastaPath = testInfo.outputPath('genome64.fa');
    const fastaBytes = writeSyntheticFasta(fastaPath, 64 * MIB);
    expect(fastaBytes).toBeGreaterThanOrEqual(64 * MIB);

    const readHeap = async (): Promise<number> => {
      const value = await page.evaluate(
        () => (performance as Performance & { memory?: { usedJSHeapSize: number } }).memory?.usedJSHeapSize ?? -1,
      );
      expect(value, 'performance.memory.usedJSHeapSize available').toBeGreaterThan(0);
      return value;
    };

    await openWindowMasker(page);
    await setFastaFilePath(page, fastaPath);
    await setSformat(page, 'ascii');

    const heapBeforeMk = await readHeap();
    await submitWm(page);
    await waitStatus(page, 1, 'completed', 480_000);
    await pickRun(page, 1);
    await expect(page.getByTestId('tool-output-text')).toBeVisible({ timeout: 60_000 });
    const heapAfterMk = await readHeap();
    const mkDelta = heapAfterMk - heapBeforeMk;

    // Mask the same 64 MiB FASTA with the counts made by run 1, writing -outfmt fasta to OPFS.
    await openOutputTab(page, 'file-0');
    await clickMaskWithCounts(page);
    await page.getByTestId('wm-param-outfmt').selectOption('fasta');
    await setOutName(page, 'masked64.fa');

    const heapBeforeMask = await readHeap();
    await submitWm(page);
    await waitStatus(page, 2, 'completed', 480_000);
    await pickRun(page, 2);
    await expect(page.getByTestId('tool-output-truncated')).toBeVisible({ timeout: 60_000 });
    const heapAfterMask = await readHeap();
    const maskDelta = heapAfterMask - heapBeforeMask;

    const toMiB = (bytes: number) => `${(bytes / MIB).toFixed(2)} MiB`;
    console.log(
      `[64 MiB heap] mk_counts: ${toMiB(heapBeforeMk)} -> ${toMiB(heapAfterMk)} (delta ${toMiB(mkDelta)}); ` +
        `mask: ${toMiB(heapBeforeMask)} -> ${toMiB(heapAfterMask)} (delta ${toMiB(maskDelta)})`,
    );

    expect(mkDelta, `mk_counts JS heap growth ${toMiB(mkDelta)}`).toBeLessThan(32 * MIB);
    expect(maskDelta, `mask JS heap growth ${toMiB(maskDelta)}`).toBeLessThan(32 * MIB);
  });
});
