// LOSAT Web WW4: database masking (-db_soft_mask / -db_hard_mask) and BLASTN WindowMasker + -db
// (-window_masker_taxid / -window_masker_db / -no_taxid_expansion) form controls, responsive
// narrow-width checks, and real-engine byte parity against native CLI `losat`.
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { copyFileSync, mkdirSync, mkdtempSync, readdirSync, readFileSync, rmSync, statSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { axeViolations } from './support/axe';
import { ENGINE, NO_ENGINE_REASON } from './support/engine';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { chooseRun, openParameters, openSearch, paste, program, result, submit, task, waitStatus } from './support/search';
import { test } from './support/test';

const LOSAT = join(REPOSITORY, 'LOSAT');
const WM_COUNTS_PATH = join(LOSAT, 'tests/winmask/counts/mg1655_u9.obinary');
const WM_QUERY_PATH = join(LOSAT, 'tests/fixtures/blastn_regression/inputs/wm_query.fa');
const WM_TAXID = 562;
const WM_URL = `https://ftp.ncbi.nlm.nih.gov/blast/windowmasker_files/${WM_TAXID}/wmasker.obinary`;
const LAST_MODIFIED = 'Thu, 01 Oct 2026 08:00:00 GMT';

const chooseDatabaseSource = (page: Page) => page.getByTestId('subject-source-database').check();
const chooseFastaSource = (page: Page) => page.getByTestId('subject-source-fasta').check();
const expectFilterCount = (page: Page, count: number) =>
  expect(page.getByTestId('subject-db-filters-count')).toHaveText(`(${count} set)`);
const checkNoTaxidExpansion = (page: Page) => page.getByTestId('db-filter-no_taxid_expansion').check();
const selectDbMaskMode = (page: Page, mode: 'none' | 'soft' | 'hard') => page.getByTestId('db-mask-mode').selectOption(mode);
const expectValidationState = (page: Page, state: 'ok' | 'invalid') =>
  expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', state);
async function withTempWorkDir(prefix: string, fn: (work: string) => Promise<void>): Promise<void> {
  const work = mkdtempSync(join(tmpdir(), prefix));
  try {
    await fn(work);
  } finally {
    rmSync(work, { recursive: true, force: true });
  }
}
async function submitReadyRun(page: Page, runNumber: number): Promise<void> {
  await expectValidationState(page, 'ok');
  await submit(page);
  await waitStatus(page, runNumber, 'completed');
}
async function enableWmTaxid(page: Page, id: number): Promise<void> {
  await page.getByTestId('param-window_masker_db').check();
  await page.getByTestId('param-window_masker_db-taxid').fill(String(id));
}
const expectWmConflict = (page: Page, visible: boolean) =>
  visible
    ? expect(page.getByTestId('param-window_masker_db-no-taxid-conflict')).toBeVisible()
    : expect(page.getByTestId('param-window_masker_db-no-taxid-conflict')).toHaveCount(0);

async function openDbFilters(page: Page): Promise<void> {
  const details = page.getByTestId('subject-db-filters');
  if ((await details.getAttribute('open')) === null) {
    await page.getByTestId('subject-db-filters-toggle').click();
  }
}

const setSubjectDbFiles = (page: Page, files: Parameters<ReturnType<Page['getByTestId']>['setInputFiles']>[0]) =>
  page.getByTestId('subject-db-files-input').setInputFiles(files);
const sha256 = (bytes: Uint8Array | string): string => createHash('sha256').update(bytes).digest('hex');

function copyMatching(srcDir: string, dstDir: string, predicate: (name: string) => boolean): string[] {
  mkdirSync(dstDir, { recursive: true });
  const copied: string[] = [];
  for (const name of readdirSync(srcDir)) {
    const full = join(srcDir, name);
    if (!statSync(full).isFile() || !predicate(name)) continue;
    const target = join(dstDir, name);
    copyFileSync(full, target);
    copied.push(target);
  }
  return copied;
}

function runNative(
  args: readonly string[],
  cwd: string,
  format: 0 | 6 | 7,
  env: Readonly<Record<string, string>> = {},
): string {
  const full = [...args, '-outfmt', String(format), '-num_threads', '1'];
  const proc = spawnSync(NATIVE!, full, { cwd, maxBuffer: 1 << 28, env: { ...process.env, ...env } });
  if (proc.status !== 0) {
    throw new Error(`native losat ${full.join(' ')} exited ${proc.status}: ${proc.stderr.toString()}`);
  }
  return proc.stdout.toString('utf8');
}

/** Verifies byte-identical outfmt 0/6/7 output vs native `losat` and returns the outfmt 6 bytes. */
async function expectRunParity(
  page: Page,
  runNumber: number,
  formats: readonly (0 | 6 | 7)[],
  opts: {
    readonly nativeArgs: readonly string[];
    readonly cwd: string;
    readonly env?: Readonly<Record<string, string>>;
    readonly mustContain?: readonly string[];
    readonly mustOmit?: readonly string[];
    readonly label: string;
  },
): Promise<string> {
  let out6 = '';
  for (const fmt of formats) {
    const res = await result(page, runNumber, fmt);
    for (const token of opts.mustContain ?? []) expect(res.command).toContain(token);
    for (const token of opts.mustOmit ?? []) expect(res.command).not.toContain(token);
    const expected = runNative(opts.nativeArgs, opts.cwd, fmt, opts.env ?? {});
    expect(res.output, `${opts.label} outfmt ${fmt}`).toBe(expected);
    if (fmt === 6) out6 = res.output;
  }
  return out6;
}

/** Verifies zero horizontal overflow and that no visible "?" HelpTip starts a line alone. */
async function expectNoOrphanHelpTips(page: Page, context: string): Promise<void> {
  const measured = await page.evaluate(() => {
    const isVisibleElement = (el: Element): boolean => {
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
      while (cur !== null && !isVisibleElement(cur)) cur = cur.previousElementSibling;
      if (cur !== null) return cur;
      const parent = node.parentElement;
      if (parent !== null && parent.classList.contains('control-tip')) return findPrev(parent);
      return null;
    };
    const tips: Array<{ id: string; hasPrev: boolean; overlapsY: boolean; followsX: boolean }> = [];
    for (const btn of document.querySelectorAll('.help-toggle')) {
      const br = btn.getBoundingClientRect();
      if (br.width < 0.5 || br.height < 0.5) continue;
      const tip = btn.closest('.help-tip') ?? btn;
      const wrapper = tip.closest('.hint-line, .param-hint-line');
      let prevRect: DOMRect | undefined;
      if (wrapper !== null && wrapper.firstChild !== null && wrapper.firstChild !== tip) {
        prevRect = lastRangeRect(wrapper.firstChild, tip);
      } else {
        const prev = findPrev(tip);
        if (prev !== null) {
          prevRect = prev.matches('.param-hint, .hint, .file-missing, .file-name')
            ? lastRangeRect(prev)
            : prev.getBoundingClientRect();
        }
      }
      const id = btn.getAttribute('data-testid') ?? btn.getAttribute('aria-label') ?? 'help-tip';
      if (prevRect === undefined) {
        tips.push({ id, hasPrev: false, overlapsY: false, followsX: false });
      } else {
        tips.push({
          id,
          hasPrev: true,
          overlapsY: Math.min(br.bottom, prevRect.bottom) > Math.max(br.top, prevRect.top),
          followsX: br.left >= prevRect.right - 1,
        });
      }
    }
    return {
      scrollWidth: document.documentElement.scrollWidth,
      clientWidth: document.documentElement.clientWidth,
      tips,
    };
  });

  expect(measured.scrollWidth, `${context}: horizontal overflow`).toBeLessThanOrEqual(measured.clientWidth);
  expect(measured.tips.length, `${context}: expected visible HelpTips`).toBeGreaterThan(0);
  for (const tip of measured.tips) {
    expect(tip.hasPrev, `${context}: ${tip.id} has preceding content`).toBe(true);
    expect(tip.overlapsY, `${context}: ${tip.id} shares a line vertically`).toBe(true);
    expect(tip.followsX, `${context}: ${tip.id} sits after preceding content`).toBe(true);
  }
}

test.describe('WW4 database masking & WindowMasker -db form controls', () => {
  test('enforces exclusive -db_soft_mask/-db_hard_mask controls, -no_taxid_expansion interactions, responsive layout, and axe a11y', async ({
    page,
  }) => {
    await openSearch(page);
    await program(page, 'blastn');
    await paste(page, 'query', '>q\nACGTACGTACGTACGTACGTACGTACGT\n');
    await chooseDatabaseSource(page);
    await setSubjectDbFiles(page, [
      { name: 'a5.nin', mimeType: 'application/octet-stream', buffer: Buffer.from('0') },
      { name: 'a5.nhr', mimeType: 'application/octet-stream', buffer: Buffer.from('1') },
      { name: 'a5.nsq', mimeType: 'application/octet-stream', buffer: Buffer.from('2') },
      { name: 'ps5.pin', mimeType: 'application/octet-stream', buffer: Buffer.from('p') },
    ]);
    await openDbFilters(page);

    // 1. Database masking select: None -> Soft -> Hard -> None, preserving algorithm across Soft/Hard.
    const maskMode = page.getByTestId('db-mask-mode');
    await expect(maskMode).toHaveValue('none');
    await expect(page.getByTestId('db-mask-algorithm')).toHaveCount(0);

    await selectDbMaskMode(page, 'soft');
    const maskAlg = page.getByTestId('db-mask-algorithm');
    await expect(maskAlg).toBeVisible();
    await expect(maskAlg).toHaveAttribute('placeholder', '11 or dust');
    await expect(page.getByTestId('db-mask-missing')).toHaveText('Enter a database soft mask algorithm ID or name.');
    await submit(page);
    await expect(page.getByTestId('search-message')).toHaveText('Enter a database soft mask algorithm ID or name.');

    await maskAlg.fill('30');
    await expect(page.getByTestId('db-mask-missing')).toHaveCount(0);
    await expectFilterCount(page, 1);

    // Switching Soft -> Hard keeps the algorithm value and enforces mutual exclusion.
    await selectDbMaskMode(page, 'hard');
    await expect(maskAlg).toHaveValue('30');
    await expectFilterCount(page, 1);

    // 2. -no_taxid_expansion checked -> enabling WindowMasker defaults to Counts file (-window_masker_db).
    await checkNoTaxidExpansion(page);
    await expectFilterCount(page, 2);

    await openParameters(page);
    const wmCheck = page.getByTestId('param-window_masker_db');
    await wmCheck.check();
    const wmBy = page.getByTestId('param-window_masker_db-by');
    await expect(wmBy).toHaveValue('file');
    await expect(wmBy.locator('option[value="taxid"]')).toBeDisabled();
    await expect(page.getByTestId('param-window_masker_db-taxid-disabled-note')).toContainText('-no_taxid_expansion');

    // 3. Responsive widths (320, 390, 768, 1280): zero overflow, zero orphaned "?" HelpTips, zero axe issues.
    for (const width of [320, 390, 768, 1280]) {
      await page.setViewportSize({ width, height: 844 });
      await expectNoOrphanHelpTips(page, `WW4 blastn -db at ${width}px`);
      if (width === 390 || width === 1280) {
        const violations = await axeViolations(page, ['main']);
        expect(violations, `axe violations at ${width}px`).toEqual([]);
      }
    }
    await page.setViewportSize({ width: 1280, height: 900 });

    // 4. Conflict state via FASTA -> Database switch when both -no_taxid_expansion and -window_masker_taxid are set.
    await chooseFastaSource(page);
    await wmBy.selectOption('taxid');
    await page.getByTestId('param-window_masker_db-taxid').fill('562');
    await chooseDatabaseSource(page);

    // Conflict banner appears in WindowMasker and Database masking & filters shows the conflict note.
    await expectWmConflict(page, true);
    await expect(page.getByTestId('db-filter-no_taxid_expansion-conflict')).toContainText('-window_masker_taxid');

    // One-click "Clear -no_taxid_expansion" resolves the conflict cleanly.
    await page.getByTestId('param-window_masker_db-clear-noexp').click();
    await expectWmConflict(page, false);
    await expect(page.getByTestId('db-filter-no_taxid_expansion')).not.toBeChecked();
    await expect(page.getByTestId('db-filter-no_taxid_expansion')).toBeDisabled();

    // 5. Switching to BLASTP updates the database mask placeholder to protein ("21 or seg").
    await program(page, 'blastp');
    await openDbFilters(page);
    await selectDbMaskMode(page, 'soft');
    await expect(page.getByTestId('db-mask-algorithm')).toHaveAttribute('placeholder', '21 or seg');
  });
});

test.describe('WW4 real engine byte parity vs native losat', () => {
  test.skip(ENGINE === undefined, NO_ENGINE_REASON);
  test.skip(NATIVE === undefined, NO_NATIVE_REASON);
  test.setTimeout(300_000);

  test('matches native losat bytes for -db_soft_mask and -db_hard_mask on nucleotide (blastn) and protein (blastp) databases', async ({
    page,
  }) => {
    await withTempWorkDir('losat-ww4-dbmask-', async (work) => {
      // 1. Nucleotide DB (`a5`) with `-task blastn`: soft mask 30, hard mask 30, and soft mask by name `windowmasker`.
      const nucDir = join(work, 'nucdb');
      const nucFiles = copyMatching(join(LOSAT, 'tests/db_search/dbmask'), nucDir, (name) => name.startsWith('a5.'));
      const nucQueryPath = join(LOSAT, 'tests/db_search/dbmask/q.fna');
      const nucQueryText = readFileSync(nucQueryPath, 'utf8');
      const nucBase = ['blastn', '-query', nucQueryPath, '-db', 'a5', '-task', 'blastn'];

      await openSearch(page);
      await program(page, 'blastn');
      await paste(page, 'query', nucQueryText);
      await chooseDatabaseSource(page);
      await setSubjectDbFiles(page, nucFiles);
      await expect(page.getByTestId('subject-db-selected')).toContainText('a5');
      await openParameters(page);
      await task(page, 'blastn');

      // Run 1: -db_soft_mask 30
      await openDbFilters(page);
      await selectDbMaskMode(page, 'soft');
      await page.getByTestId('db-mask-algorithm').fill('30');
      await submit(page);
      await waitStatus(page, 1, 'completed');

      const unmaskedNuc6 = runNative(nucBase, nucDir, 6);
      const softNuc6 = await expectRunParity(page, 1, [0, 6, 7], {
        nativeArgs: [...nucBase, '-db_soft_mask', '30'],
        cwd: nucDir,
        mustContain: ['-db_soft_mask 30'],
        mustOmit: ['-db_hard_mask'],
        label: 'blastn -db_soft_mask 30',
      });
      expect(softNuc6).not.toBe(unmaskedNuc6);

      // Run 2: switch to -db_hard_mask 30 (value preserved)
      await selectDbMaskMode(page, 'hard');
      await expect(page.getByTestId('db-mask-algorithm')).toHaveValue('30');
      await submit(page);
      await waitStatus(page, 2, 'completed');

      const hardNuc6 = await expectRunParity(page, 2, [0, 6, 7], {
        nativeArgs: [...nucBase, '-db_hard_mask', '30'],
        cwd: nucDir,
        mustContain: ['-db_hard_mask 30'],
        mustOmit: ['-db_soft_mask'],
        label: 'blastn -db_hard_mask 30',
      });
      expect(hardNuc6).not.toBe(softNuc6);
      expect(hardNuc6).not.toBe(unmaskedNuc6);

      // Run 3: -db_soft_mask by algorithm name (`windowmasker`)
      await selectDbMaskMode(page, 'soft');
      await page.getByTestId('db-mask-algorithm').fill('windowmasker');
      await submit(page);
      await waitStatus(page, 3, 'completed');
      await expectRunParity(page, 3, [6], {
        nativeArgs: [...nucBase, '-db_soft_mask', 'windowmasker'],
        cwd: nucDir,
        mustContain: ['-db_soft_mask windowmasker'],
        label: 'blastn -db_soft_mask windowmasker',
      });

      // 2. Protein DB (`ps5`): -db_soft_mask 21 and -db_hard_mask 21 (outfmt 0, 6, 7)
      const protDir = join(work, 'protdb');
      const protFiles = copyMatching(join(LOSAT, 'tests/db_search/dbmaskp'), protDir, (name) => name.startsWith('ps5.'));
      const protQueryPath = join(LOSAT, 'tests/db_search/dbmaskp/qp.faa');
      const protQueryText = readFileSync(protQueryPath, 'utf8');
      const protBase = ['blastp', '-query', protQueryPath, '-db', 'ps5'];

      await program(page, 'blastp');
      await paste(page, 'query', protQueryText);
      await setSubjectDbFiles(page, protFiles);
      await expect(page.getByTestId('subject-db-selected')).toContainText('ps5');
      await openDbFilters(page);
      await selectDbMaskMode(page, 'soft');
      await page.getByTestId('db-mask-algorithm').fill('21');
      await submit(page);
      await waitStatus(page, 4, 'completed');

      const unmaskedProt6 = runNative(protBase, protDir, 6);
      const softProt6 = await expectRunParity(page, 4, [0, 6, 7], {
        nativeArgs: [...protBase, '-db_soft_mask', '21'],
        cwd: protDir,
        mustContain: ['-db_soft_mask 21'],
        mustOmit: ['-db_hard_mask'],
        label: 'blastp -db_soft_mask 21',
      });
      expect(softProt6).not.toBe(unmaskedProt6);

      await selectDbMaskMode(page, 'hard');
      await expect(page.getByTestId('db-mask-algorithm')).toHaveValue('21');
      await submit(page);
      await waitStatus(page, 5, 'completed');

      const hardProt6 = await expectRunParity(page, 5, [0, 6, 7], {
        nativeArgs: [...protBase, '-db_hard_mask', '21'],
        cwd: protDir,
        mustContain: ['-db_hard_mask 21'],
        mustOmit: ['-db_soft_mask'],
        label: 'blastp -db_hard_mask 21',
      });
      expect(hardProt6).not.toBe(softProt6);
      expect(hardProt6).not.toBe(unmaskedProt6);
    });
  });

  test('matches native losat bytes for BLASTN -db with -window_masker_taxid, -window_masker_db + -no_taxid_expansion, and real refusal on conflict', async ({
    page,
  }) => {
    await withTempWorkDir('losat-ww4-wmdb-', async (work) => {
      const countsBytes = readFileSync(WM_COUNTS_PATH);
      const wmPath = join(work, 'wm_path');
      mkdirSync(join(wmPath, String(WM_TAXID)), { recursive: true });
      copyFileSync(WM_COUNTS_PATH, join(wmPath, String(WM_TAXID), 'wmasker.obinary'));

      const dbDir = join(work, 'wmdb');
      const dbFiles = copyMatching(join(LOSAT, 'tests/db_search'), dbDir, (name) => name.startsWith('wm_s60k.'));
      const queryText = readFileSync(WM_QUERY_PATH, 'utf8');
      const wmBase = ['blastn', '-query', WM_QUERY_PATH, '-db', 'wm_s60k', '-task', 'blastn'];

      // Stub NCBI FTP for taxid 562 (zero live network requests).
      await page.route(/^https:\/\/ftp\.ncbi\.nlm\.nih\.gov\//, async (route) => {
        const cors = { 'Access-Control-Allow-Origin': '*', 'Access-Control-Expose-Headers': 'Content-Range, Content-Length, Last-Modified' };
        if (route.request().url() !== WM_URL) {
          await route.fulfill({ status: 404, headers: cors, body: 'Not Found' });
          return;
        }
        await route.fulfill({
          status: 200,
          headers: {
            ...cors,
            'Content-Type': 'application/octet-stream',
            'Content-Length': String(countsBytes.length),
            'Last-Modified': LAST_MODIFIED,
          },
          body: countsBytes,
        });
      });

      await openSearch(page);
      await program(page, 'blastn');
      await paste(page, 'query', queryText);
      await chooseDatabaseSource(page);
      await setSubjectDbFiles(page, dbFiles);
      await expect(page.getByTestId('subject-db-selected')).toContainText('wm_s60k');
      await openParameters(page);
      await task(page, 'blastn');

      // 1. BLASTN -db + -window_masker_taxid 562 (OPFS download + search)
      await enableWmTaxid(page, WM_TAXID);
      await page.getByTestId('counts-download').click();
      await expect(page.getByTestId('counts-delete')).toBeVisible();
      await submitReadyRun(page, 1);
      const unmasked6 = runNative(wmBase, dbDir, 6);
      const taxid6 = await expectRunParity(page, 1, [0, 6, 7], {
        nativeArgs: [...wmBase, '-window_masker_taxid', String(WM_TAXID)],
        cwd: dbDir,
        env: { WINDOW_MASKER_PATH: wmPath },
        mustContain: ['-db wm_s60k', `-window_masker_taxid ${WM_TAXID}`],
        label: 'blastn -db -window_masker_taxid 562',
      });
      expect(taxid6).not.toBe(unmasked6);

      // Run details records the OPFS/downloaded WindowMasker provenance even for -db runs.
      await chooseRun(page, 1);
      await page.getByTestId('results-view-details').click();
      await expect(page.getByTestId('run-window-masker')).toContainText(`SHA-256 ${sha256(countsBytes)}`);
      await expect(page.getByTestId('run-window-masker-url')).toHaveText(WM_URL);
      await page.getByTestId('tab-search').click();

      // 2. Real engine refusal when -window_masker_taxid and -no_taxid_expansion are both active.
      await page.getByTestId('param-window_masker_db').uncheck();
      await openDbFilters(page);
      await checkNoTaxidExpansion(page);
      await chooseFastaSource(page);
      await enableWmTaxid(page, WM_TAXID);
      await chooseDatabaseSource(page);

      await expectValidationState(page, 'invalid');
      await expect(page.getByTestId('argv-validation')).toContainText("Incompatible with argument:  `no_taxid_expansion'");
      await expect(page.getByTestId('param-window_masker_db-taxid')).toHaveAttribute('aria-invalid', 'true');
      await expect(page.getByTestId('db-filter-no_taxid_expansion')).toHaveAttribute('aria-invalid', 'true');
      await expectWmConflict(page, true);

      // 3. One-click "Use Counts file" switches to -window_masker_db while keeping -no_taxid_expansion -> valid and byte-identical!
      await page.getByTestId('param-window_masker_db-use-file').click();
      await expect(page.getByTestId('param-window_masker_db-by')).toHaveValue('file');
      await page.getByTestId('param-window_masker_db-file').setInputFiles(WM_COUNTS_PATH);
      await submitReadyRun(page, 2);
      await expectRunParity(page, 2, [0, 6, 7], {
        nativeArgs: [...wmBase, '-window_masker_db', WM_COUNTS_PATH, '-no_taxid_expansion'],
        cwd: dbDir,
        mustContain: ['-db wm_s60k', '-window_masker_db mg1655_u9.obinary', '-no_taxid_expansion'],
        mustOmit: ['-window_masker_taxid'],
        label: 'blastn -db -window_masker_db -no_taxid_expansion',
      });
    });
  });
});
