// LOSAT Web local BLAST database (-db) & BLASTN repeats database (-filtering_db) UI and
// byte-parity E2E tests (WQA D12b + RG-D4 W2).
import { spawnSync } from 'node:child_process';
import { copyFileSync, mkdirSync, mkdtempSync, readdirSync, readFileSync, rmSync, statSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { axeViolations } from './support/axe';
import { ENGINE, NO_ENGINE_REASON } from './support/engine';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { openParameters, openSearch, paste, program, result, submit, waitStatus } from './support/search';
import { test } from './support/test';

const LOSAT = join(REPOSITORY, 'LOSAT');
const PRIVACY_NOTICE = 'LOSAT Web reads the database from this device; nothing is uploaded.';

const chooseDatabaseSource = (page: Page) => page.getByTestId('subject-source-database').check();
const toggleDbFilters = (page: Page) => page.getByTestId('subject-db-filters-toggle').click();
const enableFilteringDb = (page: Page) => page.getByTestId('param-filtering_db').check();
const setSubjectDbFiles = (page: Page, files: Parameters<ReturnType<Page['getByTestId']>['setInputFiles']>[0]) =>
  page.getByTestId('subject-db-files-input').setInputFiles(files);

async function expectFilteringDbParity(
  page: Page,
  runIndex: number,
  subjectArgs: readonly string[],
  queryPath: string,
  cwd: string,
  label: string,
): Promise<void> {
  const { output, command } = await result(page, runIndex, 6);
  expect(command).toContain('-filtering_db rpt_d4');
  if (subjectArgs[0] === '-db') expect(command).toContain(`-db ${subjectArgs[1]}`);
  const expected = runNative(
    ['blastn', '-query', queryPath, ...subjectArgs, '-filtering_db', 'rpt_d4', '-soft_masking', 'false'],
    cwd,
    6,
  );
  expect(output, label).toBe(expected);
}

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
  format: 0 | 6 | 7 | string,
  env: Readonly<Record<string, string>> = {},
): string {
  const full = [...args, '-outfmt', String(format), '-num_threads', '1'];
  const proc = spawnSync(NATIVE!, full, { cwd, maxBuffer: 1 << 28, env: { ...process.env, ...env } });
  if (proc.status !== 0) {
    throw new Error(`native losat ${full.join(' ')} exited ${proc.status}: ${proc.stderr.toString()}`);
  }
  return proc.stdout.toString('utf8');
}

test.describe('local BLAST database (-db) and repeats database (-filtering_db) UI', () => {
  test('switches Subject to Database, discovers multi-volume & molecule types, shows privacy and reload notices, and passes axe at 1280px and 390px', async ({
    page,
  }) => {
    await openSearch(page);

    // Switch Subject to Database (-db)
    await chooseDatabaseSource(page);
    await expect(page.getByTestId('subject-db-privacy')).toHaveText(PRIVACY_NOTICE);

    // Pick a mixed multi-volume nucleotide + protein set via multi-file fallback
    await setSubjectDbFiles(page, [
      { name: 'multi.nal', mimeType: 'text/plain', buffer: Buffer.from('TITLE Multi Nuc\nDBLIST multi.00 multi.01\n') },
      { name: 'multi.00.nin', mimeType: 'application/octet-stream', buffer: Buffer.from('0') },
      { name: 'multi.01.nin', mimeType: 'application/octet-stream', buffer: Buffer.from('1') },
      { name: 'extra.nin', mimeType: 'application/octet-stream', buffer: Buffer.from('2') },
      { name: 'prot.pin', mimeType: 'application/octet-stream', buffer: Buffer.from('p') },
      { name: 'taxdb.bti', mimeType: 'application/octet-stream', buffer: Buffer.from('i') },
      { name: 'taxdb.btd', mimeType: 'application/octet-stream', buffer: Buffer.from('d') },
      { name: 'taxonomy4blast.sqlite3', mimeType: 'application/octet-stream', buffer: Buffer.from('s') },
    ]);

    // Multi-volume alias hides .00/.01 components; both nucleotide DBs appear in select
    const select = page.getByTestId('subject-db-select');
    await expect(select).toBeVisible();
    const options = await select.locator('option').evaluateAll((els) => els.map((el) => (el as HTMLOptionElement).value));
    expect(options).toEqual(['extra', 'multi']);
    await expect(page.getByTestId('subject-db-meta')).toContainText('taxonomy: taxdb, taxonomy4blast.sqlite3');

    await select.selectOption('multi');
    await expect(page.getByTestId('search-summary-line')).toContainText('database multi');

    // Open Database filters and attach -seqidlist by name
    await toggleDbFilters(page);
    await page.getByTestId('db-filter-seqidlist-input').setInputFiles({
      name: 'subset.ids',
      mimeType: 'text/plain',
      buffer: Buffer.from('seq1\n'),
    });
    await expect(page.getByTestId('db-filter-seqidlist-name')).toContainText('subset.ids');
    await expect(page.getByTestId('subject-db-filters-count')).toHaveText('(1 set)');

    // Enable BLASTN Repeats database (-filtering_db) using the shared DatabasePicker
    await openParameters(page);
    await enableFilteringDb(page);
    await expect(page.getByTestId('param-filtering_db-privacy')).toHaveText(PRIVACY_NOTICE);
    await page.getByTestId('param-filtering_db-files-input').setInputFiles([
      { name: 'repeats.nin', mimeType: 'application/octet-stream', buffer: Buffer.from('r') },
    ]);
    await expect(page.getByTestId('param-filtering_db-selected')).toContainText('repeats');

    // Desktop 1280px & Mobile 390px: zero horizontal overflow and zero axe accessibility violations
    for (const width of [1280, 390]) {
      await page.setViewportSize({ width, height: 844 });
      const overflow = await page.evaluate(() => document.documentElement.scrollWidth > window.innerWidth);
      expect(overflow, `horizontal overflow at ${width}px`).toBe(false);
      const violations = await axeViolations(page, ['main']);
      expect(violations, `axe violations at ${width}px`).toEqual([]);
    }

    // Switching to BLASTP filters to the protein DB (`prot`) automatically
    await page.setViewportSize({ width: 1280, height: 900 });
    await program(page, 'blastp');
    await expect(page.getByTestId('subject-db-selected')).toContainText('prot');

    // Picking only a nucleotide DB in BLASTP shows a clear molecule-type mismatch error
    await setSubjectDbFiles(page, [{ name: 'only_nuc.nin', mimeType: 'application/octet-stream', buffer: Buffer.from('n') }]);
    await expect(page.getByTestId('subject-db-mismatch')).toContainText(
      'Found nucleotide database (only_nuc), not a protein database for BLASTP.',
    );

    // Reloading the page after referencing a DB shows the reload prompt ("Choose the folder again.")
    await page.reload({ waitUntil: 'domcontentloaded' });
    await expect(page.getByTestId('search-panel')).toHaveAttribute('data-describe', /^(ready|failed)$/, { timeout: 30_000 });
    await expect(page.getByTestId('subject-source-database')).toBeChecked();
    await expect(page.getByTestId('subject-db-reload')).toContainText('Choose the folder again.');
  });
});

test.describe('real engine byte-parity for -db, filter lists, taxonomy, and BLASTN -filtering_db', () => {
  test.skip(ENGINE === undefined, NO_ENGINE_REASON);
  test.skip(NATIVE === undefined, NO_NATIVE_REASON);
  test.setTimeout(300_000);

  test('matches native losat bytes for nucleotide -db, protein -db, -seqidlist, -negative_seqidlist, taxonomy/taxdb, and BLASTN -filtering_db', async ({
    page,
    browserName,
  }) => {
    const work = mkdtempSync(join(tmpdir(), 'losat-d12b-e2e-'));
    try {
      // 1. Nucleotide -db (`blastn` against `nslice`): folder pick (`webkitdirectory`) on Chromium gives `-db nucdb/nslice`, multi-file fallback on others gives `-db nslice`.
      const nucDir = join(work, 'nucdb');
      const nucFiles = copyMatching(join(LOSAT, 'tests/db_search'), nucDir, (name) => name.startsWith('nslice.'));
      const nucQueryText = readFileSync(join(LOSAT, 'tests/db_search/nslice_q.fna'), 'utf8');
      const nucQueryPath = join(work, 'nslice_q.fna');
      writeFileSync(nucQueryPath, nucQueryText);

      await openSearch(page);
      await program(page, 'blastn');
      await paste(page, 'query', nucQueryText);
      await chooseDatabaseSource(page);
      const useFolderPick = browserName === 'chromium';
      const nucDbArg = useFolderPick ? 'nucdb/nslice' : 'nslice';
      const nucCwd = useFolderPick ? work : nucDir;
      if (useFolderPick) {
        await page.getByTestId('subject-db-folder-input').setInputFiles(nucDir);
      } else {
        await setSubjectDbFiles(page, nucFiles);
      }
      await expect(page.getByTestId('subject-db-selected')).toContainText(nucDbArg);
      await submit(page);
      await waitStatus(page, 1, 'completed');

      for (const fmt of [0, 6, 7] as const) {
        const { output, command } = await result(page, 1, fmt);
        expect(command).toContain(`-db ${nucDbArg}`);
        const expected = runNative(['blastn', '-query', nucQueryPath, '-db', nucDbArg], nucCwd, fmt);
        expect(output, `blastn -db outfmt ${fmt}`).toBe(expected);
      }

      // 2. Protein -db + -seqidlist + -negative_seqidlist + taxonomy4blast.sqlite3 + taxdb (`blastp` against `taxfilter/tfp`)
      const protDir = join(work, 'protdb');
      const protFiles = copyMatching(join(LOSAT, 'tests/db_search/taxfilter'), protDir, (name) =>
        name.startsWith('tfp.') || name === 'taxonomy4blast.sqlite3',
      );
      copyFileSync(join(LOSAT, 'tests/makembindex/in/taxt5/taxdb.bti'), join(protDir, 'taxdb.bti'));
      copyFileSync(join(LOSAT, 'tests/makembindex/in/taxt5/taxdb.btd'), join(protDir, 'taxdb.btd'));
      protFiles.push(join(protDir, 'taxdb.bti'), join(protDir, 'taxdb.btd'));

      const protQueryPath = join(LOSAT, 'tests/db_search/taxfilter/q.faa');
      const protQueryText = readFileSync(protQueryPath, 'utf8');

      await program(page, 'blastp');
      await paste(page, 'query', protQueryText);
      await setSubjectDbFiles(page, protFiles);
      await expect(page.getByTestId('subject-db-selected')).toContainText('tfp');
      await expect(page.getByTestId('subject-db-meta')).toContainText('taxdb, taxonomy4blast.sqlite3');

      // Plain protein -db search + verify non-zero subject length in SubjectTable
      await submit(page);
      await waitStatus(page, 2, 'completed');
      const firstLen = await page.locator('[data-testid="subject-table"] [data-field="length"]').first().textContent();
      expect(Number((firstLen ?? '0').replace(/,/g, ''))).toBeGreaterThan(0);
      for (const fmt of [0, 6, 7] as const) {
        const { output } = await result(page, 2, fmt);
        const expected = runNative(['blastp', '-query', protQueryPath, '-db', 'tfp'], protDir, fmt, { BLASTDB: protDir });
        expect(output, `blastp -db outfmt ${fmt}`).toBe(expected);
      }

      // -seqidlist mounted by name
      const slCopy = join(protDir, 'sl_tfp.txt');
      copyFileSync(join(LOSAT, 'tests/db_search/seqidlist/sl_tfp.txt'), slCopy);
      await toggleDbFilters(page);
      await page.getByTestId('db-filter-seqidlist-input').setInputFiles(slCopy);
      await submit(page);
      await waitStatus(page, 3, 'completed');
      {
        const { output, command } = await result(page, 3, 6);
        expect(command).toContain('-seqidlist sl_tfp.txt');
        const expected = runNative(['blastp', '-query', protQueryPath, '-db', 'tfp', '-seqidlist', 'sl_tfp.txt'], protDir, 6, {
          BLASTDB: protDir,
        });
        expect(output, 'blastp -db -seqidlist outfmt 6').toBe(expected);
      }

      // -negative_seqidlist mounted by name
      const negCopy = join(protDir, 'sl_tfp_neg.txt');
      copyFileSync(join(LOSAT, 'tests/db_search/seqidlist/sl_tfp_neg.txt'), negCopy);
      await page.getByTestId('db-filter-seqidlist-remove').click();
      await page.getByTestId('db-filter-negative_seqidlist-input').setInputFiles(negCopy);
      await submit(page);
      await waitStatus(page, 4, 'completed');
      {
        const { output, command } = await result(page, 4, 6);
        expect(command).toContain('-negative_seqidlist sl_tfp_neg.txt');
        const expected = runNative(
          ['blastp', '-query', protQueryPath, '-db', 'tfp', '-negative_seqidlist', 'sl_tfp_neg.txt'],
          protDir,
          6,
          { BLASTDB: protDir },
        );
        expect(output, 'blastp -db -negative_seqidlist outfmt 6').toBe(expected);
      }

      // Taxonomy filter (-taxids 40674 Mammalia) expanded with taxonomy4blast.sqlite3 + taxdb in folder
      await page.getByTestId('db-filter-negative_seqidlist-remove').click();
      await page.getByTestId('db-filter-taxids').fill('40674');
      await submit(page);
      await waitStatus(page, 5, 'completed');
      {
        const { output, command } = await result(page, 5, 6);
        expect(command).toContain('-taxids 40674');
        const expected = runNative(['blastp', '-query', protQueryPath, '-db', 'tfp', '-taxids', '40674'], protDir, 6, {
          BLASTDB: protDir,
        });
        expect(output, 'blastp -db -taxids 40674 with taxonomy4blast.sqlite3 outfmt 6').toBe(expected);
      }
      await page.getByTestId('db-filter-taxids').fill('');

      // 3. BLASTN -filtering_db (`rpt_d4`) with `-db` (`wm_s60k`)
      const rptDir = join(work, 'rptdb');
      const wmFiles = copyMatching(join(LOSAT, 'tests/db_search'), rptDir, (name) => name.startsWith('wm_s60k.'));
      const rptFiles = copyMatching(join(LOSAT, 'tests/db_search'), rptDir, (name) => name.startsWith('rpt_d4.'));
      const wmQueryPath = join(LOSAT, 'tests/fixtures/blastn_regression/inputs/wm_query.fa');
      const wmQueryText = readFileSync(wmQueryPath, 'utf8');

      await program(page, 'blastn');
      await paste(page, 'query', wmQueryText);
      await setSubjectDbFiles(page, wmFiles);
      await expect(page.getByTestId('subject-db-selected')).toContainText('wm_s60k');
      await openParameters(page);
      await enableFilteringDb(page);
      await page.getByTestId('param-filtering_db-files-input').setInputFiles(rptFiles);
      await expect(page.getByTestId('param-filtering_db-selected')).toContainText('rpt_d4');
      await page.getByTestId('param-soft_masking').uncheck();
      await submit(page);
      await waitStatus(page, 6, 'completed');
      await expectFilteringDbParity(page, 6, ['-db', 'wm_s60k'], wmQueryPath, rptDir, 'blastn -db -filtering_db outfmt 6');

      // 4. BLASTN -filtering_db (`rpt_d4`) with `-subject` FASTA
      const wmSubjectPath = join(LOSAT, 'tests/fixtures/blastn_regression/inputs/s60k.fa');
      const wmSubjectText = readFileSync(wmSubjectPath, 'utf8');
      await page.getByTestId('subject-source-fasta').check();
      await paste(page, 'subject', wmSubjectText);
      await submit(page);
      await waitStatus(page, 7, 'completed');
      await expectFilteringDbParity(page, 7, ['-subject', wmSubjectPath], wmQueryPath, rptDir, 'blastn -subject -filtering_db outfmt 6');
    } finally {
      rmSync(work, { recursive: true, force: true });
    }
  });
});
