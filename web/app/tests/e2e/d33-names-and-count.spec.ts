// Line WQA's D33 (2026-11-28), two low defects of sweep 4, through the real application (the engine build, LOSAT_WEB_REACTORS):
//  - S4-3: a list file of the Database filters whose name starts or ends with a space is run, and shown, under that name,
//    quoted as a query file is (the shown command runs on the command line with the same bytes);
//  - S4-8: Run LOSAT's "(N runs)" never shows the count of an older argv after a help or version word is typed or taken out;
//    the WindowMasker and Masking tabs' buttons say one run for a help or version word, as the press starts one (R1 M3, R2 H1),
//    and the line under the Search tab's button changes with its label (R2 L-c);
//  - R2 H2: in the tool tabs a file whose name starts or ends with a space is named so in the argv, as its run directory holds
//    it, so the run opens it and its output is the native command line's for the command shown (WindowMasker -in, -ustat,
//    -ids; DustMasker -in; Format Archive -archive); R2 L-d: Run details and the Results summary quote a query file's name.
import { readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { nativeRun } from './support/makeblastdb-native';
import { NATIVE } from './support/native';
import { described, openFiles, openParameters, openSearch, paste, pickRun, program, submit, waitStatus } from './support/search';

test.setTimeout(240_000);
test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine runs no database search and prints no help or version: the engine build runs these tests');

const TESTS = join(REPOSITORY, 'LOSAT/tests/db_search');

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

/**
 * Two list files of the Database filters, each run by itself (the command line refuses -seqidlist with -negative_seqidlist,
 * -gilist and the others together; the unit tests name every list file): the flag, the box's test id, and a name that starts
 * or ends with a space.
 */
const LISTS = [
  ['-seqidlist', 'seqidlist', 'ids.txt ', 'seqidlist/sl_tfp.txt'],
  ['-negative_seqidlist', 'negative_seqidlist', ' neg ids.txt', 'seqidlist/sl_tfp_neg.txt'],
] as const;

test('a list file whose name starts or ends with a space is run and shown under its name, quoted (S4-3)', async ({ page }) => {
  const database = readdirSync(join(TESTS, 'taxfilter'))
    .filter((name) => name.startsWith('tfp.'))
    .map((name) => ({ name, mimeType: 'application/octet-stream', buffer: readFileSync(join(TESTS, 'taxfilter', name)) }));
  await program(page, 'blastp');
  await paste(page, 'query', readFileSync(join(TESTS, 'taxfilter/q.faa'), 'utf8'));
  await page.getByTestId('subject-source-database').check();
  await page.getByTestId('subject-db-files-input').setInputFiles(database);
  await expect(page.getByTestId('subject-db-selected')).toContainText('tfp');
  await page.getByTestId('subject-db-filters-toggle').click();
  for (const [run, [flag, id, name, source]] of LISTS.entries()) {
    if (run > 0) {
      await page.getByTestId('tab-search').click();
      await page.getByTestId(`db-filter-${LISTS[run - 1]![1]}-remove`).click();
    }
    await page.getByTestId(`db-filter-${id}-input`).setInputFiles({ name, mimeType: 'text/plain', buffer: readFileSync(join(TESTS, source)) });
    await expect(page.getByTestId(`db-filter-${id}-name`)).toContainText(name.trim());
    await submit(page);
    await waitStatus(page, run + 1, 'completed');
    await page.getByTestId('tab-results').click();
    await page.getByTestId('results-view-details').click();
    // The command names the file as the user's file is called, quoted: it runs on the command line next to the file.
    const shown = (await page.getByTestId('run-argv').textContent()) ?? '';
    expect(shown, flag).toContain(`${flag} '${name}'`);
    // The line of the database's filter lists names the file the same way (R1, L1).
    await expect(page.getByTestId('run-database-filters')).toContainText(`${flag} '${name}'`);
  }
});

test('the button names no count of an older argv after a help or version word is typed or taken out (S4-8)', async ({ page }) => {
  await program(page, 'blastn');
  await openFiles(page, 'query', [{ name: 'a.fa', text: '>a\nACGTACGTACGTACGTACGTAAAC\n' }]);
  await openFiles(page, 'subject', [
    { name: 's1.fa', text: '>s1\nACGTACGTACGTACGTACGTAAAC\n' },
    { name: 's2.fa', text: '>s2\nTTTACGTACGTACGTACGTACGT\n' },
  ]);
  await page.getByTestId('subject-mode-separate').check();
  await openParameters(page);
  const button = page.getByTestId('run-losat');
  await expect(button).toHaveText('Run LOSAT (2 runs)');
  // A sampler of the button's text and of the line under it in every frame, from the edit to the engine's answer.
  // `from` is where the samples start that come after the edit's input event.
  await page.evaluate(() => {
    const w = window as unknown as { __labels: string[]; __lines: string[]; __widths: number[]; __from: number };
    w.__labels = [];
    w.__lines = [];
    w.__widths = [];
    w.__from = 0;
    document.addEventListener('input', () => void (w.__from = w.__labels.length), true);
    const sample = () => {
      const button = document.querySelector('[data-testid="run-losat"]');
      w.__labels.push((button?.textContent ?? '').trim());
      w.__lines.push((document.querySelector('[data-testid="search-summary-line"]')?.textContent ?? '').trim());
      w.__widths.push(button?.getBoundingClientRect().width ?? 0);
      requestAnimationFrame(sample);
    };
    sample();
  });
  const sampled = async (edit: () => Promise<void>, settle: string): Promise<string[]> => {
    await edit();
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
    await expect(button).toHaveText(settle);
    return page.evaluate(() => {
      const w = window as unknown as { __labels: string[]; __lines: string[]; __widths: number[]; __from: number };
      // Nothing moves under a press: the button is as wide with every label.
      const widths = new Set(w.__widths.slice(Math.max(0, w.__from - 1)).map((width) => Math.round(width)));
      // The line under the button and the label read one answer (R2, L-c): no frame counts the runs of a search while the
      // line says that the run prints.
      const apart = w.__labels
        .map((label, index) => [label, w.__lines[index] ?? ''] as const)
        .slice(w.__from)
        .filter(([label, line]) => label.includes('runs') && line.startsWith("Print the command line's"));
      return [
        ...w.__labels.slice(w.__from),
        ...(widths.size === 1 ? [] : [`WIDTHS ${[...widths].join(',')}`]),
        ...apart.map(([label, line]) => `APART ${label} | ${line}`),
      ];
    });
  };
  for (const word of ['-version', '-h']) {
    const typed = await sampled(() => page.getByTestId('param-extra').fill(word), 'Run LOSAT');
    expect(typed.filter((label) => label.includes('runs') || label.startsWith('WIDTHS')), `after ${word}`).toEqual([]);
    await expect(page.getByTestId('search-summary-line')).toContainText("Print the command line's");
    // Taken out: the label says no count until the answer, then the 2 runs again, in a button as wide as before, and the line
    // under it changes in the same frame as the label.
    const taken = await sampled(() => page.getByTestId('param-extra').fill(''), 'Run LOSAT (2 runs)');
    expect(taken.filter((label) => label.startsWith('WIDTHS') || label.startsWith('APART')), `after ${word} was taken out`).toEqual([]);
  }
  // An edit that does not print: the count is back after a round trip, and the button never changes its width.
  const ordinary = await sampled(() => page.getByTestId('param-extra').fill('-outfmt 17'), 'Run LOSAT (2 runs)');
  expect(ordinary.filter((label) => label.startsWith('WIDTHS')), 'after -outfmt 17').toEqual([]);
  await page.getByTestId('param-extra').fill('');
  await expect(button).toHaveText('Run LOSAT (2 runs)');
  // A press at once after the word starts one run, as the label that waits for the answer says.
  await page.getByTestId('param-extra').fill('-version');
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await expect(page.getByTestId('run-2-option')).toHaveCount(0);
});

test('the WindowMasker button names no count for a help or version word, which starts one print request (S4-8, R1 M3)', async ({ page }) => {
  await page.getByTestId('tab-windowmasker').click();
  await expect(page.getByTestId('wm-panel')).toHaveAttribute('data-describe', 'ready', { timeout: 60_000 });
  const files = ['a', 'b'].map((name) => ({
    name: `${name}.fa`,
    mimeType: 'application/octet-stream',
    buffer: Buffer.from(`>${name}\nACGTACGTACGTACGTACGTAAAC\n`),
  }));
  await page.getByTestId('wm-fasta-files').setInputFiles(files);
  await expect(page.locator('[data-testid="wm-fasta-list"] > li')).toHaveCount(2);
  await page.getByTestId('wm-fasta-mode-separate').check();
  const button = page.getByTestId('wm-run-losat');
  await expect(button).toHaveText('Run LOSAT (2 runs)');
  const toggle = page.getByTestId('wm-options-toggle');
  if ((await toggle.getAttribute('aria-expanded')) !== 'true') await toggle.click();
  await page.getByTestId('wm-param-extra').fill('-version');
  await expect(button).toHaveText('Run LOSAT');
  await page.getByTestId('wm-param-extra').fill('');
  await expect(button).toHaveText('Run LOSAT (2 runs)');
  await page.getByTestId('wm-param-extra').fill('-version');
  await expect(page.getByTestId('wm-argv-validation')).not.toHaveAttribute('data-state', 'checking');
  await button.click();
  await waitStatus(page, 1, 'completed');
  await expect(page.getByTestId('run-2-option')).toHaveCount(0);
});

test('the Masking button names no count for a help or version word, which starts one print request (S4-8, R2 H1)', async ({ page }) => {
  await page.getByTestId('tab-masking').click();
  await described(page, 'ms-panel', 'ready');
  const files = ['a', 'b'].map((name) => ({
    name: `${name}.fa`,
    mimeType: 'application/octet-stream',
    buffer: Buffer.from(`>${name}\nACGTACGTACGTACGTACGTAAAC\n`),
  }));
  await page.getByTestId('ms-fasta-files').setInputFiles(files);
  await expect(page.locator('[data-testid="ms-fasta-list"] > li')).toHaveCount(2);
  await page.getByTestId('ms-fasta-mode-separate').check();
  const button = page.getByTestId('ms-run-losat');
  await expect(button).toHaveText('Run LOSAT (2 runs)');
  const toggle = page.getByTestId('ms-options-toggle');
  if ((await toggle.getAttribute('aria-expanded')) !== 'true') await toggle.click();
  for (const word of ['-version-full', '-h']) {
    await page.getByTestId('ms-param-extra').fill(word);
    await expect(button, word).toHaveText('Run LOSAT');
    await page.getByTestId('ms-param-extra').fill('');
    await expect(button, word).toHaveText('Run LOSAT (2 runs)');
  }
  await page.getByTestId('ms-param-extra').fill('-version-full');
  await expect(page.getByTestId('ms-argv-validation')).not.toHaveAttribute('data-state', 'checking');
  await button.click();
  await waitStatus(page, 1, 'completed');
  await expect(page.getByTestId('run-2-option')).toHaveCount(0);
});

/** Runs a tool tab's form (its Run LOSAT button) as run `run` of the session, until it exits 0. */
async function runTool(page: Page, prefix: 'wm' | 'ms' | 'bf', run: number): Promise<void> {
  await expect(page.getByTestId(`${prefix}-argv-validation`)).toHaveAttribute('data-state', 'ok', { timeout: 30_000 });
  await page.getByTestId(`${prefix}-run-losat`).click();
  await waitStatus(page, run, 'completed');
  await expect(page.getByTestId('tool-run-view')).toHaveAttribute('data-run', String(run));
}

/** The bytes of the standard output of the tool run shown (its first output). */
async function standardOutput(page: Page): Promise<Buffer> {
  const [download] = await Promise.all([page.waitForEvent('download'), page.getByTestId('tool-output-download').click()]);
  const path = await download.path();
  if (path === null) throw new Error('download path was null');
  return readFileSync(path);
}

test('a tool tab names a file whose name starts or ends with a space as the run directory holds it, and runs it as the command line does (R2 H2)', async ({
  page,
}) => {
  test.skip(NATIVE === undefined, 'compares with the native command line of the same commit (LOSAT_WEB_NATIVE)');
  const fasta = Buffer.from('>a\nACGTACGTACGTACGTACGTAAACGGTTACCAGGTTAACCGGTTAAGG\n>b\nTTTACGTACGTACGTACGTACGTCCAAGGTTCC\n');
  const ids = Buffer.from('a\n');
  const counts = nativeRun(['windowmasker', '-mk_counts', '-in', 'c.fa', '-out', 'c.ascii'], new Map([['c.fa', fasta]]), 'UTC').files.get('c.ascii')!;
  const file = (name: string, buffer: Buffer) => ({ name, mimeType: 'application/octet-stream', buffer });

  // WindowMasker: Make counts of 'a.fa ', then Mask with the counts file ' c.ascii' and the -ids file ' ids.txt'.
  await page.getByTestId('tab-windowmasker').click();
  await described(page, 'wm-panel', 'ready');
  await page.getByTestId('wm-fasta-files').setInputFiles([file('a.fa ', fasta)]);
  await expect(page.getByTestId('wm-command')).toHaveText("losat windowmasker -mk_counts -in 'a.fa ' -out wmasker.ascii");
  await runTool(page, 'wm', 1);
  await page.getByTestId('tab-windowmasker').click();
  await page.getByTestId('wm-stage-mask').check();
  await page.getByTestId('wm-ustat-file').setInputFiles([file(' c.ascii', counts)]);
  const toggle = page.getByTestId('wm-options-toggle');
  if ((await toggle.getAttribute('aria-expanded')) !== 'true') await toggle.click();
  await page.getByTestId('wm-param-ids-file').setInputFiles([file(' ids.txt', ids)]);
  const masked = ['windowmasker', '-ustat', ' c.ascii', '-in', 'a.fa ', '-ids', ' ids.txt'];
  await expect(page.getByTestId('wm-command')).toHaveText("losat windowmasker -ustat ' c.ascii' -in 'a.fa ' -ids ' ids.txt'");
  await runTool(page, 'wm', 2);
  const wm = await standardOutput(page);
  const wmNative = nativeRun(masked, new Map([[' c.ascii', counts], ['a.fa ', fasta], [' ids.txt', ids]]), 'UTC');
  expect(wmNative.status).toBe(0);
  expect(wm.toString(), 'windowmasker -ustat -in -ids').toBe(wmNative.stdout.toString());

  // DustMasker: -in ' b.fa'.
  await page.getByTestId('tab-masking').click();
  await described(page, 'ms-panel', 'ready');
  await page.getByTestId('ms-fasta-files').setInputFiles([file(' b.fa', fasta)]);
  await expect(page.getByTestId('ms-command')).toHaveText("losat dustmasker -in ' b.fa'");
  await runTool(page, 'ms', 3);
  const dust = await standardOutput(page);
  const dustNative = nativeRun(['dustmasker', '-in', ' b.fa'], new Map([[' b.fa', fasta]]), 'UTC');
  expect(dustNative.status).toBe(0);
  expect(dust.toString(), 'dustmasker -in').toBe(dustNative.stdout.toString());

  // Format Archive: -archive 'bn.asn '.
  const archive = readFileSync(join(REPOSITORY, 'LOSAT/tests/blast_formatter/in/bn.asn'));
  await page.getByTestId('tab-blast-formatter').click();
  await described(page, 'bf-panel', 'ready');
  await page.getByTestId('bf-archive-file-input').setInputFiles([file('bn.asn ', archive)]);
  await expect(page.getByTestId('bf-command')).toHaveText("losat blast_formatter -archive 'bn.asn ' -outfmt 8");
  await runTool(page, 'bf', 4);
  const formatted = await standardOutput(page);
  const formattedNative = nativeRun(['blast_formatter', '-archive', 'bn.asn ', '-outfmt', '8'], new Map([['bn.asn ', archive]]), 'UTC');
  expect(formattedNative.status).toBe(0);
  expect(formatted.toString(), 'blast_formatter -archive').toBe(formattedNative.stdout.toString());
});

test("Run details and the Results summary quote a query file's name as the command does (R2 L-d)", async ({ page }) => {
  await program(page, 'psiblast');
  await page.getByTestId('query-source-msa').check();
  await page.getByTestId('query-file-input').setInputFiles({
    name: 'seed.aln ',
    mimeType: 'text/plain',
    buffer: readFileSync(join(REPOSITORY, 'LOSAT/tests/fixtures/msa/ps6_b_db_0.aln')),
  });
  await page.getByTestId('query-msa-master-idx').fill('1');
  await paste(
    page,
    'subject',
    '>s1\nTMAGRVELVTGPMFAGKSTYLKNIYQQENGGNKHCLFVKHSLETRYGCGTGTIVTHAGEV\n>s2\nTNIGKIHLVIGPMFAGKSSTIIKIYQKKIKMQKRCVFIKHVKDNRYGCELTETVTHSGQT\n',
  );
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await pickRun(page, 1);
  await expect(page.getByTestId('results-query-file')).toHaveText("-in_msa 'seed.aln '");
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-query-file')).toContainText("-in_msa 'seed.aln ' (");
  await expect(page.getByTestId('run-argv')).toContainText("-in_msa 'seed.aln '");
});
