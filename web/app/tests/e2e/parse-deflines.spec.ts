// -parse_deflines typed in Additional options runs as on the command line (line WQA's sweep, D7,
// 2026-11-19; docs/web/abi_v2.md §4 `losat_web2_register_inputs`): the engine reads the query and the
// subject of the run together for its argv, with the ids of the deflines in one scope, so the saved
// outputs and the warnings are the native CLI's bytes, the results name the records by the IDs that the
// outputs show (`AB000001.1` of `>gb|AB000001.1| ...`), and an input that the CLI refuses under the
// option fails the run with the CLI's message. The FakeEngine reads no ids, so the tests run in the
// engine build (LOSAT_WEB_REACTORS); with LOSAT_WEB_NATIVE the outputs are compared with the native
// CLI of the same commit.
import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { chooseRun, inputNames, openFiles, openParameters, openSearch, paste, pickRun, program, showOutput, storedOutputs, submit, waitStatus } from './support/search';

test.setTimeout(180_000);
test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine reads no ids of the deflines: the engine build runs these tests');

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

/** Pseudo-random letters of `alphabet` (a linear congruential generator): the same letters in every run. */
function letters(alphabet: string, seed: number, length: number): string {
  let state = seed;
  let out = '';
  for (let i = 0; i < length; i++) {
    state = (Math.imul(state, 1103515245) + 12345) >>> 0;
    out += alphabet[(state >>> 16) % alphabet.length];
  }
  return out;
}

const DNA = letters('ACGT', 7, 300);
const PROTEIN = letters('ACDEFGHIKLMNPQRSTVWY', 11, 160);
const CASES = [
  {
    program: 'blastn',
    query: `>gb|AB000001.1| first query\n${DNA}\n`,
    subject: `>ref|NM_000006.1| first subject\n${letters('ACGT', 9, 120)}${DNA.slice(30, 270)}${letters('ACGT', 10, 120)}\n`,
    ids: ['AB000001.1', 'NM_000006.1'],
  },
  {
    program: 'blastp',
    query: `>sp|P12345.1|ABC_X first query\n${PROTEIN}\n`,
    subject: `>gi|6|ref|NP_000006.1| first subject\n${letters('ACDEFGHIKLMNPQRSTVWY', 13, 40)}${PROTEIN.slice(10, 150)}\n`,
    ids: ['P12345.1', 'NP_000006.1'],
  },
] as const;

/** The native CLI's run of `words` (program first) on the inputs under the names that the run's command gives them, in `-outfmt format`, one thread. */
function nativeRun(
  words: readonly string[],
  query: string,
  subject: string,
  format: 0 | 6 | 7,
  names: { readonly query: string; readonly subject: string } = { query: 'query.fa', subject: 'subject.fa' },
): { status: number | null; stdout: Buffer; stderr: Buffer } {
  mkdirSync(test.info().outputPath(), { recursive: true });
  const directory = mkdtempSync(join(test.info().outputPath(), 'native-'));
  writeFileSync(join(directory, names.query), query);
  writeFileSync(join(directory, names.subject), subject);
  const args = [words[0]!, '-query', names.query, '-subject', names.subject, ...words.slice(1), '-outfmt', String(format), '-num_threads', '1'];
  const result = spawnSync(NATIVE!, args, { cwd: directory, maxBuffer: 1 << 28 });
  return { status: result.status, stdout: result.stdout, stderr: result.stderr };
}

test('-parse_deflines runs: the outputs are the CLI\'s and the results name the records as the outputs do', async ({ page }, testInfo) => {
  if (NATIVE === undefined) testInfo.annotations.push({ type: 'skip', description: `The bytes of the outputs: ${NO_NATIVE_REASON}` });
  for (const [index, { program: id, query, subject, ids }] of CASES.entries()) {
    const run = index + 1;
    if (run > 1) await page.getByTestId('tab-search').click();
    await program(page, id);
    await paste(page, 'query', query);
    await paste(page, 'subject', subject);
    await openParameters(page);
    await page.getByTestId('param-extra').fill('-parse_deflines');
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
    await submit(page);
    await waitStatus(page, run, 'completed');
    await chooseRun(page, run);
    await expect(page.getByTestId('results-query-id')).toHaveText(ids[0]);
    await expect(page.getByTestId('results-subject-id')).toHaveText(ids[1]);
    const { outputs, warnings } = await storedOutputs(page, run);
    expect(outputs[6].toString('utf8').split('\t').slice(0, 2), `${id}: the first row's IDs`).toEqual(ids);
    if (NATIVE === undefined) continue;
    // The command names the pasted inputs by the run's own names (outfmt 0 and 7 show the subject's).
    const names = inputNames((await page.getByTestId('run-argv').textContent()) ?? '');
    for (const format of [0, 6, 7] as const) {
      const cli = nativeRun([id, '-parse_deflines'], query, subject, format, names);
      expect(cli.status, `${id} -outfmt ${format}: ${cli.stderr}`).toBe(0);
      expect(outputs[format].equals(cli.stdout), `${id} -outfmt ${format}`).toBe(true);
      expect(warnings, `${id} -outfmt ${format}: the warnings`).toBe(cli.stderr.toString('utf8'));
    }
  }
});

// -parse_deflines typed in Additional options has the engine read the inputs again, with the other reader (line RG-XT
// X7-5): until they are read, their blocks and cards keep their place (line WQA's D25, 2026-11-27). Their lines went and
// came back, so "Run LOSAT" below them moved up (143 px with these inputs) and back, and a press whose release came
// while they were read was released beside the button: no run started, and the page said nothing. Here the press
// starts, the option comes in the field, and the press is released once the query is read again.
test('Run LOSAT stays under a press while -parse_deflines has the inputs read again, and the press starts the run', async ({ page }) => {
  const { program: id, query, subject } = CASES[1];
  await program(page, id);
  await paste(page, 'query', query);
  await paste(page, 'subject', subject);
  await openParameters(page);
  const button = page.getByTestId('run-losat');
  await button.scrollIntoViewIfNeeded();
  const box = (await button.boundingBox())!;
  const point = { x: box.x + box.width / 2, y: box.y + box.height / 2 };
  await page.mouse.move(point.x, point.y);
  await page.mouse.down();
  // The button's place in the page's layout when the press holds, and at each change of an input's state until its
  // click: its row's place and its own box in the row. (Its pressed look, scale(0.97), is not in its layout and is not
  // measured.) And whether the press's point is still on the button.
  await page.evaluate((at) => {
    const place = () => {
      const target = document.querySelector<HTMLElement>('[data-testid="run-losat"]')!;
      const row = target.parentElement!.getBoundingClientRect();
      return `row ${row.left},${row.top} button ${target.offsetLeft},${target.offsetTop} ${target.offsetWidth}x${target.offsetHeight}`;
    };
    const on = () => document.elementFromPoint(at.x, at.y)?.closest('[data-testid="run-losat"]') !== null;
    const places = [place()];
    const hits = [on()];
    const reading = () => document.querySelector('[data-testid="query-source-0"]')?.getAttribute('data-status') === 'indexing';
    let readAgain: () => void = () => undefined;
    const queryReadAgain = new Promise<void>((resolve) => (readAgain = resolve));
    const observer = new MutationObserver(() => {
      places.push(place());
      hits.push(on());
      if (reading()) readAgain();
    });
    observer.observe(document.body, { subtree: true, attributes: true, attributeFilter: ['data-status', 'data-check'] });
    window.addEventListener('click', () => observer.disconnect(), { capture: true, once: true });
    Object.assign(window, { readAgainTest: { places, hits, queryReadAgain } });
  }, point);
  // The option comes in the field as typed just before the press: not by a focus of the field, which (in Chromium and
  // Firefox) scrolls the field, far below the button, into view, so that the test itself would move the button.
  await page.getByTestId('param-extra').evaluate((field: HTMLInputElement) => {
    field.value = '-parse_deflines';
    field.dispatchEvent(new Event('input', { bubbles: true }));
  });
  await page.evaluate(() => (window as unknown as { readAgainTest: { queryReadAgain: Promise<void> } }).readAgainTest.queryReadAgain);
  await page.mouse.up();
  const { places, hits } = await page.evaluate(
    () => (window as unknown as { readAgainTest: { places: string[]; hits: boolean[] } }).readAgainTest,
  );
  expect(places.length, 'the button\'s place at the press and at the changes of the inputs').toBeGreaterThan(1);
  expect(new Set(places), `the button's places: ${places.join('; ')}`).toEqual(new Set([places[0]]));
  expect(new Set(hits), 'the press\'s point on the button').toEqual(new Set([true]));
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await expect(page.getByTestId('results-run-options')).toContainText('-parse_deflines');
});

// The CLI reads the subjects first, and two subject records of one id stop it under -parse_deflines
// (LOSAT's refusal: NCBI BLAST+ shows the first for both). When both records are in one subject
// source, the form's input check (losat_web2_register with -parse_deflines) refuses that source with
// the CLI's message; when they come from two combined sources (each accepted on its own), the run
// fails with that message.
test('an input that the CLI refuses under -parse_deflines is refused in the form or fails the run with the CLI\'s message', async ({ page }, testInfo) => {
  const query = `>q1 query\n${DNA}\n`;
  const first = `>a one\n${DNA.slice(0, 150)}\n`;
  const second = `>a two\n${DNA.slice(150)}\n`;
  const subject = `${first}${second}`;
  const message = "two subject records have the Seq-id lcl|a under -parse_deflines (NCBI BLAST+ shows the first for both); this is not supported by LOSAT's BLASTN";
  await program(page, 'blastn');
  await paste(page, 'query', query);
  await paste(page, 'subject', subject);
  await openParameters(page);
  await page.getByTestId('param-extra').fill('-parse_deflines');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  const singleCheck = page.getByTestId('subject-source-0-check');
  await expect(singleCheck).toHaveAttribute('data-check', 'refused', { timeout: 30_000 });
  await expect(singleCheck).toContainText(message);
  await page.getByTestId('subject-source-0').getByRole('button', { name: 'Remove' }).click();
  await openFiles(page, 'subject', [
    { name: 's1.fa', text: first },
    { name: 's2.fa', text: second },
  ]);
  await expect(page.getByTestId('subject-source-0-check')).toHaveAttribute('data-check', 'ok');
  await expect(page.getByTestId('subject-source-1-check')).toHaveAttribute('data-check', 'ok');
  await submit(page);
  await waitStatus(page, 1, 'failed');
  await pickRun(page, 1);
  const status = page.getByTestId('results-status');
  await expect(status).toContainText('Run 1 failed: ');
  await expect(status).toContainText(message);
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The CLI's error: ${NO_NATIVE_REASON}` });
    return;
  }
  const cli = nativeRun(['blastn', '-parse_deflines'], query, subject, 0);
  expect(cli.status).toBe(1);
  expect(cli.stderr.toString('utf8')).toBe(`Error: ${message}\n`);
});

// One record as the query and the subject, with lower-case letters and -lcase_masking: LOSAT refuses
// outfmt 0 of a record that is both a query and a subject under -parse_deflines and has masks, that
// format's alone (round 2 of line WQA's sweep, D7, 2026-11-19). With -outfmt 6 the run writes the
// CLI's outfmt 6 and 7, outfmt 0 stays empty and the warnings end with the CLI's error of outfmt 0;
// without -outfmt (outfmt 0) the run fails with that error, as the CLI fails.
test('a record that is both a query and a subject and has masks: outfmt 0 alone is refused, as the CLI refuses it', async ({ page }, testInfo) => {
  const record = `>gb|AB000001.1| shared record\n${DNA.slice(0, 40).toLowerCase()}${DNA.slice(40)}\n`;
  const words = ['blastn', '-parse_deflines', '-lcase_masking'] as const;
  const message = "the record gb|AB000001.1| is both a query and a subject under -parse_deflines and has masks (NCBI BLAST+ draws them at the other record's places in outfmt 0); this is not supported by LOSAT's BLASTN";
  for (const [run, extra] of [[1, ' -outfmt 6'], [2, '']] as const) {
    if (run > 1) await page.getByTestId('tab-search').click();
    await program(page, 'blastn');
    await paste(page, 'query', record);
    await paste(page, 'subject', record);
    await openParameters(page);
    // -lcase_masking is a control of the form ("Mask lower case letters"), not an additional option.
    const mask = page.getByTestId('param-lcase_masking');
    if (!(await mask.isChecked())) await mask.check();
    await page.getByTestId('param-extra').fill(`-parse_deflines${extra}`);
    await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
    await submit(page);
    await waitStatus(page, run, run === 1 ? 'completed' : 'failed');
  }
  await pickRun(page, 2);
  await expect(page.getByTestId('results-status')).toContainText(message);
  // The results of run 1 say that outfmt 0 could not be written, from the engine's report of the outputs that
  // failed (ABI v2 §5, stream 9; round 3 of D7), not that it shows the alignments of a first 0 subjects.
  await chooseRun(page, 1);
  const notices = page.getByTestId('results-notice');
  await expect(notices.and(page.locator('[data-kind="outfmt0-failed"]'))).toHaveText(
    "outfmt 0 could not be written (its error is under Warnings in Run details), so no alignment is linked to this query's 1 subject. outfmt 6 and 7 hold the HSPs.",
  );
  await expect(notices.and(page.locator('[data-kind="outfmt0-partial"]'))).toHaveCount(0);
  await page.getByTestId('results-view-hits').click();
  await expect(page.getByTestId('detail-not-in-outfmt0')).toHaveText(
    'No outfmt 0 alignment is linked to this HSP: outfmt 0 could not be written, and its error is under Warnings in Run details. The outfmt 6 row below and outfmt 7 hold the HSP.',
  );
  await showOutput(page, 1, 0, false);
  await expect(page.getByTestId('outputs-not-written')).toHaveText('outfmt 0 could not be written: its error is under Warnings in Run details.');
  await showOutput(page, 1, 6, false);
  await expect(page.getByTestId('outputs-not-written')).toHaveCount(0);
  const { outputs, warnings } = await storedOutputs(page, 1);
  expect(outputs[0].length, 'outfmt 0 is refused').toBe(0);
  expect(outputs[6].toString('utf8').split('\t').slice(0, 2)).toEqual(['AB000001.1', 'AB000001.1']);
  if (NATIVE === undefined) {
    testInfo.annotations.push({ type: 'skip', description: `The bytes of the outputs: ${NO_NATIVE_REASON}` });
    return;
  }
  const names = inputNames((await page.getByTestId('run-argv').textContent()) ?? '');
  const zero = nativeRun(words, record, record, 0, names);
  expect(zero.status).toBe(1);
  expect(zero.stderr.toString('utf8')).toBe(`Error: ${message}\n`);
  for (const format of [6, 7] as const) {
    const cli = nativeRun(words, record, record, format, names);
    expect(cli.status, `-outfmt ${format}: ${cli.stderr}`).toBe(0);
    expect(outputs[format].equals(cli.stdout), `-outfmt ${format}`).toBe(true);
    if (format === 6) expect(warnings, 'the warnings: outfmt 6\'s, then the error of outfmt 0').toBe(`${cli.stderr.toString('utf8')}Error: ${message}\n`);
  }
});
