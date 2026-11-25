// Browser searches with -parse_deflines where a subject's Seq-id names the other molecule
// (line RG-XT item X7-5, ABI v2 version 14; docs/web/abi_v2.md §4 `losat_web2_scan_begin`, §9
// *validation* and *scan*): validation reports `parse_deflines`, the draft rescans the inputs
// with it so the record table matches the engine's `register_inputs` lengths, and the search
// completes with the native CLI's and NCBI's frozen outputs (`s1.bp.subj.*`, `s1.bn.subjx.*`).
import { spawnSync } from 'node:child_process';
import { mkdirSync, mkdtempSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { inputNames, openFiles, openParameters, openSearch, program, storedOutputs, submit, task, waitStatus } from './support/search';

test.setTimeout(180_000);
test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine reads no ids of the deflines: the engine build runs this test');

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

const PD_MOL = 'LOSAT/tests/fasta/outfmt0/pd_mol';
const FROZEN = 'LOSAT/tests/fixtures/outfmt0';
const pdMol = (name: string): Buffer => readFileSync(join(REPOSITORY, PD_MOL, name));

/** Frozen NCBI fixture bytes mapped into LOSAT's words (`losat_wording.py map-output`). */
function frozenInLosatWords(file: string): Buffer {
  const path = join(REPOSITORY, FROZEN, file);
  const mapped = spawnSync('python3', ['-I', join(REPOSITORY, 'LOSAT/tests/losat_wording.py'), 'map-output', path], {
    maxBuffer: 1 << 28,
  });
  if (mapped.error) throw mapped.error;
  if (mapped.status !== 0) throw new Error(`losat_wording.py map-output ${file}: ${mapped.stderr.toString('utf8')}`);
  return mapped.stdout;
}

/**
 * Replaces the manifest's wrapped subject path (`(Input:\ntests/fasta/outfmt0/pd_mol/<name>).`)
 * in a frozen outfmt 0 report with the page's uploaded file name (`(Input: <name>).`).
 */
function withUploadedSubjectName(outfmt0: Buffer, fileName: string): Buffer {
  const manifestPath = `(Input:\ntests/fasta/outfmt0/pd_mol/${fileName}).`;
  const uploaded = `(Input: ${fileName}).`;
  const text = outfmt0.toString('utf8');
  expect(text).toContain(manifestPath);
  return Buffer.from(text.replaceAll(manifestPath, uploaded), 'utf8');
}

function nativeRun(
  words: readonly string[],
  query: { readonly name: string; readonly bytes: Buffer },
  subject: { readonly name: string; readonly bytes: Buffer },
  format: 0 | 6 | 7,
): { status: number | null; stdout: Buffer; stderr: Buffer } {
  mkdirSync(test.info().outputPath(), { recursive: true });
  const directory = mkdtempSync(join(test.info().outputPath(), 'native-'));
  writeFileSync(join(directory, query.name), query.bytes);
  writeFileSync(join(directory, subject.name), subject.bytes);
  const args = [words[0]!, '-query', query.name, '-subject', subject.name, ...words.slice(1), '-outfmt', String(format), '-num_threads', '1'];
  const result = spawnSync(NATIVE!, args, { cwd: directory, maxBuffer: 1 << 28 });
  return { status: result.status, stdout: result.stdout, stderr: result.stderr };
}

async function checkRun(
  page: Page,
  words: readonly string[],
  query: { readonly name: string; readonly bytes: Buffer },
  subject: { readonly name: string; readonly bytes: Buffer },
  frozen: { readonly f0: string; readonly f6: string; readonly err?: string },
): Promise<void> {
  await submit(page);
  await waitStatus(page, 1, 'completed');
  const { outputs, warnings } = await storedOutputs(page, 1);
  const names = inputNames((await page.getByTestId('run-argv').textContent()) ?? '');
  expect(names).toEqual({ query: query.name, subject: subject.name });

  expect(outputs[6].equals(frozenInLosatWords(frozen.f6)), `${frozen.f6} (outfmt 6)`).toBe(true);
  expect(outputs[0].equals(withUploadedSubjectName(frozenInLosatWords(frozen.f0), subject.name)), `${frozen.f0} (outfmt 0)`).toBe(true);
  const expectedErr = frozen.err === undefined ? '' : frozenInLosatWords(frozen.err).toString('utf8');
  expect(warnings, 'warnings').toBe(expectedErr);

  if (NATIVE === undefined) {
    test.info().annotations.push({ type: 'skip', description: `The bytes of the outputs: ${NO_NATIVE_REASON}` });
    return;
  }
  for (const format of [0, 6, 7] as const) {
    const cli = nativeRun(words, query, subject, format);
    expect(cli.status, `-outfmt ${format}: ${cli.stderr.toString('utf8')}`).toBe(0);
    expect(outputs[format].equals(cli.stdout), `native -outfmt ${format}`).toBe(true);
    expect(warnings, `native -outfmt ${format} warnings`).toBe(cli.stderr.toString('utf8'));
  }
}

test('blastp -parse_deflines rescans a gnl|TI| protein subject and matches native CLI and frozen s1.bp.subj.*', async ({ page }) => {
  const query = { name: 'pq.fa', bytes: pdMol('pq.fa') };
  const subject = { name: 's_ti_prot.fa', bytes: pdMol('s_ti_prot.fa') };
  await program(page, 'blastp');
  await openFiles(page, 'query', [{ name: query.name, text: query.bytes }]);
  await openFiles(page, 'subject', [{ name: subject.name, text: subject.bytes }]);
  await expect(page.getByTestId('subject-source-0-summary')).toHaveText('2 records · 260 aa');
  await openParameters(page);
  await page.getByTestId('param-extra').fill('-parse_deflines');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await expect(page.getByTestId('subject-source-0-summary')).toHaveText('2 records · 217 aa', { timeout: 30_000 });
  await checkRun(page, ['blastp', '-parse_deflines'], query, subject, {
    f0: 's1.bp.subj.f0.out',
    f6: 's1.bp.subj.f6.out',
    err: 's1.bp.subj.f0.err',
  });
});

test('blastn -parse_deflines rescans an sp| nucleotide subject and matches native CLI and frozen s1.bn.subjx.*', async ({ page }) => {
  const query = { name: 'bnq.fa', bytes: pdMol('bnq.fa') };
  const subject = { name: 's_sp_nucx_bn.fa', bytes: pdMol('s_sp_nucx_bn.fa') };
  await program(page, 'blastn');
  await task(page, 'blastn');
  await openFiles(page, 'query', [{ name: query.name, text: query.bytes }]);
  await openFiles(page, 'subject', [{ name: subject.name, text: subject.bytes }]);
  await expect(page.getByTestId('subject-source-0-summary')).toHaveText('1 record · 287 nt');
  await openParameters(page);
  await page.getByTestId('param-extra').fill('-parse_deflines');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await expect(page.getByTestId('subject-source-0-summary')).toHaveText('1 record · 300 nt', { timeout: 30_000 });
  await checkRun(page, ['blastn', '-task', 'blastn', '-parse_deflines'], query, subject, {
    f0: 's1.bn.subjx.f0.out',
    f6: 's1.bn.subjx.f6.out',
  });
});

test('blastn -parse_deflines accepts an sp| subject with protein letters and junk line and matches native CLI', async ({ page }) => {
  const query = { name: 'bnq.fa', bytes: pdMol('bnq.fa') };
  const subject = { name: 's_sp_junk_bn.fa', bytes: pdMol('s_sp_junk_bn.fa') };
  await program(page, 'blastn');
  await task(page, 'blastn');
  await openFiles(page, 'query', [{ name: query.name, text: query.bytes }]);
  await openFiles(page, 'subject', [{ name: subject.name, text: subject.bytes }]);
  // Without -parse_deflines, the subject's E-only line is dropped and line 3 fails CheckDataLine.
  await expect(page.getByTestId('subject-source-0-check')).toHaveAttribute('data-check', 'refused', { timeout: 30_000 });
  await openParameters(page);
  await page.getByTestId('param-extra').fill('-parse_deflines');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await expect(page.getByTestId('subject-source-0-summary')).toHaveText('1 record · 6 nt', { timeout: 30_000 });
  await expect(page.getByTestId('subject-source-0-check')).toHaveAttribute('data-check', 'ok', { timeout: 30_000 });
  await submit(page);
  await waitStatus(page, 1, 'completed');
  const { outputs, warnings } = await storedOutputs(page, 1);
  const names = inputNames((await page.getByTestId('run-argv').textContent()) ?? '');
  expect(names).toEqual({ query: query.name, subject: subject.name });
  if (NATIVE === undefined) {
    test.info().annotations.push({ type: 'skip', description: `The bytes of the outputs: ${NO_NATIVE_REASON}` });
    return;
  }
  for (const format of [0, 6, 7] as const) {
    const cli = nativeRun(['blastn', '-task', 'blastn', '-parse_deflines'], query, subject, format);
    expect(cli.status, `-outfmt ${format}: ${cli.stderr.toString('utf8')}`).toBe(0);
    expect(outputs[format].equals(cli.stdout), `native -outfmt ${format}`).toBe(true);
    expect(warnings, `native -outfmt ${format} warnings`).toBe(cli.stderr.toString('utf8'));
  }
});
