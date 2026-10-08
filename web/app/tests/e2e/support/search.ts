// Helpers of the E2E tests that search through the real application (S12, S13): inputs
// pasted or opened as files, the program, the queue, and the stored outputs of a run as the
// results screen's Outputs view shows them.
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { REPOSITORY } from './harness-server';

export type Role = 'query' | 'subject';

/** A FASTA file of LOSAT/tests/fasta. */
export const fasta = (path: string): Buffer => readFileSync(join(REPOSITORY, 'LOSAT/tests/fasta', path));

/** Waits until input source `index` of a role is read, and its engine check (if any) has answered. */
export async function settled(page: Page, role: Role, index = 0): Promise<void> {
  const source = page.getByTestId(`${role}-source-${index}`);
  await expect(source).toHaveAttribute('data-status', /ready|failed/, { timeout: 30_000 });
  if ((await source.getAttribute('data-status')) === 'ready') {
    await expect(page.getByTestId(`${role}-source-${index}-check`)).not.toHaveAttribute('data-check', 'pending', {
      timeout: 30_000,
    });
  }
}

export async function paste(page: Page, role: Role, text: string): Promise<void> {
  await page.getByTestId(`${role}-input`).fill(text);
  await settled(page, role);
}

export async function openFiles(
  page: Page,
  role: Role,
  files: ReadonlyArray<{ name: string; text: string | Buffer }>,
): Promise<void> {
  const before = await page.locator(`[data-testid^="${role}-source-"][data-status]`).count();
  await page.getByTestId(`${role}-files`).setInputFiles(
    files.map((file) => ({ name: file.name, mimeType: 'text/plain', buffer: Buffer.from(file.text) })),
  );
  for (let i = 0; i < files.length; i++) await settled(page, role, before + i);
}

/** Opens the record list of a source (it starts open for a few records). */
export async function showRecords(page: Page, role: Role, index = 0): Promise<void> {
  await page
    .getByTestId(`${role}-source-${index}`)
    .locator('details.records')
    .evaluate((details) => ((details as HTMLDetailsElement).open = true));
}

export async function program(page: Page, id: string): Promise<void> {
  await page.getByTestId(`program-${id}`).check();
}

export async function submit(page: Page): Promise<void> {
  await expect(page.getByTestId('argv-validation')).not.toHaveAttribute('data-state', 'checking');
  await page.getByTestId('add-to-queue').click();
}

/** Waits until run `run` ends, and fails with its error if it ends otherwise than `status`. */
export async function waitStatus(
  page: Page,
  run: number,
  status: 'completed' | 'cancelled' | 'failed',
  timeout = 120_000,
): Promise<void> {
  const element = page.getByTestId(`run-${run}-status`);
  await expect(element).toHaveText(/^(completed|cancelled|failed)$/, { timeout });
  const actual = await element.textContent();
  if (actual !== status) {
    const error = (await page.getByTestId(`run-${run}-error`).textContent({ timeout: 1000 }).catch(() => null)) ?? '';
    throw new Error(`run ${run} ended ${actual}, not ${status}: ${error}`);
  }
}

/** Chooses run `run` in the results tab's run picker and waits until its results are read. */
export async function chooseRun(page: Page, run: number): Promise<void> {
  await page.getByTestId('tab-results').click();
  const select = page.getByTestId('results-run');
  const label = await select.locator('option', { hasText: new RegExp(`^Run ${run} · `) }).textContent();
  await select.selectOption({ label: label!.trim() });
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', String(run));
}

/**
 * Shows the stored output of one format of run `run` in the Outputs view (the run already
 * shown in the results tab, or chosen there), and waits until it is read.
 */
export async function showOutput(page: Page, run: number, format: 0 | 6 | 7, choose = true): Promise<void> {
  if (choose) await chooseRun(page, run);
  await page.getByTestId('results-view-outputs').click();
  await page.getByTestId(`format-${format}`).click();
  await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', `${run}:${format}`);
}

/** The CLI command of a completed run, and its stored output of one format. */
export async function result(page: Page, run: number, format: 0 | 6 | 7): Promise<{ command: string; output: string }> {
  await showOutput(page, run, format);
  const command = (await page.getByTestId('result-command').textContent()) ?? '';
  const output = (await page.getByTestId('result-output').textContent()) ?? '';
  await page.getByTestId('tab-search').click();
  return { command, output };
}
