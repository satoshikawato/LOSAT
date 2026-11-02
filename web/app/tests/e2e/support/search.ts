// Helpers of the E2E tests that search through the real application (S12, S13): inputs
// pasted or opened as files, the program, the runs, and the stored outputs of a run as the
// results screen's Outputs view shows them.
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { REPOSITORY } from './harness-server';

export type Role = 'query' | 'subject';

/** A FASTA file of LOSAT/tests/fasta. */
export const fasta = (path: string): Buffer => readFileSync(join(REPOSITORY, 'LOSAT/tests/fasta', path));

/**
 * Waits until the engine's describe of the chosen program has come (or failed). It comes after
 * the page opens (the Data worker starts its session and loads the serial reactor first: up to
 * a second in Firefox) or the program changes, and its help and task move every control below
 * "Program Selection" down. A click on such a control while it moves is lost: Firefox and
 * WebKit press and release some 20-300 ms apart, Playwright checks the target at the press
 * only, and the release lands on the panel (WF1, losat-baselines/web-flaky-20261016).
 */
export async function described(page: Page): Promise<void> {
  await expect(page.getByTestId('search-panel')).toHaveAttribute('data-describe', /^(ready|failed)$/, { timeout: 30_000 });
}

/** Opens the application at its search form, once the form no longer moves (`described`). */
export async function openSearch(page: Page): Promise<void> {
  await page.goto('/');
  await described(page);
}

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

/** Chooses a program, and waits for its describe (`described`). */
export async function program(page: Page, id: string): Promise<void> {
  await page.getByTestId(`program-${id}`).check();
  await described(page);
}

/** Opens "Algorithm parameters" (W4b: it starts closed), where most options are set. */
export async function openParameters(page: Page): Promise<void> {
  const toggle = page.getByTestId('algorithm-parameters-toggle');
  if ((await toggle.getAttribute('aria-expanded')) !== 'true') await toggle.click();
  await expect(page.getByTestId('parameter-form')).toBeVisible();
}

/** Chooses a task in "Program Selection" (W4b: radio buttons). */
export async function task(page: Page, value: string): Promise<void> {
  await page.getByTestId(`param-task-${value}`).check();
}

export async function submit(page: Page): Promise<void> {
  await expect(page.getByTestId('argv-validation')).not.toHaveAttribute('data-state', 'checking');
  await page.getByTestId('run-losat').click();
}

/** Waits until run `run` ends, and fails with its error if it ends otherwise than `status`. */
export async function waitStatus(
  page: Page,
  run: number,
  status: 'completed' | 'cancelled' | 'failed',
  timeout = 120_000,
): Promise<void> {
  const option = page.getByTestId(`run-${run}-option`);
  await expect(option).toHaveAttribute('data-status', /^(completed|cancelled|failed)$/, { timeout });
  const actual = await option.getAttribute('data-status');
  if (actual !== status) {
    const error = await runMessage(page, run);
    throw new Error(`run ${run} ended ${actual}, not ${status}: ${error}`);
  }
}

/** Selects run `run` in the results page's run picker (the page may show no results for it). */
export async function pickRun(page: Page, run: number): Promise<void> {
  await page.getByTestId('tab-results').click();
  const select = page.getByTestId('results-run');
  const label = await select.locator('option', { hasText: new RegExp(`^Run ${run} · `) }).textContent();
  await select.selectOption({ label: label!.trim() });
}

/** The text the results page gives for run `run` when it has no results (a failed or cancelled run). */
export async function runMessage(page: Page, run: number): Promise<string> {
  await pickRun(page, run);
  const status = page.getByTestId('results-status');
  await expect(status).toHaveAttribute('data-phase', 'unavailable');
  await expect(status).toContainText(`Run ${run} `);
  return ((await status.textContent()) ?? '').trim();
}

/** Chooses run `run` in the results tab's run picker and waits until its results are read. */
export async function chooseRun(page: Page, run: number): Promise<void> {
  await pickRun(page, run);
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
