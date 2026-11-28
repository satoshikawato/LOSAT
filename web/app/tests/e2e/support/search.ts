// Helpers of the E2E tests that search through the real application (S12, S13): inputs
// pasted or opened as files, the program, the runs, and the stored outputs of a run as the
// results screen's Outputs view shows them.
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { describeBrowserThreads } from './browser-threads';
import { REPOSITORY } from './harness-server';
import { navigate } from './navigate';
import { parseShellWords } from '../walk/model';

export type Role = 'query' | 'subject';

/** A FASTA file of LOSAT/tests/fasta. */
export const fasta = (path: string): Buffer => readFileSync(join(REPOSITORY, 'LOSAT/tests/fasta', path));

/** The panels whose `data-describe` tells whether the engine answered their describe: the search form and the tool tabs. */
export type DescribedPanel = 'search-panel' | 'wm-panel' | 'ms-panel' | 'mk-panel' | 'bf-panel';

/**
 * How long a test waits for the engine to answer a panel's describe (`described`). The first describe of a page starts
 * the engine (the Data worker's session and the serial reactor fetched and compiled), so its time follows the machine's
 * load: 30 s ran out for the search form under a local load average of 36-97 (masking-tool.spec.ts, WW4c R1 Firefox and
 * WW4c final WebKit) where the same test passed alone, and the tool tabs had waited 60 s.
 */
export const ENGINE_DESCRIBE_TIMEOUT = 60_000;

/**
 * Waits until the engine has answered the describe of `panel` (the search form's program, or a tool tab's `-xmlhelp`):
 * `answered` accepts either answer, `ready` or `failed` waits for that one. The one wait for the engine's describe in the
 * e2e tests (`ENGINE_DESCRIBE_TIMEOUT`). The search form's describe comes after the page opens (the Data worker starts
 * its session and loads the serial reactor first: up to a second in Firefox) or the mode changes, and its help and task
 * move every control below "Task Selection" down. A click on such a control while it moves is lost: Firefox and WebKit
 * press and release some 20-300 ms apart, Playwright checks the target at the press only, and the release lands on the
 * panel (WF1, losat-baselines/web-flaky-20261016).
 */
export async function described(
  page: Page,
  panel: DescribedPanel = 'search-panel',
  answer: 'answered' | 'ready' | 'failed' = 'answered',
): Promise<void> {
  await expect(page.getByTestId(panel)).toHaveAttribute('data-describe', answer === 'answered' ? /^(ready|failed)$/ : answer, {
    timeout: ENGINE_DESCRIBE_TIMEOUT,
  });
}

/** Opens the application at its search form, once the form no longer moves (`described`). */
export async function openSearch(page: Page): Promise<void> {
  await navigate(page, '/');
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

/** Opens the record list of a source (it starts closed). */
export async function showRecords(page: Page, role: Role, index = 0): Promise<void> {
  await page
    .getByTestId(`${role}-source-${index}`)
    .locator('details.records')
    .evaluate((details) => ((details as HTMLDetailsElement).open = true));
}

/** Chooses a program, and waits for its describe (`described`). */
export async function program(page: Page, id: string): Promise<void> {
  if (id === 'psiblast') {
    await page.getByTestId('program-blastp').check();
    await described(page);
    await page.getByTestId('param-task-psiblast').check();
    await described(page);
    return;
  }
  await page.getByTestId(`program-${id}`).check();
  await described(page);
}

/** Opens "Algorithm parameters" (W4b: it starts closed), where most options are set. */
export async function openParameters(page: Page): Promise<void> {
  const toggle = page.getByTestId('algorithm-parameters-toggle');
  if ((await toggle.getAttribute('aria-expanded')) !== 'true') await toggle.click();
  await expect(page.getByTestId('parameter-form')).toBeVisible();
}

/** Chooses a task in "Task Selection" (W4b: radio buttons). */
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
  try {
    await expect(option).toHaveAttribute('data-status', /^(completed|cancelled|failed)$/, { timeout });
  } catch (error) {
    // A run that did not end: whether the browser's threads run (a slow search) or wait, and where (browser-threads.ts).
    const threads = await describeBrowserThreads();
    if (threads === undefined || !(error instanceof Error)) throw error;
    error.message = `${error.message}\n\n${threads}`;
    throw error;
  }
  const actual = await option.getAttribute('data-status');
  if (actual !== status) {
    const error = await runMessage(page, run);
    throw new Error(`run ${run} ended ${actual}, not ${status}: ${error}`);
  }
}

/** Selects run `run` in the results page's run picker (the page may show no results for it). */
export async function pickRun(page: Page, run: number): Promise<void> {
  await page.getByTestId('tab-results').click();
  // By the option's value (the run's ID), not its label: the label ends in " (running)" until the
  // run ends, so a label read just before then matched no option, and "several runs" timed out
  // (CI of #731 and #732; found by the SX line).
  const value = await page.getByTestId(`run-${run}-option`).getAttribute('value');
  await page.getByTestId('results-run').selectOption(value!);
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

/**
 * The stored outputs of run `run` (outfmt 0, 6 and 7), saved from "Download All" as the bytes the
 * application stores, and the run's warnings (empty when the engine wrote none).
 */
export async function storedOutputs(page: Page, run: number): Promise<{ outputs: Record<0 | 6 | 7, Buffer>; warnings: string }> {
  await chooseRun(page, run);
  const outputs: Partial<Record<0 | 6 | 7, Buffer>> = {};
  for (const format of [0, 6, 7] as const) {
    await page.getByTestId('results-download-all').click();
    const [download] = await Promise.all([page.waitForEvent('download'), page.getByTestId('results-download-menu').getByTestId(`results-download-${format}`).click()]);
    outputs[format] = readFileSync((await download.path())!);
  }
  await page.getByTestId('results-view-details').click();
  await expect(page.getByTestId('run-details')).toBeVisible();
  const shown = await page.getByTestId('run-diagnostics').count();
  const warnings = shown === 0 ? '' : ((await page.getByTestId('run-diagnostics').textContent()) ?? '');
  return { outputs: outputs as Record<0 | 6 | 7, Buffer>, warnings };
}

/** The part of a file's name that is the UTC time of a run's "Run LOSAT" (`20261123T101530Z`). */
export const RUN_TIME = String.raw`\d{8}T\d{6}Z`;

/**
 * The source of a regular expression for the name that run `run` of `program` gives a file that the page makes for it (pasted text,
 * several files, a file with records left out; line WQA's D22, DECISIONS 44): the run's prefix of the names of its outputs, with the task
 * where the program has one, and the role (`losat-run3-blastn-megablast-20261123T101530Z.query.fa`).
 */
export function madeName(run: number, program: string, role: Role): string {
  return String.raw`losat-run${run}-${program}(?:-[a-z0-9-]+)?-${RUN_TIME}\.${role}\.fa`;
}

/** The names that a shown command (`losat PROGRAM -query A -subject B ...`, or the argv without `losat`) gives its inputs. */
export function inputNames(command: string): Record<Role, string> {
  const words = parseShellWords(command);
  return { query: words[words.indexOf('-query') + 1]!, subject: words[words.indexOf('-subject') + 1]! };
}

/** A regular expression that a whole command matches: `parts` are literal text, or `{ source }` for a part that is a regular expression. */
export function commandPattern(...parts: Array<string | { readonly source: string }>): RegExp {
  return new RegExp(`^${parts.map((part) => (typeof part === 'string' ? part.replace(/[.*+?^${}()|[\]\\]/g, '\\$&') : part.source)).join('')}$`);
}

/** The names that run `run`'s command (Run details) gives its query and subject: the user's file names, or the run's own for files that the page makes. */
export async function runInputNames(page: Page, run: number): Promise<Record<Role, string>> {
  await chooseRun(page, run);
  await page.getByTestId('results-view-details').click();
  return inputNames((await page.getByTestId('run-argv').textContent()) ?? '');
}
