// S16 onboarding through the real application: "Try an example" (the examples of
// docs/web/examples.md put in the search form, BLASTX disabled until SX, the question before
// inputs are replaced), "Getting started" (open on the first visit, closed once closed, and
// with a localStorage that throws), and the browser self-test. The FakeEngine build checks the
// form and the panel; the engine build (LOSAT_WEB_REACTORS) runs each example and the
// self-test, which must pass every case in every browser. With LOSAT_WEB_SCREENS it records the
// onboarding states for the screen review.
import { mkdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, test, type BrowserContext, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { paste, program, result, settled, submit, waitStatus } from './support/search';

test.setTimeout(300_000);

const RUNNABLE = ['blastn', 'blastp', 'tblastn', 'tblastx'] as const;
const LABELS: Readonly<Record<string, string>> = { blastn: 'BLASTN', blastp: 'BLASTP', tblastn: 'TBLASTN', tblastx: 'TBLASTX' };
const GETTING_STARTED_KEY = 'losat-web:getting-started';

/** The record IDs of an example file, as its definition lines give them. */
function recordIds(id: string, role: 'query' | 'subject'): string[] {
  const text = readFileSync(join(REPOSITORY, `web/app/src/examples/example-${id}-${role}.fasta`), 'utf8');
  return text.split('\n').filter((line) => line.startsWith('>')).map((line) => line.slice(1).split(' ')[0]!);
}

async function tryExample(page: Page, id: string): Promise<void> {
  await program(page, id);
  await page.getByTestId('try-example').click();
  await expect(page.getByTestId('example-loaded')).toBeVisible();
  await settled(page, 'query');
  await settled(page, 'subject');
}

/** Later visits start with "Getting started" closed (tests that are not about it). */
async function closedGettingStarted(context: BrowserContext): Promise<void> {
  await context.addInitScript((key) => {
    try {
      localStorage.setItem(key, 'closed');
    } catch {
      // A test of a storage that throws installs its own script.
    }
  }, GETTING_STARTED_KEY);
}

test.describe('Try an example', () => {
  test.beforeEach(async ({ context, page }) => {
    await closedGettingStarted(context);
    await page.goto('/');
  });

  for (const id of RUNNABLE) {
    test(`${LABELS[id]}: the example's files are the inputs, with their records, and nothing is searched`, async ({ page }) => {
      await tryExample(page, id);
      await expect(page.getByTestId('example-loaded')).toContainText(
        `Example loaded: ${LABELS[id]} of example-${id}-query.fasta against example-${id}-subject.fasta`,
      );
      await expect(page.getByTestId('example-loaded')).toContainText('Choose Run LOSAT to search.');
      for (const role of ['query', 'subject'] as const) {
        const ids = recordIds(id, role);
        await expect(page.getByTestId(`${role}-source-0`)).toContainText(`example-${id}-${role}.fasta`);
        await expect(page.getByTestId(`${role}-source-1`)).toHaveCount(0);
        await expect(page.getByTestId(`${role}-input`)).toHaveValue('');
        await expect(page.getByTestId(`${role}-source-0-summary`)).toContainText(`${ids.length} ${ids.length === 1 ? 'record' : 'records'}`);
        for (const [index, recordId] of ids.entries()) {
          await expect(page.getByTestId(`${role}-source-0-record-${index}`)).toHaveAttribute('aria-label', `Include record ${index + 1}, ${recordId}`);
        }
      }
      await expect(page.getByTestId(`program-${id}`)).toBeChecked();
      // Default settings, and no search: the queue is empty.
      await expect(page.getByTestId('threads')).toHaveValue('auto');
      await expect(page.getByTestId('job-title')).toHaveValue('');
      await expect(page.getByTestId('run-1')).toHaveCount(0);
      await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
    });
  }

  test('BLASTX: "Try an example" is disabled, with the reason', async ({ page }) => {
    await program(page, 'blastx');
    const button = page.getByTestId('try-example');
    await expect(button).toBeDisabled();
    await expect(button).toHaveAccessibleDescription(/BLASTX joins LOSAT Web after its certification/);
    await expect(page.getByTestId('program-unavailable')).toContainText('BLASTX joins LOSAT Web after its certification');
  });

  test('asks before it replaces inputs, and keeps them when the user says so', async ({ page }) => {
    await paste(page, 'query', '>mine my own query\nACGTACGTACGTTTGACCATGG\n');
    await page.getByTestId('threads').selectOption('2');
    await page.getByTestId('job-title').fill('My search');
    await page.getByTestId('try-example').click();
    const question = page.getByTestId('try-example-question');
    await expect(question).toContainText('Replace the query and subject in the form with the BLASTN example?');
    await expect(page.getByTestId('try-example-confirm')).toBeFocused();
    await page.getByTestId('try-example-cancel').click();
    await expect(question).toHaveCount(0);
    await expect(page.getByTestId('query-input')).toHaveValue('>mine my own query\nACGTACGTACGTTTGACCATGG\n');
    await expect(page.getByTestId('example-loaded')).toHaveCount(0);

    // A file of the subject alone is an input too.
    await page.getByTestId('query-clear').click();
    await page.getByTestId('subject-files').setInputFiles({ name: 'mine.fa', mimeType: 'text/plain', buffer: Buffer.from('>s\nACGT\n') });
    await settled(page, 'subject');
    await page.getByTestId('try-example').click();
    await expect(question).toBeVisible();
    await page.getByTestId('try-example-confirm').click();
    await expect(page.getByTestId('example-loaded')).toBeVisible();
    await settled(page, 'subject');
    await expect(page.getByTestId('subject-source-0')).toContainText('example-blastn-subject.fasta');
    await expect(page.getByTestId('subject-source-1')).toHaveCount(0);
    await expect(page.getByTestId('threads')).toHaveValue('auto');
    await expect(page.getByTestId('job-title')).toHaveValue('');
    // The user's own edit ends the line that says the example is loaded.
    await page.getByTestId('subject-source-0-remove').click();
    await expect(page.getByTestId('example-loaded')).toHaveCount(0);
  });
});

test.describe('Getting started', () => {
  test('is open on the very first visit, closed after the user closes it, and on later visits', async ({ page }) => {
    await page.goto('/');
    const panel = page.getByTestId('getting-started');
    const toggle = page.getByTestId('getting-started-toggle');
    await expect(panel).toBeVisible();
    await expect(toggle).toHaveAttribute('aria-expanded', 'true');
    await expect(panel.locator('ol > li')).toHaveCount(6);
    await expect(panel).toContainText('Choose Run LOSAT.');
    await expect(panel).toContainText('Try an example');
    await expect(panel).toContainText('stay in this browser');
    await page.getByTestId('getting-started-close').click();
    await expect(panel).toHaveCount(0);
    await expect(toggle).toHaveAttribute('aria-expanded', 'false');
    expect(await page.evaluate((key) => localStorage.getItem(key), GETTING_STARTED_KEY)).toBe('closed');
    await page.reload();
    await expect(page.getByTestId('storage-status')).toBeVisible();
    await expect(panel).toHaveCount(0);
    // The header opens it again; that does not change what later visits start with.
    await toggle.click();
    await expect(panel).toBeVisible();
    await page.reload();
    await expect(page.getByTestId('storage-status')).toBeVisible();
    await expect(panel).toHaveCount(0);
  });

  test('works with a localStorage that throws: open on each visit, closed in the tab after a close', async ({ context, page }) => {
    await context.addInitScript(() => {
      Object.defineProperty(window, 'localStorage', {
        configurable: true,
        get() {
          throw new DOMException('The operation is insecure.', 'SecurityError');
        },
      });
    });
    await page.goto('/');
    const panel = page.getByTestId('getting-started');
    await expect(panel).toBeVisible();
    await page.getByTestId('getting-started-close').click();
    await expect(panel).toHaveCount(0);
    // The page keeps working without storage.
    await paste(page, 'query', '>q\nACGTACGTACGT\n');
    await expect(page.getByTestId('query-source-0-summary')).toContainText('1 record');
    await page.reload();
    await expect(panel).toBeVisible();
  });
});

test.describe('self-test in the FakeEngine build', () => {
  test.skip(BUILD_HAS_ENGINE, 'the engine build runs the self-test');

  test('says that it needs the engine, and cannot be started', async ({ page }) => {
    await page.goto('/');
    await page.getByTestId('self-test').locator('summary').click();
    await expect(page.getByTestId('self-test-unavailable')).toHaveText(
      'The self-test needs the engine. This development build uses a fake engine, whose outputs are not search results.',
    );
    await expect(page.getByTestId('self-test-run')).toBeDisabled();
    await expect(page.getByTestId('self-test')).toContainText('BLASTX is not tested');
  });
});

// --- engine build --------------------------------------------------------------------------------

test.describe('engine build', () => {
  test.skip(!BUILD_HAS_ENGINE, 'real searches need the engine (LOSAT_WEB_REACTORS)');

  for (const id of RUNNABLE) {
    test(`${LABELS[id]}: "Try an example", then "Run LOSAT", gives results with hits`, async ({ context, page }) => {
      await closedGettingStarted(context);
      await page.goto('/');
      await tryExample(page, id);
      await submit(page);
      await waitStatus(page, 1, 'completed', 60_000);
      const { command, output } = await result(page, 1, 6);
      expect(command).toContain(`LOSAT ${id} -query example-${id}-query.fasta -subject example-${id}-subject.fasta`);
      const rows = output.split('\n').filter((line) => line.trim() !== '');
      expect(rows.length).toBeGreaterThanOrEqual(1);
      expect(rows[0]!.split('\t')[0]).toBe(recordIds(id, 'query')[0]);
    });
  }

  test('the self-test passes every case', async ({ page }) => {
    await page.goto('/');
    await page.getByTestId('self-test').locator('summary').click();
    await page.getByTestId('self-test-run').click();
    const resultBox = page.getByTestId('self-test-result');
    await expect(page.getByTestId('self-test-summary')).toHaveText(/^All 8 cases passed/, { timeout: 180_000 });
    await expect(resultBox).toHaveAttribute('data-status', 'finished');
    const cases = resultBox.locator('[data-testid^="self-test-case-"]');
    await expect(cases).toHaveCount(8);
    await expect(resultBox.locator('[data-outcome="pass"]')).toHaveCount(8);
    await expect(page.getByTestId('self-test-facts')).toContainText('Cross-origin isolated: yes.');
    // The threaded engine ran the cases of more than one thread (first), the serial engine those of one.
    await expect(resultBox.locator('[data-testid^="self-test-case-blastn-"]').first()).toContainText(/threaded engine, \d+ threads/);
    await expect(page.getByTestId('self-test-case-blastn-1')).toContainText('serial engine, 1 thread');
    // Nothing of it in the queue.
    await expect(page.getByTestId('run-1')).toHaveCount(0);
  });
});

// --- screen records ------------------------------------------------------------------------------

const SCREENS = process.env['LOSAT_WEB_SCREENS'] || undefined;

test.describe('screen records', () => {
  test.skip(SCREENS === undefined || !BUILD_HAS_ENGINE, 'screen records are taken only with LOSAT_WEB_SCREENS in the engine build');

  for (const size of [
    { name: 'desktop', width: 1280, height: 900 },
    { name: 'phone', width: 390, height: 844 },
  ] as const) {
    test(`onboarding at the ${size.name} size`, async ({ page, browserName }) => {
      const shoot = async (name: string) => {
        const directory = join(SCREENS!, browserName);
        mkdirSync(directory, { recursive: true });
        await page.waitForTimeout(300);
        await page.screenshot({ path: join(directory, `${size.name}-${name}.png`), fullPage: true });
      };
      await page.setViewportSize({ width: size.width, height: size.height });
      await page.goto('/');
      await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', /^(ok|invalid)$/, { timeout: 60_000 });
      await shoot('41-getting-started-first-visit');
      await page.getByTestId('self-test').locator('summary').click();
      await page.getByTestId('self-test-run').click();
      await expect(page.getByTestId('self-test-summary')).toHaveText(/^All 8 cases passed/, { timeout: 180_000 });
      await shoot('42-self-test-passed');
      await page.getByTestId('getting-started-close').click();
      await paste(page, 'query', '>mine\nMKVLAAGG\n');
      await program(page, 'tblastn');
      await page.getByTestId('try-example').click();
      await shoot('43-try-example-question');
      await page.getByTestId('try-example-confirm').click();
      await expect(page.getByTestId('example-loaded')).toBeVisible();
      await settled(page, 'query');
      await settled(page, 'subject');
      await shoot('44-example-loaded');
      await program(page, 'blastx');
      await shoot('45-blastx-no-example');
    });
  }
});
