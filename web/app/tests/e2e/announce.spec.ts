// W7 (line WA, stage 2): what a screen reader is told and where the keyboard's focus goes.
// WA2-AX-03: a run's status, its end (completed, failed, cancelled), an Example that loads and
// the check of an input are told by live regions that were in the page before (support/live.ts).
// WA2-AX-04: a button that leaves the page, or that cannot be pressed for a moment, does not drop
// the focus to the page's body. WA2-AX-06: n and p on the Dot Plot are told. WA2-AX-08: the
// output texts that scroll are Tab stops, also in WebKit. WA2-AX-13: the hidden file input is not
// a second button.
import { gzipSync } from 'node:zlib';
import { expect, type Locator, type Page } from '@playwright/test';
import { test } from './support/test';
import { exampleById } from '../../src/domain/examples';
import { axeViolations } from './support/axe';
import { BUILD_HAS_ENGINE } from './support/browser';
import { focused, forget, listen, told } from './support/live';
import { chooseRun, fasta, openFiles, openSearch, paste, program, settled, submit, waitStatus } from './support/search';

test.setTimeout(180_000);

/** Pseudo-random letters (a linear congruential generator): the same in every run. */
function letters(alphabet: string, seed: number, length: number): string {
  let state = seed;
  let out = '';
  for (let i = 0; i < length; i++) {
    state = (Math.imul(state, 1103515245) + 12345) >>> 0;
    out += alphabet[(state >>> 16) % alphabet.length];
  }
  return out;
}
const dna = (seed: number, length: number) => letters('ACGT', seed, length);
const PROTEIN = letters('ACDEFGHIKLMNPQRSTVWY', 7, 90);
/** Three subjects that hold the query: three subjects with HSPs, for the engine and the FakeEngine. */
const THREE_SUBJECTS = [1, 2, 3].map((n) => `>s${n} subject ${n}\n${letters('ACDEFGHIKLMNPQRSTVWY', 100 + n, 20)}${PROTEIN}\n`).join('');

async function start(page: Page): Promise<void> {
  await listen(page);
  await openSearch(page);
}

/** Presses Tab until the element with test ID `testid` has the focus; false if it never does. */
async function tabTo(page: Page, testid: string, most = 60): Promise<boolean> {
  for (let i = 0; i < most; i++) {
    await page.keyboard.press('Tab');
    if ((await focused(page)) === testid) return true;
  }
  return false;
}

/**
 * Presses `key` until the focused `box` has scrolled, as a person presses again when the text did not move.
 * In CI's Chromium (#904), an ArrowRight 17 ms after the Tab that focused the box, before the
 * focus was drawn, left it unscrolled; the machines here scroll on the first press.
 */
async function pressUntilScrolled(page: Page, key: string, box: Locator, side: 'scrollTop' | 'scrollLeft'): Promise<void> {
  await expect(async () => {
    await page.keyboard.press(key);
    expect(await box.evaluate((element, s) => element[s], side)).toBeGreaterThan(0);
  }).toPass({ timeout: 10_000 });
}

test('a run is told as it goes and when it completes, without the elapsed time; the focus stays on Results', async ({ page }) => {
  await start(page);
  await program(page, 'blastp');
  await paste(page, 'query', `>q1\n${PROTEIN}\n`);
  await paste(page, 'subject', THREE_SUBJECTS);
  await forget(page);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
  await expect.poll(() => told(page, 'status')).toContainEqual(expect.stringMatching(/^Run 1 completed: \d+ HSPs? on 3 subjects\.$/));
  const texts = await told(page);
  expect(texts).toContainEqual(expect.stringMatching(/^Status: (Preparing|Searching|Organizing results)$/));
  expect(texts.filter((text) => /\d:\d\d/.test(text))).toEqual([]);
  expect(await focused(page)).toBe('results-heading');
  // A run that the page only shows again (another tab, then back) is not told again.
  await forget(page);
  await page.getByTestId('tab-search').click();
  await page.getByTestId('tab-results').click();
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
  expect(await told(page)).toEqual([]);
});

test('a failed run is told as an alert, with its reason', async ({ page }) => {
  await start(page);
  if (BUILD_HAS_ENGINE) {
    // NCBI's "The average subject length is too short" (results.spec.ts has the same input).
    await program(page, 'tblastx');
    await paste(page, 'query', `>q1\n${dna(7, 300)}\n`);
    await paste(page, 'subject', '>s1 alpha\n');
  } else {
    // The FakeEngine fails a run of an empty query.
    await program(page, 'blastp');
    await openFiles(page, 'query', [{ name: 'empty.fa', text: '' }]);
    await paste(page, 'subject', THREE_SUBJECTS);
  }
  await forget(page);
  await submit(page);
  await waitStatus(page, 1, 'failed');
  await expect.poll(() => told(page, 'alert')).toContainEqual(expect.stringMatching(/^Run 1 failed: .+ A failed run keeps no results\.$/));
});

test('a cancelled run is told; Cancel leaves the focus on the Results heading, and Tab goes on from there', async ({ page }) => {
  await listen(page);
  if (!BUILD_HAS_ENGINE) await page.clock.install();
  await openSearch(page);
  await program(page, 'blastp');
  if (BUILD_HAS_ENGINE) {
    // A proteome against itself: many seconds (search.spec.ts).
    await openFiles(page, 'query', [{ name: 'NZ_CP006932.faa', text: fasta('NZ_CP006932.faa') }]);
    await openFiles(page, 'subject', [{ name: 'NZ_CP006932.faa', text: fasta('NZ_CP006932.faa') }]);
  } else {
    await paste(page, 'query', `>q1\n${PROTEIN}\n`);
    await paste(page, 'subject', THREE_SUBJECTS);
    // The FakeEngine's phases wait on timers: the run stays in progress until the clock goes on.
    await page.clock.pauseAt((await page.evaluate(() => Date.now())) + 1000);
  }
  await forget(page);
  await submit(page);
  const cancel = page.getByTestId('run-status-cancel');
  await expect(cancel).toBeVisible();
  await cancel.focus();
  await page.keyboard.press('Enter');
  if (!BUILD_HAS_ENGINE) await page.clock.resume();
  await waitStatus(page, 1, 'cancelled');
  await expect(page.getByTestId('results-status')).toHaveAttribute('data-run-status', 'cancelled');
  await expect.poll(() => told(page, 'status')).toContainEqual('Run 1 was cancelled. A cancelled run keeps no results.');
  expect((await told(page)).filter((text) => /\d:\d\d/.test(text))).toEqual([]);
  expect(await focused(page)).toBe('results-heading');
  await page.keyboard.press('Tab');
  expect(await focused(page)).toBe('results-run');
});

test('a run that waits for another run says so when its page opens', async ({ page }) => {
  test.skip(BUILD_HAS_ENGINE, "the FakeEngine's runs wait on the page's timers, which the test holds");
  await listen(page);
  await page.clock.install();
  await openSearch(page);
  await program(page, 'blastp');
  await paste(page, 'query', `>q1\n${PROTEIN}\n`);
  await paste(page, 'subject', THREE_SUBJECTS);
  // Run 1 stays in progress until the clock goes on; run 2 waits for it.
  await page.clock.pauseAt((await page.evaluate(() => Date.now())) + 1000);
  await submit(page);
  await expect(page.getByTestId('run-status')).toHaveAttribute('data-run', '1');
  await page.getByTestId('tab-search').click();
  await forget(page);
  await submit(page);
  await expect(page.getByTestId('run-status')).toHaveAttribute('data-run', '2');
  await expect(page.getByTestId('run-status-phase')).toHaveText('Waiting for run 1');
  await expect.poll(() => told(page, 'status')).toContainEqual('Status: Waiting for run 1');
  await page.clock.resume();
  await waitStatus(page, 1, 'completed');
  await waitStatus(page, 2, 'completed');
});

test('"Run LOSAT" that cannot start says why in an alert, and keeps the focus', async ({ page }) => {
  await start(page);
  await program(page, 'blastn');
  await forget(page);
  await page.getByTestId('run-losat').focus();
  await page.keyboard.press('Enter');
  await expect(page.getByTestId('search-message')).toHaveText('Add the query sequences: paste them or open FASTA files.');
  await expect.poll(() => told(page, 'alert')).toContainEqual('Add the query sequences: paste them or open FASTA files.');
  expect(await focused(page)).toBe('run-losat');
});

test("the engine's check of an input, its warnings and a refusal to read it are told", async ({ page }) => {
  await start(page);
  await program(page, 'blastn');
  await forget(page);
  await paste(page, 'query', `>q1\n${dna(3, 60)}\n`);
  await expect.poll(() => told(page, 'status')).toContainEqual(expect.stringContaining('Accepted by the BLASTN engine: 1 record'));
  // Another check that finds the same (the text edited, then another nucleotide mode) is not told again.
  await page.getByTestId('query-input').fill(`>q1\n${dna(3, 90)}\n`);
  await expect(page.getByTestId('query-source-0-summary')).toContainText('90 nt');
  await settled(page, 'query');
  await program(page, 'tblastx');
  await expect(page.getByTestId('query-source-0-check')).toContainText('Accepted by the TBLASTX engine: 1 record');
  expect((await told(page, 'status')).filter((text) => text.includes('Accepted by the'))).toHaveLength(1);
  await program(page, 'blastn');
  await paste(page, 'query', `>p1\n${PROTEIN}\n`);
  await expect
    .poll(() => told(page, 'status'))
    .toContainEqual(expect.stringContaining('1 record looks like protein sequences, but BLASTN reads the query as nucleotide.'));
  await page.getByTestId('query-input').fill(`\uFEFF>q1\n${dna(3, 60)}\n`);
  await expect(page.getByTestId('query-source-0')).toHaveAttribute('data-status', 'failed');
  await expect.poll(() => told(page, 'status')).toContainEqual(expect.stringContaining('This input cannot be read: '));
});

const SIXTEEN_S = exampleById('refseq-16s-bacteria');

/** The 16S set's host, stubbed at its real address (examples.spec.ts). */
async function stubSixteenS(page: Page): Promise<void> {
  const gz = gzipSync('>NR_102783.2 Bacillus subtilis 16S\nACGTACGTACGTACGTACGT\n>NR_114042.1 Escherichia coli 16S\nACGTTCGTACGTACGTTCGT\n');
  await page.route(/^https:\/\/ftp\.ncbi\.nlm\.nih\.gov\//, (route) =>
    route.fulfill({
      status: 200,
      contentType: 'application/x-gzip',
      headers: { 'Access-Control-Allow-Origin': '*', 'Content-Length': String(gz.length), 'Last-Modified': 'Tue, 27 Oct 2026 12:04:14 GMT' },
      body: gz,
    }),
  );
}

test('an Example that loads is told, and its button keeps the focus', async ({ page }) => {
  await stubSixteenS(page);
  await start(page);
  await forget(page);
  const button = page.getByTestId(`example-${SIXTEEN_S.id}`);
  await button.focus();
  await page.keyboard.press('Enter');
  await expect(page.getByTestId('subject-source-0')).toContainText(SIXTEEN_S.fileName);
  await settled(page, 'subject');
  await expect.poll(() => told(page, 'status')).toContainEqual(expect.stringContaining(`Example ${SIXTEEN_S.label} loaded.`));
  expect(await focused(page)).toBe(`example-${SIXTEEN_S.id}`);
});

test("Delete of an Example's saved copy leaves with the copy; the focus goes to the set's button", async ({ page }) => {
  await stubSixteenS(page);
  await start(page);
  const backend = await page.getByTestId('storage-status').getAttribute('data-backend');
  test.skip(backend !== 'opfs', `Delete is offered for a copy saved in this browser (OPFS); this browser keeps the files in ${backend}`);
  await page.getByTestId(`example-${SIXTEEN_S.id}`).click();
  await settled(page, 'subject');
  const remove = page.getByTestId(`example-${SIXTEEN_S.id}-delete`);
  await remove.focus();
  await page.keyboard.press('Enter');
  await expect(remove).toHaveCount(0);
  expect(await focused(page)).toBe(`example-${SIXTEEN_S.id}`);
});

test("an Example's Cancel that leaves as the load ends gives the focus back to the set's button", async ({ page }) => {
  const gz = gzipSync('>NR_102783.2 Bacillus subtilis 16S\nACGTACGTACGTACGTACGT\n');
  let release: () => void = () => undefined;
  const held = new Promise<void>((resolve) => (release = resolve));
  await page.route(/^https:\/\/ftp\.ncbi\.nlm\.nih\.gov\//, async (route) => {
    await held;
    await route.fulfill({ status: 200, contentType: 'application/x-gzip', headers: { 'Access-Control-Allow-Origin': '*' }, body: gz });
  });
  await start(page);
  await page.getByTestId(`example-${SIXTEEN_S.id}`).click();
  const cancel = page.getByTestId(`example-${SIXTEEN_S.id}-cancel`);
  await cancel.focus();
  release();
  await expect(cancel).toHaveCount(0);
  await expect(page.getByTestId('subject-source-0')).toContainText(SIXTEEN_S.fileName);
  expect(await focused(page)).toBe(`example-${SIXTEEN_S.id}`);
});

test('"Run a search" of the empty Results tab moves the focus to the Search tab', async ({ page }) => {
  await start(page);
  await page.getByTestId('tab-results').click();
  await page.getByTestId('results-run-a-search').focus();
  await page.keyboard.press('Enter');
  await expect(page.getByTestId('search-panel')).toBeVisible();
  expect(await focused(page)).toBe('tab-search');
});

test('Remove moves the focus to the next input, then to "Open FASTA files…"', async ({ page }) => {
  await start(page);
  await openFiles(page, 'subject', [
    { name: 'a.fa', text: `>a\n${dna(4, 40)}\n` },
    { name: 'b.fa', text: `>b\n${dna(5, 40)}\n` },
  ]);
  await page.getByTestId('subject-source-0-remove').focus();
  await page.keyboard.press('Enter');
  await expect(page.getByTestId('subject-source-1')).toHaveCount(0);
  expect(await focused(page)).toBe('subject-source-0-remove');
  await expect(page.getByTestId('subject-source-0')).toContainText('b.fa');
  await page.keyboard.press('Enter');
  await expect(page.getByTestId('subject-source-0')).toHaveCount(0);
  expect(await focused(page)).toBe('subject-open');
});

test('Dismiss of the notice after the tab was hidden moves the focus to the tab shown', async ({ page }) => {
  await start(page);
  await program(page, 'blastp');
  await paste(page, 'query', `>q1\n${PROTEIN}\n`);
  await paste(page, 'subject', THREE_SUBJECTS);
  // Hidden as soon as a run is in progress (search.spec.ts), shown again after the report's 1 s.
  await page.evaluate(() => {
    const busy = ['queued', 'preparing', 'running', 'finalizing'].map((status) => `[data-testid^="run-"][data-status="${status}"]`).join(', ');
    const observer = new MutationObserver(() => {
      if (document.querySelector(busy) === null) return;
      observer.disconnect();
      Object.defineProperty(document, 'visibilityState', { get: () => 'hidden', configurable: true });
      document.dispatchEvent(new Event('visibilitychange'));
    });
    observer.observe(document.body, { subtree: true, childList: true, attributes: true, attributeFilter: ['data-status'] });
  });
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await expect.poll(() => page.evaluate(() => document.visibilityState)).toBe('hidden');
  await page.waitForTimeout(1100);
  await page.evaluate(() => {
    Object.defineProperty(document, 'visibilityState', { get: () => 'visible', configurable: true });
    document.dispatchEvent(new Event('visibilitychange'));
  });
  const dismiss = page.getByTestId('resume-dismiss');
  await dismiss.focus();
  await page.keyboard.press('Enter');
  await expect(page.getByTestId('resume-notice')).toHaveCount(0);
  expect(await focused(page)).toBe('tab-results');
});

test('Next on the last subject of the Alignments keeps the focus', async ({ page }) => {
  await start(page);
  await program(page, 'blastp');
  await paste(page, 'query', `>q1\n${PROTEIN}\n`);
  await paste(page, 'subject', THREE_SUBJECTS);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  const next = page.getByTestId('alignments-next-subject');
  await next.focus();
  await page.keyboard.press('Enter');
  await expect(page.getByTestId('alignments-subject')).toContainText('Subject 2 of 3');
  await page.keyboard.press('Enter');
  await expect(page.getByTestId('alignments-subject')).toContainText('Subject 3 of 3');
  await expect(next).toBeDisabled();
  expect(await focused(page)).toBe('alignments-next-subject');
  // A press of a button that cannot be pressed does nothing.
  await page.keyboard.press('Enter');
  await expect(page.getByTestId('alignments-subject')).toContainText('Subject 3 of 3');
});

test('n and p on the Dot Plot, with its popup closed, tell the HSP selected', async ({ page }) => {
  await start(page);
  // q2 is A, which the subject holds twice: two HSPs (the FakeEngine also writes two HSPs for its second query).
  const a = dna(11, 60);
  await program(page, 'blastn');
  await paste(page, 'query', `>q1\n${dna(21, 20)}\n>q2\n${a}\n`);
  await paste(page, 'subject', `>s1\n${a}${dna(12, 60)}${a}\n`);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await page.getByTestId('query-row-1').click();
  await page.getByTestId('pane-dotplot').click();
  const canvas = page.getByTestId('dotplot-canvas');
  await expect(canvas).toHaveAttribute('data-segments', '2');
  await canvas.focus();
  await forget(page);
  const selectedText = async () => ((await page.getByTestId('dotplot-selected').textContent()) ?? '').replace(/\s+/g, ' ').trim();
  for (const key of ['n', 'p']) {
    const before = (await canvas.getAttribute('data-selected')) ?? '';
    await page.keyboard.press(key);
    await expect(canvas).not.toHaveAttribute('data-selected', before);
    await expect(page.getByTestId('dotplot-popup')).toHaveCount(0);
    const selected = await selectedText();
    expect(selected).toMatch(/^Selected: HSP \d+, query /);
    await expect.poll(async () => (await told(page, 'status')).at(-1)).toBe(selected);
  }
});

test('the output text of the Outputs view is a Tab stop and scrolls with the keyboard', async ({ page }) => {
  await start(page);
  await program(page, 'blastp');
  const queries = Array.from({ length: 60 }, (_, i) => `>q${i + 1}\n${letters('ACDEFGHIKLMNPQRSTVWY', 300 + i, 40)}${PROTEIN}\n`).join('');
  await paste(page, 'query', queries);
  await paste(page, 'subject', THREE_SUBJECTS);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  await page.getByTestId('results-view-outputs').click();
  await page.getByTestId('format-0').click();
  const output = page.getByTestId('result-output');
  await expect(output).toHaveAttribute('data-shown', '1:0');
  await expect(page.getByRole('region', { name: 'outfmt 0 output' })).toBeVisible();
  await expect(output).toHaveAttribute('tabindex', '0');
  await page.getByTestId('format-7').focus();
  expect(await tabTo(page, 'result-output', 1)).toBe(true);
  await pressUntilScrolled(page, 'PageDown', output, 'scrollTop');
  expect(await focused(page)).toBe('result-output');
});

test.describe('on a phone', () => {
  test.use({ viewport: { width: 390, height: 844 } });

  test("an HSP's alignment text that scrolls sideways is a Tab stop and scrolls with the keyboard", async ({ page }) => {
    await start(page);
    await program(page, 'blastp');
    await paste(page, 'query', `>q1\n${PROTEIN}\n`);
    await paste(page, 'subject', THREE_SUBJECTS);
    await submit(page);
    await waitStatus(page, 1, 'completed');
    await chooseRun(page, 1);
    const section = page.getByTestId('detail-section');
    await expect(section).toBeVisible();
    expect(await section.evaluate((element) => element.scrollWidth > element.clientWidth)).toBe(true);
    await expect(section).toHaveAttribute('tabindex', '0');
    await page.getByTestId('results-heading').focus();
    expect(await tabTo(page, 'detail-section')).toBe(true);
    await pressUntilScrolled(page, 'ArrowRight', section, 'scrollLeft');
    expect(await focused(page)).toBe('detail-section');
  });
});

test('the hidden file input is not a second button beside "Open FASTA files…"; a file still opens', async ({ page }) => {
  await start(page);
  await expect(page.getByRole('button', { name: 'Open query FASTA files' })).toHaveCount(0);
  const tree = await page.getByTestId('query-dropzone').ariaSnapshot();
  expect(tree.split('\n').filter((line) => line.includes('button'))).toEqual(['- button "Open FASTA files…"']);
  await openFiles(page, 'query', [{ name: 'q.fa', text: `>q\n${dna(6, 40)}\n` }]);
  await expect(page.getByTestId('query-source-0')).toContainText('q.fa');
});

test('"Show all" of the Graphic Summary moves the focus to the first subject that it adds', async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, 'the FakeEngine writes HSPs on at most three subjects; "Show all" needs more than 100');
  await start(page);
  const a = dna(11, 60);
  await program(page, 'blastn');
  await paste(page, 'query', `>q1\n${a}\n`);
  await paste(page, 'subject', Array.from({ length: 120 }, (_, i) => `>s${i + 1}\n${dna(500 + i, 40)}${a}${dna(900 + i, 40)}\n`).join(''));
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  const canvas = page.getByTestId('graphic-canvas');
  await expect(canvas).toHaveAttribute('data-rows', '100');
  // The keyboard reaches the button, as a keyboard user's does. Firefox gives the plot :focus-visible
  // (and so tells its subject) only when the button's focus came from a key, not from focus().
  await canvas.focus();
  await page.keyboard.press('Tab');
  expect(await focused(page)).toBe('graphic-show-all');
  await forget(page);
  await page.keyboard.press('Enter');
  await expect(canvas).toHaveAttribute('data-rows', '120');
  expect(await focused(page)).toBe('graphic-canvas');
  await expect.poll(() => told(page, 'status')).toContainEqual(expect.stringMatching(/^Subject 101 of 120: /));
});

/** The development build's banner (the FakeEngine) is outside the page's landmarks; the engine build has none. */
const NOT_THE_APP = ['[data-testid="fake-engine-banner"]'];

/** axe's violations in the page, once the colours of a tab just chosen have ended their short transition. */
async function violations(page: Page) {
  await page.evaluate(() =>
    Promise.all(document.getAnimations().filter((animation) => animation instanceof CSSTransition).map((animation) => animation.finished)),
  );
  return axeViolations(page, ['body'], NOT_THE_APP);
}

test('axe finds nothing on the pages that these changes touch', async ({ page }) => {
  await start(page);
  await program(page, 'blastp');
  await paste(page, 'query', `>q1\n${PROTEIN}\n`);
  await openFiles(page, 'subject', [
    { name: 'a.fa', text: THREE_SUBJECTS },
    { name: 'b.fa', text: `>b\n${dna(5, 40)}\n` },
  ]);
  expect(await violations(page)).toEqual([]);
  const queries = Array.from({ length: 60 }, (_, i) => `>q${i + 1}\n${letters('ACDEFGHIKLMNPQRSTVWY', 300 + i, 40)}${PROTEIN}\n`).join('');
  await paste(page, 'query', queries);
  await submit(page);
  await waitStatus(page, 1, 'completed');
  await chooseRun(page, 1);
  expect(await violations(page)).toEqual([]);
  for (const [tab, shown] of [
    ['pane-dotplot', 'dotplot-canvas'],
    ['results-view-details', 'run-record'],
    ['results-view-outputs', 'result-output'],
  ] as const) {
    await page.getByTestId(tab).click();
    await expect(page.getByTestId(shown)).toBeVisible();
    if (tab === 'results-view-outputs') {
      await page.getByTestId('format-0').click();
      await expect(page.getByTestId('result-output')).toHaveAttribute('data-shown', '1:0');
    }
    expect(await violations(page), tab).toEqual([]);
  }
  // A phone's width: the alignment texts scroll sideways.
  await page.setViewportSize({ width: 390, height: 844 });
  await page.getByTestId('results-view-hits').click();
  await expect(page.getByTestId('detail-section')).toBeVisible();
  expect(await violations(page)).toEqual([]);
});
