// The steps of the random walk (walk/plan.ts names their kinds) and what every step must leave
// true. Each step acts as a user does, through the page's controls, and draws its arguments
// (which field, which value, which record) as numbers in [0, 1) that it records, so a replay
// repeats them: a number picks among what the page shows at that moment, never by an ID that
// changes from one page load to the next (the lesson of gbdraw's walk, whose replays keyed steps
// by per-load IDs).
//
// After every step: no page error, console error, dialog, CSP refusal or request to another origin
// (web/AGENTS.md rule 4: research data never leaves the page), no script of the inputs ran
// (window.__walkXss), the page is no wider than its window, axe finds nothing, the form shows the
// mode and the inputs that the walk gave.
// After every press of Run LOSAT: the runs it starts are as many as the inputs the walk gave make
// (one, or one per input of a role searched apart), or it starts none and says why, and a reason
// that the walk can judge (no input, every record left out) agrees with the inputs it gave.
// After every run: it ends, and not as failed (a search that the form let through must run); a
// completed run's stored outfmt 0, 6 and 7 (and its own -outfmt) are the bytes that its command
// line, as the run's details show it, writes when the native CLI of the same commit runs it on the
// inputs the walk gave (LOSAT_WEB_NATIVE); the options of its argv are the fields marked changed,
// Additional options and the subranges, and a changed text box's or select's value is its option's
// value; two runs of the same form give the same bytes. A walk of the engine's build completes
// runs, and compares them with the native CLI when it is given (runner.ts).
import { spawnSync } from 'node:child_process';
import { mkdirSync, readFileSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { axeViolations } from '../support/axe';
import { NATIVE } from '../support/native';
import { navigate, reload } from '../support/navigate';
import { described, storedOutputs } from '../support/search';
import { alphabetOf, pickFor, type FixtureLibrary } from './fixtures';
import {
  emptyModel,
  expectedInputs,
  parseShellWords,
  withFiles,
  withIncluded,
  withMode,
  withoutSource,
  withPaste,
  type ExpectedInput,
  type InputModel,
  type Role,
} from './model';
import type { Random, WalkKind } from './plan';

export const PROGRAMS = Object.freeze(['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'] as const);
const ROLES: readonly Role[] = ['query', 'subject'];
const WIDTHS = Object.freeze([320, 390, 768, 1024, 1280, 1920]);
const FORMATS = Object.freeze([0, 6, 7] as const);
const RUN_TIMEOUT = 120_000;
const SETTLE_TIMEOUT = 30_000;

/** A step as a walk ran it: its kind, the draws it made, and what it did. */
export interface RecordedStep {
  readonly kind: WalkKind;
  readonly args: readonly number[];
  readonly note?: string;
}

/** A step's failure: `signature` names the kind of defect, so a replay reproduces the same one. */
export class StepFailure extends Error {
  constructor(
    readonly signature: string,
    message: string,
  ) {
    super(message);
  }
}

/** The draws of a step: those recorded (a replay), then new ones from the walk's generator. */
export class Draws {
  readonly taken: number[] = [];
  constructor(
    private readonly recorded: readonly number[] | undefined,
    private readonly random: Random,
  ) {}
  next(): number {
    const index = this.taken.length;
    const value = this.recorded !== undefined && index < this.recorded.length ? this.recorded[index]! : this.random.next();
    this.taken.push(value);
    return value;
  }
  /** An index in [0, count). */
  index(count: number): number {
    if (count <= 0) throw new Error('no item to pick');
    return Math.min(count - 1, Math.floor(this.next() * count));
  }
  pick<T>(items: readonly T[]): T {
    return items[this.index(items.length)]!;
  }
}

/** What the page has met since it opened, and the walk's view of it. */
export interface Walker {
  readonly page: Page;
  readonly fixtures: FixtureLibrary;
  /** Where the native CLI's files go. */
  readonly outDir: string;
  model: InputModel;
  program: string;
  /** The runs that this page started. */
  runs: number;
  /** Problems the page's events reported since the last step's check. */
  readonly problems: string[];
  notes: string[];
  /** The outputs compared with the native CLI, the runs that ended each way, and the completed runs whose inputs the walk could not tell apart. */
  readonly counts: { compared: number; completed: number; failed: number; cancelled: number; refused: number; unmatched: number; offered: number; offeredLeftOut: number };
}

/** WebKit's report of a <select>'s own inline style (see `watch`). */
const WEBKIT_SELECT_CSP =
  "Refused to apply a stylesheet because its hash, its nonce, or 'unsafe-inline' appears in neither the style-src directive nor the default-src directive of the Content Security Policy.";

/** Starts watching a page: call before it opens. */
export async function watch(page: Page, fixtures: FixtureLibrary, outDir: string): Promise<Walker> {
  const walker: Walker = {
    page,
    fixtures,
    outDir,
    model: emptyModel(),
    program: 'blastn',
    runs: 0,
    problems: [],
    notes: [],
    counts: { compared: 0, completed: 0, failed: 0, cancelled: 0, refused: 0, unmatched: 0, offered: 0, offeredLeftOut: 0 },
  };
  const engine = page.context().browser()?.browserType().name();
  page.on('pageerror', (error) => walker.problems.push(`page-error: ${error.message}`));
  page.on('console', (message) => {
    if (message.type() !== 'error') return;
    // Playwright's WebKit (WPE, GTK) reports this for every <select> that a script makes: the select's
    // own shadow tree has inline style that the page's CSP (no 'unsafe-inline') refuses. A page of
    // the site's CSP with only `document.body.append(document.createElement('select'))` gives it,
    // Chromium and Firefox do not (losat-baselines/web-qa-20261119/probe/csp8.mjs, 2026-11-19).
    // The walk's Chromium and Firefox runs still catch a style that the application itself adds,
    // and the violations below are checked in every engine.
    if (engine === 'webkit' && message.text() === WEBKIT_SELECT_CSP) return;
    walker.problems.push(`console-error: ${message.text()}`);
  });
  await page.addInitScript(() => {
    const seen: string[] = [];
    (window as unknown as { __walkCsp: string[] }).__walkCsp = seen;
    document.addEventListener('securitypolicyviolation', (event) =>
      seen.push(`${event.effectiveDirective} ${event.blockedURI} ${event.sourceFile}:${event.lineNumber}:${event.columnNumber}`),
    );
  });
  page.on('crash', () => walker.problems.push('crash: the page crashed'));
  page.on('dialog', (dialog) => {
    if (dialog.type() === 'beforeunload') {
      void dialog.accept();
      return;
    }
    walker.problems.push(`dialog: ${dialog.type()} ${dialog.message()}`);
    void dialog.dismiss();
  });
  // The context's requests: those of the page's workers (Data, Engine) as well as the page's own.
  page.context().on('request', (request) => {
    const url = new URL(request.url());
    if (url.protocol !== 'http:' && url.protocol !== 'https:') return;
    const own = new URL(page.url() === 'about:blank' ? request.url() : page.url()).origin;
    if (url.origin !== own) walker.problems.push(`other-origin: a request to ${request.url()}`);
  });
  return walker;
}

/**
 * Opens the application, reads its mode, and pastes a query and a subject of that mode, so that a
 * walk's first runs search (its steps change them).
 */
export async function open(walker: Walker): Promise<void> {
  await navigate(walker.page, '/');
  await start(walker);
}

/** The start of a walk and of the page after a reload: the mode's query and subject pasted. */
async function start(walker: Walker): Promise<void> {
  const { page, fixtures } = walker;
  await described(page);
  walker.model = emptyModel();
  walker.program = await checkedProgram(page);
  walker.runs = 0;
  for (const role of ROLES) {
    const text = pickFor(fixtures.texts, [], alphabetOf(walker.program, role), 0, role === 'query' ? 0 : 0.99).text;
    await page.getByTestId(`${role}-input`).fill(text);
    walker.model = withPaste(walker.model, role, text);
    await expect.poll(async () => (await sourceCards(page, role)).length, { timeout: SETTLE_TIMEOUT }).toBe(1);
  }
  await settle(page);
}

async function checkedProgram(page: Page): Promise<string> {
  for (const id of PROGRAMS) if (await page.getByTestId(`program-${id}`).isChecked()) return id;
  throw new StepFailure('model', 'no mode is chosen');
}

// --- waits ----------------------------------------------------------------------------------

async function tabShown(page: Page, tab: 'search' | 'results'): Promise<boolean> {
  return (await page.getByTestId(`tab-${tab}`).getAttribute('aria-selected')) === 'true';
}

async function showTab(page: Page, tab: 'search' | 'results'): Promise<void> {
  if (!(await tabShown(page, tab))) await page.getByTestId(`tab-${tab}`).click();
  await expect(page.getByTestId(`tab-${tab}`)).toHaveAttribute('aria-selected', 'true');
}

/** The form's state once nothing in it is still being read or checked. */
async function settle(page: Page): Promise<void> {
  try {
    await expect
      .poll(
        () =>
          page.evaluate(() => {
            const busy: string[] = [];
            const panel = document.querySelector('[data-testid="search-panel"]');
            const describe = panel?.getAttribute('data-describe');
            if (describe === 'loading') busy.push('describe');
            for (const element of document.querySelectorAll('[data-testid][data-status]')) {
              const id = element.getAttribute('data-testid') ?? '';
              if (/^(query|subject)-source-\d+$/.test(id) && element.getAttribute('data-status') === 'indexing') busy.push(id);
            }
            for (const element of document.querySelectorAll('[data-check="pending"]')) busy.push(element.getAttribute('data-testid') ?? 'check');
            if (document.querySelector('[data-testid="argv-validation"]')?.getAttribute('data-state') === 'checking') busy.push('validation');
            return busy.join(' ');
          }),
        { timeout: SETTLE_TIMEOUT },
      )
      .toBe('');
  } catch {
    throw new StepFailure('stuck', 'the form was still reading or checking after 30 s');
  }
}

/** The sources of a role that the page shows, in order, with their status. */
async function sourceCards(page: Page, role: Role): Promise<Array<{ status: string; check: string | null }>> {
  return page.evaluate((role) => {
    const cards: Array<{ index: number; status: string; check: string | null }> = [];
    for (const element of document.querySelectorAll(`[data-testid^="${role}-source-"][data-status]`)) {
      const match = new RegExp(`^${role}-source-(\\d+)$`).exec(element.getAttribute('data-testid') ?? '');
      if (match === null) continue;
      const check = element.querySelector(`[data-testid="${role}-source-${match[1]}-check"]`);
      cards.push({ index: Number(match[1]), status: element.getAttribute('data-status') ?? '', check: check?.getAttribute('data-check') ?? null });
    }
    return cards.sort((a, b) => a.index - b.index).map(({ status, check }) => ({ status, check }));
  }, role);
}

async function openParameters(page: Page): Promise<void> {
  const toggle = page.getByTestId('algorithm-parameters-toggle');
  if ((await toggle.getAttribute('aria-expanded')) !== 'true') await toggle.click();
  await expect(page.getByTestId('parameter-form')).toBeVisible();
}

// --- the checks after every step ------------------------------------------------------------

/** What must hold after every step; throws a StepFailure naming the first broken rule. */
export async function check(walker: Walker): Promise<void> {
  const { page } = walker;
  const problems = walker.problems.splice(0);
  if (problems.length > 0) throw new StepFailure(problems[0]!.split(':', 1)[0]!, problems.join('\n'));
  // The page's CSP refused nothing, but WebKit's own style of a <select> (see `watch`).
  const engine = page.context().browser()?.browserType().name();
  const refused = (await page.evaluate(() => (window as unknown as { __walkCsp?: string[] }).__walkCsp?.splice(0) ?? [])).filter(
    (line) => !(engine === 'webkit' && line.startsWith('style-src-attr inline ')),
  );
  if (refused.length > 0) throw new StepFailure('csp', `the page's CSP refused: ${refused.join('; ')}`);
  const xss = await page.evaluate(() => (window as unknown as { __walkXss?: number }).__walkXss);
  if (xss !== undefined) throw new StepFailure('xss', `a script of the inputs ran (${xss} times)`);
  const overflow = await page.evaluate(() => ({ scroll: document.documentElement.scrollWidth, client: document.documentElement.clientWidth }));
  if (overflow.scroll > overflow.client + 1) {
    throw new StepFailure('overflow', `the page is ${overflow.scroll} px wide in a window of ${overflow.client} px`);
  }
  // axe's rules hold in the page as the step left it, once the colours of a tab just chosen have
  // ended their transition (as announce.spec.ts checks the pages it touches). The FakeEngine's
  // banner is outside the page's landmarks.
  // A transition that another change cancels never finishes (its promise rejects): 2 s at most.
  await page.evaluate(() =>
    Promise.race([
      Promise.allSettled(document.getAnimations().filter((animation) => animation instanceof CSSTransition).map((animation) => animation.finished)),
      new Promise((resolve) => setTimeout(resolve, 2000)),
    ]),
  );
  const violations = await axeViolations(page, ['body'], ['[data-testid="fake-engine-banner"]']);
  if (violations.length > 0) {
    throw new StepFailure(
      `axe:${violations.map((violation) => violation.id).join(',')}`,
      `axe: ${violations.map((violation) => `${violation.id} (${violation.impact}) at ${violation.targets.slice(0, 3).join(', ')}`).join('; ')}`,
    );
  }
  // The form shows the walk's mode and inputs (the search panel stays mounted under the results).
  const program = await checkedProgram(page);
  if (program !== walker.program) throw new StepFailure('model', `the form shows mode ${program}; the walk chose ${walker.program}`);
  for (const role of ROLES) {
    const cards = await sourceCards(page, role);
    const want = walker.model[role].sources.length;
    if (cards.length !== want) throw new StepFailure('model', `the ${role} has ${cards.length} inputs; the walk gave ${want}`);
  }
}

// --- the steps ------------------------------------------------------------------------------

type StepBody = (walker: Walker, draws: Draws) => Promise<string>;

const mode: StepBody = async (walker, draws) => {
  const { page } = walker;
  await showTab(page, 'search');
  const id = draws.pick(PROGRAMS);
  await page.getByTestId(`program-${id}`).check();
  await described(page);
  walker.program = id;
  await settle(page);
  return id;
};

const task: StepBody = async (walker, draws) => {
  const { page } = walker;
  await showTab(page, 'search');
  await settle(page);
  const radios = page.locator('[data-testid="program-selection"] input[type="radio"]');
  const count = await radios.count();
  if (count === 0) return 'no task';
  const radio = radios.nth(draws.index(count));
  const id = (await radio.getAttribute('data-testid')) ?? '';
  await radio.check();
  await settle(page);
  return id;
};

/** Values typed into a field: valid ones mostly, and some that the engine refuses. */
const VALUES: Readonly<Record<string, readonly string[]>> = Object.freeze({
  numeric: ['0', '1', '2', '3', '5', '7', '11', '28', '100', '-1', 'abc', ''],
  decimal: ['10', '1', '0.001', '1e-5', '1e-30', '0', '2.5', '-2', 'abc', ''],
  text: ['', 'yes', 'no', '1', 'abc', '20 64 1', '12 2.2 2.5', 'm L; R -d rodents.lib'],
});

function valueFor(inputMode: string | null, draws: Draws): string {
  return draws.pick(VALUES[inputMode ?? 'text'] ?? VALUES['text']!);
}

const field: StepBody = async (walker, draws) => {
  const { page } = walker;
  await showTab(page, 'search');
  await settle(page);
  await openParameters(page);
  // Every field of the form that the page shows, but the task (its own step) and the files of
  // WindowMasker (which ask for a download). A row of two fields is one by its select; without one
  // (its pair is outside NCBI's lists: ParameterPair) its two fields are the fields.
  const ids = await page.evaluate(() =>
    [...document.querySelectorAll('[data-testid="search-panel"] [data-testid$="-field"]')]
      .filter((element) => (element as HTMLElement).offsetParent !== null)
      .filter((element) => element.getAttribute('data-kind') !== 'file' && element.getAttribute('data-testid') !== 'param-task-field')
      .map((element) => ({ id: element.getAttribute('data-testid')!.replace(/-field$/, ''), kind: element.getAttribute('data-kind') ?? '' }))
      .filter(({ id }) => document.querySelector(`[data-testid="${id}"]`) !== null),
  );
  if (ids.length === 0) return 'no field';
  const { id, kind } = draws.pick(ids);
  const note = await setField(page, id, kind, draws);
  await settle(page);
  return note;
};

/** Sets one field as a user does: a radio, a value of its select (and what the select's "Other" shows), its checkbox, or a value typed. */
async function setField(page: Page, id: string, kind: string, draws: Draws): Promise<string> {
  const control = page.getByTestId(id);
  let note = `${id} (${kind})`;
  if (kind === 'radio') {
    const radios = control.locator('input[type="radio"]');
    const radio = radios.nth(draws.index(await radios.count()));
    note += ` = ${await radio.getAttribute('value')}`;
    await radio.check();
  } else if ((await control.evaluate((element) => element.tagName)) === 'SELECT') {
    const options = await control.locator('option').evaluateAll((list) => list.map((option) => (option as HTMLOptionElement).value));
    const value = draws.pick(options);
    await control.selectOption(value);
    note += ` = ${value}`;
    if (value === 'other') note += await fillOther(page, id, draws);
  } else if ((await control.getAttribute('type')) === 'checkbox') {
    await control.click();
    note += ` ${(await control.isChecked()) ? 'on' : 'off'}`;
    const value = page.getByTestId(`${id}-value`);
    if (kind === 'filter' && (await value.count()) > 0 && draws.next() < 0.5) {
      const typed = valueFor('text', draws);
      await value.fill(typed);
      note += ` ${JSON.stringify(typed)}`;
    }
  } else {
    const typed = valueFor(await control.getAttribute('inputmode'), draws);
    await control.fill(typed);
    note += ` = ${JSON.stringify(typed)}`;
  }
  return note;
}

/**
 * Fills what a select's "Other" shows in its row: the field's own text box (ParameterField, `<id>-other`), or the two
 * fields of a row with NCBI's select of their pairs (ParameterPair: Gap Costs shows the fields of -gapopen and
 * -gapextend, Match/Mismatch Scores those of -reward and -penalty), each set as the field it is. The walk expected a
 * text box of its own for every select and waited for `param-gapcosts-other` (line WQA's sweep 2, W-1, 2026-11-20). A
 * field that is gone when its turn comes (a row of two goes back to its select once both values are emptied) is left.
 */
async function fillOther(page: Page, id: string, draws: Draws): Promise<string> {
  const own = page.getByTestId(`${id}-other`);
  const fields = page.getByTestId(`${id}-field`).locator('[data-testid$="-field"]');
  await expect(own.or(fields).first()).toBeVisible({ timeout: SETTLE_TIMEOUT });
  if ((await own.count()) > 0) {
    const typed = valueFor(await own.getAttribute('inputmode'), draws);
    await own.fill(typed);
    return ` ${JSON.stringify(typed)}`;
  }
  const shown = await fields.evaluateAll((list) =>
    list.map((element) => ({ id: element.getAttribute('data-testid')!.replace(/-field$/, ''), kind: element.getAttribute('data-kind') ?? '' })),
  );
  const notes: string[] = [];
  for (const { id: inner, kind } of shown) notes.push((await page.getByTestId(inner).isVisible()) ? await setField(page, inner, kind, draws) : `${inner} gone`);
  return `: ${notes.join(', ')}`;
}

/** Additional options: options of the CLI that the form has no field for, a few it has, and words it refuses. */
const EXTRAS = Object.freeze([
  '',
  '-max_hsps 1',
  '-max_hsps 2 -culling_limit 1',
  '-outfmt "6 qseqid sseqid pident length evalue bitscore"',
  '-outfmt 10',
  '-outfmt \'7 std qlen slen\'',
  '-qcov_hsp_perc 50',
  '-num_descriptions 5 -num_alignments 3',
  '-subject_besthit',
  '-parse_deflines',
  '-evalue 1e-3',
  '-not_an_option',
  '-max_target_seqs',
  '"unclosed',
]);

const extra: StepBody = async (walker, draws) => {
  const { page } = walker;
  await showTab(page, 'search');
  await settle(page);
  await openParameters(page);
  const text = draws.pick(EXTRAS);
  await page.getByTestId('param-extra').fill(text);
  await settle(page);
  return JSON.stringify(text);
};

const restore: StepBody = async (walker) => {
  const { page } = walker;
  await showTab(page, 'search');
  await settle(page);
  await openParameters(page);
  await page.getByTestId('restore-defaults').click();
  await settle(page);
  return '';
};

const TITLES = Object.freeze(['', 'walk', 'walk <b>bold</b> & "quotes"', `<img src=x onerror="window.__walkXss=1">`, 'タイトル', 'x'.repeat(300)]);

const title: StepBody = async (walker, draws) => {
  const { page } = walker;
  await showTab(page, 'search');
  const text = draws.pick(TITLES);
  await page.getByTestId('job-title').fill(text);
  return JSON.stringify(text.length > 40 ? `${text.slice(0, 40)}…` : text);
};

const threads: StepBody = async (walker, draws) => {
  const { page } = walker;
  await showTab(page, 'search');
  await openParameters(page);
  const select = page.getByTestId('threads');
  const options = await select.locator('option').evaluateAll((list) => list.map((option) => (option as HTMLOptionElement).value));
  const value = draws.pick(options);
  await select.selectOption(value);
  return value;
};

const paste: StepBody = async (walker, draws) => {
  const { page, fixtures } = walker;
  await showTab(page, 'search');
  const role = draws.pick(ROLES);
  const kind = draws.next();
  const which = draws.next();
  // Sometimes the box is emptied.
  const text = draws.next() < 0.1 ? '' : pickFor(fixtures.texts, fixtures.odd, alphabetOf(walker.program, role), kind, which).text;
  await page.getByTestId(`${role}-input`).fill(text);
  walker.model = withPaste(walker.model, role, text);
  // The text is read after a pause: wait until its source is there (or gone) and read.
  await expect
    .poll(async () => (await sourceCards(page, role)).length, { timeout: SETTLE_TIMEOUT })
    .toBe(walker.model[role].sources.length);
  await settle(page);
  return `${role} ${text === '' ? 'emptied' : JSON.stringify(text.split('\n', 1)[0]!.slice(0, 50))}`;
};

const openFiles: StepBody = async (walker, draws) => {
  const { page, fixtures } = walker;
  await showTab(page, 'search');
  const role = draws.pick(ROLES);
  const alphabet = alphabetOf(walker.program, role);
  const count = draws.next() < 0.7 ? 1 : 2;
  const files = [];
  for (let i = 0; i < count; i++) files.push(pickFor(fixtures.files, [], alphabet, draws.next(), draws.next()));
  await page.getByTestId(`${role}-files`).setInputFiles(files.map((file) => ({ name: file.name, mimeType: 'text/plain', buffer: Buffer.from(file.text) })));
  walker.model = withFiles(
    walker.model,
    role,
    files.map((file) => ({ name: file.name, bytes: new TextEncoder().encode(file.text) })),
  );
  await expect
    .poll(async () => (await sourceCards(page, role)).length, { timeout: SETTLE_TIMEOUT })
    .toBe(walker.model[role].sources.length);
  await settle(page);
  return `${role} ${files.map((file) => JSON.stringify(file.name)).join(' ')}`;
};

const remove: StepBody = async (walker, draws) => {
  const { page } = walker;
  await showTab(page, 'search');
  await settle(page);
  const role = draws.pick(ROLES);
  const count = walker.model[role].sources.length;
  if (count === 0) return `${role}: no input`;
  const index = draws.index(count);
  await page.getByTestId(`${role}-source-${index}-remove`).click();
  walker.model = withoutSource(walker.model, role, index);
  // A removed pasted text empties the box.
  await settle(page);
  return `${role} input ${index}`;
};

const records: StepBody = async (walker, draws) => {
  const { page } = walker;
  await showTab(page, 'search');
  await settle(page);
  const role = draws.pick(ROLES);
  const cards = await sourceCards(page, role);
  const ready = cards.flatMap((card, index) => (card.status === 'ready' ? [index] : []));
  if (ready.length === 0) return `${role}: no input read`;
  const index = draws.pick(ready);
  const details = page.getByTestId(`${role}-source-${index}`).locator('details.records');
  if ((await details.count()) === 0) return `${role} input ${index}: no record list`;
  await details.evaluate((element) => ((element as HTMLDetailsElement).open = true));
  const boxes = page.locator(`[data-testid^="${role}-source-${index}-record-"]`);
  const shown = await boxes.count();
  if (shown === 0) return `${role} input ${index}: no record shown`;
  const box = boxes.nth(draws.index(shown));
  const record = Number(/-record-(\d+)$/.exec((await box.getAttribute('data-testid')) ?? '')?.[1]);
  const wasIncluded = await box.isChecked();
  await box.click();
  walker.model = withIncluded(walker.model, role, alphabetOf(walker.program, role), index, record, !wasIncluded);
  await settle(page);
  return `${role} input ${index} record ${record} ${wasIncluded ? 'left out' : 'included'}`;
};

const combine: StepBody = async (walker, draws) => {
  const { page, fixtures } = walker;
  await showTab(page, 'search');
  await settle(page);
  const role = draws.pick(ROLES);
  let note = '';
  if ((await page.getByTestId(`${role}-mode`).count()) === 0) {
    // Two inputs read: open files of the role's alphabet.
    const read = (await sourceCards(page, role)).filter((card) => card.status === 'ready').length;
    const files = Array.from({ length: Math.max(1, 2 - read) }, () => pickFor(fixtures.files, [], alphabetOf(walker.program, role), 0, draws.next()));
    await page.getByTestId(`${role}-files`).setInputFiles(files.map((file) => ({ name: file.name, mimeType: 'text/plain', buffer: Buffer.from(file.text) })));
    walker.model = withFiles(
      walker.model,
      role,
      files.map((file) => ({ name: file.name, bytes: new TextEncoder().encode(file.text) })),
    );
    await expect.poll(async () => (await sourceCards(page, role)).length, { timeout: SETTLE_TIMEOUT }).toBe(walker.model[role].sources.length);
    await settle(page);
    note = `opened ${files.map((file) => JSON.stringify(file.name)).join(' ')}, `;
    if ((await page.getByTestId(`${role}-mode`).count()) === 0) return `${note}${role}: fewer than two inputs read`;
  }
  const next = walker.model[role].mode === 'combined' ? 'separate' : 'combined';
  await page.getByTestId(`${role}-mode-${next}`).check();
  walker.model = withMode(walker.model, role, next);
  await settle(page);
  return `${note}${role} ${next}`;
};

const region: StepBody = async (walker, draws) => {
  const { page } = walker;
  await showTab(page, 'search');
  await settle(page);
  const role = draws.pick(ROLES);
  const picker = page.getByTestId(`${role}-region`);
  if ((await picker.count()) === 0) return `${role}: no single record`;
  const length = Number(((await picker.locator('.region-scale span').last().textContent()) ?? '').replace(/[^0-9]/g, ''));
  const action = draws.next();
  if (action < 0.15 && (await page.getByTestId(`${role}-region-clear`).count()) > 0) {
    await page.getByTestId(`${role}-region-clear`).click();
    await settle(page);
    return `${role} cleared`;
  }
  // Mostly a range in the record; sometimes one past its end or the wrong way round.
  let start = 1 + Math.floor(draws.next() * Math.max(1, length - 1));
  let stop = start + Math.floor(draws.next() * (length - start + 1));
  const odd = draws.next();
  if (odd < 0.1) stop = length + 5;
  else if (odd < 0.2) [start, stop] = [stop, start];
  await page.getByTestId(`${role}-region-start`).fill(String(start));
  await page.getByTestId(`${role}-region-stop`).fill(String(stop));
  await settle(page);
  return `${role} ${start}-${stop} of ${length}`;
};

// --- runs -----------------------------------------------------------------------------------

/** What a run must have used, fixed when Run LOSAT is pressed. */
interface Submitted {
  readonly inputs: Readonly<Record<Role, readonly ExpectedInput[]>>;
  /** The flags of each field marked changed (♦): a field of two options (scores, gap costs) has two. */
  readonly changed: ReadonlyArray<readonly string[]>;
  /** The flags of Additional options' words. */
  readonly extraFlags: readonly string[];
  /** The value that a changed field of one option shows in its text box or select (its "other" box). */
  readonly values: ReadonlyMap<string, string>;
  /** The runs that the inputs make: one per input of each role searched apart. */
  readonly expectedRuns: number;
}

async function submitted(walker: Walker): Promise<Submitted> {
  const { page } = walker;
  await openParameters(page);
  const fields = await page.evaluate(() =>
    [...document.querySelectorAll('[data-testid="search-panel"] .changed[data-testid$="-field"]')].map((field) => {
      const flags = [...field.querySelectorAll('.flag')].flatMap((element) => (element.textContent ?? '').split(/\s+/).filter((word) => word.startsWith('-')));
      const id = (field.getAttribute('data-testid') ?? '').replace(/-field$/, '');
      const control = document.querySelector(`[data-testid="${id}"]`);
      let value: string | undefined;
      if (control instanceof HTMLSelectElement) {
        value = control.value === 'other' ? (document.querySelector<HTMLInputElement>(`[data-testid="${id}-other"]`)?.value ?? undefined) : control.value;
      } else if (control instanceof HTMLInputElement && (control.type === 'text' || control.type === 'number')) {
        value = control.value;
      }
      return { flags, value };
    }),
  );
  // A field with no option of its own (Job Title, Additional options) is no flag of the argv.
  const changed = fields.map((field) => field.flags).filter((flags) => flags.length > 0);
  const values = new Map(fields.flatMap((field) => (field.flags.length === 1 && field.value !== undefined ? [[field.flags[0]!, field.value.trim()] as const] : [])));
  const extraText = await page.getByTestId('param-extra').inputValue();
  let extraFlags: string[] = [];
  try {
    extraFlags = parseOptions(extraText).filter((word) => word.startsWith('-') && !/^-\d/.test(word)).map((word) => word.split('=', 1)[0]!);
  } catch {
    extraFlags = [];
  }
  const inputs = {} as Record<Role, readonly ExpectedInput[]>;
  for (const role of ROLES) {
    const cards = await sourceCards(page, role);
    const ready = cards.flatMap((card, index) => (card.status === 'ready' ? [index] : []));
    inputs[role] = expectedInputs(walker.model[role], role, alphabetOf(walker.program, role), ready);
  }
  const expectedRuns = ROLES.reduce((product, role) => product * Math.max(1, inputs[role].length), 1);
  return { inputs, changed, extraFlags, values, expectedRuns };
}

/** The words of Additional options as a POSIX shell reads them (single and double quotes). */
function parseOptions(text: string): string[] {
  const words: string[] = [];
  let word: string | undefined;
  for (let i = 0; i < text.length; i++) {
    const c = text[i]!;
    if (/\s/.test(c)) {
      if (word !== undefined) words.push(word);
      word = undefined;
    } else if (c === "'" || c === '"') {
      const end = text.indexOf(c, i + 1);
      if (end < 0) throw new Error('unclosed quote');
      word = (word ?? '') + text.slice(i + 1, end);
      i = end;
    } else word = (word ?? '') + c;
  }
  if (word !== undefined) words.push(word);
  return words;
}

/**
 * Presses Run LOSAT: the runs it starts (their numbers), or none with the reason the page gives.
 * A press that starts nothing must say why.
 */
async function pressRun(walker: Walker, before: Submitted): Promise<{ runs: number[]; refused?: string }> {
  const { page } = walker;
  await showTab(page, 'search');
  const label = (await page.getByTestId('run-losat').textContent()) ?? '';
  const count = Number(/\((\d+) runs\)/.exec(label)?.[1] ?? '1');
  if (count !== before.expectedRuns) {
    throw new StepFailure('run-count', `Run LOSAT says ${count} runs ("${label.trim()}"); the inputs the walk gave make ${before.expectedRuns}`);
  }
  await page.getByTestId('run-losat').click();
  const outcome = await expect
    .poll(
      async () => {
        if (await tabShown(page, 'results')) return 'started';
        const message = page.getByTestId('search-message');
        if ((await message.count()) > 0 && (await message.isVisible())) return 'refused';
        return 'waiting';
      },
      { timeout: 15_000 },
    )
    .not.toBe('waiting')
    .then(
      () => undefined,
      () => 'none',
    );
  if (outcome === 'none' && !(await tabShown(page, 'results'))) {
    const message = page.getByTestId('search-message');
    if ((await message.count()) === 0) throw new StepFailure('silent-refusal', 'Run LOSAT started no run and said nothing');
  }
  if (await tabShown(page, 'results')) {
    for (const role of ROLES) {
      if (before.inputs[role].every((input) => input.bytes.length === 0)) {
        throw new StepFailure('run-without-input', `Run LOSAT started a search, but the walk gave the ${role} no record to search`);
      }
    }
    const runs = Array.from({ length: count }, (_, i) => walker.runs + 1 + i);
    walker.runs += count;
    return { runs };
  }
  const refused = ((await page.getByTestId('search-message').textContent()) ?? '').trim();
  if (refused === '') throw new StepFailure('silent-refusal', 'Run LOSAT started no run and its message is empty');
  judgeRefusal(walker, refused);
  walker.counts.refused += 1;
  return { runs: [], refused };
}

/** A refusal whose reason the walk knows the truth of must agree with the inputs that it gave. */
function judgeRefusal(walker: Walker, message: string): void {
  for (const role of ROLES) {
    const title = role === 'query' ? 'query' : 'subject';
    const sources = walker.model[role].sources;
    if (message.startsWith(`Add the ${title} sequences`) && sources.length > 0) {
      throw new StepFailure('refusal-wrong', `Run LOSAT asked for the ${role} ("${message}"), but the walk gave it ${sources.length} inputs`);
    }
    if (message === `Every ${title} record is excluded.`) {
      const reading = alphabetOf(walker.program, role);
      const left = sources.filter((source) => (source.excluded[reading]?.size ?? 0) === 0);
      if (left.length > 0) {
        throw new StepFailure('refusal-wrong', `Run LOSAT said every ${role} record is excluded, but the walk left out none of ${left.map((source) => source.name).join(', ')}`);
      }
    }
  }
}

async function runStatus(page: Page, run: number, timeout = RUN_TIMEOUT): Promise<string> {
  await showTab(page, 'results');
  const option = page.getByTestId(`run-${run}-option`);
  try {
    await expect(option).toHaveAttribute('data-status', /^(completed|cancelled|failed)$/, { timeout });
  } catch {
    throw new StepFailure('run-hang', `run ${run} did not end in ${timeout / 1000} s (${await option.getAttribute('data-status')})`);
  }
  return (await option.getAttribute('data-status'))!;
}

async function pickRun(page: Page, run: number): Promise<void> {
  await showTab(page, 'results');
  const value = await page.getByTestId(`run-${run}-option`).getAttribute('value');
  await page.getByTestId('results-run').selectOption(value!);
}

/** What the run's details show: its inputs (name), its argv's option names, and the command line of each output. */
async function runDetails(page: Page, run: number): Promise<{ inputs: Record<Role, string>; options: string[]; commands: Map<string, string[]> }> {
  await pickRun(page, run);
  await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', String(run), { timeout: SETTLE_TIMEOUT });
  await page.getByTestId('results-view-details').click();
  const inputs = {} as Record<Role, string>;
  for (const role of ROLES) inputs[role] = ((await page.getByTestId(`run-input-${role}`).textContent()) ?? '').trim();
  const options = await page.getByTestId('run-argv').locator('.option-name').allTextContents();
  const commands = new Map<string, string[]>();
  for (const format of [...FORMATS.map(String), 'custom']) {
    const command = page.getByTestId(`run-command-${format}`);
    if ((await command.count()) > 0) commands.set(format, parseShellWords(((await command.textContent()) ?? '').trim()));
  }
  return { inputs, options, commands };
}

/** The bytes of the run's own -outfmt, saved from Download All. */
async function customOutput(page: Page): Promise<Buffer | undefined> {
  await page.getByTestId('results-download-all').click();
  const item = page.getByTestId('results-download-menu').getByTestId('results-download-custom');
  if ((await item.count()) === 0) {
    await page.keyboard.press('Escape');
    return undefined;
  }
  const [download] = await Promise.all([page.waitForEvent('download'), item.click()]);
  return readFileSync((await download.path())!);
}

/**
 * The file that the page offers for a role of the run chosen (Download All: the query or subject that is a file of LOSAT's, line
 * WQA's D22), saved as the user gets it: its name and bytes. Undefined where the input is the user's file as given.
 */
async function offeredInput(page: Page, role: Role): Promise<{ name: string; bytes: Buffer } | undefined> {
  await page.getByTestId('results-download-all').click();
  const item = page.getByTestId('results-download-menu').getByTestId(`results-download-${role}`);
  if ((await item.count()) === 0) {
    await page.keyboard.press('Escape');
    return undefined;
  }
  const [download] = await Promise.all([page.waitForEvent('download'), item.click()]);
  return { name: download.suggestedFilename(), bytes: readFileSync((await download.path())!) };
}

/** A file that the native CLI is given: its name, its bytes, and whether the page made it. */
interface CliFile {
  readonly name: string;
  readonly bytes: Uint8Array;
  readonly made: boolean;
}

/**
 * The file of each role that the command line names, for the native CLI, for run number `at` of the press's group (the draft makes
 * the group with the query outer and the subject inner: src/application/draft.ts `prepareAndEnqueue`). The walk's model says, from what
 * the walk gave, what the run's input is (`expectedInputs`) and whether it is made or the user's file as given. The page must offer a
 * file exactly for a made input (Download All), under a name that is none of the user's file names (it would name other bytes) and that
 * Run details shows, holding exactly the run's own expected bytes; it must offer none for a file as given, which the CLI gets by its
 * own name. The CLI then runs the command that the page shows on the files that the page offers, as the user would. The names of the
 * made files are not modelled here, but `query.fa` and `combined_*.fa`, whose role may be put before them.
 */
async function filesOf(walker: Walker, run: number, before: Submitted, at: number, shown: Record<Role, string>): Promise<Record<Role, CliFile> | undefined> {
  const subjects = Math.max(1, before.inputs.subject.length);
  const expected: Record<Role, ExpectedInput | undefined> = {
    query: before.inputs.query[Math.floor(at / subjects)],
    subject: before.inputs.subject[at % subjects],
  };
  const userNames = new Set(ROLES.flatMap((role) => walker.model[role].sources.filter((source) => source.origin === 'file').map((source) => source.name)));
  const files = {} as Record<Role, CliFile>;
  for (const role of ROLES) {
    const input = expected[role];
    if (input === undefined) throw new StepFailure('input-missing', `run ${run} (place ${at} of its group) has no ${role} among the walk's inputs`);
    const offered = await offeredInput(walker.page, role);
    if (!input.made) {
      if (offered !== undefined) throw new StepFailure('given-offered', `run ${run}: the ${role} is the user's file ${input.name}, but the page offers a file ${offered.name} for it`);
      if (!shown[role].startsWith(`${input.name}: `)) {
        throw new StepFailure('input-unmatched', `run ${run}: the ${role} the walk gave is the file ${input.name}; Run details say "${shown[role]}"`);
      }
      files[role] = { name: input.name!, bytes: input.bytes, made: false };
      continue;
    }
    if (offered === undefined) {
      throw new StepFailure('made-not-offered', `run ${run}: the ${role} is not the user's file as given (pasted, several files, or records left out), but the page offers no file for it ("${shown[role]}")`);
    }
    if (userNames.has(offered.name)) {
      throw new StepFailure('offered-user-name', `run ${run}: the page offers ${offered.name} as the ${role} file, a name that one of the user's files has`);
    }
    if (input.name !== undefined && offered.name !== input.name && offered.name !== `${role}_${input.name}`) {
      throw new StepFailure('offered-name', `run ${run}: the page offers ${offered.name} as the ${role} file; the walk expects ${input.name}`);
    }
    if (!shown[role].startsWith(`${offered.name}: `)) {
      throw new StepFailure('offered-name', `run ${run}: the page offers ${offered.name} as the ${role} file, but the run's details say "${shown[role]}"`);
    }
    if (Buffer.compare(Buffer.from(input.bytes), offered.bytes) !== 0) {
      throw new StepFailure('offered-bytes', `run ${run}: the ${role} file ${offered.name} that the page offers (${offered.bytes.length} bytes) is not what the walk gave this run (${input.bytes.length} bytes)`);
    }
    walker.counts.offered += 1;
    // Of those, the files with records left out (Run details says "without").
    if (((await walker.page.getByTestId(`run-input-${role}-from`).textContent().catch(() => null)) ?? '').includes(' without ')) walker.counts.offeredLeftOut += 1;
    files[role] = { name: offered.name, bytes: offered.bytes, made: true };
  }
  if (!files.query.made && !files.subject.made && files.query.name === files.subject.name && Buffer.compare(Buffer.from(files.query.bytes), Buffer.from(files.subject.bytes)) !== 0) {
    // Two of the user's own files of one name with different records: the command names both by one name (the user's own naming, D21).
    walker.notes.push(`run ${run}: the user's files ${files.query.name} differ and are named alike; native check skipped`);
    walker.counts.unmatched += 1;
    return undefined;
  }
  return files;
}

/**
 * Checks a run that ended: a completed one against the native CLI and the form; a failed one must
 * be one that the native CLI refuses too. Returns its outputs (format → bytes) when it completed.
 */
async function verifyRun(walker: Walker, run: number, before: Submitted, at: number): Promise<Map<string, Buffer> | undefined> {
  const { page } = walker;
  const status = await runStatus(page, run);
  walker.counts[status as 'completed' | 'failed' | 'cancelled'] += 1;
  if (status === 'cancelled') return undefined;
  if (status === 'failed') {
    // The form checks its options and inputs with the engine before a search starts: a search that
    // it let through and that then fails is a defect (or a check the form lacks).
    await pickRun(page, run);
    const text = ((await page.getByTestId('results-status').textContent()) ?? '').trim();
    throw new StepFailure('run-failed', `run ${run} failed: ${text}`);
  }
  const details = await runDetails(page, run);
  // The options of the argv: the changed fields, the words of Additional options and the subranges.
  // A changed field of two options needs one of them: the other may keep its default (scores
  // "1,-3": -penalty only).
  const options = new Set(details.options.filter((flag) => flag !== '-query' && flag !== '-subject'));
  const allowed = new Set([...before.changed.flat(), ...before.extraFlags, '-query_loc', '-subject_loc']);
  const unexpected = [...options].filter((flag) => !allowed.has(flag));
  const missing = before.changed.filter((flags) => !flags.some((flag) => options.has(flag))).map((flags) => flags.join('/'));
  if (unexpected.length > 0 || missing.length > 0) {
    throw new StepFailure(
      'form-argv',
      `run ${run}: the argv's options ${[...options].join(' ')} differ from the fields marked changed (${before.changed.map((flags) => flags.join('/')).join(' ')})` +
        ` and Additional options (${before.extraFlags.join(' ')}): unexpected ${unexpected.join(' ') || 'none'}, missing ${missing.join(' ') || 'none'}`,
    );
  }
  // A changed field's value, as its box or select shows it, is its option's value in the argv.
  const argv = parseShellWords(((await page.getByTestId('run-argv').textContent()) ?? '').trim());
  for (const [flag, value] of before.values) {
    const at = argv.indexOf(flag);
    if (at < 0 || argv[at + 1] !== value) {
      throw new StepFailure('form-argv-value', `run ${run}: the form shows ${flag} ${JSON.stringify(value)}; its argv has ${at < 0 ? 'no such option' : JSON.stringify(argv[at + 1])}`);
    }
  }
  const { outputs } = await storedOutputs(page, run);
  const stored = new Map<string, Buffer>(FORMATS.map((format) => [String(format), outputs[format]]));
  if (details.commands.has('custom')) {
    const custom = await customOutput(page);
    if (custom === undefined) throw new StepFailure('outputs', `run ${run}: its own -outfmt has a command but no download`);
    stored.set('custom', custom);
  }
  if (NATIVE === undefined) return stored;
  const files = await filesOf(walker, run, before, at, details.inputs);
  if (files === undefined) return stored;
  const { query, subject } = files;
  if (query.name === subject.name && Buffer.compare(Buffer.from(query.bytes), Buffer.from(subject.bytes)) !== 0) {
    throw new StepFailure('command-names', `run ${run}: the command line names ${query.name} for both inputs, which differ`);
  }
  const directory = join(walker.outDir, `run-${run}-${Date.now()}`);
  mkdirSync(directory, { recursive: true });
  writeFileSync(join(directory, query.name), query.bytes);
  writeFileSync(join(directory, subject.name), subject.bytes);
  for (const [format, words] of details.commands) {
    if (words[0] !== 'losat') throw new StepFailure('command', `run ${run}: the command of ${format} does not start with losat: ${words.join(' ')}`);
    const cli = spawnSync(NATIVE, words.slice(1), { cwd: directory, maxBuffer: 1 << 28, timeout: 60_000 });
    if (cli.status !== 0) {
      throw new StepFailure('native-refuses', `run ${run}: the run completed, but its command of outfmt ${format} fails in the CLI (${cli.status}): ${String(cli.stderr).slice(0, 400)}\n${words.join(' ')}`);
    }
    const web = stored.get(format);
    if (web === undefined || !web.equals(cli.stdout)) {
      writeFileSync(join(directory, `web.${format}.out`), web ?? '');
      writeFileSync(join(directory, `native.${format}.out`), cli.stdout);
      throw new StepFailure(
        `native-bytes:${format}`,
        `run ${run}: the stored outfmt ${format} (${web?.length ?? 'none'} bytes) differs from the CLI's (${cli.stdout.length} bytes) for ${words.join(' ')} (files in ${directory})`,
      );
    }
    walker.counts.compared += 1;
  }
  return stored;
}

const runOnce: StepBody = async (walker) => {
  const { page } = walker;
  await showTab(page, 'search');
  await settle(page);
  const before = await submitted(walker);
  const { runs, refused } = await pressRun(walker, before);
  if (refused !== undefined) return `refused: ${refused}`;
  const ended = [];
  for (const [at, run] of runs.entries()) {
    await verifyRun(walker, run, before, at);
    ended.push(`${run} ${await page.getByTestId(`run-${run}-option`).getAttribute('data-status')}`);
  }
  return `runs ${ended.join(', ')}`;
};

const twice: StepBody = async (walker) => {
  const { page } = walker;
  await showTab(page, 'search');
  await settle(page);
  const before = await submitted(walker);
  const first = await pressRun(walker, before);
  if (first.refused !== undefined) return `refused: ${first.refused}`;
  const second = await pressRun(walker, before);
  if (second.refused !== undefined) throw new StepFailure('twice', `the second Run LOSAT of the same form was refused: ${second.refused}`);
  if (second.runs.length !== first.runs.length) throw new StepFailure('twice', `the same form started ${first.runs.length} then ${second.runs.length} runs`);
  const notes = [];
  for (const [i, run] of first.runs.entries()) {
    const a = await verifyRun(walker, run, before, i);
    const b = await verifyRun(walker, second.runs[i]!, before, i);
    if (a !== undefined && b !== undefined) {
      for (const [format, bytes] of a) {
        if (!bytes.equals(b.get(format) ?? Buffer.alloc(0))) {
          throw new StepFailure(`twice-bytes:${format}`, `runs ${run} and ${second.runs[i]} of the same form differ in outfmt ${format}`);
        }
      }
    }
    notes.push(`${run}=${second.runs[i]}`);
  }
  return `runs ${notes.join(', ')}`;
};

const cancel: StepBody = async (walker) => {
  const { page } = walker;
  await showTab(page, 'search');
  await settle(page);
  const before = await submitted(walker);
  const { runs, refused } = await pressRun(walker, before);
  if (refused !== undefined) return `refused: ${refused}`;
  const group = page.getByTestId('run-status-cancel-group');
  const one = page.getByTestId('run-status-cancel');
  if ((await group.count()) > 0) await group.click().catch(() => undefined);
  else if ((await one.count()) > 0) await one.click().catch(() => undefined);
  const ended = [];
  for (const [at, run] of runs.entries()) {
    await verifyRun(walker, run, before, at);
    ended.push(`${run} ${await page.getByTestId(`run-${run}-option`).getAttribute('data-status')}`);
  }
  return `runs ${ended.join(', ')}`;
};

// --- the page -------------------------------------------------------------------------------

const tab: StepBody = async (walker) => {
  const { page } = walker;
  const next = (await tabShown(page, 'search')) ? 'results' : 'search';
  await showTab(page, next);
  return next;
};

const view: StepBody = async (walker, draws) => {
  const { page } = walker;
  await showTab(page, 'results');
  if (walker.runs === 0) {
    await expect(page.getByTestId('results-empty')).toBeVisible();
    return 'no run';
  }
  const run = 1 + draws.index(walker.runs);
  const status = await runStatus(page, run);
  await pickRun(page, run);
  if (status !== 'completed') {
    await expect(page.getByTestId('results-status')).toHaveAttribute('data-phase', 'unavailable');
    return `run ${run} ${status}`;
  }
  try {
    await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', String(run), { timeout: SETTLE_TIMEOUT });
  } catch {
    const text = ((await page.getByTestId('results-status').textContent().catch(() => '')) ?? '').trim();
    throw new StepFailure('results-unread', `run ${run} completed, but its results were not shown: ${text}`);
  }
  const views = ['hits', 'dotplot', 'details', 'outputs'] as const;
  const chosen = draws.pick(views);
  const testid = { hits: 'results-view-hits', dotplot: 'pane-dotplot', details: 'results-view-details', outputs: 'results-view-outputs' }[chosen];
  await page.getByTestId(testid).click();
  let note = `run ${run} ${chosen}`;
  if (chosen === 'outputs') {
    const formats = page.locator('[data-testid^="format-"]');
    const count = await formats.count();
    if (count > 0) {
      const button = formats.nth(draws.index(count));
      note += ` ${await button.getAttribute('data-testid')}`;
      await button.click();
    }
  } else if (chosen === 'hits') {
    const sorts = page.locator('[data-testid*="-sort-"]:visible');
    const count = await sorts.count();
    if (count > 0 && draws.next() < 0.6) {
      const button = sorts.nth(draws.index(count));
      note += ` ${await button.getAttribute('data-testid')}`;
      await button.click();
    }
  }
  return note;
};

const width: StepBody = async (walker, draws) => {
  const { page } = walker;
  const value = draws.pick(WIDTHS);
  await page.setViewportSize({ width: value, height: 900 });
  // The layout settles in a frame.
  await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));
  return `${value} px`;
};

const reloadPage: StepBody = async (walker) => {
  const { page } = walker;
  await reload(page);
  await described(page);
  // LOSAT Web keeps nothing for the next visit (REQ-20): the page starts again, with no input and no run.
  for (const role of ROLES) {
    if ((await sourceCards(page, role)).length > 0) throw new StepFailure('reload', `the ${role} kept its inputs over a reload`);
  }
  await showTab(page, 'results');
  await expect(page.getByTestId('results-empty')).toBeVisible();
  await showTab(page, 'search');
  await start(walker);
  return walker.program;
};

export const STEPS: Readonly<Record<WalkKind, StepBody>> = Object.freeze({
  mode,
  task,
  field,
  extra,
  restore,
  title,
  threads,
  paste,
  open: openFiles,
  remove,
  records,
  combine,
  region,
  run: runOnce,
  twice,
  cancel,
  tab,
  view,
  width,
  reload: reloadPage,
});
