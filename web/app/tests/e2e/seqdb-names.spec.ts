// Database names with a space or a double quote in every option that LOSAT Web writes a database name into (line WW,
// WW4d): the command line reads each of them as NCBI's SeqDB splits a list (a space separates two names, a double quote
// starts or ends a name with spaces), measured against NCBI BLAST+ 2.17.0 for every option (WW/ww4d/probe). The Web
// writes a name with a space between double quotes (`"my dbs/ndb"`), and refuses a name with a double quote where the
// database is chosen, with the reason.
//
// With either engine: the refusal in the Database picker, the Repeats database, WindowMasker's and makeblastdb's Input
// blocks and the Masking tab's, the check and Run LOSAT (for as long as the refused database is chosen or its files are
// in the run), and a
// database named with a space shown in Run details as a command that works when pasted. With the real engine
// (LOSAT_WEB_REACTORS) and the native binary (LOSAT_WEB_NATIVE): blastn -db on a folder `my dbs` (and an alias `my
// nslice`), BLASTN's Repeats database in a folder `rep dbs`, WindowMasker and dustmasker -infmt blastdb on a folder
// `my db`, and
// makeblastdb's -out and -input_type blastdb with a space write what the native command line writes on the same files.
//
// LOSAT_E2E_SCREENSHOT_DIR=<dir> saves the changed screens at 1280 and 390 px.
import { spawnSync } from 'node:child_process';
import { copyFileSync, mkdirSync, mkdtempSync, readdirSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { nativeRun } from './support/makeblastdb-native';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { checked, described, openParameters, openSearch, paste, pickRun, program, result, submit, waitStatus } from './support/search';
import { test } from './support/test';

test.setTimeout(300_000);

const LOSAT = join(REPOSITORY, 'LOSAT');
const SCREENSHOT_DIR = process.env['LOSAT_E2E_SCREENSHOT_DIR'];
const DESKTOP = { width: 1280, height: 900 };
const MOBILE = { width: 390, height: 844 };
const REASON = 'a double quote there starts or ends a name';
const ZONE = 'Asia/Tokyo';
test.use({ timezoneId: ZONE });

/** Saves the page at 1280 and 390 px, back at 1280 px after. */
async function shots(page: Page, name: string): Promise<void> {
  if (!SCREENSHOT_DIR) return;
  mkdirSync(SCREENSHOT_DIR, { recursive: true });
  const project = test.info().project.name;
  for (const viewport of [DESKTOP, MOBILE]) {
    await page.setViewportSize(viewport);
    await page.screenshot({ path: join(SCREENSHOT_DIR, `${project}-${name}-${viewport.width}.png`), fullPage: true });
  }
  await page.setViewportSize(DESKTOP);
}

/** A file whose bytes do not matter (the picker reads the names, and an alias's text). */
const named = (name: string, text = name) => ({ name, mimeType: 'application/octet-stream', buffer: Buffer.from(text) });

/** Copies the files of `names` from `from` into the folder `to` (made). */
function copyInto(from: string, names: readonly string[], to: string): void {
  mkdirSync(to, { recursive: true });
  for (const name of names) copyFileSync(join(from, name), join(to, name));
}

const filesOf = (dir: string, prefix: string) => readdirSync(dir).filter((name) => name.startsWith(`${prefix}.`) && /\.n[a-z]{2}$/.test(name));

/** The native command line's run of `args` in `cwd`: its exit status and streams. */
function native(args: readonly string[], cwd: string): { readonly status: number; readonly stdout: string; readonly stderr: string } {
  const run = spawnSync(NATIVE!, args, { cwd, maxBuffer: 1 << 28, env: { ...process.env, TZ: ZONE } });
  if (run.error) throw run.error;
  return { status: run.status ?? -1, stdout: run.stdout.toString('utf8'), stderr: run.stderr.toString('utf8') };
}

/** The native search's outfmt `format` (one thread), which must succeed. */
function nativeSearch(args: readonly string[], cwd: string, format: 0 | 6 | 7): string {
  const run = native([...args, '-outfmt', String(format), '-num_threads', '1'], cwd);
  if (run.status !== 0) throw new Error(`native losat ${args.join(' ')} exited ${run.status}: ${run.stderr}`);
  return run.stdout;
}

/** The exit status, standard output and standard error of tool run `run`, from its Outputs. */
async function toolRun(page: Page, run: number, status: 'completed' | 'failed'): Promise<{ exit: number; stdout: string; stderr: string }> {
  await waitStatus(page, run, status, 300_000);
  await pickRun(page, run);
  await expect(page.getByTestId('tool-run-view')).toHaveAttribute('data-run', String(run));
  const download = async (tab: string): Promise<string> => {
    if ((await page.getByTestId(tab).count()) === 0) return '';
    await page.getByTestId(tab).click();
    const [file] = await Promise.all([page.waitForEvent('download'), page.getByTestId('tool-output-download').click()]);
    return readFileSync((await file.path())!).toString('utf8');
  };
  const exitText = status === 'failed' ? ((await page.getByTestId('tool-exit-status').textContent()) ?? '') : 'Exit status 0';
  const exit = Number(/Exit status (\d+)/.exec(exitText)?.[1] ?? -1);
  return { exit, stdout: await download('tool-tab-stdout'), stderr: await download('tool-tab-stderr') };
}

test('a database named with a double quote is refused where it is chosen, with the reason; one named with a space is searched as "my ndb" (line WW, WW4d)', async ({
  page,
}) => {
  await openSearch(page);
  await program(page, 'blastn');
  await paste(page, 'query', '>q1\nACGTACGTACGTACGTACGTACGTACGTACGT\n');

  // Search's Database (-db): the picker and Run LOSAT say why.
  await page.getByTestId('subject-source-database').check();
  await page.getByTestId('subject-db-files-input').setInputFiles([named('q"db.nin'), named('q"db.nhr'), named('q"db.nsq')]);
  await expect(page.getByTestId('subject-db-selected')).toContainText('q"db');
  const dbRefused = `A database named "q"db" cannot be in -db: ${REASON}.`;
  await expect(page.getByTestId('subject-db-refused')).toHaveText(dbRefused);
  // The check's line, which is live, says it; the picker's line is not a second alert of the same sentence (WW4d R2, N3).
  await expect(page.getByTestId('subject-db-refused')).not.toHaveAttribute('role', 'alert');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'invalid');
  await expect(page.getByTestId('argv-validation')).toHaveText(dbRefused);
  await shots(page, 'db-refused');

  // BLASTN's Repeats database (-filtering_db).
  await page.getByTestId('subject-source-fasta').check();
  await paste(page, 'subject', '>s1\nACGTACGTACGTACGTACGTACGTACGTACGT\n');
  await openParameters(page);
  await page.getByTestId('param-filtering_db').check();
  await page.getByTestId('param-filtering_db-files-input').setInputFiles([named('q"rep.nin'), named('q"rep.nhr'), named('q"rep.nsq')]);
  const repRefused = `A database named "q"rep" cannot be in -filtering_db: ${REASON}.`;
  await expect(page.getByTestId('param-filtering_db-refused')).toHaveText(repRefused);
  await expect(page.getByTestId('param-filtering_db-refused')).not.toHaveAttribute('role', 'alert');
  // The check says it too, as it says -db's (WW4d R1).
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'invalid');
  await expect(page.getByTestId('argv-validation')).toHaveText(repRefused);
  await page.getByTestId('param-filtering_db-refused').scrollIntoViewIfNeeded();
  await shots(page, 'repeats-refused');
  await page.getByTestId('run-losat').click();
  await expect(page.getByTestId('search-message')).toHaveText(repRefused);
  await page.getByTestId('param-filtering_db').uncheck();

  // A database named with a space: chosen as it is called, written as the command line reads it, and the command of
  // Run details pastes.
  await page.getByTestId('subject-source-database').check();
  await page.getByTestId('subject-db-files-input').setInputFiles([named('my ndb.nin'), named('my ndb.nhr'), named('my ndb.nsq')]);
  await expect(page.getByTestId('subject-db-selected')).toContainText('my ndb');
  await expect(page.getByTestId('subject-db-refused')).toHaveCount(0);
  if (!BUILD_HAS_ENGINE) {
    // The fake engine searches nothing; the real engine's search of such a database is the test below.
    await submit(page);
    await waitStatus(page, 1, 'completed');
    await pickRun(page, 1);
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-details')).toContainText('my ndb');
    await expect(page.getByTestId('run-command-6')).toContainText(`-db '"my ndb"'`);
    await page.getByTestId('run-command-6').scrollIntoViewIfNeeded();
    await shots(page, 'run-details-space');
  }

  // WindowMasker's -in of -infmt blastdb.
  await page.getByTestId('tab-windowmasker').click();
  await described(page, 'wm-panel', 'ready');
  await page.getByTestId('wm-param-infmt').selectOption('blastdb');
  await page.getByTestId('wm-in-files').setInputFiles([named('q"wm.nin'), named('q"wm.nsq')]);
  const wmRefused = `A database named "q"wm" cannot be in -in: ${REASON}.`;
  await expect(page.getByTestId('wm-in-refused')).toHaveText(wmRefused);
  await expect(page.getByTestId('wm-in-refused')).not.toHaveAttribute('role', 'alert');
  await expect(page.getByTestId('wm-in-name')).toHaveValue('');
  await checked(page, 'wm-panel', 'invalid');
  await expect(page.getByTestId('wm-argv-validation')).toHaveText(wmRefused);
  await shots(page, 'wm-refused');
  // The refusal lasts as long as its cause (WW4d R1): other edits keep it, and Run LOSAT says it.
  await page.getByTestId('wm-stage-mask').check();
  await expect(page.getByTestId('wm-in-refused')).toHaveCount(0);
  await page.getByTestId('wm-stage-mk_counts').check();
  await expect(page.getByTestId('wm-in-refused')).toHaveText(wmRefused);
  await expect(page.getByTestId('wm-argv-validation')).toHaveText(wmRefused);
  await page.getByTestId('wm-run-losat').click();
  await expect(page.getByTestId('wm-message')).toHaveText(wmRefused);
  // The refused database's files stay in the run's files until they are removed, as makeblastdb's do.
  const toggle = page.getByTestId('wm-options-toggle');
  if ((await toggle.getAttribute('aria-expanded')) !== 'true') await toggle.click();
  const runFiles = page.locator('[data-testid="wm-run-files-list"] > li');
  await expect(runFiles).toHaveCount(2);
  while ((await runFiles.count()) > 0) await page.getByTestId('wm-run-file-remove-0').click();
  await page.getByTestId('wm-in-files').setInputFiles([named('my wm.nin'), named('my wm.nsq')]);
  await expect(page.getByTestId('wm-in-name')).toHaveValue('"my wm"');
  await expect(page.getByTestId('wm-command')).toContainText(`-in '"my wm"'`);
  await expect(page.getByTestId('wm-in-refused')).toHaveCount(0);
  // A name typed and then cleared is none: the files' database names -in again, shown as the field's placeholder (WW4d
  // R2, N1).
  await page.getByTestId('wm-in-name').fill('typed/db');
  await expect(page.getByTestId('wm-command')).toContainText('-in typed/db');
  await page.getByTestId('wm-in-name').fill('');
  await expect(page.getByTestId('wm-in-name')).toHaveValue('');
  await expect(page.getByTestId('wm-in-name')).toHaveAttribute('placeholder', '"my wm"');
  await expect(page.getByTestId('wm-command')).toContainText(`-in '"my wm"'`);
  // The one-shot error of an add (a file whose name is in the run already) is an alert.
  await page.getByTestId('wm-in-files').setInputFiles([named('my wm.nin')]);
  await expect(page.getByTestId('wm-in-error')).toHaveText('A file named "my wm.nin" is already in the run directory.');
  await expect(page.getByTestId('wm-in-error')).toHaveAttribute('role', 'alert');

  // makeblastdb's -in of -input_type blastdb: the Input block, the check and Run LOSAT (WW4d R1).
  await page.getByTestId('tab-makeblastdb').click();
  await described(page, 'mk-panel', 'ready');
  await page.getByTestId('mk-param-input_type').selectOption('blastdb');
  await page.getByTestId('mk-db-files').setInputFiles([named('q"mk.nin'), named('q"mk.nsq')]);
  const mkRefused = `A database named "q"mk" cannot be in -in: ${REASON}.`;
  await expect(page.getByTestId('mk-db-in-refused')).toHaveText(mkRefused);
  await expect(page.getByTestId('mk-db-in-refused')).not.toHaveAttribute('role', 'alert');
  await expect(page.getByTestId('mk-db-in')).toHaveValue('');
  await checked(page, 'mk-panel', 'invalid');
  await expect(page.getByTestId('mk-argv-validation')).toHaveText(mkRefused);
  await shots(page, 'mk-refused');
  await page.getByTestId('mk-run-losat').click();
  await expect(page.getByTestId('mk-message')).toHaveText(mkRefused);
  await page.getByTestId('mk-db-in').fill('mydb');
  await expect(page.getByTestId('mk-db-in-refused')).toHaveCount(0);

  // The Masking tab's dustmasker -in of -infmt blastdb (WW4c), a SeqDB name list as WindowMasker's (WW4d R1 probe).
  await page.getByTestId('tab-masking').click();
  await described(page, 'ms-panel', 'ready');
  await page.getByTestId('ms-param-infmt').fill('blastdb');
  await page.getByTestId('ms-in-files').setInputFiles([named('q"ms.nin'), named('q"ms.nsq')]);
  const msRefused = `A database named "q"ms" cannot be in -in: ${REASON}.`;
  await expect(page.getByTestId('ms-in-refused')).toHaveText(msRefused);
  await expect(page.getByTestId('ms-in-refused')).not.toHaveAttribute('role', 'alert');
  await expect(page.getByTestId('ms-in-name')).toHaveValue('');
  await checked(page, 'ms-panel', 'invalid');
  await expect(page.getByTestId('ms-argv-validation')).toHaveText(msRefused);
  await shots(page, 'ms-refused');
  await page.getByTestId('ms-run-losat').click();
  await expect(page.getByTestId('ms-message')).toHaveText(msRefused);
  const msToggle = page.getByTestId('ms-options-toggle');
  if ((await msToggle.getAttribute('aria-expanded')) !== 'true') await msToggle.click();
  const msRunFiles = page.locator('[data-testid="ms-run-files-list"] > li');
  await expect(msRunFiles).toHaveCount(2);
  while ((await msRunFiles.count()) > 0) await page.getByTestId('ms-run-file-remove-0').click();
  await page.getByTestId('ms-in-files').setInputFiles([named('my ms.nin'), named('my ms.nsq')]);
  await expect(page.getByTestId('ms-in-name')).toHaveValue('"my ms"');
  await expect(page.getByTestId('ms-command')).toContainText(`-in '"my ms"'`);
  await expect(page.getByTestId('ms-in-refused')).toHaveCount(0);
});

test.describe('real engine: database names with a space are the native command line\'s (line WW, WW4d)', () => {
  test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
  test.skip(NATIVE === undefined, NO_NATIVE_REASON);

  test('blastn -db on a folder "my dbs" and on an alias "my nslice", and BLASTN\'s Repeats database in a folder "rep dbs"', async ({ page }) => {
    const work = mkdtempSync(join(tmpdir(), 'losat-ww4d-'));
    try {
      const dbSearch = join(LOSAT, 'tests/db_search');
      const nslice = filesOf(dbSearch, 'nslice');
      copyInto(dbSearch, nslice, join(work, 'my dbs'));
      const queryText = readFileSync(join(dbSearch, 'nslice_q.fna'), 'utf8');
      writeFileSync(join(work, 'q.fna'), queryText);

      // 1. A picked folder "my dbs": -db "my dbs/nslice".
      await openSearch(page);
      await program(page, 'blastn');
      await paste(page, 'query', queryText);
      await page.getByTestId('subject-source-database').check();
      await page.getByTestId('subject-db-folder-input').setInputFiles(join(work, 'my dbs'));
      await expect(page.getByTestId('subject-db-selected')).toContainText('my dbs/nslice');
      await submit(page);
      await waitStatus(page, 1, 'completed');
      for (const format of [0, 6, 7] as const) {
        const shown = await result(page, 1, format);
        expect(shown.command).toContain(`-db '"my dbs/nslice"'`);
        expect(shown.output, `blastn -db "my dbs/nslice" outfmt ${format}`).toBe(
          nativeSearch(['blastn', '-query', 'q.fna', '-db', '"my dbs/nslice"'], work, format),
        );
      }

      // 2. An alias whose name holds a space, opened as files: -db "my nslice".
      const aliasDir = join(work, 'alias');
      copyInto(dbSearch, nslice, aliasDir);
      writeFileSync(join(aliasDir, 'my nslice.nal'), 'TITLE my nslice\nDBLIST nslice\n');
      await page.getByTestId('subject-db-files-input').setInputFiles([join(aliasDir, 'my nslice.nal'), ...nslice.map((name) => join(aliasDir, name))]);
      await page.getByTestId('subject-db-select').selectOption('my nslice');
      await submit(page);
      await waitStatus(page, 2, 'completed');
      for (const format of [0, 6] as const) {
        const shown = await result(page, 2, format);
        expect(shown.command).toContain(`-db '"my nslice"'`);
        expect(shown.output, `blastn -db "my nslice" outfmt ${format}`).toBe(
          nativeSearch(['blastn', '-query', join(work, 'q.fna'), '-db', '"my nslice"'], aliasDir, format),
        );
      }

      // 3. BLASTN's Repeats database in a folder "rep dbs": -filtering_db "rep dbs/rpt_d4", with -subject.
      copyInto(dbSearch, filesOf(dbSearch, 'rpt_d4'), join(work, 'rep dbs'));
      const inputs = join(LOSAT, 'tests/fixtures/blastn_regression/inputs');
      const wmQuery = readFileSync(join(inputs, 'wm_query.fa'), 'utf8');
      const subject = readFileSync(join(inputs, 's60k.fa'), 'utf8');
      writeFileSync(join(work, 'wm_query.fa'), wmQuery);
      await paste(page, 'query', wmQuery);
      await page.getByTestId('subject-source-fasta').check();
      await paste(page, 'subject', subject);
      await openParameters(page);
      await page.getByTestId('param-filtering_db').check();
      await page.getByTestId('param-filtering_db-folder-input').setInputFiles(join(work, 'rep dbs'));
      await expect(page.getByTestId('param-filtering_db-selected')).toContainText('rep dbs/rpt_d4');
      await page.getByTestId('param-soft_masking').uncheck();
      await submit(page);
      await waitStatus(page, 3, 'completed');
      for (const format of [0, 6] as const) {
        const shown = await result(page, 3, format);
        expect(shown.command).toContain(`-filtering_db '"rep dbs/rpt_d4"'`);
        // The report names the subject's file: the native run reads it under the run's own name.
        const subjectName = / -subject (\S+)/.exec(shown.command)![1]!;
        writeFileSync(join(work, subjectName), subject);
        expect(shown.output, `blastn -filtering_db "rep dbs/rpt_d4" outfmt ${format}`).toBe(
          nativeSearch(['blastn', '-query', 'wm_query.fa', '-subject', subjectName, '-filtering_db', '"rep dbs/rpt_d4"', '-soft_masking', 'false'], work, format),
        );
      }

      // 4. Help with a database that -db cannot name (WW4d R2, L1): the words print the command line's help, which reads
      // no database, as `losat blastn -db 'q"d/nslice' -help` does.
      await page.getByTestId('param-filtering_db').uncheck();
      copyInto(dbSearch, nslice, join(work, 'q"d'));
      await page.getByTestId('subject-source-database').check();
      await page.getByTestId('subject-db-folder-input').setInputFiles(join(work, 'q"d'));
      const dbRefused = `A database named "q"d/nslice" cannot be in -db: ${REASON}.`;
      await expect(page.getByTestId('subject-db-refused')).toHaveText(dbRefused);
      await expect(page.getByTestId('argv-validation')).toHaveText(dbRefused);
      await page.getByTestId('param-extra').fill('-help');
      await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
      await submit(page);
      await waitStatus(page, 4, 'completed');
      await pickRun(page, 4);
      await page.getByTestId('results-view-outputs').click();
      const [printed] = await Promise.all([page.waitForEvent('download'), page.getByTestId('export-printed').click()]);
      const helpNative = native(['blastn', '-db', 'q"d/nslice', '-help'], work);
      expect(helpNative.status).toBe(0);
      expect(readFileSync((await printed.path())!).toString('utf8'), 'blastn -db \'q"d/nslice\' -help').toBe(helpNative.stdout);
    } finally {
      rmSync(work, { recursive: true, force: true });
    }
  });

  test('WindowMasker and dustmasker -infmt blastdb -in "my db/wm", and makeblastdb -out "my db" and -input_type blastdb -in "my dbs/ndb"', async ({ page }) => {
    const work = mkdtempSync(join(tmpdir(), 'losat-ww4d-'));
    try {
      // 1. WindowMasker -mk_counts of a database in a folder "my db".
      const wmDb = join(LOSAT, 'tests/windowmasker/db');
      copyInto(wmDb, filesOf(wmDb, 'wm'), join(work, 'my db'));
      await openSearch(page);
      await page.getByTestId('tab-windowmasker').click();
      await described(page, 'wm-panel', 'ready');
      await page.getByTestId('wm-param-infmt').selectOption('blastdb');
      await page.getByTestId('wm-in-folder').setInputFiles(join(work, 'my db'));
      await expect(page.getByTestId('wm-in-name')).toHaveValue('"my db/wm"');
      await page.getByTestId('wm-param-sformat').selectOption('binary');
      await page.getByTestId('wm-out').fill('-');
      await expect(page.getByTestId('wm-command')).toContainText(`-in '"my db/wm"'`);
      await checked(page, 'wm-panel', 'ok');
      await page.getByTestId('wm-run-losat').click();
      const counts = await toolRun(page, 1, 'completed');
      const countsNative = native(['windowmasker', '-mk_counts', '-infmt', 'blastdb', '-in', '"my db/wm"', '-sformat', 'binary'], work);
      expect(countsNative.status).toBe(0);
      expect(counts.stdout, 'windowmasker -in "my db/wm" stdout').toBe(countsNative.stdout);
      expect(counts.stderr, 'windowmasker -in "my db/wm" stderr').toBe(countsNative.stderr);

      // 2. makeblastdb -out "my db": the command line takes no name with a space there (NCBI's message, exit status 1).
      const nucl = readFileSync(join(LOSAT, 'tests/dropin_corpus/ns.fna'));
      await page.getByTestId('tab-makeblastdb').click();
      await described(page, 'mk-panel', 'ready');
      await page.getByTestId('mk-fasta-input').fill(nucl.toString('utf8'));
      await page.getByTestId('mk-param-dbtype-nucl').check();
      await page.getByTestId('mk-param-out').fill('my db');
      await expect(page.getByTestId('mk-command')).toHaveText(`losat makeblastdb -in input.fa -dbtype nucl -out 'my db'`);
      const outNative = nativeRun(['makeblastdb', '-in', 'input.fa', '-dbtype', 'nucl', '-out', 'my db'], new Map([['input.fa', nucl]]), ZONE);
      expect(outNative.status).toBe(1);
      await checked(page, 'mk-panel');
      await page.getByTestId('mk-run-losat').click();
      const out = await toolRun(page, 2, 'failed');
      expect(out.exit).toBe(outNative.status);
      expect(out.stderr, 'makeblastdb -out "my db" stderr').toBe(outNative.stderr.toString('utf8'));

      // 3. makeblastdb -input_type blastdb on a database in a folder "my dbs": -in "my dbs/ndb", NCBI's result.
      const corpus = join(LOSAT, 'tests/dropin_corpus');
      const ndb = filesOf(corpus, 'ndb');
      await page.getByTestId('tab-makeblastdb').click();
      await page.getByTestId('mk-param-out').fill('copy');
      await page.getByTestId('mk-param-input_type').selectOption('blastdb');
      copyInto(corpus, ndb, join(work, 'my dbs'));
      await page.getByTestId('mk-db-folder').setInputFiles(join(work, 'my dbs'));
      await expect(page.getByTestId('mk-db-in')).toHaveValue('"my dbs/ndb"');
      await expect(page.getByTestId('mk-command')).toHaveText(`losat makeblastdb -in '"my dbs/ndb"' -input_type blastdb -dbtype nucl -out copy`);
      const inputs = new Map(ndb.map((name) => [`my dbs/${name}`, readFileSync(join(corpus, name))] as const));
      const dbNative = nativeRun(['makeblastdb', '-in', '"my dbs/ndb"', '-input_type', 'blastdb', '-dbtype', 'nucl', '-out', 'copy'], inputs, ZONE);
      await checked(page, 'mk-panel');
      await page.getByTestId('mk-run-losat').click();
      const db = await toolRun(page, 3, dbNative.status === 0 ? 'completed' : 'failed');
      expect(db.exit).toBe(dbNative.status);
      expect(db.stdout, 'makeblastdb -input_type blastdb -in "my dbs/ndb" stdout').toBe(dbNative.stdout.toString('utf8'));
      expect(db.stderr, 'makeblastdb -input_type blastdb -in "my dbs/ndb" stderr').toBe(dbNative.stderr.toString('utf8'));

      // 4. The Masking tab's dustmasker -infmt blastdb on the database in the folder "my db" (WW4c; WW4d R1).
      await page.getByTestId('tab-masking').click();
      await described(page, 'ms-panel', 'ready');
      await page.getByTestId('ms-param-infmt').fill('blastdb');
      await page.getByTestId('ms-in-folder').setInputFiles(join(work, 'my db'));
      await expect(page.getByTestId('ms-in-name')).toHaveValue('"my db/wm"');
      await expect(page.getByTestId('ms-command')).toContainText(`-in '"my db/wm"'`);
      await checked(page, 'ms-panel', 'ok');
      await page.getByTestId('ms-run-losat').click();
      const dust = await toolRun(page, 4, 'completed');
      const dustNative = native(['dustmasker', '-in', '"my db/wm"', '-infmt', 'blastdb'], work);
      expect(dustNative.status).toBe(0);
      expect(dust.stdout, 'dustmasker -in "my db/wm" stdout').toBe(dustNative.stdout);
      expect(dust.stderr, 'dustmasker -in "my db/wm" stderr').toBe(dustNative.stderr);
    } finally {
      rmSync(work, { recursive: true, force: true });
    }
  });
});
