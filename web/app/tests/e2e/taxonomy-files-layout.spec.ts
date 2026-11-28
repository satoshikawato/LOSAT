// The run's directory is the user's directory (line WQA D35, 2026-11-28): with a database folder picked, the
// folder's `taxdb.*` and `taxonomy4blast.sqlite3` are where the user has them (`protdb/`), the engine finds them
// through the shown command's `BLASTDB=protdb`, and no copy of them lies at the run's root. So a `-logfile` named
// like one of them that the run does not look for is a new file there, and runs as the shown command runs on the
// command line, with its bytes; one that the run looks for (`taxdb.bti` in a search of a database, the run's
// working directory first; `taxonomy4blast.sqlite3` with `-taxids`) is refused, as the command line refuses it
// (WQA DECISIONS 62).
import { spawnSync } from 'node:child_process';
import { copyFileSync, existsSync, mkdirSync, mkdtempSync, readFileSync, readdirSync, rmSync, statSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { expect } from '@playwright/test';
import { ENGINE, NO_ENGINE_REASON } from './support/engine';
import { REPOSITORY } from './support/harness-server';
import { NATIVE, NO_NATIVE_REASON } from './support/native';
import { openParameters, openSearch, paste, program, result, runMessage, submit, waitStatus } from './support/search';
import { test } from './support/test';

const LOSAT = join(REPOSITORY, 'LOSAT');
const RUNNING_NAMES = ['taxdb.btd', 'taxonomy4blast.sqlite3'] as const;

/** `work/protdb/`: the protein database `tfp` with its taxonomy files, as a user's folder holds them. */
function userFolder(work: string): string {
  const folder = join(work, 'protdb');
  mkdirSync(folder, { recursive: true });
  const source = join(LOSAT, 'tests/db_search/taxfilter');
  for (const name of readdirSync(source)) {
    if (statSync(join(source, name)).isFile() && (name.startsWith('tfp.') || name === 'taxonomy4blast.sqlite3')) copyFileSync(join(source, name), join(folder, name));
  }
  copyFileSync(join(LOSAT, 'tests/seqdb/taxdb/good/taxdb.bti'), join(folder, 'taxdb.bti'));
  copyFileSync(join(LOSAT, 'tests/seqdb/taxdb/good/taxdb.btd'), join(folder, 'taxdb.btd'));
  return folder;
}

function native(args: readonly string[], cwd: string, format: number | string): { status: number | null; stdout: string; stderr: string } {
  const proc = spawnSync(NATIVE!, [...args, '-outfmt', String(format), '-num_threads', '1'], {
    cwd,
    maxBuffer: 1 << 28,
    env: { ...process.env, BLASTDB: 'protdb' },
  });
  return { status: proc.status, stdout: proc.stdout.toString('utf8'), stderr: proc.stderr.toString('utf8') };
}

test.describe('the folder picked for -db: its taxonomy files stay where the user has them', () => {
  test.skip(ENGINE === undefined, NO_ENGINE_REASON);
  test.skip(NATIVE === undefined, NO_NATIVE_REASON);
  test.setTimeout(300_000);

  for (const name of RUNNING_NAMES) {
    test(`-logfile ${name} runs and gives the bytes of the shown command`, async ({ page }) => {
      const work = mkdtempSync(join(tmpdir(), 'losat-d35-e2e-'));
      try {
        const folder = userFolder(work);
        const queryPath = join(LOSAT, 'tests/db_search/taxfilter/q.faa');

        await openSearch(page);
        await program(page, 'blastp');
        await paste(page, 'query', readFileSync(queryPath, 'utf8'));
        await page.getByTestId('subject-source-database').check();
        await page.getByTestId('subject-db-folder-input').setInputFiles(folder);
        await expect(page.getByTestId('subject-db-selected')).toContainText('protdb/tfp');
        await openParameters(page);
        await page.getByTestId('param-extra').fill(`-logfile ${name}`);
        await submit(page);
        await waitStatus(page, 1, 'completed');

        const { command, output } = await result(page, 1, 6);
        expect(command.startsWith('BLASTDB=protdb losat blastp '), command).toBe(true);
        expect(command).toContain(`-logfile ${name}`);
        // The shown command, run from the folder's parent as shown: the same bytes, and the CLI makes no complaint.
        const expected = native(['blastp', '-query', queryPath, '-db', 'protdb/tfp', '-logfile', name], work, 6);
        expect(expected.status, expected.stderr).toBe(0);
        expect(output).toBe(expected.stdout);
      } finally {
        rmSync(work, { recursive: true, force: true });
      }
    });
  }

  for (const [label, extra, name] of [
    ['taxdb.bti, which the search looks for first in its working directory', '', 'taxdb.bti'],
    ['taxonomy4blast.sqlite3 with -taxids', '-taxids 9606 ', 'taxonomy4blast.sqlite3'],
  ] as const) {
    test(`-logfile ${label} is refused as the shown command is`, async ({ page }) => {
      const work = mkdtempSync(join(tmpdir(), 'losat-d35-e2e-'));
      try {
        const folder = userFolder(work);
        const queryPath = join(LOSAT, 'tests/db_search/taxfilter/q.faa');
        await openSearch(page);
        await program(page, 'blastp');
        await paste(page, 'query', readFileSync(queryPath, 'utf8'));
        await page.getByTestId('subject-source-database').check();
        await page.getByTestId('subject-db-folder-input').setInputFiles(folder);
        await expect(page.getByTestId('subject-db-selected')).toContainText('protdb/tfp');
        await openParameters(page);
        await page.getByTestId('param-extra').fill(`${extra}-logfile ${name}`);
        await submit(page);
        await waitStatus(page, 1, 'failed');
        expect(await runMessage(page, 1)).toContain(`a -logfile file '${name}' that is`);
        // The shown command on the command line: refused before the file is made.
        const args = ['blastp', '-query', queryPath, '-db', 'protdb/tfp', ...extra.trim().split(' ').filter(Boolean), '-logfile', name];
        const cli = native(args, work, 6);
        expect(cli.status, cli.stderr).toBe(1);
        expect(cli.stderr).toContain(`a -logfile file '${name}' that is also an input file read by this run`);
        expect(existsSync(join(work, name))).toBe(false);
      } finally {
        rmSync(work, { recursive: true, force: true });
      }
    });
  }

  test('a -logfile that names the folder\'s own taxonomy file is refused, as the file is the user\'s', async ({ page }) => {
    const work = mkdtempSync(join(tmpdir(), 'losat-d35-e2e-'));
    try {
      const folder = userFolder(work);
      await openSearch(page);
      await program(page, 'blastp');
      await paste(page, 'query', readFileSync(join(LOSAT, 'tests/db_search/taxfilter/q.faa'), 'utf8'));
      await page.getByTestId('subject-source-database').check();
      await page.getByTestId('subject-db-folder-input').setInputFiles(folder);
      await expect(page.getByTestId('subject-db-selected')).toContainText('protdb/tfp');
      await openParameters(page);
      await page.getByTestId('param-extra').fill('-logfile protdb/taxdb.bti');
      await submit(page);
      await waitStatus(page, 1, 'failed');
      expect(await runMessage(page, 1)).toContain("a -logfile file 'protdb/taxdb.bti' that is a file of the run's directory");
    } finally {
      rmSync(work, { recursive: true, force: true });
    }
  });
});
