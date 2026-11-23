// A run of -db through the package's host in Node (line WQA D12a, 2026-11-20; -filtering_db, line RG item D4): the serial reactor of
// LOSAT_WEB_REACTORS (web/adapter/tools/build_reactors.py's output), the run's directory (RunDirectory) holding the
// files that a folder pick of the database gives (every file of the database's folder) and the list of a filter, and
// a reader that reads each Blob's bytes from its file as FileReaderSync reads a File in a worker. The engine opens the
// files by the names of the argv and reads them in place; outfmt 6 is NCBI BLAST+ 2.17.0's frozen bytes for the same
// search (LOSAT/tests/fixtures/db_search) in LOSAT's words (LOSAT/tests/losat_wording.py map-output, line NW2). Skips
// without LOSAT_WEB_REACTORS.
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { closeSync, openAsBlob, openSync, readdirSync, readFileSync, readSync, statSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import type { DirectoryFile } from '../src/domain/run-files.ts';
import { DISCARD_STREAM, ROLE_QUERY, ROLE_SUBJECT } from '../src/infra/reactor/abi.ts';
import { compileReactor, instantiateSerial } from '../src/infra/reactor/instance.ts';
import { RunDirectory, type BlobReader } from '../src/infra/reactor/run-directory.ts';

const LOSAT = fileURLToPath(new URL('../../../../LOSAT', import.meta.url));
const REACTORS = process.env['LOSAT_WEB_REACTORS'];

interface Row {
  readonly id: string;
  readonly program: string;
  readonly db: string;
  readonly query: string;
  readonly extra: readonly string[];
}

interface RegressionRow {
  readonly id: string;
  readonly argv: readonly string[];
  readonly stdoutSha256: string;
}

/** A frozen NCBI output file in LOSAT's words, as LOSAT writes it. */
function frozenFileInLosatWords(fixturePath: string, label: string): Buffer {
  const mapped = spawnSync('python3', ['-I', join(LOSAT, 'tests/losat_wording.py'), 'map-output', fixturePath]);
  if (mapped.error) throw mapped.error;
  if (mapped.status !== 0) throw new Error(`losat_wording.py map-output ${label}: ${mapped.stderr.toString()}`);
  return mapped.stdout;
}

/** A row's frozen NCBI outfmt 6 (LOSAT/tests/fixtures/db_search) in LOSAT's words, as LOSAT writes it. */
function frozenInLosatWords(id: string): Buffer {
  return frozenFileInLosatWords(join(LOSAT, 'tests/fixtures/db_search', `${id}.out`), id);
}

/** A row of LOSAT/tests/fixtures/blastn_regression/manifest.tsv by its case_id. */
function blastnRegressionRow(caseId: string): RegressionRow {
  const lines = readFileSync(join(LOSAT, 'tests/fixtures/blastn_regression/manifest.tsv'), 'utf8')
    .split('\n')
    .filter((line) => line !== '' && !line.startsWith('#'));
  const header = lines[0]!.split('\t');
  for (const line of lines.slice(1)) {
    const cells = Object.fromEntries(line.split('\t').map((cell, i) => [header[i], cell])) as Record<string, string>;
    if (cells['case_id'] === caseId) {
      return {
        id: caseId,
        argv: cells['argv']!.split(' '),
        stdoutSha256: cells['stdout_sha256']!,
      };
    }
  }
  throw new Error(`blastn_regression manifest row not found: ${caseId}`);
}

/** The rows of LOSAT/tests/db_search_manifest by their ID. */
function manifestRows(ids: readonly string[]): Row[] {
  const rows = new Map<string, Row>();
  const directory = join(LOSAT, 'tests/db_search_manifest');
  for (const name of readdirSync(directory).filter((file) => file.endsWith('.tsv'))) {
    const lines = readFileSync(join(directory, name), 'utf8')
      .split('\n')
      .filter((line) => line !== '' && !line.startsWith('#'));
    const header = lines[0]!.split('\t');
    for (const line of lines.slice(1)) {
      const cells = Object.fromEntries(line.split('\t').map((cell, i) => [header[i], cell])) as Record<string, string>;
      rows.set(cells['row_id']!, {
        id: cells['row_id']!,
        program: cells['program']!,
        db: cells['db']!,
        query: cells['query']!,
        extra: cells['extra_args'] === '' ? [] : cells['extra_args']!.split(' '),
      });
    }
  }
  return ids.map((id) => rows.get(id)!);
}

/** The files of a folder, named by their paths under LOSAT/ (the CLI's directory), each a Blob of its file. */
async function folder(relative: string, paths: Map<Blob, string>): Promise<DirectoryFile[]> {
  const files: DirectoryFile[] = [];
  for (const name of readdirSync(join(LOSAT, relative))) {
    const path = join(LOSAT, relative, name);
    if (!statSync(path).isFile()) continue;
    const file = await openAsBlob(path);
    paths.set(file, path);
    files.push({ name: `${relative}/${name}`, file });
  }
  return files;
}

/** Reads a Blob's bytes from its file, synchronously; the bytes read in all, and by file. */
function osReader(paths: Map<Blob, string>): { read: BlobReader; bytes: () => number; byFile: Map<string, number> } {
  let total = 0;
  const byFile = new Map<string, number>();
  return {
    bytes: () => total,
    byFile,
    read: (blob, start, end) => {
      const bytes = new Uint8Array(end - start);
      const path = paths.get(blob)!;
      const descriptor = openSync(path, 'r');
      try {
        const got = readSync(descriptor, bytes, 0, bytes.length, start);
        total += got;
        byFile.set(path, (byFile.get(path) ?? 0) + got);
        return bytes.subarray(0, got);
      } finally {
        closeSync(descriptor);
      }
    },
  };
}

describe.skipIf(REACTORS === undefined)('a run of -db through the host', () => {
  it("reads the database's files in place, and writes NCBI's frozen outfmt 6", async () => {
    const { module } = await compileReactor(readFileSync(join(REACTORS!, 'losat-web-serial.wasm')));
    const paths = new Map<Blob, string>();
    const reader = osReader(paths);
    const directory = new RunDirectory(undefined, reader.read);
    const { abi } = await instantiateSerial(module, directory);
    const rows = manifestRows(['db.blastn.mega.6', 'db.d1.6', 'db.blastx.default.6', 'db.tblastn.default.6', 'db.tblastx.default.6', 'db.sl.bp.bsl.6', 'db.tf.bp.nosql.6', 'db.rpt.mb.6', 'db.rpt.bn.6', 'db.rpt.plus.6']);
    for (const row of rows) {
      const files = await folder(dirname(row.db), paths);
      const list = row.extra.indexOf('-seqidlist');
      if (list >= 0) files.push(...(await folder(dirname(row.extra[list + 1]!), paths)).filter(({ name }) => name === row.extra[list + 1]));
      // BLASTN's -filtering_db is a second database in the run's directory (line RG item D4): its files, by the
      // names that start with the option's value and a dot, join the folder's.
      const repeats = row.extra.indexOf('-filtering_db');
      const repeatsName = repeats >= 0 ? row.extra[repeats + 1]! : undefined;
      if (repeatsName !== undefined) {
        const names = new Set(files.map(({ name }) => name));
        files.push(...(await folder(dirname(repeatsName), paths)).filter(({ name }) => name.startsWith(`${repeatsName}.`) && !names.has(name)));
      }
      directory.open(files);
      const argv = [row.program, '-db', row.db, '-query', row.query, ...row.extra, '-num_threads', '1'];
      // The query is read with the data loaders that the verdict on the argv gives (ABI v2 §9).
      const validation = abi.validate(argv);
      expect(validation, row.id).toEqual({ ok: true, dataLoaders: true });
      const started = abi.start(argv);
      expect(started.ok && started.start.database, row.id).toBe(true);
      const query = abi.register(row.program, ROLE_QUERY, validation.ok && validation.dataLoaders, readFileSync(join(LOSAT, row.query)));
      const before = reader.bytes();
      reader.byFile.clear();
      let outfmt6: Uint8Array[] = [];
      abi.run(argv, query.handle, 0, (stream, bytes) => {
        if (stream === 6) outfmt6.push(bytes.slice());
        else if (stream === DISCARD_STREAM) outfmt6 = [];
      });
      abi.release(query.handle);
      expect(Buffer.concat(outfmt6).equals(frozenInLosatWords(row.id)), row.id).toBe(true);
      // Only what the search reads: the database's own files (`<db>.*`) and the filter's list, by range, never
      // another file of the folder (2.8 MB in all). A file may be read more than once in a run (its index at each
      // open of the database, its sequences by the search and again by the traceback and the report), as the CLI
      // reads it (line WQA D12a review, low 5).
      const own = new Set(files.filter(({ name }) => name.startsWith(`${row.db}.`) || (repeatsName !== undefined && name.startsWith(`${repeatsName}.`)) || row.extra.includes(name)).map(({ name }) => join(LOSAT, name)));
      const read = [...reader.byFile].map(([path, bytes]) => `${path.slice(LOSAT.length + 1)} ${bytes}/${statSync(path).size}`);
      console.log(`${row.id}: read ${reader.bytes() - before} bytes: ${read.join(', ')}`);
      expect([...reader.byFile.keys()].filter((path) => !own.has(path)), row.id).toEqual([]);
    }
    // A -subject search with BLASTN -filtering_db reads the repeat database from the run's directory while taking
    // the query and subject through registered handles (line RG item D4).
    const repeatSubject = blastnRegressionRow('rpt.megablast');
    const repeatQueryPath = repeatSubject.argv[repeatSubject.argv.indexOf('-query') + 1]!;
    const repeatSubjectPath = repeatSubject.argv[repeatSubject.argv.indexOf('-subject') + 1]!;
    const repeatDbName = repeatSubject.argv[repeatSubject.argv.indexOf('-filtering_db') + 1]!;
    const repeatFiles = await folder(dirname(repeatDbName), paths);
    directory.open(repeatFiles);
    const repeatArgv = ['blastn', ...repeatSubject.argv, '-num_threads', '1'];
    const repeatValidation = abi.validate(repeatArgv);
    expect(repeatValidation, repeatSubject.id).toEqual({ ok: true, dataLoaders: true });
    const repeatStarted = abi.start(repeatArgv);
    expect(repeatStarted.ok && !repeatStarted.start.database, repeatSubject.id).toBe(true);
    const repeatQuery = abi.register('blastn', ROLE_QUERY, true, readFileSync(join(LOSAT, repeatQueryPath)));
    const repeatTarget = abi.register('blastn', ROLE_SUBJECT, true, readFileSync(join(LOSAT, repeatSubjectPath)));
    const repeatBefore = reader.bytes();
    reader.byFile.clear();
    let repeatOutfmt6: Uint8Array[] = [];
    abi.run(repeatArgv, repeatQuery.handle, repeatTarget.handle, (stream, bytes) => {
      if (stream === 6) repeatOutfmt6.push(bytes.slice());
      else if (stream === DISCARD_STREAM) repeatOutfmt6 = [];
    });
    abi.release(repeatQuery.handle);
    abi.release(repeatTarget.handle);
    const repeatOutput = Buffer.concat(repeatOutfmt6);
    const repeatFrozen = frozenFileInLosatWords(join(LOSAT, 'tests/fixtures/blastn_regression', `${repeatSubject.id}.out`), repeatSubject.id);
    expect(repeatOutput.equals(repeatFrozen), repeatSubject.id).toBe(true);
    expect(createHash('sha256').update(repeatOutput).digest('hex'), repeatSubject.id).toBe(repeatSubject.stdoutSha256);
    const repeatOwn = new Set(repeatFiles.filter(({ name }) => name.startsWith(`${repeatDbName}.`)).map(({ name }) => join(LOSAT, name)));
    const repeatRead = [...reader.byFile].map(([path, bytes]) => `${path.slice(LOSAT.length + 1)} ${bytes}/${statSync(path).size}`);
    console.log(`${repeatSubject.id}: read ${reader.bytes() - repeatBefore} bytes: ${repeatRead.join(', ')}`);
    expect(reader.byFile.size, repeatSubject.id).toBeGreaterThan(0);
    expect([...reader.byFile.keys()].filter((path) => !repeatOwn.has(path)), repeatSubject.id).toEqual([]);
    // A repeat database whose files are omitted from the directory: the CLI's batch warning, with `/work::` as its
    // search path (abi_v2.md §3; line RG item D4).
    directory.open([]);
    const nucQuery = abi.register('blastn', ROLE_QUERY, true, new TextEncoder().encode('>q\nACGTACGTACGTACGTACGTACGTACGTACGT\n'));
    const nucSubject = abi.register('blastn', ROLE_SUBJECT, true, new TextEncoder().encode('>s\nACGTACGTACGTACGTACGTACGTACGTACGT\n'));
    const missingRepeats = ['blastn', '-query', 'q.fna', '-subject', 's.fna', '-filtering_db', repeatDbName, '-num_threads', '1'];
    expect(() => abi.run(missingRepeats, nucQuery.handle, nucSubject.handle, () => undefined)).toThrow(
      `LOSAT engine error: Warning: No alias or index file found for nucleotide database [${repeatDbName}] in search path [/work::].`,
    );
    abi.release(nucQuery.handle);
    abi.release(nucSubject.handle);
    // A database that the directory does not hold: the CLI's message.
    directory.open([]);
    const argv = ['blastp', '-db', 'tests/db_search/paja', '-query', 'q.faa', '-num_threads', '1'];
    const query = abi.register('blastp', ROLE_QUERY, true, new TextEncoder().encode('>q\nMKVLA\n'));
    expect(() => abi.run(argv, query.handle, 0, () => undefined)).toThrow(/LOSAT database error: No alias or index file found for protein database \[tests\/db_search\/paja\]/);
    // With the database in the directory, the names that the CLI run from LOSAT/ does not find either: an absolute
    // name (the CLI's root holds no tests/) and `..` after a name that is no folder. The CLI's message, whose
    // search path names the engine's current directory, /work (abi_v2.md §3; review of D12a, low 4).
    directory.open(await folder('tests/db_search', paths));
    for (const name of ['/tests/db_search/paja', 'nofolder/../tests/db_search/paja']) {
      const named = ['blastp', '-db', name, '-query', 'q.faa', '-num_threads', '1'];
      expect(() => abi.run(named, query.handle, 0, () => undefined), name).toThrow(
        `LOSAT database error: No alias or index file found for protein database [${name}] in search path [/work::]. ` +
          'Please verify the spelling of the BLAST database and its molecule type.',
      );
    }
    abi.release(query.handle);
  }, 120_000);
});
