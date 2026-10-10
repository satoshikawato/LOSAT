// Session files (plan §5.8, design §12.2, S15 item 4; the format: docs/web/session_file.md). A
// session file holds the completed runs of a working session: what each run searched (its
// RunSnapshot without the bytes of its inputs), what happened (its RunRecord), its outputs 0/6/7,
// HSP records and diagnostics byte for byte, and the identity of its inputs (the record table, the
// SHA-256 of the engine input and of each record, the sources and their exclusions); and, when
// chosen, the candidate tray's candidates and notes (S15 decision 1). It never holds the input
// FASTA, storage paths, tokens, or the IDs of runs, revisions or sources: references inside the
// file are positions in it.
//
// The container, before gzip: a first line naming it and its version, then blocks
// "<name> <length>\n" + <length bytes> + "\n" in a fixed order (the manifest, five per run, the
// candidates when included), then an end line. This module writes the lines, checks a manifest
// and a candidates block against the limits (the same checks on saving and on loading, so that a
// saved file always loads), and reads a container incrementally: the run blocks pass through in
// pieces and are never held; only the manifest and the candidates are collected, in buffers that
// the limits bound.
import type { FastaParserKind } from './dataset';
import { indexParser, PROGRAMS, type InputRole, type ProgramId } from './programs';
import type { RunRecord } from './run';

export const SESSION_CONTAINER = 'LOSAT-WEB-SESSION';
export const SESSION_CONTAINER_VERSION = 1;
export const SESSION_END = 'LOSAT-WEB-SESSION-END';
export const SESSION_FORMAT = 'losat-web-session';
export const SESSION_SCHEMA = 1;
export const SESSION_MIME = 'application/gzip';
export const MANIFEST_BLOCK = 'manifest';
export const CANDIDATES_BLOCK = 'candidates';

/** The blocks of a run, in the file's order: outfmt 0, 6 and 7, the HSP records (ABI stream 1) and the diagnostics (stream 3). */
export const SESSION_STREAMS = Object.freeze(['out0', 'out6', 'out7', 'hits', 'diagnostics'] as const);
export type SessionStream = (typeof SESSION_STREAMS)[number];

const MiB = 1024 * 1024;

/** What a session file may hold; a file over any limit is refused, and a session over one is not saved. */
export interface SessionLimits {
  /** Bytes of the manifest block and of the candidates block. */
  readonly manifestBytes: number;
  readonly candidatesBytes: number;
  /** Bytes of the first line, a block header line and the end line, without the line end. */
  readonly lineBytes: number;
  readonly runs: number;
  readonly argvWords: number;
  readonly argvWordChars: number;
  readonly titleChars: number;
  /** File names and the -query / -subject names. */
  readonly nameChars: number;
  /** Record IDs. */
  readonly idChars: number;
  readonly noteChars: number;
  /** The app's version and build, an engine build, a fallback reason. */
  readonly textChars: number;
  /** Sources of one input (files joined into one search input). */
  readonly sources: number;
  readonly candidates: number;
}

export const SESSION_LIMITS: SessionLimits = Object.freeze({
  manifestBytes: 64 * MiB,
  candidatesBytes: 64 * MiB,
  lineBytes: 64,
  runs: 1000,
  argvWords: 1000,
  argvWordChars: 10_000,
  titleChars: 10_000,
  nameChars: 1000,
  idChars: 10_000,
  noteChars: 100_000,
  textChars: 10_000,
  sources: 10_000,
  candidates: 1_000_000,
});

// --- the manifest and the candidates -------------------------------------------------------------

export interface SessionApp {
  readonly version: string;
  /** The git commit of the build (short SHA), or "unknown". */
  readonly build: string;
}

/** The record table of a run input, as columns: element `k` is the record at position `k` (`q_idx` / `s_idx`). */
export interface SessionRecordTable {
  readonly id: readonly string[];
  /** Residues, in the record's letters. */
  readonly length: readonly number[];
  /** Lower-case hex SHA-256 of the record's original bytes (`DatasetRecord.sha256`). */
  readonly sha256: readonly string[];
}

/** One source of a run input: a file (or pasted text), in the order the input joined them. */
export interface SessionSource {
  /** The file's name (`query.fa` or `subject.fa` for pasted text). */
  readonly name: string;
  /** Bytes of the file. */
  readonly size: number;
  /** Records of the file's record table. */
  readonly records: number;
  /** 0-based indices of the records that the run left out, ascending. */
  readonly excluded: readonly number[];
}

/** The identity of one input of a run: what the re-attachment of its original FASTA checks (REQ-23). */
export interface SessionInput {
  /** The name passed as -query / -subject. */
  readonly name: string;
  /** Lower-case hex SHA-256 of the engine input (`InputSnapshot.sha256`). */
  readonly sha256: string;
  /** Bytes of the engine input. */
  readonly length: number;
  /** The reader kind that the record table was made with (`indexParser(program, role)`). */
  readonly reader: FastaParserKind;
  readonly records: SessionRecordTable;
  readonly sources: readonly SessionSource[];
}

/** The RunRecord of a completed run (a completed run has no error). */
export type SessionRunRecord = Omit<RunRecord, 'error'>;

export interface SessionGroup {
  /** 1-based number of the group in the file: runs with the same number were queued together. */
  readonly index: number;
  readonly position: number;
  readonly size: number;
}

export interface SessionRun {
  /** The run's number in the working session that saved it. */
  readonly number: number;
  readonly title?: string;
  readonly program: ProgramId;
  readonly argv: readonly string[];
  readonly requestedThreads: number | 'auto';
  readonly group?: SessionGroup;
  readonly queuedAt: number;
  readonly record: SessionRunRecord;
  /** HSP records of the run (`ResultSetRef.hitCount`). */
  readonly hitCount: number;
  /** Bytes of each block of the run. */
  readonly blocks: Readonly<Record<SessionStream, number>>;
  readonly query: SessionInput;
  readonly subject: SessionInput;
}

export interface SessionManifest {
  readonly format: typeof SESSION_FORMAT;
  readonly schema: typeof SESSION_SCHEMA;
  readonly app: SessionApp;
  /** When the file was saved (ms since the epoch). */
  readonly savedAt: number;
  /** The file has a candidates block. */
  readonly candidates: boolean;
  readonly runs: readonly SessionRun[];
}

/** A candidate of the tray (an HSP and its note), in tray order. */
export interface SessionCandidate {
  /** 1-based position of the run in the file's `runs`. */
  readonly run: number;
  /** The HSP's `index` (its HSP record), and its query record and rank, which must agree with the record. */
  readonly index: number;
  readonly qIdx: number;
  readonly rank: number;
  readonly note: string;
  /** When it was added to the tray (ms since the epoch). */
  readonly addedAt: number;
}

/** A session file that is refused: the message says what is wrong and where (English, for the screen). */
export class SessionFileError extends Error {
  override readonly name = 'SessionFileError';
}

const NOT_SESSION = 'This is not a LOSAT Web session file';
const DAMAGED = 'The session file is damaged';
const INCOMPLETE = 'The session file is incomplete';

const PROGRAM_IDS: readonly ProgramId[] = PROGRAMS.map((program) => program.id);
const SHA256 = /^[0-9a-f]{64}$/;

/** Checks fields of one JSON block and names the field that fails ("runs[1].query.sha256"). */
class Fields {
  /**
   * The long arrays already frozen, which deepFreeze passes by without asking: WebKit answers
   * `Object.isFrozen` of an array of 100,000 elements in tens of milliseconds (fix round 2).
   */
  readonly frozen = new WeakSet<object>();

  constructor(
    private readonly block: string,
    readonly limits: SessionLimits,
  ) {}

  fail(where: string, problem: string): never {
    throw new SessionFileError(`The session file's ${this.block} is not valid: ${where} ${problem}.`);
  }

  /** A JSON object with the `required` fields, and no field outside `required` and `optional`. */
  object(value: unknown, where: string, required: readonly string[], optional: readonly string[] = []): Record<string, unknown> {
    if (typeof value !== 'object' || value === null || Array.isArray(value)) this.fail(where, 'is not a JSON object');
    const fields = value as Record<string, unknown>;
    const known = new Set([...required, ...optional]);
    for (const key of Object.keys(fields)) {
      if (!known.has(key)) this.fail(`${where}.${key}`, `is not a field of a session file (schema ${SESSION_SCHEMA})`);
    }
    for (const key of required) if (!Object.hasOwn(fields, key)) this.fail(`${where}.${key}`, 'is missing');
    return fields;
  }

  /** A whole number from 0 to 2^53 - 1. */
  count(value: unknown, where: string): number {
    if (typeof value !== 'number' || !Number.isSafeInteger(value) || value < 0) this.fail(where, 'is not a whole number of 0 or more');
    return value;
  }

  positive(value: unknown, where: string): number {
    const n = this.count(value, where);
    if (n === 0) this.fail(where, 'is 0, not 1 or more');
    return n;
  }

  text(value: unknown, where: string, maxChars: number, allowEmpty = true): string {
    if (typeof value !== 'string') this.fail(where, 'is not text');
    if (value.length > maxChars) this.fail(where, `is longer than ${maxChars} characters`);
    if (!allowEmpty && value.trim() === '') this.fail(where, 'is empty');
    return value;
  }

  sha256(value: unknown, where: string): string {
    if (typeof value !== 'string' || !SHA256.test(value)) this.fail(where, 'is not a SHA-256 (64 lower-case hex digits)');
    return value;
  }

  boolean(value: unknown, where: string): boolean {
    if (typeof value !== 'boolean') this.fail(where, 'is not true or false');
    return value;
  }

  array(value: unknown, where: string, maxItems: number): readonly unknown[] {
    if (!Array.isArray(value)) this.fail(where, 'is not a list');
    if (value.length > maxItems) this.fail(where, `has more than ${maxItems} items`);
    return value;
  }

  oneOf<T>(value: unknown, where: string, choices: readonly T[]): T {
    if (!choices.includes(value as T)) this.fail(where, `is not one of ${choices.map((choice) => JSON.stringify(choice)).join(', ')}`);
    return value as T;
  }
}

/** Records of a record table checked in one step of the manifest's checks. */
const RECORDS_PER_STEP = 20_000;

/** The work of a check in steps: the caller may let the page draw between two steps. */
type Steps<T> = Generator<void, T, void>;

function runSteps<T>(steps: Steps<T>): T {
  for (;;) {
    const step = steps.next();
    if (step.done === true) return step.value;
  }
}

async function runStepsPaced<T>(steps: Steps<T>, pause: () => Promise<void>): Promise<T> {
  for (;;) {
    const step = steps.next();
    if (step.done === true) return step.value;
    await pause();
  }
}

/**
 * Checks a parsed manifest and returns it with only its known fields. Throws a SessionFileError
 * that names the field: a format or schema other than this one's (a newer schema is "saved by a
 * newer LOSAT Web"), a missing or unknown field, a value of the wrong type or over the limits, an
 * argv that does not name the inputs, a reader kind other than the program's, or a record table
 * that disagrees with its sources and exclusions.
 */
export function checkManifest(value: unknown, limits: SessionLimits = SESSION_LIMITS): SessionManifest {
  return runSteps(manifestSteps(value, limits));
}

/**
 * `checkManifest` with `pause` between its steps, so that the page draws while a large manifest is
 * checked: a step checks at most RECORDS_PER_STEP records or freezes one column of a record table
 * (fix round 2: with 100,000 queries the checks held WebKit's page for about 260 ms, most of it
 * freezing the columns, about 35 ms each there).
 */
export function checkManifestPaced(value: unknown, limits: SessionLimits, pause: () => Promise<void>): Promise<SessionManifest> {
  return runStepsPaced(manifestSteps(value, limits), pause);
}

function* manifestSteps(value: unknown, limits: SessionLimits): Steps<SessionManifest> {
  const f: Fields = new Fields('manifest', limits);
  if (typeof value !== 'object' || value === null || Array.isArray(value)) f.fail('the manifest', 'is not a JSON object');
  const top = value as Record<string, unknown>;
  if (top.format !== SESSION_FORMAT) f.fail('format', `is not "${SESSION_FORMAT}"`);
  const schema = top.schema;
  if (typeof schema !== 'number' || !Number.isSafeInteger(schema)) f.fail('schema', 'is not a whole number');
  if (schema > SESSION_SCHEMA) {
    throw new SessionFileError(
      `The session file was saved by a newer LOSAT Web (its schema is ${schema}; this LOSAT Web reads schema ${SESSION_SCHEMA}). Open it with a newer LOSAT Web.`,
    );
  }
  if (schema !== SESSION_SCHEMA) f.fail('schema', `is ${schema}, not ${SESSION_SCHEMA}`);
  f.object(top, 'the manifest', ['format', 'schema', 'app', 'savedAt', 'candidates', 'runs']);
  const app = f.object(top.app, 'app', ['version', 'build']);
  const runs = f.array(top.runs, 'runs', limits.runs);
  if (runs.length === 0) f.fail('runs', 'is empty: the file holds no runs');
  const checkedApp = { version: f.text(app.version, 'app.version', limits.textChars), build: f.text(app.build, 'app.build', limits.textChars) };
  const savedAt = f.count(top.savedAt, 'savedAt');
  const candidates = f.boolean(top.candidates, 'candidates');
  const checkedRuns: SessionRun[] = [];
  for (const [k, run] of runs.entries()) checkedRuns.push(yield* checkRun(f, run, `runs[${k}]`));
  // The record tables' columns are frozen already (checkInput), so this does not walk them again.
  return deepFreeze({ format: SESSION_FORMAT, schema: SESSION_SCHEMA, app: checkedApp, savedAt, candidates, runs: checkedRuns }, f.frozen);
}

function* checkRun(f: Fields, value: unknown, where: string): Steps<SessionRun> {
  const run = f.object(
    value,
    where,
    ['number', 'program', 'argv', 'requestedThreads', 'queuedAt', 'record', 'hitCount', 'blocks', 'query', 'subject'],
    ['title', 'group'],
  );
  const program = f.oneOf(run.program, `${where}.program`, PROGRAM_IDS);
  const argv = f.array(run.argv, `${where}.argv`, f.limits.argvWords).map((word, i) => f.text(word, `${where}.argv[${i}]`, f.limits.argvWordChars));
  if (argv.length < 5 || argv[0] !== program || argv[1] !== '-query' || argv[3] !== '-subject') {
    f.fail(`${where}.argv`, `does not begin with "${program} -query <name> -subject <name>"`);
  }
  const query = yield* checkInput(f, run.query, `${where}.query`, program, 'query');
  const subject = yield* checkInput(f, run.subject, `${where}.subject`, program, 'subject');
  if (argv[2] !== query.name) f.fail(`${where}.argv[2]`, `is not the query's name "${query.name}"`);
  if (argv[4] !== subject.name) f.fail(`${where}.argv[4]`, `is not the subject's name "${subject.name}"`);
  const blocks = f.object(run.blocks, `${where}.blocks`, SESSION_STREAMS);
  return {
    number: f.positive(run.number, `${where}.number`),
    ...(run.title === undefined ? {} : { title: f.text(run.title, `${where}.title`, f.limits.titleChars, false) }),
    program,
    argv,
    requestedThreads: run.requestedThreads === 'auto' ? 'auto' : f.positive(run.requestedThreads, `${where}.requestedThreads`),
    ...(run.group === undefined ? {} : { group: checkGroup(f, run.group, `${where}.group`) }),
    queuedAt: f.count(run.queuedAt, `${where}.queuedAt`),
    record: checkRecord(f, run.record, `${where}.record`),
    hitCount: f.count(run.hitCount, `${where}.hitCount`),
    blocks: Object.fromEntries(SESSION_STREAMS.map((stream) => [stream, f.count(blocks[stream], `${where}.blocks.${stream}`)])) as Record<
      SessionStream,
      number
    >,
    query,
    subject,
  };
}

function checkGroup(f: Fields, value: unknown, where: string): SessionGroup {
  const group = f.object(value, where, ['index', 'position', 'size']);
  const position = f.positive(group.position, `${where}.position`);
  const size = f.positive(group.size, `${where}.size`);
  if (position > size) f.fail(`${where}.position`, `is ${position}, past the group's size ${size}`);
  return { index: f.positive(group.index, `${where}.index`), position, size };
}

function checkRecord(f: Fields, value: unknown, where: string): SessionRunRecord {
  const record = f.object(value, where, [], [
    'runtimePath',
    'threads',
    'fallbackReason',
    'engineBuild',
    'runtimeGeneration',
    'memory',
    'subjectRetained',
    'startedAt',
    'phaseTimes',
    'endedAt',
  ]);
  const optional = <T>(key: string, read: (value: unknown, where: string) => T): Record<string, T> =>
    record[key] === undefined ? {} : { [key]: read(record[key], `${where}.${key}`) };
  const text = (value: unknown, at: string) => f.text(value, at, f.limits.textChars);
  const memory = (value: unknown, at: string) => {
    const m = f.object(value, at, ['linearBytesBefore', 'linearBytesAfter', 'instanceRuns']);
    return {
      linearBytesBefore: f.count(m.linearBytesBefore, `${at}.linearBytesBefore`),
      linearBytesAfter: f.count(m.linearBytesAfter, `${at}.linearBytesAfter`),
      instanceRuns: f.count(m.instanceRuns, `${at}.instanceRuns`),
    };
  };
  const phaseTimes = (value: unknown, at: string) => {
    const times = f.object(value, at, [], ['preparing', 'running', 'finalizing']);
    return Object.fromEntries(
      (['preparing', 'running', 'finalizing'] as const).filter((phase) => times[phase] !== undefined).map((phase) => [phase, f.count(times[phase], `${at}.${phase}`)]),
    );
  };
  return {
    ...optional('runtimePath', (v, at) => f.oneOf(v, at, ['threaded', 'serial', 'fake'] as const)),
    ...optional('threads', (v, at) => f.positive(v, at)),
    ...optional('fallbackReason', text),
    ...optional('engineBuild', text),
    ...optional('runtimeGeneration', (v, at) => f.count(v, at)),
    ...optional('memory', memory),
    ...optional('subjectRetained', (v, at) => f.boolean(v, at)),
    ...optional('startedAt', (v, at) => f.count(v, at)),
    ...optional('phaseTimes', phaseTimes),
    ...optional('endedAt', (v, at) => f.count(v, at)),
  } as SessionRunRecord;
}

function* checkInput(f: Fields, value: unknown, where: string, program: ProgramId, role: InputRole): Steps<SessionInput> {
  const input = f.object(value, where, ['name', 'sha256', 'length', 'reader', 'records', 'sources']);
  const reader = f.oneOf(input.reader, `${where}.reader`, [1, 2] as const);
  if (reader !== indexParser(program, role)) {
    f.fail(`${where}.reader`, `is ${reader}, but ${program} reads its ${role} with reader ${indexParser(program, role)}`);
  }
  const table = f.object(input.records, `${where}.records`, ['id', 'length', 'sha256']);
  const ids = f.array(table.id, `${where}.records.id`, Number.MAX_SAFE_INTEGER);
  const lengths = f.array(table.length, `${where}.records.length`, ids.length);
  const hashes = f.array(table.sha256, `${where}.records.sha256`, ids.length);
  if (lengths.length !== ids.length || hashes.length !== ids.length) {
    f.fail(`${where}.records`, `has ${ids.length} IDs, ${lengths.length} lengths and ${hashes.length} SHA-256s`);
  }
  // One by one, the path built only on a failure: a table may have a million records.
  for (let k = 0; k < ids.length; k++) {
    if (k > 0 && k % RECORDS_PER_STEP === 0) yield;
    if (typeof ids[k] !== 'string' || (ids[k] as string).length > f.limits.idChars) f.text(ids[k], `${where}.records.id[${k}]`, f.limits.idChars);
    const length = lengths[k];
    if (typeof length !== 'number' || !Number.isSafeInteger(length) || length < 0) f.count(length, `${where}.records.length[${k}]`);
    if (typeof hashes[k] !== 'string' || !SHA256.test(hashes[k] as string)) f.sha256(hashes[k], `${where}.records.sha256[${k}]`);
  }
  const sources = f.array(input.sources, `${where}.sources`, f.limits.sources).map((source, i) => checkSource(f, source, `${where}.sources[${i}]`));
  if (sources.length === 0) f.fail(`${where}.sources`, 'is empty');
  const included = sources.reduce((sum, source) => sum + source.records - source.excluded.length, 0);
  if (included !== ids.length) f.fail(`${where}.sources`, `include ${included} records, but ${where}.records lists ${ids.length}`);
  // The long arrays are frozen here, one in each step (deepFreeze then passes them by).
  for (const column of [ids, lengths, hashes, ...sources.map((source) => source.excluded)]) {
    if (column.length >= RECORDS_PER_STEP) {
      yield;
      f.frozen.add(column);
    }
    Object.freeze(column);
  }
  return {
    name: f.text(input.name, `${where}.name`, f.limits.nameChars, false),
    sha256: f.sha256(input.sha256, `${where}.sha256`),
    length: f.count(input.length, `${where}.length`),
    reader,
    records: { id: ids as string[], length: lengths as number[], sha256: hashes as string[] },
    sources,
  };
}

function checkSource(f: Fields, value: unknown, where: string): SessionSource {
  const source = f.object(value, where, ['name', 'size', 'records', 'excluded']);
  const records = f.count(source.records, `${where}.records`);
  const excluded = f.array(source.excluded, `${where}.excluded`, records).map((index, i) => f.count(index, `${where}.excluded[${i}]`));
  excluded.forEach((index, i) => {
    if (index >= records) f.fail(`${where}.excluded[${i}]`, `is ${index}, not a record of the source's ${records}`);
    if (i > 0 && index <= excluded[i - 1]!) f.fail(`${where}.excluded`, 'is not in ascending order without repeats');
  });
  return {
    name: f.text(source.name, `${where}.name`, f.limits.nameChars, false),
    size: f.count(source.size, `${where}.size`),
    records,
    excluded,
  };
}

/**
 * Checks a parsed candidates block against the manifest and returns the candidates in tray order:
 * each names a run of the file, an HSP index within the run's HSP records, a query record of the
 * run, and appears once. Whether the HSP record at `index` has that query and rank is checked
 * when the run's records are read (application/session.ts).
 */
export function checkCandidates(value: unknown, manifest: SessionManifest, limits: SessionLimits = SESSION_LIMITS): readonly SessionCandidate[] {
  const f: Fields = new Fields('candidates block', limits);
  const block = f.object(value, 'the candidates block', ['candidates']);
  const seen = new Map<string, number>();
  return Object.freeze(
    f.array(block.candidates, 'candidates', limits.candidates).map((item, i): SessionCandidate => {
      const where = `candidates[${i}]`;
      const candidate = f.object(item, where, ['run', 'index', 'qIdx', 'rank', 'note', 'addedAt']);
      const position = f.positive(candidate.run, `${where}.run`);
      const run = manifest.runs[position - 1];
      if (run === undefined) f.fail(`${where}.run`, `is ${position}, but the file holds ${manifest.runs.length} runs`);
      const index = f.count(candidate.index, `${where}.index`);
      if (index >= run.hitCount) f.fail(`${where}.index`, `is ${index}, beyond the ${run.hitCount} HSP records of run ${position} in the file`);
      const qIdx = f.count(candidate.qIdx, `${where}.qIdx`);
      if (qIdx >= run.query.records.id.length) {
        f.fail(`${where}.qIdx`, `is ${qIdx}, beyond the ${run.query.records.id.length} query records of run ${position} in the file`);
      }
      const rank = f.count(candidate.rank, `${where}.rank`);
      const key = `${position}/${qIdx}/${rank}`;
      const twice = seen.get(key);
      if (twice !== undefined) f.fail(where, `is the same HSP as candidates[${twice}]`);
      seen.set(key, i);
      return Object.freeze({
        run: position,
        index,
        qIdx,
        rank,
        note: f.text(candidate.note, `${where}.note`, limits.noteChars),
        addedAt: f.count(candidate.addedAt, `${where}.addedAt`),
      });
    }),
  );
}

function deepFreeze<T>(value: T, frozen?: WeakSet<object>): T {
  if (typeof value === 'object' && value !== null && !frozen?.has(value) && !Object.isFrozen(value)) {
    Object.freeze(value);
    for (const item of Object.values(value)) deepFreeze(item, frozen);
  }
  return value;
}

// --- writing --------------------------------------------------------------------------------------

/** The container's first line. */
export const containerHeader = (): string => `${SESSION_CONTAINER} ${SESSION_CONTAINER_VERSION}\n`;
/** The line after the last block. */
export const containerEnd = (): string => `${SESSION_END}\n`;
/** The line before a block's bytes; a line end follows the bytes. */
export const blockHeader = (name: string, length: number): string => `${name} ${length}\n`;
export const runBlockName = (position: number, stream: SessionStream): string => `run${position}.${stream}`;

/** The name of a session file saved at `time` (local time); it never names an input. */
export function sessionFileName(time: number): string {
  const date = new Date(time);
  const two = (n: number) => String(n).padStart(2, '0');
  const day = `${date.getFullYear()}${two(date.getMonth() + 1)}${two(date.getDate())}`;
  return `losat-session-${day}-${two(date.getHours())}${two(date.getMinutes())}${two(date.getSeconds())}.losat-session.gz`;
}

/** Whether bytes begin as gzip data does (RFC 1952: 0x1f 0x8b). */
export const isGzip = (head: Uint8Array): boolean => head.length >= 2 && head[0] === 0x1f && head[1] === 0x8b;

// --- reading --------------------------------------------------------------------------------------

/** What the reader found, in the file's order. */
export type SessionEvent =
  | { readonly type: 'manifest'; readonly manifest: SessionManifest }
  /** The first block of run `run` (1-based position in the file) begins. */
  | { readonly type: 'run-start'; readonly run: number }
  /** Bytes of a run's block, in order: a view of the bytes given to `push`, valid until the next `push`. */
  | { readonly type: 'data'; readonly run: number; readonly stream: SessionStream; readonly bytes: Uint8Array }
  /** The last block of the run ended. */
  | { readonly type: 'run-end'; readonly run: number }
  | { readonly type: 'candidates'; readonly candidates: readonly SessionCandidate[] }
  | { readonly type: 'end' };

interface Expected {
  readonly name: string;
  readonly run?: number;
  readonly stream?: SessionStream;
  /** The block's length, which its header must state. */
  readonly length?: number;
  /** The largest length its header may state, for the manifest and the candidates. */
  readonly max?: number;
}

const LF = 0x0a;

/**
 * Reads a container (the bytes after gzip) incrementally: `push` each piece in order, then
 * `finish`. Every refusal is a SessionFileError that says what is wrong and where; after one, the
 * reader takes nothing more. Only the first line, the block header lines, the manifest and the
 * candidates are buffered, and each within its limit.
 */
export class SessionFileReader {
  private state: 'first' | 'header' | 'body' | 'body-end' | 'done' | 'failed' = 'first';
  private line: number[] = [];
  /** The blocks in the order they must come: the manifest, then those that it lists. */
  private readonly expected: Expected[];
  private next = 0;
  private block: { readonly expected: Expected; readonly length: number; remaining: number; readonly parts: Uint8Array[] } | undefined;
  private manifest: SessionManifest | undefined;
  private started = false;

  constructor(private readonly limits: SessionLimits = SESSION_LIMITS) {
    this.expected = [{ name: MANIFEST_BLOCK, max: limits.manifestBytes }];
  }

  push(bytes: Uint8Array): SessionEvent[] {
    const events: SessionEvent[] = [];
    try {
      runSteps(this.read(bytes, events));
    } catch (error) {
      this.state = 'failed';
      throw error;
    }
    return events;
  }

  /**
   * `push` with `pause` between the steps of reading the manifest (its parse, then the steps of
   * `checkManifestPaced`), so that the page draws while a large manifest is read. Call it again
   * only once it has settled.
   */
  async pushPaced(bytes: Uint8Array, pause: () => Promise<void>): Promise<SessionEvent[]> {
    const events: SessionEvent[] = [];
    try {
      await runStepsPaced(this.read(bytes, events), pause);
    } catch (error) {
      this.state = 'failed';
      throw error;
    }
    return events;
  }

  /** The end of the bytes: refuses a container that ended early. */
  finish(): void {
    const where = this.block === undefined ? '' : `block "${this.block.expected.name}"`;
    switch (this.state) {
      case 'done':
        return;
      case 'failed':
        throw new SessionFileError(`${DAMAGED}.`);
      case 'first':
        this.state = 'failed';
        throw new SessionFileError(this.started ? `${NOT_SESSION}: it ends inside its first line.` : `${NOT_SESSION}: it is empty.`);
      case 'header': {
        this.state = 'failed';
        const after = this.next === 0 ? 'its first line' : `block "${this.expected[this.next - 1]!.name}"`;
        const missing = this.next < this.expected.length ? `block "${this.expected[this.next]!.name}"` : 'its end line';
        throw new SessionFileError(
          this.line.length > 0 ? `${INCOMPLETE}: it ends inside the line after ${after}.` : `${INCOMPLETE}: it ends after ${after}, before ${missing}.`,
        );
      }
      case 'body':
        this.state = 'failed';
        throw new SessionFileError(
          `${INCOMPLETE}: it ends inside ${where} (${this.block!.length - this.block!.remaining} of its ${this.block!.length} bytes).`,
        );
      case 'body-end':
        this.state = 'failed';
        throw new SessionFileError(`${INCOMPLETE}: it ends after the bytes of ${where}, before the line end that closes it.`);
    }
  }

  private *read(bytes: Uint8Array, events: SessionEvent[]): Steps<void> {
    if (this.state === 'failed') throw new SessionFileError(`${DAMAGED}.`);
    let at = 0;
    if (bytes.length > 0) this.started = true;
    while (at < bytes.length) {
      switch (this.state) {
        case 'first':
        case 'header': {
          const lf = bytes.indexOf(LF, at);
          const end = lf < 0 ? bytes.length : lf;
          if (this.line.length + (end - at) > this.limits.lineBytes) this.longLine();
          for (let i = at; i < end; i++) this.line.push(bytes[i]!);
          at = end;
          if (lf < 0) break;
          at++;
          const text = String.fromCharCode(...this.line);
          this.line = [];
          if (this.state === 'first') this.firstLine(text);
          else this.headerLine(text, events);
          break;
        }
        case 'body': {
          const block = this.block!;
          const take = Math.min(block.remaining, bytes.length - at);
          const piece = bytes.subarray(at, at + take);
          if (block.expected.run !== undefined) events.push({ type: 'data', run: block.expected.run, stream: block.expected.stream!, bytes: piece });
          else block.parts.push(piece.slice());
          block.remaining -= take;
          at += take;
          if (block.remaining === 0) this.state = 'body-end';
          break;
        }
        case 'body-end': {
          const block = this.block!;
          if (bytes[at] !== LF) {
            throw new SessionFileError(`${DAMAGED}: block "${block.expected.name}" is longer than the ${block.length} bytes that its header states.`);
          }
          at++;
          yield* this.endBlock(events);
          break;
        }
        case 'done':
          throw new SessionFileError(`${DAMAGED}: there are bytes after its end line.`);
      }
    }
  }

  private longLine(): never {
    if (this.state === 'first') throw new SessionFileError(`${NOT_SESSION} (its first line is not "${SESSION_CONTAINER} ${SESSION_CONTAINER_VERSION}").`);
    const after = this.next === 0 ? 'its first line' : `block "${this.expected[this.next - 1]!.name}"`;
    throw new SessionFileError(`${DAMAGED}: the line after ${after} is not a block header (it is longer than ${this.limits.lineBytes} bytes).`);
  }

  private firstLine(text: string): void {
    const version = /^LOSAT-WEB-SESSION (\d{1,9})$/.exec(text);
    if (version !== null && Number(version[1]) > SESSION_CONTAINER_VERSION) {
      throw new SessionFileError(
        `The session file was saved by a newer LOSAT Web (its container is version ${Number(version[1])}; this LOSAT Web reads version ${SESSION_CONTAINER_VERSION}). Open it with a newer LOSAT Web.`,
      );
    }
    if (text !== `${SESSION_CONTAINER} ${SESSION_CONTAINER_VERSION}`) {
      throw new SessionFileError(`${NOT_SESSION} (its first line is ${printable(text)}, not "${SESSION_CONTAINER} ${SESSION_CONTAINER_VERSION}").`);
    }
    this.state = 'header';
  }

  private headerLine(text: string, events: SessionEvent[]): void {
    const wanted = this.expected[this.next];
    const after = this.next === 0 ? 'its first line' : `block "${this.expected[this.next - 1]!.name}"`;
    if (text === SESSION_END) {
      if (wanted !== undefined) throw new SessionFileError(`${INCOMPLETE}: its end line comes after ${after}, where block "${wanted.name}" belongs.`);
      this.state = 'done';
      events.push({ type: 'end' });
      return;
    }
    const header = /^([a-z0-9.]{1,32}) (0|[1-9]\d{0,15})$/.exec(text);
    if (header === null) throw new SessionFileError(`${DAMAGED}: the line after ${after} is not a block header (${printable(text)}).`);
    const [, name, digits] = header as unknown as [string, string, string];
    if (wanted === undefined) throw new SessionFileError(`${DAMAGED}: block "${name}" comes after ${after}, where its end line belongs.`);
    if (name !== wanted.name) throw new SessionFileError(`${DAMAGED}: block "${name}" comes after ${after}, where block "${wanted.name}" belongs.`);
    const length = Number(digits);
    if (!Number.isSafeInteger(length)) throw new SessionFileError(`${DAMAGED}: the length of block "${name}" is too large.`);
    if (wanted.length !== undefined && length !== wanted.length) {
      throw new SessionFileError(`${DAMAGED}: block "${name}" has ${length} bytes, but the manifest gives ${wanted.length}.`);
    }
    if (wanted.max !== undefined && length > wanted.max) {
      throw new SessionFileError(`The session file's block "${name}" has ${length} bytes, more than the ${wanted.max} that LOSAT Web reads.`);
    }
    this.next++;
    this.block = { expected: wanted, length, remaining: length, parts: [] };
    if (wanted.run !== undefined && wanted.stream === SESSION_STREAMS[0]) events.push({ type: 'run-start', run: wanted.run });
    this.state = length === 0 ? 'body-end' : 'body';
  }

  private *endBlock(events: SessionEvent[]): Steps<void> {
    const { expected, parts } = this.block!;
    this.block = undefined;
    this.state = 'header';
    if (expected.run !== undefined) {
      if (expected.stream === SESSION_STREAMS[SESSION_STREAMS.length - 1]) events.push({ type: 'run-end', run: expected.run });
      return;
    }
    const json = parseJson(parts, expected.name);
    if (expected.name === MANIFEST_BLOCK) {
      yield;
      const manifest = yield* manifestSteps(json, this.limits);
      this.manifest = manifest;
      manifest.runs.forEach((run, k) => {
        for (const stream of SESSION_STREAMS) this.expected.push({ name: runBlockName(k + 1, stream), run: k + 1, stream, length: run.blocks[stream] });
      });
      if (manifest.candidates) this.expected.push({ name: CANDIDATES_BLOCK, max: this.limits.candidatesBytes });
      events.push({ type: 'manifest', manifest });
    } else {
      events.push({ type: 'candidates', candidates: checkCandidates(json, this.manifest!, this.limits) });
    }
  }
}

function parseJson(parts: readonly Uint8Array[], name: string): unknown {
  let text: string;
  try {
    const decoder = new TextDecoder('utf-8', { fatal: true });
    text = parts.map((part, i) => decoder.decode(part, { stream: i < parts.length - 1 })).join('') + decoder.decode();
  } catch {
    throw new SessionFileError(`${DAMAGED}: block "${name}" is not UTF-8 text.`);
  }
  try {
    return JSON.parse(text) as unknown;
  } catch (error) {
    throw new SessionFileError(`${DAMAGED}: block "${name}" is not JSON (${error instanceof Error ? error.message : String(error)}).`);
  }
}

/** A line of the file in a message: quoted, its first 40 characters, controls escaped. */
function printable(text: string): string {
  const shown = text.length > 40 ? `${text.slice(0, 40)}…` : text;
  return JSON.stringify(shown);
}

// --- what a loaded run must agree with -------------------------------------------------------------

/** What the HSP records of a loaded run may name: its record tables and outputs, as the manifest gives them. */
export interface HspRecordBounds {
  /** HSP records of the run (`hitCount`). */
  readonly count: number;
  /** Records of the run's query and subject record tables. */
  readonly queries: number;
  readonly subjects: number;
  /** Bytes of the run's outfmt 0 and outfmt 6. */
  readonly out0: number;
  readonly out6: number;
}

export const hspRecordBounds = (run: SessionRun): HspRecordBounds => ({
  count: run.hitCount,
  queries: run.query.records.id.length,
  subjects: run.subject.records.id.length,
  out0: run.blocks.out0,
  out6: run.blocks.out6,
});

const isCount = (value: unknown): value is number => typeof value === 'number' && Number.isSafeInteger(value) && value >= 0;
const isFrame = (value: unknown) => value === null || (typeof value === 'number' && Number.isInteger(value) && value !== 0 && Math.abs(value) <= 3);
const isRange = (value: unknown, length: number) =>
  value === null || (Array.isArray(value) && value.length === 2 && isCount(value[0]) && isCount(value[1]) && value[0] <= value[1] && value[1] <= length);

/**
 * Checks the HSP records of a run loaded from a session file one by one, in the order of their
 * lines, as the JSON values that the file holds (docs/web/abi_v2.md §8), before anything coerces
 * them (domain/hsp-table.ts makes typed arrays of them, and the exports write them as they are).
 * A record is refused when a field is missing or of the wrong type: `index`, `q_idx`, `s_idx` and
 * `rank` whole numbers of 0 or more, coordinates whole numbers of 1 or more, frames null or -3 to
 * 3 other than 0, scores numbers, `subject_length` null or a whole number, aligned rows null or
 * text, byte ranges null or [start, end] within their output; or when an index is outside 0 to
 * count - 1 or repeated, a record index is beyond the record tables, or a query has a rank twice.
 * Fields that the ABI may add later are left alone.
 */
export class HspRecordCheck {
  private records = 0;
  private readonly seen: Uint8Array;
  /** The ranks of each query record seen so far. */
  private readonly ranks = new Map<number, Set<number>>();

  constructor(private readonly bounds: HspRecordBounds) {
    this.seen = new Uint8Array(bounds.count);
  }

  /** Why the next record (the parsed JSON of its line) is not one that the run can have, or undefined. */
  next(value: unknown): string | undefined {
    const { count, queries, subjects, out0, out6 } = this.bounds;
    const where = `HSP record ${++this.records}`;
    if (this.records > count) return `there are more HSP records than the ${count} that the manifest gives`;
    if (typeof value !== 'object' || value === null || Array.isArray(value)) return `${where} is not a JSON object`;
    const record = value as Record<string, unknown>;
    const wrong = (field: string, expected: string) =>
      record[field] === undefined ? `${where} has no ${field}` : `${where} has ${field} ${shown(record[field])}, not ${expected}`;
    for (const field of ['index', 'q_idx', 's_idx', 'rank']) if (!isCount(record[field])) return wrong(field, 'a whole number of 0 or more');
    const [index, q, s, rank] = [record.index, record.q_idx, record.s_idx, record.rank] as number[];
    if (index! >= count || this.seen[index!] === 1) return `${where} has an index that is outside 0 to ${count - 1} or repeated`;
    this.seen[index!] = 1;
    if (q! >= queries) return `${where} names query record ${q}, but the run has ${queries} query records`;
    if (s! >= subjects) return `${where} names subject record ${s}, but the run has ${subjects} subject records`;
    let ranks = this.ranks.get(q!);
    if (ranks === undefined) this.ranks.set(q!, (ranks = new Set()));
    if (ranks.has(rank!)) return `${where} has a rank that another HSP of its query has`;
    ranks.add(rank!);
    for (const field of ['q_start', 'q_end', 's_start', 's_end']) {
      if (!isCount(record[field]) || record[field] === 0) return wrong(field, 'a coordinate (a whole number of 1 or more)');
    }
    for (const field of ['query_frame', 'subject_frame']) if (!isFrame(record[field])) return wrong(field, 'null or a frame of -3 to 3 other than 0');
    // The adapter writes null for a value that is not finite (web/adapter/src/json.rs `number`), so
    // a run that the engine wrote can hold it; refusing it would make a saved session unopenable.
    for (const field of ['raw_score', 'bit_score', 'e_value']) {
      if (record[field] !== null && (typeof record[field] !== 'number' || !Number.isFinite(record[field]))) return wrong(field, 'null or a number');
    }
    if (record.subject_length !== null && !isCount(record.subject_length)) return wrong('subject_length', 'null or a whole number of 0 or more');
    for (const field of ['query_aligned', 'subject_aligned']) {
      if (record[field] !== null && typeof record[field] !== 'string') return wrong(field, 'null or text');
    }
    if (!isRange(record.out6, out6)) return wrong('out6', `null or a byte range within the run's ${out6} bytes of outfmt 6`);
    for (const field of ['out0', 'out0_subject']) {
      if (!isRange(record[field], out0)) return wrong(field, `null or a byte range within the run's ${out0} bytes of outfmt 0`);
    }
    return undefined;
  }

  /** After the last record: why their count is not the manifest's, or undefined. */
  finish(): string | undefined {
    const { count } = this.bounds;
    return this.records === count ? undefined : `there are ${this.records} HSP records, but the manifest gives ${count}`;
  }
}

/** A JSON value in a message: its first 40 characters. */
function shown(value: unknown): string {
  // JSON.parse makes a number too large for a double Infinity, which JSON.stringify would write as null.
  const text = typeof value === 'number' && !Number.isFinite(value) ? String(value) : (JSON.stringify(value) ?? String(value));
  return text.length > 40 ? `${text.slice(0, 40)}…` : text;
}

/** A record of a run input as the re-attachment compares it: its ID, length and SHA-256. */
export interface RecordIdentity {
  readonly id: string;
  readonly length: number;
  readonly sha256: string;
}

/** A file chosen as (a part of) a loaded run's original FASTA, indexed with the recorded reader kind: every record of it, in order. */
export interface ChosenFile {
  readonly name: string;
  readonly records: readonly RecordIdentity[];
}

/** The chosen file of each recorded source (`files[i]` is the index of source i's file), or why the files are not the sources. */
export type SourceMatch = { readonly ok: true; readonly files: readonly number[] } | { readonly ok: false; readonly message: string };

/** The records of a file that a source's exclusions (0-based indices, ascending) leave in. */
function includedOf(records: readonly RecordIdentity[], excluded: readonly number[]): RecordIdentity[] {
  const included: RecordIdentity[] = [];
  let next = 0;
  for (let k = 0; k < records.length; k++) {
    if (excluded[next] === k) next++;
    else included.push(records[k]!);
  }
  return included;
}

/**
 * Matches files chosen as a loaded run's original FASTA to the sources that the session file
 * recorded, whatever order they were chosen in (REQ-23; file dialogs rarely let the user order a
 * selection): a file is a source's when it has as many records as the source had and, with the
 * source's exclusions applied, the records that the run input has from that source (IDs, lengths
 * and SHA-256s). As many files as sources are expected. When some source has no file, the message
 * names the first problem of a pairing that keeps the matches found and gives each other source a
 * file left over (one with its record count first): a file's record count (`sourceMismatch`), or
 * the first record of the input that differs (`recordsMismatch`).
 */
export function matchSources(saved: SessionInput, chosen: readonly ChosenFile[]): SourceMatch {
  const { sources } = saved;
  const table = saved.records;
  const starts: number[] = [];
  let at = 0;
  for (const source of sources) {
    starts.push(at);
    at += source.records - source.excluded.length;
  }
  const fits = (i: number, j: number): boolean => {
    const source = sources[i]!;
    if (chosen[j]!.records.length !== source.records) return false;
    return includedOf(chosen[j]!.records, source.excluded).every((record, k) => {
      const position = starts[i]! + k;
      return record.id === table.id[position] && record.length === table.length[position] && record.sha256 === table.sha256[position];
    });
  };
  const fitting = sources.map((_, i) => chosen.map((_, j) => j).filter((j) => fits(i, j)));
  // Augmenting paths (Kuhn): a file that fits two sources (the same file joined twice) goes where it is needed.
  const fileOf = sources.map(() => -1);
  const sourceOf = chosen.map(() => -1);
  const assign = (i: number, visited: boolean[]): boolean => {
    for (const j of fitting[i]!) {
      if (visited[j]) continue;
      visited[j] = true;
      if (sourceOf[j]! < 0 || assign(sourceOf[j]!, visited)) {
        sourceOf[j] = i;
        fileOf[i] = j;
        return true;
      }
    }
    return false;
  };
  sources.forEach((_, i) => assign(i, chosen.map(() => false)));
  if (fileOf.every((j) => j >= 0)) return { ok: true, files: fileOf };
  const left = chosen.map((_, j) => j).filter((j) => sourceOf[j]! < 0);
  const pairing = fileOf.map((j, i) => {
    if (j >= 0 || left.length === 0) return j;
    const same = left.findIndex((k) => chosen[k]!.records.length === sources[i]!.records);
    return left.splice(Math.max(0, same), 1)[0]!;
  });
  for (const [i, j] of pairing.entries()) {
    if (j < 0) continue;
    const mismatch = sourceMismatch(sources[i]!, i, chosen[j]!.name, chosen[j]!.records.length);
    if (mismatch !== undefined) return { ok: false, message: mismatch };
  }
  const records = pairing.flatMap((j, i) => (j < 0 ? [] : includedOf(chosen[j]!.records, sources[i]!.excluded)));
  return { ok: false, message: recordsMismatch(saved, records) ?? `the chosen files do not have the records of the run's input` };
}

/**
 * Why a chosen file cannot be source `position` (0-based) of a loaded run's input, or undefined:
 * its record table does not have as many records as the saved source had.
 */
export function sourceMismatch(saved: SessionSource, position: number, fileName: string, records: number): string | undefined {
  if (records === saved.records) return undefined;
  return (
    `${JSON.stringify(fileName)} has ${records} ${records === 1 ? 'record' : 'records'}, but file ${position + 1} of the saved input ` +
    `(${JSON.stringify(saved.name)}) had ${saved.records}`
  );
}

/**
 * Why the included records of chosen files are not the records of a loaded run's input, or
 * undefined: the first record whose ID, length or SHA-256 differs (the record SHA-256s tell which
 * record changed), or a different count.
 */
export function recordsMismatch(saved: SessionInput, chosen: readonly RecordIdentity[]): string | undefined {
  const { id, length, sha256 } = saved.records;
  const count = Math.min(id.length, chosen.length);
  for (let k = 0; k < count; k++) {
    const got = chosen[k]!;
    if (got.id !== id[k] || got.length !== length[k]) {
      return `record ${k + 1} is ${JSON.stringify(id[k])} (length ${length[k]}) in the saved run, but the chosen files give ${JSON.stringify(got.id)} (length ${got.length})`;
    }
    if (got.sha256 !== sha256[k]) return `record ${k + 1} (${JSON.stringify(id[k])}) differs from the saved run's record: its bytes have another SHA-256`;
  }
  if (chosen.length !== id.length) return `the chosen files give ${chosen.length} records after the exclusions, but the saved run searched ${id.length}`;
  return undefined;
}
