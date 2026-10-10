// The session file's container and checks (domain/session-file.ts; docs/web/session_file.md): a
// container written here is read back the same in pieces of any size, and every kind of damage is
// refused with a message that says what and where: the first line, the schema, the manifest's
// fields and limits, the order and lengths of the blocks, truncation at any byte, data after the
// end, and candidates or HSP records that the runs cannot have.
import { describe, expect, it } from 'vitest';
import {
  blockHeader,
  checkCandidates,
  checkManifest,
  checkManifestPaced,
  containerEnd,
  containerHeader,
  hspRecordBounds,
  HspRecordCheck,
  isGzip,
  matchSources,
  recordsMismatch,
  runBlockName,
  SESSION_LIMITS,
  SESSION_STREAMS,
  SessionFileError,
  SessionFileReader,
  sessionFileName,
  sourceMismatch,
  type SessionCandidate,
  type SessionEvent,
  type SessionInput,
  type SessionLimits,
  type SessionManifest,
  type SessionRun,
  type SessionStream,
} from '../../src/domain/session-file';

const encoder = new TextEncoder();
const decoder = new TextDecoder();
const hex = (c: string) => c.repeat(64);

function run(overrides: Partial<SessionRun> = {}): SessionRun {
  return {
    number: 3,
    title: 'Two subjects',
    program: 'blastn',
    argv: ['blastn', '-query', 'query.fa', '-subject', 'combined_subject.fa', '-evalue', '1e-5'],
    requestedThreads: 'auto',
    group: { index: 1, position: 1, size: 2 },
    queuedAt: 1000,
    record: { runtimePath: 'fake', threads: 1, engineBuild: 'fake-engine', startedAt: 1001, phaseTimes: { preparing: 1002 }, endedAt: 1100 },
    hitCount: 2,
    blocks: { out0: 12, out6: 30, out7: 0, hits: 40, diagnostics: 5 },
    query: {
      name: 'query.fa',
      sha256: hex('a'),
      length: 20,
      reader: 1,
      records: { id: ['q1', ''], length: [8, 4], sha256: [hex('1'), hex('2')] },
      sources: [{ name: 'query.fa', size: 20, records: 2, excluded: [] }],
    },
    subject: {
      name: 'combined_subject.fa',
      sha256: hex('b'),
      length: 50,
      reader: 1,
      records: { id: ['s1', 's3'], length: [10, 12], sha256: [hex('3'), hex('5')] },
      sources: [
        { name: 'a.fa', size: 30, records: 2, excluded: [1] },
        { name: 'b.fa', size: 20, records: 1, excluded: [] },
      ],
    },
    ...overrides,
  };
}

function manifest(overrides: Partial<SessionManifest> = {}): SessionManifest {
  return {
    format: 'losat-web-session',
    schema: 1,
    app: { version: '0.0.0', build: 'abc1234' },
    savedAt: 1_700_000_000_000,
    candidates: true,
    runs: [run(), run({ number: 4, group: { index: 1, position: 2, size: 2 }, program: 'tblastn', argv: tblastnArgv(), query: proteinQuery() })],
    ...overrides,
  };
}

const tblastnArgv = () => ['tblastn', '-query', 'p.fa', '-subject', 'combined_subject.fa'];
const proteinQuery = () => ({
  name: 'p.fa',
  sha256: hex('c'),
  length: 9,
  reader: 2 as const,
  records: { id: ['p1'], length: [5], sha256: [hex('6')] },
  sources: [{ name: 'p.fa', size: 9, records: 1, excluded: [] }],
});

const CANDIDATES: readonly SessionCandidate[] = [
  { run: 1, index: 1, qIdx: 0, rank: 1, note: 'look <b>here</b>', addedAt: 5 },
  { run: 2, index: 0, qIdx: 0, rank: 0, note: '', addedAt: 6 },
];

/** Bytes of each block of a run: distinct letters, so that a misplaced byte shows. */
function blockBytes(position: number, stream: SessionStream, length: number): Uint8Array {
  const letter = 'abcdefghij'[(position * 5 + SESSION_STREAMS.indexOf(stream)) % 10]!;
  return encoder.encode(letter.repeat(length));
}

interface Parts {
  readonly manifest?: unknown;
  readonly candidates?: unknown;
  /** Text appended or put in place of parts, for damaged files. */
  readonly firstLine?: string;
  readonly end?: string;
  /** Replaces a block's header line (its name and length) by this text, keeping its bytes. */
  readonly headers?: Readonly<Record<string, string>>;
  /** Leaves out a block. */
  readonly without?: string;
}

/** A container as the session writes one (without gzip). */
function container(m: SessionManifest = manifest(), candidates: readonly SessionCandidate[] = CANDIDATES, parts: Parts = {}): Uint8Array {
  const chunks: Uint8Array[] = [];
  const text = (t: string) => chunks.push(encoder.encode(t));
  const block = (name: string, bytes: Uint8Array) => {
    if (parts.without === name) return;
    text(parts.headers?.[name] ?? blockHeader(name, bytes.length));
    chunks.push(bytes);
    text('\n');
  };
  text(parts.firstLine ?? containerHeader());
  block('manifest', encoder.encode(JSON.stringify(parts.manifest ?? m)));
  m.runs.forEach((r, k) => {
    for (const stream of SESSION_STREAMS) block(runBlockName(k + 1, stream), blockBytes(k + 1, stream, r.blocks[stream]));
  });
  if (m.candidates) block('candidates', encoder.encode(JSON.stringify(parts.candidates ?? { candidates })));
  text(parts.end ?? containerEnd());
  const out = new Uint8Array(chunks.reduce((sum, c) => sum + c.length, 0));
  let at = 0;
  for (const c of chunks) {
    out.set(c, at);
    at += c.length;
  }
  return out;
}

/** Reads a container in pieces of `size` bytes; the data events are joined per block. */
function read(bytes: Uint8Array, size = bytes.length, limits?: SessionLimits) {
  const reader = new SessionFileReader(limits);
  const events: SessionEvent[] = [];
  const blocks = new Map<string, number[]>();
  for (let at = 0; at < bytes.length; at += Math.max(1, size)) {
    for (const event of reader.push(bytes.slice(at, at + Math.max(1, size)))) {
      if (event.type === 'data') {
        const key = runBlockName(event.run, event.stream);
        blocks.set(key, [...(blocks.get(key) ?? []), ...event.bytes]);
        // Views are valid until the next push only; keep a copy above.
      } else {
        events.push(event);
      }
    }
  }
  reader.finish();
  return { events, blocks };
}

function refusal(action: () => unknown): string {
  try {
    action();
  } catch (error) {
    expect(error).toBeInstanceOf(SessionFileError);
    return (error as Error).message;
  }
  throw new Error('not refused');
}

describe('session file container', () => {
  it('reads back the manifest, the run blocks and the candidates, in pieces of any size', () => {
    const m = manifest();
    const bytes = container(m);
    for (const size of [1, 2, 7, 64, bytes.length]) {
      const { events, blocks } = read(bytes, size);
      expect(events.map((e) => e.type)).toEqual(['manifest', 'run-start', 'run-end', 'run-start', 'run-end', 'candidates', 'end']);
      expect((events[0] as { manifest: SessionManifest }).manifest).toEqual(m);
      expect((events[5] as { candidates: readonly SessionCandidate[] }).candidates).toEqual(CANDIDATES);
      m.runs.forEach((r, k) => {
        for (const stream of SESSION_STREAMS) {
          const got = blocks.get(runBlockName(k + 1, stream)) ?? [];
          expect(decoder.decode(new Uint8Array(got))).toBe(decoder.decode(blockBytes(k + 1, stream, r.blocks[stream])));
        }
      });
    }
  });

  it('reads a session without candidates, and empty blocks', () => {
    const m = manifest({ candidates: false });
    const { events } = read(container(m), 3);
    expect(events.map((e) => e.type)).toEqual(['manifest', 'run-start', 'run-end', 'run-start', 'run-end', 'end']);
  });

  it('freezes the manifest it returns and keeps only its known fields', () => {
    const checked = checkManifest(JSON.parse(JSON.stringify(manifest())));
    expect(Object.isFrozen(checked)).toBe(true);
    expect(Object.isFrozen(checked.runs[0]!.subject.records.id)).toBe(true);
  });

  it('names its files by the time of saving, never by an input', () => {
    const name = sessionFileName(new Date(2026, 9, 10, 8, 5, 9).getTime());
    expect(name).toBe('losat-session-20261010-080509.losat-session.gz');
    expect(isGzip(new Uint8Array([0x1f, 0x8b, 8]))).toBe(true);
    expect(isGzip(encoder.encode('LOSAT'))).toBe(false);
  });
});

describe('session file refusals', () => {
  it('refuses a file that is not a session file, or of a newer container', () => {
    expect(refusal(() => read(encoder.encode('>q1\nACGT\n')))).toMatch(/^This is not a LOSAT Web session file \(its first line is ">q1"/);
    expect(refusal(() => read(encoder.encode('x'.repeat(200))))).toMatch(/^This is not a LOSAT Web session file/);
    expect(refusal(() => read(new Uint8Array()))).toBe('This is not a LOSAT Web session file: it is empty.');
    expect(refusal(() => read(container(manifest(), CANDIDATES, { firstLine: 'LOSAT-WEB-SESSION 2\n' })))).toMatch(
      /saved by a newer LOSAT Web \(its container is version 2/,
    );
  });

  it('refuses another format or schema, a newer schema as saved by a newer LOSAT Web', () => {
    expect(refusal(() => read(container(manifest(), CANDIDATES, { manifest: { ...manifest(), schema: 2 } })))).toMatch(
      /^The session file was saved by a newer LOSAT Web \(its schema is 2; this LOSAT Web reads schema 1\)/,
    );
    expect(refusal(() => checkManifest({ ...manifest(), schema: 0 }))).toMatch(/manifest is not valid: schema is 0, not 1/);
    expect(refusal(() => checkManifest({ ...manifest(), schema: '1' }))).toMatch(/schema is not a whole number/);
    expect(refusal(() => checkManifest({ ...manifest(), format: 'other' }))).toMatch(/format is not "losat-web-session"/);
    expect(refusal(() => checkManifest([]))).toMatch(/the manifest is not a JSON object/);
  });

  it('refuses a manifest that is not JSON or not UTF-8', () => {
    const notJson = container(manifest(), CANDIDATES, { manifest: 'x' });
    // The manifest block of `container` is JSON.stringify of the value: a string is JSON; break it.
    const broken = encoder.encode(decoder.decode(notJson).replace('manifest 3\n"x"', 'manifest 3\n{"x'));
    expect(refusal(() => read(broken))).toMatch(/block "manifest" is not JSON/);
    const bytes = container();
    const at = decoder.decode(bytes).indexOf('{"format"');
    const invalid = bytes.slice();
    invalid[at + 1] = 0xff;
    expect(refusal(() => read(invalid))).toMatch(/block "manifest" is not UTF-8 text/);
  });

  it('names a missing field, an unknown field and a value of the wrong type', () => {
    const m = JSON.parse(JSON.stringify(manifest())) as Record<string, unknown> & { runs: Record<string, unknown>[] };
    for (const field of ['app', 'savedAt', 'candidates', 'runs']) {
      const copy = { ...m };
      delete copy[field];
      expect(refusal(() => checkManifest(copy))).toContain(`${field} is missing`);
    }
    for (const field of ['number', 'program', 'argv', 'requestedThreads', 'queuedAt', 'record', 'hitCount', 'blocks', 'query', 'subject']) {
      const r = { ...m.runs[0]! };
      delete r[field];
      expect(refusal(() => checkManifest({ ...m, runs: [r] }))).toContain(`runs[0].${field} is missing`);
    }
    expect(refusal(() => checkManifest({ ...m, opfsPath: 'tmp/x' }))).toContain('the manifest.opfsPath is not a field of a session file');
    expect(refusal(() => checkManifest({ ...m, runs: [{ ...m.runs[0]!, runId: 'r' }] }))).toContain('runs[0].runId is not a field');
    const cases: Array<[string, (r: Record<string, unknown>) => void, RegExp]> = [
      ['number', (r) => (r.number = '3'), /runs\[0\]\.number is not a whole number/],
      ['number', (r) => (r.number = 1.5), /runs\[0\]\.number is not a whole number/],
      ['number', (r) => (r.number = 0), /runs\[0\]\.number is 0/],
      ['queuedAt', (r) => (r.queuedAt = -1), /runs\[0\]\.queuedAt is not a whole number of 0 or more/],
      ['hitCount', (r) => (r.hitCount = 2 ** 53), /runs\[0\]\.hitCount is not a whole number/],
      ['program', (r) => (r.program = 'megablast'), /runs\[0\]\.program is not one of/],
      ['argv', (r) => (r.argv = ['blastn', '-query', 'query.fa', '-subject', 7]), /runs\[0\]\.argv\[4\] is not text/],
      ['title', (r) => (r.title = '  '), /runs\[0\]\.title is empty/],
      ['requestedThreads', (r) => (r.requestedThreads = 'many'), /runs\[0\]\.requestedThreads is not a whole number/],
      ['record', (r) => (r.record = { runtimePath: 'gpu' }), /runs\[0\]\.record\.runtimePath is not one of/],
      ['record', (r) => (r.record = { error: 'x' }), /runs\[0\]\.record\.error is not a field/],
      ['record', (r) => (r.record = { phaseTimes: { running: 'soon' } }), /runs\[0\]\.record\.phaseTimes\.running is not a whole number/],
      ['blocks', (r) => (r.blocks = { out0: 1, out6: 1, out7: 1, hits: 1 }), /runs\[0\]\.blocks\.diagnostics is missing/],
      ['group', (r) => (r.group = { index: 1, position: 3, size: 2 }), /runs\[0\]\.group\.position is 3, past the group's size 2/],
      ['query', (r) => ((r.query as Record<string, unknown>).sha256 = 'A'.repeat(64)), /runs\[0\]\.query\.sha256 is not a SHA-256/],
      ['query', (r) => ((r.query as Record<string, unknown>).reader = 2), /runs\[0\]\.query\.reader is 2, but blastn reads its query with reader 1/],
      [
        'records',
        (r) => (((r.subject as Record<string, unknown>).records as Record<string, unknown>).length = [10]),
        /runs\[0\]\.subject\.records has 2 IDs, 1 lengths and 2 SHA-256s/,
      ],
      [
        'records',
        (r) => (((r.subject as Record<string, unknown>).records as Record<string, unknown>).id = ['s1', 5]),
        /runs\[0\]\.subject\.records\.id\[1\] is not text/,
      ],
      [
        'excluded',
        (r) => ((r.subject as { sources: Record<string, unknown>[] }).sources[0]!.excluded = [2]),
        /runs\[0\]\.subject\.sources\[0\]\.excluded\[0\] is 2, not a record of the source's 2/,
      ],
      [
        'excluded',
        (r) => {
          const source = (r.subject as { sources: Record<string, unknown>[] }).sources[0]!;
          source.records = 3;
          source.excluded = [1, 1];
        },
        /runs\[0\]\.subject\.sources\[0\]\.excluded is not in ascending order without repeats/,
      ],
      ['argv', (r) => (r.argv = ['blastp', '-query', 'query.fa', '-subject', 'combined_subject.fa']), /runs\[0\]\.argv does not begin with "blastn -query/],
      ['argv', (r) => (r.argv = ['blastn', '-query', 'other.fa', '-subject', 'combined_subject.fa']), /runs\[0\]\.argv\[2\] is not the query's name "query\.fa"/],
    ];
    for (const [, change, message] of cases) {
      const r = JSON.parse(JSON.stringify(m.runs[0])) as Record<string, unknown>;
      change(r);
      expect(refusal(() => checkManifest({ ...m, runs: [r] }))).toMatch(message);
    }
  });

  it('checks that the sources and their exclusions make the record table', () => {
    const m = JSON.parse(JSON.stringify(manifest())) as SessionManifest & { runs: Array<{ subject: { sources: unknown[] } }> };
    m.runs[0]!.subject.sources = [{ name: 'a.fa', size: 30, records: 2, excluded: [1] }];
    expect(refusal(() => checkManifest(m))).toMatch(/runs\[0\]\.subject\.sources include 1 records, but runs\[0\]\.subject\.records lists 2/);
  });

  it('refuses values over the limits', () => {
    const small: SessionLimits = { ...SESSION_LIMITS, runs: 1, argvWords: 6, argvWordChars: 8, titleChars: 5, nameChars: 10, idChars: 2, noteChars: 3 };
    const m = manifest();
    expect(refusal(() => checkManifest(m, small))).toMatch(/runs has more than 1 items/);
    const one = { ...m, runs: [m.runs[0]!] };
    expect(refusal(() => checkManifest(one, small))).toMatch(/runs\[0\]\.argv has more than 6 items/);
    const r = { ...m.runs[0]!, argv: ['blastn', '-query', 'query.fa', '-subject', 'combined_subject.fa'] };
    expect(refusal(() => checkManifest({ ...m, runs: [r] }, { ...small, argvWords: 9 }))).toMatch(/runs\[0\]\.argv\[4\] is longer than 8 characters/);
    expect(refusal(() => checkManifest({ ...m, runs: [r] }, { ...small, argvWords: 9, argvWordChars: 100 }))).toMatch(
      /runs\[0\]\.subject\.name is longer than 10 characters/,
    );
    expect(refusal(() => checkManifest({ ...m, runs: [r] }, { ...small, argvWords: 9, argvWordChars: 100, nameChars: 100, idChars: 1 }))).toMatch(
      /runs\[0\]\.query\.records\.id\[0\] is longer than 1 characters/,
    );
    expect(refusal(() => checkManifest({ ...m, runs: [{ ...r, title: 'x'.repeat(6) }] }, { ...SESSION_LIMITS, titleChars: 5 }))).toMatch(/runs\[0\]\.title is longer than 5 characters/);
    const notes = { candidates: [{ ...CANDIDATES[0]!, note: 'four' }] };
    expect(refusal(() => checkCandidates(notes, m, small))).toMatch(/candidates\[0\]\.note is longer than 3 characters/);
    // A manifest block over the limit is refused from its header, before it is read.
    expect(refusal(() => read(container(), undefined, { ...SESSION_LIMITS, manifestBytes: 100 }))).toMatch(
      /block "manifest" has \d+ bytes, more than the 100 that LOSAT Web reads/,
    );
    expect(refusal(() => read(container(), undefined, { ...SESSION_LIMITS, candidatesBytes: 10 }))).toMatch(/block "candidates" has \d+ bytes, more than the 10/);
  });

  it('refuses blocks that are missing, extra, out of order, or of another length', () => {
    expect(refusal(() => read(container(manifest(), CANDIDATES, { without: 'run1.out7' })))).toMatch(
      /block "run1\.hits" comes after block "run1\.out6", where block "run1\.out7" belongs/,
    );
    expect(refusal(() => read(container(manifest(), CANDIDATES, { without: 'candidates' })))).toMatch(
      /its end line comes after block "run2\.diagnostics", where block "candidates" belongs/,
    );
    expect(refusal(() => read(container(manifest({ candidates: false }), CANDIDATES, { end: 'candidates 2\n{}\nLOSAT-WEB-SESSION-END\n' })))).toMatch(
      /block "candidates" comes after block "run2\.diagnostics", where its end line belongs/,
    );
    expect(refusal(() => read(container(manifest(), CANDIDATES, { headers: { 'run1.out6': 'run1.out0 30\n' } })))).toMatch(
      /block "run1\.out0" comes after block "run1\.out0", where block "run1\.out6" belongs/,
    );
    expect(refusal(() => read(container(manifest(), CANDIDATES, { headers: { 'run2.hits': 'run2.hits 39\n' } })))).toMatch(
      /block "run2\.hits" has 39 bytes, but the manifest gives 40/,
    );
    expect(refusal(() => read(container(manifest(), CANDIDATES, { headers: { 'run1.out0': 'run1.out0 012\n' } })))).toMatch(
      /the line after block "manifest" is not a block header \("run1\.out0 012"\)/,
    );
    expect(refusal(() => read(container(manifest(), CANDIDATES, { headers: { 'run1.out0': `run1.out0 ${'9'.repeat(16)}\n` } })))).toMatch(
      /the length of block "run1\.out0" is too large|has 9999999999999999 bytes/,
    );
  });

  it('refuses a block longer than its header states, and bytes after the end line', () => {
    const bytes = container();
    const text = decoder.decode(bytes);
    // One byte more in run1.out0 than its header (and the manifest) states.
    const longer = text.replace(`run1.out0 12\n${'f'.repeat(12)}\n`, `run1.out0 12\n${'f'.repeat(13)}\n`);
    expect(longer).not.toBe(text);
    expect(refusal(() => read(encoder.encode(longer)))).toMatch(/block "run1\.out0" is longer than the 12 bytes that its header states/);
    expect(refusal(() => read(encoder.encode(`${text}x`)))).toMatch(/there are bytes after its end line/);
  });

  it('refuses a container cut at any byte: at every block boundary and inside every block', () => {
    const bytes = container();
    const boundaries: number[] = [];
    for (let cut = 0; cut < bytes.length; cut++) {
      const message = refusal(() => read(bytes.slice(0, cut), 5));
      if (bytes[cut - 1] === 0x0a) boundaries.push(cut);
      expect(message).toMatch(/^(The session file is incomplete|This is not a LOSAT Web session file)/);
    }
    // After the first line, each block's header line, its bytes and its line end.
    expect(boundaries.length).toBeGreaterThan(2 + 11 * 2);
    expect(refusal(() => read(bytes.slice(0, decoder.decode(bytes).indexOf('run1.out6'))))).toMatch(
      /it ends after block "run1\.out0", before block "run1\.out6"/,
    );
    const inside = decoder.decode(bytes).indexOf('run2.hits 40\n') + 'run2.hits 40\n'.length + 7;
    expect(refusal(() => read(bytes.slice(0, inside)))).toMatch(/it ends inside block "run2\.hits" \(7 of its 40 bytes\)/);
    expect(refusal(() => read(bytes.slice(0, bytes.length - 1)))).toMatch(/it ends inside the line after block "candidates"/);
  });

  it('takes nothing more after a refusal', () => {
    const reader = new SessionFileReader();
    expect(() => reader.push(encoder.encode('nope\n'))).toThrow(SessionFileError);
    expect(() => reader.push(container())).toThrow(/damaged/);
    expect(() => reader.finish()).toThrow(/damaged/);
  });

  it('refuses candidates that name a run, an HSP or a query record that the file does not have, or one HSP twice', () => {
    const m = manifest();
    const check = (candidate: Partial<SessionCandidate>) => () => checkCandidates({ candidates: [{ ...CANDIDATES[0]!, ...candidate }] }, m);
    expect(refusal(check({ run: 3 }))).toMatch(/candidates\[0\]\.run is 3, but the file holds 2 runs/);
    expect(refusal(check({ index: 2 }))).toMatch(/candidates\[0\]\.index is 2, beyond the 2 HSP records of run 1 in the file/);
    expect(refusal(check({ qIdx: 2 }))).toMatch(/candidates\[0\]\.qIdx is 2, beyond the 2 query records of run 1 in the file/);
    expect(refusal(check({ rank: -1 }))).toMatch(/candidates\[0\]\.rank is not a whole number/);
    expect(refusal(check({ rank: 4294967296 }))).toMatch(/candidates\[0\]\.rank is 4294967296, beyond the 2 HSP records of run 1 in the file/);
    expect(refusal(check({ note: 5 as unknown as string }))).toMatch(/candidates\[0\]\.note is not text/);
    expect(refusal(() => checkCandidates({ candidates: [CANDIDATES[0], CANDIDATES[1], CANDIDATES[0]] }, m))).toMatch(
      /candidates\[2\] is the same HSP as candidates\[0\]/,
    );
    expect(refusal(() => checkCandidates({ list: [] }, m))).toMatch(/candidates is missing|list is not a field/);
    expect(refusal(() => read(container(m, CANDIDATES, { candidates: { candidates: [{ ...CANDIDATES[0]!, run: 9 }] } })))).toMatch(
      /The session file's candidates block is not valid: candidates\[0\]\.run is 9/,
    );
  });
});

describe('what a loaded run must agree with', () => {
  const r = run();
  /** An HSP record as the engine writes it (docs/web/abi_v2.md §8). */
  const hsp = (overrides: Record<string, unknown> = {}): Record<string, unknown> => ({
    index: 0,
    q_idx: 0,
    s_idx: 1,
    rank: 0,
    raw_score: 10,
    bit_score: 20.5,
    e_value: 1e-5,
    q_start: 1,
    q_end: 8,
    s_start: 12,
    s_end: 5,
    query_frame: null,
    subject_frame: null,
    subject_length: 20,
    query_aligned: 'ACGTACGT',
    subject_aligned: 'ACGTACGT',
    out6: [0, 15],
    out0: [0, 6],
    out0_subject: [0, 3],
    ...overrides,
  });
  const second = (overrides: Record<string, unknown> = {}) =>
    hsp({ index: 1, q_idx: 1, rank: 0, query_aligned: null, subject_aligned: null, out6: [15, 30], out0: null, out0_subject: null, ...overrides });
  /** The first problem of records checked in order, as their JSON lines would be (the JSON text, parsed). */
  const problem = (records: readonly unknown[]): string | undefined => {
    const check = new HspRecordCheck(hspRecordBounds(r));
    for (const record of records) {
      const found = check.next(JSON.parse(JSON.stringify(record)) as unknown);
      if (found !== undefined) return found;
    }
    return check.finish();
  };

  it('accepts HSP records within the run, in any order of their indices, with frames and fields that the ABI may add', () => {
    expect(problem([hsp(), second()])).toBeUndefined();
    expect(problem([second(), hsp()])).toBeUndefined();
    expect(problem([hsp({ query_frame: -3, subject_frame: 2, later_field: 'x' }), second({ subject_length: null })])).toBeUndefined();
    // The adapter writes a value that is not finite as null (web/adapter/src/json.rs).
    expect(problem([hsp({ bit_score: null, e_value: null }), second({ raw_score: null })])).toBeUndefined();
  });

  it('names the HSP record that the run cannot have', () => {
    expect(problem([hsp()])).toBe('there are 1 HSP records, but the manifest gives 2');
    expect(problem([hsp(), second(), hsp({ index: 2 })])).toBe('there are more HSP records than the 2 that the manifest gives');
    expect(problem([hsp(), hsp({ index: 0, rank: 1 })])).toBe('HSP record 2 has an index that is outside 0 to 1 or repeated');
    expect(problem([hsp(), second({ q_idx: 2 })])).toBe('HSP record 2 names query record 2, but the run has 2 query records');
    expect(problem([hsp({ s_idx: 2 }), second()])).toBe('HSP record 1 names subject record 2, but the run has 2 subject records');
    expect(problem([hsp(), second({ q_idx: 0 })])).toBe('HSP record 2 has a rank that another HSP of its query has');
    expect(problem([hsp({ q_start: 0 }), second()])).toBe('HSP record 1 has q_start 0, not a coordinate (a whole number of 1 or more)');
    expect(problem([hsp({ s_end: 2.5 }), second()])).toMatch(/has s_end 2\.5, not a coordinate/);
    expect(problem([hsp({ query_frame: 4 }), second()])).toBe('HSP record 1 has query_frame 4, not null or a frame of -3 to 3 other than 0');
    expect(problem([hsp({ subject_frame: 0 }), second()])).toMatch(/has subject_frame 0, not null or a frame/);
    expect(problem([hsp({ out6: [0, 31] }), second()])).toBe("HSP record 1 has out6 [0,31], not null or a byte range within the run's 30 bytes of outfmt 6");
    expect(problem([hsp({ out6: [5, 4] }), second()])).toMatch(/has out6 \[5,4\]/);
    expect(problem([hsp({ out0: [10, 13] }), second()])).toMatch(/has out0 \[10,13\], not null or a byte range within the run's 12 bytes of outfmt 0/);
    expect(problem([hsp({ out0_subject: [0, 13] }), second()])).toMatch(/has out0_subject \[0,13\]/);
  });

  it('refuses ranks that are not 0 to n - 1 within each query, which the typed table would wrap or leave unfound (code review 2 L1)', () => {
    // 2^32 becomes 0 in the table's Int32Array, giving two HSPs "rank 0".
    expect(problem([hsp({ rank: 4294967296 }), second()])).toBe('HSP record 1 has rank 4294967296, but the run has 2 HSP records');
    // A gap: the query's one HSP has rank 1; a rank of the count or more is never found.
    expect(problem([hsp({ rank: 1 }), second()])).toBe('the ranks of the HSP records of query record 0 are not 0 to 0 (there is a rank 1)');
    // A repeat.
    expect(problem([hsp({ q_idx: 1 }), second()])).toBe('HSP record 2 has a rank that another HSP of its query has');
    // Two HSPs of one query with ranks 0 and 2 (a rank of the run's count or more is refused at once).
    expect(problem([hsp(), second({ q_idx: 0, rank: 2 })])).toBe('HSP record 2 has rank 2, but the run has 2 HSP records');
    // Ranks in any order are 0 to n - 1 as a set.
    expect(problem([hsp({ rank: 1 }), second({ q_idx: 0, rank: 0 })])).toBeUndefined();
  });

  it('refuses fields of the wrong type that a typed array or a coercion would have hidden (code review M1)', () => {
    // Each would become another value in the table's typed arrays: null and true become 0 or 1,
    // 2^32 becomes 0 in an Int32Array, 259 becomes 3 in an Int8Array, "0" becomes 0.
    const cases: Array<[Record<string, unknown>, string]> = [
      [{ s_idx: null }, 'HSP record 1 has s_idx null, not a whole number of 0 or more'],
      [{ s_idx: undefined }, 'HSP record 1 has no s_idx'],
      [{ s_idx: 0.5 }, 'HSP record 1 has s_idx 0.5, not a whole number of 0 or more'],
      [{ q_idx: 4294967296 }, 'HSP record 1 names query record 4294967296, but the run has 2 query records'],
      [{ q_idx: '0' }, 'HSP record 1 has q_idx "0", not a whole number of 0 or more'],
      [{ index: -1 }, 'HSP record 1 has index -1, not a whole number of 0 or more'],
      [{ rank: 2.5 }, 'HSP record 1 has rank 2.5, not a whole number of 0 or more'],
      [{ q_start: true }, 'HSP record 1 has q_start true, not a coordinate (a whole number of 1 or more)'],
      [{ q_end: 2 ** 53 }, 'HSP record 1 has q_end 9007199254740992, not a coordinate (a whole number of 1 or more)'],
      [{ query_frame: 259 }, 'HSP record 1 has query_frame 259, not null or a frame of -3 to 3 other than 0'],
      [{ subject_frame: '1' }, 'HSP record 1 has subject_frame "1", not null or a frame of -3 to 3 other than 0'],
      [{ bit_score: true }, 'HSP record 1 has bit_score true, not null or a number'],
      [{ e_value: '1e-5' }, 'HSP record 1 has e_value "1e-5", not null or a number'],
      [{ subject_length: -1 }, 'HSP record 1 has subject_length -1, not null or a whole number of 0 or more'],
      [{ query_aligned: 5 }, 'HSP record 1 has query_aligned 5, not null or text'],
      [{ out6: [0, '15'] }, 'HSP record 1 has out6 [0,"15"], not null or a byte range within the run\'s 30 bytes of outfmt 6'],
      [{ out0: [0, 6, 7] }, "HSP record 1 has out0 [0,6,7], not null or a byte range within the run's 12 bytes of outfmt 0"],
      [{ out0_subject: 'x'.repeat(100) }, `HSP record 1 has out0_subject "${'x'.repeat(39)}…, not null or a byte range within the run's 12 bytes of outfmt 0`],
    ];
    for (const [change, message] of cases) expect(problem([hsp(change), second()])).toBe(message);
    expect(problem([[0], second()])).toBe('HSP record 1 is not a JSON object');
    expect(problem([null, second()])).toBe('HSP record 1 is not a JSON object');
    // JSON numbers too large for a double are Infinity once parsed.
    const check = new HspRecordCheck(hspRecordBounds(r));
    expect(check.next(JSON.parse(JSON.stringify(hsp()).replace('"raw_score":10', '"raw_score":1e999')))).toBe('HSP record 1 has raw_score Infinity, not null or a number');
  });

  it('names the first record of chosen files that differs from the saved input, or the count', () => {
    const saved = r.subject;
    const same = saved.records.id.map((id, k) => ({ id, length: saved.records.length[k]!, sha256: saved.records.sha256[k]! }));
    expect(recordsMismatch(saved, same)).toBeUndefined();
    expect(recordsMismatch(saved, [same[0]!, { ...same[1]!, sha256: hex('9') }])).toBe(
      'record 2 ("s3") differs from the saved run\'s record: its bytes have another SHA-256',
    );
    expect(recordsMismatch(saved, [{ ...same[0]!, id: 's9' }, same[1]!])).toBe(
      'record 1 is "s1" (length 10) in the saved run, but the chosen files give "s9" (length 10)',
    );
    expect(recordsMismatch(saved, [...same, same[0]!])).toBe('the chosen files give 3 records after the exclusions, but the saved run searched 2');
    expect(sourceMismatch(saved.sources[0]!, 0, 'a.fa', 2)).toBeUndefined();
    expect(sourceMismatch(saved.sources[0]!, 0, 'a2.fa', 1)).toBe('"a2.fa" has 1 record, but file 1 of the saved input ("a.fa") had 2');
  });

  it('matches chosen files to the recorded sources by their records, whatever order they were chosen in (code review L3)', () => {
    const saved = r.subject;
    const record = (id: string, length: number, c: string) => ({ id, length, sha256: hex(c) });
    const a = { name: 'a.fa', records: [record('s1', 10, '3'), record('s2', 7, '4')] };
    const b = { name: 'b.fa', records: [record('s3', 12, '5')] };
    expect(matchSources(saved, [a, b])).toEqual({ ok: true, files: [0, 1] });
    expect(matchSources(saved, [b, a])).toEqual({ ok: true, files: [1, 0] });
    // The names do not matter, nor the record that the run left out.
    expect(matchSources(saved, [{ ...b, name: 'x.fa' }, { name: 'y.fa', records: [record('s1', 10, '3'), record('other', 1, '9')] }])).toEqual({
      ok: true,
      files: [1, 0],
    });
    // A file that fits two sources goes where it is needed: x fits both, y only the first.
    const twice: SessionInput = {
      ...saved,
      records: { id: ['s1', 's1'], length: [10, 10], sha256: [hex('3'), hex('3')] },
      sources: [
        { name: 'a.fa', size: 30, records: 2, excluded: [1] },
        { name: 'a.fa', size: 30, records: 2, excluded: [0] },
      ],
    };
    const x = { name: 'x.fa', records: [record('s1', 10, '3'), record('s1', 10, '3')] };
    const y = { name: 'y.fa', records: [record('s1', 10, '3'), record('s9', 10, '8')] };
    expect(matchSources(twice, [x, y])).toEqual({ ok: true, files: [1, 0] });
    // The refusals name a file's record count, or the first record of the input that differs.
    expect(matchSources(saved, [b, { name: 'a2.fa', records: [record('s1', 10, '3')] }])).toEqual({
      ok: false,
      message: '"a2.fa" has 1 record, but file 1 of the saved input ("a.fa") had 2',
    });
    expect(matchSources(saved, [b, { name: 'a.fa', records: [record('s1', 10, '9'), record('s2', 7, '4')] }])).toEqual({
      ok: false,
      message: 'record 1 ("s1") differs from the saved run\'s record: its bytes have another SHA-256',
    });
    expect(matchSources(saved, [a, { name: 'b.fa', records: [record('s4', 12, '5')] }])).toEqual({
      ok: false,
      message: 'record 2 is "s3" (length 12) in the saved run, but the chosen files give "s4" (length 12)',
    });
    expect(matchSources(twice, [y, y])).toEqual({ ok: false, message: 'record 2 is "s1" (length 10) in the saved run, but the chosen files give "s9" (length 10)' });
  });
});

describe('the paced checks of a large manifest (fix round 2)', () => {
  /** A run whose query has `n` records, as a session of `n` queries saves it. */
  function large(n: number): SessionManifest {
    const query = {
      name: 'query.fa',
      sha256: hex('a'),
      length: 20 * n,
      reader: 1 as const,
      records: {
        id: Array.from({ length: n }, (_, k) => `q${k}`),
        length: Array.from({ length: n }, () => 8),
        sha256: Array.from({ length: n }, () => hex('1')),
      },
      sources: [{ name: 'query.fa', size: 20 * n, records: n + 2, excluded: [3, n] }],
    };
    return manifest({ runs: [run({ query })] });
  }

  it('gives the manifest that checkManifest gives, frozen and written the same, pausing between its steps', async () => {
    const value = large(50_000);
    let pauses = 0;
    const paced = await checkManifestPaced(JSON.parse(JSON.stringify(value)), SESSION_LIMITS, async () => void pauses++);
    const whole = checkManifest(JSON.parse(JSON.stringify(value)));
    expect(paced).toEqual(whole);
    expect(JSON.stringify(paced)).toBe(JSON.stringify(value));
    expect(Object.isFrozen(paced.runs[0]!.query.records.id)).toBe(true);
    expect(Object.isFrozen(paced.runs[0]!.query.records.sha256)).toBe(true);
    expect(Object.isFrozen(paced.runs[0]!.query.sources[0]!.excluded)).toBe(true);
    // Steps of at most 20,000 records, and a step for each long column frozen.
    expect(pauses).toBeGreaterThanOrEqual(2 + 3);
  });

  it('refuses what checkManifest refuses, with the same message', async () => {
    const value = JSON.parse(JSON.stringify(large(30_000)));
    value.runs[0].query.records.sha256[25_000] = 'not a hash';
    const message = 'runs[0].query.records.sha256[25000] is not a SHA-256 (64 lower-case hex digits)';
    expect(() => checkManifest(JSON.parse(JSON.stringify(value)))).toThrow(message);
    await expect(checkManifestPaced(value, SESSION_LIMITS, () => Promise.resolve())).rejects.toThrow(message);
  });

  it('reads a container with pushPaced as push reads it, pausing while it checks the manifest', async () => {
    const m = large(25_000);
    const bytes = container(m, [CANDIDATES[0]!]);
    const reader = new SessionFileReader();
    let pauses = 0;
    const events: SessionEvent[] = [];
    for (let at = 0; at < bytes.length; at += 65_536) {
      for (const event of await reader.pushPaced(bytes.slice(at, at + 65_536), async () => void pauses++)) {
        if (event.type !== 'data') events.push(event);
      }
    }
    reader.finish();
    expect(events).toEqual(read(bytes, 65_536).events);
    expect(events[0]).toMatchObject({ type: 'manifest' });
    expect(pauses).toBeGreaterThanOrEqual(2);
    // A refusal fails the reader as push does.
    const damaged = new SessionFileReader();
    await expect(damaged.pushPaced(encoder.encode('LOSAT-WEB-SESSION 1\nmanifest 2\n{}\n'), () => Promise.resolve())).rejects.toThrow(SessionFileError);
    expect(() => damaged.push(new Uint8Array(1))).toThrow('The session file is damaged.');
  });
});
