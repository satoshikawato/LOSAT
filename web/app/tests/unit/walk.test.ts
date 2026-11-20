// The pure parts of the random walk (tests/e2e/random-walk.spec.ts): its plan covers what it
// says, its shrinking finds a minimal reproduction within its budget, its model of the inputs
// joins records as the data layer does, and its inputs have hits in every mode.
import { describe, expect, it } from 'vitest';
import { alphabetOf, makeFixtures, pickFor, XSS_PAYLOAD } from '../e2e/walk/fixtures';
import {
  emptyModel,
  expectedInputs,
  includedBytes,
  joinedBytes,
  parseShellWords,
  recordRanges,
  withFiles,
  withIncluded,
  withMode,
  withPaste,
  type ModelSource,
} from '../e2e/walk/model';
import { BARRIER_KINDS, coverageOf, coverageTargets, createRandom, OPENING, planWalk, PR_SEED, stepCount, WALK_KINDS } from '../e2e/walk/plan';
import { pairPhaseBudget, shrinkWalk } from '../e2e/walk/shrink';
import { shellWords } from '../../src/domain/argv';

describe('walk plan', () => {
  it('covers every kind and every ordered pair at the pull-request level, in a few walks', () => {
    const plan = planWalk(PR_SEED, { level: 'pr' });
    expect(coverageOf(plan, 'pr').missing).toEqual([]);
    expect(plan.walks.length).toBeLessThanOrEqual(3);
    expect(stepCount(plan)).toBeLessThanOrEqual(90);
  });

  it('covers every ordered pair as adjacent steps at the full level', () => {
    const plan = planWalk(PR_SEED, { level: 'full' });
    expect(coverageOf(plan, 'full').missing).toEqual([]);
    expect(coverageTargets('full').adjacent!.size).toBe(WALK_KINDS.length ** 2);
  });

  it('starts every walk by running the inputs it starts with', () => {
    for (const level of ['pr', 'full'] as const) for (const walk of planWalk(PR_SEED, { level }).walks) expect(walk.steps[0]!.kind).toBe(OPENING);
  });

  it('is the same for the same seed, and differs for another', () => {
    expect(planWalk(7, { level: 'pr' })).toEqual(planWalk(7, { level: 'pr' }));
    expect(planWalk(7, { level: 'pr' })).not.toEqual(planWalk(8, { level: 'pr' }));
  });

  it('counts no pair across a barrier', () => {
    const [barrier] = BARRIER_KINDS;
    const plan = { walks: [{ steps: [{ kind: 'mode' as const }, { kind: barrier! }, { kind: 'run' as const }] }] };
    const { missing } = coverageOf(plan, 'pr');
    expect(missing).toContain('sequence mode>run');
    expect(missing).not.toContain(`sequence mode>${barrier}`);
    expect(missing).not.toContain(`sequence ${barrier}>run`);
  });

  it('draws numbers in [0, 1) from a seed', () => {
    const a = createRandom(1);
    const b = createRandom(1);
    const draws = Array.from({ length: 1000 }, () => a.next());
    expect(draws).toEqual(Array.from({ length: 1000 }, () => b.next()));
    expect(Math.min(...draws)).toBeGreaterThanOrEqual(0);
    expect(Math.max(...draws)).toBeLessThan(1);
  });
});

describe('walk shrink', () => {
  /** A predicate that reproduces when every step of `needed` is in the candidate, and counts its calls. */
  const needs = (needed: readonly number[]) => {
    const calls: number[][] = [];
    return {
      calls,
      reproduces: async (indices: readonly number[]) => {
        calls.push([...indices]);
        return needed.every((index) => indices.includes(index));
      },
    };
  };
  const prefix = [...Array(20).keys()];

  it('finds the failing step alone', async () => {
    const { reproduces, calls } = needs([20]);
    const result = await shrinkWalk({ prefix, failing: 20, budget: 10, reproduces });
    expect(result).toEqual({ indices: [20], replays: 1, finished: true, confirmed: true });
    expect(calls).toEqual([[20]]);
  });

  it('finds a pair, nearest step first', async () => {
    const { reproduces } = needs([17, 20]);
    const result = await shrinkWalk({ prefix, failing: 20, budget: 10, reproduces });
    expect(result.indices).toEqual([17, 20]);
    expect(result.replays).toBe(4);
  });

  it('finds three steps with ddmin after the pairs', async () => {
    const { reproduces } = needs([2, 9, 20]);
    const result = await shrinkWalk({ prefix, failing: 20, budget: 60, reproduces });
    expect(result.indices).toEqual([2, 9, 20]);
    expect(result.confirmed).toBe(true);
  });

  it('stops at its budget with the shortest reproduction found', async () => {
    const { reproduces, calls } = needs([2, 9, 20]);
    const result = await shrinkWalk({ prefix, failing: 20, budget: 6, reproduces });
    expect(calls.length).toBe(6);
    expect(result.finished).toBe(false);
    expect(result.indices.at(-1)).toBe(20);
    expect(pairPhaseBudget(6)).toBe(3);
  });

  it('starts no replay after its deadline', async () => {
    const { reproduces, calls } = needs([2, 9, 20]);
    const result = await shrinkWalk({ prefix, failing: 20, budget: 60, reproduces, deadline: Date.now() - 1 });
    expect(calls.length).toBe(0);
    expect(result).toEqual({ indices: [...prefix, 20], replays: 0, finished: false, confirmed: false });
  });

  it('keeps the whole walk when nothing smaller reproduces', async () => {
    const result = await shrinkWalk({ prefix: [0, 1, 2], failing: 3, budget: 20, reproduces: async (indices) => indices.length === 4 });
    expect(result.indices).toEqual([0, 1, 2, 3]);
    expect(result.confirmed).toBe(false);
  });
});

describe('walk model', () => {
  const text = (value: string) => new TextEncoder().encode(value);
  const decode = (bytes: Uint8Array) => new TextDecoder().decode(bytes);

  it('cuts records from one ">" at a line start to the next', () => {
    const bytes = text('>a x>y\nAC\n>b\nGT\n>c\nTT');
    expect(recordRanges(bytes)).toEqual([
      [0, 10],
      [10, 16],
      [16, 21],
    ]);
  });

  it('searches the included records, and joins inputs with a line end where one has none', () => {
    const source = (value: string, excluded: number[] = []): ModelSource => ({ origin: 'file', name: 'f', bytes: text(value), excluded: { nucleotide: new Set(excluded) } });
    expect(decode(includedBytes(source('>a\nAC\n>b\nGT\n>c\nTT\n', [1]), 'nucleotide'))).toBe('>a\nAC\n>c\nTT\n');
    // Another reading of the same source has its own records left out (none here).
    expect(decode(includedBytes(source('>a\nAC\n>b\nGT\n', [1]), 'protein'))).toBe('>a\nAC\n>b\nGT\n');
    expect(decode(joinedBytes([source('>a\nAC'), source('>b\nGT\n'), source('>c\nTT\n', [0]), source('>d\nCC\n')], 'nucleotide'))).toBe('>a\nAC\n>b\nGT\n>d\nCC\n');
  });

  it('names the inputs as the command line does: the pasted text, one file, or several together or apart', () => {
    let model = withPaste(emptyModel(), 'query', '>q\nAC\n');
    expect(expectedInputs(model.query, 'query', 'nucleotide', [0]).map((input) => input.name)).toEqual(['query.fa']);
    model = withFiles(model, 'query', [{ name: 'a.fa', bytes: text('>a\nGG\n') }]);
    expect(expectedInputs(model.query, 'query', 'nucleotide', [0, 1]).map((input) => input.name)).toEqual(['combined_query.fa']);
    expect(expectedInputs(model.query, 'query', 'nucleotide', [1]).map((input) => input.name)).toEqual(['a.fa']);
    model = withMode(model, 'query', 'separate');
    expect(expectedInputs(model.query, 'query', 'nucleotide', [0, 1]).map((input) => input.name)).toEqual(['query.fa', 'a.fa']);
    // Apart, an input whose records are all left out is no run of its own.
    const apart = withIncluded(model, 'query', 'nucleotide', 1, 0, false);
    expect(expectedInputs(apart.query, 'query', 'nucleotide', [0, 1]).map((input) => input.name)).toEqual(['query.fa']);
    expect(expectedInputs(withMode(apart, 'query', 'combined').query, 'query', 'nucleotide', [0, 1]).map((input) => input.name)).toEqual(['combined_query.fa']);
    model = withIncluded(withPaste(model, 'query', '>q\nAC\n>r\nGT\n'), 'query', 'nucleotide', 0, 0, false);
    expect(decode(expectedInputs(model.query, 'query', 'nucleotide', [0])[0]!.bytes)).toBe('>r\nGT\n');
    expect(withPaste(model, 'query', '').query.sources.map((source) => source.name)).toEqual(['a.fa']);
  });

  it('reads a pasted text as its text box gives it: line ends as LF', () => {
    const model = withPaste(emptyModel(), 'query', '>q\r\nAC\r\nGT\r>r\nTT\n');
    expect(decode(model.query.sources[0]!.bytes)).toBe('>q\nAC\nGT\n>r\nTT\n');
  });

  it('reads back the words that the page quotes', () => {
    const words = ['losat', 'blastn', '-query', "it's.faa", '-outfmt', '6 qseqid sseqid', '-subject', '<svg onload="x">.fa', 'a\\b'];
    expect(parseShellWords(shellWords(words))).toEqual(words);
  });
});

describe('walk fixtures', () => {
  const fixtures = makeFixtures(PR_SEED);

  it('has texts and files of both alphabets, and a defline that carries a script', () => {
    for (const alphabet of ['nucleotide', 'protein'] as const) {
      expect(fixtures.texts.some((item) => item.alphabet === alphabet)).toBe(true);
      expect(fixtures.files.some((item) => item.alphabet === alphabet)).toBe(true);
    }
    expect(fixtures.texts.some((item) => item.text.includes(XSS_PAYLOAD))).toBe(true);
    expect(makeFixtures(PR_SEED)).toEqual(fixtures);
  });

  it('writes nucleotide records in ACGT and proteins without stops', () => {
    for (const item of [...fixtures.texts, ...fixtures.files]) {
      const sequence = item.text
        .split('\n')
        .filter((line) => !line.startsWith('>'))
        .join('')
        .toUpperCase();
      expect(sequence).toMatch(item.alphabet === 'nucleotide' ? /^[ACGT]+$/ : /^[ACDEFGHIKLMNPQRSTVWY]+$/);
    }
  });

  it('picks the mode alphabet for a low draw, another or an odd text for a high one', () => {
    expect(alphabetOf('blastx', 'query')).toBe('nucleotide');
    expect(alphabetOf('blastx', 'subject')).toBe('protein');
    expect(pickFor(fixtures.texts, fixtures.odd, 'protein', 0.1, 0.5).alphabet).toBe('protein');
    expect(pickFor(fixtures.texts, fixtures.odd, 'protein', 0.75, 0.5).alphabet).toBe('nucleotide');
    expect(fixtures.odd).toContain(pickFor(fixtures.texts, fixtures.odd, 'protein', 0.9, 0.5));
  });
});
