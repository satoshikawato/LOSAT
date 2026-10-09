// Deterministic synthetic FASTA for the results screen's measurements (results-measure.spec.ts).
// Nothing here is a search result: these are only the inputs of real searches.

/** mulberry32: a small generator with a full-length stream (a linear congruential generator's bits repeat within 2^18 draws). */
class Random {
  private state: number;
  constructor(seed: number) {
    this.state = seed >>> 0;
  }
  next(): number {
    this.state = (this.state + 0x6d2b79f5) >>> 0;
    let t = this.state;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return (t ^ (t >>> 14)) >>> 0;
  }
  below(n: number): number {
    return this.next() % n;
  }
}

const dna = (random: Random, length: number): string => {
  let out = '';
  for (let i = 0; i < length; i++) out += 'ACGT'[random.below(4)];
  return out;
};

/** The sequence with about one letter in a hundred replaced by another. */
const mutate = (random: Random, sequence: string): string =>
  [...sequence].map((letter) => (random.below(100) === 0 ? 'ACGT'[('ACGT'.indexOf(letter) + 1 + random.below(3)) % 4]! : letter)).join('');

const COMPLEMENT: Readonly<Record<string, string>> = {
  A: 'T',
  C: 'G',
  G: 'C',
  T: 'A',
};
const reverseComplement = (sequence: string): string =>
  [...sequence]
    .reverse()
    .map((letter) => COMPLEMENT[letter]!)
    .join('');

export const QUERY_LENGTH = 100;
/** The kinds of query of `manyQueries`, by index. */
export type QueryKind = 'unique' | 'shared' | 'wide' | 'random';
export const WIDE_SUBJECTS = 200;
const BIG_SUBJECTS = 3;
const BIG_LENGTH = 60_000;
const SHARED_LENGTH = 4_000;
const WIDE_BLOCK = 120;

export const queryKind = (index: number): QueryKind =>
  index % 1000 === 500 ? 'wide' : index % 10 === 9 ? 'random' : index % 10 >= 7 ? 'shared' : 'unique';
export const queryId = (index: number): string => `${queryKind(index) === 'wide' ? 'wide' : 'q'}${String(index).padStart(6, '0')}`;

export interface ManyQueries {
  readonly queries: Buffer;
  readonly subjects: Buffer;
  /** How many queries of each kind. */
  readonly kinds: Readonly<Record<QueryKind, number>>;
}

/**
 * `count` BLASTN queries of 100 letters, cut from the subjects (every fifth on the other
 * strand) or random, and 3 + 200 subject records: 70% of the queries come from one of three
 * long subjects (one subject, one HSP), 20% from a block that the three long subjects share
 * with about 1% differences (three subjects), every thousandth ("wide") from a block that 200
 * short subjects share (200 subjects), and 10% are random (no hits).
 */
export function manyQueries(count: number): ManyQueries {
  const random = new Random(20261008);
  const shared = dna(random, SHARED_LENGTH);
  const wideBlock = dna(random, WIDE_BLOCK);
  const unique: string[] = [];
  const subjects: string[] = [];
  for (let k = 0; k < BIG_SUBJECTS; k++) {
    const left = dna(random, BIG_LENGTH / 2);
    const right = dna(random, BIG_LENGTH / 2);
    unique.push(left, right);
    subjects.push(`>big${k + 1} long subject ${k + 1}\n${left}${mutate(random, shared)}${right}\n`);
  }
  for (let k = 0; k < WIDE_SUBJECTS; k++) {
    subjects.push(
      `>short${String(k + 1).padStart(3, '0')} short subject ${k + 1}\n${dna(random, 90)}${mutate(random, wideBlock)}${dna(random, 90)}\n`,
    );
  }
  const kinds: Record<QueryKind, number> = {
    unique: 0,
    shared: 0,
    wide: 0,
    random: 0,
  };
  const parts: string[] = [];
  for (let i = 0; i < count; i++) {
    const kind = queryKind(i);
    kinds[kind]++;
    let sequence: string;
    if (kind === 'random') sequence = dna(random, QUERY_LENGTH);
    else if (kind === 'wide') {
      const start = random.below(WIDE_BLOCK - QUERY_LENGTH + 1);
      sequence = wideBlock.slice(start, start + QUERY_LENGTH);
    } else if (kind === 'shared') {
      const start = random.below(SHARED_LENGTH - QUERY_LENGTH + 1);
      sequence = shared.slice(start, start + QUERY_LENGTH);
    } else {
      const source = unique[random.below(unique.length)]!;
      const start = random.below(source.length - QUERY_LENGTH + 1);
      sequence = source.slice(start, start + QUERY_LENGTH);
    }
    if (kind !== 'random' && i % 5 === 2) sequence = reverseComplement(sequence);
    parts.push(`>${queryId(i)} synthetic ${kind} query ${i + 1}\n${sequence}\n`);
  }
  return {
    queries: Buffer.from(parts.join('')),
    subjects: Buffer.from(subjects.join('')),
    kinds,
  };
}

/**
 * One record that repeats a unit of `unit` letters `copies` times, with about 1% differences
 * between the copies: searched against itself, every pair of copies gives a diagonal of the
 * dot plot, so one query-subject pair has many HSPs.
 */
export function repeats(unit: number, copies: number): Buffer {
  const random = new Random(77);
  const base = dna(random, unit);
  let sequence = '';
  for (let i = 0; i < copies; i++) sequence += mutate(random, base);
  return Buffer.from(`>repeat unit ${unit} x ${copies}\n${sequence.replace(/(.{70})/g, '$1\n').replace(/\n$/, '')}\n`);
}
