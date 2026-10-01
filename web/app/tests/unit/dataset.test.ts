import { describe, expect, it } from 'vitest';
import {
  includedRecords,
  normalizeExclusion,
  recordMismatch,
  type DatasetRecord,
  type DatasetRevision,
} from '../../src/domain/dataset';

const record = (index: number, id: string): DatasetRecord => ({
  index,
  id,
  header_offset: index * 10,
  sequence_offset: index * 10 + 3,
  end_offset: index * 10 + 10,
  length: 6,
  line_layout: { kind: 'uniform', width: 6, eol: 1 },
  residue_counts: { A: 6 },
  sha256: `hash-${id}`,
});

describe('dataset revisions', () => {
  const revision: DatasetRevision = {
    revisionId: 'r',
    sourceId: 's',
    parser: 0,
    records: [record(0, 'a'), record(1, 'b'), record(2, 'c')],
    excluded: [1],
  };

  it('lists the included records in the original order', () => {
    expect(includedRecords(revision).map((r) => r.id)).toEqual(['a', 'c']);
  });

  it('sorts and deduplicates an exclusion and rejects indices outside the table', () => {
    expect(normalizeExclusion([2, 0, 2], 3)).toEqual([0, 2]);
    expect(() => normalizeExclusion([3], 3)).toThrow(RangeError);
    expect(() => normalizeExclusion([-1], 3)).toThrow(RangeError);
    expect(() => normalizeExclusion([0.5], 3)).toThrow(RangeError);
  });
});

describe('record tables at register', () => {
  const table = [
    { id: 'q1', length: 8 },
    { id: 'q2', length: 4 },
  ];

  it('agrees when every ID and length is the same', () => {
    expect(recordMismatch(table, [...table])).toBeUndefined();
  });

  it('names the first record that differs', () => {
    expect(recordMismatch(table, [table[0]!, { id: 'q2', length: 5 }])).toBe(
      'record 2 is "q2" (length 4) in the record table, but the engine read "q2" (length 5)',
    );
    expect(recordMismatch(table, [{ id: 'x', length: 8 }, table[1]!])).toContain('record 1 is "q1"');
  });

  it('reports a different number of records', () => {
    expect(recordMismatch(table, table.slice(0, 1))).toBe('the record table has 2 records, but the engine read 1');
  });
});
