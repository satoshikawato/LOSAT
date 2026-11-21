import { describe, expect, it } from 'vitest';
import { recordMismatch } from 'losat/domain';
import { includedRecords, normalizeExclusion, recordName, type DatasetRecord, type DatasetRevision } from '../../src/domain/dataset';

const record = (index: number, id: string): DatasetRecord => ({
  index,
  id,
  header_offset: index * 10,
  sequence_offset: index * 10 + 3,
  end_offset: index * 10 + 10,
  length: 6,
  residue_counts: { A: 6 },
  sha256: `hash-${id}`,
});

describe('dataset revisions', () => {
  const revision: DatasetRevision = {
    revisionId: 'r',
    sourceId: 's',
    parser: 1,
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

describe('record names in the results', () => {
  it("names a record without a title as NCBI's tabular outputs do", () => {
    expect(recordName('q1', 'query', 0)).toBe('q1');
    expect(recordName('', 'query', 0)).toBe('Query_1');
    expect(recordName('', 'subject', 1)).toBe('Subject_2');
  });

  // line WQA's sweep, D7, 2026-11-19: a run of -parse_deflines names its records as its outputs do.
  it('names a record by the ID that the outputs show, where the engine read the ids of the deflines', () => {
    const shown = { query: ['AB000001.1'], subject: ['s1', 'NP_000006.1'] };
    expect(recordName('gb|AB000001.1|', 'query', 0, shown)).toBe('AB000001.1');
    expect(recordName('gi|6|ref|NP_000006.1|', 'subject', 1, shown)).toBe('NP_000006.1');
    expect(recordName('gb|AB000001.1|', 'query', 0)).toBe('gb|AB000001.1|');
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
