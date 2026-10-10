import { describe, expect, it } from 'vitest';
import {
  includedRecords,
  isFirstLineRefusal,
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

describe('isFirstLineRefusal', () => {
  it('is true for the refusal of a first line that may be a sequence identifier, and for no other refusal', () => {
    // The adapter's wording (web/adapter/src/scan/ncbi.rs), also with a prefix and an odd first line.
    const refusal = (line: string) =>
      `the first line (${JSON.stringify(line)}) is not a defline and may be a sequence identifier that NCBI BLAST+ fetches through a data loader (from GenBank or a BLAST database), which is not supported by LOSAT Web (start the input with a '>' defline)`;
    expect(isFirstLineRefusal(refusal('AB123456'))).toBe(true);
    expect(isFirstLineRefusal(`BLAST query error: ${refusal('lcl|a)b')}`)).toBe(true);
    // A gap line and NCBI's "Near line N" refusal are not fixed by a defline in front of the text.
    expect(isFirstLineRefusal("line 2 is a gap line ('>?'), which NCBI BLAST+ reads ... not supported by LOSAT Web")).toBe(false);
    expect(isFirstLineRefusal("BLAST query error: CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.")).toBe(false);
    expect(isFirstLineRefusal('')).toBe(false);
  });
});
