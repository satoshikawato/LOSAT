import { describe, expect, it } from 'vitest';
import { fromText, leftOutIds, leftOutTexts, partText } from '../../src/domain/input-parts';

// Line WQA's D22 (2026-11-23): the words that say which records an input leaves out.
describe('input parts', () => {
  it('names the IDs of the records left out, the first three', () => {
    expect(leftOutIds(['r2'])).toBe('r2');
    expect(leftOutIds(['r2', 'r5'])).toBe('r2 and r5');
    expect(leftOutIds(['r2', 'r5', 'r7'])).toBe('r2, r5 and r7');
    expect(leftOutIds(['r2', 'r5', 'r7', 'r8', 'r9'])).toBe('r2, r5, r7 and 2 more');
  });

  it('says the records that a file holds of those it has', () => {
    expect(partText({ name: 'q.fa', records: 2, leftOut: ['r2'] })).toBe('q.fa: 2 of 3 records (r2 left out)');
    expect(partText({ name: 'q.fa', records: 1, leftOut: ['r2'] })).toBe('q.fa: 1 of 2 records (r2 left out)');
    expect(partText({ name: 'b.fa', records: 3, leftOut: [] })).toBe('b.fa: 3 records');
    expect(partText({ name: 'pasted text', records: 1, leftOut: [] })).toBe('pasted text: 1 record');
  });

  it('tells what an input is made of, and only the parts that leave records out as left out', () => {
    const parts = [
      { name: 'a.fa', records: 2, leftOut: [] },
      { name: 'b.fa', records: 1, leftOut: ['b2', 'b3'] },
    ];
    expect(fromText(parts)).toBe('a.fa: 2 records; b.fa: 1 of 3 records (b2 and b3 left out)');
    expect(leftOutTexts(parts)).toEqual(['b.fa: 1 of 3 records (b2 and b3 left out)']);
    expect(leftOutTexts(undefined)).toEqual([]);
  });
});
