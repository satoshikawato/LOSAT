import { describe, expect, it } from 'vitest';
import { fromText, leavesOut, leftOutIds, partText, worthSaying } from '../../src/domain/input-parts';

// Line WQA's D22 (2026-11-23): the words that say which records an input leaves out.
const part = (name: string, records: number, leftOut: string[], leftOutCount = leftOut.length) => ({ name, records, leftOut, leftOutCount });

describe('input parts', () => {
  it('names the IDs of the records left out, the first three', () => {
    expect(leftOutIds(['r2'])).toBe('r2');
    expect(leftOutIds(['r2', 'r5'])).toBe('r2 and r5');
    expect(leftOutIds(['r2', 'r5', 'r7'])).toBe('r2, r5 and r7');
    expect(leftOutIds(['r2', 'r5', 'r7', 'r8', 'r9'])).toBe('r2, r5, r7 and 2 more');
    // A part keeps the first IDs and their number.
    expect(leftOutIds(['r2', 'r5', 'r7'], 1_000_000)).toBe('r2, r5, r7 and 999,997 more');
  });

  it('says the records that a file holds of those it has', () => {
    expect(partText(part('q.fa', 2, ['r2']))).toBe('q.fa without r2 (2 of 3 records)');
    expect(partText(part('q.fa.gz', 1, ['r1 (#2)']))).toBe('q.fa.gz without r1 (#2) (1 of 2 records)');
    expect(partText(part('big.fa', 10, ['x1', 'x2', 'x3'], 1_000_000))).toBe('big.fa without x1, x2, x3 and 999,997 more (10 of 1,000,010 records)');
    expect(partText(part('b.fa', 3, []))).toBe('b.fa (3 records)');
    expect(partText(part('pasted text', 1, []))).toBe('pasted text (1 record)');
  });

  it('tells what an input is made of', () => {
    const parts = [
      part('a.fa', 2, []),
      part('b.fa', 1, ['b2', 'b3']),
    ];
    expect(fromText(parts)).toBe('a.fa (2 records) and b.fa without b2 and b3 (1 of 3 records)');
    expect(leavesOut(parts)).toBe(true);
    expect(leavesOut(undefined)).toBe(false);
  });

  it('has something to tell for records left out or several parts, not for pasted text alone', () => {
    expect(worthSaying(undefined)).toBe(false);
    expect(worthSaying([part('pasted text', 3, [])])).toBe(false);
    expect(worthSaying([part('pasted text', 2, ['r2'])])).toBe(true);
    expect(worthSaying([part('a.fa', 1, []), part('b.fa', 1, [])])).toBe(true);
  });
});
