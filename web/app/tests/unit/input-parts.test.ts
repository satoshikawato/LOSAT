import { describe, expect, it } from 'vitest';
import { fromText, inputLabel, leavesOut, leftOutIds, partText, readDifferenceText, worthSaying } from '../../src/domain/input-parts';

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

  it('has something to tell for every file of LOSAT\'s (the name is the run\'s own), not for the user\'s file as given', () => {
    expect(worthSaying(undefined)).toBe(false);
    expect(worthSaying([part('pasted text', 3, [])])).toBe(true);
    expect(worthSaying([part('pasted text', 2, ['r2'])])).toBe(true);
    expect(worthSaying([part('a.fa', 1, []), part('b.fa', 1, [])])).toBe(true);
  });

  it('calls a file of LOSAT\'s by what it is made of in a list of runs, and the user\'s file by its name', () => {
    expect(inputLabel({ name: 'q.fa' })).toBe('q.fa');
    expect(inputLabel({ name: 'losat-run1-blastn-20261123T101530Z.query.fa', from: [part('pasted text', 1, [])] })).toBe('pasted text');
    expect(inputLabel({ name: 'x', from: [part('a.fa', 2, []), part('b.fa', 1, [])] })).toBe('a.fa + b.fa');
    expect(inputLabel({ name: 'x', from: [part('q.fa', 2, ['r2'])] })).toBe('q.fa, some records left out');
  });

  it('explains when reading an assembled input as one file differs from the per-file tables (D17)', () => {
    const table = [
      { id: 'q1', length: 180 },
      { id: 'q3', length: 180 },
    ];
    expect(readDifferenceText(table, undefined, 2)).toBeUndefined();
    expect(readDifferenceText(table, [...table], 2)).toBeUndefined();
    // Combined files with different record count:
    expect(
      readDifferenceText(
        Array.from({ length: 5 }, (_, i) => ({ id: `q${i + 1}`, length: 180 })),
        Array.from({ length: 6 }, (_, i) => ({ id: `q${i + 1}`, length: 180 })),
        2,
      ),
    ).toBe('Read as one file, as the command line reads it, these give 6 records (5 when each file is read alone).');
    // Single-file selection with different record count:
    expect(readDifferenceText(table, [table[0]!], 1)).toBe(
      'Read as one file, as the command line reads it, this selection gives 1 record (2 in the record table).',
    );
    // Single-file selection with same count, one record length differs:
    expect(readDifferenceText(table, [table[0]!, { id: 'q3', length: 120 }], 1)).toBe(
      'Read as one file, as the command line reads it, 1 record differs from the record table (record 2 (q3): length 120 instead of 180).',
    );
    // Combined files with same count, multiple records differ (including ID):
    expect(
      readDifferenceText(
        table,
        [
          { id: 'q1cr', length: 120 },
          { id: 'q3', length: 150 },
        ],
        2,
      ),
    ).toBe(
      'Read as one file, as the command line reads it, 2 records differ from reading each file alone (first: record 1: "q1cr" (length 120) instead of "q1" (length 180)).',
    );
  });
});
