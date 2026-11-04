// The UI's display helpers (src/ui/format.ts).
import { describe, expect, it } from 'vitest';
import {
  formatCounted,
  formatDateTime,
  formatRuns,
  framesLabel,
  framesPhrase,
  framesText,
  lostResultsText,
  optionHelp,
} from '../../src/ui/format';

describe('formatDateTime', () => {
  it('writes a local time in ISO 8601 form, day and month never swapped', () => {
    expect(formatDateTime(new Date(2026, 9, 9, 21, 49, 30).getTime())).toBe('2026-10-09 21:49:30');
    expect(formatDateTime(new Date(2026, 0, 2, 3, 4, 5).getTime())).toBe('2026-01-02 03:04:05');
  });

  it('pads every field and drops the milliseconds', () => {
    expect(formatDateTime(new Date(999, 11, 31, 0, 0, 0, 999).getTime())).toBe('0999-12-31 00:00:00');
  });
});

describe('optionHelp', () => {
  it("keeps the help's first paragraph without the references to NCBI's sources", () => {
    expect(optionHelp('Cost to open a gap [default: 11]\n\nMore text')).toBe('Cost to open a gap [default: 11]');
    expect(optionHelp('Reward for a match. NCBI reference: blast_options.c:12')).toBe('Reward for a match.');
    expect(optionHelp(undefined)).toBe('');
  });
});

describe('formatCounted', () => {
  it('writes the noun in the singular for one only', () => {
    expect(formatCounted(1, 'subject')).toBe('1 subject');
    expect(formatCounted(3, 'subject')).toBe('3 subjects');
    expect(formatCounted(1, 'HSP')).toBe('1 HSP');
    expect(formatCounted(12345, 'HSP')).toBe('12,345 HSPs');
    expect(formatCounted(0, 'HSP')).toBe('0 HSPs');
  });
});

describe('formatRuns', () => {
  it('names one run, two runs, and a list with "and"', () => {
    expect(formatRuns([2])).toBe('run 2');
    expect(formatRuns([3, 1])).toBe('runs 1 and 3');
    expect(formatRuns([1, 2])).toBe('runs 1 and 2');
    expect(formatRuns([5, 1, 3])).toBe('runs 1, 3 and 5');
  });

  it('writes three or more runs in a row as a range, so that hundreds of runs stay short', () => {
    expect(formatRuns([1, 2, 3])).toBe('runs 1–3');
    expect(formatRuns([1, 2, 3, 5])).toBe('runs 1–3 and 5');
    expect(formatRuns([7, 1, 3, 5, 6, 3])).toBe('runs 1, 3 and 5–7');
    expect(formatRuns(Array.from({ length: 300 }, (_, i) => i + 1))).toBe('runs 1–300');
  });
});

describe('lostResultsText', () => {
  it('is empty when every stored result is there', () => {
    expect(lostResultsText([], [])).toBe('');
  });

  it('names the runs the browser removed the stored results of, with "it" for one run and "them" for more', () => {
    expect(lostResultsText([1], [])).toBe('The browser removed the stored results of run 1. Run it again.');
    expect(lostResultsText([1, 2, 3], [])).toBe('The browser removed the stored results of runs 1–3. Run them again.');
  });

  it('names the runs whose stored results could not be read otherwise', () => {
    expect(lostResultsText([], [4])).toBe('The stored results of run 4 could not be read. Run it again.');
    expect(lostResultsText([], [4, 6])).toBe('The stored results of runs 4 and 6 could not be read. Run them again.');
  });

  it('says both, removed first, and "them" for one run of each', () => {
    expect(lostResultsText([2], [4])).toBe(
      'The browser removed the stored results of run 2. The stored results of run 4 could not be read. Run them again.',
    );
  });
});

describe('frames', () => {
  const tblastn = { subjectFrame: 2 };
  const blastx = { queryFrame: -1 };
  const tblastx = { queryFrame: -2, subjectFrame: 2 };

  it('names only the translated sequence’s frame where one sequence is translated (W4b screen review L9)', () => {
    expect(framesLabel(tblastn)).toBe('Subject frame');
    expect(framesLabel(blastx)).toBe('Query frame');
    expect(framesLabel(tblastx)).toBe('Frames (q/s)');
    expect(framesLabel({})).toBeUndefined();
  });

  it('writes the frames with their signs, the query’s first', () => {
    expect(framesText(tblastn)).toBe('+2');
    expect(framesText(blastx)).toBe('-1');
    expect(framesText(tblastx)).toBe('-2/+2');
    expect(framesText(tblastx, ' / ')).toBe('-2 / +2');
    expect(framesText({})).toBe('');
  });

  it('puts them in a sentence', () => {
    expect(framesPhrase(tblastn)).toBe('subject frame +2');
    expect(framesPhrase(blastx)).toBe('query frame -1');
    expect(framesPhrase(tblastx)).toBe('frames -2 / +2');
    expect(framesPhrase({})).toBe('');
  });
});
