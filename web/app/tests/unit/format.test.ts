// The UI's display helpers (src/ui/format.ts).
import { describe, expect, it } from 'vitest';
import { MEMORY_RESULTS_CAPACITY } from '../../src/ports/data';
import { MEMORY_FULL_MESSAGE } from '../../src/infra/data/block-store';
import { NBSP, withUnit } from '../../src/domain/units';
import {
  formatBytes,
  formatCounted,
  formatDateTime,
  formatRuns,
  framesLabel,
  framesPhrase,
  framesText,
  formatMillis,
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
  it("keeps the help's first paragraph", () => {
    expect(optionHelp('Cost to open a gap [default: 11]\n\nMore text')).toBe('Cost to open a gap [default: 11]');
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

// A number and its unit never part at a line break (line WQA's sweep, D9, 2026-11-19): the storage
// note broke "at most 512" / "MB" at 390 px. Every size and unit text of the UI joins them with U+00A0.
describe('a number and its unit', () => {
  it('are joined by a no-break space', () => {
    expect(withUnit(512, 'MB')).toBe(`512${NBSP}MB`);
    expect(withUnit('1,200', 'aa')).toBe(`1,200${NBSP}aa`);
  });

  it('are joined by it in byte sizes, whatever the unit', () => {
    expect(formatBytes(111)).toBe(`111${NBSP}bytes`);
    expect(formatBytes(981_200)).toBe(`981.2${NBSP}kB`);
    expect(formatBytes(36_800_000)).toBe(`36.8${NBSP}MB`);
    expect(formatBytes(3_400_000_000)).toBe(`3.4${NBSP}GB`);
  });

  it('are joined by it in the time spans of a run', () => {
    expect(formatMillis(420)).toBe(`420${NBSP}ms`);
    expect(formatMillis(2600)).toBe(`2.6${NBSP}s`);
  });

  it('are joined by it in the storage note and its message, and have no plain space between them', () => {
    expect(MEMORY_RESULTS_CAPACITY).toBe(`512${NBSP}MB`);
    expect(MEMORY_FULL_MESSAGE).toContain(`at most 512${NBSP}MB in a tab`);
    for (const text of [formatBytes(5_000), formatMillis(5), MEMORY_RESULTS_CAPACITY, MEMORY_FULL_MESSAGE]) {
      expect(text).not.toMatch(/\d (bytes|kB|MB|GB|ms|s)\b/);
    }
  });
});
