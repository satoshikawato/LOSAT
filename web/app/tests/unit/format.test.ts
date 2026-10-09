// The UI's display helpers (src/ui/format.ts).
import { describe, expect, it } from 'vitest';
import { formatDateTime, framesLabel, framesPhrase, framesText } from '../../src/ui/format';

describe('formatDateTime', () => {
  it('writes a local time in ISO 8601 form, day and month never swapped', () => {
    expect(formatDateTime(new Date(2026, 9, 9, 21, 49, 30).getTime())).toBe('2026-10-09 21:49:30');
    expect(formatDateTime(new Date(2026, 0, 2, 3, 4, 5).getTime())).toBe('2026-01-02 03:04:05');
  });

  it('pads every field and drops the milliseconds', () => {
    expect(formatDateTime(new Date(999, 11, 31, 0, 0, 0, 999).getTime())).toBe('0999-12-31 00:00:00');
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
