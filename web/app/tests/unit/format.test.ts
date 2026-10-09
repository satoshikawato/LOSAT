// The UI's display helpers (src/ui/format.ts).
import { describe, expect, it } from 'vitest';
import { formatDateTime } from '../../src/ui/format';

describe('formatDateTime', () => {
  it('writes a local time in ISO 8601 form, day and month never swapped', () => {
    expect(formatDateTime(new Date(2026, 9, 9, 21, 49, 30).getTime())).toBe('2026-10-09 21:49:30');
    expect(formatDateTime(new Date(2026, 0, 2, 3, 4, 5).getTime())).toBe('2026-01-02 03:04:05');
  });

  it('pads every field and drops the milliseconds', () => {
    expect(formatDateTime(new Date(999, 11, 31, 0, 0, 0, 999).getTime())).toBe('0999-12-31 00:00:00');
  });
});
