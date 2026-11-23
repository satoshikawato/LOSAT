import { describe, expect, it } from 'vitest';
import type { RunStatus } from '../../src/domain/run';
import { phaseLabel, tabWarning } from '../../src/ui/run-phase';

// Line WQA's D20 round 3, F1 (2026-11-23): a run that prints the command line's help or version searches nothing, so its status
// block and the phase told to a screen reader do not say "Searching", and its tab warning does not say "the search".
describe('run phase words', () => {
  const STATUSES: readonly RunStatus[] = ['queued', 'preparing', 'running', 'finalizing', 'completed', 'cancelled', 'failed'];
  it('say Searching only for a run that searches', () => {
    expect(phaseLabel('running', false)).toBe('Searching');
    expect(phaseLabel('running', true)).toBe('Printing');
    expect(tabWarning(false)).toBe('Closing or reloading this tab ends the search');
    expect(tabWarning(true)).toBe('Closing or reloading this tab ends the run');
  });
  it('are the same for the other phases', () => {
    for (const status of STATUSES.filter((s) => s !== 'running')) expect(phaseLabel(status, true), status).toBe(phaseLabel(status, false));
    expect(STATUSES.map((s) => phaseLabel(s, true))).not.toContain('Searching');
  });
});
