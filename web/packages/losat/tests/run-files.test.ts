// The options that name a run's files (stream 8), matched by the command line's file names of the engine's *start*
// (`namedRunFiles`; line RG item F1w-r2, 2026-11-28).
import { describe, expect, it } from 'vitest';
import { namedRunFiles, runFileOption, type RunFile } from '../src/domain/index.ts';

/** The start as runtime.ts gives it to `namedRunFiles`: the whole `StartInfo`, the `-logfile` value as written too. */
type Start = Parameters<typeof namedRunFiles>[1] & { readonly logFile: string | null };

const bytes = new Uint8Array([1]);
const labels = (names: readonly string[], start: Start): string[][] =>
  namedRunFiles(
    names.map((name): RunFile => ({ name, bytes })),
    start,
  ).map((file) => [file.name, runFileOption(file)]);

describe('namedRunFiles', () => {
  it('labels a file by the names of the start, not by the -logfile value as written', () => {
    // `-logfile - -export_search_strategy ./-`: `-logfile -` makes no file, and the strategy file is `-`.
    expect(labels(['-'], { logFile: '-', logFileName: null, strategyFile: '-' })).toEqual([['-', 'export_search_strategy']]);
    // `-logfile " s.asn" -export_search_strategy "./ s.asn"`: the log `s.asn` and the strategy file ` s.asn`.
    expect(labels(['s.asn', ' s.asn'], { logFile: ' s.asn', logFileName: 's.asn', strategyFile: ' s.asn' })).toEqual([
      ['s.asn', 'logfile'],
      [' s.asn', 'export_search_strategy'],
    ]);
    // psiblast `-out_pssm ./- -logfile -`: the PSSM file `-`.
    expect(labels(['-'], { logFile: '-', logFileName: null, strategyFile: null, outPssm: '-' })).toEqual([['-', 'out_pssm']]);
  });

  it('labels the two files of one name after a link apart, and the file of -logfile first', () => {
    // `-logfile /var/run/../tmp/X -export_search_strategy /var/run/../tmp/X`: two files (NCBI 2.17.0 blastp).
    expect(labels(['/var/tmp/X', '/tmp/X'], { logFile: '/var/run/../tmp/X', logFileName: '/var/tmp/X', strategyFile: '/tmp/X' })).toEqual([
      ['/var/tmp/X', 'logfile'],
      ['/tmp/X', 'export_search_strategy'],
    ]);
    // A name of both is the file of -logfile (the engine refuses a strategy file that is the log's).
    expect(labels(['logs'], { logFile: 'logs/', logFileName: 'logs', strategyFile: 'logs', outPssm: 'logs' })).toEqual([['logs', 'logfile']]);
    expect(labels(['s.asn'], { logFile: null, logFileName: null, strategyFile: 's.asn', outPssm: 's.asn' })).toEqual([
      ['s.asn', 'export_search_strategy'],
    ]);
  });

  it('throws on a file that no name of the start gives', () => {
    expect(() => labels(['run.log'], { logFile: null, logFileName: null, strategyFile: null })).toThrow(
      'has a file that the argv does not name: "run.log"',
    );
  });
});
