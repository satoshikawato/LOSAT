// The options that name a run's files (stream 8): the engine's label of each, checked against the command line's file
// names of the engine's *start*, and those names in order where the archive labels no file (`namedRunFiles`; line RG
// item F1w-r2, 2026-11-28).
import { describe, expect, it } from 'vitest';
import { namedRunFiles, runFileOption, type RunFileOption, type StoredRunFile } from '../src/domain/index.ts';

/** The start as runtime.ts gives it to `namedRunFiles`: the whole `StartInfo`, the `-logfile` value as written too. */
type Start = Parameters<typeof namedRunFiles>[1] & { readonly logFile: string | null };

const bytes = new Uint8Array([1]);
/** The labels of the files `entries`: a name alone is an entry without a label (an earlier module's archive). */
const labels = (entries: ReadonlyArray<string | readonly [string, RunFileOption]>, start: Start): string[][] =>
  namedRunFiles(
    entries.map((entry): StoredRunFile => (typeof entry === 'string' ? { name: entry, bytes } : { name: entry[0], bytes, label: entry[1] })),
    start,
  ).map((file) => [file.name, runFileOption(file)]);

describe('namedRunFiles', () => {
  it("labels each file by the engine's label, which the start's name of that option gives", () => {
    // `-logfile /data/../tmp/X -export_search_strategy /tmp/X` without a folder `/data` (NCBI 2.17.0 blastp makes only
    // the strategy file /tmp/X; the start names the log /tmp/X too).
    expect(
      labels([['/tmp/X', 'export_search_strategy']], { logFile: '/data/../tmp/X', logFileName: '/tmp/X', strategyFile: '/tmp/X' }),
    ).toEqual([['/tmp/X', 'export_search_strategy']]);
    // `-logfile sub/../X -export_search_strategy X` without a folder `sub`.
    expect(labels([['X', 'export_search_strategy']], { logFile: 'sub/../X', logFileName: 'X', strategyFile: 'X' })).toEqual([
      ['X', 'export_search_strategy'],
    ]);
    // psiblast `-logfile /data/../tmp/P -out_pssm /tmp/P` and `-logfile sub/../P -out_ascii_pssm P`: only the PSSM file.
    expect(
      labels([['/tmp/P', 'out_pssm']], { logFile: '/data/../tmp/P', logFileName: '/tmp/P', strategyFile: null, outPssm: '/tmp/P' }),
    ).toEqual([['/tmp/P', 'out_pssm']]);
    expect(
      labels([['P', 'out_ascii_pssm']], { logFile: 'sub/../P', logFileName: 'P', strategyFile: null, outAsciiPssm: 'P' }),
    ).toEqual([['P', 'out_ascii_pssm']]);
    // Both files, and the PSSM of each round (`-save_each_pssm`).
    expect(
      labels(
        [
          ['s.asn', 'export_search_strategy'],
          ['run.log', 'logfile'],
        ],
        { logFile: 'run.log', logFileName: 'run.log', strategyFile: 's.asn' },
      ),
    ).toEqual([
      ['s.asn', 'export_search_strategy'],
      ['run.log', 'logfile'],
    ]);
    expect(
      labels([['out.asn.2', 'out_pssm']], { logFile: null, logFileName: null, strategyFile: null, outPssm: 'out.asn', saveEachPssm: true }),
    ).toEqual([['out.asn.2', 'out_pssm']]);
  });

  it("throws on a labeled file that the start's name of its option does not give", () => {
    expect(() => labels([['run.log', 'export_search_strategy']], { logFile: 'run.log', logFileName: 'run.log', strategyFile: 's.asn' })).toThrow(
      'has a file of -export_search_strategy that the argv names otherwise: "run.log"',
    );
    expect(() => labels([['s.asn', 'logfile']], { logFile: '-', logFileName: null, strategyFile: 's.asn' })).toThrow(
      'has a file of -logfile that the argv names otherwise: "s.asn"',
    );
  });

  it('labels a file without a label by the names of the start, not by the -logfile value as written', () => {
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

  it('labels the two files of one name after a link apart, and a file without a label of both names as the log', () => {
    // `-logfile /var/run/../tmp/X -export_search_strategy /var/run/../tmp/X`: two files (NCBI 2.17.0 blastp).
    expect(
      labels(['/var/tmp/X', '/tmp/X'], { logFile: '/var/run/../tmp/X', logFileName: '/var/tmp/X', strategyFile: '/tmp/X' }),
    ).toEqual([
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
