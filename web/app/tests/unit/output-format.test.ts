import { describe, expect, it } from 'vitest';
import {
  CUSTOM_OUTPUT,
  customOutfmt,
  OUTPUT_FORMATS,
  outfmtOf,
  outfmtOption,
  outputFileName,
  runOutputs,
} from '../../src/domain/output-format';

const argv = (...words: string[]) => ['blastn', '-query', 'q.fa', '-subject', 's.fa', ...words];

// The output of the argv's own -outfmt (WA2-OP-03, Owner 2026-11-05), read as the adapter reads it
// (`run::parse_run`, web/adapter/src/run.rs: the same values in `outfmt_is_read_as_on_the_command_line`).
describe('customOutfmt', () => {
  it.each(['0', '6', '7', ' 7 ', '0 qseqid', '+6'])('adds no output for -outfmt %j', (value) => {
    expect(customOutfmt(argv('-outfmt', value))).toBeUndefined();
    expect(runOutputs(argv('-outfmt', value))).toEqual(OUTPUT_FORMATS);
  });

  it.each(['6 qseqid sseqid pident', '7 std qcovs', '6 std', '6 delim=,', ' 6 -qseqid ', '+7 sseqid', '06 qseqid'])(
    'adds an output for -outfmt %j, without the white space around it',
    (value) => {
      expect(customOutfmt(argv('-outfmt', value))).toBe(value.trim());
      expect(customOutfmt(argv(`-outfmt=${value}`))).toBe(value.trim());
      expect(runOutputs(argv('-outfmt', value))).toEqual([0, 6, 7, CUSTOM_OUTPUT]);
    },
  );

  it('adds none without -outfmt, and for the values that the engine refuses', () => {
    expect(customOutfmt(argv())).toBeUndefined();
    expect(customOutfmt(argv('-outfmt', '10 qseqid'))).toBeUndefined();
    expect(customOutfmt(argv('-outfmt', '6\tqseqid'))).toBeUndefined();
    expect(customOutfmt(argv('-evalue', '1'))).toBeUndefined();
  });

  it('reads -outfmt where the argv gives it, in either form', () => {
    expect(outfmtOption(argv('-evalue', '1', '-outfmt', '7 std'))).toEqual({ at: 7, words: 2, value: '7 std' });
    expect(outfmtOption(argv('-outfmt=6 delim== qseqid'))).toEqual({ at: 5, words: 1, value: '6 delim== qseqid' });
    expect(outfmtOption(argv('-outfmt'))).toBeUndefined();
  });

  it('gives the -outfmt value of each output', () => {
    const custom = argv('-outfmt', '7 std qcovs');
    expect(OUTPUT_FORMATS.map((format) => outfmtOf(custom, format))).toEqual(['0', '6', '7']);
    expect(outfmtOf(custom, CUSTOM_OUTPUT)).toBe('7 std qcovs');
  });
});

describe('outputFileName', () => {
  const run = { number: 3, program: 'blastn', task: 'megablast', queuedAt: 1_762_255_000_000 };

  it("names the output of the argv's -outfmt after its format, as BLAST+'s custom format", () => {
    expect(outputFileName({ ...run, argv: argv() }, 6)).toBe('losat-run3-blastn-megablast-20251104T111640Z.outfmt6.txt');
    expect(outputFileName({ ...run, argv: argv('-outfmt', '6 qseqid') }, CUSTOM_OUTPUT)).toBe(
      'losat-run3-blastn-megablast-20251104T111640Z.outfmt6-custom.txt',
    );
    expect(outputFileName({ ...run, argv: argv('-outfmt', '+07 std qcovs') }, CUSTOM_OUTPUT)).toBe(
      'losat-run3-blastn-megablast-20251104T111640Z.outfmt7-custom.txt',
    );
    expect(outputFileName({ ...run, argv: argv('-outfmt', '7 std') }, 7)).toBe('losat-run3-blastn-megablast-20251104T111640Z.outfmt7.txt');
  });
});

describe('outfmtOption', () => {
  it('reads the names of the inputs as values', () => {
    expect(outfmtOption(['blastn', '-query', '-outfmt=7 x.fa', '-subject', 's.fa'])).toBeUndefined();
    expect(customOutfmt(['blastn', '-query', 'q.fa', '-subject', '-outfmt=6 x.fa'])).toBeUndefined();
  });
});
