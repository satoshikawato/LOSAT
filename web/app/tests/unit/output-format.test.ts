import { describe, expect, it } from 'vitest';
import { CUSTOM_OUTPUT, OUTPUT_FORMATS } from 'losat/domain';
import {
  customOutfmt,
  FORMAT_NAMES,
  formatNumber,
  isBinaryOutput,
  outfmtOf,
  outfmtOption,
  outputFileName,
  outputMediaType,
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

  it.each([
    '6 qseqid sseqid pident',
    '7 std qcovs',
    '6 std',
    '6 delim=,',
    ' 6 -qseqid ',
    '+7 sseqid',
    '06 qseqid',
    // Every other format (W11), with or without words.
    '5',
    ' 11 ',
    '+010',
    '10 qseqid',
    '20 delim=| std',
    '17 SQ',
    '5 qseqid',
  ])(
    'adds an output for -outfmt %j, without the white space around it',
    (value) => {
      expect(customOutfmt(argv('-outfmt', value))).toBe(value.trim());
      expect(customOutfmt(argv(`-outfmt=${value}`))).toBe(value.trim());
      expect(runOutputs(argv('-outfmt', value))).toEqual([0, 6, 7, CUSTOM_OUTPUT]);
    },
  );

  it('adds none without -outfmt, and for the values that are not numbers', () => {
    expect(customOutfmt(argv())).toBeUndefined();
    expect(customOutfmt(argv('-outfmt', '6\tqseqid'))).toBeUndefined();
    expect(customOutfmt(argv('-outfmt', 'abc'))).toBeUndefined();
    expect(customOutfmt(argv('-evalue', '1'))).toBeUndefined();
  });

  it('reads the format number as NCBI does', () => {
    expect(['5', ' +05 ', '10 qseqid', '0', '6\tqseqid', 'x', ''].map(formatNumber)).toEqual([5, 5, 10, 0, undefined, undefined, undefined]);
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

  // W11: each format's extension (NCBI's Web BLAST downloads), "-custom" where words follow the number.
  it.each([
    ['5', 'outfmt5.xml', 'application/xml'],
    ['8', 'outfmt8.asn', 'text/plain'],
    ['9', 'outfmt9.asnb', 'application/octet-stream'],
    ['10', 'outfmt10.csv', 'text/csv'],
    ['10 qseqid sseqid', 'outfmt10-custom.csv', 'text/csv'],
    ['11', 'outfmt11.asn', 'text/plain'],
    ['12', 'outfmt12.json', 'application/json'],
    ['13', 'outfmt13.zip', 'application/zip'],
    ['14', 'outfmt14.zip', 'application/zip'],
    ['15', 'outfmt15.json', 'application/json'],
    ['16', 'outfmt16.xml', 'application/xml'],
    ['17', 'outfmt17.sam', 'text/plain'],
    ['17 SQ', 'outfmt17-custom.sam', 'text/plain'],
    ['+020', 'outfmt20.csv', 'text/csv'],
    // Words that the format ignores (NCBI clears them) configure nothing.
    ['5 qseqid', 'outfmt5.xml', 'application/xml'],
  ])('names the output of -outfmt %j …%s (%s)', (value, name, type) => {
    const custom = argv('-outfmt', value);
    expect(outputFileName({ ...run, argv: custom }, CUSTOM_OUTPUT)).toBe(`losat-run3-blastn-megablast-20251104T111640Z.${name}`);
    expect(outputMediaType(custom, CUSTOM_OUTPUT)).toBe(type);
    expect(outputMediaType(custom, 0)).toBe('text/plain');
    expect(isBinaryOutput(custom, CUSTOM_OUTPUT)).toBe(['9', '13', '14'].includes(value));
    expect(isBinaryOutput(custom, 6)).toBe(false);
  });
});

describe('FORMAT_NAMES', () => {
  it("has BLAST+'s name of every format that a run writes", () => {
    for (const number of [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 20]) expect(FORMAT_NAMES[number]).toBeTruthy();
    expect(FORMAT_NAMES[20]).toBe('Comma-separated values with header lines');
  });
});

describe('outfmtOption', () => {
  it('reads the names of the inputs as values', () => {
    expect(outfmtOption(['blastn', '-query', '-outfmt=7 x.fa', '-subject', 's.fa'])).toBeUndefined();
    expect(customOutfmt(['blastn', '-query', 'q.fa', '-subject', '-outfmt=6 x.fa'])).toBeUndefined();
  });
});
