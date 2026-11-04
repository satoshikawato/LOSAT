import { describe, expect, it } from 'vitest';
import { inputLead } from '../../src/domain/input-format';

const bytes = (text: string) => new TextEncoder().encode(text);

describe('inputLead', () => {
  it.each([
    ['UTF-8', [0xef, 0xbb, 0xbf, 0x3e, 0x71]],
    ['UTF-16', [0xff, 0xfe, 0x3e, 0x00]],
    ['UTF-16', [0xfe, 0xff, 0x00, 0x3e]],
  ] as const)('names a %s byte order mark', (mark, lead) => {
    expect(inputLead(new Uint8Array(lead), true)).toEqual({ mark });
  });

  it('names the mark of a pasted text that starts with U+FEFF', () => {
    expect(inputLead(bytes('\uFEFF>q\nACGT\n'), true)).toEqual({ mark: 'UTF-8' });
  });

  it.each([
    ['GenBank', 'LOCUS       NC_045512              29903 bp ss-RNA     linear   VRL 18-JUL-2020\nDEFINITION  x\n'],
    ['GenBank', '\n\nLOCUS       AB000001 10 bp DNA\n'],
    ['EMBL', 'ID   X56734; SV 1; linear; mRNA; STD; PLN; 1859 BP.\nXX\n'],
    ['EMBL', 'ID   INS_HUMAN               Reviewed;         110 AA.\n'],
    ['FASTQ', '@r1 extra\nACGTN\n+\nIIII#\n@r2\nAC\n+r2\nII\n'],
    ['FASTQ', '@r1\r\nACGT\r\n+\r\nIIII\r\n'],
  ] as const)('names a %s file', (format, text) => {
    expect(inputLead(bytes(text), true)).toEqual({ format });
  });

  it('names a FASTQ file whose preview ends before the third line, but not a whole input of two lines', () => {
    // A long read whose sequence line runs past the preview, and a preview that ends after the sequence line.
    expect(inputLead(bytes(`@read/1\n${'ACGT'.repeat(600)}`), false)).toEqual({ format: 'FASTQ' });
    expect(inputLead(bytes('@r1\nACGT\n'), false)).toEqual({ format: 'FASTQ' });
    expect(inputLead(bytes('@r1\nACGT\n'), true)).toBeUndefined();
    expect(inputLead(bytes('@r1\nACGT'), true)).toBeUndefined();
    expect(inputLead(bytes('@r1'), false)).toBeUndefined();
  });

  it.each([
    ['FASTA', '>q1 LOCUS\nACGT\n'],
    ['FASTA after blank lines', '\n>q\nACGT\n'],
    ['a bare sequence', 'ACGTACGT\n'],
    ['an empty input', ''],
    ['a word that starts with LOCUS', 'LOCUSX\n'],
    ['ID without the three spaces', 'ID X\n'],
    ['@ followed by a space', '@ r1\nACGT\n+\nIIII\n'],
    ['@ with a third line that is not +', '@r1\nACGT\nACGT\n'],
    ['@ with a second line that is not a sequence', '@r1 x\n12 34\n+\n'],
  ])('names nothing for %s', (_, text) => {
    expect(inputLead(bytes(text), true)).toBeUndefined();
  });
});
