// The blocks of lines of the Outputs view (src/domain/output-chunks.ts).
import { describe, expect, it } from 'vitest';
import { outputChunks } from '../../src/domain/output-chunks';

describe('outputChunks', () => {
  it('joins back to the output, byte for byte', () => {
    for (const text of ['', 'a', 'a\n', '\n', '\n\n\n', 'a\nb', 'a\nb\n', 'a\n\nb\nc\nd\ne\n\n', 'a\r\nb\r\n']) {
      expect(outputChunks(text, 2).map((chunk) => chunk.text).join('')).toBe(text);
    }
  });

  it('cuts after every full block of lines', () => {
    expect(outputChunks('a\nb\nc\nd\ne', 2)).toEqual([
      { text: 'a\nb\n', lines: 2 },
      { text: 'c\nd\n', lines: 2 },
      { text: 'e', lines: 1 },
    ]);
  });

  it('has no last block when the output ends with a full block', () => {
    expect(outputChunks('a\nb\n', 2)).toEqual([{ text: 'a\nb\n', lines: 2 }]);
  });

  it('counts the lines of a last block that ends with a line break', () => {
    expect(outputChunks('a\nb\nc\n', 2)).toEqual([
      { text: 'a\nb\n', lines: 2 },
      { text: 'c\n', lines: 1 },
    ]);
  });

  it('has no block for an empty output', () => {
    expect(outputChunks('')).toEqual([]);
  });
});
