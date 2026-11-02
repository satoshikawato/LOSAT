// The blocks of lines of the Outputs view. The browser lays out a block (CSS `content-visibility:
// auto`) only when it comes near the view. One text node of a whole output keeps the shaped text
// of every line: about 250 bytes for each character in Chromium 149, 3.9 GiB for the outfmt 6 of
// 15 MB of a TBLASTX search (v4-W results/s29browser).

/** Lines of one block: a few screens of the view (60vh). */
export const OUTPUT_CHUNK_LINES = 128;

export interface OutputChunk {
  readonly text: string;
  readonly lines: number;
}

/** Cuts `text` after every `linesPerChunk` lines; the blocks joined are `text`. */
export function outputChunks(text: string, linesPerChunk: number = OUTPUT_CHUNK_LINES): OutputChunk[] {
  const chunks: OutputChunk[] = [];
  let start = 0;
  let lines = 0;
  for (let newline = text.indexOf('\n'); newline !== -1; newline = text.indexOf('\n', newline + 1)) {
    lines++;
    if (lines === linesPerChunk) {
      chunks.push({ text: text.slice(start, newline + 1), lines });
      start = newline + 1;
      lines = 0;
    }
  }
  if (start < text.length) chunks.push({ text: text.slice(start), lines: text.endsWith('\n') ? lines : lines + 1 });
  return chunks;
}
