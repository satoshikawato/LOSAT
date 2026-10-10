// The record that an engine message points at by its line (plan §5.4). Since session SF the
// engine reads every input with its port of NCBI BLAST+'s FASTA reader, whose errors name a
// line of the input ("BLAST query error: CFastaReader: Near line 7, there's a line that
// doesn't look like plausible data, ...") as do LOSAT Web's rejections ("line 3 is a gap
// line ('>?'), ... not supported by LOSAT Web"; "the first line (...) is not a defline ...").
// The Data worker has the bytes that the engine read and where each record lies in them, so
// it finds the record; the engine's message stays as it is.
//
// Line numbers are those of NCBI's line reader (CStreamLineReader), which the engine ports
// and the adapter's index scan follows (web/adapter/src/scan/ncbi.rs `Lines`; abi_v2.md §9:
// every line counts, CR LF is one line end). CR, LF and CR LF each end a line, a lone CR
// included, except where the reader reads past an end of line and joins two lines of the
// file: after a lone CR in a file of LF (or CR LF) line ends, the LF that ends that line of
// the file is not a line end; after an LF in a file of CR line ends, the CR that ends that
// line of the file is not one either.

const CR = 0x0d;
const LF = 0x0a;

/** "Near line 7," and "line 3 is a gap line": NCBI's line number in a message. */
const LINE_IN_MESSAGE = /\bline (\d+)\b/;
/** LOSAT Web's rejection of a first line that NCBI may read as a sequence identifier. */
const FIRST_LINE_IN_MESSAGE = /\bthe first line\b/;

/** The 1-based line that a message names, or undefined. */
export function lineInMessage(message: string): number | undefined {
  const match = LINE_IN_MESSAGE.exec(message);
  if (match !== null) return Number(match[1]);
  return FIRST_LINE_IN_MESSAGE.test(message) ? 1 : undefined;
}

/** The end-of-line style of NCBI's line reader: unknown until the first line end. */
type Eol = 'unknown' | 'cr' | 'lf' | 'crlf' | 'mixed';

/**
 * The offset of the first byte of line `line` (1-based) as NCBI's line reader numbers the
 * lines of `bytes`, or undefined when the input has fewer lines.
 */
export function ncbiLineStart(bytes: Uint8Array, line: number): number | undefined {
  if (!Number.isSafeInteger(line) || line < 1) return undefined;
  let eol: Eol = 'unknown';
  // A CR ended the line in the CR styles: an LF next is part of its end of line.
  let crEnd = false;
  // A CR in the LF styles: an LF next makes it a CR LF; anything else ends the line there.
  let lfCr = false;
  // The end of line of the file that the reader consumed while it read past a line end.
  let dropLf = false;
  let dropCr = false;
  let inLine = false;
  let number = 0;
  for (let offset = 0; offset < bytes.length; offset++) {
    const byte = bytes[offset]!;
    if (byte === LF && dropLf) {
      dropLf = false;
      continue;
    }
    if (byte === CR && dropCr) {
      dropCr = false;
      continue;
    }
    if (crEnd) {
      crEnd = false;
      if (eol === 'unknown') eol = byte === LF ? 'crlf' : 'cr';
      inLine = false;
      if (byte === LF) continue;
    }
    if (lfCr) {
      lfCr = false;
      inLine = false;
      if (byte === LF) continue;
      // A lone CR: the rest of the file's line is read again, joined to the next one.
      dropLf = true;
      eol = eol === 'crlf' ? 'cr' : 'mixed';
    }
    if (!inLine) {
      inLine = true;
      number++;
      if (number === line) return offset;
    }
    if (byte === CR) {
      if (eol === 'lf' || eol === 'crlf') lfCr = true;
      else crEnd = true;
    } else if (byte === LF) {
      inLine = false;
      if (eol === 'unknown') eol = 'crlf';
      else if (eol === 'crlf') eol = 'lf';
      else if (eol === 'cr') {
        // An LF in a file of CR line ends: the rest of the file's line is read again.
        eol = 'mixed';
        dropCr = true;
      }
    }
  }
  return undefined;
}

/**
 * The position, among `spans` (the [start, end) byte ranges of the records of `bytes`, in
 * order), of the record that holds the line that `message` names; undefined when the
 * message names no line, or the line lies outside every record (before the first record,
 * or between the inputs of a combined run input).
 */
export function recordOfMessageLine(
  message: string,
  bytes: Uint8Array,
  spans: ReadonlyArray<readonly [number, number]>,
): number | undefined {
  const line = lineInMessage(message);
  if (line === undefined) return undefined;
  const offset = ncbiLineStart(bytes, line);
  if (offset === undefined) return undefined;
  let low = 0;
  let high = spans.length;
  while (low < high) {
    const middle = (low + high) >>> 1;
    if (spans[middle]![1] <= offset) low = middle + 1;
    else high = middle;
  }
  const span = spans[low];
  return span !== undefined && span[0] <= offset ? low : undefined;
}
