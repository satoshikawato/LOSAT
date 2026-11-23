// Reading the outfmt 0 text that the engine wrote (plan §4.5). Only the subject heading's
// title is cut out for the subject list; the detail shows the heading and the HSP's
// section as written.
import { LEGACY_NAMES, LEGACY_NAME_LENGTH, NAMED_REFERENCES } from './html-named-references';

/**
 * The title of a subject heading (the bytes of an HSP's `out0_subject` range): the text
 * after the leading "> " up to the "Length=" line, as NCBI wrapped it. A list cell shows
 * it on one line (HTML white space), so the wrapping needs no change.
 */
export function headingTitle(heading: string): string {
  const end = heading.lastIndexOf('\nLength=');
  const title = end < 0 ? heading : heading.slice(0, end);
  return title.replace(/^>\s?/, '');
}

/**
 * The text that a browser shows for outfmt 0 bytes written with `-html`, as the HTML tokenizer of the
 * WHATWG HTML Living Standard reads them (section "Tokenization"): the tags, comments and doctypes removed
 * and every character reference decoded, in one pass from left to right (`&amp;lt;` shows `&lt;`), after
 * the input stream's line ends (CR LF and CR) became LF. NCBI encodes a title's `< > & "` but writes a
 * local id raw, before and inside its anchor (`>x<1y<a name=x<1y></a>`), so the text holds what a
 * browser reads as text: a "<" that starts no tag, a reference without ";". The stored bytes and the
 * export stay as written; only what the viewer shows changes. Do not use it for a run without `-html`,
 * whose titles may hold `&lt;` as written. (What the tree builder and the style sheet then do with an
 * element that a raw id writes, such as `<script>` or `</pre>`, is not the tokenizer's: the viewer shows
 * that id's letters as the tokenizer reads them.)
 */
export function htmlShownText(html: string): string {
  const input = html.includes('\r') ? html.replace(/\r\n?/g, '\n') : html;
  const special = /[<&]/g;
  let out = '';
  let at = 0;
  for (;;) {
    special.lastIndex = at;
    const found = special.exec(input);
    if (found === null) return out + input.slice(at);
    out += input.slice(at, found.index);
    const [text, end] = input[found.index] === '<' ? markup(input, found.index + 1) : characterReference(input, found.index + 1);
    out += text;
    at = end;
  }
}

/** The tokenizer's white space in a tag (tab, LF, FF and space; a CR became LF before). */
const isSpace = (char: string | undefined): boolean => char === '\t' || char === '\n' || char === '\f' || char === ' ';
const isAsciiAlpha = (char: string | undefined): boolean => char !== undefined && /^[A-Za-z]$/.test(char);
const isAsciiAlphanumeric = (char: string | undefined): boolean => char !== undefined && /^[A-Za-z0-9]$/.test(char);

/**
 * What a "<" at `at - 1` gives (the tag open state): the text it shows and where the text after it
 * starts. A tag starts at an ASCII letter, an end tag at "/" and one, a comment at "!--", a doctype at
 * "!DOCTYPE" (any case), and any other "!" or "?", or "/" before anything else than a letter or ">",
 * opens a bogus comment to the next ">"; "</>" is nothing. A "<" before anything else is text, and one
 * at the end is "<" ("</" at the end, "</").
 */
function markup(input: string, at: number): [string, number] {
  const char = input[at];
  if (char === undefined) return ['<', at];
  if (isAsciiAlpha(char)) return ['', tagEnd(input, at)];
  if (char === '!') {
    if (input.startsWith('--', at + 1)) return ['', commentEnd(input, at + 3)];
    // A doctype and a bogus comment both end at the next ">" (a doctype's quoted identifier too).
    return ['', toGreaterThan(input, input.slice(at + 1, at + 8).toUpperCase() === 'DOCTYPE' ? at + 8 : at + 1)];
  }
  if (char === '?') return ['', toGreaterThan(input, at)];
  if (char === '/') {
    const next = input[at + 1];
    if (next === undefined) return ['</', at + 1];
    if (isAsciiAlpha(next)) return ['', tagEnd(input, at + 1)];
    if (next === '>') return ['', at + 2];
    return ['', toGreaterThan(input, at + 1)];
  }
  return ['<', at];
}

/** Where a bogus comment or a doctype that starts at `at` ends: after the next ">", or at the end. */
function toGreaterThan(input: string, at: number): number {
  const end = input.indexOf('>', at);
  return end < 0 ? input.length : end + 1;
}

/**
 * Where a comment whose text starts at `at` (after "<!--") ends: "<!-->" and "<!--->" end at once;
 * else after the first "-->" or "--!>" in its text, or at the end.
 */
function commentEnd(input: string, at: number): number {
  if (input[at] === '>') return at + 1;
  if (input.startsWith('->', at)) return at + 2;
  const ends = [input.indexOf('-->', at), input.indexOf('--!>', at)].map((end, kind) => (end < 0 ? Infinity : end + 3 + kind));
  const end = Math.min(...ends);
  return end === Infinity ? input.length : end;
}

/**
 * Where a start or end tag whose name starts at `at` ends: after the ">" that the tokenizer's states of
 * the tag name and its attributes reach outside a quoted value, or at the end (a tag that the text ends
 * in shows nothing). A quote opens a value only where a value starts (after "=", white space allowed);
 * an unquoted value (NCBI writes a local id's anchor so: `<a name=a'b>`) runs to white space or ">" and
 * keeps any quote in it as a letter; a quote or "=" that starts an attribute is a letter of its name.
 */
function tagEnd(input: string, at: number): number {
  let state: 'name' | 'beforeName' | 'attributeName' | 'afterName' | 'beforeValue' | 'double' | 'single' | 'unquoted' | 'afterValue' | 'selfClosing' = 'name';
  for (let i = at; i < input.length; i++) {
    const char = input[i]!;
    switch (state) {
      case 'name':
        if (isSpace(char)) state = 'beforeName';
        else if (char === '/') state = 'selfClosing';
        else if (char === '>') return i + 1;
        break;
      case 'beforeName':
        if (char === '/') state = 'selfClosing';
        else if (char === '>') return i + 1;
        else if (!isSpace(char)) state = 'attributeName';
        break;
      case 'attributeName':
        if (isSpace(char)) state = 'afterName';
        else if (char === '/') state = 'selfClosing';
        else if (char === '=') state = 'beforeValue';
        else if (char === '>') return i + 1;
        break;
      case 'afterName':
        if (char === '/') state = 'selfClosing';
        else if (char === '=') state = 'beforeValue';
        else if (char === '>') return i + 1;
        else if (!isSpace(char)) state = 'attributeName';
        break;
      case 'beforeValue':
        if (char === '"') state = 'double';
        else if (char === "'") state = 'single';
        else if (char === '>') return i + 1;
        else if (!isSpace(char)) state = 'unquoted';
        break;
      case 'double':
        if (char === '"') state = 'afterValue';
        break;
      case 'single':
        if (char === "'") state = 'afterValue';
        break;
      case 'unquoted':
        if (isSpace(char)) state = 'beforeName';
        else if (char === '>') return i + 1;
        break;
      case 'afterValue':
      case 'selfClosing':
        // A letter right after a quoted value or a "/" starts the next attribute's name.
        if (isSpace(char)) state = 'beforeName';
        else if (char === '/') state = 'selfClosing';
        else if (char === '>') return i + 1;
        else state = 'attributeName';
        break;
    }
  }
  return input.length;
}

/**
 * The code points that a numeric character reference of 0x80-0x9F gives (the tokenizer's table of
 * windows-1252's letters); the other numbers of that range stay.
 */
const C1_REPLACEMENTS: ReadonlyMap<number, number> = new Map([
  [0x80, 0x20ac], [0x82, 0x201a], [0x83, 0x0192], [0x84, 0x201e], [0x85, 0x2026], [0x86, 0x2020], [0x87, 0x2021],
  [0x88, 0x02c6], [0x89, 0x2030], [0x8a, 0x0160], [0x8b, 0x2039], [0x8c, 0x0152], [0x8e, 0x017d], [0x91, 0x2018],
  [0x92, 0x2019], [0x93, 0x201c], [0x94, 0x201d], [0x95, 0x2022], [0x96, 0x2013], [0x97, 0x2014], [0x98, 0x02dc],
  [0x99, 0x2122], [0x9a, 0x0161], [0x9b, 0x203a], [0x9c, 0x0153], [0x9e, 0x017e], [0x9f, 0x0178],
]);

/** The digits of a numeric character reference, read where `lastIndex` is set. */
const DECIMAL_DIGITS = /[0-9]+/y;
const HEX_DIGITS = /[0-9A-Fa-f]+/y;

/**
 * What an "&" at `at - 1` in text gives (the character reference state): the text it shows and where the
 * text after it starts.
 * - "&#" and decimal digits, or "&#x" (or "&#X") and hexadecimal ones, with or without ";": the number's
 *   character; 0, a surrogate or a number past 0x10FFFF is U+FFFD, and 0x80-0x9F as C1_REPLACEMENTS.
 *   Without a digit, "&#" or "&#x" is text.
 * - A name: the longest one of the table (`html-named-references.ts`) that the text starts with: a name
 *   with ";" where the letters and digits after "&" are that name and ";" follows, else the longest of the
 *   names without ";" (LEGACY_NAMES) that starts them (`&ltz` shows `<z`, `&notin` shows `¬in`).
 * - Anything else: "&" is text.
 */
function characterReference(input: string, at: number): [string, number] {
  if (input[at] === '#') {
    const hex = input[at + 1] === 'x' || input[at + 1] === 'X';
    const digits = hex ? HEX_DIGITS : DECIMAL_DIGITS;
    const start = at + (hex ? 2 : 1);
    digits.lastIndex = start;
    const run = digits.exec(input)?.[0];
    if (run === undefined) return [input.slice(at - 1, start), start];
    let value = 0;
    for (const digit of run) value = Math.min(value * (hex ? 16 : 10) + parseInt(digit, 16), 0x110000);
    const end = start + run.length + (input[start + run.length] === ';' ? 1 : 0);
    return [numericCharacter(value), end];
  }
  if (!isAsciiAlphanumeric(input[at])) return ['&', at];
  let runEnd = at;
  while (isAsciiAlphanumeric(input[runEnd])) runEnd++;
  const name = input.slice(at, runEnd);
  if (input[runEnd] === ';') {
    const characters = NAMED_REFERENCES.get(name);
    if (characters !== undefined) return [characters, runEnd + 1];
  }
  for (let length = Math.min(name.length, LEGACY_NAME_LENGTH); length > 1; length--) {
    const legacy = name.slice(0, length);
    if (LEGACY_NAMES.has(legacy)) return [NAMED_REFERENCES.get(legacy)!, at + length];
  }
  return ['&', at];
}

/** The character of a numeric character reference of `value` (the numeric character reference end state). */
function numericCharacter(value: number): string {
  if (value === 0 || value > 0x10ffff || (value >= 0xd800 && value <= 0xdfff)) return '\u{fffd}';
  return String.fromCodePoint(C1_REPLACEMENTS.get(value) ?? value);
}
