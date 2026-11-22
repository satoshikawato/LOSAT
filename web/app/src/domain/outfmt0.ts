// Reading the outfmt 0 text that the engine wrote (plan §4.5). Only the subject heading's
// title is cut out for the subject list; the detail shows the heading and the HSP's
// section as written.

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

/** The four entities that NCBI's HTML encoding of a title writes (`< > & "`). */
const HTML_ENTITIES: readonly (readonly [string, string])[] = [
  ['&quot;', '"'],
  ['&amp;', '&'],
  ['&lt;', '<'],
  ['&gt;', '>'],
];

/**
 * The text that a browser shows for outfmt 0 bytes written with `-html`: the tags removed (a tag runs
 * from "<" to the next ">" outside a quoted attribute value; a raw "<" in the HTML text always starts
 * one) and the four entities decoded, in one pass from left to right (`&amp;lt;` shows `&lt;`). The
 * stored bytes and the export stay as written; only what the viewer shows changes. Do not use it for
 * a run without `-html`, whose titles may hold `&lt;` as written.
 */
export function htmlShownText(html: string): string {
  const special = /[<&]/g;
  let out = '';
  let at = 0;
  for (;;) {
    special.lastIndex = at;
    const found = special.exec(html);
    if (found === null) return out + html.slice(at);
    out += html.slice(at, found.index);
    at = found.index;
    if (html[at] === '<') {
      // A browser's reading of a tag: a quote opens a quoted value only where a value starts (after "=",
      // spaces allowed); an unquoted value (NCBI writes a local id's anchor so: `<a name=a'b>`) runs to
      // a space or ">" and keeps any quote in it as a letter.
      let quote = '';
      let value: 'none' | 'expected' | 'unquoted' = 'none';
      at += 1;
      while (at < html.length) {
        const char = html[at]!;
        at += 1;
        if (quote !== '') {
          if (char === quote) quote = '';
        } else if (char === '>') {
          break;
        } else if (/\s/.test(char)) {
          if (value === 'unquoted') value = 'none';
        } else if (value === 'unquoted') {
          // part of the value
        } else if (value === 'expected') {
          if (char === '"' || char === "'") {
            quote = char;
            value = 'none';
          } else {
            value = 'unquoted';
          }
        } else if (char === '=') {
          value = 'expected';
        }
      }
    } else {
      const entity = HTML_ENTITIES.find(([name]) => html.startsWith(name, at));
      out += entity === undefined ? '&' : entity[1];
      at += entity === undefined ? 1 : entity[0].length;
    }
  }
}
