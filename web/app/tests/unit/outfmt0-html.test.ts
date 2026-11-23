// The text that the viewer shows for outfmt 0 of a run made with -html (line RG T4a2, decisions
// T4a2-5 and T4a2-17): NCBI's HTML as a browser's HTML tokenizer reads it, the tags, comments and
// doctypes removed and the character references decoded. The expected strings are written by hand
// from the engine's fixtures (LOSAT/tests/fixtures/outfmt0), NCBI BLAST+ 2.17.0's bytes of local ids
// (line RG T4a2 review e S3) and the tokenizer's rules (WHATWG HTML, section "Tokenization").
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { LEGACY_NAMES, NAMED_REFERENCES } from '../../src/domain/html-named-references';
import { headingTitle } from '../../src/domain/outfmt0';
import { htmlShownText } from '../../src/domain/outfmt0-html';

const fixture = (name: string) => readFileSync(resolve(import.meta.dirname, '../../../../LOSAT/tests/fixtures/outfmt0', name), 'utf8');

/** The lines of a fixture from the first that starts with `first` up to and including the first one that starts with `last` after it. */
function lines(text: string, first: string, last: string): string {
  const all = text.split('\n');
  const from = all.findIndex((line) => line.startsWith(first));
  expect(from, first).toBeGreaterThanOrEqual(0);
  const to = all.findIndex((line, at) => at >= from && line.startsWith(last));
  expect(to, last).toBeGreaterThanOrEqual(from);
  return all.slice(from, to + 1).join('\n') + '\n';
}

describe('htmlShownText', () => {
  const protein = fixture('html.bp.sids.f0.pd.out');

  it('shows a subject heading with a report link and an anchor as its plain text', () => {
    const heading = lines(protein, '><a title="Show report for XP_000101.1"', 'Length=');
    expect(heading).toContain('<a name=XP_000101></a>');
    expect(htmlShownText(heading)).toBe('>XP_000101.1 subject one two HSPs\nLength=450\n');
    expect(headingTitle(htmlShownText(heading))).toBe('XP_000101.1 subject one two HSPs');
  });

  it('shows a local identifier, which has an anchor and no link, as its plain text', () => {
    const heading = lines(protein, '>s3_local_identifier', 'Length=');
    expect(heading.startsWith('>s3_local_identifier<a name=s3_local_identifier></a> subject three')).toBe(true);
    expect(htmlShownText(heading)).toBe(
      '>s3_local_identifier subject three with a title that is longer than the description \n' +
        'column of the organism report table so that NCBI cuts it and \nappends three dots\nLength=150\n',
    );
    expect(headingTitle(htmlShownText(heading))).toBe(
      's3_local_identifier subject three with a title that is longer than the description \n' +
        'column of the organism report table so that NCBI cuts it and \nappends three dots',
    );
  });

  it('decodes the entities of a title that holds the characters of a template: &lt;@score@&gt;', () => {
    const heading = lines(protein, '><a title="Show report for XP_000105.1"', 'Length=');
    expect(heading).toContain('&lt;@score@&gt; s &lt;@evalue@&gt; e &lt;@acc_hd@&gt; \na &lt;@descr_text@&gt; d &lt;@table_rows@&gt; end');
    expect(headingTitle(htmlShownText(heading))).toBe('XP_000105.1 subject five <@score@> s <@evalue@> e <@acc_hd@> \na <@descr_text@> d <@table_rows@> end');
  });

  it('shows a masked nucleotide row without its font tags', () => {
    const rows = lines(fixture('html.bn.mk.lc.f0.out'), 'Query  61   <font', 'Query  61   <font');
    expect(rows).toBe(
      'Query  61   <font color="#000000">atgagacaaagtcga</font>GACAAGACCATTGTGGTAAAAGTTC<font color="#000000">gcagttgggtattgtccgag</font>  120\n',
    );
    expect(htmlShownText(rows)).toBe('Query  61   atgagacaaagtcgaGACAAGACCATTGTGGTAAAAGTTCgcagttgggtattgtccgag  120\n');
  });

  it('shows a tblastn heading with an anchor only, and a masked protein row with the frame line before it', () => {
    const section = lines(fixture('html.tn.mk.segtouch.f0.out'), '><a name=Subject_1></a>', 'Query  1    SYVP');
    expect(htmlShownText(section)).toBe(
      '> ns back translation\nLength=468\n\n Score = 171 bits (434),  Expect = 3e-61, Method: Compositional matrix adjust.\n' +
        ' Identities = 116/116 (100%), Positives = 116/116 (100%), Gaps = 0/116 (0%)\n Frame = +1\n\n' +
        'Query  1    SYVPLHFRERFFLHSNMGYWLDLCENYPNMPNNAWYINIRTCRNCqqqqqqqqqqqqqqq  60\n',
    );
  });

  it('removes a tag up to the ">" outside its quoted attribute values, and decodes in one pass from left to right', () => {
    expect(htmlShownText('a<a title="x > y" href=\'z>\'>b</a>c')).toBe('abc');
    expect(htmlShownText('&amp;lt;&amp;&lt;&quot;&gt;&#38;&auml;&')).toBe('&lt;&<">&\u00e4&');
    expect(htmlShownText('no markup\n')).toBe('no markup\n');
    // A "<" that closes nowhere starts a tag to the end of the text.
    expect(htmlShownText('x<b')).toBe('x');
  });

  it('opens a quoted value only where a value starts: an unquoted value keeps a quote as a letter', () => {
    // NCBI writes the anchor of a local id unquoted (the bytes of `blastp -subject -parse_deflines -html` for the ids a'b, c"d, e>f).
    expect(htmlShownText(">a'b<a name=a'b></a> apostrophe id title\nLength=450\n")).toBe(">a'b apostrophe id title\nLength=450\n");
    expect(htmlShownText('>c"d<a name=c"d></a> dquote id title\nLength=450\n')).toBe('>c"d dquote id title\nLength=450\n');
    // An unquoted value ends at ">": the "f>" after it is text, as a browser shows it.
    expect(htmlShownText('>e>f<a name=e>f></a> gt id title\nLength=450\n')).toBe('>e>ff> gt id title\nLength=450\n');
    // A quoted value after "=" (spaces allowed) holds ">" and the other quote, and ends at its own quote.
    expect(htmlShownText('<a title = "x > \'y\'" name=z>t</a>')).toBe('t');
    expect(htmlShownText("<a title='say \"> hi\"' href=u>t</a>")).toBe('t');
    // After an unquoted value and a space the next attribute may be quoted.
    expect(htmlShownText('<a name=a\'b title="p>q">t</a>')).toBe('t');
    // Only tab, LF, FF, CR and space end an unquoted value (review e N2): with U+00A0 the value goes on,
    // "=" and the quote are its letters, and the first ">" ends the tag.
    expect(htmlShownText('<a name=x\u00a0b="p>q">t</a>')).toBe('q">t');
    expect(htmlShownText('<a name=x\fb="p>q">t</a>')).toBe('t');
    expect(htmlShownText('<a name=x\rb="p>q">t</a>')).toBe('t');
    // A quote or "=" that starts an attribute is a letter of its name, not the start of a value.
    expect(htmlShownText('<a ="x>y" b>z')).toBe('y" b>z');
    expect(htmlShownText('<a b"c="d>e">f')).toBe('f');
    expect(htmlShownText('<a b="c"d="e>f">g')).toBe('g');
    expect(htmlShownText('<br/>a<a / b="c>d">e')).toBe('ae');
  });

  it('shows the local ids that NCBI writes raw as a browser shows them (review e S3)', () => {
    // NCBI BLAST+ 2.17.0, `blastp -subject -parse_deflines -html -outfmt 0` (T/t4/review_t4a2_e/out/web.ids.pd.ncbi.out).
    expect(htmlShownText('>x<1y<a name=x<1y></a> title 0\nLength=150\n')).toBe('>x<1y title 0\nLength=150\n');
    expect(htmlShownText('>m&#62;n<a name=m&#62;n></a> title 2\nLength=150\n')).toBe('>m>n title 2\nLength=150\n');
    expect(htmlShownText('>k&ltz<a name=k&ltz></a> title 3\nLength=150\n')).toBe('>k<z title 3\nLength=150\n');
    expect(headingTitle(htmlShownText('>k&ltz<a name=k&ltz></a> title 3\nLength=150\n'))).toBe('k<z title 3');
  });

  it('starts a tag only at "<" before an ASCII letter, "/", "!" or "?"', () => {
    for (const text of ['a<1b', 'a< b', 'a<=b', 'a<<1', 'a<\u00e9b', 'a<']) expect(htmlShownText(text), text).toBe(text);
    expect(htmlShownText('a</')).toBe('a</');
    expect(htmlShownText('a</>b')).toBe('ab');
    expect(htmlShownText('a</ x>b')).toBe('ab');
    expect(htmlShownText('a</1>b')).toBe('ab');
    expect(htmlShownText('a</b c="d>e">f')).toBe('af');
    expect(htmlShownText('a<?x y="z>w">b')).toBe('aw">b');
    expect(htmlShownText('a<!x>b<!>c')).toBe('abc');
    expect(htmlShownText('a<!doctype html>b<!DOCTYPE x "y>z">c')).toBe('abz">c');
    expect(htmlShownText('a<!-x>b')).toBe('ab');
  });

  it('ends a comment as the tokenizer does', () => {
    expect(htmlShownText('a<!-- b > c -->d')).toBe('ad');
    expect(htmlShownText('a<!-->b')).toBe('ab');
    expect(htmlShownText('a<!--->b')).toBe('ab');
    expect(htmlShownText('a<!---->b')).toBe('ab');
    expect(htmlShownText('a<!-- x --!> y')).toBe('a y');
    expect(htmlShownText('a<!--!>b-->c')).toBe('ac');
    expect(htmlShownText('a<!-- -- --->b')).toBe('ab');
    expect(htmlShownText('a<!--b')).toBe('a');
  });

  it('decodes every character reference as a browser does', () => {
    // Numeric: decimal and hexadecimal, with and without ";".
    expect(htmlShownText('&#62;&#62&#x3e;&#X3E&#x3Eg')).toBe('>>>>>g');
    // No digit: the text stays.
    expect(htmlShownText('&#;&#x;&#xg&#a')).toBe('&#;&#x;&#xg&#a');
    // 0, a surrogate and a number past U+10FFFF give U+FFFD; 0x80-0x9F windows-1252's letters, the
    // numbers of that range that it does not use stay; other controls and noncharacters stay.
    expect(htmlShownText('&#0;&#xD800;&#x110000;&#99999999999999999999;')).toBe('\ufffd'.repeat(4));
    expect(htmlShownText('&#128;&#x9F;&#x81;&#1;&#xFFFF;&#13;')).toBe('\u20ac\u0178\u0081\u0001\uffff\r');
    expect(htmlShownText('&#x1F600;')).toBe('\u{1f600}');
    // Named, with ";" (every name of the table) and without it (the legacy names only, the longest
    // that starts the letters).
    expect(htmlShownText('&hearts;&hearts&CounterClockwiseContourIntegral;')).toBe('\u2665&hearts\u2233');
    expect(htmlShownText('&amp&ampx&AMP;&LT&notin;&notin&notit;&ltz;')).toBe('&&x&<\u2209\u00acin\u00acit;<z;');
    expect(htmlShownText('&frac12&frac123&iexcl!')).toBe('\u00bd\u00bd3\u00a1!');
    expect(htmlShownText('&nbsp;&NewLine;&Tab;')).toBe('\u00a0\n\t');
    // An "&" before anything else, or before a name that the table does not have, is text.
    expect(htmlShownText('& &; &x; &unknown; &')).toBe('& &; &x; &unknown; &');
  });

  it('reads CR LF and CR as LF, as the input stream does', () => {
    expect(htmlShownText('a\r\nb\rc\n')).toBe('a\nb\nc\n');
    expect(htmlShownText('<a\rname=x\r>t</a>')).toBe('t');
  });
});

describe('the named character references', () => {
  it('are the table of the WHATWG HTML standard (2,231 identifiers, 106 of them without ";")', () => {
    expect(NAMED_REFERENCES.size).toBe(2125);
    expect(LEGACY_NAMES.size).toBe(106);
    for (const name of LEGACY_NAMES) expect(NAMED_REFERENCES.has(name), name).toBe(true);
    expect([NAMED_REFERENCES.get('amp'), NAMED_REFERENCES.get('nGt'), NAMED_REFERENCES.get('zwnj')]).toEqual(['&', '\u226b\u20d2', '\u200c']);
  });
});
