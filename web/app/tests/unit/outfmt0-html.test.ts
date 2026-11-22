// The text that the viewer shows for outfmt 0 of a run made with -html (line RG T4a2, decision
// T4a2-Web): NCBI's HTML with the tags removed and the four entities decoded. The expected strings
// are written by hand from the engine's fixtures (LOSAT/tests/fixtures/outfmt0).
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { headingTitle, htmlShownText } from '../../src/domain/outfmt0';

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
    expect(htmlShownText('&amp;lt;&amp;&lt;&quot;&gt;&#38;&auml;&')).toBe('&lt;&<">&#38;&auml;&');
    expect(htmlShownText('no markup\n')).toBe('no markup\n');
    // A "<" that closes nowhere starts a tag to the end of the text.
    expect(htmlShownText('x<b')).toBe('x');
  });
});
