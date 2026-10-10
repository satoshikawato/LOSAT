// The dot plot's SVG file (src/domain/plot-svg.ts, application/plot-export.ts): the same axes,
// ticks, colours, opacity classes and line rules as the canvas, the view as shown, every text
// escaped, nothing that could run or load, and a file written in blocks through the Writer.
import { describe, expect, it } from 'vitest';
import { exportDotPlotSvg } from '../../src/application/plot-export';
import { EQUAL_WEIGHTS, fullView, type View } from '../../src/domain/plot-geometry';
import { labelledTicks, ORIENTATIONS, type MeasureText } from '../../src/domain/plot-layout';
import { axisTicks, identityClass, IDENTITY_CLASSES } from '../../src/domain/plot-scale';
import { dotPlotFileName, dotPlotSvg, SVG_MIME, xmlText, type DotPlotSvgInput } from '../../src/domain/plot-svg';
import { memoryDownloader, type SavedFile } from './support/memory-downloader';

/** 0.55 em a letter, as the browsers' sans-serif fonts come to. */
const measure: MeasureText = (text, font) => text.length * 0.55 * Number.parseFloat(/(\d+(?:\.\d+)?)px/.exec(font)![1]!);

interface Hsp {
  readonly q: [number, number];
  readonly s: [number, number];
  readonly orientation?: 'forward' | 'reverse' | 'unknown';
  readonly pident?: string;
}

function lines(hsps: readonly Hsp[]): DotPlotSvgInput['lines'] {
  const pick = (f: (h: Hsp) => number) => Float64Array.from(hsps, f);
  return {
    count: hsps.length,
    x0: pick((h) => h.q[0]),
    y0: pick((h) => h.s[0]),
    x1: pick((h) => h.q[1]),
    y1: pick((h) => h.s[1]),
    batch: Uint8Array.from(hsps, (h) => ORIENTATIONS.indexOf(h.orientation ?? 'forward') * 4 + identityClass(h.pident ?? '99')),
  };
}

const BOX = { left: 100, top: 50, width: 400, height: 200 };

function plot(hsps: readonly Hsp[], over: Partial<DotPlotSvgInput> = {}): DotPlotSvgInput {
  return {
    run: 3,
    queryId: 'q1',
    subjectId: 's1',
    queryLength: 1000,
    subjectLength: 500,
    units: { query: 'nt', subject: 'nt' },
    width: 520,
    height: 260,
    box: BOX,
    view: fullView({ x: 1000, y: 500 }),
    weights: EQUAL_WEIGHTS,
    toScale: true,
    lines: lines(hsps),
    measure,
    ...over,
  };
}

const svg = (input: DotPlotSvgInput): string => [...dotPlotSvg(input)].join('');
const lineElements = (text: string): string[] => text.match(/<line [^>]*\/>/g) ?? [];

interface Tag {
  readonly name: string;
  readonly attrs: Readonly<Record<string, string>>;
}

/**
 * Reads the file as the XML it must be: only tags with quoted attributes and text in which no `<`
 * and no bare `&` appear. Returns the tags in order; throws where the file is not so, or where the
 * tags do not nest.
 */
function tags(text: string): Tag[] {
  const body = text.replace(/^<\?xml[^>]*\?>\n/, '');
  const out: Tag[] = [];
  const open: string[] = [];
  let last = 0;
  for (const m of body.matchAll(/<(\/?)([A-Za-z][\w:-]*)((?:\s+[\w:-]+="[^"<]*")*)\s*(\/?)>/g)) {
    const between = body.slice(last, m.index);
    if (/[<>]|&(?!(?:amp|lt|gt|quot|apos);)/.test(between)) throw new Error(`unescaped text: ${between.slice(0, 80)}`);
    last = m.index + m[0].length;
    const [, closing, name, attributes, selfClosing] = m;
    if (closing === '/') {
      if (open.pop() !== name) throw new Error(`tags do not nest at </${name}>`);
      continue;
    }
    const attrs: Record<string, string> = {};
    for (const a of attributes!.matchAll(/([\w:-]+)="([^"<]*)"/g)) attrs[a[1]!] = a[2]!;
    out.push({ name: name!, attrs });
    if (selfClosing !== '/') open.push(name!);
  }
  if (/[<>]/.test(body.slice(last)) || open.length > 0) throw new Error('the file is cut or not closed');
  return out;
}

describe('dotPlotSvg: the file', () => {
  it('is an SVG of the plot\'s size in CSS pixels, with a title and a description that name the run and the format', () => {
    const text = svg(plot([{ q: [100, 300], s: [50, 150] }]));
    expect(text.startsWith('<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="520" height="260" viewBox="0 0 520 260"')).toBe(true);
    expect(text).toContain('<title>LOSAT Web dot plot of run 3: query q1 against subject s1</title>');
    expect(text).toMatch(/<desc>[^<]*run 3[^<]*query q1 \(1000 nt\)[^<]*subject s1 \(500 nt\)[^<]*LOSAT Web format[^<]*not an NCBI graphic\.<\/desc>/);
    expect(tags(text)[0]!.name).toBe('svg');
    expect(text.trimEnd().endsWith('</svg>')).toBe(true);
  });

  it('draws an HSP from its start to its end in the frame\'s pixels, as the canvas does', () => {
    // Full view of 1000 × 500 in 400 × 200 px: 0.4 px a letter on both axes.
    const text = svg(plot([{ q: [100, 300], s: [50, 150] }]));
    expect(lineElements(text)).toEqual(['<line x1="140" y1="70" x2="220" y2="110"/>']);
  });

  it('groups lines by colour and opacity class: the script\'s colours, the identity classes', () => {
    const hsps: Hsp[] = [
      { q: [0, 100], s: [0, 100], pident: '95.2' },
      { q: [200, 300], s: [100, 0], orientation: 'reverse', pident: '75' },
      { q: [400, 500], s: [200, 300], orientation: 'unknown', pident: '61.5' },
      { q: [600, 700], s: [300, 400], pident: '30' },
    ];
    const groups = tags(svg(plot(hsps))).filter((t) => t.attrs['stroke-opacity'] !== undefined);
    expect(groups.map((g) => [g.attrs['stroke'], g.attrs['stroke-opacity']])).toEqual([
      ['#1f77b4', String(IDENTITY_CLASSES[0]!.opacity)],
      ['#1f77b4', '1'],
      ['#ff7f0e', String(IDENTITY_CLASSES[2]!.opacity)],
      ['#7f7f7f', String(IDENTITY_CLASSES[1]!.opacity)],
    ]);
  });

  it('writes only the lines the view can show, and none of the groups that stay empty', () => {
    const hsps: Hsp[] = [
      { q: [100, 200], s: [100, 200] },
      { q: [800, 900], s: [400, 450], orientation: 'reverse' },
    ];
    const zoomed: View = { x0: 0, x1: 500, y0: 0, y1: 250 };
    const text = svg(plot(hsps, { view: zoomed }));
    expect(lineElements(text)).toHaveLength(1);
    expect(text).not.toContain('#ff7f0e');
  });

  it('leaves out opaque lines that fall on the pixels of one before them, and no translucent one', () => {
    const same = (pident: string): Hsp[] => [
      { q: [100, 300], s: [50, 150], pident },
      { q: [100.2, 300.1], s: [50.1, 150], pident },
    ];
    expect(lineElements(svg(plot(same('99'))))).toHaveLength(1);
    expect(lineElements(svg(plot(same('50'))))).toHaveLength(2);
  });

  it('draws an HSP shorter than a pixel or two as a 3 px dash', () => {
    const text = svg(plot([{ q: [500, 500], s: [250, 250] }]));
    // (500, 250) is at (300, 150).
    expect(lineElements(text)).toEqual(['<line x1="298.5" y1="150" x2="301.5" y2="150"/>']);
  });

  it('clips the lines to the frame by a nested svg, with no clip path', () => {
    const nested = tags(svg(plot([{ q: [0, 1000], s: [0, 500] }]))).filter((t) => t.name === 'svg')[1]!;
    expect(nested.attrs).toMatchObject({ x: '100', y: '50', width: '400', height: '200', viewBox: '100 50 400 200', overflow: 'hidden' });
  });

  it('writes the grid, the frame and the ticks of the canvas, and the labels it would label', () => {
    const input = plot([]);
    const text = svg(input);
    const xTicks = axisTicks(0, 1000, 'nt', 0.4);
    const yTicks = axisTicks(0, 500, 'nt', 0.4);
    const grid = tags(text).find((t) => t.attrs['stroke'] === '#d3d3d3')!;
    expect(grid.attrs['d']!.match(/V/g)).toHaveLength(xTicks.major.length + xTicks.minor.length);
    expect(grid.attrs['d']!.match(/H/g)).toHaveLength(yTicks.major.length + yTicks.minor.length);
    const frame = tags(text).find((t) => t.name === 'rect' && t.attrs['stroke'] === '#000000')!;
    expect(frame.attrs).toMatchObject({ x: '99.5', y: '49.5', width: '401', height: '201' });
    const labels = [...text.matchAll(/<text [^>]*>([^<]*)<\/text>/g)].map((m) => m[1]!);
    const expected = [...labelledTicks(xTicks, 0.4, measure), ...labelledTicks(yTicks, 0.4, measure)].map((t) => t.toLocaleString('en-US'));
    expect(labels.slice(0, expected.length)).toEqual(expected);
    expect(labels.slice(expected.length)).toEqual(['Query q1 (bp)', 'Subject s1 (bp)']);
  });

  it('titles the axes with the IDs and the unit of the view\'s span (kbp from 5,000 letters)', () => {
    const input = plot([], { queryLength: 20000, subjectLength: 20000, view: fullView({ x: 20000, y: 20000 }), units: { query: 'nt', subject: 'aa' } });
    const text = svg(input);
    expect(text).toContain('>Query q1 (kbp)</text>');
    expect(text).toContain('>Subject s1 (kaa)</text>');
  });

  it('says in its description when the axes are not to scale, and the nt per aa of a protein axis', () => {
    expect(svg(plot([], { toScale: false }))).toContain('The axes are not to scale');
    expect(svg(plot([]))).not.toContain('not to scale');
    expect(svg(plot([], { weights: { x: 3, y: 1 } }))).toContain('>Query q1 (bp; drawn at 3 nt per aa)</text>');
  });

  it('cuts a long ID of the title so that the unit stays, as the canvas does', () => {
    const text = svg(plot([], { subjectId: 'S'.repeat(300) }));
    expect(text).toMatch(/>Subject S+… \(bp\)<\/text>/);
  });

  it('leaves out the selected HSP\'s halo and has no legend or other drawing but the plot', () => {
    const text = svg(plot([{ q: [100, 300], s: [50, 150] }]));
    expect(text).not.toContain('rgba(');
    expect(text).not.toContain('#ffc400');
  });
});

describe('dotPlotSvg: text from the inputs', () => {
  const hostile = `a"><script>alert(1)</script><x onload='alert(2)' href="http://evil.example/?a=1&b=2">`;

  it('escapes & < > " \' in IDs, in the title, the description and the axis titles', () => {
    const text = svg(plot([{ q: [100, 300], s: [50, 150] }], { queryId: hostile, subjectId: hostile }));
    const all = tags(text); // throws where any text is unescaped
    expect(new Set(all.map((t) => t.name))).toEqual(new Set(['svg', 'title', 'desc', 'rect', 'path', 'g', 'line', 'text']));
    expect(all.some((t) => t.name === 'script' || t.name === 'x')).toBe(false);
    expect(all.flatMap((t) => Object.keys(t.attrs)).some((a) => a.startsWith('on') || a.includes('href'))).toBe(false);
    expect(text).toContain('a&quot;&gt;&lt;script&gt;alert(1)&lt;/script&gt;&lt;x onload=&apos;alert(2)&apos; href=&quot;http://evil.example/?a=1&amp;b=2&quot;&gt;');
  });

  it('has no script, no event attribute and no reference for ordinary IDs either', () => {
    const text = svg(plot([{ q: [100, 300], s: [50, 150] }]));
    const all = tags(text);
    expect(all.map((t) => t.name)).not.toContain('script');
    expect(all.map((t) => t.name)).not.toContain('use');
    expect(all.map((t) => t.name)).not.toContain('image');
    expect(all.flatMap((t) => Object.keys(t.attrs)).filter((a) => a !== 'xmlns' && (a.startsWith('on') || a.includes('href')))).toEqual([]);
    expect(text).not.toMatch(/href|url\(|<script|<style|<foreignObject/i);
  });

  it('replaces characters that XML cannot hold, so that the file still parses', () => {
    expect(xmlText('a\u0001b\u0000c￾d\uD800e\u{1F9EC}f\tg')).toBe('a�b�c�d�e\u{1F9EC}f\tg');
    const text = svg(plot([], { queryId: 'x\u0007y' }));
    expect(() => tags(text)).not.toThrow();
    expect(text).toContain('x�y');
  });
});

describe('exportDotPlotSvg: the file', () => {
  it('is named by the run and the records\' positions, never by an ID, and saved as image/svg+xml', async () => {
    const saved: SavedFile[] = [];
    const written = await exportDotPlotSvg(memoryDownloader((file) => saved.push(file)), plot([{ q: [100, 300], s: [50, 150] }], { queryId: 'secret/name' }), 2, 7);
    expect(saved).toHaveLength(1);
    expect(saved[0]!.name).toBe('losat-run3-dotplot-q2-s7.svg');
    expect(saved[0]!.mime).toBe(SVG_MIME);
    expect(saved[0]!.mime).toBe('image/svg+xml');
    expect(saved[0]!.name).not.toContain('secret');
    expect(written).toBe(saved[0]!.bytes.length);
    const text = new TextDecoder().decode(saved[0]!.bytes);
    expect(text).toBe(svg(plot([{ q: [100, 300], s: [50, 150] }], { queryId: 'secret/name' })));
  });

  it('writes a pair of 5,993 HSPs one line each, in order, without building the file at once', async () => {
    const hsps: Hsp[] = Array.from({ length: 5993 }, (_, i) => ({ q: [i * 150, i * 150 + 120], s: [i * 70, i * 70 + 60], orientation: i % 2 === 0 ? 'forward' : 'reverse', pident: String(40 + (i % 30)) }));
    const input = plot(hsps, { view: fullView({ x: 1_000_000, y: 500_000 }) });
    const blocks = [...dotPlotSvg(input)];
    expect(blocks.length).toBeGreaterThan(5993);
    expect(Math.max(...blocks.map((b) => b.length))).toBeLessThan(1500);
    const saved: SavedFile[] = [];
    await exportDotPlotSvg(memoryDownloader((file) => saved.push(file)), input, 1, 1);
    const text = new TextDecoder().decode(saved[0]!.bytes);
    expect(() => tags(text)).not.toThrow();
    // All of these are translucent (an identity of 40 to 69), so none is left out as a repeat.
    expect(lineElements(text)).toHaveLength(5993);
  });

  it('saves nothing when the building fails', async () => {
    const saved: SavedFile[] = [];
    const broken: MeasureText = () => {
      throw new Error('no font');
    };
    await expect(exportDotPlotSvg(memoryDownloader((file) => saved.push(file)), plot([], { measure: broken }), 1, 1)).rejects.toThrow('no font');
    expect(saved).toEqual([]);
  });

  it('names the file by positions', () => {
    expect(dotPlotFileName(12, 1, 30)).toBe('losat-run12-dotplot-q1-s30.svg');
  });
});
