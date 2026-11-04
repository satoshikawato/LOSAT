import { describe, expect, it } from 'vitest';
import { SourceRefs } from '../../src/application/source-refs';

function refs() {
  const dropped: string[] = [];
  return { refs: new SourceRefs((sourceId) => dropped.push(sourceId)), dropped };
}

describe('SourceRefs', () => {
  it('drops a source when its last holder lets it go, not before', () => {
    const { refs: r, dropped } = refs();
    r.set('draft', ['a', 'b']);
    r.set('run:1', ['a']);
    r.set('draft', ['b']);
    expect(dropped).toEqual([]);
    expect(r.holds('a')).toBe(true);
    r.release('run:1');
    expect(dropped).toEqual(['a']);
    expect(r.holds('a')).toBe(false);
    r.release('draft');
    expect(dropped).toEqual(['a', 'b']);
  });

  it('keeps what a holder names again and drops each source once', () => {
    const { refs: r, dropped } = refs();
    r.set('draft', ['a']);
    r.set('draft', ['a', 'a']);
    r.set('index:query:paste', ['a']);
    r.release('index:query:paste');
    expect(dropped).toEqual([]);
    r.release('draft');
    r.release('draft');
    r.release('nobody');
    expect(dropped).toEqual(['a']);
  });

  it('holds a source again after it was dropped (a new hold of the same id)', () => {
    const { refs: r, dropped } = refs();
    r.set('draft', ['a']);
    r.release('draft');
    r.set('run:2', ['a']);
    expect(r.holds('a')).toBe(true);
    r.release('run:2');
    expect(dropped).toEqual(['a', 'a']);
  });
});
