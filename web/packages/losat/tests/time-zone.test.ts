// The engine's time zone (src/infra/reactor/time-zone.ts; line NP WW4b1, D-WW4b-2, 2026-11-27): the POSIX rule that the
// host puts into WASI's `TZ`, made from the browser's zone for the current year. The same rules are read by jiff and
// compared with the zone files hour by hour in LOSAT (writedb/local_offset_tests.rs,
// `the_rules_of_the_web_host_are_their_zones`).
import { describe, expect, it } from 'vitest';
import { posixTimeZone } from '../src/infra/reactor/time-zone.ts';
import { createWasi } from '../src/infra/reactor/wasi.ts';

const IN_2026 = new Date(Date.UTC(2026, 5, 15, 12));
const IN_2028 = new Date(Date.UTC(2028, 0, 10, 12));

describe('the POSIX rule of a time zone', () => {
  it('is the offset of a zone without summer time, named by the offset', () => {
    expect(posixTimeZone('Asia/Tokyo', IN_2026)).toBe('<+09>-9');
    expect(posixTimeZone('Asia/Kolkata', IN_2026)).toBe('<+0530>-5:30');
    expect(posixTimeZone('UTC', IN_2026)).toBe('<+00>0');
    expect(posixTimeZone('America/Phoenix', IN_2026)).toBe('<-07>7');
  });

  it("holds the year's changes of a zone with summer time: month, week, weekday and the local time before each", () => {
    expect(posixTimeZone('America/New_York', IN_2026)).toBe('<-05>5<-04>4,M3.2.0/2,M11.1.0/2');
    expect(posixTimeZone('Europe/London', IN_2026)).toBe('<+00>0<+01>-1,M3.5.0/1,M10.5.0/2');
    expect(posixTimeZone('Europe/Berlin', IN_2026)).toBe('<+01>-1<+02>-2,M3.5.0/2,M10.5.0/3');
    // The southern hemisphere: summer time at the start and the end of the year.
    expect(posixTimeZone('Australia/Sydney', IN_2026)).toBe('<+10>-10<+11>-11,M10.1.0/2,M4.1.0/3');
    // A half-hour change.
    expect(posixTimeZone('Australia/Lord_Howe', IN_2026)).toBe('<+1030>-10:30<+11>-11,M10.1.0/2,M4.1.0/2');
  });

  it('is made for the year of the date it is asked at', () => {
    // Israel's summer time starts on the Friday before the last Sunday of March: the last Friday of March 2026 (the
    // 27th, M3.5.5), but the fourth of March 2028 (the 24th; the 31st is a Friday after that Sunday, the 26th), so a
    // rule made for another year than the date's misses the change (line NP WW4b1 R2, finding 4, 2026-11-28).
    expect(posixTimeZone('Asia/Jerusalem', IN_2026)).toBe('<+02>-2<+03>-3,M3.5.5/2,M10.5.0/2');
    expect(posixTimeZone('Asia/Jerusalem', IN_2028)).toBe('<+02>-2<+03>-3,M3.4.5/2,M10.5.0/2');
  });

  it("is the browser's zone by default, and UTC for a zone that Intl does not know", () => {
    const zone = Intl.DateTimeFormat().resolvedOptions().timeZone;
    expect(posixTimeZone()).toBe(posixTimeZone(zone));
    expect(posixTimeZone('Not/AZone', IN_2026)).toBe('<+00>0');
  });
});

describe("the engine's environment", () => {
  it("is made each time the engine reads it: a zone changed between two reads is the second's (every instance's WASI, the threaded module's too)", () => {
    // Line NP WW4b1 R2, finding 5, 2026-11-28: the engine reads its environment at each run's start (abi_v2.md §3).
    const memory = new WebAssembly.Memory({ initial: 1 });
    const wasi = createWasi();
    wasi.bind({ exports: { memory } } as unknown as WebAssembly.Instance);
    const imports = wasi.imports as Record<string, (...args: number[]) => number>;
    const read = () => {
      expect(imports['environ_sizes_get']!(0, 4)).toBe(0);
      const view = new DataView(memory.buffer);
      const [count, size] = [view.getUint32(0, true), view.getUint32(4, true)];
      expect(imports['environ_get']!(16, 64)).toBe(0);
      expect(view.getUint32(16, true)).toBe(64);
      return { count, text: new TextDecoder().decode(new Uint8Array(memory.buffer, 64, size)) };
    };
    const saved = process.env['TZ'];
    try {
      process.env['TZ'] = 'Asia/Tokyo';
      expect(read()).toEqual({ count: 1, text: 'TZ=<+09>-9\0' });
      process.env['TZ'] = 'America/Phoenix';
      expect(read()).toEqual({ count: 1, text: 'TZ=<-07>7\0' });
    } finally {
      if (saved === undefined) delete process.env['TZ'];
      else process.env['TZ'] = saved;
    }
  });
});
