// The engine's local time zone (docs/web/abi_v2.md §3; line NP WW4b1, D-WW4b-2, 2026-11-27): the host puts `TZ` into WASI's
// environment, so that the dates that makeblastdb writes into a database (`.pin`/`.nin`, the alias files, the mask and
// GI-mask files) are the user's local time, as when the user runs the command line on their machine. The engine reads
// `TZ` as glibc does (LOSAT/src/writedb/local_offset.rs): it finds no zone file on WASI and reads the value as a POSIX
// rule, with jiff. The browser gives the zone (`Intl`) but not its rules, so the host makes the rule of the current year
// from the zone's offsets: every day of the year is looked at, each change of the offset is found to the second, and the
// two changes of a zone with summer time become the rule's start and end, as month, week and weekday (`Mm.w.d`, w = 5 for
// the last such weekday of the month) and the local time before the change. A zone of one offset all year is that
// offset (`<+09>-9`); the names are the offsets (`<+0530>`), as the dates carry no zone name. A zone whose offset
// changes more than twice in the year (Africa/Casablanca, around Ramadan) gets the offset that it has now.

/** One day, in milliseconds. */
const DAY_MS = 86_400_000;

/** The offset from UTC (seconds east) of `zone` at each instant, read from the browser's `Intl`. */
function offsetsOf(zone: string): (ms: number) => number {
  const format = new Intl.DateTimeFormat('en-US', {
    timeZone: zone,
    hourCycle: 'h23',
    year: 'numeric',
    month: 'numeric',
    day: 'numeric',
    hour: 'numeric',
    minute: 'numeric',
    second: 'numeric',
  });
  return (ms) => {
    const parts: Record<string, number> = {};
    for (const part of format.formatToParts(new Date(ms))) if (part.type !== 'literal') parts[part.type] = Number(part.value);
    const local = Date.UTC(parts['year']!, parts['month']! - 1, parts['day']!, parts['hour']! % 24, parts['minute']!, parts['second']!);
    return Math.round((local - Math.floor(ms / 1000) * 1000) / 1000);
  };
}

/** `seconds` east of UTC as POSIX writes an offset: west of UTC is positive, `[-]h[:mm[:ss]]`. */
function posixOffset(seconds: number): string {
  const west = -seconds;
  const sign = west < 0 ? '-' : '';
  const value = Math.abs(west);
  const hours = Math.floor(value / 3600);
  const minutes = Math.floor((value % 3600) / 60);
  const rest = value % 60;
  const pad = (n: number) => String(n).padStart(2, '0');
  return `${sign}${hours}${minutes === 0 && rest === 0 ? '' : `:${pad(minutes)}`}${rest === 0 ? '' : `:${pad(rest)}`}`;
}

/** The name of an offset: `<+09>`, `<-05>`, `<+0530>` (POSIX's quoted form, which jiff and glibc read). */
function offsetName(seconds: number): string {
  const sign = seconds < 0 ? '-' : '+';
  const value = Math.abs(seconds);
  const pad = (n: number) => String(n).padStart(2, '0');
  const minutes = Math.floor((value % 3600) / 60);
  return `<${sign}${pad(Math.floor(value / 3600))}${minutes === 0 ? '' : pad(minutes)}>`;
}

/** The POSIX date and time of a change at `ms` (UTC) from the offset `before`: `Mm.w.d/time` in the local time before it. */
function changeRule(ms: number, before: number): string {
  const local = new Date(ms + before * 1000);
  const month = local.getUTCMonth();
  const day = local.getUTCDate();
  const daysInMonth = new Date(Date.UTC(local.getUTCFullYear(), month + 1, 0)).getUTCDate();
  const week = day + 7 > daysInMonth ? 5 : Math.ceil(day / 7);
  const seconds = local.getUTCHours() * 3600 + local.getUTCMinutes() * 60 + local.getUTCSeconds();
  const pad = (n: number) => String(n).padStart(2, '0');
  const time = `${Math.floor(seconds / 3600)}${seconds % 3600 === 0 ? '' : `:${pad(Math.floor((seconds % 3600) / 60))}`}${seconds % 60 === 0 ? '' : `:${pad(seconds % 60)}`}`;
  return `M${month + 1}.${week}.${local.getUTCDay()}/${time}`;
}

/**
 * The POSIX rule (`TZ`) of the time zone `zone` in the year of `now`: `<+09>-9`, `<-05>5<-04>4,M3.2.0/2,M11.1.0/2`. The
 * zone is the browser's by default; a zone that `Intl` does not know is UTC.
 */
export function posixTimeZone(zone?: string, now: Date = new Date()): string {
  let offsetAt: (ms: number) => number;
  try {
    offsetAt = offsetsOf(zone ?? Intl.DateTimeFormat().resolvedOptions().timeZone);
  } catch {
    offsetAt = () => 0;
  }
  const current = offsetAt(now.getTime());
  const year = new Date(now.getTime() + current * 1000).getUTCFullYear();
  // The changes whose local time before them is in the year: each found to the second between two days.
  const changes: Array<{ readonly ms: number; readonly before: number; readonly after: number }> = [];
  const first = Date.UTC(year, 0, 1) - DAY_MS;
  const last = Date.UTC(year + 1, 0, 1) + DAY_MS;
  let before = offsetAt(first);
  for (let day = first + DAY_MS; day <= last; day += DAY_MS) {
    const after = offsetAt(day);
    if (after === before) continue;
    let low = day - DAY_MS;
    let high = day;
    while (high - low > 1000) {
      const middle = low + Math.floor((high - low) / 2000) * 1000;
      if (offsetAt(middle) === before) low = middle;
      else high = middle;
    }
    if (new Date(high + before * 1000).getUTCFullYear() === year) changes.push({ ms: high, before, after });
    before = after;
  }
  if (changes.length !== 2 || changes[0]!.after !== changes[1]!.before || changes[1]!.after !== changes[0]!.before) {
    return `${offsetName(current)}${posixOffset(current)}`;
  }
  const standard = Math.min(changes[0]!.before, changes[0]!.after);
  const summer = Math.max(changes[0]!.before, changes[0]!.after);
  const start = changes.find((change) => change.after === summer)!;
  const end = changes.find((change) => change.after === standard)!;
  return `${offsetName(standard)}${posixOffset(standard)}${offsetName(summer)}${posixOffset(summer)},${changeRule(start.ms, start.before)},${changeRule(end.ms, end.before)}`;
}
