// The words for what an input that is not the user's file as given is made of (InputSnapshot.from; line WQA's D22,
// 2026-11-23): the run's Run details, its header and its Outputs say which records a file left out.
import type { InputPart } from './run';

/** The most IDs a sentence names; the rest are counted. */
const MOST_IDS = 3;

/** "r2", "r2 and r5", "r2, r5, r7 and 4 more". */
export function leftOutIds(ids: readonly string[]): string {
  const named = ids.slice(0, MOST_IDS);
  const more = ids.length - named.length;
  if (more > 0) return `${named.join(', ')} and ${more} more`;
  return named.length < 2 ? named.join('') : `${named.slice(0, -1).join(', ')} and ${named[named.length - 1]}`;
}

const records = (count: number): string => `${count} ${count === 1 ? 'record' : 'records'}`;

/** "q.fa: 2 of 3 records (r2 left out)"; "b.fa: 3 records" where nothing is left out. */
export function partText(part: InputPart): string {
  if (part.leftOut.length === 0) return `${part.name}: ${records(part.records)}`;
  return `${part.name}: ${part.records} of ${records(part.records + part.leftOut.length)} (${leftOutIds(part.leftOut)} left out)`;
}

/** The parts that leave records out, as `partText` says them. */
export function leftOutTexts(parts: readonly InputPart[] | undefined): readonly string[] {
  return (parts ?? []).filter((part) => part.leftOut.length > 0).map(partText);
}

/** What the input is made of, as `partText` says each part, for an input that is not the user's file as given. */
export function fromText(parts: readonly InputPart[]): string {
  return parts.map(partText).join('; ');
}
