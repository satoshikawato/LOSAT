// The words for what an input that is not the user's file as given is made of (InputSnapshot.from; line WQA's D22,
// 2026-11-23): the run's Run details, its header and its Outputs say which records a file left out.
import type { InputPart } from './run';

/** The most IDs a sentence names, and a part keeps; the rest are counted. */
export const MOST_LEFT_OUT_IDS = 3;

/** "a", "a and b", "a, b and c". */
function listOf(words: readonly string[]): string {
  return words.length < 2 ? words.join('') : `${words.slice(0, -1).join(', ')} and ${words[words.length - 1]}`;
}

/** "r2", "r2 and r5", "r2, r5, r7 and 4 more": the first IDs of `count` records. */
export function leftOutIds(ids: readonly string[], total: number = ids.length): string {
  const named = ids.slice(0, MOST_LEFT_OUT_IDS);
  const more = total - named.length;
  return more > 0 ? `${named.join(', ')} and ${count(more)} more` : listOf(named);
}

/** A count as the source card writes it: 1,000,010. */
const count = (n: number): string => n.toLocaleString('en-US');

const records = (n: number): string => `${count(n)} ${n === 1 ? 'record' : 'records'}`;

/** "q.fa without r2 (2 of 3 records)"; "b.fa (3 records)" where nothing is left out. */
export function partText(part: InputPart): string {
  if (part.leftOutCount === 0) return `${part.name} (${records(part.records)})`;
  return `${part.name} without ${leftOutIds(part.leftOut, part.leftOutCount)} (${count(part.records)} of ${records(part.records + part.leftOutCount)})`;
}

/** What the input is made of, as `partText` says each part: "a.fa (2 records) and b.fa without b2 (1 of 2 records)". */
export function fromText(parts: readonly InputPart[]): string {
  return listOf(parts.map(partText));
}

/** Whether a part leaves records out. */
export function leavesOut(parts: readonly InputPart[] | undefined): boolean {
  return (parts ?? []).some((part) => part.leftOutCount > 0);
}

/**
 * Whether the page tells what an input is made of: for every file that LOSAT makes (InputSnapshot.from: pasted text, several
 * files, a file with records left out), whose name is the run's own and says nothing of it (DECISIONS 44).
 */
export function worthSaying(parts: readonly InputPart[] | undefined): parts is readonly InputPart[] {
  return parts !== undefined;
}

/**
 * What an input is called in a list of runs (the run picker): the user's file by its name; a file that LOSAT made by what it is
 * made of, not by its long name (`losat-run3-blastn-20261123T101530Z.query.fa`): "pasted text", "a.fa + b.fa", "q.fa, some records
 * left out".
 */
export function inputLabel(input: { readonly name: string; readonly from?: readonly InputPart[] | undefined }): string {
  if (input.from === undefined) return input.name;
  const names = input.from.map((part) => part.name).join(' + ');
  return leavesOut(input.from) ? `${names}, some records left out` : names;
}
