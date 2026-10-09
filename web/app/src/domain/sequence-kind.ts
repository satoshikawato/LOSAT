// The application's estimate of whether a record is a nucleotide or a protein sequence
// (plan §5.4). It is a warning only: the program and its inputs never change because of it
// (design §3.1), and whether the engine can read a record is the engine's decision (the
// input check). The estimate counts the letters of the index scan's `residue_counts`.
import type { SequenceKind } from './programs';

/** Nucleotide letters, A C G T U and N, as the share of all letters that marks a nucleotide sequence. */
const NUCLEOTIDE_LETTERS = new Set(['A', 'C', 'G', 'T', 'U', 'N']);
const NUCLEOTIDE_SHARE = 0.9;

/** The estimated kind, or undefined for a record without letters. */
export function estimateKind(residueCounts: Readonly<Record<string, number>>): SequenceKind | undefined {
  let letters = 0;
  let nucleotides = 0;
  for (const [key, count] of Object.entries(residueCounts)) {
    if (!/^[A-Za-z]$/.test(key)) continue;
    letters += count;
    if (NUCLEOTIDE_LETTERS.has(key.toUpperCase())) nucleotides += count;
  }
  if (letters === 0) return undefined;
  return nucleotides >= NUCLEOTIDE_SHARE * letters ? 'nucleotide' : 'protein';
}

/** Whether the estimate of a record differs from the kind the program reads. */
export function looksLikeOtherKind(residueCounts: Readonly<Record<string, number>>, expected: SequenceKind): boolean {
  const kind = estimateKind(residueCounts);
  return kind !== undefined && kind !== expected;
}
