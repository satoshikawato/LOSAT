// The engine's programs and the kind of sequence that each reads as its query and its subject.
// What the search form shows for each program is the application's (web/app/src/domain/programs.ts).

export type ProgramId = 'blastn' | 'blastp' | 'blastx' | 'tblastn' | 'tblastx';
export type SequenceKind = 'nucleotide' | 'protein';
export type InputRole = 'query' | 'subject';

/** The kinds of sequence that a program reads. */
export interface ProgramInputs {
  readonly query: SequenceKind;
  readonly subject: SequenceKind;
}

export const PROGRAM_INPUTS: Readonly<Record<ProgramId, ProgramInputs>> = Object.freeze({
  blastn: { query: 'nucleotide', subject: 'nucleotide' },
  blastp: { query: 'protein', subject: 'protein' },
  blastx: { query: 'nucleotide', subject: 'protein' },
  tblastn: { query: 'protein', subject: 'nucleotide' },
  tblastx: { query: 'nucleotide', subject: 'nucleotide' },
});
