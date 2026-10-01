// Program descriptors: the only place that lists what differs between programs in
// the UI. Defaults, choices and help text come from the engine (`describe`), not here.
import type { FastaParserKind } from './dataset';

export type ProgramId = 'blastn' | 'blastp' | 'blastx' | 'tblastn' | 'tblastx';
export type SequenceKind = 'nucleotide' | 'protein';

export interface ProgramDescriptor {
  readonly id: ProgramId;
  /** Tab label shown to users. The final display names are still to be decided. */
  readonly label: string;
  readonly query: SequenceKind;
  readonly subject: SequenceKind;
  /** The engine's FASTA reader for this program, which the index scan follows (plan TD-8). */
  readonly fastaParser: FastaParserKind;
}

export const PROGRAMS: readonly ProgramDescriptor[] = Object.freeze([
  // Kind 0 is the `bio::io::fasta` reader; kind 1 is BLASTX's NCBI-style reader (ABI v2 §4).
  { id: 'blastn', label: 'BLASTN', query: 'nucleotide', subject: 'nucleotide', fastaParser: 0 },
  { id: 'blastp', label: 'BLASTP', query: 'protein', subject: 'protein', fastaParser: 0 },
  { id: 'blastx', label: 'BLASTX', query: 'nucleotide', subject: 'protein', fastaParser: 1 },
  { id: 'tblastn', label: 'TBLASTN', query: 'protein', subject: 'nucleotide', fastaParser: 0 },
  { id: 'tblastx', label: 'TBLASTX', query: 'nucleotide', subject: 'nucleotide', fastaParser: 0 },
]);

export function programById(id: ProgramId): ProgramDescriptor {
  const program = PROGRAMS.find((p) => p.id === id);
  if (program === undefined) throw new Error(`unknown program: ${id}`);
  return program;
}
