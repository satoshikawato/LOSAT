// Program descriptors: the only place that lists what differs between programs in
// the UI. Defaults, choices and help text come from the engine (`describe`), not here.

export type ProgramId = 'blastn' | 'blastp' | 'blastx' | 'tblastn' | 'tblastx';
export type SequenceKind = 'nucleotide' | 'protein';

export interface ProgramDescriptor {
  readonly id: ProgramId;
  /** Tab label shown to users. The final display names are still to be decided. */
  readonly label: string;
  readonly query: SequenceKind;
  readonly subject: SequenceKind;
}

export const PROGRAMS: readonly ProgramDescriptor[] = Object.freeze([
  { id: 'blastn', label: 'BLASTN', query: 'nucleotide', subject: 'nucleotide' },
  { id: 'blastp', label: 'BLASTP', query: 'protein', subject: 'protein' },
  { id: 'blastx', label: 'BLASTX', query: 'nucleotide', subject: 'protein' },
  { id: 'tblastn', label: 'TBLASTN', query: 'protein', subject: 'nucleotide' },
  { id: 'tblastx', label: 'TBLASTX', query: 'nucleotide', subject: 'nucleotide' },
]);

export function programById(id: ProgramId): ProgramDescriptor {
  const program = PROGRAMS.find((p) => p.id === id);
  if (program === undefined) throw new Error(`unknown program: ${id}`);
  return program;
}
