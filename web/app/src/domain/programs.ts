// Program descriptors: the only place that lists what differs between programs in the UI
// (plan §2.1 "O"). A descriptor says which engine options the search form shows, in which
// section and under which label. Defaults, choices and help text come from the engine
// (`describe`), not from here (plan §5.3, DW-4); a field whose flag the engine does not
// describe is not shown.
import type { FastaParserKind } from './dataset';

export type ProgramId = 'blastn' | 'blastp' | 'blastx' | 'tblastn' | 'tblastx';
export type SequenceKind = 'nucleotide' | 'protein';
export type InputRole = 'query' | 'subject';

/**
 * How the form edits an option's value. Every kind writes the text as the user chose it
 * (the engine's parser reads it, as on the command line); none checks BLAST's rules.
 * - `text`: free text (numbers too: the engine's parser decides what it accepts).
 * - `choice`: a list; the choices are `describe`'s, or the field's own `choices` where
 *   `describe` writes them only in the help text.
 * - `boolean`: an option that takes `true` or `false`.
 * - `flag`: an option without a value, written alone.
 * - `gencode`: a genetic code; the engine's list (`describe` `query_gencodes` /
 *   `subject_gencodes`) gives the choices.
 */
export type FieldKind = 'text' | 'choice' | 'boolean' | 'flag' | 'gencode';

export interface ParameterField {
  readonly flag: string;
  readonly label: string;
  readonly kind: FieldKind;
  /** For `choice` fields whose choices `describe` gives only in its help text. */
  readonly choices?: readonly string[];
  /** For `text` fields: a hint for the on-screen keyboard. */
  readonly inputMode?: 'numeric' | 'decimal' | 'text';
  /** Suggestions offered with a `text` field (they do not limit what can be typed). */
  readonly suggestions?: readonly string[];
}

export interface ParameterSection {
  readonly title: string;
  readonly fields: readonly ParameterField[];
}

export interface ProgramDescriptor {
  readonly id: ProgramId;
  /** Tab label (the maintainer's decision before S12: the NCBI program names). */
  readonly label: string;
  /** One line under the tabs. */
  readonly summary: string;
  readonly query: SequenceKind;
  readonly subject: SequenceKind;
  /** The engine's FASTA reader for this program, which the index scan follows (plan TD-8). */
  readonly fastaParser: FastaParserKind;
  /** Why the program cannot be searched in this release; undefined when it can. */
  readonly unavailable?: string;
  /** The search form, in display order (as NCBI Web BLAST's "Algorithm parameters"). */
  readonly sections: readonly ParameterSection[];
}

const integer = { kind: 'text', inputMode: 'numeric' } as const;
const decimal = { kind: 'text', inputMode: 'decimal' } as const;
/** Matrix names of NCBI BLAST+; the engine says which ones it searches with. */
const MATRICES = ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250'];

const evalue: ParameterField = { flag: '-evalue', label: 'Expect threshold', ...decimal };
const maxTargetSeqs: ParameterField = { flag: '-max_target_seqs', label: 'Max target sequences', ...integer };
const wordSize: ParameterField = { flag: '-word_size', label: 'Word size', ...integer };
const maxHsps: ParameterField = { flag: '-max_hsps', label: 'Max HSPs per subject', ...integer };
const threshold: ParameterField = { flag: '-threshold', label: 'Neighboring words threshold', ...decimal };
const windowSize: ParameterField = { flag: '-window_size', label: 'Two-hit window size', ...integer };
const matrix: ParameterField = { flag: '-matrix', label: 'Matrix', kind: 'text', suggestions: MATRICES };
const gapOpen: ParameterField = { flag: '-gapopen', label: 'Gap open cost', ...integer };
const gapExtend: ParameterField = { flag: '-gapextend', label: 'Gap extend cost', ...integer };
// `describe` explains the values of -comp_based_stats only in its help text.
const compBasedStats: ParameterField = {
  flag: '-comp_based_stats',
  label: 'Compositional adjustments',
  kind: 'choice',
  choices: ['0', '1', '2', '3'],
};
const seg: ParameterField = { flag: '-seg', label: 'Low-complexity filter (SEG)', kind: 'text', suggestions: ['no', 'yes'] };
const lcaseMasking: ParameterField = { flag: '-lcase_masking', label: 'Mask lower-case letters', kind: 'flag' };

export const PROGRAMS: readonly ProgramDescriptor[] = Object.freeze([
  {
    id: 'blastn',
    label: 'BLASTN',
    summary: 'Nucleotide query against nucleotide subjects',
    query: 'nucleotide',
    subject: 'nucleotide',
    // Kind 0 is the `bio::io::fasta` reader; kind 1 is BLASTX's NCBI-style reader (ABI v2 §4).
    fastaParser: 0,
    sections: [
      {
        title: 'Program selection',
        fields: [
          // `describe` lists the tasks and the template values only in its help text.
          {
            flag: '-task',
            label: 'Task',
            kind: 'choice',
            choices: ['megablast', 'blastn', 'dc-megablast', 'blastn-short'],
          },
          {
            flag: '-template_type',
            label: 'Template type (discontiguous)',
            kind: 'choice',
            choices: ['coding', 'optimal', 'coding_and_optimal'],
          },
          {
            flag: '-template_length',
            label: 'Template length (discontiguous)',
            kind: 'choice',
            choices: ['16', '18', '21'],
          },
        ],
      },
      { title: 'General parameters', fields: [maxTargetSeqs, evalue, wordSize, maxHsps] },
      {
        title: 'Scoring parameters',
        fields: [
          { flag: '-reward', label: 'Match reward', ...integer },
          { flag: '-penalty', label: 'Mismatch penalty', ...integer },
          gapOpen,
          gapExtend,
        ],
      },
      {
        title: 'Filters and masking',
        fields: [
          { flag: '-dust', label: 'Low-complexity filter (DUST)', kind: 'text', suggestions: ['no', 'yes'] },
          lcaseMasking,
        ],
      },
      {
        title: 'Restrict results',
        fields: [
          { flag: '-perc_identity', label: 'Percent identity', ...decimal },
          { flag: '-subject_besthit', label: 'Subject best hit', kind: 'flag' },
        ],
      },
    ],
  },
  {
    id: 'blastp',
    label: 'BLASTP',
    summary: 'Protein query against protein subjects',
    query: 'protein',
    subject: 'protein',
    fastaParser: 0,
    sections: [
      { title: 'Program selection', fields: [{ flag: '-task', label: 'Task', kind: 'choice' }] },
      { title: 'General parameters', fields: [maxTargetSeqs, evalue, wordSize, threshold, windowSize, maxHsps] },
      { title: 'Scoring parameters', fields: [matrix, gapOpen, gapExtend, compBasedStats] },
      { title: 'Filters and masking', fields: [seg] },
    ],
  },
  {
    id: 'blastx',
    label: 'BLASTX',
    summary: 'Translated nucleotide query against protein subjects',
    query: 'nucleotide',
    subject: 'protein',
    fastaParser: 1,
    unavailable: 'BLASTX joins LOSAT Web after its certification. Until then, use BLASTX of the LOSAT command line.',
    sections: [],
  },
  {
    id: 'tblastn',
    label: 'TBLASTN',
    summary: 'Protein query against translated nucleotide subjects',
    query: 'protein',
    subject: 'nucleotide',
    fastaParser: 0,
    sections: [
      { title: 'Program selection', fields: [{ flag: '-task', label: 'Task', kind: 'choice' }] },
      { title: 'General parameters', fields: [maxTargetSeqs, evalue, wordSize, threshold, windowSize] },
      { title: 'Scoring parameters', fields: [matrix, gapOpen, gapExtend, compBasedStats] },
      {
        title: 'Filters and masking',
        fields: [seg, { flag: '-soft_masking', label: 'Mask for lookup table only', kind: 'boolean' }, lcaseMasking],
      },
      { title: 'Genetic code', fields: [{ flag: '-db_gencode', label: 'Subject genetic code', kind: 'gencode' }] },
      {
        title: 'Extension and statistics',
        fields: [
          { flag: '-xdrop_gap', label: 'X-dropoff, preliminary gapped (bits)', ...decimal },
          { flag: '-xdrop_gap_final', label: 'X-dropoff, final gapped (bits)', ...decimal },
          { flag: '-sum_stats', label: 'Sum statistics', kind: 'boolean' },
        ],
      },
    ],
  },
  {
    id: 'tblastx',
    label: 'TBLASTX',
    summary: 'Translated nucleotide query against translated nucleotide subjects',
    query: 'nucleotide',
    subject: 'nucleotide',
    fastaParser: 0,
    sections: [
      {
        title: 'General parameters',
        fields: [
          maxTargetSeqs,
          evalue,
          wordSize,
          threshold,
          windowSize,
          { flag: '-culling_limit', label: 'Culling limit', ...integer },
        ],
      },
      { title: 'Filters and masking', fields: [seg] },
      {
        title: 'Genetic codes',
        fields: [
          { flag: '-query_gencode', label: 'Query genetic code', kind: 'gencode' },
          { flag: '-db_gencode', label: 'Subject genetic code', kind: 'gencode' },
        ],
      },
    ],
  },
] satisfies ProgramDescriptor[]);

export function programById(id: ProgramId): ProgramDescriptor {
  const program = PROGRAMS.find((p) => p.id === id);
  if (program === undefined) throw new Error(`unknown program: ${id}`);
  return program;
}

/** The program's sequence kind for an input role. */
export function sequenceKind(program: ProgramDescriptor, role: InputRole): SequenceKind {
  return role === 'query' ? program.query : program.subject;
}

/** The unit of a position in an input of this kind. */
export function residueUnit(kind: SequenceKind): string {
  return kind === 'nucleotide' ? 'nt' : 'aa';
}
