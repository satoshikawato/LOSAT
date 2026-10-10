// Program descriptors: the only place that lists what differs between programs in the UI
// (plan §2.1 "O"). A descriptor says which engine options the search form shows, in which
// section, where on the screen and under which label; the sections and labels are those of
// NCBI BLAST's two-sequence page (docs/web/ncbi_ui_mapping.md §1). Defaults, choices and
// help text come from the engine (`describe`), not from here (plan §5.3, DW-4); a field whose
// flag the engine does not describe is not shown.
import type { Unit } from './coordinates';
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
 * - `radio`: a choice shown as radio buttons (the task, NCBI's "Program Selection"); the
 *   engine's default is checked while nothing is set.
 * - `boolean`: an option that takes `true` or `false`.
 * - `flag`: an option without a value, written alone.
 * - `gencode`: a genetic code; the engine's list (`describe` `query_gencodes` /
 *   `subject_gencodes`) gives the choices.
 */
export type FieldKind = 'text' | 'choice' | 'radio' | 'boolean' | 'flag' | 'gencode';

export interface ParameterField {
  readonly flag: string;
  readonly label: string;
  readonly kind: FieldKind;
  /** For `choice` and `radio` fields whose choices `describe` gives only in its help text. */
  readonly choices?: readonly string[];
  /** Names shown for choices (NCBI's words); the value written is the choice itself. */
  readonly choiceLabels?: Readonly<Record<string, string>>;
  /** For `text` fields: a hint for the on-screen keyboard. */
  readonly inputMode?: 'numeric' | 'decimal' | 'text';
  /** Suggestions offered with a `text` field (they do not limit what can be typed). */
  readonly suggestions?: readonly string[];
  /**
   * The name of a row that shows this field with its neighbours of the same row (NCBI's
   * "Match/Mismatch Scores" and "Gap Costs", which are one select on NCBI's page).
   */
  readonly row?: string;
}

/**
 * Where a section is shown: in "Program Selection", in the query's or the subject's block
 * (the genetic codes), or under "Algorithm parameters".
 */
export type SectionPlacement = 'program' | 'query' | 'subject' | 'algorithm';

export interface ParameterSection {
  readonly title: string;
  readonly placement: SectionPlacement;
  readonly fields: readonly ParameterField[];
  /**
   * A section shown only while this option's value (set, or else the engine's default) is
   * one of `values`, or while one of the section's own fields has a value, so that a value
   * is never written to the argv from a hidden field (NCBI's "Discontiguous Word Options").
   */
  readonly shownFor?: { readonly flag: string; readonly values: readonly string[] };
}

export interface ProgramDescriptor {
  readonly id: ProgramId;
  /** Tab label (the maintainer's decision before S12: the NCBI program names). */
  readonly label: string;
  /** The sentence under the tabs (NCBI's "BLASTN programs search ..."). */
  readonly summary: string;
  /** What the program searches, for the line beside the "Run LOSAT" button ("nucleotide subjects"). */
  readonly searches: string;
  /** The query, where the line beside the "Run LOSAT" button names it ("protein query"). */
  readonly queryNote?: string;
  /** What each task is for, after the task's name in the line beside the "Run LOSAT" button. */
  readonly taskNotes?: Readonly<Record<string, string>>;
  readonly query: SequenceKind;
  readonly subject: SequenceKind;
  /** Why the program cannot be searched in this release; undefined when it can. */
  readonly unavailable?: string;
  /**
   * The search form. The order of the sections and fields is the order of the argv: the
   * task first, as on a command line, then the genetic codes and the algorithm parameters
   * in the order of the page.
   */
  readonly sections: readonly ParameterSection[];
}

const integer = { kind: 'text', inputMode: 'numeric' } as const;
const decimal = { kind: 'text', inputMode: 'decimal' } as const;
/** Matrix names of NCBI BLAST+; the engine says which ones it searches with. */
const MATRICES = ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250'];

// General Parameters.
const maxTargetSeqs: ParameterField = { flag: '-max_target_seqs', label: 'Max target sequences', ...integer };
const evalue: ParameterField = { flag: '-evalue', label: 'Expect threshold', ...decimal };
const wordSize: ParameterField = { flag: '-word_size', label: 'Word size', ...integer };
// Scoring Parameters: NCBI's "Gap Costs" select is one row of the two costs here.
const matrix: ParameterField = { flag: '-matrix', label: 'Matrix', kind: 'text', suggestions: MATRICES };
const gapOpen: ParameterField = { flag: '-gapopen', label: 'Existence', row: 'Gap Costs', ...integer };
const gapExtend: ParameterField = { flag: '-gapextend', label: 'Extension', row: 'Gap Costs', ...integer };
// `describe` explains the values of -comp_based_stats only in its help text; the names are NCBI's.
const compBasedStats: ParameterField = {
  flag: '-comp_based_stats',
  label: 'Compositional adjustments',
  kind: 'choice',
  choices: ['0', '1', '2', '3'],
  choiceLabels: {
    '0': '0: No adjustment',
    '1': '1: Composition-based statistics',
    '2': '2: Conditional compositional score matrix adjustment',
    '3': '3: Universal compositional score matrix adjustment',
  },
};
// Filters and Masking.
const seg: ParameterField = { flag: '-seg', label: 'Low complexity regions filter (SEG)', kind: 'text', suggestions: ['no', 'yes'] };
const lcaseMasking: ParameterField = { flag: '-lcase_masking', label: 'Mask lower case letters', kind: 'flag' };
// Other Parameters: LOSAT's options that NCBI's page does not show.
const maxHsps: ParameterField = { flag: '-max_hsps', label: 'Max HSPs per subject', ...integer };
const threshold: ParameterField = { flag: '-threshold', label: 'Neighboring words threshold', ...decimal };
const windowSize: ParameterField = { flag: '-window_size', label: 'Two-hit window size', ...integer };

const OTHER = 'Other Parameters';

export const PROGRAMS: readonly ProgramDescriptor[] = Object.freeze([
  {
    id: 'blastn',
    label: 'BLASTN',
    summary: 'BLASTN searches nucleotide subjects using a nucleotide query.',
    searches: 'nucleotide subjects',
    taskNotes: {
      megablast: 'optimized for highly similar sequences',
      'dc-megablast': 'optimized for more dissimilar sequences',
      blastn: 'optimized for somewhat similar sequences',
      'blastn-short': 'optimized for short sequences',
    },
    query: 'nucleotide',
    subject: 'nucleotide',
    sections: [
      {
        title: 'Program Selection',
        placement: 'program',
        fields: [
          // `describe` lists the tasks only in its help text. LOSAT adds blastn-short, which
          // NCBI's page replaces with its "Short queries" adjustment.
          {
            flag: '-task',
            label: 'Optimize for',
            kind: 'radio',
            choices: ['megablast', 'dc-megablast', 'blastn', 'blastn-short'],
            choiceLabels: {
              megablast: 'Highly similar sequences (megablast)',
              'dc-megablast': 'More dissimilar sequences (discontiguous megablast)',
              blastn: 'Somewhat similar sequences (blastn)',
              'blastn-short': 'Short sequences (blastn-short)',
            },
          },
        ],
      },
      { title: 'General Parameters', placement: 'algorithm', fields: [maxTargetSeqs, evalue, wordSize] },
      {
        title: 'Scoring Parameters',
        placement: 'algorithm',
        fields: [
          { flag: '-reward', label: 'Match', row: 'Match/Mismatch Scores', ...integer },
          { flag: '-penalty', label: 'Mismatch', row: 'Match/Mismatch Scores', ...integer },
          gapOpen,
          gapExtend,
        ],
      },
      {
        title: 'Filters and Masking',
        placement: 'algorithm',
        fields: [
          { flag: '-dust', label: 'Low complexity regions filter (DUST)', kind: 'text', suggestions: ['no', 'yes'] },
          lcaseMasking,
        ],
      },
      {
        title: 'Discontiguous Word Options',
        placement: 'algorithm',
        shownFor: { flag: '-task', values: ['dc-megablast'] },
        // `describe` lists the template values only in its help text.
        fields: [
          { flag: '-template_length', label: 'Template length', kind: 'choice', choices: ['16', '18', '21'] },
          { flag: '-template_type', label: 'Template type', kind: 'choice', choices: ['coding', 'optimal', 'coding_and_optimal'] },
        ],
      },
      {
        title: OTHER,
        placement: 'algorithm',
        fields: [
          maxHsps,
          { flag: '-perc_identity', label: 'Percent identity', ...decimal },
          { flag: '-subject_besthit', label: 'Subject best hit', kind: 'flag' },
        ],
      },
    ],
  },
  {
    id: 'blastp',
    label: 'BLASTP',
    summary: 'BLASTP searches protein subjects using a protein query.',
    searches: 'protein subjects',
    taskNotes: { blastp: 'protein-protein BLAST', 'blastp-fast': 'Quick BLASTP', 'blastp-short': 'short queries' },
    query: 'protein',
    subject: 'protein',
    sections: [
      {
        title: 'Program Selection',
        placement: 'program',
        fields: [
          {
            flag: '-task',
            label: 'Algorithm',
            kind: 'radio',
            choiceLabels: {
              blastp: 'blastp (protein-protein BLAST)',
              'blastp-fast': 'Quick BLASTP (blastp-fast)',
              'blastp-short': 'blastp-short (short queries)',
            },
          },
        ],
      },
      { title: 'General Parameters', placement: 'algorithm', fields: [maxTargetSeqs, evalue, wordSize] },
      { title: 'Scoring Parameters', placement: 'algorithm', fields: [matrix, gapOpen, gapExtend, compBasedStats] },
      { title: 'Filters and Masking', placement: 'algorithm', fields: [seg] },
      { title: OTHER, placement: 'algorithm', fields: [threshold, windowSize, maxHsps] },
    ],
  },
  {
    id: 'blastx',
    label: 'BLASTX',
    summary: 'BLASTX searches protein subjects using a translated nucleotide query.',
    searches: 'protein subjects',
    queryNote: 'translated nucleotide query',
    query: 'nucleotide',
    subject: 'protein',
    unavailable: 'BLASTX joins LOSAT Web after its certification. Until then, use BLASTX of the LOSAT command line.',
    sections: [],
  },
  {
    id: 'tblastn',
    label: 'TBLASTN',
    summary: 'TBLASTN searches translated nucleotide subjects using a protein query.',
    searches: 'translated nucleotide subjects',
    queryNote: 'protein query',
    query: 'protein',
    subject: 'nucleotide',
    sections: [
      { title: 'Program Selection', placement: 'program', fields: [{ flag: '-task', label: 'Algorithm', kind: 'radio' }] },
      { title: 'Genetic code', placement: 'subject', fields: [{ flag: '-db_gencode', label: 'Genetic code', kind: 'gencode' }] },
      { title: 'General Parameters', placement: 'algorithm', fields: [maxTargetSeqs, evalue, wordSize] },
      { title: 'Scoring Parameters', placement: 'algorithm', fields: [matrix, gapOpen, gapExtend, compBasedStats] },
      {
        title: 'Filters and Masking',
        placement: 'algorithm',
        fields: [seg, { flag: '-soft_masking', label: 'Mask for lookup table only', kind: 'boolean' }, lcaseMasking],
      },
      {
        title: OTHER,
        placement: 'algorithm',
        fields: [
          threshold,
          windowSize,
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
    summary: 'TBLASTX searches translated nucleotide subjects using a translated nucleotide query.',
    searches: 'translated nucleotide subjects',
    queryNote: 'translated nucleotide query',
    query: 'nucleotide',
    subject: 'nucleotide',
    sections: [
      { title: 'Genetic code', placement: 'query', fields: [{ flag: '-query_gencode', label: 'Genetic code', kind: 'gencode' }] },
      { title: 'Genetic code', placement: 'subject', fields: [{ flag: '-db_gencode', label: 'Genetic code', kind: 'gencode' }] },
      {
        title: 'General Parameters',
        placement: 'algorithm',
        // NCBI's "Max matches in a query range" (HSP_RANGE_MAX) is BLAST+'s -culling_limit.
        fields: [maxTargetSeqs, evalue, wordSize, { flag: '-culling_limit', label: 'Max matches in a query range', ...integer }],
      },
      { title: 'Filters and Masking', placement: 'algorithm', fields: [seg] },
      { title: OTHER, placement: 'algorithm', fields: [threshold, windowSize] },
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

/**
 * The index scan's reader for an input role of a program: the engine's reader with the flags
 * of the role's sequence kind (abi_v2.md §4, §9), kind 1 for nucleotide input and kind 2 for
 * protein input. BLASTX, which cannot be searched yet, has its kinds as well.
 */
export function indexParser(program: ProgramId, role: InputRole): FastaParserKind {
  return sequenceKind(programById(program), role) === 'nucleotide' ? 1 : 2;
}

/** The unit of a position in an input of this kind (domain/coordinates.ts `Unit`). */
export function residueUnit(kind: SequenceKind): Unit {
  return kind === 'nucleotide' ? 'nt' : 'aa';
}

/**
 * The line beside the "Run LOSAT" button (NCBI's "Search nucleotide sequence using Megablast
 * (Optimize for highly similar sequences)"): what is searched, with which program and task,
 * and where. `task` is the task the search uses (set, or the engine's default), if known.
 */
export function searchSummary(program: ProgramDescriptor, task: string | undefined): string {
  const query = program.queryNote === undefined ? '' : ` (${program.queryNote})`;
  const note = task === undefined ? undefined : program.taskNotes?.[task];
  const taskText = task === undefined ? '' : `, task ${task}${note === undefined ? '' : ` (${note})`}`;
  return `Search ${program.searches} with ${program.label}${query}${taskText}. Runs in this browser.`;
}
