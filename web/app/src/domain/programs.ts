// Program descriptors: the only place that lists what differs between programs in the UI
// (plan §2.1 "O"). A descriptor says which engine options the search form shows, in which
// section, where on the screen and under which label; the sections and labels are those of
// NCBI BLAST's two-sequence page (docs/web/ncbi_ui_mapping.md §1). Defaults, choices and
// help text come from the engine (`describe`), not from here (plan §5.3, DW-4); a field whose
// flag the engine does not describe is not shown.
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
 * - `radio`: a choice shown as radio buttons (the task, in "Task Selection"); the
 *   engine's default is checked while nothing is set.
 * - `boolean`: an option that takes `true` or `false`.
 * - `flag`: an option without a value, written alone.
 * - `gencode`: a genetic code; the engine's list (`describe` `query_gencodes` /
 *   `subject_gencodes`) gives the choices.
 * - `filter`: a low complexity filter (`no`, `yes` or the filter's arguments): NCBI's "Low
 *   complexity regions" checkbox, with the value in a text box while it is checked.
 */
export type FieldKind = 'text' | 'choice' | 'radio' | 'boolean' | 'flag' | 'gencode' | 'filter';

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
  /**
   * For `text` fields: the values of NCBI's select, offered in a select (one list, or a list per
   * task); its "Other" takes any other value in a text box, which the engine's parser reads.
   */
  readonly list?: readonly string[] | Readonly<Record<string, readonly string[]>>;
  /** NCBI's sentence beside the control ("Select the maximum number of aligned sequences to display"). */
  readonly hint?: string;
  /**
   * The name of a row that shows this field with its neighbours of the same row (NCBI's
   * "Match/Mismatch Scores", "Gap Costs", "Filter" and "Mask").
   */
  readonly row?: string;
}

/**
 * NCBI's select of a row of two fields, in NCBI's order: "Match/Mismatch Scores" (1,-2 writes
 * -reward 1 -penalty -2) and "Gap Costs" ("Existence: 11 Extension: 1" writes -gapopen 11
 * -gapextend 1). Its "Other" shows the two fields. The values are those of NCBI's page; which
 * pairs the engine reads is still its `validate`'s and its search's to say.
 */
export interface PairSelect {
  /** NCBI's id of the select (`matchscores`, `gapcosts`). */
  readonly id: string;
  /** What separates the two values in NCBI's values ("1,-2", "11 1"). */
  readonly separator: ',' | ' ';
  /**
   * The options whose values (as the search uses them, joined with ",") choose the list: the
   * matrix, or BLASTN's scores. One list (key '') without.
   */
  readonly by?: readonly string[];
  readonly lists: Readonly<Record<string, readonly string[]>>;
  /** NCBI's default of each list (its "Deflt"), which a change of the list chooses. */
  readonly defaults?: Readonly<Record<string, string>>;
  /** Values offered for these tasks only, and then the default of their list (NCBI's "Linear" gap costs of megablast). */
  readonly only?: { readonly tasks: readonly string[]; readonly values: readonly string[] };
  /** NCBI's words for a value ("Existence: 11 Extension: 1"); the value itself without. */
  readonly labels?: (value: string) => string;
}

/**
 * Where a section is shown: in "Task Selection", in the query's or the subject's block
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
  /** NCBI's selects of the section's rows of two fields, by row name. */
  readonly pairs?: Readonly<Record<string, PairSelect>>;
  /**
   * Options that BLAST reads only together (each requires the others: NCBI's "Discontiguous
   * Word Options"), with the values of NCBI's page, which are the engine's for the task that
   * shows the section: a value set in one writes the others too, at their values set or else
   * these, so that the argv never holds one alone.
   */
  readonly together?: Readonly<Record<string, string>>;
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
  /**
   * The search form. The order of the sections and fields is the order of the argv: the
   * task first, as on a command line, then the genetic codes and the algorithm parameters
   * in the order of the page.
   */
  readonly sections: readonly ParameterSection[];
}

const integer = { kind: 'text', inputMode: 'numeric' } as const;
const decimal = { kind: 'text', inputMode: 'decimal' } as const;
/** NCBI's matrix select, in its order; "Other" takes the other names that the engine reads. */
export const MATRICES: readonly string[] = ['PAM30', 'PAM70', 'BLOSUM90', 'BLOSUM80', 'BLOSUM62', 'BLOSUM50', 'BLOSUM45', 'PAM250'];
/** NCBI's word size selects of the protein programs. */
export const PROTEIN_WORD_SIZES: readonly string[] = ['2', '3', '5', '6'];
/** NCBI's compositional adjustments. */
export const COMPOSITION_MODES: readonly string[] = ['0', '1', '2', '3'];

// General Parameters.
const maxTargetSeqs: ParameterField = {
  flag: '-max_target_seqs',
  label: 'Max target sequences',
  ...integer,
  list: ['10', '50', '100', '250', '500', '1000', '5000'],
  hint: 'Select the maximum number of aligned sequences to display',
};
const evalue: ParameterField = { flag: '-evalue', label: 'Expect threshold', ...decimal };
const wordSize: ParameterField = { flag: '-word_size', label: 'Word size', ...integer, list: PROTEIN_WORD_SIZES };
// NCBI's "Max matches in a query range" (HSP_RANGE_MAX) is BLAST+'s -culling_limit.
const cullingLimit: ParameterField = { flag: '-culling_limit', label: 'Max matches in a query range', ...integer };
// Scoring Parameters: NCBI's "Gap Costs" select writes the two costs of its row.
const matrix: ParameterField = { flag: '-matrix', label: 'Matrix', kind: 'text', list: MATRICES };
const gapOpen: ParameterField = { flag: '-gapopen', label: 'Existence', row: 'Gap Costs', ...integer };
const gapExtend: ParameterField = { flag: '-gapextend', label: 'Extension', row: 'Gap Costs', ...integer };
const gapCostsLabel = (value: string): string => {
  const [existence, extension] = value.split(' ');
  return value === '0 0' ? 'Linear' : `Existence: ${existence} Extension: ${extension}`;
};
// NCBI's lists of BLASTN's two-sequence page (gapcosts_1-2 ...): each pair of scores has its own,
// and megablast's "Linear" (0 0) is the default where it is listed.
const MATCH_SCORES: PairSelect = { id: 'matchscores', separator: ',', lists: { '': ['1,-2', '1,-3', '1,-4', '2,-3', '4,-5', '1,-1'] } };
const NUCLEOTIDE_GAP_COSTS: PairSelect = {
  id: 'gapcosts',
  separator: ' ',
  by: ['-reward', '-penalty'],
  lists: {
    '1,-2': ['0 0', '5 2', '2 2', '1 2', '0 2', '3 1', '2 1', '1 1'],
    '1,-3': ['0 0', '5 2', '2 2', '1 2', '0 2', '2 1', '1 1'],
    '1,-4': ['0 0', '5 2', '1 2', '0 2', '2 1', '1 1'],
    '2,-3': ['0 0', '4 4', '2 4', '0 4', '3 3', '6 2', '5 2', '4 2', '2 2'],
    '4,-5': ['0 0', '12 8', '6 5', '5 5', '4 5', '3 5'],
    '1,-1': ['5 2', '3 2', '2 2', '1 2', '0 2', '4 1', '3 1', '2 1'],
  },
  defaults: { '1,-2': '5 2', '1,-3': '5 2', '1,-4': '5 2', '2,-3': '5 2', '4,-5': '12 8', '1,-1': '5 2' },
  only: { tasks: ['megablast'], values: ['0 0'] },
  labels: gapCostsLabel,
};
// NCBI's gap costs of BLOSUM62, all of which LOSAT's BLASTP and TBLASTN read.
const BLOSUM62_GAP_COSTS = ['11 2', '10 2', '9 2', '8 2', '7 2', '6 2', '13 1', '12 1', '11 1', '10 1', '9 1'];
// NCBI's lists of the protein programs' pages (the same on BLASTP's, TBLASTN's and TBLASTX's), by matrix.
export const PROTEIN_GAP_COSTS: PairSelect = {
  id: 'gapcosts',
  separator: ' ',
  by: ['-matrix'],
  lists: {
    BLOSUM62: BLOSUM62_GAP_COSTS,
    PAM30: ['7 2', '6 2', '5 2', '10 1', '9 1', '8 1', '13 3', '15 3', '14 1', '14 2'],
    PAM70: ['8 2', '7 2', '6 2', '11 1', '10 1', '9 1', '12 3', '11 2'],
    BLOSUM80: ['8 2', '7 2', '6 2', '11 1', '10 1', '9 1'],
    BLOSUM45: ['13 3', '12 3', '11 3', '10 3', '15 2', '14 2', '13 2', '12 2', '19 1', '18 1', '17 1', '16 1'],
    BLOSUM50: ['13 3', '12 3', '11 3', '10 3', '9 3', '16 2', '15 2', '14 2', '13 2', '12 2', '19 1', '18 1', '17 1', '16 1', '15 1'],
    PAM250: ['15 3', '14 3', '13 3', '12 3', '11 3', '17 2', '16 2', '15 2', '14 2', '13 2', '21 1', '20 1', '19 1', '18 1', '17 1'],
    BLOSUM90: ['9 2', '8 2', '7 2', '6 2', '11 1', '10 1', '9 1'],
  },
  defaults: {
    BLOSUM62: '11 1',
    PAM30: '9 1',
    PAM70: '10 1',
    BLOSUM80: '10 1',
    BLOSUM45: '14 2',
    BLOSUM50: '13 2',
    PAM250: '15 2',
    BLOSUM90: '10 1',
  },
  labels: gapCostsLabel,
};
const PROTEIN_SCORING: Pick<ParameterSection, 'pairs'> = { pairs: { 'Gap Costs': PROTEIN_GAP_COSTS } };
// LOSAT's BLASTP scores with NCBI's matrices and the gap costs of their lists (IDENTITY, which
// NCBI's page does not list, is refused); LOSAT's BLASTX with BLOSUM62 and gap costs 11/1 only;
// both adjust with modes 0 and 2: their lists keep these of NCBI's values ("Other" still takes any
// matrix and costs, which the engine refuses with its message; Compositional adjustments has no
// "Other"). A series that makes the engine read more widens them (the engine's e2e of the selects
// fails until then).
const BLOSUM62_ONLY: Pick<ParameterSection, 'pairs'> = {
  pairs: { 'Gap Costs': { id: 'gapcosts', separator: ' ', lists: { '': ['11 1'] }, labels: gapCostsLabel } },
};
// `describe` explains the values of -comp_based_stats only in its help text; the names are NCBI's.
const compBasedStats: ParameterField = {
  flag: '-comp_based_stats',
  label: 'Compositional adjustments',
  kind: 'choice',
  choices: COMPOSITION_MODES,
  choiceLabels: {
    '0': '0: No adjustment',
    '1': '1: Composition-based statistics',
    '2': '2: Conditional compositional score matrix adjustment',
    '3': '3: Universal compositional score matrix adjustment',
  },
};
// Filters and Masking: NCBI's "Filter" row, its "Low complexity regions" checkbox.
const seg: ParameterField = { flag: '-seg', label: 'Low complexity regions', kind: 'filter', row: 'Filter' };
// NCBI's "Mask" row: two checkboxes.
const softMasking: ParameterField = { flag: '-soft_masking', label: 'Mask for lookup table only', kind: 'boolean', row: 'Mask' };
const lcaseMasking: ParameterField = { flag: '-lcase_masking', label: 'Mask lower case letters', kind: 'flag', row: 'Mask' };
// Other Parameters: LOSAT's options that NCBI's page does not show.
const maxHsps: ParameterField = { flag: '-max_hsps', label: 'Max HSPs per subject', ...integer };
const threshold: ParameterField = { flag: '-threshold', label: 'Neighboring words threshold', ...decimal };
const windowSize: ParameterField = { flag: '-window_size', label: 'Two-hit window size', ...integer };

// BLAST+'s -num_descriptions / -num_alignments (NCBI's page has them as "Descriptions" / "Alignments" of its
// display settings; they also set the hit list sizes, so a run may search twice, abi_v2.md §7).
const numDescriptions: ParameterField = { flag: '-num_descriptions', label: 'Number of descriptions', ...integer };
const numAlignments: ParameterField = { flag: '-num_alignments', label: 'Number of alignments', ...integer };

const OTHER = 'Other Parameters';

export const PROGRAMS: readonly ProgramDescriptor[] = Object.freeze([
  {
    id: 'blastn',
    label: 'BLASTN',
    summary: 'BLASTN mode searches nucleotide subjects using a nucleotide query.',
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
        title: 'Task Selection',
        placement: 'program',
        fields: [
          // `describe` lists the tasks in the order of NCBI's page. LOSAT adds blastn-short,
          // which NCBI's page replaces with its "Short queries" adjustment.
          {
            flag: '-task',
            label: 'Optimize for',
            kind: 'radio',
            choiceLabels: {
              megablast: 'Highly similar sequences (megablast)',
              'dc-megablast': 'More dissimilar sequences (discontiguous megablast)',
              blastn: 'Somewhat similar sequences (blastn)',
              'blastn-short': 'Short sequences (blastn-short)',
            },
          },
        ],
      },
      {
        title: 'General Parameters',
        placement: 'algorithm',
        fields: [
          maxTargetSeqs,
          evalue,
          {
            ...wordSize,
            list: {
              megablast: ['16', '20', '24', '28', '32', '48', '64'],
              'dc-megablast': ['11', '12'],
              blastn: ['7', '11', '15'],
              'blastn-short': ['7', '11', '15'],
            },
          },
          cullingLimit,
        ],
      },
      {
        title: 'Scoring Parameters',
        placement: 'algorithm',
        fields: [
          { flag: '-reward', label: 'Match', row: 'Match/Mismatch Scores', ...integer },
          { flag: '-penalty', label: 'Mismatch', row: 'Match/Mismatch Scores', ...integer },
          gapOpen,
          gapExtend,
        ],
        pairs: { 'Match/Mismatch Scores': MATCH_SCORES, 'Gap Costs': NUCLEOTIDE_GAP_COSTS },
      },
      {
        title: 'Filters and Masking',
        placement: 'algorithm',
        fields: [
          { flag: '-dust', label: 'Low complexity regions', kind: 'filter', row: 'Filter' },
          softMasking,
          lcaseMasking,
        ],
      },
      {
        title: 'Discontiguous Word Options',
        placement: 'algorithm',
        shownFor: { flag: '-task', values: ['dc-megablast'] },
        // BLAST+ requires the two together; dc-megablast's template is 18, coding.
        together: { '-template_length': '18', '-template_type': 'coding' },
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
          numDescriptions,
          numAlignments,
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
    summary: 'BLASTP mode searches protein subjects using a protein query.',
    searches: 'protein subjects',
    taskNotes: { blastp: 'protein-protein BLAST', 'blastp-fast': 'Quick BLASTP', 'blastp-short': 'optimized for short sequences' },
    query: 'protein',
    subject: 'protein',
    sections: [
      {
        title: 'Task Selection',
        placement: 'program',
        fields: [
          {
            flag: '-task',
            label: 'Algorithm',
            kind: 'radio',
            choiceLabels: {
              blastp: 'blastp (protein-protein BLAST)',
              'blastp-fast': 'Quick BLASTP (blastp-fast)',
              'blastp-short': 'Short sequences (blastp-short)',
            },
          },
        ],
      },
      // LOSAT's BLASTP reads word sizes 2, 3 and 5 to 7 (7, which NCBI's select lacks, through "Other").
      { title: 'General Parameters', placement: 'algorithm', fields: [maxTargetSeqs, evalue, wordSize, cullingLimit] },
      {
        title: 'Scoring Parameters',
        placement: 'algorithm',
        fields: [matrix, gapOpen, gapExtend, { ...compBasedStats, choices: ['0', '2'] }],
        ...PROTEIN_SCORING,
      },
      { title: 'Filters and Masking', placement: 'algorithm', fields: [seg, softMasking, lcaseMasking] },
      { title: OTHER, placement: 'algorithm', fields: [numDescriptions, numAlignments, threshold, windowSize, maxHsps] },
    ],
  },
  {
    id: 'blastx',
    label: 'BLASTX',
    summary: 'BLASTX mode searches protein subjects using a translated nucleotide query.',
    searches: 'protein subjects',
    queryNote: 'translated nucleotide query',
    query: 'nucleotide',
    subject: 'protein',
    sections: [
      { title: 'Task Selection', placement: 'program', fields: [{ flag: '-task', label: 'Algorithm', kind: 'radio' }] },
      { title: 'Genetic code', placement: 'query', fields: [{ flag: '-query_gencode', label: 'Genetic code', kind: 'gencode' }] },
      {
        title: 'General Parameters',
        placement: 'algorithm',
        // LOSAT's BLASTX reads word sizes 3 and 5 only.
        fields: [
          maxTargetSeqs,
          evalue,
          { ...wordSize, list: ['3', '5'] },
          cullingLimit,
        ],
      },
      {
        title: 'Scoring Parameters',
        placement: 'algorithm',
        fields: [{ ...matrix, list: ['BLOSUM62'] }, gapOpen, gapExtend, { ...compBasedStats, choices: ['0', '2'] }],
        ...BLOSUM62_ONLY,
      },
      { title: 'Filters and Masking', placement: 'algorithm', fields: [seg, softMasking, lcaseMasking] },
      {
        title: OTHER,
        placement: 'algorithm',
        fields: [
          numDescriptions,
          numAlignments,
          threshold,
          windowSize,
          maxHsps,
          { flag: '-strand', label: 'Query strand', kind: 'choice' },
          { flag: '-sum_stats', label: 'Sum statistics', kind: 'boolean' },
        ],
      },
    ],
  },
  {
    id: 'tblastn',
    label: 'TBLASTN',
    summary: 'TBLASTN mode searches translated nucleotide subjects using a protein query.',
    searches: 'translated nucleotide subjects',
    queryNote: 'protein query',
    query: 'protein',
    subject: 'nucleotide',
    sections: [
      { title: 'Task Selection', placement: 'program', fields: [{ flag: '-task', label: 'Algorithm', kind: 'radio' }] },
      { title: 'Genetic code', placement: 'subject', fields: [{ flag: '-db_gencode', label: 'Genetic code', kind: 'gencode' }] },
      { title: 'General Parameters', placement: 'algorithm', fields: [maxTargetSeqs, evalue, wordSize, cullingLimit] },
      {
        title: 'Scoring Parameters',
        placement: 'algorithm',
        fields: [matrix, gapOpen, gapExtend, compBasedStats],
        ...PROTEIN_SCORING,
      },
      {
        title: 'Filters and Masking',
        placement: 'algorithm',
        fields: [seg, softMasking, lcaseMasking],
      },
      {
        title: OTHER,
        placement: 'algorithm',
        fields: [
          numDescriptions,
          numAlignments,
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
    summary: 'TBLASTX mode searches translated nucleotide subjects using a translated nucleotide query.',
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
        // NCBI's select (2 and 3); LOSAT's TBLASTX reads word sizes 2 to 4, as NCBI's (4 after "Other").
        fields: [
          maxTargetSeqs,
          evalue,
          { ...wordSize, list: ['2', '3'] },
          cullingLimit,
        ],
      },
      // TBLASTX is ungapped: NCBI's page has the matrix and no gap costs.
      { title: 'Scoring Parameters', placement: 'algorithm', fields: [matrix] },
      { title: 'Filters and Masking', placement: 'algorithm', fields: [seg, softMasking, lcaseMasking] },
      { title: OTHER, placement: 'algorithm', fields: [numDescriptions, numAlignments, threshold, windowSize] },
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
 * The FASTA reader kind of the program's input (abi_v2.md §9, plan TD-8): the engine reads a
 * nucleotide input with NCBI's reader and the nucleotide flags (kind 1), a protein input with
 * the protein flags (kind 2), and the index scan reads it with the same kind.
 */
export function fastaParser(program: ProgramDescriptor, role: InputRole): FastaParserKind {
  return sequenceKind(program, role) === 'nucleotide' ? 1 : 2;
}

/** The unit of a position in an input of this kind. */
export function residueUnit(kind: SequenceKind): string {
  return kind === 'nucleotide' ? 'nt' : 'aa';
}

/**
 * The line beside the "Run LOSAT" button (NCBI's "Search nucleotide sequence using Megablast
 * (Optimize for highly similar sequences)"): what is searched, with which program and task,
 * and where. `task` is the task the search uses (set, or the engine's default), if known.
 */
export function searchSummary(program: ProgramDescriptor, task: string | undefined): string {
  return searchSummaryParts(program, task)
    .map((part) => part.text)
    .join('');
}

/**
 * The line beside the "Run LOSAT" button in parts: what is searched and the program with its
 * task are set off (`strong`), as NCBI sets off "database core_nt" and "Megablast (Optimize for
 * highly similar sequences)".
 */
export function searchSummaryParts(
  program: ProgramDescriptor,
  task: string | undefined,
): ReadonlyArray<{ readonly text: string; readonly strong?: true }> {
  const query = program.queryNote === undefined ? '' : ` (${program.queryNote})`;
  const note = task === undefined ? undefined : program.taskNotes?.[task];
  const taskText = task === undefined ? '' : `, task ${task}${note === undefined ? '' : ` (${note})`}`;
  return [
    { text: 'Search ' },
    { text: program.searches, strong: true },
    { text: ' with ' },
    { text: `${program.label}${query}${taskText}`, strong: true },
    { text: '.' },
  ];
}
