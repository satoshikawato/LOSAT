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
   * For `text` fields: the values of NCBI's select, offered in a select (one list, or a list per
   * task); its "Other" takes any other value in a text box, which the engine's parser reads.
   */
  readonly list?: readonly string[] | Readonly<Record<string, readonly string[]>>;
  /** NCBI's sentence beside the control ("Select the maximum number of aligned sequences to display"). */
  readonly hint?: string;
  /**
   * The name of a row that shows this field with its neighbours of the same row (NCBI's
   * "Match/Mismatch Scores" and "Gap Costs", which are one select on NCBI's page).
   */
  readonly row?: string;
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
}

/**
 * A default that the engine resolves itself and `describe` does not report (the option has no
 * clap default): by task, or from other options. The values are those of the CLI's help
 * (`losat <program> -help`; describe.json carries the same text, and a unit test checks each
 * value against it). This table is a temporary source (line WR, 2026-11-02): the engine's route
 * is `describe` reporting the defaults it resolves (docs/web/ncbi_ui_mapping.md §2).
 */
export interface ResolvedDefault {
  /** The value, or the value of each task. */
  readonly value: string | Readonly<Record<string, string>>;
  /**
   * Options whose value decides this one ("11, or the matrix's best value"): the value holds
   * while each of them keeps its own default.
   */
  readonly unless?: readonly string[];
  /**
   * BLAST+ tells the option given from not given (the hit list options, abi_v2.md §7): a value
   * set is written even when it equals the default.
   */
  readonly given?: true;
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
  /** The engine's FASTA reader for this program, which the index scan follows (plan TD-8). */
  readonly fastaParser: FastaParserKind;
  /**
   * False when the engine takes no -query_loc / -subject_loc for the program: the form shows
   * no subrange and writes none (DW-9).
   */
  readonly subrange?: false;
  /**
   * The search form. The order of the sections and fields is the order of the argv: the
   * task first, as on a command line, then the genetic codes and the algorithm parameters
   * in the order of the page.
   */
  readonly sections: readonly ParameterSection[];
  /** The defaults that the engine resolves and `describe` does not report, by flag. */
  readonly defaults?: Readonly<Record<string, ResolvedDefault>>;
}

const integer = { kind: 'text', inputMode: 'numeric' } as const;
const decimal = { kind: 'text', inputMode: 'decimal' } as const;
/** NCBI's matrix select, in its order; "Other" takes the other names that the engine reads. */
const MATRICES = ['PAM30', 'PAM70', 'BLOSUM90', 'BLOSUM80', 'BLOSUM62', 'BLOSUM50', 'BLOSUM45', 'PAM250'];
/** NCBI's word size selects of the protein programs. */
const PROTEIN_WORD_SIZES = ['2', '3', '5', '6'];

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
// Scoring Parameters: NCBI's "Gap Costs" select is one row of the two costs here.
const matrix: ParameterField = { flag: '-matrix', label: 'Matrix', kind: 'text', list: MATRICES };
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

/**
 * The hit list options: BLAST+ tells each given from not given, and -max_target_seqs excludes
 * the other two (abi_v2.md §7), so a default is shown only while the others are not set.
 */
const HIT_LIST_DEFAULTS: Readonly<Record<string, ResolvedDefault>> = {
  '-max_target_seqs': { value: '500', given: true, unless: ['-num_descriptions', '-num_alignments'] },
  '-num_descriptions': { value: '500', given: true, unless: ['-max_target_seqs'] },
  '-num_alignments': { value: '250', given: true, unless: ['-max_target_seqs'] },
};
/** The values of BLASTP and TBLASTN that the matrix decides ("or the matrix's best value"). */
const GAP_COST_DEFAULTS: Readonly<Record<string, ResolvedDefault>> = {
  '-gapopen': { value: '11', unless: ['-matrix'] },
  '-gapextend': { value: '1', unless: ['-matrix'] },
  '-window_size': { value: '40', unless: ['-matrix'] },
};

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
    // Kind 0 is the `bio::io::fasta` reader; kind 1 is BLASTX's NCBI-style reader (ABI v2 §4).
    fastaParser: 0,
    sections: [
      {
        title: 'Task Selection',
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
      },
      {
        title: 'Filters and Masking',
        placement: 'algorithm',
        fields: [
          { flag: '-dust', label: 'Low complexity regions filter (DUST)', kind: 'text', suggestions: ['no', 'yes'] },
          softMasking,
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
          numDescriptions,
          numAlignments,
          maxHsps,
          { flag: '-perc_identity', label: 'Percent identity', ...decimal },
          { flag: '-subject_besthit', label: 'Subject best hit', kind: 'flag' },
        ],
      },
    ],
    defaults: {
      ...HIT_LIST_DEFAULTS,
      '-evalue': { value: { megablast: '10', 'dc-megablast': '10', blastn: '10', 'blastn-short': '1000' } },
      '-word_size': { value: { megablast: '28', 'dc-megablast': '11', blastn: '11', 'blastn-short': '7' } },
      '-reward': { value: { megablast: '1', 'dc-megablast': '2', blastn: '2', 'blastn-short': '1' } },
      '-penalty': { value: { megablast: '-2', 'dc-megablast': '-3', blastn: '-3', 'blastn-short': '-3' } },
      '-gapopen': { value: { megablast: '0', 'dc-megablast': '5', blastn: '5', 'blastn-short': '5' } },
      '-gapextend': { value: { megablast: '0', 'dc-megablast': '2', blastn: '2', 'blastn-short': '2' } },
      '-dust': { value: { megablast: '20 64 1', 'dc-megablast': '20 64 1', blastn: '20 64 1', 'blastn-short': 'no' } },
      '-soft_masking': { value: 'true' },
    },
  },
  {
    id: 'blastp',
    label: 'BLASTP',
    summary: 'BLASTP mode searches protein subjects using a protein query.',
    searches: 'protein subjects',
    taskNotes: { blastp: 'protein-protein BLAST', 'blastp-fast': 'Quick BLASTP' },
    query: 'protein',
    subject: 'protein',
    fastaParser: 0,
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
            },
          },
        ],
      },
      { title: 'General Parameters', placement: 'algorithm', fields: [maxTargetSeqs, evalue, wordSize] },
      { title: 'Scoring Parameters', placement: 'algorithm', fields: [matrix, gapOpen, gapExtend, compBasedStats] },
      { title: 'Filters and Masking', placement: 'algorithm', fields: [seg, softMasking, lcaseMasking] },
      { title: OTHER, placement: 'algorithm', fields: [numDescriptions, numAlignments, threshold, windowSize, maxHsps] },
    ],
    defaults: {
      ...HIT_LIST_DEFAULTS,
      ...GAP_COST_DEFAULTS,
      '-evalue': { value: '10' },
      '-word_size': { value: { blastp: '3', 'blastp-fast': '5' } },
      '-matrix': { value: 'BLOSUM62' },
      '-seg': { value: 'no' },
      '-soft_masking': { value: 'false' },
      '-threshold': { value: { blastp: '11', 'blastp-fast': '20' }, unless: ['-matrix', '-word_size'] },
    },
  },
  {
    id: 'blastx',
    label: 'BLASTX',
    summary: 'BLASTX mode searches protein subjects using a translated nucleotide query.',
    searches: 'protein subjects',
    queryNote: 'translated nucleotide query',
    query: 'nucleotide',
    subject: 'protein',
    fastaParser: 0,
    // The engine refuses -query_loc and -subject_loc for BLASTX (line SX, DW-11).
    subrange: false,
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
          { flag: '-culling_limit', label: 'Max matches in a query range', ...integer },
        ],
      },
      // LOSAT's BLASTX scores with BLOSUM62 and gap costs 11/1 only.
      { title: 'Scoring Parameters', placement: 'algorithm', fields: [{ ...matrix, list: ['BLOSUM62'] }, gapOpen, gapExtend, compBasedStats] },
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
    defaults: {
      ...HIT_LIST_DEFAULTS,
      '-word_size': { value: { blastx: '3', 'blastx-fast': '5' } },
      '-culling_limit': { value: '0' },
      '-matrix': { value: 'BLOSUM62' },
      '-gapopen': { value: '11' },
      '-gapextend': { value: '1' },
      '-threshold': { value: { blastx: '12', 'blastx-fast': '20' }, unless: ['-word_size'] },
      '-window_size': { value: '40' },
      '-sum_stats': { value: 'true' },
    },
  },
  {
    id: 'tblastn',
    label: 'TBLASTN',
    summary: 'TBLASTN mode searches translated nucleotide subjects using a protein query.',
    searches: 'translated nucleotide subjects',
    queryNote: 'protein query',
    query: 'protein',
    subject: 'nucleotide',
    fastaParser: 0,
    sections: [
      { title: 'Task Selection', placement: 'program', fields: [{ flag: '-task', label: 'Algorithm', kind: 'radio' }] },
      { title: 'Genetic code', placement: 'subject', fields: [{ flag: '-db_gencode', label: 'Genetic code', kind: 'gencode' }] },
      { title: 'General Parameters', placement: 'algorithm', fields: [maxTargetSeqs, evalue, wordSize] },
      { title: 'Scoring Parameters', placement: 'algorithm', fields: [matrix, gapOpen, gapExtend, compBasedStats] },
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
    defaults: {
      ...HIT_LIST_DEFAULTS,
      ...GAP_COST_DEFAULTS,
      '-evalue': { value: '10' },
      '-word_size': { value: { tblastn: '3', 'tblastn-fast': '5' } },
      '-matrix': { value: 'BLOSUM62' },
      '-seg': { value: '12 2.2 2.5' },
      '-soft_masking': { value: 'false' },
      '-threshold': { value: { tblastn: '13', 'tblastn-fast': '20' }, unless: ['-matrix', '-word_size'] },
      '-xdrop_gap': { value: '15' },
      '-xdrop_gap_final': { value: '25' },
      '-sum_stats': { value: 'true' },
    },
  },
  {
    id: 'tblastx',
    label: 'TBLASTX',
    summary: 'TBLASTX mode searches translated nucleotide subjects using a translated nucleotide query.',
    searches: 'translated nucleotide subjects',
    queryNote: 'translated nucleotide query',
    query: 'nucleotide',
    subject: 'nucleotide',
    fastaParser: 0,
    sections: [
      { title: 'Genetic code', placement: 'query', fields: [{ flag: '-query_gencode', label: 'Genetic code', kind: 'gencode' }] },
      { title: 'Genetic code', placement: 'subject', fields: [{ flag: '-db_gencode', label: 'Genetic code', kind: 'gencode' }] },
      {
        title: 'General Parameters',
        placement: 'algorithm',
        // NCBI's "Max matches in a query range" (HSP_RANGE_MAX) is BLAST+'s -culling_limit.
        fields: [
          maxTargetSeqs,
          evalue,
          { ...wordSize, list: ['2', '3'] },
          { flag: '-culling_limit', label: 'Max matches in a query range', ...integer },
        ],
      },
      // TBLASTX is ungapped: NCBI's page has the matrix and no gap costs.
      { title: 'Scoring Parameters', placement: 'algorithm', fields: [matrix] },
      { title: 'Filters and Masking', placement: 'algorithm', fields: [seg, softMasking, lcaseMasking] },
      { title: OTHER, placement: 'algorithm', fields: [numDescriptions, numAlignments, threshold, windowSize] },
    ],
    defaults: {
      ...HIT_LIST_DEFAULTS,
      '-word_size': { value: '3' },
      '-culling_limit': { value: '0' },
      '-matrix': { value: 'BLOSUM62' },
      '-soft_masking': { value: 'false' },
      '-threshold': { value: '13', unless: ['-matrix'] },
      '-window_size': { value: '40', unless: ['-matrix'] },
    },
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
