// The reproduction of a run (S15 instructions, item 5; design §12.3): the LOSAT command of each
// output format and the NCBI BLAST+ command to compare it with, both made from the run's argv
// (the RunSnapshot), never from the search form; and the notes that say what a comparison can
// show. The NCBI command is not a mere change of the program's name: every option is one that
// LOSAT Web checked against NCBI BLAST+ 2.17.0, and a run that NCBI cannot run as LOSAT does
// gets no NCBI command but the reason (docs/evidence/losat_web_w6/check_commands.py runs the
// commands of fixed argvs with both programs and compares the outputs byte for byte).
import { toShellCommand } from './argv';
import type { OutputFormat } from './output-format';
import { programById, type InputRole, type ProgramId } from './programs';
import { approvedExceptions } from './verification';

/** The NCBI BLAST+ release that LOSAT's outputs are compared with (the gate records' oracle). */
export const NCBI_BLAST_VERSION = '2.17.0';

const VALUE = true;
const FLAG = false;
const REGIONS = { '-query_loc': VALUE, '-subject_loc': VALUE };
const GENERAL = { '-max_target_seqs': VALUE, '-evalue': VALUE, '-word_size': VALUE };
const PROTEIN_SCORING = { '-matrix': VALUE, '-gapopen': VALUE, '-gapextend': VALUE, '-comp_based_stats': VALUE };

/**
 * The options that the search form writes for each program (domain/programs.ts, and the
 * regions), with whether each takes a value, as `<program> -help` of NCBI BLAST+ 2.17.0 lists
 * them (checked 2026-10-10): NCBI's program has every one of them. BLASTX, which the form cannot
 * search yet, has the options of LOSAT's `blastx` that NCBI's has, without the regions (LOSAT's
 * `blastx` has none).
 */
const NCBI_OPTIONS: Readonly<Record<ProgramId, Readonly<Record<string, boolean>>>> = {
  blastn: {
    '-task': VALUE,
    ...GENERAL,
    '-reward': VALUE,
    '-penalty': VALUE,
    '-gapopen': VALUE,
    '-gapextend': VALUE,
    '-dust': VALUE,
    '-lcase_masking': FLAG,
    '-template_length': VALUE,
    '-template_type': VALUE,
    '-max_hsps': VALUE,
    '-perc_identity': VALUE,
    '-subject_besthit': FLAG,
    ...REGIONS,
  },
  blastp: {
    '-task': VALUE,
    ...GENERAL,
    ...PROTEIN_SCORING,
    '-seg': VALUE,
    '-threshold': VALUE,
    '-window_size': VALUE,
    '-max_hsps': VALUE,
    ...REGIONS,
  },
  blastx: {
    ...GENERAL,
    ...PROTEIN_SCORING,
    '-query_gencode': VALUE,
    '-seg': VALUE,
    '-threshold': VALUE,
    '-window_size': VALUE,
    '-max_hsps': VALUE,
  },
  tblastn: {
    '-task': VALUE,
    '-db_gencode': VALUE,
    ...GENERAL,
    ...PROTEIN_SCORING,
    '-seg': VALUE,
    '-soft_masking': VALUE,
    '-lcase_masking': FLAG,
    '-threshold': VALUE,
    '-window_size': VALUE,
    '-xdrop_gap': VALUE,
    '-xdrop_gap_final': VALUE,
    '-sum_stats': VALUE,
    ...REGIONS,
  },
  tblastx: {
    '-query_gencode': VALUE,
    '-db_gencode': VALUE,
    ...GENERAL,
    '-culling_limit': VALUE,
    '-seg': VALUE,
    '-threshold': VALUE,
    '-window_size': VALUE,
    ...REGIONS,
  },
};

/**
 * The genetic codes that NCBI BLAST+ 2.17.0's command line takes for -query_gencode and
 * -db_gencode ("values between: 1-6, 9-16, 21-31, 33"; blast_args.cpp:997-1056). LOSAT's
 * TBLASTN also takes 32 (PD-TLOSAN-LOCAL-GENCODE-32), which NCBI's command line refuses.
 */
const NCBI_GENETIC_CODES: ReadonlySet<string> = new Set(
  [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33].map(String),
);
const GENETIC_CODE_FLAGS = ['-query_gencode', '-db_gencode'];

/** The LOSAT command that writes one output format of a run (domain/argv.ts). */
export function losatCommand(argv: readonly string[], format: OutputFormat): string {
  return toShellCommand(argv, format);
}

/**
 * The NCBI BLAST+ command with the same program, inputs and options, and `-outfmt`: the LOSAT
 * command without its first word, so that both quote every word alike. No `-num_threads`:
 * the outputs do not depend on it, and NCBI searches a `-subject` with one thread.
 */
export function ncbiCommand(argv: readonly string[], format: OutputFormat): string {
  const losat = toShellCommand(argv, format);
  return losat.slice(losat.indexOf(' ') + 1);
}

export interface NcbiComparison {
  /** Why NCBI BLAST+ cannot run the run's options as LOSAT did; empty when it can. */
  readonly refused: readonly string[];
  /** Approved exceptions that apply (domain/verification.ts): LOSAT's outputs can differ from NCBI's for them. */
  readonly exceptions: readonly string[];
}

/** Whether the NCBI command of a run compares with LOSAT's, and what to know when it does. */
export function ncbiComparison(argv: readonly string[]): NcbiComparison {
  const program = argv[0] as ProgramId;
  const known = NCBI_OPTIONS[program];
  if (known === undefined) return { refused: [`LOSAT Web does not know the program ${argv[0] ?? ''}.`], exceptions: [] };
  const label = programById(program).label;
  const refused: string[] = [];
  let i = 5;
  while (i < argv.length) {
    const flag = argv[i]!;
    const takesValue = known[flag];
    if (takesValue === undefined) {
      refused.push(
        flag.startsWith('-')
          ? `LOSAT Web has not checked ${flag} against ${label} of NCBI BLAST+ ${NCBI_BLAST_VERSION}, so it gives no command to compare with.`
          : `The run's arguments have "${flag}" where an option was expected.`,
      );
      break;
    }
    const value = takesValue ? argv[i + 1] : undefined;
    if (takesValue && value === undefined) {
      refused.push(`The run's arguments end with ${flag} without its value.`);
      break;
    }
    if (value !== undefined && GENETIC_CODE_FLAGS.includes(flag) && !NCBI_GENETIC_CODES.has(value.trim())) {
      const decision = program === 'tblastn' && flag === '-db_gencode' ? ' (PD-TLOSAN-LOCAL-GENCODE-32)' : '';
      refused.push(
        `NCBI BLAST+ ${NCBI_BLAST_VERSION} does not accept ${flag} ${value}: its command line takes the genetic codes ` +
          `1-6, 9-16, 21-31 and 33. LOSAT searched with it${decision}, so there is no NCBI command to compare with.`,
      );
    }
    i += takesValue ? 2 : 1;
  }
  return { refused, exceptions: approvedExceptions(program, argv) };
}

/**
 * The notes of the commands: where the files go, and the threads. The browser does not know the
 * folders of the files, so the commands name them as the run did (design §12.3: no invented
 * paths). `sameBytes` says whether the query and the subject are the same bytes, which matters
 * only when they have the same name. `saved` are the roles whose input is not the chosen file as
 * it is (`inputIsChosenFile`): the commands must run on the file saved from the run for them.
 */
export function commandNotes(argv: readonly string[], sameBytes: boolean, saved: readonly InputRole[] = []): readonly string[] {
  const query = argv[2] ?? '';
  const subject = argv[4] ?? '';
  const notes: string[] = [
    query !== subject
      ? `The commands name the inputs as the run did. Put the files ${query} and ${subject} in one folder and run the commands there; ` +
        'the browser does not know the folders of your files.'
      : sameBytes
        ? `The commands name the input as the run did. Put the file ${query} in a folder and run the commands there; ` +
          'the browser does not know the folders of your files.'
        : `The query and the subject are both named ${query}, but they differ: save them under two names and change the ` +
          'names in the commands to match. The browser does not know the folders of your files.',
  ];
  if (saved.length > 0) {
    const names = [...new Set(saved.map((role) => (role === 'query' ? query : subject)))];
    notes.push(
      `Run the commands with the input FASTA saved from this run (under "The input FASTA of this run" below), not with the file you chose: ` +
        `${names.join(' and ')} there ${names.length === 1 ? 'is' : 'are'} what the run searched.`,
    );
  }
  notes.push('The outputs do not depend on the threads, so the commands do not set -num_threads.');
  return notes;
}

/** Where one part of a run input comes from, as far as the application knows. */
export interface InputPart {
  readonly origin: 'file' | 'paste';
  readonly name: string;
  /** The records of the part's whole source; the run searched all of them but `excluded`. */
  readonly records: number;
  /** Records of the source that the run left out (a session file records them); none when absent. */
  readonly excluded?: number;
}

/**
 * Whether a role's run input is the file that was chosen, byte for byte: one whole file (or
 * nothing known about the parts), not pasted text, a joined input, or a file with records left out.
 */
export function inputIsChosenFile(parts: ReadonlyArray<InputPart | undefined>): boolean {
  if (parts.length === 0) return true;
  const part = parts[0];
  return parts.length === 1 && part !== undefined && part.origin === 'file' && (part.excluded ?? 0) === 0;
}

/**
 * How the bytes that the engine searched relate to what was chosen (S15 item 5): the same bytes
 * as one file or the pasted text, or files joined, or a file with records left out. A part is
 * undefined when it is a selection of a source's records (some left out).
 */
export function inputRelation(role: InputRole, name: string, records: number, parts: ReadonlyArray<InputPart | undefined>): string {
  const count = `${records} ${records === 1 ? 'record' : 'records'}`;
  const whole = (part: InputPart | undefined) => part !== undefined && (part.excluded ?? 0) === 0;
  // One record left out reads in the singular; a part that is unknown (undefined) may have left out several.
  const only = (left: ReadonlyArray<InputPart | undefined>) => !left.some((part) => part === undefined) && left.reduce((sum, part) => sum + (part?.excluded ?? 0), 0) === 1;
  if (parts.length === 0) return `${name} has the ${count} that the run searched.`;
  if (parts.length > 1) {
    const names = parts.map((part) => part?.name).filter((part): part is string => part !== undefined);
    const leftOut = parts.some((part) => !whole(part)) ? `, without the ${only(parts) ? 'record' : 'records'} left out` : '';
    return (
      `${name} joins the ${parts.length} ${role} inputs${names.length > 0 ? ` (${names.join(', ')})` : ''} in the order chosen${leftOut}: ` +
      `${count}. It is no single file you chose.`
    );
  }
  const part = parts[0];
  if (part === undefined) return `${name} has the ${count} that the run searched; the records left out of the chosen ${role} are not in it.`;
  if (!whole(part)) {
    return `${name} has the ${count} that the run searched from the file ${part.name}; the ${only(parts) ? 'record' : 'records'} left out of it ${only(parts) ? 'is' : 'are'} not in it.`;
  }
  if (part.origin === 'paste') return `${name} is the pasted ${role} text, as the run searched it (${count}).`;
  return `${name} has the same bytes as the file ${part.name} (${count}).`;
}
