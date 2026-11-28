// The options of the command lines whose value is a list that the command line splits (line WW, WW4b2 R3 and WW4d): one
// place for the forms (`listNameRefusal`, `seqdbNameList`, `seqdbNames`) and the Supported mark (`toolOptionKey` in
// verification.ts). No imports: the build's generator of the verification table (build/verification.ts) reads it too.

/**
 * makeblastdb's options whose value the command line splits at commas (`NStr::Split` of the value with the spaces at
 * its ends cut: NCBI makeblastdb.cpp `x_ProcessMaskData`; LOSAT writedb `mask_data/process.rs` `split_commas`).
 */
export const COMMA_LIST_FLAGS: ReadonlySet<string> = new Set(['-mask_data', '-mask_id', '-mask_desc', '-gi_mask_name']);

/** The number of entries of a comma list as the command line splits it: none for an empty value, else one more than its commas. */
export function commaListLength(value: string): number {
  const cut = value.trim();
  return cut === '' ? 0 : cut.split(',').length;
}

/**
 * The options whose value the command line reads as a list of BLAST database (or file) names, as NCBI's SeqDB splits
 * one (`SeqDB_SplitQuoted`): a space separates two names, a double quote starts or ends a name that holds spaces, so no
 * name can hold a double quote. Measured against NCBI BLAST+ 2.17.0 for every such option that the Web writes (line WW,
 * WW4d probe, 2026-11-28): the search programs' `-db`, BLASTN's `-filtering_db`, makeblastdb's `-in` (its FASTA files
 * and its `-input_type blastdb` databases), and the `-in` of WindowMasker, dustmasker and segmasker with `-infmt blastdb`.
 * The same flag means other things in other commands (`isSeqdbList`): their `-in` is one file's name with another
 * `-infmt`.
 */
export const SEQDB_LIST_FLAGS: ReadonlySet<string> = new Set(['-db', '-filtering_db', '-in']);

/** An option of a command: the flag, the command's tool where the flag depends on it, and the tool's `-infmt` value. */
export interface ListOption {
  /**
   * The tool of the command line (`makeblastdb`, `windowmasker`, …) for a flag whose meaning depends on it (`-in`, the
   * comma lists); omitted for `-db` and `-filtering_db`, which every search program that has them reads alike.
   */
  readonly command?: string;
  readonly flag: string;
  /** The `-infmt` value of WindowMasker, dustmasker or segmasker, which decides what its `-in` is (`fasta` when not given). */
  readonly inputFormat?: string;
}

/**
 * The tools whose `-in` is a SeqDB name list with `-infmt blastdb` (one file's name with another `-infmt`): WindowMasker,
 * and the Masking tab's dustmasker and segmasker (measured against NCBI BLAST+ 2.17.0, WW4d R1 probe, 2026-11-28).
 */
const SEQDB_INFMT_COMMANDS: ReadonlySet<string> = new Set(['windowmasker', 'dustmasker', 'segmasker']);

/** Whether the option's value is a SeqDB name list (`SEQDB_LIST_FLAGS`), for its command. */
export function isSeqdbList(option: ListOption): boolean {
  if (!SEQDB_LIST_FLAGS.has(option.flag)) return false;
  if (option.flag !== '-in') return true;
  if (option.command === 'makeblastdb') return true;
  return SEQDB_INFMT_COMMANDS.has(option.command ?? '') && option.inputFormat?.trim() === 'blastdb';
}

/** What a list option's splitter takes out of a name (`listNameRefusal`): the mark, and why. */
interface ListSplit {
  readonly mark: string;
  readonly reason: string;
}

const SEQDB_LIST_SPLIT: ListSplit = { mark: '"', reason: 'a double quote there starts or ends a name' };
const COMMA_LIST_SPLIT: ListSplit = { mark: ',', reason: 'a comma separates its files' };

function listSplit(option: ListOption): ListSplit | undefined {
  if (isSeqdbList(option)) return SEQDB_LIST_SPLIT;
  if (option.command === 'makeblastdb' && COMMA_LIST_FLAGS.has(option.flag)) return COMMA_LIST_SPLIT;
  return undefined;
}

/**
 * Why a name cannot be written in the list of names that an option takes, as the command line splits the list, or
 * undefined when it can (and for an option that takes one name). The one rule for every list option: the SeqDB lists
 * (`isSeqdbList`; no name can hold a double quote) and makeblastdb's options of `COMMA_LIST_FLAGS` (no entry can hold a
 * comma). `what` says what the name is of: a file opened, or a database.
 */
export function listNameRefusal(option: ListOption, name: string, what: 'file' | 'database' = 'file'): string | undefined {
  const split = listSplit(option);
  if (split === undefined || !name.includes(split.mark)) return undefined;
  return `A ${what} named "${name}" cannot be in ${option.flag}: ${split.reason}.`;
}

/**
 * The one writer of a SeqDB name list (`isSeqdbList`), as the command line reads it back (`seqdbNames`): the names
 * separated by spaces, a name that holds a space between double quotes (`SeqDB_CombineAndQuote`); an empty name is left
 * out (the reader drops it), a name of spaces is quoted as any name with a space. A name with a double quote cannot be
 * written: it throws the refusal (`listNameRefusal`), which the place that takes the name shows before it comes here.
 */
export function seqdbNameList(flag: string, names: readonly string[], what: 'file' | 'database' = 'database'): string {
  return names
    .filter((name) => name !== '')
    .map((name) => {
      if (name.includes(SEQDB_LIST_SPLIT.mark)) throw new Error(`A ${what} named "${name}" cannot be in ${flag}: ${SEQDB_LIST_SPLIT.reason}.`);
      return name.includes(' ') ? `"${name}"` : name;
    })
    .join(' ');
}

/**
 * The names of a SeqDB name list, as the command line reads them (`SeqDB_SplitQuoted`): a space separates two names, a
 * double quote starts or ends a name (the quotes are not part of it; an unclosed one runs to the end), and empty names
 * are left out. The one reader for the places that need a list's names back from an argv.
 */
export function seqdbNames(value: string): readonly string[] {
  const names: string[] = [];
  let current = '';
  let quoted = false;
  for (const char of value) {
    if (char === '"') {
      if (current !== '') names.push(current);
      current = '';
      quoted = !quoted;
    } else if (char === ' ' && !quoted) {
      if (current !== '') names.push(current);
      current = '';
    } else {
      current += char;
    }
  }
  if (current !== '') names.push(current);
  return Object.freeze(names);
}
