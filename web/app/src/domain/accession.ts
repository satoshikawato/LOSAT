// Accessions typed in a paste box (line WF-1, docs/web/fetch_design.md §4): NCBI's box takes
// "accession number(s) … or FASTA sequence(s)"; LOSAT fetches an accession's FASTA only when the
// user asks (ports/fetch.ts) and writes it into the box.
import type { SequenceKind } from './programs';

export type FetchSourceId = 'ncbi' | 'ebi' | 'cncb';

/**
 * The sources, in the order of the choice ("Fetch from Fastest / NCBI / EBI / CNCB"; Fastest is
 * the default, FetchChoice). EBI is ENA for nucleotides and UniProt for proteins. DDBJ sends no
 * CORS header (fetch_design.md §2).
 */
export const FETCH_SOURCES: ReadonlyArray<{ readonly id: FetchSourceId; readonly label: string }> = [
  { id: 'ncbi', label: 'NCBI' },
  { id: 'ebi', label: 'EBI' },
  { id: 'cncb', label: 'CNCB' },
];

/**
 * Where "Fetch" asks: a source, or 'fastest' (the default, line WF-3), the first source that
 * answers this browser among those that serve the accessions (application/source-race.ts).
 */
export type FetchChoice = FetchSourceId | 'fastest';

/**
 * The fixed public record that the race asks each source for, never the user's accession
 * (fetch_design.md §10): nucleotide MH011443.1 (TP53 exon 5, 192 bytes) at all three; protein
 * the E. coli thr operon leader peptide, which NCBI and CNCB name NP_414542.1 and UniProt P0AD86
 * (CNCB has no UniProt accessions).
 */
export const RACE_RECORDS: Readonly<Record<SequenceKind, Readonly<Record<FetchSourceId, string>>>> = {
  nucleotide: { ncbi: 'MH011443.1', ebi: 'MH011443.1', cncb: 'MH011443.1' },
  protein: { ncbi: 'NP_414542.1', ebi: 'P0AD86', cncb: 'NP_414542.1' },
};

/** Bytes of FASTA that one fetch writes into a paste box at most. */
export const FETCH_BYTE_LIMIT = 10_000_000;

/** Accessions that one fetch takes at most (NCBI answers 3 requests a second without a key). */
export const FETCH_ACCESSION_LIMIT = 100;

/**
 * INSDC (A12345, AB123456, AAA12345, WGS AAAA01000001), RefSeq (NR_074540, NZ_CP006932,
 * WP_000001), CNCB GenBase (C_AA004835) and UniProt (P12345, Q9Y261, A0A023GPI8) accessions,
 * each with an optional version. A sequence line has no digits, so it never matches.
 */
const ACCESSION =
  /^(?:[A-Z]{1,6}\d{5,}|[A-Z]{1,2}_[A-Z]{0,6}\d{5,}|[OPQ]\d[A-Z\d]{3}\d|[A-NR-Z]\d(?:[A-Z][A-Z\d]{2}\d){1,2})(?:\.\d+)?$/;
/** RefSeq accessions (two letters and "_"), which ENA does not serve. */
const REFSEQ = /^[A-Z]{2}_/;
/** CNCB GenBase's own accessions (C_AA004835), which only GenBase serves (NCBI and ENA refuse them, 2026-11-03). */
const GENBASE = /^C_/;
/** UniProt accessions (uniprot.org/help/accession_numbers), the only ones UniProt's REST serves. */
const UNIPROT = /^(?:[OPQ]\d[A-Z\d]{3}\d|[A-NR-Z]\d(?:[A-Z][A-Z\d]{2}\d){1,2})(?:\.\d+)?$/;

/**
 * The accessions of a text that has only accessions, one per line (blank lines aside), in
 * upper case; undefined for any other text (FASTA, a bare sequence, an empty box).
 */
export function accessionLines(text: string): readonly string[] | undefined {
  const lines = text
    .split(/\r?\n/)
    .map((line) => line.trim().toUpperCase())
    .filter((line) => line !== '');
  return lines.length > 0 && lines.every((line) => ACCESSION.test(line)) ? lines : undefined;
}

/**
 * The accessions of a text that has only accessions but not one per line (separated by white
 * space or commas, as lists are often copied), as typed; undefined for any other text, and for a
 * text that accessionLines takes. NCBI's reader cannot read such a line (WA2-OP-04): the source's
 * card offers to put each accession on its own line, where Fetch takes them.
 */
export function accessionsOnALine(text: string): readonly string[] | undefined {
  if (accessionLines(text) !== undefined) return undefined;
  const words = text.split(/[\s,]+/).filter((word) => word !== '');
  return words.length > 0 && words.every((word) => ACCESSION.test(word.toUpperCase())) ? words : undefined;
}

/**
 * The part of the engine's reason that names a first line NCBI BLAST+ would fetch as a Seq-id
 * (adapter, scan/ncbi.rs). The adapter's answer has no kind or code, only this text, so the
 * app tells that reason by this part of it (WA2-OP-04); the CLI's message is not the adapter's.
 */
const SEQ_ID_REASON = 'may be a sequence identifier that NCBI BLAST+ fetches through a data loader';

/** True when a source's error is the reason for a first line that NCBI BLAST+ would fetch as a Seq-id. */
export function isSeqIdReason(error: string | undefined): boolean {
  return error !== undefined && error.includes(SEQ_ID_REASON);
}

/**
 * Whether the card of a pasted text that failed with `error` says how Fetch takes accessions
 * (WA2-OP-04): the reason is the Seq-id one and the box's first line is an accession, which is
 * what Fetch takes. A GI number, a PDB id (`4HHB`, `1ABC_A`) and a tagged id (`sp|P04637|...`)
 * are Seq-ids to the engine but not to Fetch, and the user has put them as they are.
 */
export function pointsToFetch(error: string | undefined, paste: string): boolean {
  if (!isSeqIdReason(error)) return false;
  const first = paste.split(/\r?\n/).find((line) => line.trim() !== '');
  return first !== undefined && ACCESSION.test(first.trim().toUpperCase());
}

/**
 * Why a source cannot serve an accession of a kind, without asking it, or undefined (measured on
 * 2026-11-03): ENA has no RefSeq records, UniProt has only its own accessions, only GenBase has
 * GenBase's, and GenBase has no UniProt accessions (an empty answer).
 */
export function unserved(source: FetchSourceId, kind: SequenceKind, accession: string): string | undefined {
  switch (source) {
    case 'ncbi':
      return GENBASE.test(accession) ? 'NCBI has no GenBase records' : undefined;
    case 'ebi':
      if (kind === 'nucleotide' && REFSEQ.test(accession)) return 'EBI (ENA) has no RefSeq records';
      if (kind === 'nucleotide' && GENBASE.test(accession)) return 'EBI (ENA) has no GenBase records';
      if (kind === 'protein' && !UNIPROT.test(accession)) return 'EBI (UniProt) has UniProt accessions only';
      return undefined;
    case 'cncb':
      return kind === 'protein' && UNIPROT.test(accession) ? 'CNCB has no UniProt records' : undefined;
  }
}

/** The sources that serve every one of the accessions (unserved), in the order of FETCH_SOURCES; it may be none. */
export function servingSources(kind: SequenceKind, accessions: readonly string[]): readonly FetchSourceId[] {
  return FETCH_SOURCES.map((source) => source.id).filter((source) =>
    accessions.every((accession) => unserved(source, kind, accession) === undefined),
  );
}

/** A part of a box that no one source serves: the accessions that the same sources serve (fetchGroups). */
export interface FetchGroup {
  readonly sources: readonly FetchSourceId[];
  readonly accessions: readonly string[];
  /** The lines of the box (among its accessions, 0-based) that the accessions are, ascending. */
  readonly lines: readonly number[];
}

/**
 * The box's accessions by the set of sources that serve each, in the order of each group's first
 * accession in the box (line WF-5, fetch_design.md §11). Each group has a "Fetch" of its own,
 * which sends its accessions to one of its sources; fetchAccessions sends a box that one source
 * serves whole.
 */
export function fetchGroups(kind: SequenceKind, accessions: readonly string[]): readonly FetchGroup[] {
  const groups = new Map<string, { sources: readonly FetchSourceId[]; accessions: string[]; lines: number[] }>();
  for (const [line, accession] of accessions.entries()) {
    const sources = servingSources(kind, [accession]);
    const key = sources.join(',');
    const group = groups.get(key) ?? { sources, accessions: [], lines: [] };
    group.accessions.push(accession);
    group.lines.push(line);
    groups.set(key, group);
  }
  return [...groups.values()];
}

export function sourceLabel(source: FetchSourceId): string {
  return FETCH_SOURCES.find((candidate) => candidate.id === source)!.label;
}

/** "NCBI", "NCBI or CNCB", "NCBI, EBI or CNCB" (`word` "and" for a list of all of them). */
export function sourceList(sources: readonly FetchSourceId[], word: 'or' | 'and' = 'or'): string {
  const labels = sources.map(sourceLabel);
  return labels.length < 2 ? (labels[0] ?? '') : `${labels.slice(0, -1).join(', ')} ${word} ${labels[labels.length - 1]}`;
}
