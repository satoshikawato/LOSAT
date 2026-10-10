// Extraction (instruction items 2-4, design §11.4, plan §5.7): which intervals of which original
// records to read for the HSPs of the candidate tray, and the FASTA that the application writes
// from them. Two outputs, kept apart:
//
// - Sequences: the original residues of a record, read from the source File by the Data worker
//   (`readResidues`), in the record's own direction (the only orientation; aligning a sequence
//   to the hit's direction is later work), case and letters as the file has them. A region is
//   the HSP's interval (`hit`), that interval with left and right flanks on the record's axis
//   (`flanked`, left toward position 1), or the whole record (`whole`). Several HSPs on the same
//   record of the same run and role give one sequence each (`separate`, overlaps kept) or one
//   interval from their smallest to their largest position (`spanning`); nothing else is merged,
//   and records of different runs are never joined, even when they come from the same file.
//   `whole` gives one sequence per record, listing every HSP of the tray on it. Ends cut a region
//   and the header keeps the requested interval beside the actual one. The header (an
//   application format, not NCBI's):
//
//     >{id}:{from}-{to} run={N} {role}_record={K} length={L} unit={nt|aa} hsps={Q.R,...} hit_strand={plus|minus|unknown|mixed}[ requested={a}-{b}]
//
//   {from}-{to} is the extracted interval (1-based, both ends included, from <= to, whatever the
//   strand); N the run's number; role `query` or `subject`; K the record's position in the run's
//   input plus 1 (q_idx + 1 or s_idx + 1, which tells records with the same ID apart); L the
//   record's length in its unit; Q.R each HSP (query record q_idx + 1, HSP rank + 1); hit_strand
//   the strand of the HSPs on this record (`unknown` for a BLASTN HSP of one letter, whose record
//   does not give it; `mixed` where they differ); `requested` only where the record's ends cut
//   the region (a may be 0 or less, b more than L). A record without an ID is named as the
//   engine names it in outfmt 6, `Query_K` or `Subject_K` (`recordLabel`).
//
// - Gapped alignments: for each HSP, the query and subject rows of its HSP record
//   (`query_aligned`, `subject_aligned`) exactly as the engine wrote them (gaps `-`, case as
//   written), never re-derived from the original residues, with the record's coordinates (start
//   > end kept: the direction of the search, never reversed):
//
//     >{qid}:{q_start}-{q_end} run={N} query_record={K} hsp={Q.R} aligned[ frame={f}]
//     >{sid}:{s_start}-{s_end} run={N} subject_record={K} hsp={Q.R} aligned[ frame={f}]
//
//   `frame` only for a translated sequence. An HSP whose record has no aligned rows is reported.
//
// Residues are wrapped at 60 per line. Output is bytes (headers UTF-8), so that a writer can
// stream it record by record (S15). Nothing here reads a file or computes a BLAST value.
import {
  clipToRecord,
  commonStrand,
  flanked,
  intervalLength,
  intervalText,
  isClipped,
  spanning,
  spanOn,
  translates,
  unitOf,
  withinRecord,
  type ClippedInterval,
  type Flanks,
  type HspCoordinates,
  type HspSpan,
  type Interval,
  type Strand,
  type Unit,
} from './coordinates';
import { programById, type InputRole, type ProgramId } from './programs';

/** Residues per line of the FASTA that extraction writes (an application choice). */
export const FASTA_LINE_WIDTH = 60;

/** The HSP record fields that extraction uses (docs/web/abi_v2.md §8; `HspRecord` has them all). */
export interface ExtractedHsp extends HspCoordinates {
  readonly q_idx: number;
  readonly s_idx: number;
  readonly rank: number;
}

/** The HSP record fields of a gapped alignment. */
export interface AlignedHsp extends ExtractedHsp {
  readonly query_aligned: string | null;
  readonly subject_aligned: string | null;
}

/** The run of an HSP, as extraction names it. */
export interface RunRef {
  readonly runId: string;
  /** The run's number in the working session (`RunSnapshot.number`). */
  readonly number: number;
  readonly program: ProgramId;
}

/** One HSP of the tray on one of its records (query or subject). */
export interface ExtractionTarget {
  readonly runId: string;
  readonly runNumber: number;
  readonly role: InputRole;
  /** 0-based position of the record in the run's input of `role`: the HSP's `q_idx` or `s_idx`. */
  readonly position: number;
  readonly recordId: string;
  /** The record's length in `unit` (the record table's `length`). */
  readonly recordLength: number;
  readonly unit: Unit;
  /** The HSP's interval and strand on the record. */
  readonly span: HspSpan;
  /** The HSP's label `Q.R` (`hspLabel`). */
  readonly hsp: string;
}

export type ExtractionRegion =
  | { readonly kind: 'hit' }
  | { readonly kind: 'flanked'; readonly flanks: Flanks }
  | { readonly kind: 'whole' };

/** How several HSPs on one record of one run and role become sequences. */
export type ExtractionJoin = 'separate' | 'spanning';

export interface ExtractionOptions {
  readonly region: ExtractionRegion;
  readonly join: ExtractionJoin;
}

/** One sequence to extract: an interval of one record of a run's input. */
export interface ExtractionPiece {
  readonly runId: string;
  readonly runNumber: number;
  readonly role: InputRole;
  readonly position: number;
  readonly recordId: string;
  readonly recordLength: number;
  readonly unit: Unit;
  /** The interval asked for and the one the record has (`actual` is what is read). */
  readonly interval: ClippedInterval;
  /** The labels of the HSPs that the piece is for, in tray order. */
  readonly hsps: readonly string[];
  readonly strand: Strand | 'mixed';
}

export interface ExtractionPlan {
  /** In tray order: a joined piece takes the place of its first HSP. */
  readonly pieces: readonly ExtractionPiece[];
  /** What extraction did not decide, and why (English, one sentence each). */
  readonly notes: readonly string[];
}

/** One `readResidues` call: a record of a run's input and the intervals of its pieces. */
export interface ReadRequest {
  readonly runId: string;
  readonly role: InputRole;
  readonly position: number;
  /** The actual interval of each piece, in the order of `pieces`. */
  readonly intervals: readonly Interval[];
  /** Indices into the plan's pieces. */
  readonly pieces: readonly number[];
}

/** "3.2": the query record q_idx + 1 and the HSP rank + 1 (the HSP's place in the results). */
export const hspLabel = (qIdx: number, rank: number): string => `${qIdx + 1}.${rank + 1}`;

/**
 * The name of a record in extraction headers: its ID, or for a record without one (no title)
 * the name that the engine gives it in outfmt 6 (its local ID `Query_K` or `Subject_K`, K the
 * record's position in the input plus 1; LOSAT/src/blastinput/fasta_reader `shown_id`).
 */
export function recordLabel(id: string, role: InputRole, position: number): string {
  return id !== '' ? id : `${role === 'query' ? 'Query' : 'Subject'}_${position + 1}`;
}

/** The target of an HSP on its query or subject record (the record table's ID and length). */
export function extractionTarget(
  run: RunRef,
  hsp: ExtractedHsp,
  role: InputRole,
  record: { readonly id: string; readonly length: number },
): ExtractionTarget {
  const program = programById(run.program);
  const kind = role === 'query' ? program.query : program.subject;
  return {
    runId: run.runId,
    runNumber: run.number,
    role,
    position: role === 'query' ? hsp.q_idx : hsp.s_idx,
    recordId: record.id,
    recordLength: record.length,
    unit: unitOf(kind),
    span: spanOn(hsp, role, kind),
    hsp: hspLabel(hsp.q_idx, hsp.rank),
  };
}

const groupKey = (target: ExtractionTarget): string => JSON.stringify([target.runId, target.role, target.position]);

/** The pieces to extract for the targets, in tray order (see the file comment). */
export function planExtraction(targets: readonly ExtractionTarget[], options: ExtractionOptions): ExtractionPlan {
  const groups = new Map<string, ExtractionTarget[]>();
  for (const target of targets) {
    if (!withinRecord(target.span, target.recordLength)) {
      throw new RangeError(`HSP ${target.hsp} at ${intervalText(target.span)} is not within ${target.role} record ${target.position + 1} (${target.recordLength} ${target.unit})`);
    }
    const key = groupKey(target);
    const group = groups.get(key);
    if (group === undefined) {
      groups.set(key, [target]);
    } else {
      const first = group[0]!;
      if (first.recordId !== target.recordId || first.recordLength !== target.recordLength || first.unit !== target.unit) {
        throw new Error(`the HSPs on ${target.role} record ${target.position + 1} of run ${target.runNumber} disagree about the record`);
      }
      group.push(target);
    }
  }
  const { region, join } = options;
  const cut = (wanted: Interval, length: number): ClippedInterval =>
    region.kind === 'flanked' ? flanked(wanted, region.flanks, length) : clipToRecord(wanted, length);
  const piece = (members: readonly ExtractionTarget[], interval: ClippedInterval): ExtractionPiece => {
    const first = members[0]!;
    return {
      runId: first.runId,
      runNumber: first.runNumber,
      role: first.role,
      position: first.position,
      recordId: first.recordId,
      recordLength: first.recordLength,
      unit: first.unit,
      interval,
      hsps: [...new Set(members.map((member) => member.hsp))],
      strand: commonStrand(members.map((member) => member.span)),
    };
  };
  const pieces: ExtractionPiece[] = [];
  if (region.kind !== 'whole' && join === 'separate') {
    for (const target of targets) pieces.push(piece([target], cut(target.span, target.recordLength)));
  } else {
    for (const members of groups.values()) {
      const { recordLength } = members[0]!;
      const interval =
        region.kind === 'whole'
          ? clipToRecord({ from: 1, to: recordLength }, recordLength)
          : cut(spanning(members.map((member) => member.span)), recordLength);
      pieces.push(piece(members, interval));
    }
  }
  const notes = [...new Set(targets.filter((target) => target.span.strand === 'unknown').map(unknownStrandNote))];
  return { pieces, notes };
}

function unknownStrandNote(target: ExtractionTarget): string {
  return (
    `HSP ${target.hsp} of run ${target.runNumber} covers one letter of ${target.role} record ${target.position + 1}, ` +
    'so its HSP record does not give its strand and the extraction does not decide it (hit_strand=unknown); ' +
    'the Strand= line of its outfmt 0 section shows it.'
  );
}

/** The `readResidues` calls of a plan: one per record, with the intervals of its pieces. */
export function readRequests(plan: ExtractionPlan): readonly ReadRequest[] {
  const requests = new Map<string, { runId: string; role: InputRole; position: number; intervals: Interval[]; pieces: number[] }>();
  plan.pieces.forEach((piece, i) => {
    const key = JSON.stringify([piece.runId, piece.role, piece.position]);
    let request = requests.get(key);
    if (request === undefined) {
      request = { runId: piece.runId, role: piece.role, position: piece.position, intervals: [], pieces: [] };
      requests.set(key, request);
    }
    request.intervals.push(piece.interval.actual);
    request.pieces.push(i);
  });
  return [...requests.values()];
}

/** The header line of a sequence (without `>` and the line end); see the file comment. */
export function sequenceHeader(piece: ExtractionPiece): string {
  const { actual, requested } = piece.interval;
  return (
    `${recordLabel(piece.recordId, piece.role, piece.position)}:${intervalText(actual)} run=${piece.runNumber} ` +
    `${piece.role}_record=${piece.position + 1} length=${piece.recordLength} unit=${piece.unit} ` +
    `hsps=${piece.hsps.join(',')} hit_strand=${piece.strand}` +
    (isClipped(piece.interval) ? ` requested=${requested.from}-${requested.to}` : '')
  );
}

/** The FASTA record of a piece from the residues that `readResidues` read for its actual interval. */
export function sequenceFasta(piece: ExtractionPiece, residues: Uint8Array): Uint8Array {
  const expected = intervalLength(piece.interval.actual);
  if (residues.length !== expected) {
    throw new RangeError(`${residues.length} residues were given for ${intervalText(piece.interval.actual)} (${expected} ${piece.unit})`);
  }
  return fastaRecord(sequenceHeader(piece), residues);
}

/** The header lines (without `>`) of an HSP's two aligned rows; see the file comment. */
export function alignmentHeaders(
  run: RunRef,
  hsp: AlignedHsp,
  ids: { readonly query: string; readonly subject: string },
): { readonly query: string; readonly subject: string } {
  const program = programById(run.program);
  const label = hspLabel(hsp.q_idx, hsp.rank);
  const frame = (value: number | null, translated: boolean) => (translated && value !== null && value !== 0 ? ` frame=${value}` : '');
  const line = (role: InputRole, id: string, position: number, start: number, end: number, value: number | null, translated: boolean) =>
    `${recordLabel(id, role, position)}:${start}-${end} run=${run.number} ${role}_record=${position + 1} hsp=${label} aligned${frame(value, translated)}`;
  return {
    query: line('query', ids.query, hsp.q_idx, hsp.q_start, hsp.q_end, hsp.query_frame, translates(run.program, program.query)),
    subject: line('subject', ids.subject, hsp.s_idx, hsp.s_start, hsp.s_end, hsp.subject_frame, translates(run.program, program.subject)),
  };
}

/**
 * The two FASTA records of an HSP's gapped alignment, from its record's aligned rows as written;
 * undefined when the record has no aligned rows (nothing is invented).
 */
export function alignmentFasta(
  run: RunRef,
  hsp: AlignedHsp,
  ids: { readonly query: string; readonly subject: string },
): Uint8Array | undefined {
  if (hsp.query_aligned === null || hsp.subject_aligned === null) return undefined;
  const headers = alignmentHeaders(run, hsp, ids);
  const query = fastaRecord(headers.query, encoder.encode(hsp.query_aligned));
  const subject = fastaRecord(headers.subject, encoder.encode(hsp.subject_aligned));
  const out = new Uint8Array(query.length + subject.length);
  out.set(query);
  out.set(subject, query.length);
  return out;
}

/** Why an HSP has no gapped alignment to write (one sentence, for the extraction's notes). */
export function missingAlignmentNote(run: RunRef, hsp: AlignedHsp): string {
  return `HSP ${hspLabel(hsp.q_idx, hsp.rank)} of run ${run.number} has no aligned sequences in its HSP record, so it has no gapped alignment to write.`;
}

const encoder = new TextEncoder();
const LF = 0x0a;

/** `>header`, then the letters in lines of `width`, every line ended with LF. */
export function fastaRecord(header: string, letters: Uint8Array, width = FASTA_LINE_WIDTH): Uint8Array {
  if (!Number.isSafeInteger(width) || width < 1) throw new RangeError(`a line of ${width} letters`);
  const head = encoder.encode(`>${header}\n`);
  const lines = Math.ceil(letters.length / width);
  const out = new Uint8Array(head.length + letters.length + lines);
  out.set(head);
  let at = head.length;
  for (let start = 0; start < letters.length; start += width) {
    const line = letters.subarray(start, Math.min(letters.length, start + width));
    out.set(line, at);
    at += line.length;
    out[at++] = LF;
  }
  return out;
}
