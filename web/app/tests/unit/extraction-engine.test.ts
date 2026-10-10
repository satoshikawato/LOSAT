// Extraction against real searches (S14 instruction item 6): the serial reactor runs in Node.
// Generated FASTA files are indexed through DataService with the engine's scan (reader kinds 1
// and 2), run inputs are built from their revisions, the searches run, their outputs are
// committed, and for every HSP record (read back with readHspRecords) the residues that
// readResidues reads for the hit interval, flanked regions, spanning intervals and whole records
// (planned by domain/extraction.ts) are the generated sequence's letters, case and U kept. The
// coordinates are cross-checked against the engine's own letters: the aligned rows without gaps
// are the extracted hit interval (case-insensitive; reverse complement for a BLASTN minus hit;
// for TBLASTN and TBLASTX the extracted nucleotides translated here in the HSP's frame with the
// standard code - test-only code - where the aligned letter is not masked, X or *).
// BLASTX cannot run until SX (its coordinates are covered by extraction.test.ts), and a first
// record without a defline cannot be indexed in this tree (checkRecordTable), so neither is here.
// It needs LOSAT_WEB_REACTORS.
import { beforeAll, describe, expect, it } from 'vitest';
import { findReactors } from '../../build/reactors';
import { interval } from '../../src/domain/coordinates';
import { recordMismatch, type DatasetRevision, type FastaParserKind, type RecordKey } from '../../src/domain/dataset';
import {
  extractionTarget,
  planExtraction,
  readRequests,
  recordLabel,
  sequenceFasta,
  type ExtractionOptions,
  type RunRef,
} from '../../src/domain/extraction';
import { splitOutfmt6Row } from '../../src/domain/outfmt6';
import type { InputRole, ProgramId } from '../../src/domain/programs';
import { sha256Hex } from '../../src/infra/browser/platform';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { ROLE_QUERY, ROLE_SUBJECT, type ReactorAbi } from '../../src/infra/reactor/abi';
import { ReactorInputChecker } from '../../src/infra/reactor/checker';
import { instantiateSerial } from '../../src/infra/reactor/instance';
import { ReactorScanner } from '../../src/infra/reactor/scanner';
import { RunOutputWriter } from '../../src/infra/run-output/writer';
import type { HspRecord } from '../../src/ports/engine';
import type { OutputStream } from '../../src/ports/run-output';
import { FastaWriter, residues, seeded, type Lines, type ReaderKind } from './support/fasta-writer';

const latin1 = new TextDecoder('latin1');
const AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY';

// --- test-only sequence helpers ---------------------------------------------------------------------

const COMPLEMENT: Readonly<Record<string, string>> = { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N', R: 'Y', Y: 'R', K: 'M', M: 'K', S: 'S', W: 'W', B: 'V', V: 'B', D: 'H', H: 'D' };
/** Upper case, U as T: how the engine reads a nucleotide. */
const normalized = (letters: string) => letters.toUpperCase().replaceAll('U', 'T');
const reverseComplement = (letters: string) => [...normalized(letters)].reverse().map((c) => COMPLEMENT[c] ?? 'N').join('');

/** The standard genetic code (NCBI table 1), codons in TCAG order. */
const CODE = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
const BASE: Readonly<Record<string, number>> = { T: 0, C: 1, A: 2, G: 3 };
function translate(nucleotides: string): string {
  let out = '';
  for (let i = 0; i + 3 <= nucleotides.length; i += 3) {
    const [a, b, c] = [BASE[nucleotides[i]!], BASE[nucleotides[i + 1]!], BASE[nucleotides[i + 2]!]];
    out += a === undefined || b === undefined || c === undefined ? 'X' : CODE[a * 16 + b * 4 + c]!;
  }
  return out;
}
/** A random coding sequence of a protein (no stop codons for its residues). */
function backTranslate(random: () => number, protein: string): string {
  const codons = new Map<string, string[]>();
  for (let i = 0; i < 64; i++) {
    const codon = ['T', 'C', 'A', 'G'][i >> 4]! + ['T', 'C', 'A', 'G'][(i >> 2) & 3]! + ['T', 'C', 'A', 'G'][i & 3]!;
    codons.set(CODE[i]!, [...(codons.get(CODE[i]!) ?? []), codon]);
  }
  return [...protein].map((aa) => {
    const choices = codons.get(aa.toUpperCase())!;
    return choices[Math.floor(random() * choices.length)]!;
  }).join('');
}

/** `letters` with `insert` written over it from 0-based `at`. */
const plant = (letters: string, at: number, insert: string) => letters.slice(0, at) + insert + letters.slice(at + insert.length);

// --- the service and the searches --------------------------------------------------------------------

interface Source {
  readonly name: string;
  readonly kind: ReaderKind;
  readonly records: ReadonlyArray<{ readonly title: string; readonly letters: string; readonly lines: Lines }>;
  /** Text before the first record (comments, blank lines). */
  readonly lead?: string;
}

interface Indexed {
  readonly revision: DatasetRevision;
  readonly letters: readonly string[];
}

interface Run {
  readonly ref: RunRef;
  readonly revisions: Readonly<Record<InputRole, readonly string[]>>;
  /** The generated letters of each record of the run input, by position. */
  readonly letters: Readonly<Record<InputRole, readonly string[]>>;
  readonly records: Readonly<Record<InputRole, readonly RecordKey[]>>;
  readonly hsps: readonly HspRecord[];
}

const reactors = findReactors();

describe.skipIf(reactors === undefined)('extraction reads the residues that the engine searched (real searches)', () => {
  let abi: ReactorAbi;
  let data: DataService;
  let runs = 0;
  const summary = { runs: 0, hsps: 0, pieces: 0, clippedLeft: 0, clippedRight: 0, alignedLetters: 0 };

  beforeAll(async () => {
    abi = (await instantiateSerial(new WebAssembly.Module(reactors!.serial.bytes as BufferSource))).abi;
    let token = 0;
    data = new DataService({
      store: new MemoryBlockStore(),
      scanner: new ReactorScanner(async () => abi),
      checker: new ReactorInputChecker(async () => abi),
      digest: sha256Hex,
      newToken: () => `token-${++token}`,
      cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
      // Small reads, so that every record is read in several slices.
      readChunkBytes: 4096,
    });
  }, 60_000);

  async function index(source: Source): Promise<Indexed> {
    const writer = new FastaWriter(seeded(source.name.length), source.kind);
    if (source.lead !== undefined) writer.raw(source.lead);
    for (const record of source.records) writer.record(record.title, record.letters, record.lines);
    const ref = await data.addSource(new File([writer.toString()], source.name));
    const revision = await data.indexSource(ref.sourceId, source.kind as FastaParserKind);
    expect(revision.records.map((record) => record.length), `${source.name}: the scan's lengths`).toEqual(source.records.map((record) => record.letters.length));
    return { revision, letters: source.records.map((record) => record.letters) };
  }

  /** The letters of the included records of revisions, in run-input order. */
  function included(indexed: readonly Indexed[], revisions: readonly DatasetRevision[]): string[] {
    return revisions.flatMap((revision, i) => indexed[i]!.letters.filter((_, k) => !revision.excluded.includes(k)));
  }

  async function search(
    program: ProgramId,
    words: readonly string[],
    query: { readonly indexed: readonly Indexed[]; readonly revisions?: readonly DatasetRevision[] },
    subject: { readonly indexed: readonly Indexed[]; readonly revisions?: readonly DatasetRevision[] },
  ): Promise<Run> {
    const runId = `run-${++runs}`;
    const ref: RunRef = { runId, number: runs, program };
    const roles = { query, subject };
    const revisions = {} as Record<InputRole, readonly string[]>;
    const letters = {} as Record<InputRole, readonly string[]>;
    const records = {} as Record<InputRole, readonly RecordKey[]>;
    const handles: number[] = [];
    const port = await data.openRun(runId);
    const writer = new RunOutputWriter(port);
    try {
      for (const role of ['subject', 'query'] as const) {
        const chosen = roles[role].revisions ?? roles[role].indexed.map((each) => each.revision);
        revisions[role] = chosen.map((revision) => revision.revisionId);
        letters[role] = included([...roles[role].indexed], [...chosen]);
        const input = await data.buildRunInput(revisions[role]);
        records[role] = input.records;
        const registered = abi.register(program, role === 'query' ? ROLE_QUERY : ROLE_SUBJECT, input.bytes);
        handles.push(registered.handle);
        expect(recordMismatch(input.records, registered.records), `${runId} ${role}: the engine read the record table's records`).toBeUndefined();
      }
      abi.run([program, '-query', 'query.fa', '-subject', 'subject.fa', ...words, '-num_threads', '1'], handles[1]!, handles[0]!, (stream, bytes) =>
        writer.write(stream as OutputStream, bytes),
      );
      writer.end();
    } finally {
      for (const handle of handles) abi.release(handle);
    }
    const ref6 = await data.commitRun(runId);
    const all = await data.readHits(runId);
    const hsps = await data.readHspRecords(runId, all.map((hsp) => hsp.index).reverse());
    expect([...hsps].reverse(), `${runId}: readHspRecords gives the records that readHits gives`).toEqual(all);
    expect(ref6.hitCount).toBe(all.length);
    summary.runs++;
    summary.hsps += all.length;
    return { ref, revisions, letters, records, hsps: all };
  }

  /** Every region and join for every HSP on both of its records reads the generated letters. */
  async function checkExtraction(run: Run): Promise<void> {
    const options: readonly ExtractionOptions[] = [
      { region: { kind: 'hit' }, join: 'separate' },
      { region: { kind: 'flanked', flanks: { left: 37, right: 250 } }, join: 'separate' },
      { region: { kind: 'flanked', flanks: { left: 400, right: 3 } }, join: 'spanning' },
      { region: { kind: 'hit' }, join: 'spanning' },
      { region: { kind: 'whole' }, join: 'separate' },
    ];
    for (const role of ['query', 'subject'] as const) {
      const targets = run.hsps.map((hsp) => {
        const position = role === 'query' ? hsp.q_idx : hsp.s_idx;
        return extractionTarget(run.ref, hsp, role, run.records[role][position]!);
      });
      for (const option of options) {
        const plan = planExtraction(targets, option);
        for (const request of readRequests(plan)) {
          const got = await data.readResidues(run.revisions[role], request.position, request.intervals);
          const generated = run.letters[role][request.position]!;
          expect(got.origin.length).toBe(generated.length);
          request.pieces.forEach((p, i) => {
            const piece = plan.pieces[p]!;
            const { from, to } = piece.interval.actual;
            expect(latin1.decode(got.residues[i]), `${run.ref.runId} ${role} ${request.position + 1} ${from}-${to}`).toBe(generated.slice(from - 1, to));
            sequenceFasta(piece, got.residues[i]!);
            summary.pieces++;
            if (piece.interval.clippedLeft) summary.clippedLeft++;
            if (piece.interval.clippedRight) summary.clippedRight++;
          });
        }
      }
    }
  }

  /** The aligned rows without gaps are the extracted hit intervals, in the HSP's direction or frame. */
  async function checkAligned(run: Run): Promise<void> {
    for (const hsp of run.hsps) {
      for (const role of ['query', 'subject'] as const) {
        const [start, end, frame, aligned] =
          role === 'query' ? [hsp.q_start, hsp.q_end, hsp.query_frame, hsp.query_aligned] : [hsp.s_start, hsp.s_end, hsp.subject_frame, hsp.subject_aligned];
        const position = role === 'query' ? hsp.q_idx : hsp.s_idx;
        const where = `${run.ref.runId} HSP ${hsp.index} ${role}`;
        const read = latin1.decode((await data.readResidues(run.revisions[role], position, [interval(start, end)])).residues[0]);
        const row = aligned!.replaceAll('-', '');
        const program = run.ref.program;
        const translated = (program === 'tblastn' && role === 'subject') || program === 'tblastx';
        if (program === 'blastn') {
          expect(row.toUpperCase(), where).toBe(start <= end ? normalized(read) : reverseComplement(read));
          summary.alignedLetters += row.length;
          continue;
        }
        const expected = translated ? translate(frame! < 0 ? reverseComplement(read) : normalized(read)) : read.toUpperCase();
        if (translated) expect(read.length % 3, `${where}: whole codons`).toBe(0);
        expect(row.length, where).toBe(expected.length);
        for (let i = 0; i < row.length; i++) {
          const letter = row[i]!;
          if (letter !== letter.toUpperCase() || letter === 'X' || letter === '*' || expected[i] === 'X') continue;
          expect(letter, `${where} letter ${i + 1}`).toBe(expected[i]);
          summary.alignedLetters++;
        }
      }
    }
  }

  it(
    'BLASTN: plus and minus strands, flanks cut at both ends, lower case, U and comment lines, an untitled subject, a record over 65536 residues (checkpoints) and a uniform one with CR LF',
    async () => {
      const random = seeded(101);
      const query = residues(random, 700, 'ACGT', 0.1);
      const near = plant(residues(random, 900, 'ACGT'), 2, query.slice(10, 260).replaceAll('T', 'U'));
      const far = plant(residues(random, 1500, 'acgtACGT'), 1500 - 205, reverseComplement(query.slice(400, 600)).toLowerCase());
      const untitled = plant(residues(random, 600, 'ACGT'), 100, query.slice(300, 420));
      const ragged = plant(residues(random, 140_000, 'ACGT', 0.2), 100_000, reverseComplement(query.slice(100, 400)));
      const uniform = plant(residues(random, 70_000, 'ACGT'), 69_500, query.slice(450, 690));
      const q = await index({
        name: 'query.fa',
        kind: 1,
        lead: '; generated query\n\n',
        records: [{ title: 'q1 the query', letters: query, lines: { kind: 'ragged', minWidth: 50, maxWidth: 70, eol: '\n', comments: true } }],
      });
      const s = await index({
        name: 'subjects.fa',
        kind: 1,
        records: [
          { title: 'near plus hit at the start', letters: near, lines: { kind: 'uniform', width: 60, eol: '\n' } },
          { title: 'far minus hit at the end', letters: far, lines: { kind: 'ragged', minWidth: 30, maxWidth: 90, eol: '\n', comments: true } },
          { title: '', letters: untitled, lines: { kind: 'uniform', width: 80, eol: '  \n' } },
          { title: 'ragged', letters: ragged, lines: { kind: 'ragged', minWidth: 55, maxWidth: 65, eol: '\n', comments: true } },
        ],
      });
      const crlf = await index({ name: 'crlf.fa', kind: 1, records: [{ title: 'uniform crlf', letters: uniform, lines: { kind: 'uniform', width: 60, eol: '\r\n' } }] });
      expect(s.revision.records.map((record) => record.line_layout.kind)).toEqual(['uniform', 'checkpoints', 'uniform', 'checkpoints']);
      expect(s.revision.records[2]!.line_layout).toEqual({ kind: 'uniform', width: 80, eol: 3 });
      expect(crlf.revision.records[0]!.line_layout).toEqual({ kind: 'uniform', width: 60, eol: 2 });
      const run = await search('blastn', [], { indexed: [q] }, { indexed: [s, crlf] });
      const strands = (sIdx: number) => run.hsps.filter((hsp) => hsp.s_idx === sIdx).map((hsp) => (hsp.s_start <= hsp.s_end ? 'plus' : 'minus'));
      expect(strands(0)).toContain('plus');
      expect(strands(1)).toContain('minus');
      expect(strands(2)).toContain('plus');
      expect(strands(3)).toContain('minus');
      expect(strands(4)).toContain('plus');
      expect(run.hsps.some((hsp) => hsp.s_idx === 3 && Math.min(hsp.s_start, hsp.s_end) > 65_536)).toBe(true);
      // The untitled record is named in outfmt 6 as recordLabel names it.
      const untitledHsp = run.hsps.find((hsp) => hsp.s_idx === 2)!;
      const row = splitOutfmt6Row(new TextDecoder().decode(await data.readOutputRange(run.ref.runId, 6, untitledHsp.out6![0], untitledHsp.out6![1])));
      expect(row.sseqid).toBe(recordLabel('', 'subject', 2));
      const clipped = { ...summary };
      await checkExtraction(run);
      expect(summary.clippedLeft - clipped.clippedLeft, 'flanks cut at position 1').toBeGreaterThan(0);
      expect(summary.clippedRight - clipped.clippedRight, "flanks cut at the record's end").toBeGreaterThan(0);
      await checkAligned(run);
    },
    60_000,
  );

  it('BLASTN with -subject_loc and with -query_loc: coordinates stay those of the whole record', async () => {
    const random = seeded(202);
    const query = residues(random, 800, 'ACGT');
    const subject = plant(plant(residues(random, 5000, 'ACGT'), 500, query.slice(0, 200)), 3000, reverseComplement(query.slice(500, 750)));
    const q = await index({ name: 'loc-query.fa', kind: 1, records: [{ title: 'lq', letters: query, lines: { kind: 'uniform', width: 70, eol: '\n' } }] });
    const s = await index({ name: 'loc-subject.fa', kind: 1, records: [{ title: 'ls', letters: subject, lines: { kind: 'uniform', width: 60, eol: '\n' } }] });
    const subjectLoc = await search('blastn', ['-subject_loc', '2001-4500'], { indexed: [q] }, { indexed: [s] });
    expect(subjectLoc.hsps.length).toBeGreaterThan(0);
    for (const hsp of subjectLoc.hsps) expect(Math.min(hsp.s_start, hsp.s_end)).toBeGreaterThan(2000);
    const queryLoc = await search('blastn', ['-query_loc', '451-800'], { indexed: [q] }, { indexed: [s] });
    expect(queryLoc.hsps.length).toBeGreaterThan(0);
    for (const hsp of queryLoc.hsps) expect(hsp.q_start).toBeGreaterThan(450);
    for (const run of [subjectLoc, queryLoc]) {
      await checkExtraction(run);
      await checkAligned(run);
    }
  }, 60_000);

  it('two subject records with the same ID: the hit on the second, then the first left out (positions shift)', async () => {
    const random = seeded(303);
    const query = residues(random, 400, 'ACGT');
    const second = plant(residues(random, 1200, 'ACGT'), 600, query.slice(50, 350));
    const q = await index({ name: 'dup-query.fa', kind: 1, records: [{ title: 'dq', letters: query, lines: { kind: 'uniform', width: 60, eol: '\n' } }] });
    const s = await index({
      name: 'dup-subject.fa',
      kind: 1,
      records: [
        { title: 'dup first', letters: residues(random, 1000, 'ACGT'), lines: { kind: 'uniform', width: 60, eol: '\n' } },
        { title: 'dup second', letters: second, lines: { kind: 'uniform', width: 60, eol: '\n' } },
      ],
    });
    const both = await search('blastn', [], { indexed: [q] }, { indexed: [s] });
    expect(new Set(both.hsps.map((hsp) => hsp.s_idx))).toEqual(new Set([1]));
    await checkExtraction(both);
    await checkAligned(both);
    const withoutFirst = await data.reviseDataset(s.revision.revisionId, [0]);
    const shifted = await search('blastn', [], { indexed: [q] }, { indexed: [s], revisions: [withoutFirst] });
    expect(new Set(shifted.hsps.map((hsp) => hsp.s_idx))).toEqual(new Set([0]));
    const origin = (await data.readResidues(shifted.revisions.subject, 0, [])).origin;
    expect([origin.id, origin.recordIndex]).toEqual(['dup', 1]);
    await checkExtraction(shifted);
    await checkAligned(shifted);
  }, 60_000);

  it('BLASTP, TBLASTN with a minus frame, and TBLASTX', async () => {
    const random = seeded(404);
    const protein = residues(random, 220, AMINO_ACIDS, 0.15);
    const subjectProtein = plant(residues(random, 500, AMINO_ACIDS), 300, protein.slice(20, 180).toUpperCase());
    const pq = await index({ name: 'p-query.fa', kind: 2, records: [{ title: 'pq protein', letters: protein, lines: { kind: 'ragged', minWidth: 40, maxWidth: 70, eol: '\r\n', comments: true } }] });
    const ps = await index({ name: 'p-subject.fa', kind: 2, records: [{ title: 'ps', letters: subjectProtein, lines: { kind: 'uniform', width: 60, eol: '\n' } }] });
    const blastp = await search('blastp', [], { indexed: [pq] }, { indexed: [ps] });
    expect(blastp.hsps.length).toBeGreaterThan(0);

    const coding = backTranslate(random, protein.slice(30, 200));
    const genome = plant(plant(residues(random, 3000, 'ACGT', 0.1), 1700, reverseComplement(coding)), 200, backTranslate(random, protein.slice(0, 60)));
    const ns = await index({ name: 'genome.fa', kind: 1, records: [{ title: 'genome', letters: genome, lines: { kind: 'ragged', minWidth: 60, maxWidth: 61, eol: '\n', comments: true } }] });
    const tblastn = await search('tblastn', [], { indexed: [pq] }, { indexed: [ns] });
    expect(tblastn.hsps.some((hsp) => (hsp.subject_frame ?? 0) < 0)).toBe(true);
    expect(tblastn.hsps.some((hsp) => (hsp.subject_frame ?? 0) > 0)).toBe(true);

    const nucleotideQuery = residues(random, 900, 'ACGT');
    const nucleotideSubject = plant(plant(residues(random, 2500, 'ACGT'), 400, nucleotideQuery.slice(100, 520)), 1800, reverseComplement(nucleotideQuery.slice(600, 840)));
    const nq = await index({ name: 'tx-query.fa', kind: 1, records: [{ title: 'txq', letters: nucleotideQuery, lines: { kind: 'uniform', width: 60, eol: '\n' } }] });
    const nsx = await index({ name: 'tx-subject.fa', kind: 1, records: [{ title: 'txs', letters: nucleotideSubject, lines: { kind: 'uniform', width: 70, eol: '\r\n' } }] });
    const tblastx = await search('tblastx', [], { indexed: [nq] }, { indexed: [nsx] });
    expect(tblastx.hsps.some((hsp) => (hsp.query_frame ?? 0) * (hsp.subject_frame ?? 0) < 0)).toBe(true);

    for (const run of [blastp, tblastn, tblastx]) {
      await checkExtraction(run);
      await checkAligned(run);
    }
    console.log(
      `extraction against real searches: ${summary.runs} runs, ${summary.hsps} HSPs, ${summary.pieces} pieces read, ` +
        `${summary.clippedLeft} cut at position 1 and ${summary.clippedRight} at the end, ${summary.alignedLetters} aligned letters compared`,
    );
  }, 60_000);
});
