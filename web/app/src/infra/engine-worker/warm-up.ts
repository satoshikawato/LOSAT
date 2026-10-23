// The warm-up of a prepared instance. The browser compiles a Wasm function when it is first
// called (V8's lazy compilation, the default in Chromium and Node), on the thread that calls
// it, so the first search of a page compiles the engine code it calls while it runs. A
// preparation that knows the argv of the coming search runs it once, on one thread, on the
// small built-in sequences below and discards the outputs. The compiled code belongs to the
// compiled module, which the thread workers share. One thread starts no thread worker: a
// failure stays on the thread that calls the engine (`abi.stoppedBy`) and cannot reach a
// search queued behind the warm-up. In Chromium 149 (line W s12) the first search of a 7 kb
// TBLASTX input on one thread took 68-273 ms without it and 29-74 ms with it (later searches
// 17-55 ms). Longer searches gain less or nothing: what remains is the time until the
// optimizing compiler's code of their hot functions is ready, which the small sequences do
// not reach.
import { programById, type ProgramId, type SequenceKind } from '../../domain/programs';
import { REGION_FLAG } from '../../domain/region';
import { ROLE_QUERY, ROLE_SUBJECT, type ReactorAbi } from '../reactor/abi';

/** Options whose values are positions in the user's inputs: the warm-up leaves them out. */
const INPUT_OPTIONS: ReadonlySet<string> = new Set(Object.values(REGION_FLAG));

// The standard genetic code, codons in TCAG order.
const BASES = 'TCAG';
const CODE = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
const SENSE_CODONS = [...CODE].flatMap((amino, i) =>
  amino === '*' ? [] : [BASES[i >> 4]! + BASES[(i >> 2) & 3]! + BASES[i & 3]!],
);

/** `count` sense codons from a linear congruential generator started at `seed`. */
function codons(seed: number, count: number): string {
  let state = seed >>> 0;
  let dna = '';
  for (let i = 0; i < count; i++) {
    state = (Math.imul(state, 1_103_515_245) + 12_345) >>> 0;
    dna += SENSE_CODONS[(state >>> 16) % SENSE_CODONS.length];
  }
  return dna;
}

function translate(dna: string): string {
  let protein = '';
  for (let i = 0; i + 3 <= dna.length; i += 3) {
    protein += CODE[16 * BASES.indexOf(dna[i]!) + 4 * BASES.indexOf(dna[i + 1]!) + BASES.indexOf(dna[i + 2]!)];
  }
  return protein;
}

// Sequences without a stop codon in frame 1 (936 and 903 nt): the query is three blocks with a low-complexity
// stretch for the masking, the subject the same blocks in the opposite order, the middle one with
// one codon more and one codon changed. Each block is a hit of its own (none lies within
// another's ranges), and the middle one has a gap.
const A = codons(1, 120);
const B = codons(3, 80);
const C = codons(2, 100);
const QUERY = A + 'CAG'.repeat(12) + B + C;
const SUBJECT = C + B.slice(0, 120) + 'GGC' + B.slice(120, 180) + 'TGG' + B.slice(183) + A;
const encoder = new TextEncoder();
const fasta = (dna: string, kind: SequenceKind) => encoder.encode(`>warm-up\n${kind === 'protein' ? translate(dna) : dna}\n`);

/** The built-in query and subject of `program`, nucleotide or protein as the program reads them. */
export function warmUpInputs(program: ProgramId): { readonly query: Uint8Array; readonly subject: Uint8Array } {
  const descriptor = programById(program);
  return { query: fasta(QUERY, descriptor.query), subject: fasta(SUBJECT, descriptor.subject) };
}

/**
 * Runs the search of `argv` (program first, without -num_threads) once on one thread on the
 * built-in sequences and discards its outputs. Returns whether the run completed. An argv
 * that the engine refuses or a run that fails is no error here (the search reports its own);
 * a stop of the instance stays in `abi.stoppedBy` for the caller.
 */
export function warmUp(abi: ReactorAbi, argv: readonly string[]): boolean {
  const runArgv: string[] = [];
  for (let i = 0; i < argv.length; i++) {
    if (INPUT_OPTIONS.has(argv[i]!)) i++;
    else runArgv.push(argv[i]!);
  }
  runArgv.push('-num_threads', '1');
  const handles: number[] = [];
  try {
    const program = argv[0] as ProgramId;
    const inputs = warmUpInputs(program);
    if (abi.validate(runArgv) !== undefined) return false;
    handles.push(abi.register(program, ROLE_QUERY, inputs.query).handle);
    handles.push(abi.register(program, ROLE_SUBJECT, inputs.subject).handle);
    abi.run(runArgv, handles[0]!, handles[1]!, () => undefined);
    return true;
  } catch {
    return false;
  } finally {
    for (const handle of handles) {
      if (abi.stoppedBy !== undefined) break;
      try {
        abi.release(handle);
      } catch {
        // A handle that cannot be released only keeps its copy of the built-in sequence.
      }
    }
  }
}
