// The expectations of V-BR (plan §6.1, §6.2, TD-5): the native CLI of the same commit, run
// with the argv that the application ran (one format at a time, one thread), as V-ABI does
// (web/adapter/tests/v_abi.js), and the NCBI-frozen bytes of LOSAT/tests/outfmt0_manifest.tsv.
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { join, resolve } from 'node:path';
import type { OutputFormat } from '../../../src/domain/output-format';
import type { SearchCase } from '../harness/engine';
import { REPOSITORY } from './harness-server';

/** The native LOSAT binary of the same commit as the reactors (LOSAT_WEB_NATIVE). */
export const NATIVE = process.env['LOSAT_WEB_NATIVE'] || undefined;
export const NO_NATIVE_REASON = 'V-BR needs LOSAT_WEB_NATIVE: the native LOSAT built from the same commit as the reactors';

const FORMATS: readonly OutputFormat[] = [0, 6, 7];

export interface NativeExpectation {
  readonly sha256: Record<OutputFormat, string>;
  /** SHA-256 of the CLI's standard error (the warnings; stream 3 of the adapter). */
  readonly stderrSha256: string;
}

function sha256(bytes: Uint8Array | string): string {
  return createHash('sha256').update(bytes).digest('hex');
}

const cache = new Map<string, NativeExpectation>();

/** Runs the native CLI for every format with `argv` (program first) from `cwd` (relative to the repository, or absolute). */
export function nativeExpectation(argv: readonly string[], cwd: string): NativeExpectation {
  const key = JSON.stringify([argv, cwd]);
  const cached = cache.get(key);
  if (cached !== undefined) return cached;
  if (NATIVE === undefined) throw new Error(NO_NATIVE_REASON);
  const formats: Partial<Record<OutputFormat, string>> = {};
  let stderr: string | undefined;
  for (const format of FORMATS) {
    const args = [...argv, '-outfmt', String(format), '-num_threads', '1'];
    const result = spawnSync(NATIVE, args, { cwd: resolve(REPOSITORY, cwd), maxBuffer: 1 << 30 });
    if (result.status !== 0) throw new Error(`native ${args.join(' ')} failed (${result.status}): ${result.stderr}`);
    formats[format] = sha256(result.stdout);
    const warnings = sha256(result.stderr);
    if (stderr !== undefined && stderr !== warnings) throw new Error(`native ${argv.join(' ')}: the warnings differ between formats`);
    stderr = warnings;
  }
  const expectation = { sha256: formats as Record<OutputFormat, string>, stderrSha256: stderr! };
  cache.set(key, expectation);
  return expectation;
}

/** A V-BR search: the inputs are repository files, named in the argv relative to `cwd`. */
export interface VbrCase {
  readonly id: string;
  readonly program: SearchCase['program'];
  readonly options: readonly string[];
  /** The directory that the argv's file names are relative to: repository-relative, or absolute. */
  readonly cwd: string;
  readonly query: string;
  readonly subject: string;
  /** NCBI BLAST+ 2.17.0 stdout SHA-256 by format (every V-BR cell is certified, S08). */
  readonly frozen: Partial<Record<OutputFormat, { readonly sha256: string }>>;
}

/** The search of a V-BR case at `threads`. */
export function searchOf(vbr: VbrCase, threads: number | 'auto'): SearchCase {
  const file = (name: string) => ({ url: `/__files/${vbr.cwd === '.' ? '' : `${vbr.cwd}/`}${name}`, name });
  return {
    id: `${vbr.id} n=${threads}`,
    program: vbr.program,
    options: vbr.options,
    query: file(vbr.query),
    subject: file(vbr.subject),
    threads,
  };
}

/**
 * Cases from LOSAT/tests/outfmt0_manifest.tsv (run from LOSAT/). Rows of the same search
 * (program, task, inputs and options) that freeze different formats become one case.
 */
export function manifestCase(id: string, rowIds: readonly string[]): VbrCase {
  const lines = readFileSync(join(REPOSITORY, 'LOSAT/tests/outfmt0_manifest.tsv'), 'utf8').split('\n');
  const header = lines.find((line) => line.startsWith('fixture_id\t'))!.split('\t');
  const rows = lines
    .filter((line) => line !== '' && !line.startsWith('#') && !line.startsWith('fixture_id\t'))
    .map((line) => Object.fromEntries(line.split('\t').map((value, i) => [header[i]!, value])) as Record<string, string>);
  let vbr: VbrCase | undefined;
  const frozen: { [F in OutputFormat]?: { sha256: string } } = {};
  for (const rowId of rowIds) {
    const row = rows.find((candidate) => candidate['fixture_id'] === rowId);
    if (row === undefined) throw new Error(`no fixture ${rowId} in outfmt0_manifest.tsv`);
    // A row under an approved deviation (for example a genetic code) is not plain NCBI bytes.
    if (row['contract']) throw new Error(`${rowId} has the contract ${row['contract']}; V-BR uses plain fixtures`);
    const program = row['program'] as SearchCase['program'];
    const options = [...(row['task'] ? ['-task', row['task']] : []), ...(row['extra_args'] ? row['extra_args'].split(' ') : [])];
    const next: VbrCase = { id, program, options, cwd: 'LOSAT', query: row['query']!, subject: row['subject']!, frozen };
    if (vbr !== undefined && JSON.stringify([vbr.program, vbr.options, vbr.query, vbr.subject]) !== JSON.stringify([program, options, next.query, next.subject])) {
      throw new Error(`${rowId} is not the same search as the other rows of ${id}`);
    }
    vbr = next;
    const format = Number(row['outfmt'] || '0') as OutputFormat;
    frozen[format] = { sha256: row['stdout_sha256']! };
  }
  return vbr!;
}

/**
 * The V-BR cells (plan §7 S09): BLASTP, TBLASTN, BLASTN and TBLASTX, every format each
 * supports (0, 6, 7), compared with the native CLI and with NCBI's frozen bytes. BLASTX
 * joins in SX.
 */
export function vbrCases(): VbrCase[] {
  return [
    manifestCase('compact.blastn', ['compact.blastn']),
    manifestCase('multi.blastn', ['multi.blastn']),
    manifestCase('mask.lcase.megablast', ['mask.lcase.megablast']),
    manifestCase('width.blastp', ['width.blastp']),
    manifestCase('width.tblastn', ['width.tblastn']),
    manifestCase('tblastx.multi', ['tblastx.multi.0', 'tblastx.multi.7']),
    {
      // A genetic code other than the standard one (TLOSAN Stage D fixtures, V-ABI quick).
      id: 'tblastn.code4',
      program: 'tblastn',
      options: ['-task', 'tblastn', '-db_gencode', '4'],
      cwd: '.',
      query: 'docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures/code4.faa',
      subject: 'docs/evidence/tlosan_stage_d/all_codes_20260925/fixtures/code4.fna',
      frozen: {},
    },
    {
      // A query of 64,069 residues: TBLASTN searches it in four batches of 20,000 residues,
      // each with its own thread pool, so threads of one pool start while those of the
      // previous one end (the ThreadHost waits for them, S09).
      id: 'tblastn.batches',
      program: 'tblastn',
      options: [],
      cwd: 'LOSAT',
      query: 'tests/fasta/SicyWSV.faa',
      subject: 'tests/fasta/outfmt0/tblastx_multi_subject.fasta',
      frozen: {},
    },
    {
      // Several subjects that the threaded module reduces in its wasm-threads-only path
      // (verification_cells.tsv, the TBLASTX threaded command-WASI row; V-ABI quick).
      id: 'tblastx.threaded-direct',
      program: 'tblastx',
      options: [],
      cwd: '.',
      query: 'docs/evidence/losat_web_e1c/run-20260928T203929Z/tblastx-threaded-direct/q.fna',
      subject: 'docs/evidence/losat_web_e1c/run-20260928T203929Z/tblastx-threaded-direct/s.fna',
      frozen: {},
    },
  ];
}
