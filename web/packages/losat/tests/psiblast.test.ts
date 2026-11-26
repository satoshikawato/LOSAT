// Unit and Wasm reactor parity tests for psiblast and tblastn -in_pssm (PSW1).
import { spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { closeSync, mkdtempSync, openAsBlob, openSync, readdirSync, readFileSync, readSync, rmSync, statSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { namedRunFiles, PROGRAM_INPUTS, readStoredZip, runFileOption, type DirectoryFile } from '../src/domain/index.ts';
import { DISCARD_STREAM, ROLE_QUERY, ROLE_SUBJECT } from '../src/infra/reactor/abi.ts';
import { compileReactor, instantiateSerial } from '../src/infra/reactor/instance.ts';
import { RunDirectory, type BlobReader } from '../src/infra/reactor/run-directory.ts';

const LOSAT = fileURLToPath(new URL('../../../../LOSAT', import.meta.url));
const REACTORS = process.env['LOSAT_WEB_REACTORS'];
const NATIVE = process.env['LOSAT_NATIVE_BIN'];

function sha256(bytes: Uint8Array): string {
  return createHash('sha256').update(bytes).digest('hex');
}

function osReader(paths: Map<Blob, string>): BlobReader {
  return (blob, start, end) => {
    const bytes = new Uint8Array(end - start);
    const descriptor = openSync(paths.get(blob)!, 'r');
    try {
      const got = readSync(descriptor, bytes, 0, bytes.length, start);
      return bytes.subarray(0, got);
    } finally {
      closeSync(descriptor);
    }
  };
}

async function mountFiles(relPaths: readonly string[], paths: Map<Blob, string>): Promise<DirectoryFile[]> {
  const out: DirectoryFile[] = [];
  for (const rel of relPaths) {
    const full = join(LOSAT, rel);
    const blob = await openAsBlob(full);
    paths.set(blob, full);
    out.push({ name: rel, file: blob });
  }
  return out;
}

async function mountFolder(relDir: string, paths: Map<Blob, string>): Promise<DirectoryFile[]> {
  const out: DirectoryFile[] = [];
  for (const name of readdirSync(join(LOSAT, relDir))) {
    const full = join(LOSAT, relDir, name);
    if (!statSync(full).isFile()) continue;
    const blob = await openAsBlob(full);
    paths.set(blob, full);
    out.push({ name: `${relDir}/${name}`, file: blob });
  }
  return out;
}

describe('psiblast package domain support', () => {
  it('lists psiblast with protein query and protein subject', () => {
    expect(PROGRAM_INPUTS.psiblast).toEqual({ query: 'protein', subject: 'protein' });
  });

  it('labels generated PSSM and ASCII PSSM stream-8 files including -save_each_pssm suffixes', () => {
    const start = {
      logFile: null,
      logFileName: null,
      strategyFile: null,
      outPssm: 'out.asn',
      outAsciiPssm: 'out.ascii',
      saveEachPssm: true,
    };
    const raw = [
      { name: 'out.asn.1', bytes: new Uint8Array([1]) },
      { name: 'out.asn.2', bytes: new Uint8Array([2]) },
      { name: 'out.ascii.1', bytes: new Uint8Array([3]) },
      { name: 'out.ascii.2', bytes: new Uint8Array([4]) },
    ];
    const labeled = namedRunFiles(raw, start);
    expect(labeled.map((f) => [f.name, runFileOption(f)])).toEqual([
      ['out.asn.1', 'out_pssm'],
      ['out.asn.2', 'out_pssm'],
      ['out.ascii.1', 'out_ascii_pssm'],
      ['out.ascii.2', 'out_ascii_pssm'],
    ]);
  });
});

describe.skipIf(REACTORS === undefined)('psiblast and tblastn -in_pssm through the serial Wasm reactor', () => {
  it('runs psiblast (-subject, -db, -in_msa, -in_pssm, PSSM outputs) and tblastn -in_pssm cleanly', async () => {
    const { module } = await compileReactor(readFileSync(join(REACTORS!, 'losat-web-serial.wasm')));
    const paths = new Map<Blob, string>();
    const directory = new RunDirectory(undefined, osReader(paths));
    const { abi } = await instantiateSerial(module, directory);

    const desc = abi.describe('psiblast') as { program: string; parameters: Array<{ flag: string }> };
    expect(desc.program).toBe('psiblast');
    expect(desc.parameters.some((p) => p.flag === '-num_iterations')).toBe(true);
    expect(desc.parameters.some((p) => p.flag === '-in_msa')).toBe(true);
    expect(desc.parameters.some((p) => p.flag === '-in_pssm')).toBe(true);
    expect(desc.parameters.some((p) => p.flag === '-out_pssm')).toBe(true);
    expect(desc.parameters.some((p) => p.flag === '-out_ascii_pssm')).toBe(true);

    // 1) psiblast -query + -subject with -num_iterations 2 and PSSM outputs
    const q1 = abi.register('psiblast', ROLE_QUERY, true, readFileSync(join(LOSAT, 'tests/fasta/psiblast/q1b_gi.faa')));
    const s1 = abi.register('psiblast', ROLE_SUBJECT, true, readFileSync(join(LOSAT, 'tests/fasta/psiblast/q1_subj_ids.faa')));
    const argv1 = [
      'psiblast',
      '-query',
      'tests/fasta/psiblast/q1b_gi.faa',
      '-subject',
      'tests/fasta/psiblast/q1_subj_ids.faa',
      '-num_iterations',
      '2',
      '-out_pssm',
      'chk.asn',
      '-out_ascii_pssm',
      'chk.ascii',
      '-save_pssm_after_last_round',
      '-num_threads',
      '1',
    ];
    const st1 = abi.start(argv1);
    expect(st1.ok).toBe(true);
    if (st1.ok) {
      expect(st1.start.outPssm).toBe('chk.asn');
      expect(st1.start.outAsciiPssm).toBe('chk.ascii');
    }
    expect(abi.validate(argv1).ok).toBe(true);
    const streams1 = new Map<number, Uint8Array[]>();
    abi.run(argv1, q1.handle, s1.handle, (stream, bytes) => {
      if (stream === DISCARD_STREAM) return;
      const list = streams1.get(stream) ?? [];
      list.push(bytes.slice());
      streams1.set(stream, list);
    });
    abi.release(q1.handle);
    abi.release(s1.handle);
    const out6_1 = Buffer.concat(streams1.get(6) ?? []);
    expect(out6_1.length).toBeGreaterThan(0);
    const files1 = readStoredZip(Buffer.concat(streams1.get(8) ?? []));
    expect(files1.map((f) => f.name).sort()).toEqual(['chk.ascii', 'chk.asn']);
    if (NATIVE) {
      const tmp = mkdtempSync(join(tmpdir(), 'losat-pkg-psi-'));
      try {
        const res = spawnSync(
          NATIVE,
          [
            'psiblast',
            '-query',
            'tests/fasta/psiblast/q1b_gi.faa',
            '-subject',
            'tests/fasta/psiblast/q1_subj_ids.faa',
            '-num_iterations',
            '2',
            '-out_pssm',
            join(tmp, 'chk.asn'),
            '-out_ascii_pssm',
            join(tmp, 'chk.ascii'),
            '-save_pssm_after_last_round',
            '-outfmt',
            '6',
            '-num_threads',
            '1',
          ],
          { cwd: LOSAT, maxBuffer: 1 << 26 },
        );
        expect(res.status).toBe(0);
        expect(sha256(out6_1)).toBe(sha256(res.stdout));
        for (const f of files1) {
          expect(sha256(f.bytes), f.name).toBe(sha256(readFileSync(join(tmp, f.name))));
        }
      } finally {
        rmSync(tmp, { recursive: true, force: true });
      }
    }

    // 2) psiblast -in_msa + -db
    const msaFiles = [
      ...(await mountFolder('tests/db_search', paths)),
      ...(await mountFiles(['tests/fixtures/msa/ps6_b_db_0.aln'], paths)),
    ];
    directory.open(msaFiles);
    const argv2 = ['psiblast', '-in_msa', 'tests/fixtures/msa/ps6_b_db_0.aln', '-db', 'tests/db_search/paja', '-num_threads', '1'];
    const st2 = abi.start(argv2);
    expect(st2.ok && st2.start.queryFile === 'tests/fixtures/msa/ps6_b_db_0.aln' && st2.start.database).toBe(true);
    expect(abi.validate(argv2).ok).toBe(true);
    const out6_2: Uint8Array[] = [];
    abi.run(argv2, 0, 0, (stream, bytes) => {
      if (stream === 6) out6_2.push(bytes.slice());
    });
    expect(Buffer.concat(out6_2).length).toBeGreaterThan(0);

    // 3) tblastn -in_pssm + -subject
    directory.open(await mountFiles(['tests/fixtures/pssm/q1b_gi.asn'], paths));
    const s3 = abi.register('tblastn', ROLE_SUBJECT, true, readFileSync(join(LOSAT, 'tests/fasta/psiblast/q1b_gi_nt.fna')));
    const argv3 = ['tblastn', '-in_pssm', 'tests/fixtures/pssm/q1b_gi.asn', '-subject', 'tests/fasta/psiblast/q1b_gi_nt.fna', '-num_threads', '1'];
    const st3 = abi.start(argv3);
    expect(st3.ok && st3.start.queryFile === 'tests/fixtures/pssm/q1b_gi.asn' && !st3.start.database).toBe(true);
    expect(abi.validate(argv3).ok).toBe(true);
    const out6_3: Uint8Array[] = [];
    abi.run(argv3, 0, s3.handle, (stream, bytes) => {
      if (stream === 6) out6_3.push(bytes.slice());
    });
    abi.release(s3.handle);
    expect(Buffer.concat(out6_3).length).toBeGreaterThan(0);
  }, 120_000);
});
