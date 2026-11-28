import { describe, expect, it } from 'vitest';
import { sha256Hex } from 'losat';
import type { EngineRunRequest } from 'losat/ports';
import { Coordinator, isToolRun } from '../../src/application/coordinator';
import { holdsInputs, SearchDraft } from '../../src/application/draft';
import { ResultsBrowser } from '../../src/application/results';
import { buildArgv } from '../../src/domain/argv';
import { assembleSearchDirectory } from '../../src/domain/blast-db';
import { hspTable } from '../../src/domain/hsp-table';
import { MODE_TABS, modeTabOf, programById, querySourcesOf } from '../../src/domain/programs';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import { runEnd } from './support/run-end';
import { NodeInflater } from './support/node-inflater';

function makeHarness() {
  const store = new MemoryBlockStore();
  let token = 0;
  const data = new DataService({
    store,
    scanner: new FakeScanner(),
    checker: new FakeInputChecker(),
    inflater: new NodeInflater(),
    digest: sha256Hex,
    newToken: () => `token-${++token}`,
    cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
  });
  const engine = new FakeEngine();
  const seenRequests: EngineRunRequest[] = [];
  const wrappedEngine = Object.assign(Object.create(Object.getPrototypeOf(engine)), engine, {
    run: (request: EngineRunRequest, output: MessagePort, onPhase: (phase: 'preparing' | 'running' | 'finalizing') => void) => {
      seenRequests.push(request);
      return engine.run(request, output, onPhase);
    },
  }) as FakeEngine;
  let idCounter = 1;
  const coordinator = new Coordinator({
    engine: wrappedEngine,
    data,
    downloader: { save: () => undefined },
    now: () => 1_764_100_000_000,
    newRunId: () => `run-${idCounter++}`,
    digest: sha256Hex,
  });
  const draft = new SearchDraft({
    engine: wrappedEngine,
    data,
    enqueueAll: (requests) => coordinator.enqueueAll(requests),
    fetcher: { fetch: async () => ({ ok: false, reason: 'not-found' }) },
    now: () => 1_764_100_000_000,
    debounceMs: 0,
  });
  const results = new ResultsBrowser({
    data,
    describe: (p) => wrappedEngine.describe(p),
    runs: coordinator.state,
    verification: { ncbi: '2.17.0', sources: [], programs: {} },
  });
  return { draft, coordinator, results, seenRequests };
}

/** Submits and waits for the run to end: the submit resolves when the run is accepted, not when it is done. */
async function submitFileQuery(draft: SearchDraft, coordinator: Coordinator, seenRequests: EngineRunRequest[]): Promise<string> {
  await draft.idle();
  const submitted = await draft.submit();
  expect(submitted.ok, submitted.ok ? '' : submitted.message).toBe(true);
  const runId = submitted.runIds![0]!;
  await runEnd(coordinator, runId);
  expect(seenRequests[0]!.query).toBeUndefined();
  return runId;
}

describe('psiblast & tblastn -in_pssm domain, draft, coordinator, and results', () => {
  it('exposes expected query sources per program, keeps 5 mode tabs with psiblast under blastp, and builds -in_pssm / -in_msa argv cleanly', () => {
    expect(MODE_TABS.map((p) => p.id)).toEqual(['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']);
    expect(modeTabOf('psiblast')).toBe('blastp');
    expect(modeTabOf('blastp')).toBe('blastp');
    expect(querySourcesOf('psiblast')).toEqual(['fasta', 'msa', 'pssm']);
    expect(querySourcesOf('tblastn')).toEqual(['fasta', 'pssm']);
    expect(querySourcesOf('blastp')).toEqual(['fasta']);
    expect(programById('psiblast').label).toBe('PSIBLAST');

    const pssmArgv = buildArgv({
      program: 'tblastn',
      queryFile: { flag: '-in_pssm', name: 'seed.asn' },
      subjectName: 'genome.fna',
      parameters: [['-evalue', '1e-3']],
    });
    expect(pssmArgv).toEqual(['tblastn', '-in_pssm', 'seed.asn', '-subject', 'genome.fna', '-evalue', '1e-3']);

    for (const [prog, flag, name] of [
      ['psiblast', '-in_pssm', 'x.asn'],
      ['psiblast', '-in_msa', 'x.aln'],
      ['tblastn', '-in_pssm', 'x.asn'],
    ] as const) {
      expect(() =>
        buildArgv({ program: prog, queryName: 'q.fa', subjectName: 's.fa', parameters: [], extra: [flag, name] }),
      ).toThrow(`${flag} is set in the form: Query`);
    }
  });

  it('mounts queryFile into run directory and stores queryLength/round in HspTable', () => {
    const pssm = new File(['PssmWithParameters ::= {}'], 'seed.asn');
    const assembled = assembleSearchDirectory({ queryFile: pssm });
    expect(assembled.files.map((f) => f.name)).toEqual(['seed.asn']);

    const table = hspTable([
      {
        index: 0,
        q_idx: 0,
        s_idx: 0,
        rank: 0,
        raw_score: 120,
        bit_score: 55.2,
        e_value: 1e-10,
        q_start: 1,
        q_end: 80,
        s_start: 5,
        s_end: 84,
        query_frame: null,
        subject_frame: null,
        query_length: 142,
        subject_length: 300,
        out6: [0, 40],
        out0: null,
        out0_subject: null,
        round: 2,
      },
    ]);
    expect(table.queryLength[0]).toBe(142);
    expect(table.round[0]).toBe(2);
  });

  it('submits psiblast with -in_msa, -num_iterations, -out_pssm, -out_ascii_pssm and opens query-file results cleanly', async () => {
    const { draft, coordinator, results, seenRequests } = makeHarness();
    draft.setProgram('psiblast');
    draft.setQuerySource('msa');
    const msa = new File(['CLUSTAL W\n\ns1 MKTAYIAKQRQISFVKSHFSRQ\ns2 MKTAYIAKQRQISFVKSHFSRQ\n'], 'align.aln');
    draft.setQueryFile(msa);
    draft.setMsaMasterIdx('2');
    draft.setPaste('subject', '>s1\nMKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQ\n');
    draft.setField('-num_iterations', '2');
    draft.setField('-out_pssm', 'out.asn');
    draft.setField('-out_ascii_pssm', 'out.ascii');
    draft.setField('-save_pssm_after_last_round', true);
    expect(holdsInputs(draft.state.get())).toBe(true);
    const runId = await submitFileQuery(draft, coordinator, seenRequests);

    const run = coordinator.state.get().runs.find((r) => r.snapshot.runId === runId);
    expect(run?.status).toBe('completed');
    expect(seenRequests).toHaveLength(1);
    expect(seenRequests[0]!.files?.map((f) => f.name)).toEqual(['align.aln']);
    expect(seenRequests[0]!.argv).toEqual(
      expect.arrayContaining(['psiblast', '-in_msa', 'align.aln', '-msa_master_idx', '2', '-num_iterations', '2', '-out_pssm', 'out.asn', '-out_ascii_pssm', 'out.ascii', '-save_pssm_after_last_round']),
    );
    if (run && !isToolRun(run)) {
      expect(run.snapshot.queryFile).toEqual({ flag: '-in_msa', name: 'align.aln', size: msa.size });
      expect(run.record.runFiles?.map((f) => [f.name, f.option])).toEqual([
        ['out.asn', 'out_pssm'],
        ['out.ascii', 'out_ascii_pssm'],
      ]);
    }

    await results.open(runId);
    expect(results.state.get().phase, results.state.get().message).toBe('ready');
    expect(results.state.get().queries).toHaveLength(1);
    expect(results.state.get().queries[0]!.length).toBe(120);
    expect(results.state.get().subjects.length).toBeGreaterThan(0);
  });

  it('switches between blastp tasks and psiblast via selectProteinAlgorithm', () => {
    const { draft } = makeHarness();
    draft.setProgram('blastp');
    draft.selectProteinAlgorithm('psiblast');
    expect(draft.state.get().program).toBe('psiblast');
    draft.selectProteinAlgorithm('blastp-fast');
    expect(draft.state.get().program).toBe('blastp');
    expect(draft.state.get().values.blastp['-task']).toBe('blastp-fast');
  });

  it('resets incompatible querySource when switching programs and supports tblastn -in_pssm', async () => {
    const { draft, coordinator, seenRequests } = makeHarness();
    draft.setProgram('psiblast');
    draft.setQuerySource('msa');
    expect(draft.state.get().querySource).toBe('msa');
    draft.setProgram('tblastn');
    expect(draft.state.get().querySource).toBe('fasta');

    draft.setQuerySource('pssm');
    draft.setQueryFile(new File(['PssmWithParameters ::= {}'], 'query.asn'));
    draft.setPaste('subject', '>nt1\nATGAAAACCGCTTATATTGCTAAACAACGTCAAATTAGCTTTGTTAAAAGCCATTTTTCTCGTCAACTG\n');
    await submitFileQuery(draft, coordinator, seenRequests);
    expect(seenRequests[0]!.argv.slice(0, 3)).toEqual(['tblastn', '-in_pssm', 'query.asn']);
    expect(coordinator.state.get().runs[0]!.status).toBe('completed');
  });
});
