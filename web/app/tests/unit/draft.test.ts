import { describe, expect, it } from 'vitest';
import { Coordinator, type EnqueueAllResult, type SearchRequest } from '../../src/application/coordinator';
import { SearchDraft, type DraftSource } from '../../src/application/draft';
import type { InputRole } from '../../src/domain/programs';
import { sha256Hex } from '../../src/infra/browser/platform';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import type { ValidationResult } from '../../src/ports/engine';

function setup(
  options: {
    validate?: (argv: readonly string[]) => ValidationResult;
    useCoordinator?: boolean;
    checkInput?: DataService['checkInput'];
  } = {},
) {
  let token = 0;
  const data = new DataService({
    store: new MemoryBlockStore(),
    scanner: new FakeScanner(),
    checker: new FakeInputChecker(),
    digest: sha256Hex,
    newToken: () => `token-${++token}`,
    cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
  });
  const checks: string[] = [];
  const checkInput = data.checkInput.bind(data);
  data.checkInput = async (program, role, revisionIds) => {
    checks.push(`${program}|${role}|${revisionIds.join(',')}`);
    return (options.checkInput ?? checkInput)(program, role, revisionIds);
  };
  const fake = new FakeEngine();
  const validated: Array<readonly string[]> = [];
  const engine = {
    describe: fake.describe.bind(fake),
    validate: async (argv: readonly string[]): Promise<ValidationResult> => {
      validated.push(argv);
      return options.validate?.(argv) ?? { ok: true };
    },
  };
  let id = 0;
  const coordinator = new Coordinator({
    engine: { ...engine, run: fake.run.bind(fake), cancel: fake.cancel.bind(fake) },
    data,
    downloader: { save: () => undefined },
    now: () => 0,
    newRunId: () => `run-${++id}`,
  });
  const requests: SearchRequest[][] = [];
  const draft = new SearchDraft({
    engine,
    data,
    enqueueAll: async (batch): Promise<EnqueueAllResult> => {
      requests.push([...batch]);
      return options.useCoordinator ? coordinator.enqueueAll(batch) : { ok: true, runIds: batch.map((_, i) => `r${i}`) };
    },
    now: () => 0,
    debounceMs: 0,
  });
  return { draft, data, requests, validated, coordinator, checks };
}

const file = (name: string, text: string) => new File([text], name, { type: 'text/plain' });
const sourceOf = (draft: SearchDraft, role: InputRole, index = 0): DraftSource => draft.state.get()[role].sources[index]!;
const datasetOf = (request: SearchRequest, role: InputRole) => {
  const input = request[role];
  if (!('dataset' in input)) throw new Error('expected a dataset input');
  return input.dataset;
};

describe('SearchDraft inputs', () => {
  it('indexes pasted text and checks it with the engine', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>q1\nACGT\n>q2\nGGCC\n');
    await draft.idle();
    const source = sourceOf(draft, 'query');
    expect(source).toMatchObject({ origin: 'paste', name: 'query.fa', status: 'ready', excluded: [] });
    expect(source.base?.records.map((r) => r.id)).toEqual(['q1', 'q2']);
    expect(source.check).toEqual({ state: 'ok', records: 2 });
    expect(source.head).toBe('>q1\nACGT\n>q2\nGGCC\n');
  });

  it('replaces the paste source when the text changes, and removes it when the text is emptied', async () => {
    const { draft } = setup();
    draft.setPaste('subject', '>s1\nACGT\n');
    await draft.idle();
    draft.setPaste('subject', '>s2\nACGT\n');
    await draft.idle();
    expect(draft.state.get().subject.sources).toHaveLength(1);
    expect(sourceOf(draft, 'subject').base?.records[0]?.id).toBe('s2');
    draft.setPaste('subject', '');
    await draft.idle();
    expect(draft.state.get().subject.sources).toHaveLength(0);
  });

  it('reports the index scan error of a text without a defline, and adds one on request', async () => {
    const { draft } = setup();
    draft.setPaste('query', 'ACGTACGT\n');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'failed', error: 'Expected > at record start.' });
    draft.addDefline('query');
    await draft.idle();
    expect(draft.state.get().query.paste).toBe('>pasted_query\nACGTACGT\n');
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready' });
  });

  it('shows the engine refusal of a record and maps it to the record table; excluding it passes the check', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>a\nACGT\n>b\nAC!GT\n>c\nACGT\n');
    await draft.idle();
    const refused = sourceOf(draft, 'query').check;
    expect(refused).toMatchObject({ state: 'refused', record: 1 });
    expect(refused?.state === 'refused' && refused.message).toContain('query record 2 (b)');
    draft.setIncluded('query', 'paste', [1], false);
    await draft.idle();
    expect(sourceOf(draft, 'query').check).toEqual({ state: 'ok', records: 2 });
  });

  it('maps a refused record among the included records only', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>a\nACGT\n>b\nACGT\n>c\nAC!GT\n');
    await draft.idle();
    draft.setIncluded('query', 'paste', [0], false);
    await draft.idle();
    // The engine read b and c (record 2 of its input is c, index 2 in the table).
    expect(sourceOf(draft, 'query').check).toMatchObject({ state: 'refused', record: 2 });
  });

  it('checks the inputs again for another program', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>a\nACGT\n');
    await draft.idle();
    draft.setProgram('blastx');
    await draft.idle();
    // BLASTX cannot be searched until SX: nothing is read again or checked.
    expect(sourceOf(draft, 'query').status).toBe('ready');
    expect(sourceOf(draft, 'query').check).toBeUndefined();
    draft.setProgram('tblastx');
    await draft.idle();
    expect(sourceOf(draft, 'query').check).toEqual({ state: 'ok', records: 1 });
  });
});

describe('SearchDraft submit', () => {
  it('refuses to queue without inputs, with every record excluded, or with a refused input', async () => {
    const { draft, requests } = setup();
    expect(await draft.submit()).toEqual({ ok: false, message: 'Add the query sequences: paste them or open FASTA files.' });
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nAC!GT\n');
    await draft.idle();
    const refused = await draft.submit();
    expect(refused.ok).toBe(false);
    expect(!refused.ok && refused.message).toMatch(/^Subject \(pasted\): subject record 1 \(s\) has/);
    draft.setIncluded('subject', 'paste', [0], false);
    expect(await draft.submit()).toEqual({ ok: false, message: 'Every subject record is excluded.' });
    expect(requests).toHaveLength(0);
    expect(draft.state.get().message).toEqual({ kind: 'error', text: 'Every subject record is excluded.' });
  });

  it('refuses BLASTX until SX', async () => {
    const { draft } = setup();
    draft.setProgram('blastx');
    const result = await draft.submit();
    expect(result.ok).toBe(false);
    expect(!result.ok && result.message).toMatch(/^BLASTX joins LOSAT Web after its certification/);
  });

  it('queues one combined search of all files, named as plan §5.3 says', async () => {
    const { draft, requests } = setup();
    draft.setPaste('query', '>q\nACGT\n');
    draft.addFiles('subject', [file('a.fa', '>a\nACGT\n'), file('b.fa', '>b\nACGT')]);
    await draft.idle();
    expect(draft.runCount()).toBe(1);
    expect(await draft.submit()).toMatchObject({ ok: true });
    expect(requests).toHaveLength(1);
    const [request] = requests[0]!;
    expect(datasetOf(request!, 'query').name).toBe('query.fa');
    expect(datasetOf(request!, 'subject')).toEqual({
      name: 'combined_subject.fa',
      revisionIds: [sourceOf(draft, 'subject', 0).base!.revisionId, sourceOf(draft, 'subject', 1).base!.revisionId],
    });
    expect(draft.state.get().message).toEqual({ kind: 'info', text: 'Added to the queue.' });
  });

  it('queues separate searches as every query input against every subject input', async () => {
    const { draft, requests } = setup();
    draft.addFiles('query', [file('q1.fa', '>q1\nACGT\n'), file('q2.fa', '>q2\nACGT\n')]);
    draft.addFiles('subject', [file('s1.fa', '>s1\nACGT\n'), file('s2.fa', '>s2\nACGT\n'), file('s3.fa', '>s3\nACGT\n')]);
    await draft.idle();
    draft.setMode('subject', 'separate');
    expect(draft.runCount()).toBe(3);
    draft.setMode('query', 'separate');
    expect(draft.runCount()).toBe(6);
    await draft.submit();
    const names = requests[0]!.map((r) => `${datasetOf(r, 'query').name}|${datasetOf(r, 'subject').name}`);
    expect(names).toEqual(['q1.fa|s1.fa', 'q1.fa|s2.fa', 'q1.fa|s3.fa', 'q2.fa|s1.fa', 'q2.fa|s2.fa', 'q2.fa|s3.fa']);
    expect(draft.state.get().message).toEqual({ kind: 'info', text: 'Added 6 runs to the queue as a group.' });
  });

  it('uses the same revision for the same selection, so that a kept subject is the same input (R1)', async () => {
    const { draft, requests } = setup();
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s1\nACGT\n>s2\nACGT\n>s3\nACGT\n');
    await draft.idle();
    const base = sourceOf(draft, 'subject').base!.revisionId;
    await draft.submit();
    draft.setPaste('query', '>other\nGGGG\n');
    await draft.submit();
    draft.setIncluded('subject', 'paste', [1], false);
    await draft.submit();
    await draft.submit();
    const subjects = requests.map((batch) => datasetOf(batch[0]!, 'subject').revisionIds[0]);
    expect(subjects[0]).toBe(base);
    expect(subjects[1]).toBe(base);
    expect(subjects[2]).not.toBe(base);
    expect(subjects[3]).toBe(subjects[2]);
    const queries = requests.map((batch) => datasetOf(batch[0]!, 'query').revisionIds[0]);
    expect(queries[1]).not.toBe(queries[0]);
  });

  it('writes the region of a role that has one record, and only then', async () => {
    const { draft, requests } = setup();
    draft.setPaste('query', '>q1\nACGTACGTAC\n>q2\nACGT\n');
    draft.setPaste('subject', '>s\nACGTACGTACGTACGTACGT\n');
    await draft.idle();
    expect(draft.regionRecord('query')).toBeUndefined();
    expect(draft.regionRecord('subject')?.id).toBe('s');
    // A role of two records has no region.
    draft.setRegion('query', { start: '2', stop: '5' });
    expect(draft.region('query')).toBeUndefined();
    draft.setRegion('subject', { start: ' 3', stop: '012 ' });
    expect(draft.parameters()).toEqual([['-subject_loc', '3-12']]);
    draft.setIncluded('query', 'paste', [1], false);
    await draft.idle();
    draft.setRegion('query', { start: '2', stop: '5' });
    expect(draft.parameters()).toEqual([
      ['-query_loc', '2-5'],
      ['-subject_loc', '3-12'],
    ]);
    draft.setRegion('subject', { start: '3', stop: '21' });
    const result = await draft.submit();
    expect(result).toEqual({
      ok: false,
      message: 'Subject region: The stop must be between 1 and 20, the length of the record.',
    });
    draft.setRegion('subject', undefined);
    await draft.submit();
    expect(requests[0]![0]!.parameters).toEqual([['-query_loc', '2-5']]);
  });

  it('writes only the values the user set and that differ from the engine default', async () => {
    const { draft, requests } = setup();
    draft.setProgram('tblastx');
    await draft.idle();
    draft.setField('-evalue', '10'); // the engine's default
    draft.setField('-word_size', '');
    draft.setField('-seg', 'no');
    draft.setField('-max_target_seqs', ' 50 ');
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT\n');
    await draft.submit();
    expect(requests[0]![0]!.parameters).toEqual([
      ['-max_target_seqs', '50'],
      ['-seg', 'no'],
    ]);
    // Each program keeps its own values.
    draft.setProgram('blastn');
    draft.setField('-task', 'dc-megablast');
    draft.setField('-template_type', 'optimal');
    draft.setField('-template_length', '21');
    draft.setField('-lcase_masking', true);
    expect(draft.parameters()).toEqual([
      ['-task', 'dc-megablast'],
      ['-template_type', 'optimal'],
      ['-template_length', '21'],
      ['-lcase_masking', true],
    ]);
    draft.setField('-task', 'blastn');
    expect(draft.parameters()).toEqual([
      ['-task', 'blastn'],
      ['-lcase_masking', true],
    ]);
    draft.setProgram('tblastx');
    expect(draft.parameters()).toEqual([
      ['-max_target_seqs', '50'],
      ['-seg', 'no'],
    ]);
  });

  it("validates the draft's argv with the engine and shows its message", async () => {
    const { draft, validated } = setup({
      validate: (argv) => (argv.includes('-word_size') ? { ok: false, message: 'BLAST query/options error: bad' } : { ok: true }),
    });
    await draft.idle();
    expect(draft.state.get().validation).toEqual({ state: 'ok' });
    draft.setField('-word_size', '3');
    await draft.idle();
    expect(draft.state.get().validation).toEqual({ state: 'invalid', message: 'BLAST query/options error: bad' });
    expect(validated.at(-1)).toEqual(['blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-word_size', '3']);
  });

  it('queues runs in the coordinator, a group with shared subject bytes', async () => {
    const { draft, coordinator } = setup({ useCoordinator: true });
    draft.addFiles('query', [file('q1.fa', '>q1\nACGT\n'), file('q2.fa', '>q2\nACGT\n')]);
    draft.setPaste('subject', '>s\nACGT\n');
    await draft.idle();
    draft.setMode('query', 'separate');
    const result = await draft.submit();
    expect(result.ok).toBe(true);
    const runs = coordinator.state.get().runs;
    expect(runs.map((run) => run.snapshot.group?.position)).toEqual([1, 2]);
    expect(runs[0]!.snapshot.group?.groupId).toBe(runs[1]!.snapshot.group?.groupId);
    expect(runs[0]!.snapshot.subject.bytes).toBe(runs[1]!.snapshot.subject.bytes);
    expect(runs.map((run) => run.snapshot.query.name)).toEqual(['q1.fa', 'q2.fa']);
  });

  it('reports a failure to queue instead of dropping it', async () => {
    const { draft } = setup({ useCoordinator: true, validate: () => { throw new Error('the Data worker stopped'); } });
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT\n');
    const result = await draft.submit();
    expect(result).toEqual({ ok: false, message: 'The options could not be checked: the Data worker stopped' });
    expect(draft.state.get()).toMatchObject({ submitting: false, message: { kind: 'error' } });
  });

  it("maps a subject record that the engine refuses, and reports it in short above the button", async () => {
    const { draft } = setup();
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s1\nACGT\n>s2\nAC!GT\n');
    await draft.idle();
    expect(sourceOf(draft, 'subject').check).toMatchObject({ state: 'refused', record: 1 });
    expect(draft.readiness()).toBe(
      'Subject (pasted): the engine refuses a record. Exclude it or correct the input to search.',
    );
    draft.setIncluded('subject', 'paste', [1], false);
    await draft.idle();
    expect(draft.readiness()).toBeUndefined();
  });

  it('does not let a check that could not run keep the search from the queue', async () => {
    const { draft, requests } = setup({
      checkInput: async () => {
        throw new Error('the engine could not allocate 9 bytes for an input');
      },
    });
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT\n');
    await draft.idle();
    expect(sourceOf(draft, 'query').check).toEqual({
      state: 'error',
      message: 'the engine could not allocate 9 bytes for an input',
    });
    expect(draft.readiness()).toBeUndefined();
    expect((await draft.submit()).ok).toBe(true);
    expect(requests).toHaveLength(1);
  });

  it('checks a selection once: going back to an earlier selection or program uses its verdict', async () => {
    const { draft, checks } = setup();
    draft.setPaste('query', '>a\nACGT\n>b\nACGT\n');
    await draft.idle();
    draft.setIncluded('query', 'paste', [1], false);
    await draft.idle();
    draft.setIncluded('query', 'paste', [1], true);
    await draft.idle();
    draft.setProgram('tblastx');
    await draft.idle();
    draft.setProgram('blastn');
    await draft.idle();
    expect(checks.map((check) => check.split('|').slice(0, 2).join('|'))).toEqual([
      'blastn|query',
      'blastn|query',
      'tblastx|query',
    ]);
  });

  it('shows a pasted text as being read at once, before it is indexed', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>a\nACGT\n');
    await draft.idle();
    draft.setPaste('query', '>b\nACGT\n');
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'indexing' });
    await draft.idle();
    expect(sourceOf(draft, 'query').base?.records[0]?.id).toBe('b');
  });

  it("keeps a region with its record: another record leaves it aside", async () => {
    const { draft } = setup();
    draft.setPaste('subject', '>s1\nACGTACGTACGTACGTACGT\n');
    await draft.idle();
    draft.setRegion('subject', { start: '3', stop: '12' });
    expect(draft.parameters()).toEqual([['-subject_loc', '3-12']]);
    draft.setPaste('subject', '>s2\nACGTACGT\n');
    await draft.idle();
    expect(draft.region('subject')).toBeUndefined();
    expect(draft.parameters()).toEqual([]);
  });
});
