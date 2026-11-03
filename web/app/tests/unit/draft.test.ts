import { describe, expect, it } from 'vitest';
import { Coordinator, type EnqueueAllResult, type SearchRequest } from '../../src/application/coordinator';
import { SearchDraft, type DraftSource } from '../../src/application/draft';
import type { InputRole } from '../../src/domain/programs';
import { sha256Hex } from '../../src/infra/browser/platform';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeAccessionFetcher } from '../../src/infra/fake/fake-fetcher';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import type { EnginePrepareHint, ValidationResult } from '../../src/ports/engine';

function setup(
  options: {
    validate?: (argv: readonly string[]) => ValidationResult;
    useCoordinator?: boolean;
    checkInput?: DataService['checkInput'];
    fetcher?: FakeAccessionFetcher;
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
  const prepared: EnginePrepareHint[] = [];
  const engine = {
    describe: fake.describe.bind(fake),
    validate: async (argv: readonly string[]): Promise<ValidationResult> => {
      validated.push(argv);
      return options.validate?.(argv) ?? { ok: true };
    },
    prepare: (hint: EnginePrepareHint) => void prepared.push(hint),
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
    fetcher: options.fetcher ?? new FakeAccessionFetcher(),
    now: () => 0,
    debounceMs: 0,
  });
  return { draft, data, requests, validated, prepared, coordinator, checks };
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

  it('reads a text without a defline as one record without a title, and adds a defline on request', async () => {
    const { draft } = setup();
    draft.setPaste('query', 'ACGTACGT\n');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready', parser: 1 });
    expect(sourceOf(draft, 'query').base?.records.map((r) => [r.id, r.header_offset, r.sequence_offset, r.length])).toEqual([['', 0, 0, 8]]);
    draft.addDefline('query');
    await draft.idle();
    expect(draft.state.get().query.paste).toBe('>pasted_query\nACGTACGT\n');
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready' });
  });

  it('shows the engine refusal of a record and maps it to the record table; excluding it passes the check', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>a\nACGT\n>b\nACGT!ACGT\n>c\nACGT\n');
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
    draft.setPaste('query', '>a\nACGT\n>b\nACGT\n>c\nACGT!ACGT\n');
    await draft.idle();
    draft.setIncluded('query', 'paste', [0], false);
    await draft.idle();
    // The engine read b and c (record 2 of its input is c, index 2 in the table).
    expect(sourceOf(draft, 'query').check).toMatchObject({ state: 'refused', record: 2 });
  });

  it('checks the inputs again for another program, and reads them again when it reads the role with another kind', async () => {
    const { draft, checks } = setup();
    draft.setPaste('query', '>a\nMKLV\n');
    draft.setPaste('subject', '>s\n$$$$\n');
    await draft.idle();
    const revision = sourceOf(draft, 'query').base!.revisionId;
    expect(sourceOf(draft, 'query').base!.records[0]!.length).toBe(3);
    expect(sourceOf(draft, 'subject')).toMatchObject({ status: 'failed', parser: 1 });
    // BLASTX reads its query as BLASTN does (kind 1): the query is checked again, not read again.
    // Its subject is protein (kind 2): the subject is read again, the failed one too.
    draft.setProgram('blastx');
    await draft.idle();
    expect(sourceOf(draft, 'query').base!.revisionId).toBe(revision);
    expect(sourceOf(draft, 'query').check).toEqual({ state: 'ok', records: 1 });
    expect(checks.at(-1)).toMatch(/^blastx\|query\|/);
    expect(sourceOf(draft, 'subject')).toMatchObject({ status: 'failed', parser: 2 });
    // BLASTP reads its query with the protein flags (kind 2): read again.
    draft.setProgram('blastp');
    await draft.idle();
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready', parser: 2 });
    expect(sourceOf(draft, 'query').base!.revisionId).not.toBe(revision);
    expect(sourceOf(draft, 'query').base!.records[0]!.length).toBe(4);
  });

  it('reads a text pasted just before a program change once, with the new kind', async () => {
    const { draft, data } = setup();
    const kinds: number[] = [];
    const indexSource = data.indexSource.bind(data);
    data.indexSource = async (sourceId, parser) => {
      kinds.push(parser);
      return indexSource(sourceId, parser);
    };
    draft.setPaste('query', '>a\nMKLV\n>b\nMKLV\n');
    draft.setProgram('blastp');
    await draft.idle();
    expect(kinds).toEqual([2]);
    expect(sourceOf(draft, 'query')).toMatchObject({ status: 'ready', parser: 2 });
  });
});

describe('SearchDraft accessions (line WF-1)', () => {
  const NR = '>NR_074540.1 Bacillus cereus 16S\nACGTACGT\n';
  const AB = '>AB000001.1 Rhizoctonia solani\nGGCCAATT\n';

  it('keeps a box of accessions without a source until it is fetched, then reads the FASTA written into it', async () => {
    const fetcher = new FakeAccessionFetcher({ 'NR_074540.1': NR, 'AB000001.1': AB });
    const { draft } = setup({ fetcher });
    draft.setPaste('query', 'nr_074540.1\n\nAB000001.1\n');
    await draft.idle();
    expect(draft.state.get().query.sources).toHaveLength(0);
    expect(draft.accessions('query')).toEqual(['NR_074540.1', 'AB000001.1']);
    expect(draft.readiness()).toBe('Fetch the query accessions first.');
    expect(fetcher.requests).toHaveLength(0);
    draft.fetchAccessions('query');
    await draft.idle();
    expect(fetcher.requests.map((r) => [r.source, r.kind, r.accession])).toEqual([
      ['ncbi', 'nucleotide', 'NR_074540.1'],
      ['ncbi', 'nucleotide', 'AB000001.1'],
    ]);
    expect(draft.state.get().query.paste).toBe(NR + AB);
    expect(draft.state.get().query.fetch).toBeUndefined();
    expect(sourceOf(draft, 'query').base?.records.map((r) => r.id)).toEqual(['NR_074540.1', 'AB000001.1']);
  });

  it('asks the chosen source for the kind the program reads, and keeps the choice', async () => {
    const fetcher = new FakeAccessionFetcher({ P12345: '>sp|P12345|AATM_RABIT\nMALLHS\n' });
    const { draft } = setup({ fetcher });
    draft.setProgram('blastp');
    draft.setFetchSource('ebi');
    draft.setPaste('subject', 'P12345');
    draft.fetchAccessions('subject');
    await draft.idle();
    expect(fetcher.requests.map((r) => [r.source, r.kind, r.accession])).toEqual([['ebi', 'protein', 'P12345']]);
    expect(draft.state.get().fetchSource).toBe('ebi');
    expect(draft.state.get().subject.paste).toBe('>sp|P12345|AATM_RABIT\nMALLHS\n');
  });

  it('says why a fetch failed, leaves the box as it was, and asks EBI nothing it cannot serve', async () => {
    const fetcher = new FakeAccessionFetcher({ 'AB000001.1': AB, 'ZZ000001.1': { ok: false, reason: 'refused', status: 429 } });
    const { draft } = setup({ fetcher });
    draft.setPaste('query', 'AB000001.1\nAB999999.1');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(draft.state.get().query.fetch).toEqual({ state: 'failed', message: 'NCBI has no nucleotide record AB999999.1.' });
    expect(draft.state.get().query.paste).toBe('AB000001.1\nAB999999.1');
    draft.setFetchSource('ebi');
    expect(draft.state.get().query.fetch).toBeUndefined();
    draft.setPaste('query', 'NR_074540.1');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(draft.state.get().query.fetch).toEqual({
      state: 'failed',
      message: 'EBI (ENA) has no RefSeq records: fetch NR_074540.1 from NCBI or CNCB.',
    });
    draft.setFetchSource('cncb');
    draft.setPaste('query', 'ZZ000001.1');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(draft.state.get().query.fetch).toEqual({
      state: 'failed',
      message: 'CNCB asks to wait (too many requests): fetch again in a moment.',
    });
    expect(fetcher.requests.map((r) => [r.source, r.accession])).toEqual([
      ['ncbi', 'AB000001.1'],
      ['ncbi', 'AB999999.1'],
      ['cncb', 'ZZ000001.1'],
    ]);
  });

  it('asks EBI nothing when one accession of the box is not for it', async () => {
    const fetcher = new FakeAccessionFetcher({ 'AB000001.1': AB });
    const { draft } = setup({ fetcher });
    draft.setFetchSource('ebi');
    draft.setPaste('query', 'AB000001.1\nNR_074540.1');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(fetcher.requests).toHaveLength(0);
    expect(draft.state.get().query.fetch).toMatchObject({ state: 'failed', message: expect.stringContaining('NR_074540.1') });
  });

  it('runs no search of accessions: Run neither fetches nor waits for a running fetch', async () => {
    const fetcher = new FakeAccessionFetcher({ 'NR_074540.1': NR });
    const { draft, requests } = setup({ fetcher });
    draft.setPaste('subject', '>s\nACGT\n');
    fetcher.hold();
    draft.setPaste('query', 'NR_074540.1');
    draft.fetchAccessions('query');
    expect(await draft.submit()).toEqual({ ok: false, message: 'The query accession is still being fetched.' });
    fetcher.release();
    await draft.idle();
    expect(requests).toHaveLength(0);
    expect(draft.state.get().query.paste).toBe(NR);
  });

  it('drops a running fetch when the box is edited, and limits the fetched text', async () => {
    const fetcher = new FakeAccessionFetcher({ 'NR_074540.1': NR, 'NC_000001.11': `>NC_000001.11\n${'A'.repeat(10_000_001)}\n` });
    const { draft } = setup({ fetcher });
    fetcher.hold();
    draft.setPaste('query', 'NR_074540.1');
    draft.fetchAccessions('query');
    await Promise.resolve();
    expect(draft.state.get().query.fetch).toEqual({ state: 'fetching', source: 'ncbi', done: 0, total: 1 });
    draft.setPaste('query', '>mine\nACGT\n');
    fetcher.release();
    await draft.idle();
    expect(draft.state.get().query.paste).toBe('>mine\nACGT\n');
    expect(draft.state.get().query.fetch).toBeUndefined();
    draft.setPaste('query', 'NC_000001.11');
    draft.fetchAccessions('query');
    await draft.idle();
    expect(draft.state.get().query.fetch).toMatchObject({ state: 'failed', message: expect.stringContaining('larger than 10 MB') });
  });
});

describe('SearchDraft submit', () => {
  it('refuses to queue without inputs, with every record excluded, or with a refused input', async () => {
    const { draft, requests } = setup();
    expect(await draft.submit()).toEqual({ ok: false, message: 'Add the query sequences: paste them or open FASTA files.' });
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT!ACGT\n');
    await draft.idle();
    const refused = await draft.submit();
    expect(refused.ok).toBe(false);
    expect(!refused.ok && refused.message).toMatch(/^Subject \(pasted\): subject record 1 \(s\) has/);
    draft.setIncluded('subject', 'paste', [0], false);
    expect(await draft.submit()).toEqual({ ok: false, message: 'Every subject record is excluded.' });
    expect(requests).toHaveLength(0);
    expect(draft.state.get().message).toEqual({ kind: 'error', text: 'Every subject record is excluded.' });
  });

  it("queues BLASTX with the regions of its nucleotide query and protein subject (DW-39)", async () => {
    const { draft, requests } = setup();
    draft.setPaste('query', '>q\nACGTACGTAC\n');
    draft.setPaste('subject', '>s\nMKLVMKLV\n');
    await draft.idle();
    draft.setProgram('blastx');
    await draft.idle();
    expect(draft.regionRecord('query')?.id).toBe('q');
    expect(draft.regionRecord('subject')?.id).toBe('s');
    draft.setRegion('query', { start: '2', stop: '8' });
    draft.setRegion('subject', { start: '3', stop: '8' });
    expect(draft.parameters()).toEqual([
      ['-query_loc', '2-8'],
      ['-subject_loc', '3-8'],
    ]);
    // The subject's region is checked against its residues, the query's against its bases.
    draft.setRegion('subject', { start: '3', stop: '9' });
    expect(await draft.submit()).toEqual({
      ok: false,
      message: 'Subject region: The stop must be between 1 and 8, the length of the record.',
    });
    draft.setRegion('subject', undefined);
    draft.setRegion('query', { start: '2', stop: '10' });
    expect(await draft.submit()).toMatchObject({ ok: true });
    expect(requests.at(-1)![0]!).toMatchObject({ program: 'blastx', parameters: [['-query_loc', '2-10']] });
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
    expect(draft.state.get().message).toBeUndefined();
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
    expect(draft.state.get().message).toBeUndefined();
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
    // The argv follows the sections of the page: Filters and Masking, then Discontiguous Word Options.
    expect(draft.parameters()).toEqual([
      ['-task', 'dc-megablast'],
      ['-lcase_masking', true],
      ['-template_length', '21'],
      ['-template_type', 'optimal'],
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

  it('restores the defaults of some fields: their values are cleared, the others kept, and the argv validated again', async () => {
    const { draft, validated } = setup();
    draft.setField('-task', 'dc-megablast');
    draft.setField('-evalue', '1e-5');
    draft.setField('-lcase_masking', true);
    draft.setField('-template_length', '18');
    draft.setProgram('tblastx');
    draft.setField('-evalue', '5');
    draft.setProgram('blastn');
    await draft.idle();
    draft.resetFields(['-evalue', '-lcase_masking', '-template_length', '-word_size']);
    expect(draft.state.get().values.blastn).toEqual({ '-task': 'dc-megablast' });
    expect(Object.isFrozen(draft.state.get().values.blastn)).toBe(true);
    // Another program's values stay.
    expect(draft.state.get().values.tblastx).toEqual({ '-evalue': '5' });
    await draft.idle();
    expect(validated.at(-1)).toEqual(['blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-task', 'dc-megablast']);
  });

  it('writes the words of Additional options after the form\'s, for each program, and refuses an unclosed quote', async () => {
    const { draft, requests, validated } = setup();
    draft.setField('-evalue', '1e-5');
    draft.setExtra(" -ungapped  -strand 'plus' ");
    await draft.idle();
    expect(validated.at(-1)).toEqual([
      'blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-evalue', '1e-5', '-ungapped', '-strand', 'plus',
    ]);
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT\n');
    expect((await draft.submit()).ok).toBe(true);
    expect(requests[0]![0]!.parameters).toEqual([['-evalue', '1e-5']]);
    expect(requests[0]![0]!.extra).toEqual(['-ungapped', '-strand', 'plus']);
    // Each program keeps its own words; a program without them adds none.
    draft.setProgram('blastp');
    await draft.idle();
    expect(draft.state.get().extra).toMatchObject({ blastn: " -ungapped  -strand 'plus' ", blastp: '' });
    expect(validated.at(-1)).toEqual(['blastp', '-query', 'query.fa', '-subject', 'subject.fa']);
    expect((await draft.submit()).ok).toBe(true);
    expect(requests[1]![0]!).not.toHaveProperty('extra');
    // An unclosed quote: the argv is not made, and nothing is queued.
    draft.setExtra("-seg '12 2.2");
    await draft.idle();
    expect(draft.state.get().validation).toEqual({ state: 'invalid', message: "Additional options: a quote (') is not closed" });
    expect(await draft.submit()).toEqual({ ok: false, message: "Additional options: a quote (') is not closed" });
    expect(requests).toHaveLength(2);
  });

  it('runs the words of Additional options in the coordinator; it refuses the options that LOSAT manages', async () => {
    const { draft, coordinator } = setup({ useCoordinator: true });
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s\nACGT\n');
    draft.setExtra('-num_threads 4');
    await draft.idle();
    expect(draft.state.get().validation).toEqual({ state: 'invalid', message: '-num_threads is managed by LOSAT' });
    expect(await draft.submit()).toEqual({ ok: false, message: '-num_threads is managed by LOSAT' });
    draft.setExtra('-qcov_hsp_perc 50 -ungapped');
    expect((await draft.submit()).ok).toBe(true);
    expect(coordinator.state.get().runs[0]!.snapshot.argv.slice(5)).toEqual(['-qcov_hsp_perc', '50', '-ungapped']);
  });

  it('refuses in Additional options the options that the form sets, and the subranges, with the control that sets them', async () => {
    const { draft } = setup();
    await draft.idle();
    const refusal = async (text: string) => {
      draft.setExtra(text);
      await draft.idle();
      const validation = draft.state.get().validation;
      return validation.state === 'invalid' ? validation.message : validation.state;
    };
    // A -task here would change the defaults that the form leaves out of the argv.
    expect(await refusal('-ungapped -task=blastn')).toBe('-task is set in the form: Optimize for');
    expect(await refusal('-evalue 1')).toBe('-evalue is set in the form: Expect threshold');
    expect(await refusal('-penalty -3')).toBe('-penalty is set in the form: Match/Mismatch Scores');
    expect(await refusal('-gapopen 5')).toBe('-gapopen is set in the form: Gap Costs');
    expect(await refusal("-dust '20 64 1'")).toBe('-dust is set in the form: Low complexity regions');
    expect(await refusal('-culling_limit 2')).toBe('-culling_limit is set in the form: Max matches in a query range');
    expect(await refusal('-query_loc 1-10')).toBe('-query_loc is set in the form: Query subrange');
    expect(await refusal('-subject_loc=1-10')).toBe('-subject_loc is set in the form: Subject subrange');
    // Another program's form: BLASTP has no -dust field, TBLASTX's genetic codes are fields.
    draft.setProgram('tblastx');
    expect(await refusal('-db_gencode 4')).toBe('-db_gencode is set in the form: Genetic code');
    draft.setProgram('blastp');
    expect(await refusal('-qcov_hsp_perc 50 -use_sw_tback')).toBe('ok');
  });

  it('gives every run of the draft its Job Title, outside the argv', async () => {
    const { draft, requests, coordinator } = setup({ useCoordinator: true });
    draft.addFiles('subject', [file('s1.fa', '>s1\nACGT\n'), file('s2.fa', '>s2\nACGT\n')]);
    draft.setPaste('query', '>q\nACGT\n');
    await draft.idle();
    draft.setMode('subject', 'separate');
    expect(draft.state.get().title).toBe('');
    draft.setTitle(' Two plasmids ');
    expect(draft.state.get().title).toBe(' Two plasmids ');
    await draft.submit();
    expect(requests[0]!.map((r) => r.title)).toEqual([' Two plasmids ', ' Two plasmids ']);
    expect(coordinator.state.get().runs.map((run) => run.snapshot.title)).toEqual(['Two plasmids', 'Two plasmids']);
    expect(coordinator.state.get().runs[0]!.snapshot.argv).not.toContain('Two plasmids');
    // An empty title is no title.
    draft.setTitle('   ');
    await draft.submit();
    expect(requests[1]!.every((r) => !('title' in r))).toBe(true);
  });

  it('asks the engine to prepare once the argv is valid and both roles have inputs, with the size and argv of the run', async () => {
    const { draft, prepared, validated } = setup();
    draft.setPaste('query', '>q1\nACGT\n');
    await draft.idle();
    expect(prepared).toEqual([]);
    draft.setPaste('subject', '>s1\nACGTACGT\n');
    await draft.idle();
    const argv = validated.at(-1);
    expect(argv?.slice(0, 5)).toEqual(['blastn', '-query', 'query.fa', '-subject', 'subject.fa']);
    expect(prepared.at(-1)).toEqual({ requestedThreads: 'auto', inputBytes: 9 + 13, argv });
    draft.setThreads(8);
    draft.setThreads(4);
    await draft.idle();
    // Stepping through the choices hints only the last one.
    expect(prepared.filter((hint) => hint.requestedThreads !== 'auto')).toEqual([{ requestedThreads: 4, inputBytes: 22, argv }]);
  });

  it('does not ask the engine to prepare for an argv that it refuses', async () => {
    const { draft, prepared } = setup({ validate: () => ({ ok: false, message: 'no' }) });
    draft.setPaste('query', '>q1\nACGT\n');
    draft.setPaste('subject', '>s1\nACGT\n');
    await draft.idle();
    draft.setThreads(2);
    await draft.idle();
    expect(prepared).toEqual([]);
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

  it('refuses to combine a source that starts without a defline after another source', async () => {
    const { draft } = setup();
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', 'ACGTACGT\n');
    draft.addFiles('subject', [new File(['>s\nACGT\n'], 's.fa')]);
    await draft.idle();
    // The pasted subject comes first: its sequence without a defline is the input's first record.
    expect(draft.readiness()).toBeUndefined();
    draft.setPaste('subject', '');
    draft.addFiles('subject', [new File(['GGCC\n'], 'bare.fa')]);
    await draft.idle();
    expect(draft.readiness()).toBe(
      'Subject (file bare.fa) starts with a sequence without a defline (">"), which cannot follow another input.',
    );
    draft.setMode('subject', 'separate');
    expect(draft.readiness()).toBeUndefined();
  });

  it("maps a subject record that the engine refuses, and reports it in short above the button", async () => {
    const { draft } = setup();
    draft.setPaste('query', '>q\nACGT\n');
    draft.setPaste('subject', '>s1\nACGT\n>s2\nACGT!ACGT\n');
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
