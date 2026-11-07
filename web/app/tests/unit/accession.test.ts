import { describe, expect, it } from 'vitest';
import { groupSite, SourceRace, type RaceStanding } from '../../src/application/source-race';
import {
  accessionLines,
  accessionsOnALine,
  fetchGroups,
  isSeqIdReason,
  pointsToFetch,
  RACE_RECORDS,
  servingSources,
  sourceList,
  unserved,
} from '../../src/domain/accession';
import { BrowserAccessionFetcher, fetchUrl } from '../../src/infra/browser/accession-fetcher';
import { FakeAccessionFetcher } from '../../src/infra/fake/fake-fetcher';

describe('accessionLines', () => {
  it('takes INSDC, RefSeq, GenBase and UniProt accessions, one per line, in upper case', () => {
    expect(accessionLines(' NR_074540.1 \n\nab000001\nAAA12345.2\nAAAA01000001\nNZ_CP006932.1\nC_AA004835.1\nP12345\nQ9Y261\nA0A023GPI8\n')).toEqual([
      'NR_074540.1',
      'AB000001',
      'AAA12345.2',
      'AAAA01000001',
      'NZ_CP006932.1',
      'C_AA004835.1',
      'P12345',
      'Q9Y261',
      'A0A023GPI8',
    ]);
  });

  it('takes nothing else: FASTA, a bare sequence, a gi, an empty box, a mix', () => {
    for (const text of ['>NR_074540.1\nACGT', 'ACGTACGT', 'MKVLAAGIV', '12345678', '', '\n', 'NR_074540.1\nACGT']) {
      expect(accessionLines(text)).toBeUndefined();
    }
  });

  it('finds accessions that are not one per line (WA2-OP-04), as typed', () => {
    expect(accessionsOnALine('AB000001.1 MH011443.1')).toEqual(['AB000001.1', 'MH011443.1']);
    expect(accessionsOnALine('ab000001.1,MH011443.1\n')).toEqual(['ab000001.1', 'MH011443.1']);
    expect(accessionsOnALine(' NR_074540.1, \tP12345 ,\r\nQ9Y261\n\nAB000001.1,')).toEqual(['NR_074540.1', 'P12345', 'Q9Y261', 'AB000001.1']);
    // No version (a sequence to the reader): accessions still.
    expect(accessionsOnALine('AB000001 MH011443')).toEqual(['AB000001', 'MH011443']);
    expect(accessionsOnALine('P12345 Q9Y261')).toEqual(['P12345', 'Q9Y261']);
    // A sequence, a GenBank ORIGIN line and a protein sequence have no accession-shaped word.
    for (const text of ['ACGTACGTAC GGTTACGTAA', '        1 gatcctccat atacaacggt', '       61 aaccggttaa cc', 'MKLVFFAEDV GSNKGAIIGL']) {
      expect(accessionsOnALine(text), text).toBeUndefined();
    }
    // One per line is accessionLines'; anything but accessions (a gi, a sequence, a defline) is FASTA's.
    for (const text of ['AB000001.1\nMH011443.1', 'AB000001.1', '', ' , ', 'AB000001.1 12345678', 'AB000001.1 ACGT', '>q AB000001.1,MH011443.1']) {
      expect(accessionsOnALine(text), text).toBeUndefined();
    }
  });

  it('points to Fetch only for an accession on the first line of a Seq-id reason (WA2-OP-04)', () => {
    const reason = (line: string) =>
      `the first line ("${line}") is not a defline and may be a sequence identifier that NCBI BLAST+ fetches through a data loader (from GenBank or a BLAST database), which is not supported in the browser (start the input with a '>' defline)`;
    const accessions = ['AB000001.1', 'ab000001', 'NR_074540.1', 'P12345', 'Q9Y261', 'C_AA004835.1', 'AAA12345'];
    for (const line of accessions) {
      expect(pointsToFetch(reason(line), `${line}\nACGT\n`), line).toBe(true);
      expect(pointsToFetch(reason(line), `\n  ${line} \r\n>x\nACGT`), line).toBe(true);
    }
    // A GI number, a PDB id and a tagged id are Seq-ids to the engine, which Fetch does not take.
    for (const line of ['12345678', '4HHB', '1ABC_A', 'sp|P04637|P53_HUMAN', 'ref|NM_000546.6|', 'gi|12345', 'AB000001.1 ACGT']) {
      expect(pointsToFetch(reason(line), `${line}\nACGT\n`), line).toBe(false);
    }
    // Another reason, a reader's error, or no reason: no sentence, whatever the first line is.
    expect(pointsToFetch("CFastaReader: Near line 2, there's a line that doesn't look like plausible data", 'AB000001.1\nACGT')).toBe(false);
    expect(pointsToFetch(undefined, 'AB000001.1\nACGT')).toBe(false);
    expect(pointsToFetch(reason('AB000001.1'), '')).toBe(false);
    expect(isSeqIdReason(reason('x'))).toBe(true);
    expect(isSeqIdReason('line 3 is a gap line')).toBe(false);
  });

  it('knows what EBI cannot serve without asking it', () => {
    expect(unserved('ebi', 'nucleotide', 'NR_074540.1')).toBe('EBI (ENA) has no RefSeq records');
    expect(unserved('ebi', 'nucleotide', 'AB000001.1')).toBeUndefined();
    expect(unserved('ebi', 'protein', 'WP_000001.1')).toBe('EBI (UniProt) has UniProt accessions only');
    expect(unserved('ebi', 'protein', 'P12345')).toBeUndefined();
    expect(unserved('ncbi', 'protein', 'WP_000001.1')).toBeUndefined();
    expect(servingSources('nucleotide', ['AB000001.1', 'NR_074540.1'])).toEqual(['ncbi', 'cncb']);
    expect(servingSources('protein', ['WP_000001.1'])).toEqual(['ncbi', 'cncb']);
  });

  it('knows that CNCB has no UniProt records and that only CNCB has GenBase records (2026-11-03)', () => {
    expect(unserved('cncb', 'protein', 'P12345')).toBe('CNCB has no UniProt records');
    expect(unserved('cncb', 'protein', 'NP_414542.1')).toBeUndefined();
    expect(servingSources('protein', ['P12345'])).toEqual(['ncbi', 'ebi']);
    expect(unserved('ncbi', 'nucleotide', 'C_AA004835.1')).toBe('NCBI has no GenBase records');
    expect(unserved('ebi', 'nucleotide', 'C_AA004835.1')).toBe('EBI (ENA) has no GenBase records');
    expect(servingSources('nucleotide', ['C_AA004835.1'])).toEqual(['cncb']);
    expect(servingSources('protein', ['C_AA004835.1', 'P12345'])).toEqual([]);
    expect(sourceList(['cncb'])).toBe('CNCB');
    expect(sourceList(['ncbi', 'cncb'])).toBe('NCBI or CNCB');
    expect(sourceList(['ncbi', 'ebi', 'cncb'], 'and')).toBe('NCBI, EBI and CNCB');
  });

  it('groups a box by the sites that serve each accession, in the order of the box (line WF-5)', () => {
    expect(fetchGroups('protein', ['P12345', 'C_AA004835.1', 'Q9Y261', 'C_AA000002.1'])).toEqual([
      { sources: ['ncbi', 'ebi'], accessions: ['P12345', 'Q9Y261'], lines: [0, 2] },
      { sources: ['cncb'], accessions: ['C_AA004835.1', 'C_AA000002.1'], lines: [1, 3] },
    ]);
    // RefSeq proteins are neither UniProt's nor only GenBase's: a third group.
    expect(fetchGroups('protein', ['C_AA004835.1', 'NP_414542.1', 'P12345', 'WP_000001.1'])).toEqual([
      { sources: ['cncb'], accessions: ['C_AA004835.1'], lines: [0] },
      { sources: ['ncbi', 'cncb'], accessions: ['NP_414542.1', 'WP_000001.1'], lines: [1, 3] },
      { sources: ['ncbi', 'ebi'], accessions: ['P12345'], lines: [2] },
    ]);
    // CNCB serves every nucleotide accession: a nucleotide box is never one that no one site serves.
    expect(servingSources('nucleotide', ['C_AA004835.1', 'NR_074540.1', 'AB000001.1', 'P12345'])).toEqual(['cncb']);
  });

  it('starts a group at the site chosen by hand if it serves the group, else at the race, else at the first site', () => {
    const standing: RaceStanding = { answered: [{ source: 'ebi', ms: 900 }], failed: [] };
    expect(groupSite(['ncbi', 'ebi'], 'fastest', undefined)).toBe('ncbi');
    expect(groupSite(['ncbi', 'ebi'], 'fastest', standing)).toBe('ebi');
    // A race still running that has not reached the group's sites: the first of them.
    expect(groupSite(['ncbi', 'cncb'], 'fastest', standing)).toBe('ncbi');
    expect(groupSite(['ncbi', 'ebi'], 'ncbi', standing)).toBe('ncbi');
    expect(groupSite(['ncbi', 'ebi'], 'cncb', standing)).toBe('ebi');
    expect(groupSite(['cncb'], 'ebi', standing)).toBe('cncb');
    // Not a site that failed the race while another of the group's has yet to answer.
    expect(groupSite(['ncbi', 'ebi'], 'fastest', { answered: [{ source: 'cncb', ms: 300 }], failed: ['ncbi'] })).toBe('ebi');
    expect(groupSite(['ncbi', 'ebi'], 'fastest', { answered: [{ source: 'cncb', ms: 300 }], failed: ['ncbi', 'ebi'] })).toBe('ncbi');
  });

  it('races with records that every source serves, of the kind raced', () => {
    for (const kind of ['nucleotide', 'protein'] as const) {
      for (const [source, accession] of Object.entries(RACE_RECORDS[kind])) {
        expect(accessionLines(accession), accession).toEqual([accession]);
        expect(unserved(source as keyof (typeof RACE_RECORDS)[typeof kind], kind, accession)).toBeUndefined();
      }
    }
  });
});

describe('SourceRace', () => {
  it('times each source from its request to its answer, and orders them by arrival', async () => {
    let clock = 0;
    const fetcher = new FakeAccessionFetcher({ 'MH011443.1': '>MH011443.1\nACGT\n' });
    const changes: Array<RaceStanding | undefined> = [];
    const race = new SourceRace({ fetcher, now: () => clock, onChange: (_, standing) => changes.push(standing) });
    fetcher.hold();
    const first = race.fastest('nucleotide', ['ncbi', 'cncb'], new AbortController().signal);
    clock = 300;
    fetcher.release('ebi');
    await new Promise((resolve) => setTimeout(resolve, 0));
    clock = 800;
    fetcher.release('cncb');
    expect(await first).toBe('cncb');
    clock = 1500;
    fetcher.release('ncbi');
    await race.idle();
    expect(changes.at(-1)).toEqual({
      answered: [
        { source: 'ebi', ms: 300 },
        { source: 'cncb', ms: 800 },
        { source: 'ncbi', ms: 1500 },
      ],
      failed: [],
    });
    // Kept: no second race.
    expect(await race.fastest('nucleotide', ['ncbi', 'ebi', 'cncb'], new AbortController().signal)).toBe('ebi');
    expect(fetcher.requests).toHaveLength(3);
  });

  it('rejects a waiting fetch that aborts, and goes on without it', async () => {
    const fetcher = new FakeAccessionFetcher({ 'MH011443.1': '>MH011443.1\nACGT\n' });
    const race = new SourceRace({ fetcher, now: () => 0, onChange: () => undefined });
    fetcher.hold();
    const controller = new AbortController();
    const waiting = race.fastest('nucleotide', ['ncbi'], controller.signal);
    const other = race.fastest('nucleotide', ['ncbi', 'cncb'], new AbortController().signal);
    controller.abort(new DOMException('edited', 'AbortError'));
    await expect(waiting).rejects.toThrow('edited');
    fetcher.release('cncb');
    expect(await other).toBe('cncb');
    fetcher.release();
    await race.idle();
    expect(await race.fastest('nucleotide', ['ncbi'], new AbortController().signal)).toBe('ncbi');
    expect(fetcher.requests).toHaveLength(3);
  });

  it('times from the moment the request goes out, not from the wait between NCBI requests', async () => {
    let clock = 0;
    const changes: Array<RaceStanding | undefined> = [];
    const fetcher = {
      // Each request waits until 1000 before it goes out, then is answered at once.
      fetch: async (request: { readonly onStart?: () => void }) => {
        clock = 1000;
        request.onStart?.();
        return { ok: true as const, fasta: '>x\nA\n' };
      },
    };
    const race = new SourceRace({ fetcher, now: () => clock, onChange: (_, standing) => changes.push(standing) });
    await race.fastest('nucleotide', ['ncbi'], new AbortController().signal);
    await race.idle();
    expect(changes.at(-1)?.answered.map((entry) => entry.ms)).toEqual([0, 0, 0]);
  });

  it('counts a source that has not answered by the deadline as failed, and races again after a race without answers', async () => {
    const fetcher = new FakeAccessionFetcher({ 'MH011443.1': '>MH011443.1\nACGT\n' });
    const changes: Array<RaceStanding | undefined> = [];
    const race = new SourceRace({ fetcher, now: () => 0, deadlineMs: 20, onChange: (_, standing) => changes.push(standing) });
    fetcher.hold();
    // Every eligible source failed: the first of them.
    expect(await race.fastest('nucleotide', ['cncb', 'ncbi'], new AbortController().signal)).toBe('cncb');
    await race.idle();
    expect(changes.at(-1)).toBeUndefined();
    fetcher.release();
    expect(await race.fastest('nucleotide', ['ncbi'], new AbortController().signal)).toBe('ncbi');
    expect(fetcher.requests).toHaveLength(6);
  });
});

describe('BrowserAccessionFetcher', () => {
  const answer = (status: number, body: string) => async () => new Response(body, { status });
  const request = { source: 'ncbi', kind: 'nucleotide', accession: 'NR_074540.1', maxBytes: 1000 } as const;

  it('asks each source for the accession only', () => {
    expect(fetchUrl('ncbi', 'nucleotide', 'NR_074540.1')).toBe(
      'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_074540.1&rettype=fasta&retmode=text',
    );
    expect(fetchUrl('ncbi', 'protein', 'P12345')).toContain('db=protein&id=P12345&');
    expect(fetchUrl('ebi', 'nucleotide', 'AB000001.1')).toBe('https://www.ebi.ac.uk/ena/browser/api/fasta/AB000001.1');
    expect(fetchUrl('ebi', 'protein', 'P12345')).toBe('https://rest.uniprot.org/uniprotkb/P12345.fasta');
    expect(fetchUrl('cncb', 'nucleotide', 'C_AA004835.1')).toBe('https://ngdc.cncb.ac.cn/genbase/api/file/fasta?acc=C_AA004835.1');
  });

  it("reads each source's way of saying that it has no such record", async () => {
    const signal = new AbortController().signal;
    for (const [status, body] of [
      [200, 'Error: F a i l e d  t o  u n d e r s t a n d  i d :  Z Z 9 9\n'],
      [200, ''],
      [400, 'Error'],
      [404, 'not found'],
    ] as const) {
      const fetcher = new BrowserAccessionFetcher({ fetch: answer(status, body) });
      expect(await fetcher.fetch(request, signal)).toEqual({ ok: false, reason: 'not-found' });
    }
    const fetcher = new BrowserAccessionFetcher({ fetch: answer(429, '') });
    expect(await fetcher.fetch(request, signal)).toEqual({ ok: false, reason: 'refused', status: 429 });
  });

  it("asks past the browser's HTTP cache only for the race, without cookies or referrer", async () => {
    const inits: RequestInit[] = [];
    const fetcher = new BrowserAccessionFetcher({
      fetch: async (_url, init) => {
        inits.push(init ?? {});
        return new Response('>x\nA\n');
      },
    });
    const signal = new AbortController().signal;
    await fetcher.fetch({ ...request, source: 'cncb' }, signal);
    await fetcher.fetch({ ...request, source: 'cncb', fresh: true }, signal);
    expect(inits.map((init) => [init.cache, init.credentials, init.referrerPolicy])).toEqual([
      [undefined, 'omit', 'no-referrer'],
      ['no-store', 'omit', 'no-referrer'],
    ]);
  });

  it('returns the FASTA with a final newline, and drops an answer over the limit', async () => {
    const signal = new AbortController().signal;
    const ok = new BrowserAccessionFetcher({ fetch: answer(200, '>NR_074540.1 x\nACGT') });
    expect(await ok.fetch(request, signal)).toEqual({ ok: true, fasta: '>NR_074540.1 x\nACGT\n' });
    const big = new BrowserAccessionFetcher({ fetch: answer(200, `>x\n${'A'.repeat(2000)}\n`) });
    expect(await big.fetch(request, signal)).toEqual({ ok: false, reason: 'too-large' });
  });

  it('reports a network failure as unreachable, and rejects when the user aborts', async () => {
    const failing = new BrowserAccessionFetcher({
      fetch: async () => {
        throw new TypeError('Failed to fetch');
      },
    });
    expect(await failing.fetch(request, new AbortController().signal)).toEqual({
      ok: false,
      reason: 'unreachable',
      detail: 'Failed to fetch',
    });
    const controller = new AbortController();
    controller.abort();
    await expect(failing.fetch(request, controller.signal)).rejects.toThrow();
  });

  it('gives up when no answer comes in time', async () => {
    const fetcher = new BrowserAccessionFetcher({
      timeoutMs: 10,
      fetch: (_url, init) =>
        new Promise((_resolve, reject) => init?.signal?.addEventListener('abort', () => reject(init.signal!.reason))),
    });
    expect(await fetcher.fetch({ ...request, source: 'cncb' }, new AbortController().signal)).toEqual({
      ok: false,
      reason: 'unreachable',
      detail: 'no answer in time',
    });
  });

  it('starts NCBI requests at least 340 ms apart', async () => {
    let clock = 0;
    const starts: number[] = [];
    const fetcher = new BrowserAccessionFetcher({
      now: () => clock,
      fetch: async () => {
        starts.push(clock);
        return new Response('>x\nA\n');
      },
    });
    const signal = new AbortController().signal;
    await fetcher.fetch(request, signal);
    clock = 100;
    const started = performance.now();
    await fetcher.fetch(request, signal);
    expect(performance.now() - started).toBeGreaterThanOrEqual(230);
    expect(starts).toEqual([0, 100]);
  });
});
