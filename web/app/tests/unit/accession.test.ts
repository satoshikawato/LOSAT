import { describe, expect, it } from 'vitest';
import { accessionLines, unserved } from '../../src/domain/accession';
import { BrowserAccessionFetcher, fetchUrl } from '../../src/infra/browser/accession-fetcher';

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

  it('knows what EBI cannot serve without asking it', () => {
    expect(unserved('ebi', 'nucleotide', 'NR_074540.1')).toBe('EBI (ENA) has no RefSeq records');
    expect(unserved('ebi', 'nucleotide', 'AB000001.1')).toBeUndefined();
    expect(unserved('ebi', 'protein', 'WP_000001.1')).toBe('EBI (UniProt) has UniProt accessions only');
    expect(unserved('ebi', 'protein', 'P12345')).toBeUndefined();
    expect(unserved('ncbi', 'protein', 'WP_000001.1')).toBeUndefined();
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
