import { describe, expect, it, vi } from 'vitest';
import { parseToolDescription } from 'losat/domain';
import { BlastFormatterDraft, holdsBlastFormatterInputs } from '../../src/application/blast-formatter-draft';
import type { EnqueueAllResult, ToolRequest } from '../../src/application/coordinator';
import { buildBlastFormatterArgv, changedBlastFormatterOptionsCount, isCountsOutputFile } from '../../src/domain/tool-argv';
import { blastFormatterOptionParameters, toolLabel, toolStageLabel } from '../../src/domain/tools';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import BLAST_FORMATTER_XMLHELP_DOC from '../../src/infra/fake/blast_formatter-xmlhelp.json';

const BF_DESC = parseToolDescription(BLAST_FORMATTER_XMLHELP_DOC.xml);

describe('blast_formatter domain helpers', () => {
  it('lists formatting options in NCBI usage order without block flags', () => {
    const flags = blastFormatterOptionParameters(BF_DESC).map((p) => p.flag);
    expect(flags).toEqual([
      '-show_gis',
      '-num_descriptions',
      '-num_alignments',
      '-line_length',
      '-html',
      '-sorthits',
      '-sorthsps',
      '-max_target_seqs',
      '-parse_deflines',
    ]);
  });

  it('builds blast_formatter argv omitting default -outfmt 0 and -out -', () => {
    expect(
      buildBlastFormatterArgv({
        archiveName: 'bn.asn',
        outfmt: '8',
        outName: '',
        values: {},
        description: BF_DESC,
      }),
    ).toEqual(['blast_formatter', '-archive', 'bn.asn', '-outfmt', '8']);

    expect(
      buildBlastFormatterArgv({
        archiveName: 'bn.asn',
        outfmt: '0',
        outName: '-',
        values: { '-max_target_seqs': '5', '-parse_deflines': true },
        description: BF_DESC,
        extra: ['-logfile', 'run.log'],
      }),
    ).toEqual(['blast_formatter', '-archive', 'bn.asn', '-max_target_seqs', '5', '-parse_deflines', '-logfile', 'run.log']);
  });

  it('counts changed formatting options and identifies labels', () => {
    expect(changedBlastFormatterOptionsCount({}, '', BF_DESC)).toBe(0);
    expect(changedBlastFormatterOptionsCount({ '-max_target_seqs': '10', '-html': true }, '-logfile a.log', BF_DESC)).toBe(3);
    expect(toolLabel('blast_formatter')).toBe('Format Archive');
    expect(toolStageLabel('blast_formatter', 'format')).toBe('blast_formatter');
    expect(isCountsOutputFile('format', ['blast_formatter', '-out', 'x'], 'x', 'blast_formatter')).toBe(false);
  });
});

describe('BlastFormatterDraft', () => {
  function setup() {
    const engine = new FakeEngine();
    const queued: ToolRequest[] = [];
    const enqueueToolAll = vi.fn(async (requests: readonly ToolRequest[]): Promise<EnqueueAllResult> => {
      queued.push(...requests);
      return { ok: true, runIds: requests.map((_, idx) => `run-${idx + 1}`) };
    });
    const getSessionArchiveFile = vi.fn(async (_runId: string, name: string) => new File(['Blast4-archive ::= {}'], name));
    const draft = new BlastFormatterDraft({
      engine,
      enqueueToolAll,
      getSessionArchiveFile,
      debounceMs: 0,
    });
    return { engine, draft, queued, enqueueToolAll, getSessionArchiveFile };
  }

  it('loads -xmlhelp, submits a picked archive file with run directory files, and tracks leave guard', async () => {
    const { draft, queued } = setup();
    expect(holdsBlastFormatterInputs(draft.state.get())).toBe(false);
    draft.ensureLoaded();
    await vi.waitFor(() => expect(draft.state.get().description).toBeDefined());

    const archive = new File(['Blast4-archive ::= {}'], 'db.asn');
    const dbFile = new File(['dbbytes'], 's3.nin');
    Object.defineProperty(dbFile, 'webkitRelativePath', { value: 'in/db/s3.nin' });
    draft.setArchiveFile(archive);
    draft.addRunFiles([dbFile]);
    draft.setOutfmt('12');
    draft.setOutName('report.json');
    draft.setField('-max_target_seqs', '2');

    expect(holdsBlastFormatterInputs(draft.state.get())).toBe(true);
    expect(draft.previewCommandWords()).toEqual([
      'losat',
      'blast_formatter',
      '-archive',
      'db.asn',
      '-outfmt',
      '12',
      '-out',
      'report.json',
      '-max_target_seqs',
      '2',
    ]);

    const res = await draft.submit();
    expect(res.ok).toBe(true);
    expect(queued).toHaveLength(1);
    expect(queued[0]!.tool).toBe('blast_formatter');
    expect(queued[0]!.stage).toBe('format');
    expect(queued[0]!.argv).toEqual([
      'blast_formatter',
      '-archive',
      'db.asn',
      '-outfmt',
      '12',
      '-out',
      'report.json',
      '-max_target_seqs',
      '2',
    ]);
    expect(queued[0]!.directory.map((d) => d.name)).toEqual(['db.asn', 'in/db/s3.nin']);
  });

  it('supports session -outfmt 11 shortcut and refuses -rid via validateTool', async () => {
    const { draft, queued, getSessionArchiveFile } = setup();
    draft.useSessionArchive({ runId: 'r1', runNumber: 1, program: 'blastn', name: 'run_1_blastn.asn' });
    await vi.waitFor(() => expect(draft.state.get().description).toBeDefined());

    const res = await draft.submit();
    expect(res.ok).toBe(true);
    expect(getSessionArchiveFile).toHaveBeenCalledWith('r1', 'run_1_blastn.asn');
    expect(queued[0]!.directory[0]).toMatchObject({
      name: 'run_1_blastn.asn',
      origin: { kind: 'archive', runNumber: 1 },
    });

    // -rid in Additional options is refused.
    draft.setArchiveSource('paste');
    draft.setArchivePaste('');
    draft.setExtra('-rid ABCDEFGH01R');
    const bad = await draft.submit();
    expect(bad.ok).toBe(false);
    expect(!bad.ok && bad.message).toMatch(/BLAST_FORMATTER|Choose|Paste/);
  });
});
