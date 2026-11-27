import { describe, expect, it, vi } from 'vitest';
import { parseToolDescription } from 'losat/domain';
import type { ToolEngine } from 'losat/ports';
import { BlastFormatterDraft } from '../../src/application/blast-formatter-draft';
import type { EnqueueAllResult, ToolRequest } from '../../src/application/coordinator';
import { buildBlastFormatterArgv, changedBlastFormatterOptionsCount, isCountsOutputFile } from '../../src/domain/tool-argv';
import {
  BLAST_FORMATTER_DEFAULT_OUTFMT,
  BLAST_FORMATTER_OUTFMT_CHOICES,
  blastFormatterOptionParameters,
  toolLabel,
  toolStageLabel,
} from '../../src/domain/tools';
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

  it('offers the formats that blast_formatter writes, 0 (its default) first', () => {
    // Each PR that makes blast_formatter write another format adds it to the choices; the e2e test runs each on the engine.
    expect(BLAST_FORMATTER_OUTFMT_CHOICES.map((choice) => choice.value)).toEqual(['0', '1', '2', '3', '4', '8', '9', '12']);
    expect(BLAST_FORMATTER_DEFAULT_OUTFMT).toBe('0');
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
  function setup(options: { readonly busy?: boolean } = {}) {
    const fake = new FakeEngine();
    const validateTool = vi.fn((...args: Parameters<ToolEngine['validateTool']>) => fake.validateTool(args[0]));
    const engine = { describeTool: (tool: string) => fake.describeTool(tool), validateTool };
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
      isBusy: () => options.busy === true,
      debounceMs: 0,
    });
    return { engine, draft, queued, enqueueToolAll, getSessionArchiveFile, validateTool };
  }

  async function waitDescribed(draft: BlastFormatterDraft): Promise<void> {
    await vi.waitFor(() => expect(draft.state.get().description).toBeDefined());
  }

  it('loads -xmlhelp, submits a picked archive file with run directory files, and tracks leave guard', async () => {
    const { draft, queued } = setup();
    expect(draft.holdsInputs()).toBe(false);
    draft.ensureLoaded();
    await waitDescribed(draft);

    const archive = new File(['Blast4-archive ::= {}'], 'db.asn');
    const dbFile = new File(['dbbytes'], 's3.nin');
    Object.defineProperty(dbFile, 'webkitRelativePath', { value: 'in/db/s3.nin' });
    draft.setArchiveFile(archive);
    draft.addRunFiles([dbFile]);
    draft.setOutfmt('12');
    draft.setOutName('report.json');
    draft.setField('-max_target_seqs', '2');

    expect(draft.holdsInputs()).toBe(true);
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
    await waitDescribed(draft);

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
  it('refuses a second file under a name already in the run directory, as WindowMasker does', async () => {
    const { draft } = setup();
    draft.addRunFiles([new File(['a'], 'app.ini')]);
    draft.addRunFiles([new File(['bb'], 'app.ini')]);
    expect(draft.state.get().runFiles.map((entry) => entry.name)).toEqual(['app.ini']);
    expect(draft.state.get().runFilesError).toBe('A file named "app.ini" is already in the run directory.');
  });

  it('checks the archive and the files of the run as one run directory', async () => {
    const { draft, validateTool } = setup();
    draft.ensureLoaded();
    await waitDescribed(draft);
    draft.setArchiveFile(new File(['Blast4-archive ::= {}'], 'bn.asn'));
    draft.addRunFiles([new File(['other'], 'bn.asn')]);
    await vi.waitFor(() => expect(draft.state.get().validation.state).toBe('invalid'));
    expect(draft.state.get().validation).toEqual({
      state: 'invalid',
      message: 'A file named "bn.asn" is already in the run directory.',
    });
    const res = await draft.submit();
    expect(res.ok).toBe(false);
    // The check refused the names before asking the engine; the run directory of a valid form names the archive once.
    draft.removeRunFile(draft.state.get().runFiles[0]!.key);
    await vi.waitFor(() => expect(draft.state.get().validation.state).toBe('ok'));
    const last = validateTool.mock.calls.at(-1)!;
    expect(last[0]).toEqual(['blast_formatter', '-archive', 'bn.asn']);
    expect(last[1]?.directory?.map((entry) => entry.name)).toEqual(['bn.asn']);
  });

  it('runs words that name -archive themselves from the run directory, as the CLI does', async () => {
    const { draft, queued, validateTool } = setup();
    draft.ensureLoaded();
    await waitDescribed(draft);
    draft.addRunFiles([new File(['Blast4-archive ::= {}'], 'res.asn')]);
    draft.setExtra('-archive res.asn');
    await vi.waitFor(() => expect(draft.state.get().validation.state).toBe('ok'));
    expect(validateTool.mock.calls.at(-1)![0]).toEqual(['blast_formatter', '-archive', 'res.asn']);

    const res = await draft.submit();
    expect(res.ok).toBe(true);
    expect(queued).toHaveLength(1);
    expect(queued[0]!.argv).toEqual(['blast_formatter', '-archive', 'res.asn']);
    // No archive of the source: the run directory is the files of the run, and the engine opens -archive there.
    expect(queued[0]!.directory.map((entry) => entry.name)).toEqual(['res.asn']);
  });

  it('checks the argv while the engine is busy, as WindowMasker does', async () => {
    const { draft, validateTool } = setup({ busy: true });
    draft.ensureLoaded();
    await waitDescribed(draft);
    draft.setArchivePaste('Blast4-archive ::= {}');
    draft.setExtra('-max_target_seqs');
    await vi.waitFor(() => expect(draft.state.get().validation.state).toBe('invalid'));
    expect(validateTool).toHaveBeenCalled();
    draft.setExtra('');
    await vi.waitFor(() => expect(draft.state.get().validation.state).toBe('ok'));
  });

  it('shows nothing to check before an archive is chosen, and the CLI\'s -rid refusal once -rid is typed', async () => {
    const { draft, validateTool } = setup();
    draft.ensureLoaded();
    await waitDescribed(draft);
    await vi.waitFor(() => expect(draft.state.get().validation.state).toBe('idle'));
    expect(validateTool).not.toHaveBeenCalled();
    draft.setExtra('-rid ABCDEFGH01R');
    await vi.waitFor(() => expect(draft.state.get().validation.state).toBe('invalid'));
    expect(validateTool.mock.calls.at(-1)![0]).toEqual(['blast_formatter', '-rid', 'ABCDEFGH01R']);
  });

  it('reports a quote left open in Additional options once', async () => {
    const { draft } = setup();
    draft.ensureLoaded();
    await waitDescribed(draft);
    draft.setArchivePaste('Blast4-archive ::= {}');
    draft.setExtra("-logfile 'a.log");
    await vi.waitFor(() => expect(draft.state.get().validation.state).toBe('invalid'));
    expect(draft.state.get().validation).toEqual({ state: 'invalid', message: "Additional options: a quote (') is not closed" });
  });
});
