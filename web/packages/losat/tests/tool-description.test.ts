// Tests for `parseToolDescription` over every XML file of LOSAT/src/ncbi_xmlhelp/ and for
// `describeTool` through `reactorControl`, `WasmEngine` and the serial reactor in Node.
import { readdirSync, readFileSync } from 'node:fs';
import { join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { decodeXmlEntities, parseToolDescription } from '../src/domain/tool-description.ts';
import { WasmEngine } from '../src/infra/engine-worker/gateway.ts';
import type { EngineCommand, EngineEvent, WorkerTool } from '../src/infra/engine-worker/protocol.ts';
import { reactorControl } from '../src/infra/reactor/control.ts';
import { compileReactor, instantiateSerial } from '../src/infra/reactor/instance.ts';

const REPOSITORY = fileURLToPath(new URL('../../../..', import.meta.url));
const XMLHELP_DIR = join(REPOSITORY, 'LOSAT/src/ncbi_xmlhelp');
const REACTORS = process.env['LOSAT_WEB_REACTORS'];

const xmlFiles = readdirSync(XMLHELP_DIR)
  .filter((name) => name.endsWith('.xml'))
  .sort();

describe('parseToolDescription', () => {
  it('parses every XML file in LOSAT/src/ncbi_xmlhelp/', () => {
    expect(xmlFiles.length).toBeGreaterThanOrEqual(17);
    for (const file of xmlFiles) {
      const xml = readFileSync(join(XMLHELP_DIR, file), 'utf8');
      const parsed = parseToolDescription(xml);
      const expectedName = file.slice(0, -'.xml'.length);
      expect(parsed.tool, file).toBe(expectedName);
      expect(parsed.parameters.length, file).toBeGreaterThan(0);
      for (const param of parsed.parameters) {
        expect(param.flag, `${file}:${param.flag}`).toMatch(/^-/);
        expect(typeof param.help, `${file}:${param.flag}`).toBe('string');
      }
      for (const [a, b] of [...parsed.excludes, ...parsed.requires]) {
        expect(a, file).toMatch(/^-/);
        expect(b, file).toMatch(/^-/);
      }
    }
  });

  it('extracts every parameter, default, value list, range and exclusion from windowmasker.xml', () => {
    const parsed = parseToolDescription(readFileSync(join(XMLHELP_DIR, 'windowmasker.xml'), 'utf8'));
    expect(parsed.tool).toBe('windowmasker');
    expect(parsed.version).toBe('1.0.0');
    expect(parsed.description).toBe('Window based sequence masker');
    expect(parsed.parameters).toHaveLength(38);
    expect(parsed.parameters.filter((p) => p.takesValue)).toHaveLength(29);
    expect(parsed.parameters.filter((p) => !p.takesValue)).toHaveLength(9);
    expect(parsed.excludes).toHaveLength(76);
    expect(parsed.requires).toHaveLength(0);

    const byFlag = new Map(parsed.parameters.map((p) => [p.flag, p]));

    expect(byFlag.get('-checkdup')).toMatchObject({
      kind: 'Boolean',
      takesValue: true,
      mandatory: false,
      defaultValue: 'false',
      group: '1',
    });
    expect(byFlag.get('-dust')).toMatchObject({ kind: 'Boolean', defaultValue: 'F' });
    expect(byFlag.get('-text_match')).toMatchObject({ kind: 'Boolean', defaultValue: 'T' });
    expect(byFlag.get('-in')).toMatchObject({ kind: 'File_In', defaultValue: '-' });
    expect(byFlag.get('-out')).toMatchObject({ kind: 'File_Out', defaultValue: '-' });
    expect(byFlag.get('-exclude_ids')).toMatchObject({ kind: 'String', defaultValue: '' });
    expect(byFlag.get('-ids')).toMatchObject({ kind: 'String', defaultValue: '' });
    expect(byFlag.get('-meta')).toMatchObject({ kind: 'String', defaultValue: '' });
    expect(byFlag.get('-mem')).toMatchObject({
      kind: 'Integer',
      defaultValue: '1536',
      range: { min: 1, max: 2147483647, kind: 'integer' },
      constraintDescription: 'greater or equal to 1',
    });
    expect(byFlag.get('-unit')).toMatchObject({
      kind: 'Integer',
      range: { min: 1, max: 16, kind: 'integer' },
      constraintDescription: '1..16',
    });
    expect(byFlag.get('-unit')?.defaultValue).toBeUndefined();
    expect(byFlag.get('-t_extend')).toMatchObject({
      kind: 'Integer',
      range: { min: 0, max: 2147483647, kind: 'integer' },
    });
    expect(byFlag.get('-t_high_pct')).toMatchObject({
      kind: 'Real',
      range: { min: 0, max: 100, kind: 'double' },
      constraintDescription: '0..100',
    });
    expect(byFlag.get('-sformat')).toMatchObject({
      kind: 'String',
      defaultValue: 'ascii',
      choices: ['ascii', 'binary', 'oascii', 'obinary'],
      caseSensitive: true,
    });
    expect(byFlag.get('-outfmt')).toMatchObject({
      kind: 'String',
      defaultValue: 'interval',
      choices: [
        'fasta',
        'interval',
        'maskinfo_asn1_bin',
        'maskinfo_asn1_text',
        'maskinfo_xml',
        'seqloc_asn1_bin',
        'seqloc_asn1_text',
        'seqloc_xml',
      ],
      caseSensitive: true,
    });
    expect(byFlag.get('-infmt')).toMatchObject({
      kind: 'String',
      defaultValue: 'fasta',
      choices: ['blastdb', 'fasta', 'seqids'],
      caseSensitive: true,
    });
    expect(byFlag.get('-mk_counts')).toMatchObject({
      kind: 'Flag',
      takesValue: false,
      mandatory: false,
      group: '1',
      help: 'generate frequency counts for a database',
    });
    expect(byFlag.get('-h')).toMatchObject({ kind: 'Flag', takesValue: false });
    expect(byFlag.get('-h')?.group).toBeUndefined();
    expect(parsed.excludes).toContainEqual(['-mk_counts', '-ustat']);
    expect(parsed.excludes).toContainEqual(['-convert', '-outfmt']);
    expect(parsed.excludes).toContainEqual(['-ustat', '-sformat']);
  });

  it('parses mandatory flags, requirements, aliases, hidden/append flags and XML entities', () => {
    const makeblastdb = parseToolDescription(readFileSync(join(XMLHELP_DIR, 'makeblastdb.xml'), 'utf8'));
    const byFlag = new Map(makeblastdb.parameters.map((p) => [p.flag, p]));
    expect(byFlag.get('-dbtype')?.mandatory).toBe(true);
    expect(byFlag.get('-in')?.mandatory).toBe(false);
    expect(byFlag.get('-logfile')?.append).toBe(true);
    expect(byFlag.get('-randomize')?.hidden).toBe(true);
    expect(byFlag.get('-taxid')?.constraintDescription).toBe('>=0');
    expect(byFlag.get('-taxid_map')?.help).toContain('\nFormat:<SequenceId> <TaxonomyId><newline>');
    expect(makeblastdb.requires).toContainEqual(['-gi_mask', '-parse_seqids']);
    expect(makeblastdb.excludes).toContainEqual(['-taxid', '-taxid_map']);

    const convert = parseToolDescription(readFileSync(join(XMLHELP_DIR, 'convert2blastmask.xml'), 'utf8'));
    expect(convert.parameters.filter((p) => p.mandatory).map((p) => p.flag)).toEqual(['-masking_algorithm', '-masking_options']);

    const rpsbproc = parseToolDescription(readFileSync(join(XMLHELP_DIR, 'rpsbproc.xml'), 'utf8'));
    const rpsFlags = new Map(rpsbproc.parameters.map((p) => [p.flag, p]));
    expect(rpsFlags.get('-d')?.alias).toBe('-data-path');
    expect(rpsFlags.get('-f')?.alias).toBe('-show-families');
    expect(rpsFlags.get('-m')?.help).toContain('"rep"');
    expect(rpsbproc.detailedDescription).toContain('RPS-BLAST');

    expect(decodeXmlEntities('&lt;&gt;&amp;&apos;&quot;&#xa;&#x9;&#65;&amp;lt;')).toBe('<>&\'"\n\tA&lt;');
  });

  it('rejects malformed XML documents', () => {
    expect(() => parseToolDescription('')).toThrow();
    expect(() => parseToolDescription('<other/>')).toThrow(/ncbi_application/);
    expect(() => parseToolDescription('<ncbi_application><program><name>x</name></program></ncbi_application>')).toThrow(/arguments/);
    expect(() => parseToolDescription('<ncbi_application><program><name>x</wrong></program><arguments/></ncbi_application>')).toThrow(/mismatched/);
  });
});

describe('WasmEngine.describeTool', () => {
  const ASSETS = {
    serial: { url: '/serial.wasm', sha256: 'a'.repeat(64), size: 1, artifactSha256: 'a'.repeat(64) },
    threads: { url: '/threads.wasm', sha256: 'b'.repeat(64), size: 1, artifactSha256: 'c'.repeat(64), transform: 'shared-memory guard' },
  };
  const WINDOWMASKER_XML = readFileSync(join(XMLHELP_DIR, 'windowmasker.xml'), 'utf8');
  const UNUSED_CONTROL = {
    describe: () => Promise.reject(new Error('unused')),
    start: () => Promise.reject(new Error('unused')),
    validate: () => Promise.reject(new Error('unused')),
  };

  it('delegates to control.describeTool when provided', async () => {
    const expected = parseToolDescription(WINDOWMASKER_XML);
    const engine = new WasmEngine({
      assets: ASSETS,
      control: {
        ...UNUSED_CONTROL,
        describeTool: async (tool) => {
          expect(tool).toBe('windowmasker');
          return expected;
        },
      },
    });
    expect(await engine.describeTool('windowmasker')).toBe(expected);
  });

  it('falls back to runTool([tool, "-xmlhelp"]) when control has no describeTool', async () => {
    const encode = (value: string) => new TextEncoder().encode(value);
    let workerRef: { onmessage: ((event: MessageEvent<EngineEvent>) => void) | null } | undefined;
    const engine = new WasmEngine({
      assets: ASSETS,
      control: { ...UNUSED_CONTROL },
      createWorker: () => {
        const worker = {
          onmessage: null as ((event: MessageEvent<EngineEvent>) => void) | null,
          onerror: null,
          onmessageerror: null,
          postMessage: (command: EngineCommand) => {
            if (command.type === 'tool') {
              const toolCmd = command as WorkerTool;
              expect(toolCmd.argv).toEqual(['windowmasker', '-xmlhelp']);
              queueMicrotask(() => {
                worker.onmessage?.({ data: { type: 'tool-chunk', generation: toolCmd.generation, runId: toolCmd.runId, stream: 'stdout', bytes: encode(WINDOWMASKER_XML) } } as MessageEvent<EngineEvent>);
                worker.onmessage?.({
                  data: {
                    type: 'tool-done',
                    generation: toolCmd.generation,
                    runId: toolCmd.runId,
                    result: { status: 0, memory: { linearBytesBefore: 1, linearBytesAfter: 2, instanceRuns: 1 } },
                  },
                } as MessageEvent<EngineEvent>);
              });
            }
          },
          terminate: () => undefined,
        };
        workerRef = worker;
        return worker as unknown as Worker;
      },
      loadModules: async () => ({ serial: { asset: ASSETS.serial }, threads: { asset: ASSETS.threads } }),
    });
    const description = await engine.describeTool('windowmasker');
    expect(workerRef).toBeDefined();
    expect(description.tool).toBe('windowmasker');
    expect(description.parameters).toHaveLength(38);
  });
});

describe.skipIf(REACTORS === undefined)('reactorControl.describeTool on the serial reactor', () => {
  it('describes windowmasker from the reactor and matches windowmasker.xml', async () => {
    const { module } = await compileReactor(readFileSync(join(REACTORS!, 'losat-web-serial.wasm')));
    const open = async () => (await instantiateSerial(module)).abi;
    const control = reactorControl(open, open);
    const live = await control.describeTool!('windowmasker');
    const frozen = parseToolDescription(readFileSync(join(XMLHELP_DIR, 'windowmasker.xml'), 'utf8'));
    expect({ ...live, version: frozen.version }).toEqual(frozen);
  }, 60_000);
});
