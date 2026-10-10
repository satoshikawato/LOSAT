// Settings files (S15 item 3): the document, and every refusal of a file that is not exactly it.
import { describe, expect, it } from 'vitest';
import {
  MAX_THREADS,
  parseSettings,
  readOptions,
  serializeSettings,
  settingsFileName,
  settingsOfRun,
  settingsText,
  SETTINGS_FORMAT,
  SETTINGS_MAX_BYTES,
  SETTINGS_NOTE,
  threadLimit,
  type SearchSettings,
} from '../../src/domain/settings-file';

const encode = (text: string) => new TextEncoder().encode(text);
const document = (fields: Record<string, unknown>) =>
  encode(JSON.stringify({ format: SETTINGS_FORMAT, schema: 1, program: 'blastn', options: [], threads: 'auto', ...fields }));
const refusal = (bytes: Uint8Array) => {
  const read = parseSettings(bytes);
  if (read.ok) throw new Error('expected a refusal');
  return read.message;
};

describe('settings file', () => {
  const settings: SearchSettings = {
    program: 'blastn',
    options: ['-task', 'blastn', '-penalty', '-3', '-lcase_masking', '-dust', '20 64 1', '-query_loc', '2-50'],
    threads: 2,
  };

  it('writes the program, the options and the threads with the format, the schema and a note, as JSON.stringify does', () => {
    const text = serializeSettings(settings);
    const expected = { format: SETTINGS_FORMAT, schema: 1, note: SETTINGS_NOTE, program: 'blastn', options: settings.options, threads: 2 };
    expect(text).toBe(`${JSON.stringify(expected, null, 2)}\n`);
    expect(serializeSettings({ program: 'tblastx', options: [], threads: 'auto' })).toBe(
      `${JSON.stringify({ ...expected, program: 'tblastx', options: [], threads: 'auto' }, null, 2)}\n`,
    );
    // One part per option word: the file is written in blocks, never built whole.
    expect([...settingsText(settings)].length).toBeGreaterThan(settings.options.length);
    expect(SETTINGS_NOTE).toMatch(/LOSAT Web application format, not an NCBI BLAST\+ file/);
    expect(settingsFileName('tblastn')).toBe('losat-settings-tblastn.json');
  });

  it('reads back what it writes', () => {
    expect(parseSettings(encode(serializeSettings(settings)))).toEqual({ ok: true, settings });
    expect(parseSettings(document({ note: undefined }))).toEqual({ ok: true, settings: { program: 'blastn', options: [], threads: 'auto' } });
  });

  it('takes the settings of a run from its argv after the inputs, and its threads', () => {
    const argv = ['tblastn', '-query', 'q.faa', '-subject', 'combined_subject.fa', '-db_gencode', '4', '-subject_loc', '10-90'];
    expect(settingsOfRun({ program: 'tblastn', argv, requestedThreads: 'auto' })).toEqual({
      program: 'tblastn',
      options: ['-db_gencode', '4', '-subject_loc', '10-90'],
      threads: 'auto',
    });
  });

  it('refuses a file that is too large, not UTF-8, or not JSON', () => {
    expect(refusal(new Uint8Array(SETTINGS_MAX_BYTES + 1))).toBe(
      `The file has ${SETTINGS_MAX_BYTES + 1} bytes; a settings file has at most ${SETTINGS_MAX_BYTES} (1 MiB).`,
    );
    expect(refusal(new Uint8Array([0x7b, 0xff, 0x7d]))).toBe('The file is not UTF-8 text, so it is not a LOSAT Web settings file.');
    expect(refusal(encode('{"format": '))).toMatch(/^The file is not JSON \(.+\)\.$/);
    expect(refusal(encode('[1, 2]'))).toBe('The file is not a LOSAT Web settings file: it is not a JSON object.');
    expect(refusal(encode('null'))).toBe('The file is not a LOSAT Web settings file: it is not a JSON object.');
  });

  it('refuses another format, a newer or a wrong schema, and fields that schema 1 does not have', () => {
    expect(refusal(document({ format: 'LOSAT Web session' }))).toBe(
      'The file is not a LOSAT Web settings file: its "format" is "LOSAT Web session", not "LOSAT Web search settings".',
    );
    expect(refusal(document({ format: undefined }))).toBe(
      'The file is not a LOSAT Web settings file: its "format" is missing, not "LOSAT Web search settings".',
    );
    expect(refusal(document({ schema: 2, extra: true }))).toBe('The file has schema 2: a newer LOSAT Web saved it. This one reads schema 1.');
    expect(refusal(document({ schema: '1' }))).toBe('"schema" must be a whole number from 1 (it is "1").');
    expect(refusal(document({ schema: 0 }))).toBe('"schema" must be a whole number from 1 (it is 0).');
    expect(refusal(document({ schema: 1.5 }))).toBe('"schema" must be a whole number from 1 (it is 1.5).');
    expect(refusal(document({ query: 'q.fa' }))).toBe('The file has a field that schema 1 does not have: "query".');
    expect(refusal(document({ note: 3 }))).toBe('"note" must be a text of at most 10000 characters.');
    expect(refusal(document({ note: 'x'.repeat(10_001) }))).toBe('"note" must be a text of at most 10000 characters.');
  });

  it('refuses an unknown program', () => {
    expect(refusal(document({ program: 'megablast' }))).toBe(
      '"program" must be one of blastn, blastp, blastx, tblastn, tblastx (it is "megablast").',
    );
    expect(refusal(document({ program: undefined }))).toBe(
      '"program" must be one of blastn, blastp, blastx, tblastn, tblastx (it is missing).',
    );
  });

  it('refuses options that are not a bounded list of bounded words, or that hold a flag LOSAT Web sets itself', () => {
    expect(refusal(document({ options: '-evalue 1' }))).toBe('"options" must be a list of words (it is "-evalue 1").');
    expect(refusal(document({ options: undefined }))).toBe('"options" must be a list of words (it is missing).');
    expect(refusal(document({ options: Array(1001).fill('-lcase_masking') }))).toBe(
      '"options" has 1001 words; a settings file has at most 1000.',
    );
    expect(parseSettings(document({ options: Array(1000).fill('-lcase_masking') })).ok).toBe(true);
    expect(refusal(document({ options: ['-evalue', 10] }))).toBe('"options" word 2 must be a text (it is 10).');
    expect(refusal(document({ options: ['-dust', 'x'.repeat(10_001)] }))).toBe(
      '"options" word 2 has 10001 characters; a word has at most 10000.',
    );
    expect(parseSettings(document({ options: ['-dust', 'x'.repeat(10_000)] })).ok).toBe(true);
    for (const flag of ['-query', '-subject', '-out', '-outfmt', '-num_threads']) {
      expect(refusal(document({ options: ['-evalue', '1', flag, 'x'] }))).toBe(
        `"options" word 3 is ${flag}, which LOSAT Web sets itself: a settings file has no inputs, outputs or thread count among its options.`,
      );
    }
  });

  it('refuses threads other than "auto" or a whole number in the form range', () => {
    for (const threads of [0, MAX_THREADS + 1, 2.5, '4', 'Auto', null, undefined]) {
      expect(refusal(document({ threads }))).toMatch(/^"threads" must be "auto" or a whole number from 1 to 16 \(it is .+\)\.$/);
    }
    expect(parseSettings(document({ threads: MAX_THREADS })).ok).toBe(true);
    expect(parseSettings(document({ threads: 1 })).ok).toBe(true);
  });

  it('cuts a long value in a message', () => {
    expect(refusal(document({ program: 'x'.repeat(100) }))).toBe(
      `"program" must be one of blastn, blastp, blastx, tblastn, tblastx (it is "${'x'.repeat(56)}...).`,
    );
  });

  it('offers the threads of the processor, at most 16, and 4 where the browser does not say', () => {
    expect(threadLimit(8)).toBe(8);
    expect(threadLimit(64)).toBe(16);
    expect(threadLimit(undefined)).toBe(4);
    expect(threadLimit(0)).toBe(4);
  });
});

describe('readOptions', () => {
  const grammar: Record<string, boolean> = { '-penalty': true, '-evalue': true, '-lcase_masking': false, '-dust': true };
  const takesValue = (flag: string) => grammar[flag];

  it('reads a value that starts with "-" after a flag that takes a value, and bare flags', () => {
    expect(readOptions(['-penalty', '-3', '-lcase_masking', '-evalue', '-1e-5'], takesValue).options).toEqual([
      { flag: '-penalty', value: '-3', words: ['-penalty', '-3'] },
      { flag: '-lcase_masking', value: true, words: ['-lcase_masking'] },
      { flag: '-evalue', value: '-1e-5', words: ['-evalue', '-1e-5'] },
    ]);
  });

  it('gives an undescribed flag the next word unless it looks like an option, and reports stray words', () => {
    const read = readOptions(['10', '-foo', 'bar', '-baz', '-3', '-qux', '-lcase_masking', '-evalue'], takesValue);
    expect(read.options.map((option) => option.words)).toEqual([['-foo', 'bar'], ['-baz', '-3'], ['-qux'], ['-lcase_masking']]);
    expect(read.stray).toEqual(['10', '-evalue']);
  });
});
