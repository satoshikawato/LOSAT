import { describe, expect, it } from 'vitest';
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { toProgramDescription } from 'losat';
import { CUSTOM_OUTPUT, type ProgramId } from 'losat/domain';
import {
  buildArgv,
  commandWords,
  DATABASE_FLAGS,
  optionNames,
  optionValue,
  refusedFields,
  RESERVED_FLAGS,
  shellWords,
  splitOptions,
  toShellCommand,
} from '../../src/domain/argv';
import { PROGRAMS } from '../../src/domain/programs';

describe('buildArgv', () => {
  it('puts the program and input names first, then parameters in the given order', () => {
    const argv = buildArgv({
      program: 'blastn',
      queryName: 'q.fa',
      subjectName: 's.fa',
      parameters: [
        ['-task', 'blastn'],
        ['-evalue', '1e-5'],
        ['-lcase_masking', true],
      ],
    });
    expect(argv).toEqual(['blastn', '-query', 'q.fa', '-subject', 's.fa', '-task', 'blastn', '-evalue', '1e-5', '-lcase_masking']);
    expect(Object.isFrozen(argv)).toBe(true);
  });

  it.each(RESERVED_FLAGS)('rejects the application-managed flag %s', (flag) => {
    expect(() =>
      buildArgv({ program: 'blastp', queryName: 'q.fa', subjectName: 's.fa', parameters: [[flag, 'x']] }),
    ).toThrow(/managed by LOSAT/);
  });

  it('rejects parameters that are not flags', () => {
    expect(() =>
      buildArgv({ program: 'blastp', queryName: 'q.fa', subjectName: 's.fa', parameters: [['evalue', '1']] }),
    ).toThrow();
  });

  it('puts the words of Additional options after the parameters, as they are', () => {
    const argv = buildArgv({
      program: 'blastn',
      queryName: 'q.fa',
      subjectName: 's.fa',
      parameters: [['-evalue', '1e-5']],
      extra: ['-ungapped', '-dust', '20 64 1', '-qcov_hsp_perc=50'],
    });
    expect(argv).toEqual(['blastn', '-query', 'q.fa', '-subject', 's.fa', '-evalue', '1e-5', '-ungapped', '-dust', '20 64 1', '-qcov_hsp_perc=50']);
  });

  it.each(RESERVED_FLAGS.filter((flag) => flag !== '-num_threads').flatMap((flag) => [flag, `${flag}=x`]))(
    'rejects the application-managed %s among the words of Additional options',
    (word) => {
      expect(() =>
        buildArgv({ program: 'blastp', queryName: 'q.fa', subjectName: 's.fa', parameters: [], extra: ['-evalue', '1', word] }),
      ).toThrow(`${word.split('=')[0]} is managed by LOSAT`);
    },
  );

  it.each(['-num_threads', '-num_threads=2'])('names the control of the threads for %s in Additional options (W5 OP)', (word) => {
    expect(() =>
      buildArgv({ program: 'blastn', queryName: 'q.fa', subjectName: 's.fa', parameters: [], extra: [word, '2'] }),
    ).toThrow('-num_threads is set in the form: Threads');
  });

  // The engine reads -outfmt as the CLI does, and adds the output of the other formats and of 6 and 7
  // with fields (WA2-OP-03, W11): the application passes its words on, and the engine refuses what a
  // run does not write.
  it.each([
    [['-outfmt', '6 qseqid sseqid pident']],
    [['-outfmt=7 std qcovs']],
    [['-outfmt', '5']],
    [['-outfmt', '6']],
  ])('passes %j of Additional options to the engine as they are', (words) => {
    const argv = buildArgv({ program: 'blastp', queryName: 'q.fa', subjectName: 's.fa', parameters: [], extra: words });
    expect(argv).toEqual(['blastp', '-query', 'q.fa', '-subject', 's.fa', ...words]);
  });

  // The engine's validation drops a last empty word of its NUL-separated argv, and a run adds
  // -num_threads after a last --: the validation and the run would read them differently.
  it.each([
    ['', "Additional options: '' is not accepted"],
    ['--', 'Additional options: -- is not accepted'],
    ['a\0b', 'Additional options: a NUL character is not accepted'],
  ])('refuses the word %j among the words of Additional options', (word, message) => {
    expect(() =>
      buildArgv({ program: 'blastn', queryName: 'q.fa', subjectName: 's.fa', parameters: [], extra: ['-ungapped', word] }),
    ).toThrow(message);
  });

  // WA-FM-04 and WA-FM-05: the CLI prints its help for these words, and its message for the
  // others names -subject or -num_threads, which the browser sets itself.
  const extra = (program: ProgramId, words: readonly string[]) => () =>
    buildArgv({ program, queryName: 'q.fa', subjectName: 's.fa', parameters: [], extra: words });
  it.each(['-h', '-help', '-help-full', '--help'])('refuses %s, which asks the CLI for its help, in a few words', (word) => {
    for (const program of PROGRAMS) expect(extra(program.id, ['-evalue', '1', word]), program.id).toThrow(`${word} is not available in the browser`);
  });

  it('refuses the database options, also written -flag=value, as the browser build cannot search a database', () => {
    for (const program of PROGRAMS) {
      expect(extra(program.id, ['-db', 'nt']), program.id).toThrow(/^-db \(a database search\) is not supported in the browser yet$/);
      for (const flag of DATABASE_FLAGS.slice(1)) {
        expect(extra(program.id, [`${flag}=x`]), `${program.id} ${flag}`).toThrow(`${flag} (a database search filter) is not supported in the browser yet`);
      }
    }
  });

  it("lists the database options of every mode's engine (describe)", () => {
    const described = JSON.parse(readFileSync(resolve(import.meta.dirname, '../../src/infra/fake/describe.json'), 'utf8')) as Record<string, unknown>;
    for (const program of PROGRAMS) {
      const flags = toProgramDescription(described[program.id]).parameters.map((option) => option.flag);
      expect(flags, program.id).toEqual(expect.arrayContaining([...DATABASE_FLAGS]));
    }
  });

  it('refuses -mt_mode, which goes with the threads that LOSAT sets, where the CLI has it', () => {
    for (const program of ['blastn', 'blastp', 'blastx', 'tblastn'] as const) {
      expect(extra(program, ['-mt_mode', '1']), program).toThrow('-mt_mode is managed by LOSAT');
    }
    // TBLASTX's CLI has no -mt_mode: the engine's message says so.
    expect(extra('tblastx', ['-mt_mode', '1'])()).toEqual(['tblastx', '-query', 'q.fa', '-subject', 's.fa', '-mt_mode', '1']);
  });
});

describe('optionNames (WA-FM-07)', () => {
  const valued = new Set(['-task', '-penalty', '-dust', '-strand', '-evalue']);
  const takesValue = (flag: string) => valued.has(flag);

  it("tells an option's name from a value that starts with '-', with the options that take a value", () => {
    const words = ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-penalty', '-3', '-lcase_masking', '-dust', '20 64 1', '-strand=minus', '-evalue', '-1'];
    const parts = optionNames(words, takesValue);
    expect(parts.filter((part) => part.name).map((part) => part.text)).toEqual(['-query', '-subject', '-penalty', '-lcase_masking', '-dust', '-strand=minus', '-evalue']);
    // The text is the command line's, so that a copy reads back as the same words.
    const text = parts.map((part) => part.text).join(' ');
    expect(text).toBe(shellWords(words));
    expect(splitOptions(text)).toEqual(words);
  });

  it('reads the options that the application manages as taking a value (the command of an output format)', () => {
    const parts = optionNames(commandWords(['blastn', '-query', '-q.fa', '-subject', 's.fa'], 6), takesValue);
    expect(parts.filter((part) => part.name).map((part) => part.text)).toEqual(['-query', '-subject', '-outfmt']);
    const custom = optionNames(commandWords(['blastn', '-query', 'q.fa', '-outfmt', '-6 qseqid', '-subject', 's.fa'], CUSTOM_OUTPUT), takesValue);
    expect(custom.filter((part) => part.name).map((part) => part.text)).toEqual(['-query', '-outfmt', '-subject']);
  });

  it("names every word that starts with '-' while the grammar is unknown", () => {
    expect(optionNames(['-penalty', '-3'], undefined).map((part) => part.name)).toEqual([true, true]);
  });
});

describe('splitOptions', () => {
  it.each([
    ['', []],
    ['  \t\n ', []],
    ['-ungapped  -qcov_hsp_perc 50', ['-ungapped', '-qcov_hsp_perc', '50']],
    ["-dust '20 64 1'", ['-dust', '20 64 1']],
    ['-dust "20 64 1" -penalty -3', ['-dust', '20 64 1', '-penalty', '-3']],
    ['-dust 20\\ 64\\ 1', ['-dust', '20 64 1']],
    ["a'b c'd", ['ab cd']],
    ["''", ['']],
    ['"a \\"b\\" \\\\ \\x"', ['a "b" \\ \\x']],
    ["'a\\b'", ['a\\b']],
  ])('reads %j as a shell does', (text, words) => {
    expect(splitOptions(text)).toEqual(words);
  });

  it('refuses an unclosed quote', () => {
    expect(() => splitOptions("-dust '20 64 1")).toThrow("Additional options: a quote (') is not closed");
    expect(() => splitOptions('-dust "20 64 1')).toThrow('Additional options: a quote (") is not closed');
  });

  it('reads back the words of the CLI command', () => {
    const words = ['-dust', '20 64 1', "it's", '', 'a"b', 'back\\slash'];
    const command = toShellCommand(['blastn', ...words], 6);
    expect(splitOptions(command)).toEqual(['losat', 'blastn', ...words, '-outfmt', '6']);
  });
});

describe('shellWords', () => {
  it('quotes the words that need it, so that they read back', () => {
    const words = ['-task', 'blastn', '-dust', '20 64 1', '-penalty', '-3', "it's"];
    expect(shellWords(words)).toBe("-task blastn -dust '20 64 1' -penalty -3 'it'\\''s'");
    expect(splitOptions(shellWords(words))).toEqual(words);
  });
});

describe('toShellCommand', () => {
  it('quotes only the words that need quoting and always writes -outfmt', () => {
    expect(toShellCommand(['blastn', '-query', 'my query.fa', '-subject', "it's.fa"], 6)).toBe(
      "losat blastn -query 'my query.fa' -subject 'it'\\''s.fa' -outfmt 6",
    );
  });

  // The argv's own -outfmt (WA2-OP-03): the command of a format that every run writes has that
  // format's -outfmt in its place; the command of the output of the argv's -outfmt is the argv.
  it("writes the format's -outfmt for the formats of every run, and the argv's own once for its output", () => {
    const argv = ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-outfmt', '6 qseqid sseqid pident', '-evalue', '1'];
    expect(toShellCommand(argv, 0)).toBe('losat blastn -query q.fa -subject s.fa -evalue 1 -outfmt 0');
    expect(toShellCommand(argv, 7)).toBe('losat blastn -query q.fa -subject s.fa -evalue 1 -outfmt 7');
    expect(toShellCommand(argv, CUSTOM_OUTPUT)).toBe("losat blastn -query q.fa -subject s.fa -outfmt '6 qseqid sseqid pident' -evalue 1");
    const joined = ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-outfmt=7 std qcovs'];
    expect(commandWords(joined, 6)).toEqual(['losat', 'blastn', '-query', 'q.fa', '-subject', 's.fa', '-outfmt', '6']);
    expect(commandWords(joined, CUSTOM_OUTPUT)).toEqual(['losat', ...joined]);
    // -outfmt 6 adds no output, and the commands of the formats leave it out too.
    expect(commandWords(['blastp', '-outfmt', '6'], 0)).toEqual(['losat', 'blastp', '-outfmt', '0']);
    // The CLI writes 13 and 14 as files: the command names the -out of the archive's files.
    expect(toShellCommand(['blastp', '-query', 'q.fa', '-outfmt', '13'], CUSTOM_OUTPUT)).toBe('losat blastp -query q.fa -outfmt 13 -out results');
    expect(toShellCommand(['blastp', '-outfmt=+014 x'], CUSTOM_OUTPUT)).toBe("losat blastp '-outfmt=+014 x' -out results");
    expect(toShellCommand(['blastp', '-outfmt', '14'], 6)).toBe('losat blastp -outfmt 6');
  });
});

describe('optionValue', () => {
  it("reads the word after an option's last occurrence, and nothing where it is not given", () => {
    const argv = ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-task', 'blastn', '-evalue', '1e-5', '-task', 'dc-megablast'];
    expect(optionValue(argv, '-task')).toBe('dc-megablast');
    expect(optionValue(argv, '-evalue')).toBe('1e-5');
    expect(optionValue(argv, '-word_size')).toBeUndefined();
    expect(optionValue(['blastn', '-task'], '-task')).toBeUndefined();
  });
});

describe('refusedFields', () => {
  const fields = (message: string, extra = '') => {
    const { extra: inExtra, flags } = refusedFields(message, extra);
    return { extra: inExtra, flags: [...flags].sort() };
  };
  it("names the form's fields of the options that a message names, and no words of its sentences", () => {
    expect(fields("error: invalid value '3' for '-word_size <WORD_SIZE>': expected an integer >= 4\n\nFor more information, try '--help'.")).toEqual({
      extra: false,
      flags: ['-word_size'],
    });
    expect(fields('BLAST query/options error: Invalid discontiguous template parameters: word size must be either 11 or 12')).toEqual({
      extra: false,
      flags: [],
    });
    expect(fields('the e-value of -penalty -3 and -2')).toEqual({ extra: false, flags: ['-penalty'] });
  });
  it('names "Additional options" for the options typed there, also those of the form', () => {
    expect(fields('-evalue is set in the form: Expect threshold', '-evalue abc')).toEqual({ extra: true, flags: [] });
    expect(fields("error: unknown option or argument '-foo'; use -help for CLI v2 syntax", '-qcov_hsp_perc 50 -foo')).toEqual({
      extra: true,
      flags: ['-help'],
    });
    expect(fields('-num_threads is set in the form: Threads', '-num_threads=2')).toEqual({ extra: true, flags: [] });
    expect(fields("error: invalid value 'x' for '-qcov_hsp_perc <QCOV_HSP_PERC>'", '-qcov_hsp_perc x')).toEqual({ extra: true, flags: [] });
  });
  it('does not read the options of clap\'s Usage line, which the argv lists, as named by the sentence', () => {
    const message = [
      "error: the argument '-max_target_seqs <MAX_TARGET_SEQS>' cannot be used with '-num_alignments <int_value>'",
      '',
      'Usage: losat blastn -query <PATH> -subject <PATH> -evalue <EVALUE> -max_target_seqs <MAX_TARGET_SEQS>',
      '',
      "For more information, try '--help'.",
    ].join('\n');
    expect(fields(message)).toEqual({ extra: false, flags: ['-max_target_seqs', '-num_alignments'] });
    // --help of the closing line is not the option typed in Additional options.
    expect(fields(message, '--help')).toEqual({ extra: false, flags: ['-max_target_seqs', '-num_alignments'] });
  });
  it('names "Additional options" for a word typed there that is not an option name', () => {
    expect(fields('--help is not available in the browser', '--help')).toEqual({ extra: true, flags: [] });
    expect(fields('--help is not available in the browser', '-strand plus')).toEqual({ extra: false, flags: [] });
    expect(fields("error: unknown option or argument 'evalue'; use -help for CLI v2 syntax", 'evalue 1e-5')).toEqual({
      extra: true,
      flags: ['-help'],
    });
    expect(fields("error: unknown option or argument 'evalue'", '-strand plus')).toEqual({ extra: false, flags: [] });
    // The engine quotes the whole word, with its =value.
    expect(fields("error: unknown option or argument '--evalue=1e-5'; use -help for CLI v2 syntax", '--evalue=1e-5')).toEqual({
      extra: true,
      flags: ['-help'],
    });
    expect(fields("error: unknown option or argument 'evalue=10'", '-strand plus evalue=10')).toEqual({ extra: true, flags: [] });
  });
  it('names "Additional options" for the messages about its text', () => {
    expect(fields("Additional options: a quote (') is not closed", "-strand 'plus")).toEqual({ extra: true, flags: [] });
    expect(fields("Additional options: '' is not accepted", "''")).toEqual({ extra: true, flags: [] });
  });
});
