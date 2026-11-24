import { describe, expect, it } from 'vitest';
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { toProgramDescription } from 'losat';
import { CUSTOM_OUTPUT, type ProgramId } from 'losat/domain';
import {
  DATABASE_MASK_FLAGS,
  DATABASE_MASK_PROGRAMS,
  buildArgv,
  commandWords,
  DATABASE_FLAGS,
  optionNames,
  optionValue,
  refusedFields,
  runTask,
  RESERVED_FLAGS,
  shellWords,
  splitOptions,
  userArgv,
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

  // WA-FM-04: the CLI's message for these words names -subject or -num_threads, which the browser sets itself.
  const extra = (program: ProgramId, words: readonly string[]) => () =>
    buildArgv({ program, queryName: 'q.fa', subjectName: 's.fa', parameters: [], extra: words });
  // Line WQA's D14, 2026-11-20: the engine prints the command line's help or version for these words, as the CLI does
  // (WA-FM-05 refused the help words).
  it.each(['-h', '-help', '-help-full', '-xmlhelp', '-version', '-version-full', '-version-full-xml', '-version-full-json'])(
    'passes %s, which asks the CLI for its help or version, to the engine',
    (word) => {
      for (const program of PROGRAMS) expect(extra(program.id, ['-evalue', '1', word])(), program.id).toContain(word);
    },
  );

  it('refuses the database options, also written -flag=value, as the browser build cannot search a database', () => {
    for (const program of PROGRAMS) {
      expect(extra(program.id, ['-db', 'nt']), program.id).toThrow(/^-db \(a database search\) is not supported in the browser yet$/);
      for (const flag of DATABASE_FLAGS.slice(1)) {
        expect(extra(program.id, [`${flag}=x`]), `${program.id} ${flag}`).toThrow(`${flag} (a database search filter) is not supported in the browser yet`);
      }
    }
  });

  // The database masks of every mode's CLI: their parser would name -subject, which the browser
  // sets itself (WA-FM-04).
  it('refuses the database masks of every mode, also written -flag=value', () => {
    expect(DATABASE_MASK_PROGRAMS).toEqual(['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']);
    for (const program of DATABASE_MASK_PROGRAMS) {
      for (const flag of DATABASE_MASK_FLAGS) {
        expect(extra(program, [flag, '11']), `${program} ${flag}`).toThrow(`${flag} (a database search filter) is not supported in the browser yet`);
        expect(extra(program, [`${flag}=dust`]), `${program} ${flag}`).toThrow(`${flag} (a database search filter) is not supported in the browser yet`);
      }
    }
  });

  // BLASTN's -filtering_db searches a BLAST database of repeats: the engine reads it from the run's directory, but the
  // form has no input for its files yet (D4-W2, with WQA D12b; Owner 2026-11-22).
  it('refuses the database of repeats of BLASTN, also written -flag=value', () => {
    for (const words of [['-filtering_db', 'repeats'], ['-filtering_db=repeats']]) {
      expect(extra('blastn', words)).toThrow(/^-filtering_db \(a BLAST database of repeats\) is not supported in the browser yet$/);
    }
  });

  it("lists the database masks in the engine's options of the modes that have them (describe)", () => {
    const described = JSON.parse(readFileSync(resolve(import.meta.dirname, '../../src/infra/fake/describe.json'), 'utf8')) as Record<string, unknown>;
    for (const program of PROGRAMS) {
      const flags = toProgramDescription(described[program.id]).parameters.map((option) => option.flag);
      const has = DATABASE_MASK_FLAGS.filter((flag) => flags.includes(flag));
      expect(has, program.id).toEqual(DATABASE_MASK_PROGRAMS.includes(program.id) ? [...DATABASE_MASK_FLAGS] : []);
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

  it('reads -flag=value as the engine does', () => {
    expect(optionValue(['blastp', '-task=blastp-short'], '-task')).toBe('blastp-short');
    expect(optionValue(['blastp', '-task', 'blastp', '-task=blastp-short'], '-task')).toBe('blastp-short');
    expect(optionValue(['blastp', '-task=blastp-short', '-task', 'blastp-fast'], '-task')).toBe('blastp-fast');
    expect(optionValue(['blastp', '-task_x=1'], '-task')).toBeUndefined();
  });
});

describe('runTask', () => {
  const options = [{ flag: '-task', defaultValue: 'blastp' }, { flag: '-evalue' }];

  it('is the task that the engine runs: the argv\'s -task in either form, else the default task', () => {
    expect(runTask(['blastp', '-task', 'blastp-short'], options)).toBe('blastp-short');
    // The engine reads -task=blastp-short as -task blastp-short (NCBI's reader).
    expect(runTask(['blastp', '-task=blastp-short'], options)).toBe('blastp-short');
    expect(runTask(['blastp', '-evalue', '1'], options)).toBe('blastp');
    // TBLASTX has no -task; without the engine's description only the argv's -task is known.
    expect(runTask(['tblastx'], [{ flag: '-evalue' }])).toBeUndefined();
    expect(runTask(['blastp', '-task=blastp-short'], [])).toBe('blastp-short');
  });
});

describe('refusedFields', () => {
  const fields = (message: string, extra = '') => {
    const { extra: inExtra, flags } = refusedFields(message, extra);
    return { extra: inExtra, flags: [...flags].sort() };
  };
  // The engine's argument error: a usage that lists every option of the program, a line of 72 "=",
  // an empty line, and two posts (the second with two spaces after "Error:").
  const USAGE = [
    'USAGE',
    '  blastn [-h] [-help] [-query input_file] [-subject subject_input_file]',
    '    [-evalue evalue] [-word_size int_value] [-culling_limit int_value]',
    '    [-best_hit_overhang float_value] [-ungapped] [-num_threads int_value]',
    '',
    'DESCRIPTION',
    '   Nucleotide-Nucleotide BLAST (LOSAT)',
    '',
    "Use '-help' to print detailed descriptions of command line arguments",
    '='.repeat(72),
    '',
  ];
  const ncbi = (post: string) => [...USAGE, `Error: ${post}`, `Error:  ${post}`, ''].join('\n');
  it("names the form's fields of the options that a message names, and no words of its sentences", () => {
    expect(fields('LOSAT query/options error: Invalid discontiguous template parameters: word size must be either 11 or 12')).toEqual({
      extra: false,
      flags: [],
    });
    expect(fields('the e-value of -penalty -3 and -2')).toEqual({ extra: false, flags: ['-penalty'] });
  });
  it("names the options of NCBI's argument error posts, and not those of its usage", () => {
    // Both options of "Incompatible with" (the first by its name, the second in backquotes).
    expect(fields(ncbi('Argument "culling_limit". Incompatible with argument:  `best_hit_overhang\''))).toEqual({
      extra: false,
      flags: ['-best_hit_overhang', '-culling_limit'],
    });
    expect(fields(ncbi('Argument "mask_data". Must be specified, as it is required by argument:  `mask_id\''))).toEqual({
      extra: false,
      flags: ['-mask_data', '-mask_id'],
    });
    // The value in backquotes is not an option.
    expect(fields(ncbi('Argument "word_size". Illegal value, expected >=4:  `3\''))).toEqual({ extra: false, flags: ['-word_size'] });
    expect(fields(ncbi('Argument "evalue". Argument cannot be converted:  `abc\''))).toEqual({ extra: false, flags: ['-evalue'] });
    // Only "Value is missing" writes the dash.
    expect(fields(ncbi('Argument "-evalue". Value is missing'))).toEqual({ extra: false, flags: ['-evalue'] });
    expect(fields(ncbi('Argument "dbtype". Mandatory value is missing:  `String, `nucl\', `prot\''))).toEqual({ extra: false, flags: ['-dbtype'] });
    expect(fields(ncbi('Argument with this name is defined already: evalue'))).toEqual({ extra: false, flags: ['-evalue'] });
  });
  it("names \"Additional options\" for the words that NCBI's posts name there, also those of the form", () => {
    expect(fields(ncbi('Argument "culling_limit". Incompatible with argument:  `best_hit_overhang\''), '-culling_limit 1 -best_hit_overhang 0.1')).toEqual({
      extra: true,
      flags: [],
    });
    // Only one of the two typed there: the other is the form's field.
    expect(fields(ncbi('Argument "culling_limit". Incompatible with argument:  `best_hit_overhang\''), '-best_hit_overhang 0.1')).toEqual({
      extra: true,
      flags: ['-culling_limit'],
    });
    expect(fields(ncbi('Argument "-evalue". Value is missing'), '-evalue')).toEqual({ extra: true, flags: [] });
    expect(fields(ncbi('Argument with this name is defined already: evalue'), '-evalue 1 -evalue 2')).toEqual({ extra: true, flags: [] });
    // The usage lists -word_size and -ungapped, which the extra box types: not named.
    expect(fields(ncbi('Argument "evalue". Argument cannot be converted:  `abc\''), '-word_size 11 -ungapped')).toEqual({
      extra: false,
      flags: ['-evalue'],
    });
  });
  it('reads the typed word that Unknown argument, Invalid argument and the offending value name', () => {
    // `-foo` is "foo": the word without its first dash; `--foo` is "-foo".
    expect(fields(ncbi('Unknown argument: "foo"'), '-qcov_hsp_perc 50 -foo')).toEqual({ extra: true, flags: [] });
    expect(fields(ncbi('Unknown argument: "foo"'), '-foo=3')).toEqual({ extra: true, flags: [] });
    expect(fields(ncbi('Unknown argument: "-foo"'), '--foo')).toEqual({ extra: true, flags: [] });
    expect(fields(ncbi('Unknown argument: "foo"'), '-strand plus')).toEqual({ extra: false, flags: [] });
    // The whole typed word.
    expect(fields(ncbi('Invalid argument: -ungapped=1'), '-ungapped=1')).toEqual({ extra: true, flags: [] });
    expect(fields(ncbi('Invalid argument: -ungapped=1'), '-strand plus')).toEqual({ extra: false, flags: [] });
    expect(fields(ncbi('Too many positional arguments (1), the offending value: evalue'), 'evalue 1e-5')).toEqual({ extra: true, flags: [] });
    expect(fields(ncbi('Too many positional arguments (1), the offending value: evalue'), '-strand plus')).toEqual({ extra: false, flags: [] });
    // CRLF line ends are read alike.
    expect(fields(ncbi('Unknown argument: "foo"').replaceAll('\n', '\r\n'), '-foo')).toEqual({ extra: true, flags: [] });
  });
  it("reads LOSAT's one-line refusals whole, with the dash of their names", () => {
    expect(fields("the NCBI BLAST+ option -remote is not supported by LOSAT's BLASTN", '-remote')).toEqual({ extra: true, flags: [] });
    expect(fields("the NCBI BLAST+ option -remote is not supported by LOSAT's BLASTN", '-strand plus')).toEqual({
      extra: false,
      flags: ['-remote'],
    });
    expect(fields("the option -dryrun runs on the command line of LOSAT's BLASTN only: use it there")).toEqual({ extra: false, flags: ['-dryrun'] });
  });
  it('names "Additional options" for the options typed there, also those of the form', () => {
    expect(fields('-evalue is set in the form: Expect threshold', '-evalue abc')).toEqual({ extra: true, flags: [] });
    expect(fields('-num_threads is set in the form: Threads', '-num_threads=2')).toEqual({ extra: true, flags: [] });
  });
  it('names "Additional options" for a word typed there that is not an option name', () => {
    expect(fields('--help is not available in the browser', '--help')).toEqual({ extra: true, flags: [] });
    expect(fields('--help is not available in the browser', '-strand plus')).toEqual({ extra: false, flags: [] });
    // A word that a message quotes, with its =value.
    expect(fields("the word 'evalue' is not read", 'evalue 1e-5')).toEqual({ extra: true, flags: [] });
    expect(fields("the word 'evalue' is not read", '-strand plus')).toEqual({ extra: false, flags: [] });
    expect(fields("the word '--evalue=1e-5' is not read", '--evalue=1e-5')).toEqual({ extra: true, flags: [] });
  });
  it('names "Additional options" for the messages about its text', () => {
    expect(fields("Additional options: a quote (') is not closed", "-strand 'plus")).toEqual({ extra: true, flags: [] });
    expect(fields("Additional options: '' is not accepted", "''")).toEqual({ extra: true, flags: [] });
  });
});

// Line WQA's D20, 2026-11-23: a run that printed the command line's help or version is shown as the command the user typed.
describe('userArgv', () => {
  const argv = ['blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-version-full', '-num_threads', '2'];
  it('leaves out the placeholder inputs of a printed run, which read none', () => {
    expect(userArgv({ argv, prints: 'version' })).toEqual(['blastn', '-version-full', '-num_threads', '2']);
    expect(userArgv({ argv: ['blastp', '-query', 'q.fa', '-subject', 's.fa', '-h'], prints: 'help' })).toEqual(['blastp', '-h']);
  });
  it('keeps the options and the words of Additional options, such as -conffile and -outfmt', () => {
    expect(userArgv({ argv: ['blastn', '-query', 'query.fa', '-subject', 'subject.fa', '-help', '-conffile', 's.ini', '-outfmt', '17'], prints: 'help' })).toEqual([
      'blastn',
      '-help',
      '-conffile',
      's.ini',
      '-outfmt',
      '17',
    ]);
  });
  it('is the argv itself for a search', () => {
    expect(userArgv({ argv })).toBe(argv);
  });
});
