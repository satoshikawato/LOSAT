import { describe, expect, it } from 'vitest';
import { buildArgv, optionValue, RESERVED_FLAGS, shellWords, splitOptions, toShellCommand } from '../../src/domain/argv';

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
    ).toThrow(/managed by LOSAT Web/);
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

  it.each(RESERVED_FLAGS.flatMap((flag) => [flag, `${flag}=x`]))('rejects the application-managed %s among the words of Additional options', (word) => {
    expect(() =>
      buildArgv({ program: 'blastp', queryName: 'q.fa', subjectName: 's.fa', parameters: [], extra: ['-evalue', '1', word] }),
    ).toThrow(`${word.split('=')[0]} is managed by LOSAT`);
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
