import { describe, expect, it } from 'vitest';
import { buildArgv, RESERVED_FLAGS, toShellCommand } from '../../src/domain/argv';

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
});

describe('toShellCommand', () => {
  it('quotes only the words that need quoting', () => {
    expect(toShellCommand(['blastn', '-query', 'my query.fa', '-subject', "it's.fa"])).toBe(
      "LOSAT blastn -query 'my query.fa' -subject 'it'\\''s.fa'",
    );
  });
});
