// The reproduction of a run (S15 item 5, design §12.3): the LOSAT and NCBI BLAST+ commands of an
// argv, why a run has no NCBI command, and the notes. With LOSAT_WEB_COMMANDS_OUT set, the
// commands of a fixed list of argvs are written to that JSON file, which
// docs/evidence/losat_web_w6/check_commands.py runs with both programs and compares byte for byte.
import { writeFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { buildArgv } from '../../src/domain/argv';
import type { OutputFormat } from '../../src/domain/output-format';
import type { ProgramId } from '../../src/domain/programs';
import {
  commandNotes,
  inputRelation,
  losatCommand,
  NCBI_BLAST_VERSION,
  ncbiCommand,
  ncbiComparison,
} from '../../src/domain/reproduce';
import { approvedExceptions } from '../../src/domain/verification';

const argvOf = (program: ProgramId, query: string, subject: string, words: readonly string[]) => [program, '-query', query, '-subject', subject, ...words];

describe('commands', () => {
  it('gives the LOSAT command and the NCBI command of the same program, inputs and options, quoted alike, without -num_threads', () => {
    const argv = argvOf('blastn', "my query (1).fa", "subject's.fa", ['-task', 'blastn', '-penalty', '-3', '-dust', '20 64 1', '-lcase_masking']);
    expect(losatCommand(argv, 6)).toBe(
      "LOSAT blastn -query 'my query (1).fa' -subject 'subject'\\''s.fa' -task blastn -penalty -3 -dust '20 64 1' -lcase_masking -outfmt 6",
    );
    expect(ncbiCommand(argv, 6)).toBe(
      "blastn -query 'my query (1).fa' -subject 'subject'\\''s.fa' -task blastn -penalty -3 -dust '20 64 1' -lcase_masking -outfmt 6",
    );
    expect(ncbiCommand(argvOf('tblastx', 'q.fa', 's.fa', []), 0)).toBe('tblastx -query q.fa -subject s.fa -outfmt 0');
    expect(ncbiCommand(argvOf('tblastx', 'q.fa', 's.fa', []), 7)).not.toMatch(/num_threads/);
  });

  it('compares every option the search form writes, regions and negative values included', () => {
    const forms: Array<[ProgramId, string[]]> = [
      ['blastn', ['-task', 'dc-megablast', '-max_target_seqs', '5', '-evalue', '1e-5', '-word_size', '11', '-reward', '1', '-penalty', '-2', '-gapopen', '2', '-gapextend', '1', '-dust', 'no', '-lcase_masking', '-template_length', '18', '-template_type', 'coding', '-max_hsps', '2', '-perc_identity', '90', '-subject_besthit', '-query_loc', '1-50', '-subject_loc', '3-40']],
      ['blastp', ['-task', 'blastp-fast', '-max_target_seqs', '5', '-evalue', '1', '-word_size', '3', '-matrix', 'BLOSUM45', '-gapopen', '15', '-gapextend', '2', '-comp_based_stats', '0', '-seg', 'yes', '-threshold', '12', '-window_size', '30', '-max_hsps', '1', '-query_loc', '1-50', '-subject_loc', '3-40']],
      ['tblastn', ['-task', 'tblastn-fast', '-db_gencode', '1', '-max_target_seqs', '5', '-evalue', '1', '-word_size', '3', '-matrix', 'PAM30', '-gapopen', '9', '-gapextend', '1', '-comp_based_stats', '0', '-seg', 'no', '-soft_masking', 'true', '-lcase_masking', '-threshold', '13', '-window_size', '40', '-xdrop_gap', '20', '-xdrop_gap_final', '30', '-sum_stats', 'false', '-query_loc', '1-50', '-subject_loc', '3-40']],
      ['tblastx', ['-query_gencode', '2', '-db_gencode', '1', '-max_target_seqs', '5', '-evalue', '1', '-word_size', '2', '-culling_limit', '1', '-seg', '12 2.2 2.5', '-threshold', '14', '-window_size', '30', '-query_loc', '1-50', '-subject_loc', '3-40']],
    ];
    for (const [program, words] of forms) {
      expect(ncbiComparison(argvOf(program, 'q.fa', 's.fa', words)), program).toEqual({ refused: [], exceptions: [] });
    }
  });

  it('gives no NCBI command for a genetic code that NCBI BLAST+ refuses, and says why', () => {
    const comparison = ncbiComparison(argvOf('tblastn', 'q.faa', 's.fna', ['-db_gencode', '32']));
    expect(comparison.refused).toEqual([
      `NCBI BLAST+ ${NCBI_BLAST_VERSION} does not accept -db_gencode 32: its command line takes the genetic codes 1-6, 9-16, 21-31 and 33. ` +
        'LOSAT searched with it (PD-TLOSAN-LOCAL-GENCODE-32), so there is no NCBI command to compare with.',
    ]);
    // The approved exception of a non-default subject code is the one of the verification badge.
    expect(comparison.exceptions).toEqual(approvedExceptions('tblastn', argvOf('tblastn', 'q.faa', 's.fna', ['-db_gencode', '32'])));
    expect(comparison.exceptions[0]).toMatch(/^Approved exception \(PD-TLOSAN-LOCAL-GENCODE-32\)/);
    expect(ncbiComparison(argvOf('tblastx', 'q', 's', ['-query_gencode', '7'])).refused[0]).toMatch(/does not accept -query_gencode 7:.*LOSAT searched with it, so/);
  });

  it('gives the NCBI command with the approved exception of a non-default subject code of TBLASTN and TBLASTX', () => {
    for (const program of ['tblastn', 'tblastx'] as const) {
      const argv = argvOf(program, 'q', 's', ['-db_gencode', '4']);
      const comparison = ncbiComparison(argv);
      expect(comparison.refused).toEqual([]);
      expect(comparison.exceptions).toEqual(approvedExceptions(program, argv));
      expect(comparison.exceptions).toHaveLength(1);
    }
    expect(ncbiComparison(argvOf('tblastx', 'q', 's', ['-query_gencode', '4'])).exceptions).toEqual([]);
  });

  it('gives no NCBI command for an option it has not checked, or for arguments that are not options', () => {
    expect(ncbiComparison(argvOf('blastp', 'q', 's', ['-ungapped'])).refused).toEqual([
      `LOSAT Web has not checked -ungapped against BLASTP of NCBI BLAST+ ${NCBI_BLAST_VERSION}, so it gives no command to compare with.`,
    ]);
    expect(ncbiComparison(argvOf('blastx', 'q', 's', ['-query_loc', '1-9'])).refused).toHaveLength(1);
    expect(ncbiComparison(argvOf('blastn', 'q', 's', ['1e-5'])).refused).toEqual(['The run\'s arguments have "1e-5" where an option was expected.']);
    expect(ncbiComparison(argvOf('blastn', 'q', 's', ['-evalue'])).refused).toEqual(["The run's arguments end with -evalue without its value."]);
  });

  it('notes where the files go and the threads, and when two inputs of one name differ', () => {
    expect(commandNotes(argvOf('blastn', 'query.fa', 'combined_subject.fa', []), false)).toEqual([
      'The commands name the inputs as the run did. Put the files query.fa and combined_subject.fa in one folder and run the commands there; ' +
        'the browser does not know the folders of your files.',
      'The outputs do not depend on the threads, so the commands do not set -num_threads.',
    ]);
    expect(commandNotes(argvOf('blastn', 'a.fa', 'a.fa', []), true)[0]).toMatch(/^The commands name the input as the run did\. Put the file a\.fa in a folder/);
    expect(commandNotes(argvOf('blastn', 'a.fa', 'a.fa', []), false)[0]).toMatch(/^The query and the subject are both named a\.fa, but they differ/);
  });

  it('says how the input of a run relates to what was chosen', () => {
    expect(inputRelation('query', 'q.fa', 3, [{ origin: 'file', name: 'q.fa', records: 3 }])).toBe('q.fa has the same bytes as the file q.fa (3 records).');
    expect(inputRelation('subject', 'subject.fa', 1, [{ origin: 'paste', name: 'subject.fa', records: 1 }])).toBe(
      'subject.fa is the pasted subject text, as the run searched it (1 record).',
    );
    expect(inputRelation('query', 'q.fa', 2, [undefined])).toBe(
      'q.fa has the 2 records that the run searched; the records left out of the chosen query are not in it.',
    );
    expect(inputRelation('subject', 'combined_subject.fa', 5, [{ origin: 'file', name: 'a.fa', records: 2 }, undefined])).toBe(
      'combined_subject.fa joins the 2 subject inputs (a.fa) in the order chosen, without the records left out: 5 records. It is no single file you chose.',
    );
    expect(inputRelation('subject', 'combined_subject.fa', 4, [{ origin: 'file', name: 'a.fa', records: 2 }, { origin: 'paste', name: 'subject.fa', records: 2 }])).toBe(
      'combined_subject.fa joins the 2 subject inputs (a.fa, subject.fa) in the order chosen: 4 records. It is no single file you chose.',
    );
  });
});

/**
 * The argvs whose commands check_commands.py runs: each program, the defaults and options the
 * form writes (a region and a negative value among them), names as the application gives them
 * (pasted, joined, a file name with a space and a quote), and the genetic codes of the exceptions.
 * The script names the FASTA files of each case.
 */
const FIXED_CASES: ReadonlyArray<{ readonly id: string; readonly program: ProgramId; readonly query: string; readonly subject: string; readonly options: readonly string[] }> = [
  { id: 'blastn.default', program: 'blastn', query: 'query.fa', subject: 'subject.fa', options: [] },
  {
    id: 'blastn.options',
    program: 'blastn',
    query: 'multi query.fa',
    subject: "subject's.fa",
    options: ['-task', 'blastn', '-reward', '1', '-penalty', '-2', '-gapopen', '2', '-gapextend', '1', '-evalue', '1e-5', '-max_target_seqs', '5', '-max_hsps', '2', '-dust', 'no', '-lcase_masking', '-subject_besthit'],
  },
  { id: 'blastn.region', program: 'blastn', query: 'e2d_n_q.fa', subject: 'combined_subject.fa', options: ['-task', 'dc-megablast', '-query_loc', '1001-6000', '-subject_loc', '501-6600'] },
  { id: 'blastp.default', program: 'blastp', query: 'query.fa', subject: 'subject.fa', options: [] },
  {
    id: 'blastp.options',
    program: 'blastp',
    query: 'e2d_p_q.faa',
    subject: 'subject.fa',
    // LOSAT's BLASTP searches BLOSUM62 11/1 with -comp_based_stats 2 only (an explicit rejection otherwise).
    options: ['-task', 'blastp-fast', '-seg', 'yes', '-threshold', '12', '-window_size', '30', '-max_hsps', '1', '-evalue', '1e-3', '-query_loc', '100-400'],
  },
  { id: 'tblastn.default', program: 'tblastn', query: 'query.fa', subject: 'subject.fa', options: [] },
  {
    id: 'tblastn.options',
    program: 'tblastn',
    query: 'e2d_t_pq.faa',
    subject: 'e2d_t_ts.fa',
    // LOSAT's TBLASTN refuses -task tblastn-fast explicitly.
    options: ['-comp_based_stats', '0', '-seg', 'yes', '-soft_masking', 'true', '-lcase_masking', '-xdrop_gap', '20', '-xdrop_gap_final', '30', '-sum_stats', 'false', '-query_loc', '20-400', '-subject_loc', '42-700'],
  },
  { id: 'tblastn.gencode4', program: 'tblastn', query: 'query.fa', subject: 'subject.fa', options: ['-db_gencode', '4'] },
  { id: 'tblastn.gencode32', program: 'tblastn', query: 'query.fa', subject: 'subject.fa', options: ['-db_gencode', '32'] },
  { id: 'tblastx.default', program: 'tblastx', query: 'query.fa', subject: 'subject.fa', options: [] },
  {
    id: 'tblastx.options',
    program: 'tblastx',
    query: 'e2d_x_q.fa',
    subject: 'e2d_t_ts.fa',
    options: ['-query_gencode', '2', '-max_target_seqs', '2', '-culling_limit', '1', '-seg', 'no', '-threshold', '14', '-window_size', '30', '-evalue', '0.01', '-query_loc', '30-5000', '-subject_loc', '42-700'],
  },
  { id: 'tblastx.gencode5', program: 'tblastx', query: 'query.fa', subject: 'subject.fa', options: ['-db_gencode', '5'] },
  { id: 'blastx.default', program: 'blastx', query: 'query.fa', subject: 'subject.fa', options: [] },
  { id: 'blastx.options', program: 'blastx', query: 'query.fa', subject: 'subject.fa', options: ['-evalue', '1e-5', '-query_gencode', '11', '-max_target_seqs', '3', '-seg', 'no'] },
];
const FORMATS: readonly OutputFormat[] = [0, 6, 7];

describe('the commands of the fixed comparison cases', () => {
  const cases = FIXED_CASES.map((c) => {
    const argv = buildArgv({ program: c.program, queryName: c.query, subjectName: c.subject, parameters: pairs(c.options) });
    const comparison = ncbiComparison(argv);
    return {
      id: c.id,
      program: c.program,
      query: c.query,
      subject: c.subject,
      argv,
      refused: comparison.refused,
      exceptions: comparison.exceptions,
      commands: FORMATS.map((format) => ({ format, losat: losatCommand(argv, format), ncbi: ncbiCommand(argv, format) })),
    };
  });

  it('has an NCBI command for every case but the code NCBI refuses, and the exceptions of the subject codes', () => {
    expect(cases.filter((c) => c.refused.length > 0).map((c) => c.id)).toEqual(['tblastn.gencode32']);
    expect(cases.filter((c) => c.exceptions.length > 0).map((c) => c.id)).toEqual(['tblastn.gencode4', 'tblastn.gencode32', 'tblastx.gencode5']);
    for (const c of cases) for (const command of c.commands) expect(command.losat).toBe(`LOSAT ${command.ncbi}`);
  });

  it.runIf(process.env.LOSAT_WEB_COMMANDS_OUT !== undefined)('writes them for check_commands.py', () => {
    writeFileSync(process.env.LOSAT_WEB_COMMANDS_OUT!, `${JSON.stringify({ ncbi: NCBI_BLAST_VERSION, formats: FORMATS, cases }, null, 2)}\n`);
  });
});

/** The parameters of an option list as the form gives them: a flag alone, or a flag and its value. */
function pairs(words: readonly string[]): Array<readonly [string, string | true]> {
  const out: Array<readonly [string, string | true]> = [];
  for (let i = 0; i < words.length; i++) {
    const next = words[i + 1];
    if (next === undefined || /^-[A-Za-z_]/.test(next)) {
      out.push([words[i]!, true]);
    } else {
      out.push([words[i]!, next]);
      i++;
    }
  }
  return out;
}
