// The example files from the site (ports/examples.ts). They are imported as URLs that are
// never inlined, so the build emits each as a hashed asset of the release (and the release
// manifest lists it), and they are fetched from the site like the engine modules.
import { exampleFileName } from '../../domain/examples';
import type { InputRole, ProgramId } from '../../domain/programs';
import type { ExampleFiles, ExampleSource } from '../../ports/examples';
import blastnQuery from '../../examples/example-blastn-query.fasta?url&no-inline';
import blastnSubject from '../../examples/example-blastn-subject.fasta?url&no-inline';
import blastpQuery from '../../examples/example-blastp-query.fasta?url&no-inline';
import blastpSubject from '../../examples/example-blastp-subject.fasta?url&no-inline';
import blastxQuery from '../../examples/example-blastx-query.fasta?url&no-inline';
import blastxSubject from '../../examples/example-blastx-subject.fasta?url&no-inline';
import tblastnQuery from '../../examples/example-tblastn-query.fasta?url&no-inline';
import tblastnSubject from '../../examples/example-tblastn-subject.fasta?url&no-inline';
import tblastxQuery from '../../examples/example-tblastx-query.fasta?url&no-inline';
import tblastxSubject from '../../examples/example-tblastx-subject.fasta?url&no-inline';

const URLS: Readonly<Record<ProgramId, Readonly<Record<InputRole, string>>>> = {
  blastn: { query: blastnQuery, subject: blastnSubject },
  blastp: { query: blastpQuery, subject: blastpSubject },
  blastx: { query: blastxQuery, subject: blastxSubject },
  tblastn: { query: tblastnQuery, subject: tblastnSubject },
  tblastx: { query: tblastxQuery, subject: tblastxSubject },
};
const TIMEOUT_MS = 30_000;

async function fetchFile(program: ProgramId, role: InputRole): Promise<File> {
  const name = exampleFileName(program, role);
  let response: Response;
  try {
    response = await fetch(URLS[program][role], { credentials: 'same-origin', signal: AbortSignal.timeout(TIMEOUT_MS) });
  } catch (error) {
    throw new Error(`The example file ${name} could not be read from the site: ${error instanceof Error ? error.message : String(error)}`);
  }
  if (!response.ok) throw new Error(`The example file ${name} could not be read from the site (HTTP ${response.status}).`);
  return new File([await response.blob()], name, { type: 'text/plain' });
}

export const browserExamples: ExampleSource = {
  async files(program): Promise<ExampleFiles> {
    const [query, subject] = await Promise.all([fetchFile(program, 'query'), fetchFile(program, 'subject')]);
    return { query, subject };
  },
};
