// "Try an example" (S16, domain/examples.ts): the program's example query and subject are
// fetched from the site and put in the search draft as if the user had opened those files,
// with default settings. Nothing is searched: the user then chooses "Run LOSAT". The UI asks
// before an example replaces anything that the user entered (`replaces`): inputs, the
// program's settings, the Job Title.
import { EXAMPLE_TOPICS, exampleFileName } from '../domain/examples';
import { programById, type InputRole, type ProgramId } from '../domain/programs';
import type { ExampleSource } from '../ports/examples';
import type { DraftState, SearchDraft } from './draft';
import { Store } from './store';

export interface ExamplesState {
  /** The program whose example is being fetched. */
  readonly loading?: ProgramId;
  /** The example put in the draft last, and what it searches. */
  readonly loaded?: { readonly program: ProgramId; readonly text: string };
  /** Why the latest example could not be loaded. */
  readonly error?: string;
}

export interface ExamplesDeps {
  readonly source: ExampleSource;
  readonly draft: Pick<SearchDraft, 'state' | 'loadExample'>;
}

const ROLES: readonly InputRole[] = ['query', 'subject'];

/** What the user entered that loading an example would replace (SearchDraft.loadExample). */
export interface Replaced {
  /** Pasted text or files of either role. */
  readonly inputs: boolean;
  /** The program's fields, the threads or how the records of a role are searched. */
  readonly settings: boolean;
  readonly title: boolean;
}

/** Whether the draft's inputs are now exactly the program's example files (the user changed nothing of them). */
export function holdsExample(draft: DraftState, program: ProgramId): boolean {
  return (
    draft.program === program &&
    ROLES.every((role) => {
      const { paste, sources } = draft[role];
      return paste === '' && sources.length === 1 && sources[0]!.origin === 'file' && sources[0]!.name === exampleFileName(program, role);
    })
  );
}

export class Examples {
  readonly state = new Store<ExamplesState>({});

  constructor(private readonly deps: ExamplesDeps) {}

  /** Why the program's example cannot be tried, or undefined: a program that cannot be searched yet (BLASTX until SX). */
  unavailable(program: ProgramId): string | undefined {
    return programById(program).unavailable;
  }

  /** What loading `program`'s example would replace that the user entered; undefined when nothing. */
  replaces(program: ProgramId): Replaced | undefined {
    const draft = this.deps.draft.state.get();
    const replaced: Replaced = {
      inputs: ROLES.some((role) => draft[role].paste !== '' || draft[role].sources.length > 0),
      settings: Object.keys(draft.values[program]).length > 0 || draft.threads !== 'auto' || ROLES.some((role) => draft[role].mode !== 'combined'),
      title: draft.title !== '',
    };
    return replaced.inputs || replaced.settings || replaced.title ? replaced : undefined;
  }

  /**
   * Fetches the program's example and puts it in the draft, replacing its inputs, the program's
   * settings and the Job Title (the UI asked first if the user had entered any). A failed fetch
   * changes nothing in the draft and says why.
   */
  async load(program: ProgramId): Promise<void> {
    const unavailable = this.unavailable(program);
    if (unavailable !== undefined) {
      this.state.set({ error: unavailable });
      return;
    }
    if (this.state.get().loading !== undefined) return;
    this.state.set({ loading: program });
    try {
      const files = await this.deps.source.files(program);
      this.deps.draft.loadExample(program, files);
      const label = programById(program).label;
      this.state.set({
        loaded: {
          program,
          text: `Example loaded: ${label} of ${files.query.name} against ${files.subject.name}, ${EXAMPLE_TOPICS[program]}. Choose Run LOSAT to search.`,
        },
      });
    } catch (error) {
      this.state.set({ error: error instanceof Error ? error.message : String(error) });
    }
  }
}
