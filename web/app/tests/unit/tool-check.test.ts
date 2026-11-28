// The check of a tool tab's words (`ToolDraftCore`: debounced, one at a time, an older answer dropped) ends with the
// answer for the form as it is when the engine answers late: the page's first check waits for the Engine worker's
// start (its two files fetched and compiled), and the edits made meanwhile are checked after it (line WW, WW4j). A
// check stays `checking` only while the engine has not answered the newest one.
import { describe, expect, it, vi } from 'vitest';
import type { DirectoryFile } from 'losat/domain';
import type { ToolValidation } from 'losat/ports';
import { ToolDraft } from '../../src/application/tool-draft';
import { FakeEngine } from '../../src/infra/fake/fake-engine';

interface Asked {
  readonly argv: readonly string[];
  readonly names: readonly string[];
  readonly answer: (verdict: ToolValidation) => void;
}

/** A WindowMasker form whose engine answers each check only when the test says so. */
function heldDraft() {
  const fake = new FakeEngine();
  const asked: Asked[] = [];
  const draft = new ToolDraft({
    engine: {
      describeTool: (tool) => fake.describeTool(tool),
      validateTool: (argv: readonly string[], options: { readonly directory?: readonly DirectoryFile[] } = {}) =>
        new Promise<ToolValidation>((resolve) => asked.push({ argv, names: (options.directory ?? []).map((entry) => entry.name), answer: resolve })),
    },
    data: {
      inflateToSessionFile: () => Promise.reject(new Error('not used')),
      dropSessionFile: () => Promise.resolve(),
      getToolOutputFile: () => Promise.reject(new Error('not used')),
    },
    enqueueToolAll: () => Promise.resolve({ ok: true, runIds: [] }),
    debounceMs: 0,
  });
  draft.ensureLoaded();
  return { draft, asked };
}

const folder = (draft: ToolDraft) =>
  draft.addRunFiles(['nin', 'nsq', 'nhr'].map((ext, i) => new File([new Uint8Array([i])], `wm.${ext}`)), ['my db/wm.nin', 'my db/wm.nsq', 'my db/wm.nhr']);

/** The form of seqdb-names.spec.ts:344 (WindowMasker -infmt blastdb on the folder "my db"), made while the first check waits. */
async function fillWhileFirstCheckWaits(draft: ToolDraft, asked: Asked[]): Promise<void> {
  await vi.waitFor(() => expect(asked).toHaveLength(1));
  expect(draft.state.get().validation).toEqual({ state: 'checking' });
  draft.setField('-infmt', 'blastdb');
  folder(draft);
  draft.setField('-sformat', 'binary');
  draft.setOutName('-');
}

/** Its argv (the form writes no `-out` for `-`, standard output). */
const LAST = ['windowmasker', '-mk_counts', '-in', '"my db/wm"', '-infmt', 'blastdb', '-sformat', 'binary'];

describe('the check of a tool tab whose engine answers late', () => {
  it('answers for the form as it is now once the first check (the engine start) answers, the edits made meanwhile checked after it', async () => {
    const { draft, asked } = heldDraft();
    await fillWhileFirstCheckWaits(draft, asked);
    // One at a time: nothing else is asked while the first check waits, and the form shows that it is checking.
    await new Promise((resolve) => setTimeout(resolve, 20));
    expect(asked).toHaveLength(1);
    expect(asked[0]!.argv).toEqual(['windowmasker', '-mk_counts', '-in', 'input.fa', '-out', 'wmasker.ascii']);
    expect(draft.state.get().validation).toEqual({ state: 'checking' });
    // The engine has started: its answer to the first form is dropped, and the form as it is now is checked.
    asked[0]!.answer({ ok: false, message: 'the first form' });
    await vi.waitFor(() => expect(asked).toHaveLength(2));
    expect(asked[1]!.argv).toEqual(LAST);
    expect(asked[1]!.names).toEqual(['my db/wm.nin', 'my db/wm.nsq', 'my db/wm.nhr']);
    expect(draft.state.get().validation).toEqual({ state: 'checking' });
    asked[1]!.answer({ ok: true });
    await vi.waitFor(() => expect(draft.state.get().validation).toEqual({ state: 'ok' }));
    expect(asked).toHaveLength(2);
  });

  it('asks again after each answer while edits keep coming, and shows only the answer for the last form', async () => {
    const { draft, asked } = heldDraft();
    await vi.waitFor(() => expect(asked).toHaveLength(1));
    draft.setField('-infmt', 'blastdb');
    asked[0]!.answer({ ok: true });
    await vi.waitFor(() => expect(asked).toHaveLength(2));
    // An edit while the second check waits: its answer (for -infmt blastdb with no -in yet) is dropped.
    folder(draft);
    asked[1]!.answer({ ok: false, message: 'no -in' });
    await vi.waitFor(() => expect(asked).toHaveLength(3));
    expect(draft.state.get().validation).toEqual({ state: 'checking' });
    expect(asked[2]!.argv).toContain('"my db/wm"');
    asked[2]!.answer({ ok: false, message: 'refused' });
    await vi.waitFor(() => expect(draft.state.get().validation).toEqual({ state: 'invalid', message: 'refused' }));
  });

  it('lets Run LOSAT wait for the answer to the form as it is now, and queues nothing for a refused form', async () => {
    const { draft, asked } = heldDraft();
    await fillWhileFirstCheckWaits(draft, asked);
    const submitted = draft.submit();
    asked[0]!.answer({ ok: true });
    await vi.waitFor(() => expect(asked).toHaveLength(2));
    expect(asked[1]!.argv).toEqual(LAST);
    asked[1]!.answer({ ok: false, message: 'Error: refused' });
    expect(await submitted).toEqual({ ok: false, message: 'Error: refused' });
  });
});
