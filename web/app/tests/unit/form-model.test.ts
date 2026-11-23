import { describe, expect, it, vi } from 'vitest';
import { toProgramDescription } from 'losat';
import DESCRIBE from '../../src/infra/fake/describe.json';
import { SearchDraft } from '../../src/application/draft';
import { normalizeBooleanDefault, writtenValue } from '../../src/domain/parameters';
import type { ParameterField } from '../../src/domain/programs';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { DataService } from '../../src/infra/data/data-service';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import { searchFormModel, type FormModel } from '../../src/ui/form-model';
import { REFUSAL_ID } from '../../src/ui/refusal';

function makeDraft() {
  const engine = new FakeEngine();
  const data = new DataService({
    store: new MemoryBlockStore(),
    scanner: new FakeScanner(),
    inflater: { open: async () => undefined },
    checker: new FakeInputChecker(),
    digest: async () => '0'.repeat(64),
    newToken: () => crypto.randomUUID(),
    cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
  });
  const draft = new SearchDraft({
    engine,
    data,
    enqueueAll: async () => ({ ok: true, runIds: [] }),
    fetcher: { fetch: async () => ({ ok: false, reason: 'not-found' }) },
    now: () => 0,
    debounceMs: 0,
  });
  return draft;
}

describe('searchFormModel', () => {
  it('exposes search sections, known defaults, pairs, taxid field and refusal through FormModel', async () => {
    const draft = makeDraft();
    await vi.waitFor(() => expect(draft.state.get().description).toBeDefined());
    const model = searchFormModel(draft, draft.state.get());
    expect(model.key).toBe('blastn');
    expect(model.idPrefix).toBe('param');
    expect(model.refusalId).toBe(REFUSAL_ID);
    expect(model.sections?.map((s) => s.title)).toContain('General Parameters');
    expect(model.knownDefault('-word_size')).toEqual({ defaultValue: '28' });
    expect(model.taxidFieldOf?.('-window_masker_taxid')?.flag).toBe('-window_masker_taxid');

    const scoring = model.sections?.find((s) => s.title === 'Scoring Parameters');
    const scoreFields = scoring?.fields.filter((f) => f.row === 'Match/Mismatch Scores') ?? [];
    const pair = model.pairOf?.('Match/Mismatch Scores', scoreFields);
    expect(pair?.select.id).toBe('matchscores');
    expect(model.pairValue?.(pair!)).toBe('1,-2');
    model.setPair?.(pair!, '1,-3');
    expect(searchFormModel(draft, draft.state.get()).pairValue?.(pair!)).toBe('1,-3');

    model.setExtra?.('-unknown_flag');
    await vi.waitFor(() => expect(draft.state.get().validation.state).toBe('invalid'));
    const refused = searchFormModel(draft, draft.state.get());
    expect(refused.refusal?.extra).toBe(true);
  });

  it('provides genetic codes for translated programs', async () => {
    const draft = makeDraft();
    draft.setProgram('tblastx');
    await vi.waitFor(() => expect(draft.state.get().description?.program).toBe('tblastx'));
    const tblastxModel = searchFormModel(draft, draft.state.get());
    expect(tblastxModel.queryGencodes).toEqual(toProgramDescription(DESCRIBE.tblastx).queryGencodes);
    expect(tblastxModel.subjectGencodes).toEqual(toProgramDescription(DESCRIBE.tblastx).subjectGencodes);
  });
});

describe('FormModel boolean select and tool-shaped models', () => {
  it('normalizes NCBI boolean defaults (F, T, false) and omits defaults when writtenValue checks them', () => {
    expect(normalizeBooleanDefault('F')).toBe('false');
    expect(normalizeBooleanDefault('T')).toBe('true');
    expect(normalizeBooleanDefault('false')).toBe('false');
    expect(normalizeBooleanDefault('true')).toBe('true');
    expect(normalizeBooleanDefault(undefined)).toBeUndefined();

    const dustField: ParameterField = { flag: '-dust', label: '', kind: 'boolean', select: true };
    expect(writtenValue(dustField, '', { defaultValue: 'F' })).toBeUndefined();
    expect(writtenValue(dustField, 'false', { defaultValue: 'F' })).toBeUndefined();
    expect(writtenValue(dustField, 'true', { defaultValue: 'F' })).toBe('true');

    const textMatchField: ParameterField = { flag: '-text_match', label: '', kind: 'boolean', select: true };
    expect(writtenValue(textMatchField, 'true', { defaultValue: 'T' })).toBeUndefined();
    expect(writtenValue(textMatchField, 'false', { defaultValue: 'T' })).toBe('false');
  });

  it('supports a minimal non-search FormModel without search-only props', () => {
    const values: Record<string, string | boolean> = {};
    const files: Record<string, File> = {};
    const toolModel: FormModel = {
      key: 'windowmasker:mk_counts',
      values,
      options: [{ flag: '-unit', help: 'number of bases in a unit', takesValue: true }],
      files,
      knownDefault: () => undefined,
      fieldList: () => [],
      refusal: undefined,
      idPrefix: 'wm-param',
      refusalId: 'wm-argv-validation-message',
      setField: (flag, value) => {
        values[flag] = value;
      },
      setFile: (flag, file) => {
        if (file === undefined) delete files[flag];
        else files[flag] = file;
      },
    };
    toolModel.setField('-unit', '15');
    expect(toolModel.values['-unit']).toBe('15');
    expect(toolModel.sections).toBeUndefined();
    expect(toolModel.pairOf).toBeUndefined();
  });
});
