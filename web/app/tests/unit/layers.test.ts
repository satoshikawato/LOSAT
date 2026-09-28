// Keeps the dependency direction of docs/losat_web_gui_plan.md §3.3 enforced: lints small
// probe modules placed (virtually) in each layer and expects the import rule to fire.
import { ESLint } from 'eslint';
import { fileURLToPath } from 'node:url';
import { beforeAll, describe, expect, it } from 'vitest';

const eslint = new ESLint({ cwd: fileURLToPath(new URL('../..', import.meta.url)) });

async function importRuleFires(filePath: string, specifier: string): Promise<boolean> {
  const code = `import * as probe from '${specifier}';\nexport const used = probe;\n`;
  const [result] = await eslint.lintText(code, { filePath });
  return (result?.messages ?? []).some((message) => message.ruleId === 'no-restricted-imports');
}

describe('layer boundaries', () => {
  // The first lint loads the parsers and the config, which can take several seconds.
  beforeAll(async () => {
    await eslint.lintText('export {};\n', { filePath: 'src/domain/warmup.ts' });
  }, 60_000);

  it.each([
    ['src/domain/probe.ts', '../application/coordinator'],
    ['src/domain/probe.ts', '../ports/engine'],
    ['src/domain/probe.ts', '../infra/fake/fake-engine'],
    ['src/domain/probe.ts', 'vue'],
    ['src/ports/probe.ts', '../application/coordinator'],
    ['src/ports/probe.ts', '../infra/fake/fake-engine'],
    ['src/application/probe.ts', '../infra/fake/fake-engine'],
    ['src/application/probe.ts', '../ui/useStore'],
    ['src/application/probe.ts', 'vue'],
    ['src/infra/probe.ts', '../application/coordinator'],
    ['src/infra/probe.ts', 'vue'],
    ['src/ui/probe.ts', '../infra/fake/fake-engine'],
  ])('%s may not import %s', async (filePath, specifier) => {
    expect(await importRuleFires(filePath, specifier)).toBe(true);
  });

  it.each([
    ['src/ports/probe.ts', '../domain/programs'],
    ['src/application/probe.ts', '../ports/engine'],
    ['src/infra/probe.ts', '../ports/engine'],
    ['src/ui/probe.ts', '../application/coordinator'],
    ['src/ui/probe.ts', '../domain/argv'],
    ['src/composition.ts', './infra/fake/fake-engine'],
  ])('%s may import %s', async (filePath, specifier) => {
    expect(await importRuleFires(filePath, specifier)).toBe(false);
  });
});
