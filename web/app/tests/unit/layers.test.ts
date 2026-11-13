// Keeps the dependency direction of docs/losat_web_gui_plan.md §3.3 enforced: lints small
// probe modules placed (virtually) in each layer and expects the import rule to fire. The
// package `losat` (web/packages/losat, web/AGENTS.md rule 5) has the same layers and its own rule.
import { ESLint } from 'eslint';
import { fileURLToPath } from 'node:url';
import { beforeAll, describe, expect, it } from 'vitest';

const eslint = new ESLint({ cwd: fileURLToPath(new URL('../..', import.meta.url)) });
const packageEslint = new ESLint({ cwd: fileURLToPath(new URL('../../../packages/losat', import.meta.url)) });

async function importRuleFires(filePath: string, specifier: string, linter = eslint): Promise<boolean> {
  const code = `import * as probe from '${specifier}';\nexport const used = probe;\n`;
  const [result] = await linter.lintText(code, { filePath });
  return (result?.messages ?? []).some((message) => message.ruleId === 'no-restricted-imports');
}

describe('layer boundaries', () => {
  // The first lint loads the parsers and the config, which can take several seconds.
  beforeAll(async () => {
    await eslint.lintText('export {};\n', { filePath: 'src/domain/warmup.ts' });
    await packageEslint.lintText('export {};\n', { filePath: 'src/domain/warmup.ts' });
  }, 60_000);

  it.each([
    ['src/domain/probe.ts', '../application/coordinator'],
    ['src/domain/probe.ts', '../ports/data'],
    ['src/domain/probe.ts', '../infra/fake/fake-engine'],
    ['src/domain/probe.ts', 'vue'],
    ['src/domain/probe.ts', 'losat'],
    ['src/domain/probe.ts', 'losat/ports'],
    ['src/ports/probe.ts', 'losat'],
    ['src/ports/probe.ts', '../application/coordinator'],
    ['src/ports/probe.ts', '../infra/fake/fake-engine'],
    ['src/application/probe.ts', '../infra/fake/fake-engine'],
    ['src/application/probe.ts', '../ui/useStore'],
    ['src/application/probe.ts', 'vue'],
    ['src/application/probe.ts', 'losat'],
    ['src/infra/probe.ts', '../application/coordinator'],
    ['src/infra/probe.ts', 'vue'],
    ['src/ui/probe.ts', '../infra/fake/fake-engine'],
    ['src/ui/probe.ts', 'losat'],
  ])('%s may not import %s', async (filePath, specifier) => {
    expect(await importRuleFires(filePath, specifier)).toBe(true);
  });

  it.each([
    ['src/ports/probe.ts', '../domain/programs'],
    ['src/domain/probe.ts', 'losat/domain'],
    ['src/ports/probe.ts', 'losat/ports'],
    ['src/application/probe.ts', '../ports/data'],
    ['src/application/probe.ts', 'losat/ports'],
    ['src/infra/probe.ts', '../ports/data'],
    ['src/infra/probe.ts', 'losat'],
    ['src/ui/probe.ts', 'losat/domain'],
    ['src/ui/probe.ts', '../application/coordinator'],
    ['src/ui/probe.ts', '../domain/argv'],
    ['src/composition.ts', './infra/fake/fake-engine'],
    ['src/composition.ts', 'losat'],
  ])('%s may import %s', async (filePath, specifier) => {
    expect(await importRuleFires(filePath, specifier)).toBe(false);
  });

  it.each([
    ['src/domain/probe.ts', '../ports/engine'],
    ['src/domain/probe.ts', '../infra/digest'],
    ['src/domain/probe.ts', 'vue'],
    ['src/ports/probe.ts', '../infra/digest'],
    ['src/infra/probe.ts', 'vue'],
    ['src/infra/probe.ts', '../../../app/src/domain/programs'],
  ])('package: %s may not import %s', async (filePath, specifier) => {
    expect(await importRuleFires(filePath, specifier, packageEslint)).toBe(true);
  });

  it.each([
    ['src/ports/probe.ts', '../domain/programs'],
    ['src/infra/probe.ts', '../ports/engine'],
    ['src/infra/reactor/probe.ts', '../../domain/programs'],
  ])('package: %s may import %s', async (filePath, specifier) => {
    expect(await importRuleFires(filePath, specifier, packageEslint)).toBe(false);
  });
});
