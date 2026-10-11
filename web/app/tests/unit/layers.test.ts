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
    ['src/domain/probe.ts', '../sw/protocol'],
    ['src/domain/probe.ts', '../boot/plan'],
    ['src/sw/probe.ts', '../application/release'],
    ['src/sw/probe.ts', '../ports/offline'],
    ['src/sw/probe.ts', '../infra/browser/offline'],
    ['src/sw/probe.ts', '../boot/plan'],
    ['src/sw/probe.ts', 'vue'],
    ['src/ports/probe.ts', '../application/coordinator'],
    ['src/ports/probe.ts', '../infra/fake/fake-engine'],
    ['src/ports/probe.ts', '../sw/protocol'],
    ['src/ports/probe.ts', '../boot/plan'],
    ['src/application/probe.ts', '../infra/fake/fake-engine'],
    ['src/application/probe.ts', '../ui/useStore'],
    ['src/application/probe.ts', 'vue'],
    ['src/application/probe.ts', '../sw/protocol'],
    ['src/application/probe.ts', '../boot/plan'],
    ['src/infra/probe.ts', '../application/coordinator'],
    ['src/infra/probe.ts', 'vue'],
    ['src/infra/browser/probe.ts', '../../sw/prepare'],
    ['src/infra/browser/probe.ts', '../../sw/cleanup'],
    ['src/infra/browser/probe.ts', '../../boot/plan'],
    ['src/ui/probe.ts', '../infra/fake/fake-engine'],
    ['src/ui/probe.ts', '../sw/protocol'],
    ['src/ui/probe.ts', '../boot/plan'],
    ['src/boot/probe.ts', '../application/coordinator'],
    ['src/boot/probe.ts', '../ports/site'],
    ['src/boot/probe.ts', '../infra/browser/site'],
    ['src/boot/probe.ts', '../ui/useStore'],
    ['src/boot/probe.ts', 'vue'],
    ['src/boot/probe.ts', '../sw/prepare'],
    ['src/boot/probe.ts', '../sw/serve'],
  ])('%s may not import %s', async (filePath, specifier) => {
    expect(await importRuleFires(filePath, specifier)).toBe(true);
  });

  it.each([
    ['src/ports/probe.ts', '../domain/programs'],
    ['src/sw/probe.ts', '../domain/site-index'],
    ['src/application/probe.ts', '../ports/engine'],
    ['src/infra/probe.ts', '../ports/engine'],
    ['src/ui/probe.ts', '../application/coordinator'],
    ['src/ui/probe.ts', '../domain/argv'],
    ['src/composition.ts', './infra/fake/fake-engine'],
    ['src/boot/probe.ts', '../domain/release-manifest'],
    ['src/boot/probe.ts', '../sw/protocol'],
    ['src/sw/probe.ts', './prepare'],
    ['src/infra/browser/probe.ts', '../../sw/protocol'],
  ])('%s may import %s', async (filePath, specifier) => {
    expect(await importRuleFires(filePath, specifier)).toBe(false);
  });
});
