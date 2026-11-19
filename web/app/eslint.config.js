// Dependency direction (docs/losat_web_gui_plan.md §3.3):
//   ui -> application -> ports <- infra, and domain depends on nothing.
// Only the composition root (src/main.ts, src/composition.ts) may import every layer.
// The package `losat` (web/packages/losat, the engine host) has the same layers: `losat/domain`
// is a domain that every layer may import, `losat/ports` a port, and `losat` itself (its infra)
// only infra and the composition root may import.
import tseslint from 'typescript-eslint';
import pluginVue from 'eslint-plugin-vue';

const layer = (name) => [`**/${name}`, `**/${name}/**`];
const restrict = (files, forbidden) => ({
  files,
  rules: {
    'no-restricted-imports': [
      'error',
      {
        // The package's infra (`losat`) belongs with infra; `**/ports` also matches `losat/ports`.
        paths: forbidden.includes('infra')
          ? [{ name: 'losat', message: 'This layer must not depend on "infra" (the package\'s is `losat`; see docs/losat_web_gui_plan.md §3.3).' }]
          : [],
        patterns: forbidden.map((name) => ({
          group: name === 'vue' ? ['vue'] : layer(name),
          message: `This layer must not depend on "${name}" (see docs/losat_web_gui_plan.md §3.3).`,
        })),
      },
    ],
  },
});

export default tseslint.config(
  { ignores: ['dist/**', 'node_modules/**', 'test-results/**', 'playwright-report/**'] },
  ...tseslint.configs.recommended,
  ...pluginVue.configs['flat/essential'],
  // The root of this configuration's tsconfig; the application's unit tests lint with both configurations.
  { languageOptions: { parserOptions: { tsconfigRootDir: import.meta.dirname } } },
  {
    files: ['**/*.vue'],
    languageOptions: { parserOptions: { parser: tseslint.parser } },
  },
  restrict(['src/domain/**'], ['application', 'ports', 'infra', 'ui', 'vue']),
  restrict(['src/ports/**'], ['application', 'infra', 'ui', 'vue']),
  restrict(['src/application/**'], ['infra', 'ui', 'vue']),
  restrict(['src/infra/**'], ['application', 'ui', 'vue']),
  restrict(['src/ui/**'], ['infra']),
  // Firefox refuses `new BroadcastChannel()` where the site's data is blocked, and the page did not start (line WQA's
  // sweep, D5, 2026-11-19): the application opens its channels with siteChannel, which tells no other tab then.
  {
    files: ['src/**'],
    ignores: ['src/infra/browser/site-channel.ts'],
    rules: {
      'no-restricted-syntax': [
        'error',
        {
          selector: "NewExpression[callee.name='BroadcastChannel']",
          message: "The browser refuses a BroadcastChannel where the site's data is blocked: open it with siteChannel (src/infra/browser/site-channel.ts).",
        },
      ],
    },
  },
  // The e2e tests take `test` from tests/e2e/support/test.ts, which puts the WebKit pages'
  // wake lock stand-in in place.
  {
    files: ['tests/e2e/**'],
    ignores: ['tests/e2e/support/test.ts'],
    rules: {
      'no-restricted-imports': [
        'error',
        {
          paths: [
            {
              name: '@playwright/test',
              importNames: ['test'],
              message: 'Import test from tests/e2e/support/test.ts.',
            },
          ],
        },
      ],
    },
  },
);
