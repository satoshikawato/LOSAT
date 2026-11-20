// Dependency direction (web/AGENTS.md rule 5), as in the application:
//   infra -> ports -> domain, and domain depends on nothing.
// The package imports nothing of the application and no UI framework.
import tseslint from 'typescript-eslint';

const layer = (name) => [`**/${name}`, `**/${name}/**`];
const restrict = (files, forbidden) => ({
  files,
  rules: {
    'no-restricted-imports': [
      'error',
      {
        patterns: [
          ...forbidden.map((name) => ({
            group: layer(name),
            message: `This layer must not depend on "${name}" (see web/AGENTS.md rule 5).`,
          })),
          { group: ['**/app', '**/app/**', 'vue'], message: 'The package must not depend on the application or Vue.' },
        ],
      },
    ],
  },
});

export default tseslint.config(
  { ignores: ['dist/**', 'node_modules/**'] },
  ...tseslint.configs.recommended,
  // The root of this configuration's tsconfig; the application's unit tests lint with both configurations.
  { languageOptions: { parserOptions: { tsconfigRootDir: import.meta.dirname } } },
  restrict(['src/domain/**'], ['ports', 'infra']),
  restrict(['src/ports/**'], ['infra']),
  restrict(['src/infra/**'], []),
  // Firefox refuses `new BroadcastChannel()` where the site's data is blocked, on the page and in every worker, and every
  // search failed (line WQA's sweep, D8, 2026-11-20): the host links its workers with MessageChannels (link.ts). Only the
  // test hooks of test builds open one.
  {
    files: ['src/**'],
    rules: {
      'no-restricted-syntax': [
        'error',
        {
          selector:
            "NewExpression[callee.name='BroadcastChannel']:not(IfStatement[test.name='__LOSAT_TEST_HOOKS__'] *, ConditionalExpression[test.name='__LOSAT_TEST_HOOKS__'] *)",
          message: 'The browser refuses a BroadcastChannel where the site\'s data is blocked: link workers with a MessageChannel (src/infra/engine-worker/link.ts).',
        },
      ],
    },
  },
);
