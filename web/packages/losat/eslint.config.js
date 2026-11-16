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
);
