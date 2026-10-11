// Dependency direction (docs/losat_web_gui_plan.md §3.3):
//   ui -> application -> ports <- infra, and domain depends on nothing.
// Only the composition root (src/main.ts, src/composition.ts) may import every layer.
// The site's service worker (src/sw) and the boot page's script (src/boot) run apart from the
// app (docs/web/release_manifest.md): sw may use the domain and its own files; boot the domain
// and sw/protocol.ts; infra's client of the worker sw/protocol.ts only. No other layer uses them.
import tseslint from 'typescript-eslint';
import pluginVue from 'eslint-plugin-vue';

const layer = (name) => [`**/${name}`, `**/${name}/**`];
// 'sw/*': the worker's files except protocol.ts, the contract between the worker and the pages
// (gitignore patterns: a file can be re-included only when its directory itself is not excluded).
const groupOf = (name) => {
  if (name === 'vue') return ['vue'];
  if (name === 'sw/*') return ['**/sw/**', '!**/sw/protocol', '!**/sw/protocol.ts'];
  return layer(name);
};
const restrict = (files, forbidden) => ({
  files,
  rules: {
    'no-restricted-imports': [
      'error',
      {
        patterns: forbidden.map((name) => ({
          group: groupOf(name),
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
  {
    files: ['**/*.vue'],
    languageOptions: { parserOptions: { parser: tseslint.parser } },
  },
  restrict(['src/domain/**'], ['application', 'ports', 'infra', 'ui', 'sw', 'boot', 'vue']),
  // The site's service worker runs apart from every page: the domain and its own files.
  restrict(['src/sw/**'], ['application', 'ports', 'infra', 'ui', 'boot', 'vue']),
  restrict(['src/ports/**'], ['application', 'infra', 'ui', 'sw', 'boot', 'vue']),
  restrict(['src/application/**'], ['infra', 'ui', 'sw', 'boot', 'vue']),
  restrict(['src/infra/**'], ['application', 'ui', 'sw/*', 'boot', 'vue']),
  restrict(['src/ui/**'], ['infra', 'sw', 'boot']),
  // The boot page's script runs before any release: the domain, and the worker's URL and scope.
  restrict(['src/boot/**'], ['application', 'ports', 'infra', 'ui', 'sw/*', 'vue']),
);
