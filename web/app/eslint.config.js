// Dependency direction (docs/losat_web_gui_plan.md §3.3):
//   ui -> application -> ports <- infra, and domain depends on nothing.
// Only the composition root (src/main.ts, src/composition.ts) may import every layer.
import tseslint from 'typescript-eslint';
import pluginVue from 'eslint-plugin-vue';

const layer = (name) => [`**/${name}`, `**/${name}/**`];
const restrict = (files, forbidden) => ({
  files,
  rules: {
    'no-restricted-imports': [
      'error',
      {
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
  {
    files: ['**/*.vue'],
    languageOptions: { parserOptions: { parser: tseslint.parser } },
  },
  restrict(['src/domain/**'], ['application', 'ports', 'infra', 'ui', 'vue']),
  restrict(['src/ports/**'], ['application', 'infra', 'ui', 'vue']),
  restrict(['src/application/**'], ['infra', 'ui', 'vue']),
  restrict(['src/infra/**'], ['application', 'ui', 'vue']),
  restrict(['src/ui/**'], ['infra']),
);
