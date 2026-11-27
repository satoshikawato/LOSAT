// After the whole run: the renaming of support/test.ts over every project's output directory, for a test
// whose worker died before its fixtures' teardown (the run has failed then; this keeps its evidence
// uploadable, support/artifact-names.ts).
import type { FullConfig } from '@playwright/test';
import { renameUnuploadable } from './artifact-names';

export default function globalTeardown(config: FullConfig): void {
  for (const directory of new Set(config.projects.map((project) => project.outputDir))) {
    const names = renameUnuploadable(directory);
    if (names.length > 0) console.error(`Renamed files whose names CI's artifact upload refuses, under ${directory}: ${names.join(', ')}`);
  }
}
