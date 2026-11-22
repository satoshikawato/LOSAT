// An input that holds no sequence (an empty file, white space only, a defline alone) is read and run as
// the command line runs it, but its card shows no "Accepted" check mark: it says that the input adds
// nothing to the search and that the run may end with an error (line WQA's sweep, D4, 2026-11-19). NCBI
// BLAST+ 2.17.0 warns of such an input or ends the run with an error, by role, program and the other
// inputs; an input with an empty record among records with residues runs with a warning, and keeps the
// mark. One whose records with residues are left out says so (D4 review).
import { expect, type Page } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { openFiles, openSearch, paste, pickRun, program, showRecords, submit, waitStatus, type Role } from './support/search';

test.setTimeout(180_000);

test.beforeEach(async ({ page }) => {
  await openSearch(page);
});

const NO_SEQUENCE = 'This input has no sequence, so it adds nothing to the search, and the run may end with an error. Add a sequence or remove it.';
const EXCLUDED =
  'Only records without a sequence are included, so this input adds nothing to the search, and the run may end with an error. Include a record that has a sequence, or remove this input.';
const SEQUENCE = '>s1\nACGTACGTACGTAAACCCGGGTTTACGTACGTACGTAAACCCGGGTTA\n';
/** The inputs, and whether each holds no sequence. */
const INPUTS = [
  { name: 'zero.fa', text: '', none: true },
  { name: 'whitespace.fa', text: '   \n\n\t\n', none: true },
  { name: 'defline-only.fa', text: '>empty record\n', none: true },
  { name: 'one-empty.fa', text: `>a\nACGTACGTACGTACGTACGTACGTAAACCCGGGTTT\n>empty record\n>b\nTTTGGGCCCAAACGTACGTACGTACGTACGT\n`, none: false },
] as const;

/** The check of source `index` of a role: the warning without the check mark, or the check mark. */
async function expectCheck(page: Page, role: Role, index: number, none: boolean): Promise<void> {
  const check = page.getByTestId(`${role}-source-${index}-check`);
  await expect(check).toHaveAttribute('data-check', 'ok');
  if (none) {
    await expect(check).toHaveAttribute('data-sequence', 'none');
    await expect(check).toHaveText(NO_SEQUENCE);
    await expect(check.locator('svg')).toHaveCount(0);
    await expect(check).not.toContainText('Accepted');
  } else {
    await expect(check).not.toHaveAttribute('data-sequence');
    await expect(check.locator('svg')).toHaveCount(1);
    await expect(check).toContainText('Accepted by the BLASTN engine: 3 records');
  }
}

for (const role of ['query', 'subject'] as const) {
  test(`a ${role} that holds no sequence has no check mark, and runs as the command line runs it`, async ({ page }) => {
    const other: Role = role === 'query' ? 'subject' : 'query';
    await program(page, 'blastn');
    await paste(page, other, SEQUENCE);
    await openFiles(page, role, INPUTS);
    for (const [index, input] of INPUTS.entries()) {
      await expect(page.getByTestId(`${role}-source-${index}`)).toContainText(input.name);
      await expectCheck(page, role, index, input.none);
    }
    // one-empty.fa with its two records with residues left out: only its empty record is included.
    const mixed = INPUTS.findIndex((input) => input.name === 'one-empty.fa');
    await showRecords(page, role, mixed);
    await page.getByTestId(`${role}-source-${mixed}-record-0`).uncheck();
    await page.getByTestId(`${role}-source-${mixed}-record-2`).uncheck();
    const check = page.getByTestId(`${role}-source-${mixed}-check`);
    await expect(check).toHaveAttribute('data-sequence', 'excluded');
    await expect(check).toHaveText(EXCLUDED);
    await expect(check.locator('svg')).toHaveCount(0);
    await page.getByTestId(`${role}-source-${mixed}-record-0`).check();
    await expect(check).not.toHaveAttribute('data-sequence');
    await expect(check).toContainText('Accepted by the BLASTN engine: 2 records');
    // A defline alone, searched alone: the form does not refuse it, and the run ends as the command
    // line's does (NCBI BLAST+ 2.17.0, exit status 3).
    for (let i = INPUTS.length - 1; i >= 0; i--) {
      if (INPUTS[i]!.name !== 'defline-only.fa') await page.getByTestId(`${role}-source-${i}-remove`).click();
    }
    await expect(page.locator(`[data-testid^="${role}-source-"][data-status]`)).toHaveCount(1);
    await expectCheck(page, role, 0, true);
    await submit(page);
    await expect(page.getByTestId('run-1-option')).toHaveAttribute('data-status', /^(completed|failed)$/, { timeout: 120_000 });
    if (BUILD_HAS_ENGINE) {
      await waitStatus(page, 1, 'failed');
      await pickRun(page, 1);
      await expect(page.getByTestId('results-status')).toContainText(
        role === 'query' ? 'LOSAT engine error: Warning: Sequence contains no data' : 'LOSAT engine error: The average subject length is too short',
      );
    }
  });
}
