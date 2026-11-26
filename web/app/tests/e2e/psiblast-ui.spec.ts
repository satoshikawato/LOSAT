// LOSAT Web PSIBLAST & TBLASTN -in_pssm UI, PSSM download, and responsive viewport E2E tests (PSW1).
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect } from '@playwright/test';
import { axeViolations } from './support/axe';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { openParameters, openSearch, paste, pickRun, program, submit, waitStatus } from './support/search';
import { test } from './support/test';

const MSA_FIXTURE = 'LOSAT/tests/fixtures/msa/ps6_b_db_0.aln';
const PSSM_FIXTURE = 'LOSAT/tests/fixtures/pssm/q1b_gi.asn';
const TBLASTN_SUBJECT_FIXTURE = 'LOSAT/tests/fasta/psiblast/q1b_gi_nt.fna';
const PSIBLAST_QUERY_FIXTURE = 'LOSAT/tests/fasta/psiblast/q1b_gi.faa';
const PSIBLAST_SUBJ_IDS_FIXTURE = 'LOSAT/tests/fasta/psiblast/q1_subj_ids.faa';

const PSIBLAST_SUBJECT_FASTA = [
  '>s1',
  'TMAGRVELVTGPMFAGKSTYLKNIYQQENGGNKHCLFVKHSLETRYGCGTGTIVTHAGEV',
  '>s2',
  'TNIGKIHLVIGPMFAGKSSTIIKIYQKKIKMQKRCVFIKHVKDNRYGCELTETVTHSGQT',
  '',
].join('\n');

async function attachQueryFile(page: import('@playwright/test').Page, name: string, fixturePath: string): Promise<void> {
  await page.getByTestId('query-file-input').setInputFiles({
    name,
    mimeType: 'text/plain',
    buffer: readFileSync(join(REPOSITORY, fixturePath)),
  });
  await expect(page.getByTestId('query-file-name')).toContainText(name);
}

test.describe('PSIBLAST and TBLASTN -in_pssm UI', () => {
  test('selects PSI-BLAST under BLASTP Algorithm radios, configures MSA query + iterations + PSSM outputs, runs, and offers PSSM downloads', async ({ page }) => {
    await openSearch(page);

    // Top-level mode tabs stay at NCBI's 5 tabs (no 6th PSIBLAST tab).
    await expect(page.locator('[data-testid="program-tabs"] label')).toHaveCount(5);
    await expect(page.getByTestId('program-psiblast')).toHaveCount(0);

    await program(page, 'psiblast');
    await expect(page.getByTestId('program-blastp')).toBeChecked();
    await expect(page.getByTestId('param-task-psiblast')).toBeChecked();

    // Query source radiogroup offers FASTA, MSA (-in_msa), and PSSM (-in_pssm)
    await expect(page.getByTestId('query-source')).toBeVisible();
    await page.getByTestId('query-source-msa').check();
    await expect(page.getByTestId('query-file-picker')).toBeVisible();
    await expect(page.getByTestId('query-file-block')).toContainText('Multiple sequence alignment');
    await expect(page.getByTestId('query-file-block')).not.toContainText('Multiple sequence alignment (FASTA)');

    await attachQueryFile(page, 'seed.aln', MSA_FIXTURE);
    await page.getByTestId('query-msa-master-idx').fill('1');

    await paste(page, 'subject', PSIBLAST_SUBJECT_FASTA);

    await openParameters(page);
    await page.getByTestId('param-num_iterations').fill('2');
    await page.getByTestId('param-out_pssm').fill('out.asn');
    await page.getByTestId('param-out_ascii_pssm').fill('out.ascii');
    await page.getByTestId('param-save_pssm_after_last_round').check();

    await expect(page.getByTestId('search-summary-line')).toContainText('PSIBLAST (MSA seed.aln)');

    await submit(page);
    await waitStatus(page, 1, 'completed');
    await pickRun(page, 1);

    await expect(page.getByTestId('results-query-file')).toContainText('-in_msa seed.aln');

    // Run details shows Query file instead of FASTA record count
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-query-file')).toContainText('-in_msa seed.aln');

    // Outputs tab shows generated PSSM files sub-tab with download buttons
    await page.getByTestId('results-view-outputs').click();
    await page.getByRole('tab', { name: /-out_pssm out\.asn, -out_ascii_pssm out\.ascii/ }).click();
    await expect(page.getByTestId('run-files')).toBeVisible();
    await expect(page.getByTestId('run-file-download-0')).toContainText('Download out.asn');
    await expect(page.getByTestId('run-file-download-1')).toContainText('Download out.ascii');
  });

  test('shows multi-round PSI-BLAST results one round at a time, defaulting to the last round and switching to Round 1', async ({ page }) => {
    await openSearch(page);
    await program(page, 'psiblast');
    await paste(page, 'query', readFileSync(join(REPOSITORY, PSIBLAST_QUERY_FIXTURE), 'utf8'));
    await paste(page, 'subject', readFileSync(join(REPOSITORY, PSIBLAST_SUBJ_IDS_FIXTURE), 'utf8'));
    await openParameters(page);
    await page.getByTestId('param-num_iterations').fill('2');

    await submit(page);
    await waitStatus(page, 1, 'completed');
    await pickRun(page, 1);

    const roundSelect = page.getByTestId('results-round');
    await expect(roundSelect).toBeVisible();
    await expect(roundSelect).toHaveValue('2');

    if (BUILD_HAS_ENGINE) {
      // Round 2: 8 subjects x 1 HSP each; NP_000549.1 (sIdx 0) has E-value 8.24e-18 (rank #2 behind gnl|db|tag7 5.59e-18).
      await expect(page.getByTestId('subject-count')).toHaveText('8 shown');
      await expect(page.getByTestId('subject-row-0')).toContainText('8.24e-18');
      await expect(page.getByTestId('subject-row-0').locator('[data-field="hsps"]')).toHaveText('1');
      await expect(page.getByTestId('subject-row-3')).toContainText('5.59e-18');
      await expect(page.getByTestId('subject-row-3').locator('[data-field="hsps"]')).toHaveText('1');

      // Switch to Round 1: 8 subjects x 1 HSP each; NP_000549.1 (sIdx 0) is #1 with E-value 1.35e-19.
      await roundSelect.selectOption('1');
      await expect(roundSelect).toHaveValue('1');
      await expect(page.getByTestId('subject-count')).toHaveText('8 shown');
      await expect(page.getByTestId('subject-row-0')).toContainText('1.35e-19');
      await expect(page.getByTestId('subject-row-0').locator('[data-field="hsps"]')).toHaveText('1');
      await expect(page.getByTestId('subject-row-0')).toHaveAttribute('data-order', '1');
    } else {
      await roundSelect.selectOption('1');
      await expect(roundSelect).toHaveValue('1');
    }
  });

  test('supports TBLASTN -in_pssm query upload, executes search, and keeps layout overflow-free from 320px to 1280px', async ({ page }) => {
    await openSearch(page);
    await program(page, 'tblastn');
    await expect(page.getByTestId('query-source-pssm')).toBeVisible();
    await expect(page.getByTestId('query-source-msa')).toHaveCount(0);

    await page.getByTestId('query-source-pssm').check();
    await attachQueryFile(page, 'profile.asn', PSSM_FIXTURE);

    for (const width of [320, 768, 900, 1280]) {
      await page.setViewportSize({ width, height: 800 });
      const overflow = await page.evaluate(() => document.documentElement.scrollWidth - document.documentElement.clientWidth);
      expect(overflow, `viewport width ${width}px`).toBeLessThanOrEqual(0);
    }

    await page.setViewportSize({ width: 1280, height: 800 });
    expect(await axeViolations(page, ['main'])).toEqual([]);

    await paste(page, 'subject', readFileSync(join(REPOSITORY, TBLASTN_SUBJECT_FIXTURE), 'utf8'));
    await submit(page);
    await waitStatus(page, 1, 'completed');
    await pickRun(page, 1);
    await expect(page.getByTestId('results-query-file')).toContainText('-in_pssm profile.asn');
    await expect(page.getByTestId('results-hits')).toBeVisible();
  });
});
