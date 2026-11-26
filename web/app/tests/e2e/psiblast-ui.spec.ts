// LOSAT Web PSIBLAST & TBLASTN -in_pssm UI, PSSM download, and responsive viewport E2E tests (PSW1).
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect } from '@playwright/test';
import { axeViolations } from './support/axe';
import { REPOSITORY } from './support/harness-server';
import { openParameters, openSearch, paste, pickRun, program, submit, waitStatus } from './support/search';
import { test } from './support/test';

const MSA_FIXTURE = 'LOSAT/tests/fixtures/msa/ps6_b_db_0.aln';
const PSSM_FIXTURE = 'LOSAT/tests/fixtures/pssm/q1b_gi.asn';
const TBLASTN_SUBJECT_FIXTURE = 'LOSAT/tests/fasta/psiblast/q1b_gi_nt.fna';

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
  test('switches to PSIBLAST, configures MSA query + iterations + PSSM outputs, runs, and offers PSSM downloads', async ({ page }) => {
    await openSearch(page);
    await program(page, 'psiblast');

    // Query source radiogroup offers FASTA, MSA (-in_msa), and PSSM (-in_pssm)
    await expect(page.getByTestId('query-source')).toBeVisible();
    await page.getByTestId('query-source-msa').check();
    await expect(page.getByTestId('query-file-picker')).toBeVisible();

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
