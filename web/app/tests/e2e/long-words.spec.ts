// Long words in the page's messages: a name the user gave (an MSA, PSSM or database file) in the line under
// Run LOSAT, and the engine's usage text for an unknown option, whose line of 72 `=` is one word, must not
// make the page wider than the window (line WQA's sweep 4, S4-4 and S4-7, 2026-11-28).
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { REPOSITORY } from './support/harness-server';
import { openParameters, openSearch, paste, program } from './support/search';
import { test } from './support/test';

const WIDTHS = [320, 390, 768, 1280] as const;
const LONG = 'Y'.repeat(400);
const MSA = readFileSync(join(REPOSITORY, 'LOSAT/tests/fixtures/msa/ps6_b_db_0.aln'));
const PSSM = readFileSync(join(REPOSITORY, 'LOSAT/tests/fixtures/pssm/q1b_gi.asn'));
const SUBJECT = '>s1\nTMAGRVELVTGPMFAGKSTYLKNIYQQENGGNKHCLFVKHSLETRYGCGTGTIVTHAGEV\n';

/** Shows that at each width the page has no sideways scroll, and that `selector` (one element) lies inside the window. */
async function expectFits(page: Page, selector: string, context: string): Promise<void> {
  for (const width of WIDTHS) {
    await page.setViewportSize({ width, height: 900 });
    const measured = await page.evaluate((css) => {
      const box = document.querySelector(css)!.getBoundingClientRect();
      return { page: document.documentElement.scrollWidth, window: document.documentElement.clientWidth, right: box.right };
    }, selector);
    expect(measured.page, `${context} at ${width} px: the page scrolls sideways`).toBeLessThanOrEqual(measured.window);
    expect(measured.right, `${context} at ${width} px: the text passes the window's edge`).toBeLessThanOrEqual(measured.window);
  }
}

const SUMMARY = '[data-testid="search-summary-line"]';

test.describe('long words stay inside the page', () => {
  for (const name of [`${LONG}.aln`, 'hemoglobin_alpha_family_seed_alignment.aln']) {
    test(`an MSA named ${name.length} characters in one word wraps in the line under Run LOSAT`, async ({ page }) => {
      await openSearch(page);
      await program(page, 'psiblast');
      await page.getByTestId('query-source-msa').check();
      await page.getByTestId('query-file-input').setInputFiles({ name, mimeType: 'text/plain', buffer: MSA });
      await expect(page.getByTestId('query-file-name')).toContainText(name.slice(0, 20));
      await paste(page, 'subject', SUBJECT);
      await expect(page.getByTestId('search-summary-line')).toContainText(`PSIBLAST (MSA ${name.slice(0, 20)}`);
      await expectFits(page, SUMMARY, 'MSA name');
    });
  }

  test('a PSSM file with a name of 400 characters in one word wraps in the line under Run LOSAT', async ({ page }) => {
    await openSearch(page);
    await program(page, 'psiblast');
    await page.getByTestId('query-source-pssm').check();
    await page.getByTestId('query-file-input').setInputFiles({ name: `${LONG}.asn`, mimeType: 'application/octet-stream', buffer: PSSM });
    await paste(page, 'subject', SUBJECT);
    await expect(page.getByTestId('search-summary-line')).toContainText('PSIBLAST (PSSM YYYY');
    await expectFits(page, SUMMARY, 'PSSM name');
  });

  test('a database with a name of 400 characters in one word wraps in the line under Run LOSAT', async ({ page }) => {
    await openSearch(page);
    await paste(page, 'query', '>q\nACGTACGTACGTACGTACGTAAAAACCCCCGGGGGTTTTTACGTACGTAC\n');
    await page.getByTestId('subject-source-database').check();
    await page.getByTestId('subject-db-files-input').setInputFiles([
      { name: `${LONG}.nin`, mimeType: 'application/octet-stream', buffer: Buffer.from('0') },
    ]);
    await expect(page.getByTestId('search-summary-line')).toContainText('database YYYY');
    await expectFits(page, SUMMARY, 'database name');
  });

  for (const mode of ['blastn', 'tblastx'] as const) {
    test(`the usage text for an unknown option in ${mode} wraps its line of 72 = inside the card`, async ({ page }) => {
      await openSearch(page);
      await program(page, mode);
      await paste(page, 'query', '>q\nACGTACGTACGTACGTACGTAAAAACCCCCGGGGGTTTTTACGTACGTAC\n');
      await paste(page, 'subject', '>s\nACGTACGTACGTACGTACGTAAAAACCCCCGGGGGTTTTTACGTACGTAC\n');
      await openParameters(page);
      await page.getByTestId('param-extra').fill('-not_an_option');
      const validation = page.getByTestId('argv-validation');
      await expect(validation).toHaveAttribute('data-state', 'invalid');
      await expect(validation).toContainText('='.repeat(72));
      await expectFits(page, '[data-testid="argv-validation"] .error', 'usage text');
    });
  }
});
