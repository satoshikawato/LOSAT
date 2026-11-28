// Long words in the page's messages: a name the user gave (an MSA, PSSM or database file) in the line under
// Run LOSAT, and the engine's usage text for an unknown option, whose line of 72 `=` is one word, must not
// make the page wider than the window (line WQA's sweep 4, S4-4 and S4-7, 2026-11-28). The same word in the other
// boxes that show a name or an ID (the review's M1): the Alignments summary, the HSP list's heading, the Dot plot's
// caption, Format Archive's chosen file, WindowMasker's output tab and download button, and the words of a failed run.
// The rule must not break the page's own words (reviews r2 and r3: "TBLASTN" and the Dot Plot's "Subject" over two
// lines at 270 to 370 px), and a help stays inside the page when the window narrows while it is open (review r2, L2).
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { openParameters, openSearch, paste, pickRun, program, runMessage, submit, waitStatus } from './support/search';
import { test } from './support/test';

const WIDTHS = [320, 390, 768, 1280] as const;
const LONG = 'Y'.repeat(400);
const MSA = readFileSync(join(REPOSITORY, 'LOSAT/tests/fixtures/msa/ps6_b_db_0.aln'));
const PSSM = readFileSync(join(REPOSITORY, 'LOSAT/tests/fixtures/pssm/q1b_gi.asn'));
const SUBJECT = '>s1\nTMAGRVELVTGPMFAGKSTYLKNIYQQENGGNKHCLFVKHSLETRYGCGTGTIVTHAGEV\n';

/** Shows that at each width the page has no sideways scroll, and that every element of `selector` lies inside the window. */
async function expectFits(page: Page, selector: string, context: string): Promise<void> {
  for (const width of WIDTHS) {
    await page.setViewportSize({ width, height: 900 });
    const measured = await page.evaluate((css) => {
      const boxes = [...document.querySelectorAll(css)].map((element) => element.getBoundingClientRect());
      return { page: document.documentElement.scrollWidth, window: document.documentElement.clientWidth, right: Math.max(...boxes.map((box) => box.right)) };
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

  // The Results page: an ID of an ordinary length (39 characters: a RefSeq accession and its species) and a very long one.
  for (const length of [39, 400]) {
    test(`the Results page wraps a sequence ID of ${length} characters (Alignments, HSP list, Dot plot caption)`, async ({ page }) => {
      test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
      test.setTimeout(180_000);
      const id = 'NC_000913.3_Escherichia_coli_K12_MG1655_'.padEnd(length, 'x').slice(0, length);
      const unit = 'ACGTTGCAAGCTTGCAGGCATCGATCGGATCGATTTAGCGCGATATCGCGATTAGCGGCTAATCGGATCGATCGTAGGCTAGCTAGGATCGGCTAGCTAG';
      const fasta = `>${id} description\n${unit}${unit}\n`;
      await openSearch(page);
      await paste(page, 'query', fasta);
      await paste(page, 'subject', fasta);
      await submit(page);
      await waitStatus(page, 1, 'completed');
      await pickRun(page, 1);
      await expect(page.getByTestId('alignments-summary')).toContainText(id.slice(0, 30));
      await expectFits(page, '[data-testid="alignments-summary"]', `Alignments summary, ID of ${length}`);
      await expectFits(page, '[data-testid="hsp-table"] h3', `HSP list heading, ID of ${length}`);
      await page.getByTestId('pane-dotplot').click();
      const caption = '.plot-caption:has([data-testid="dotplot-title"])';
      await expect(page.locator(caption)).toContainText(id.slice(0, 30));
      await expectFits(page, caption, `Dot plot caption, ID of ${length}`);
    });
  }

  for (const name of ['rid_search_strategy_archive_2026_11_28_with_a_long_name.asn', `${'a'.repeat(400)}.asn`]) {
    test(`Format Archive wraps the chosen file's name of ${name.length} characters`, async ({ page }) => {
      await openSearch(page);
      await page.getByTestId('tab-blast-formatter').click();
      await page.getByTestId('bf-archive-file-input').setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from('x') });
      await expect(page.getByTestId('bf-archive-chosen')).toContainText(name.slice(0, 30));
      await expectFits(page, '[data-testid="bf-archive-chosen"]', 'Format Archive file');
    });
  }

  test("WindowMasker's output tab and download button wrap an -out name of 200 characters", async ({ page }) => {
    test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
    test.setTimeout(180_000);
    const unit = 'ACGTTGCAAGCTTGCAGGCATCGATCGGATCGATTTAGCGCGATATCGCGATTAGCGGCTAATCGGATCGATCGTAGGCTAGCTAGGATCGGCTAGCTAG';
    await openSearch(page);
    await page.getByTestId('tab-windowmasker').click();
    await expect(page.getByTestId('wm-panel')).toHaveAttribute('data-describe', 'ready', { timeout: 60_000 });
    await page.getByTestId('wm-fasta-input').fill(`>w1\n${unit.repeat(40)}\n`);
    await page.getByTestId('wm-out').fill(`${'Z'.repeat(200)}.counts`);
    await expect(page.getByTestId('wm-argv-validation')).not.toHaveAttribute('data-state', 'checking');
    await page.getByTestId('wm-run-losat').click();
    await waitStatus(page, 1, 'completed', 120_000);
    await pickRun(page, 1);
    await expect(page.getByTestId('tool-tab-file-0')).toContainText('ZZZZ');
    await page.getByTestId('tool-tab-file-0').click();
    await expectFits(page, '[data-testid="tool-tab-file-0"]', 'WindowMasker output tab');
    await expectFits(page, '[data-testid="tool-output-download"]', 'WindowMasker download button');
  });

  // The words of a failed run (the results box) carry the name that the user gave: -out_pssm in a folder of 400 characters.
  test('the words of a failed run wrap a name of 400 characters', async ({ page }) => {
    test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
    test.setTimeout(180_000);
    const fasta = (path: string): string => readFileSync(join(REPOSITORY, 'LOSAT/tests/fasta/psiblast', path), 'utf8');
    await openSearch(page);
    await program(page, 'psiblast');
    await paste(page, 'query', fasta('q1b_gi.faa'));
    await paste(page, 'subject', fasta('q1_subj_ids.faa'));
    await openParameters(page);
    await page.getByTestId('param-num_iterations').fill('2');
    await page.getByTestId('param-out_pssm').fill(`${LONG}/x.asn`);
    await submit(page);
    await waitStatus(page, 1, 'failed');
    expect(await runMessage(page, 1)).toContain('YYYY');
    await expectFits(page, '[data-testid="results-status"]', 'failed run message');
  });

  // A source's card says which input stopped the search; with a long name in it the words wrap.
  test("an input's card names the file that stopped the search, and wraps its name of 400 characters", async ({ page }) => {
    test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
    await openSearch(page);
    await page.getByTestId('query-files').setInputFiles([
      { name: `${LONG}.fa`, mimeType: 'text/plain', buffer: Buffer.from('>a\nAC#GT\n') },
      { name: 'second.fa', mimeType: 'text/plain', buffer: Buffer.from('>q\nACGTACGTACGTACGTACGTAAAAACCCCCGGGGGTTTTTACGTACGTAC\n') },
    ]);
    await expect(page.locator('.warning').filter({ hasText: 'Not searched' }).first()).toContainText('YYYY', { timeout: 30_000 });
    await expectFits(page, '.warning', 'source card warning');
  });
});

// The page's own words: a tab's name has a line as wide as its tab, so the rule for long words never breaks it.
test.describe("the page's own words and its helps", () => {
  const NARROW_TO_WIDE = [280, 290, 300, 320, 340, 360, 370, 390, 800, 1280] as const;

  /**
   * The page's words that take more than one line inside a word: a run of letters or digits whose letters are on two
   * lines (a break at a hyphen, a slash or a space is not inside a word). Lists, text boxes and hidden text are left out.
   */
  async function wordsBrokenInside(page: Page): Promise<string[]> {
    return page.evaluate(() => {
      const broken: string[] = [];
      const walker = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
      const range = document.createRange();
      for (let node = walker.nextNode(); node !== null; node = walker.nextNode()) {
        const parent = node.parentElement;
        if (parent === null || parent.closest('select, textarea, script, style, .visually-hidden') !== null) continue;
        for (const word of (node.textContent ?? '').matchAll(/[\p{L}\p{N}]+/gu)) {
          range.setStart(node, word.index);
          range.setEnd(node, word.index + word[0].length);
          const lines = new Set([...range.getClientRects()].filter((rect) => rect.width > 0).map((rect) => Math.round(rect.top)));
          if (lines.size > 1) broken.push(`${word[0]} (${parent.className || parent.tagName.toLowerCase()})`);
        }
      }
      return broken;
    });
  }

  /** Shows that no word of the page now breaks inside at any of the widths. */
  async function expectWholeWords(page: Page, screen: string): Promise<void> {
    for (const width of NARROW_TO_WIDE) {
      await page.setViewportSize({ width, height: 900 });
      expect(await wordsBrokenInside(page), `${screen} at ${width} px: words broken over two lines`).toEqual([]);
    }
  }

  // The page's own words in every screen that the cases above open, with ordinary names and IDs (D32's reviews r2 and
  // r3: the mode tabs' "TBLASTN" and the Dot Plot's "Subject" broke inside at 270 to 370 px), with Inter and with the
  // font that the browser uses when Inter does not come (blocked web fonts; a first load before the font came).
  for (const font of ['Inter', 'the fallback font'] as const) {
    test(`the page's own words never break inside a word, from 280 to 1280 px, in ${font}`, async ({ page }) => {
      test.setTimeout(180_000);
      if (font === 'the fallback font') await page.route('**/fonts/inter/**', (route) => route.abort());
      await openSearch(page);
      const inter = await page.evaluate(async () => {
        await document.fonts.ready;
        return [...document.fonts].some((face) => face.family.replace(/["']/g, '') === 'Inter' && face.status === 'loaded');
      });
      expect(inter, 'Inter is loaded').toBe(font === 'Inter');
      const fasta = (id: string) => `>${id}\nACGTTGCAAGCTTGCAGGCATCGATCGGATCGATTTAGCGCGATATCGCGATTAGCGGCTAATCGGATCGATCGTAGGCTAGCTAG\n`;
      await paste(page, 'query', fasta('q1'));
      await paste(page, 'subject', fasta('s1'));
      await openParameters(page);
      await expectWholeWords(page, 'BLASTN form with its parameters');
      for (const mode of ['tblastx', 'psiblast'] as const) {
        await program(page, mode);
        await openParameters(page);
        await expectWholeWords(page, `${mode} form with its parameters`);
      }
      await program(page, 'blastn');
      await page.getByTestId('subject-source-database').check();
      await expectWholeWords(page, 'BLASTN form with a database as subjects');
      await page.getByTestId('subject-source-fasta').check();
      await submit(page);
      await waitStatus(page, 1, 'completed');
      await pickRun(page, 1);
      await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
      await expectWholeWords(page, 'Results');
      await page.getByTestId('pane-dotplot').click();
      await expect(page.getByTestId('dotplot-subject')).toBeVisible();
      await expectWholeWords(page, 'Results, Dot Plot');
      await page.getByTestId('tab-windowmasker').click();
      for (const stage of ['mk_counts', 'mask', 'convert'] as const) {
        await page.getByTestId(`wm-stage-${stage}`).check({ force: true });
        await expectWholeWords(page, `WindowMasker, ${stage}`);
      }
      await page.getByTestId('tab-blast-formatter').click();
      await expectWholeWords(page, 'Format Archive, file');
      await page.getByTestId('bf-source-paste').check();
      await expectWholeWords(page, 'Format Archive, paste');
    });
  }

  test('every help opened in a wide window stays inside the page when the window narrows', async ({ page }) => {
    // 24 helps opened at 1280 px and seen at three widths, then opened at two: WebKit took 33 s.
    test.setTimeout(120_000);
    await page.setViewportSize({ width: 1280, height: 900 });
    await openSearch(page);
    await openParameters(page);
    const toggles = await page
      .locator('.help-toggle:visible')
      .evaluateAll((buttons) => buttons.map((button) => button.getAttribute('data-testid') ?? ''));
    expect(toggles.length).toBeGreaterThan(20);
    /** "inside" when the open help of `id` and the page are inside the window, else where they are. */
    const where = (id: string) =>
      page.evaluate((help) => {
        const box = document.querySelector(`[data-testid="${help}"]`)!.getBoundingClientRect();
        const window = document.documentElement.clientWidth;
        const scroll = document.documentElement.scrollWidth;
        return box.left >= 0 && box.right <= window && scroll <= window
          ? 'inside'
          : `help ${box.left.toFixed(1)}..${box.right.toFixed(1)}, page ${scroll} in a window of ${window}`;
      }, id.replace(/-toggle$/, ''));
    for (const id of toggles) {
      await page.setViewportSize({ width: 1280, height: 900 });
      await page.getByTestId(id).click();
      await expect(page.getByTestId(id.replace(/-toggle$/, ''))).toBeVisible();
      for (const width of [390, 320, 280]) {
        await page.setViewportSize({ width, height: 900 });
        await expect
          .poll(() => where(id), { message: `${id} at ${width} px after opening at 1280 px`, timeout: 3_000 })
          .toBe('inside');
      }
      await page.getByTestId(id).click();
      await expect(page.getByTestId(id.replace(/-toggle$/, ''))).toHaveCount(0);
    }
    // Opened at the width it is seen at.
    for (const width of [280, 320]) {
      await page.setViewportSize({ width, height: 900 });
      for (const id of toggles) {
        await page.getByTestId(id).click();
        await expect.poll(() => where(id), { message: `${id} opened at ${width} px`, timeout: 3_000 }).toBe('inside');
        await page.getByTestId(id).click();
      }
    }
  });
});
