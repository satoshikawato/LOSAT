// Long words in the page's messages: a name the user gave (an MSA, PSSM or database file) in the line under
// Run LOSAT, and the engine's usage text for an unknown option, whose line of 72 `=` is one word, must not
// make the page wider than the window (line WQA's sweep 4, S4-4 and S4-7, 2026-11-28). The same word in the other
// boxes that show a name or an ID (the review's M1): the Alignments summary, the HSP list's heading, the Dot plot's
// caption, Format Archive's chosen file, WindowMasker's output tab and download button, and the words of a failed run.
// The rule must not break the page's own words (reviews r2 and r3: "TBLASTN" and the Dot Plot's "Subject" over two
// lines at 270 to 370 px), and a help stays inside the page when the window narrows while it is open (review r2, L2).
// The guard of the page's own words visits every tab of the page's tab bar and every stage of each (read from the page,
// so a tab added later fails it until it is visited), and the tool tabs' names of files (Masking, makeblastdb,
// WindowMasker: inputs, outputs, download buttons, Run details) wrap without breaking the page's words (D32 on dev).
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import { expect, type Page } from '@playwright/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { REPOSITORY } from './support/harness-server';
import { described, openFiles, openParameters, openSearch, paste, pickRun, program, runMessage, submit, waitStatus } from './support/search';
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

/** Widths from a small phone to a desktop, with the steps where the page's words came close to their boxes' edges. */
const NARROW_TO_WIDE = [280, 290, 300, 320, 340, 360, 370, 390, 800, 1280] as const;
/** A DNA sequence of 100 letters, for runs that need one. */
const DNA = 'ACGTTGCAAGCTTGCAGGCATCGATCGGATCGATTTAGCGCGATATCGCGATTAGCGGCTAATCGGATCGATCGTAGGCTAGCTAGGATCGGCTAGCTAG';

/**
 * The page's words that take more than one line inside a word: a run of letters or digits whose letters are on two
 * lines (a break at a hyphen, a slash or a space is not inside a word). Lists, text boxes and hidden text are left out,
 * and so are the words of `names` (names that the user gave: they wrap where they must) and the words of a string
 * without spaces that no line can hold whole (a SHA-256 of 64 digits, a run's file name, a long ID): one wider than
 * the page's column, or wider than its line where the line's width does not come from its text (a block, not a flex
 * or grid item, an inline block, a cell, a float or a placed box, which a word shrinks to and then breaks in).
 */
async function wordsBrokenInside(page: Page, names: readonly string[] = []): Promise<string[]> {
  return page.evaluate((given) => {
    const broken: string[] = [];
    const walker = document.createTreeWalker(document.body, NodeFilter.SHOW_TEXT);
    const range = document.createRange();
    const measure = document.createElement('canvas').getContext('2d')!;
    // The widest line of the page: the window less its 16 px side gutters.
    const column = document.documentElement.clientWidth - 32;
    /** The box whose lines hold an element's text: the element, or its nearest ancestor that is not inline. */
    const lineBox = (element: Element): Element => {
      let box = element;
      while (box.parentElement !== null && getComputedStyle(box).display === 'inline') box = box.parentElement;
      return box;
    };
    /** Whether a box's width comes from its content (so that a word could make it wider instead of breaking). */
    const sizedByContent = (box: Element): boolean => {
      const style = getComputedStyle(box);
      const outer = box.parentElement === null ? 'block' : getComputedStyle(box.parentElement).display;
      return (
        /^(inline-|table-cell|table-caption|flex|grid)/.test(style.display) ||
        style.float !== 'none' ||
        /^(absolute|fixed)$/.test(style.position) ||
        /(^|-)(flex|grid)$/.test(outer)
      );
    };
    for (let node = walker.nextNode(); node !== null; node = walker.nextNode()) {
      const parent = node.parentElement;
      if (parent === null || parent.closest('select, textarea, script, style, .visually-hidden') !== null) continue;
      const text = node.textContent ?? '';
      for (const word of text.matchAll(/[\p{L}\p{N}]+/gu)) {
        if (given.some((name) => name.includes(word[0]))) continue;
        range.setStart(node, word.index);
        range.setEnd(node, word.index + word[0].length);
        const lines = new Set([...range.getClientRects()].filter((rect) => rect.width > 0).map((rect) => Math.round(rect.top)));
        if (lines.size < 2) continue;
        const start = text.slice(0, word.index).search(/\S*$/);
        const end = word.index + word[0].length + (text.slice(word.index + word[0].length).match(/^\S*/)?.[0].length ?? 0);
        const style = getComputedStyle(parent);
        measure.font = `${style.fontStyle} ${style.fontWeight} ${style.fontSize} ${style.fontFamily}`;
        const width = measure.measureText(text.slice(start, end)).width;
        const box = lineBox(parent);
        const boxStyle = getComputedStyle(box);
        const line = box.clientWidth - parseFloat(boxStyle.paddingLeft) - parseFloat(boxStyle.paddingRight);
        if (width > column || (width > line && !sizedByContent(box))) continue;
        broken.push(`${word[0]} (${parent.className || parent.tagName.toLowerCase()})`);
      }
    }
    return broken;
  }, names);
}

/** Shows that no word of the page breaks inside at any of the widths. */
async function expectWholeWords(page: Page, screen: string): Promise<void> {
  for (const width of NARROW_TO_WIDE) {
    await page.setViewportSize({ width, height: 900 });
    expect(await wordsBrokenInside(page), `${screen} at ${width} px: words broken over two lines`).toEqual([]);
  }
}

/** Shows that at 320 and 390 px the page has no sideways scroll and keeps its own words whole while the user's `names` wrap. */
async function expectNamesWrap(page: Page, names: readonly string[], screen: string): Promise<void> {
  for (const width of [320, 390]) {
    await page.setViewportSize({ width, height: 900 });
    const measured = await page.evaluate(() => ({ page: document.documentElement.scrollWidth, window: document.documentElement.clientWidth }));
    expect(measured.page, `${screen} at ${width} px: the page scrolls sideways`).toBeLessThanOrEqual(measured.window);
    expect(await wordsBrokenInside(page, names), `${screen} at ${width} px: the page's words broken over two lines`).toEqual([]);
  }
}

/** Opens every closed options section ("… options", "Algorithm parameters") of the tab shown. */
async function openSections(page: Page): Promise<void> {
  const closed = page.locator('.parameters-bar[aria-expanded="false"]:visible');
  for (let left = 10; left > 0 && (await closed.count()) > 0; left--) await closed.first().click();
  await expect(closed).toHaveCount(0);
}

/** Opens a tool tab and waits for its form (the engine's `-xmlhelp`). */
async function openTool(page: Page, tab: 'windowmasker' | 'masking' | 'makeblastdb'): Promise<void> {
  await page.getByTestId(`tab-${tab}`).click();
  await described(page, tab === 'windowmasker' ? 'wm-panel' : tab === 'masking' ? 'ms-panel' : 'mk-panel', 'ready');
}

/** Presses a tool tab's Run LOSAT once its command is checked, and waits until run `run` completes on the Results tab. */
async function runTool(page: Page, prefix: 'wm' | 'ms' | 'mk', run: number): Promise<void> {
  await expect(page.getByTestId(`${prefix}-argv-validation`)).toHaveAttribute('data-state', 'ok', { timeout: 30_000 });
  await page.getByTestId(`${prefix}-run-losat`).click();
  await waitStatus(page, run, 'completed', 120_000);
  await pickRun(page, run);
  await expect(page.getByTestId('tool-run-view')).toHaveAttribute('data-run', String(run));
}

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

  // The screens that show the names of a tool tab's files (D32 on dev: the Masking and makeblastdb tabs came after the
  // page rule): names of 200 characters in one word in the inputs' lists, the outputs' tabs and download buttons, the
  // made files, the Run details, and the search that uses the made database.
  const LONG_FA = `${'Z'.repeat(200)}.fa`;
  const FASTA = `>n1 long names\n${DNA}${DNA}\n`;

  test('Masking wraps the names of its input files, its run files and its -out of 200 characters', async ({ page }) => {
    test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
    test.setTimeout(240_000);
    const second = `${'Y'.repeat(200)}.fa`;
    const conf = `${'V'.repeat(200)}.conf`;
    const out = `${'W'.repeat(200)}.txt`;
    const names = [LONG_FA, second, conf, out];
    await openSearch(page);
    await openTool(page, 'masking');
    await page.getByTestId('ms-fasta-files').setInputFiles([
      { name: LONG_FA, mimeType: 'text/plain', buffer: Buffer.from(FASTA) },
      { name: second, mimeType: 'text/plain', buffer: Buffer.from(FASTA.replace('n1', 'n2')) },
    ]);
    await expect(page.getByTestId('ms-fasta-name-1')).toContainText('YYYY');
    await page.getByTestId('ms-out').fill(out);
    await openSections(page);
    await page.getByTestId('ms-run-files-input').setInputFiles({ name: conf, mimeType: 'text/plain', buffer: Buffer.from('[BLAST]\n') });
    await expect(page.getByTestId('ms-run-file-name-0')).toContainText('VVVV');
    await expectNamesWrap(page, names, 'Masking form');
    await runTool(page, 'ms', 1);
    await expect(page.getByTestId('tool-tab-file-0')).toContainText('WWWW');
    await page.getByTestId('tool-tab-file-0').click();
    await expectNamesWrap(page, names, 'Masking run, Outputs');
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('tool-details-files')).toContainText('VVVV');
    await expectNamesWrap(page, names, 'Masking run, Run details');
  });

  test('makeblastdb wraps the names of its input file and of the database files of 200 characters, also in the search that uses them', async ({
    page,
  }) => {
    test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
    test.setTimeout(240_000);
    const out = 'X'.repeat(200);
    const names = [LONG_FA, out];
    await openSearch(page);
    await openTool(page, 'makeblastdb');
    await page.getByTestId('mk-input-files').setInputFiles({ name: LONG_FA, mimeType: 'text/plain', buffer: Buffer.from(FASTA) });
    await expect(page.getByTestId('mk-input-name-0')).toContainText('ZZZZ');
    await page.getByTestId('mk-param-dbtype-nucl').check();
    await page.getByTestId('mk-param-out').fill(out);
    await openSections(page);
    await expectNamesWrap(page, names, 'makeblastdb form');
    await runTool(page, 'mk', 1);
    await expect(page.getByTestId('tool-made-file-name-0')).toContainText('XXXX');
    await expectNamesWrap(page, names, 'makeblastdb run, Outputs (made files)');
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('tool-details-files')).toContainText('ZZZZ');
    await expectNamesWrap(page, names, 'makeblastdb run, Run details');
    await page.getByTestId('results-view-outputs').click();
    await page.getByTestId('tool-use-as-database').click();
    await expect(page.getByTestId('subject-db-selected')).toContainText('XXXX');
    await expectNamesWrap(page, names, 'Search with the made database');
  });

  test('WindowMasker wraps the names of its FASTA file and of its counts file of 200 characters', async ({ page }) => {
    const counts = `${'U'.repeat(200)}.obinary`;
    const names = [LONG_FA, counts];
    await openSearch(page);
    await openTool(page, 'windowmasker');
    await page.getByTestId('wm-fasta-files').setInputFiles({ name: LONG_FA, mimeType: 'text/plain', buffer: Buffer.from(FASTA) });
    await expect(page.getByTestId('wm-fasta-name-0')).toContainText('ZZZZ');
    await openSections(page);
    await expectNamesWrap(page, names, 'WindowMasker, counts');
    await page.getByTestId('wm-stage-mask').check({ force: true });
    await page.getByTestId('wm-ustat-file').setInputFiles({ name: counts, mimeType: 'application/octet-stream', buffer: Buffer.from('counts') });
    await expect(page.getByTestId('wm-ustat-chosen')).toContainText('UUUU');
    await expectNamesWrap(page, names, 'WindowMasker, mask');
  });

  test("a search's Run details, Outputs and Download All menu wrap a query file's name of 200 characters", async ({ page }) => {
    test.skip(!BUILD_HAS_ENGINE, 'needs the engine build (LOSAT_WEB_REACTORS)');
    test.setTimeout(180_000);
    await openSearch(page);
    await openFiles(page, 'query', [{ name: LONG_FA, text: FASTA }]);
    await paste(page, 'subject', FASTA);
    await expectNamesWrap(page, [LONG_FA], 'Search form');
    await submit(page);
    await waitStatus(page, 1, 'completed');
    await pickRun(page, 1);
    await expect(page.getByTestId('results-hits')).toHaveAttribute('data-run', '1');
    await page.getByTestId('results-view-details').click();
    await expect(page.getByTestId('run-details')).toContainText('ZZZZ');
    await expectNamesWrap(page, [LONG_FA], 'Run details');
    await page.getByTestId('results-view-outputs').click();
    await expectNamesWrap(page, [LONG_FA], 'Outputs');
    await page.getByTestId('results-download-all').click();
    await expect(page.getByTestId('results-download-menu')).toBeVisible();
    await expectNamesWrap(page, [LONG_FA], 'Download All menu');
  });
});

// The page's own words: a tab's name has a line as wide as its tab, so the rule for long words never breaks it.
test.describe("the page's own words and its helps", () => {
  /** Waits until the view shown has read what it shows. */
  async function viewRead(page: Page): Promise<void> {
    await expect(page.getByRole('tabpanel').getByText(/^Reading\b/)).toHaveCount(0, { timeout: 30_000 });
  }

  /**
   * How the guard visits each tab of the page's tab bar: every stage of the tab (read from the page: the Search modes and
   * input sources, WindowMasker's stages, the Masking tools, makeblastdb's input types, Format Archive's sources; the
   * Results views of each run), with its options sections open.
   */
  const VISITS: Readonly<Record<string, (page: Page) => Promise<void>>> = {
    'tab-search': async (page) => {
      const modes = await page.getByTestId('program-tabs').locator('input[type="radio"]').evaluateAll((inputs) =>
        inputs.map((input) => (input as HTMLInputElement).value),
      );
      expect(modes.length).toBeGreaterThan(4);
      // PSI-BLAST is a choice inside BLASTP's Task Selection (as on NCBI BLAST), with sources of its own.
      for (const mode of [...modes, 'psiblast']) {
        await program(page, mode);
        await openParameters(page);
        await expectWholeWords(page, `Search, ${mode} with its parameters`);
        for (const role of ['query', 'subject'] as const) {
          const source = (value = '') => page.locator(`input[type="radio"][name="${role}-source"]${value && `[value="${value}"]`}`);
          const values = await source().evaluateAll((inputs) => inputs.map((input) => (input as HTMLInputElement).value));
          for (const value of values.filter((value) => value !== 'fasta')) {
            await source(value).check();
            await expectWholeWords(page, `Search, ${mode}, ${role} from ${value}`);
          }
          if (values.length > 1) await source('fasta').check();
        }
      }
      await program(page, 'blastn');
    },
    'tab-results': async (page) => {
      const runs = await page.getByTestId('results-run').locator('option[data-testid]').evaluateAll((options) =>
        options.map((option) => (option as HTMLOptionElement).value),
      );
      expect(runs.length, 'runs of the session').toBe(3);
      for (const run of runs) {
        await page.getByTestId('results-run').selectOption(run);
        const name = await page.getByTestId('results-run').locator(`option[value="${run}"]`).getAttribute('data-testid');
        const views = page.locator('.results-tabs [role="tab"]');
        await expect(views.first()).toBeVisible();
        for (const view of await views.evaluateAll((tabs) => tabs.map((tab) => tab.getAttribute('data-testid') ?? ''))) {
          await page.getByTestId(view).click();
          await expect(page.getByTestId(view)).toHaveAttribute('aria-selected', 'true');
          await viewRead(page);
          await expectWholeWords(page, `Results of ${name}, ${view}`);
        }
        if (await page.getByTestId('results-download-all').isVisible()) {
          await page.getByTestId('results-download-all').click();
          await expect(page.getByTestId('results-download-menu')).toBeVisible();
          await expectWholeWords(page, `Results of ${name}, Download All menu`);
          await page.getByTestId('results-download-all').click();
        }
      }
    },
    'tab-windowmasker': async (page) => {
      await described(page, 'wm-panel', 'ready');
      for (const stage of await page.getByTestId('wm-stage-tabs').locator('input[type="radio"]').all()) {
        await stage.check({ force: true });
        await openSections(page);
        await expectWholeWords(page, `WindowMasker, ${await stage.getAttribute('value')}`);
      }
    },
    'tab-masking': async (page) => {
      for (const tool of await page.getByTestId('ms-tool-tabs').locator('input[type="radio"]').all()) {
        await tool.check();
        await described(page, 'ms-panel', 'ready');
        await openSections(page);
        await expectWholeWords(page, `Masking, ${await tool.getAttribute('value')}`);
      }
    },
    'tab-makeblastdb': async (page) => {
      await described(page, 'mk-panel', 'ready');
      const field = page.getByTestId('mk-param-input_type');
      const isSelect = (await field.evaluate((element) => element.tagName)) === 'SELECT';
      const types = await field
        .locator(isSelect ? 'option' : 'input[type="radio"]')
        .evaluateAll((choices) => choices.map((choice) => (choice as HTMLInputElement).value).filter((value) => value !== ''));
      expect(types.length, 'input types').toBeGreaterThan(1);
      for (const type of types) {
        if (isSelect) await field.selectOption(type);
        else await field.locator(`input[value="${type}"]`).check();
        await openSections(page);
        await expectWholeWords(page, `makeblastdb, -input_type ${type}`);
      }
      if (isSelect) await field.selectOption('fasta');
      else await field.locator('input[value="fasta"]').check();
    },
    'tab-blast-formatter': async (page) => {
      // A source with nothing to choose from (no archive of this session) cannot be chosen.
      for (const source of await page.locator('input[type="radio"][data-testid^="bf-source-"]:enabled').all()) {
        await source.check();
        await openSections(page);
        await expectWholeWords(page, `Format Archive, ${await source.getAttribute('value')}`);
      }
    },
  };

  // The page's own words in every tab and stage, with ordinary names and IDs (D32's reviews r2 and r3: the mode tabs'
  // "TBLASTN" and the Dot Plot's "Subject" broke inside at 270 to 370 px), with Inter and with the font that the
  // browser uses when Inter does not come (blocked web fonts; a first load before the font came). The tabs are read
  // from the page's tab bar: a tab that the guard does not know fails it.
  for (const font of ['Inter', 'the fallback font'] as const) {
    test(`the page's own words never break inside a word in any tab or stage, from 280 to 1280 px, in ${font}`, async ({ page }) => {
      test.setTimeout(900_000);
      if (font === 'the fallback font') await page.route('**/fonts/inter/**', (route) => route.abort());
      await openSearch(page);
      const inter = await page.evaluate(async () => {
        await document.fonts.ready;
        return [...document.fonts].some((face) => face.family.replace(/["']/g, '') === 'Inter' && face.status === 'loaded');
      });
      expect(inter, 'Inter is loaded').toBe(font === 'Inter');
      // A run of each kind for the Results tab: a search (Hits, Dot Plot), a tool's outputs (Masking), a tool's made
      // files (makeblastdb).
      const fasta = (id: string) => `>${id}\n${DNA}\n`;
      await paste(page, 'query', fasta('q1'));
      await paste(page, 'subject', fasta('s1'));
      await submit(page);
      await waitStatus(page, 1, 'completed');
      await openTool(page, 'masking');
      await page.getByTestId('ms-fasta-input').fill(fasta('m1'));
      await runTool(page, 'ms', 2);
      await openTool(page, 'makeblastdb');
      await page.getByTestId('mk-fasta-input').fill(fasta('d1'));
      await page.getByTestId('mk-param-dbtype-nucl').check();
      await runTool(page, 'mk', 3);
      const tabs = await page
        .getByRole('tablist', { name: 'Main' })
        .getByRole('tab')
        .evaluateAll((elements) => elements.map((element) => element.getAttribute('data-testid') ?? ''));
      expect([...tabs].sort(), "the tab bar's tabs are the tabs that the guard visits").toEqual(Object.keys(VISITS).sort());
      for (const tab of tabs) {
        await page.getByTestId(tab).click();
        await expect(page.getByTestId(tab)).toHaveAttribute('aria-selected', 'true');
        await VISITS[tab]!(page);
      }
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
