// W8 (line WA, stage 2): the marks of the form and of what shows on hover and focus.
// WA2-AX-05: a field that a refusal of the options names is marked as refused and points to the
// message. WA2-AX-07: the "?" help and the Graphic Summary's popover close with Escape, and the
// help stays while the pointer moves onto it. WA2-AX-09: text boxes and selects focused with the
// keyboard have the 2 px outline of the other controls. WA2-AX-14: in forced colours (Windows'
// high contrast), the selected mode, tab and row are marked in a system colour.
import { expect, type Locator } from '@playwright/test';
import { test } from './support/test';
import { BUILD_HAS_ENGINE } from './support/browser';
import { openParameters, openSearch, paste, program } from './support/search';

test.setTimeout(120_000);

const REFUSAL = 'argv-validation-message';

/** Expects a control marked as refused, pointing to the refusal's message. */
async function expectRefused(control: Locator, message: string | RegExp): Promise<void> {
  await expect(control).toHaveAttribute('aria-invalid', 'true');
  await expect(control).toHaveAttribute('aria-describedby', REFUSAL);
  await expect(control).toHaveAccessibleDescription(message);
}

/** Expects a control with no mark of a refusal. */
async function expectAccepted(control: Locator): Promise<void> {
  await expect(control).not.toHaveAttribute('aria-invalid');
  await expect(control).not.toHaveAttribute('aria-describedby');
}

test('a refusal of "Additional options" marks its text box, which points to the message, until it is fixed', async ({ page }) => {
  await openSearch(page);
  await openParameters(page);
  const extra = page.getByTestId('param-extra');
  const validation = page.getByTestId('argv-validation');
  await expectAccepted(extra);
  // The Web's own message about a word typed there: the form's field of the same option is not marked.
  await extra.fill('-evalue abc');
  await expect(validation).toHaveText('-evalue is set in the form: Expect threshold');
  await expectRefused(extra, '-evalue is set in the form: Expect threshold');
  await expectAccepted(page.getByTestId('param-evalue'));
  await expect(page.locator(`#${REFUSAL}`)).toHaveText('-evalue is set in the form: Expect threshold');
  // The engine's refusal of an option typed there.
  await extra.fill('-qcov_hsp_perc 50 -foo');
  await expect(validation).toHaveText("error: unknown option or argument '-foo'; use -help for CLI v2 syntax");
  await expectRefused(extra, "error: unknown option or argument '-foo'; use -help for CLI v2 syntax");
  // A message about its text.
  await extra.fill("-strand 'plus");
  await expect(validation).toHaveText("Additional options: a quote (') is not closed");
  await expectRefused(extra, "Additional options: a quote (') is not closed");
  // Fixed: no mark.
  await extra.fill('-qcov_hsp_perc 50');
  await expect(validation).toHaveAttribute('data-state', 'ok');
  await expectAccepted(extra);
  await expect(page.locator(`#${REFUSAL}`)).toHaveCount(0);
});

test("a subrange's problem marks the fields that it is about, which point to it", async ({ page }) => {
  await openSearch(page);
  await paste(page, 'subject', `>s1\n${'ACGTACGTAC'.repeat(10)}\n`);
  const start = page.getByTestId('subject-region-start');
  const stop = page.getByTestId('subject-region-stop');
  const problem = page.getByTestId('subject-region-problem');
  const marked = async (field: Locator, text: string) => {
    await expect(field).toHaveAttribute('aria-invalid', 'true');
    await expect(field).toHaveAttribute('aria-describedby', 'subject-region-problem');
    await expect(field).toHaveAccessibleDescription(text);
  };
  await start.fill('5');
  await stop.fill('101');
  await expect(problem).toHaveText('The stop must be between 1 and 100, the length of the record.');
  await marked(stop, 'The stop must be between 1 and 100, the length of the record.');
  await expectAccepted(start);
  await start.fill('x');
  await expect(problem).toHaveText('The start and the stop are whole numbers.');
  await marked(start, 'The start and the stop are whole numbers.');
  await marked(stop, 'The start and the stop are whole numbers.');
  await stop.fill('60');
  await expectAccepted(stop);
  // A start after the stop: both.
  await start.fill('70');
  await expect(problem).toHaveText('The start cannot be larger than the stop.');
  await marked(start, 'The start cannot be larger than the stop.');
  await marked(stop, 'The start cannot be larger than the stop.');
  await start.fill('5');
  await expect(problem).toHaveCount(0);
  await expectAccepted(start);
  await expectAccepted(stop);
});

test("the engine's refusal of a form's value marks the field that holds it, until it is fixed", async ({ page }) => {
  test.skip(!BUILD_HAS_ENGINE, "the engine's refusals of values need the engine build");
  await openSearch(page);
  await program(page, 'blastn');
  await openParameters(page);
  await page.getByTestId('param-word_size').selectOption('other');
  const other = page.getByTestId('param-word_size-other');
  await other.fill('3');
  const message = /^error: invalid value '3' for '-word_size <WORD_SIZE>': expected an integer >= 4/;
  await expect(page.getByTestId('argv-validation')).toHaveText(message);
  await expectRefused(other, message);
  // The text box after "Other" holds the value; the select beside it and the other fields are not marked.
  await expectAccepted(page.getByTestId('param-word_size'));
  await expectAccepted(page.getByTestId('param-evalue'));
  await expectAccepted(page.getByTestId('param-extra'));
  await other.fill('');
  await expect(page.getByTestId('argv-validation')).toHaveAttribute('data-state', 'ok');
  await expectAccepted(other);
});
