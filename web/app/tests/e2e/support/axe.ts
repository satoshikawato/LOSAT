// axe-core (a devDependency of the same version as the second audit's, W5) on parts of the page.
// The page's CSP (script-src 'self') refuses an inline script tag, so axe's source runs through
// the driver's evaluate, which is not a page script.
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import type { Page } from '@playwright/test';
import { REPOSITORY } from './harness-server';

let source: string | undefined;

/** A rule that a part of the page breaks, with the elements that break it. */
export interface AxeViolation {
  readonly id: string;
  readonly impact: string | null;
  readonly targets: readonly string[];
}

/** The violations of axe's rules (its default set) in the elements that `selectors` select, without those that `exclude` selects. */
export async function axeViolations(page: Page, selectors: readonly string[], exclude: readonly string[] = []): Promise<AxeViolation[]> {
  source ??= readFileSync(join(REPOSITORY, 'web/app/node_modules/axe-core/axe.min.js'), 'utf8');
  if (!(await page.evaluate(() => 'axe' in window))) await page.evaluate(`${source}\n;void 0`);
  return page.evaluate(async ([include, left]) => {
    type Result = { violations: { id: string; impact?: string | null; nodes: { target: unknown[] }[] }[] };
    const axe = (window as unknown as { axe: { run(context: unknown, options: unknown): Promise<Result> } }).axe;
    const result = await axe.run(left.length === 0 ? { include } : { include, exclude: left }, { resultTypes: ['violations'] });
    return result.violations.map((violation) => ({
      id: violation.id,
      impact: violation.impact ?? null,
      targets: violation.nodes.map((node) => node.target.join(' ')),
    }));
  }, [[...selectors], [...exclude]] as const);
}
