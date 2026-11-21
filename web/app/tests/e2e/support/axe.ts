// axe-core (a devDependency of the same version as the second audit's, W5) on parts of the page.
// The page's CSP (script-src 'self') refuses an inline script tag, so axe's source runs through
// the driver's evaluate, which is not a page script.
import { readFileSync } from 'node:fs';
import { createRequire } from 'node:module';
import type { Page } from '@playwright/test';

let source: string | undefined;

/** A rule that a part of the page breaks, with the elements that break it. */
export interface AxeViolation {
  readonly id: string;
  readonly impact: string | null;
  readonly targets: readonly string[];
}

/** The violations of axe's rules (its default set) in the elements that `selectors` select, without those that `exclude` selects. */
export async function axeViolations(page: Page, selectors: readonly string[], exclude: readonly string[] = []): Promise<AxeViolation[]> {
  // Where npm installed it (the workspace's node_modules, web/package.json).
  source ??= readFileSync(createRequire(import.meta.url).resolve('axe-core/axe.min.js'), 'utf8');
  if (!(await page.evaluate(() => 'axe' in window))) {
    await page.evaluate(`${source}\n;void 0`);
    await page.evaluate(readIdrefsFromAttributes);
  }
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

/**
 * Makes axe read the ID references (aria-labelledby, aria-controls, ...) from the attributes only,
 * not also from ARIA element reflection (`element.ariaLabelledByElements` and the like). In Firefox
 * 159 (Playwright 1.65.0, Firefox build 1572), the first read of `ariaLabelledByElements` of the
 * results panel that "Run LOSAT" opened killed the page's content process (SIGSEGV, "Page
 * crashed"), in every run (PR #1308, 2026-11-19); Firefox 157 (build 1555) and the next Firefox 159
 * build (1578, Playwright 1.66.0-alpha) do not. axe reads it for aria-allowed-role (the implicit
 * role of a <section> depends on its name). The same page without that read did not crash (also
 * with GNOME_ACCESSIBILITY=1). The application sets every ID reference as an attribute and uses no
 * reflected ARIA property (tests/unit/aria-references.test.ts), so axe checks the same references
 * either way, in every engine.
 */
function readIdrefsFromAttributes(): void {
  type Standards = { ariaAttrs: Record<string, { prop?: string }> };
  const axe = (window as unknown as { axe: { utils: { getStandards(): Standards }; configure(spec: unknown): void } }).axe;
  const reflected = Object.entries(axe.utils.getStandards().ariaAttrs).filter(([, attr]) => attr.prop);
  // An empty name: axe's deepMerge would turn null into {}, which axe would take for a property name.
  axe.configure({ standards: { ariaAttrs: Object.fromEntries(reflected.map(([name]) => [name, { prop: '' }])) } });
}
