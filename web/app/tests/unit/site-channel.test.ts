import { afterEach, describe, expect, it, vi } from 'vitest';
import { browserCountsChanges } from '../../src/infra/browser/counts-changes';
import { browserExampleDeletions } from '../../src/infra/browser/example-deletions';

// Review of line WQA's sweep, D5 (2026-11-19): Firefox, where the site's data is blocked (cookies
// blocked), refuses `new BroadcastChannel()` with SecurityError ("The operation is insecure."). The
// composition root opens the channels of the deletions, so the page did not start at all.
describe("the changes told to the site's other tabs", () => {
  afterEach(() => {
    vi.unstubAllGlobals();
  });

  const refusing = class {
    constructor() {
      throw new DOMException('The operation is insecure.', 'SecurityError');
    }
  };

  it.each([
    ['a browser that refuses BroadcastChannel', refusing],
    ['a browser without BroadcastChannel', undefined],
  ])('%s: told to no other tab, and nothing throws', (_case, channel) => {
    vi.stubGlobal('BroadcastChannel', channel);
    const examples = browserExampleDeletions();
    const counts = browserCountsChanges();
    examples.announce('pdbaa');
    counts.announce(9606);
    examples.subscribe(() => undefined)();
    counts.subscribe(() => undefined)();
  });

  it('a browser that opens the channel tells the other tabs', async () => {
    const heard = new Promise<string>((resolve) => browserExampleDeletions().subscribe(resolve));
    browserExampleDeletions().announce('pdbaa');
    expect(await heard).toBe('pdbaa');

    const heardCounts = new Promise<number>((resolve) => browserCountsChanges().subscribe(resolve));
    browserCountsChanges().announce(562);
    expect(await heardCounts).toBe(562);
  });
});
