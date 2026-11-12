// DIAGNOSTIC ONLY (branch diag/ph-webkit-hang, line WFL item PH): not for merge. The commonest CI stall is at page start: the
// Data worker's describe ("loading") or a source's indexing ("indexing") never comes. Each test opens the page,
// waits for the describe, and pastes a query and a subject (indexing and the engine check), and nothing else.
import { test } from './support/test';
import { openSearch, paste, program } from './support/search';

test.setTimeout(60_000);

test('page start: describe, paste query and subject', async ({ page }) => {
  await openSearch(page);
  await program(page, 'blastn');
  await paste(page, 'query', '>q1\nACGTACGTACGTAGCTAGCTAGCATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCATCGAT\n');
  await paste(page, 'subject', '>s1\nACGTACGTACGTAGCTAGCTAGCATCGATCGATCGATCGTAGCTAGCTAGCTAGCTAGCATCGAT\n');
});
