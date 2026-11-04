// Why the selected query shows no HSP, in the words of the results' notices (ResultNotices):
// each notice starts with its sentence here, and the Dot Plot tab shows it in place of a plot.
import type { NoHspReason } from '../application/results';

export const NO_HSP_TEXT: Readonly<Record<NoHspReason, string>> = {
  'empty-query': 'The query has no sequences.',
  'no-tabular-hits': '-num_alignments 0 leaves the hit list of outfmt 6 and 7 empty.',
  'no-queries': 'No queries match the query filters.',
  'no-hits': 'No significant similarity found for this query.',
  'filtered-out': 'No HSPs of this query match the view filters.',
};
