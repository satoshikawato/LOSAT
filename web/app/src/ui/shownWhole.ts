// The one page of the results is cut as NCBI's (ClassicResults.vue): the Graphic Summary draws and
// the Descriptions list the first 100 subjects until "Show all N". Which query was shown whole is
// kept here, outside the components, so that it stays while the results tab is left (Search,
// Candidates) and shown again, as the results' tab does (AppView.vue); another run or query is
// cut again. Each part has its own "Show all", as on NCBI's page.
import { ref } from 'vue';

/** The run and query (`${runId}|${qIdx}`) shown whole, by part; '' for none. */
export const shownWhole = {
  graphic: ref(''),
  descriptions: ref(''),
};

/** The key of a run's query in `shownWhole`. */
export const wholeKey = (runId: string | undefined, qIdx: number | undefined): string => `${runId}|${qIdx}`;
