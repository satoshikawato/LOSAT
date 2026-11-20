// The words of a download, shared by the lines of the Example sets and of the
// WindowMasker counts (application/examples.ts, application/counts.ts): the site's host and the
// reason, in the words of the in-place sentences (docs/web/ncbi_ui_mapping.md §4).

/** No answer: offline, blocked, timed out, or a connection that ended early. */
export function unreachableSentence(host: string, detail: string): string {
  return `${host} could not be reached (${detail}).`;
}

/** An HTTP status other than the file: a 5xx may pass on another try. */
export function refusedSentence(host: string, status: number): string {
  return status >= 500 ? `${host} answered HTTP ${status}: try again in a moment.` : `${host} answered HTTP ${status}.`;
}

/** The date (yyyy-mm-dd) of an HTTP date such as Last-Modified, if it is one. */
export function dateOf(httpDate: string | null | undefined): string | undefined {
  const time = Date.parse(httpDate ?? '');
  return Number.isNaN(time) ? undefined : new Date(time).toISOString().slice(0, 10);
}
