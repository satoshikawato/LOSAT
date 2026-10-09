// Reading the outfmt 0 text that the engine wrote (plan §4.5). Only the subject heading's
// title is cut out for the subject list; the detail shows the heading and the HSP's
// section as written.

/**
 * The title of a subject heading (the bytes of an HSP's `out0_subject` range): the text
 * after the leading "> " up to the "Length=" line, as NCBI wrapped it. A list cell shows
 * it on one line (HTML white space), so the wrapping needs no change.
 */
export function headingTitle(heading: string): string {
  const end = heading.lastIndexOf('\nLength=');
  const title = end < 0 ? heading : heading.slice(0, end);
  return title.replace(/^>\s?/, '');
}
