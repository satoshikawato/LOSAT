// The dot plot's SVG file (S15, WP-E): the plot's text blocks (domain/plot-svg.ts) written through
// the Writer, one block after another, so that a pair with thousands of HSPs is never one string.
import { dotPlotFileName, dotPlotSvg, SVG_MIME, type DotPlotSvgInput } from '../domain/plot-svg';
import type { Downloader } from '../ports/download';
import { writeFile } from './export-writer';

/**
 * Saves the plot as `losat-run{N}-dotplot-q{Q}-s{S}.svg` (Q and S are the 1-based positions of the
 * query and subject records). A failure saves nothing and is thrown again; resolves with the file's
 * length in bytes.
 */
export function exportDotPlotSvg(
  downloader: Pick<Downloader, 'open'>,
  plot: DotPlotSvgInput,
  queryPosition: number,
  subjectPosition: number,
): Promise<number> {
  return writeFile(downloader, dotPlotFileName(plot.run, queryPosition, subjectPosition), SVG_MIME, async (writer) => {
    for (const block of dotPlotSvg(plot)) await writer.text(block);
  });
}
