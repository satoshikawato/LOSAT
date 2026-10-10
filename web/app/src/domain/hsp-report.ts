// The results report of LOSAT Web (S15 item 2; design §12.3): one static, self-contained HTML page
// that shows a run's results without searching again. It has no script, no external URL, no image,
// no form and no link; its Content-Security-Policy allows nothing but its inline style, so a page
// saved with hostile IDs still loads nothing and runs nothing. It is an application format made
// from the results as the engine wrote them, not an NCBI BLAST report: the HSP tables hold the
// outfmt 6 fields and the alignments the outfmt 0 headings and sections, as written. Every text
// that comes from the run (IDs, names, titles, argv, outfmt texts, warnings) goes through
// `escapeHtml`; the page's own tags and attributes are the constants here. The functions make the
// parts of the page; application/result-export.ts writes them in order, in blocks.
import { toShellCommand } from './argv';
import type { Unit } from './coordinates';
import { hspLabel } from './extraction';
import { filterWords, SCOPE_LABELS, type ExportedHsp, type ExportInput, type ExportRun, type ExportScopeInfo } from './hsp-export';
import { OUTFMT6_FIELDS } from './outfmt6';
import type { OutputFormat } from './output-format';

const ENTITIES: Readonly<Record<string, string>> = { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' };

const MARKUP = /[&<>"']/;

/** Text as HTML text: `& < > " '` become entities, so no text of the run is read as markup. */
export function escapeHtml(text: string): string {
  // Most fields of a report hold none of these characters: they are returned as they are.
  return MARKUP.test(text) ? text.replace(/[&<>"']/g, (c) => ENTITIES[c]!) : text;
}

/** The page loads nothing and runs nothing: only its inline style applies. */
export const REPORT_CSP = "default-src 'none'; style-src 'unsafe-inline'; base-uri 'none'; form-action 'none'";

const STYLE = `
:root { color-scheme: light dark; --fg: #1f2328; --bg: #ffffff; --muted: #59636e; --line: #d1d9e0; --note: #fff8c5; }
@media (prefers-color-scheme: dark) { :root { --fg: #e6edf3; --bg: #0d1117; --muted: #9198a1; --line: #3d444d; --note: #3b2e00; } }
body { margin: 0 auto; max-width: 76rem; padding: 16px; color: var(--fg); background: var(--bg);
  font: 14px/1.45 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
h1 { font-size: 1.5rem; margin: 0 0 0.5rem; }
h2 { font-size: 1.15rem; margin: 1.5rem 0 0.5rem; border-bottom: 1px solid var(--line); padding-bottom: 0.25rem; }
h3 { font-size: 1rem; margin: 1rem 0 0.25rem; }
.statement { background: var(--note); border: 1px solid var(--line); padding: 0.5rem 0.75rem; }
.muted { color: var(--muted); }
.notice { background: var(--note); padding: 0.25rem 0.5rem; }
dl { display: grid; grid-template-columns: max-content minmax(0, 1fr); gap: 0.15rem 1rem; margin: 0; }
dt { font-weight: 600; }
dd { margin: 0; overflow-wrap: anywhere; }
code, pre { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 12px; }
code { overflow-wrap: anywhere; }
pre { margin: 0.25rem 0; overflow-x: auto; white-space: pre; }
.table { overflow-x: auto; }
table { border-collapse: collapse; font-size: 12px; }
th, td { border: 1px solid var(--line); padding: 2px 6px; text-align: left; white-space: nowrap; }
th { background: var(--note); }
.query { margin-top: 2rem; }
.hsp-label { margin: 0.5rem 0 0; font-weight: 600; }
footer { margin-top: 2rem; }
`;

/** The parts of the page's top: the run, the commands, the verification badge and the scope. */
export interface ReportHead {
  readonly run: ExportRun;
  /** The output formats of the program (the commands). */
  readonly formats: readonly OutputFormat[];
  readonly scope: ExportScopeInfo;
  readonly exportedAt: number;
  /** The outfmt 0 headings and sections are in the report (default true). */
  readonly alignments?: boolean;
}

const iso = (ms: number): string => new Date(ms).toISOString();
const item = (term: string, value: string): string => `<dt>${term}</dt><dd>${value}</dd>\n`;
const count = (n: number, one: string): string => `${n} ${n === 1 ? one : `${one}s`}`;

function inputItem(role: string, input: ExportInput): string {
  return item(
    role,
    `${escapeHtml(input.name)}: ${count(input.records, 'record')}, ${input.bytes} bytes<br>` +
      `<span class="muted">SHA-256 of the FASTA given to the engine: <code>${escapeHtml(input.sha256)}</code></span>`,
  );
}

const ALIGNMENTS_INCLUDED = 'Included: the outfmt 0 headings and sections of the HSPs that outfmt 0 shows, as written.';
const ALIGNMENTS_LEFT_OUT =
  'Not included: this report was saved without the outfmt 0 text of the alignments. Each HSP table says whether outfmt 0 shows the HSP.';

function scopeItems(scope: ExportScopeInfo): string {
  let items = item('HSPs', `${escapeHtml(SCOPE_LABELS[scope.scope])}: ${count(scope.hsps, 'HSP')}`);
  if (scope.scope === 'all') {
    items += item('View filters', 'Not applied: every HSP of the run, in the engine’s order.');
  } else {
    const words = scope.filters === undefined ? [] : filterWords(scope.filters);
    items += item('View filters', words.length === 0 ? 'None set' : words.map(escapeHtml).join('; '));
  }
  if (scope.query !== undefined) items += item('Query', `${scope.query.position + 1}: ${escapeHtml(scope.query.id)}`);
  if (scope.markedSubjects !== undefined) {
    items += item('Marked subjects', scope.markedSubjects.map((s) => `${s.position + 1}: ${escapeHtml(s.id)}`).join('; '));
  }
  return items;
}

/** The page up to the first query: the statement, the run, the commands, the verification and the scope. */
export function reportHead(head: ReportHead): string {
  const { run, scope } = head;
  const engine = [run.runtimePath, run.threads === undefined ? undefined : count(run.threads, 'thread')].filter((w) => w !== undefined).join(', ');
  const times = [
    item('Queued', iso(run.queuedAt)),
    run.startedAt === undefined ? '' : item('Started', iso(run.startedAt)),
    run.endedAt === undefined ? '' : item('Ended', iso(run.endedAt)),
  ].join('');
  const commands = head.formats.map((format) => `<li>outfmt ${format}: <code>${escapeHtml(toShellCommand(run.argv, format))}</code></li>\n`).join('');
  return (
    '<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n' +
    `<meta http-equiv="Content-Security-Policy" content="${REPORT_CSP}">\n` +
    '<meta name="referrer" content="no-referrer">\n<meta name="viewport" content="width=device-width, initial-scale=1">\n' +
    `<title>LOSAT Web report, run ${run.number}</title>\n<style>${STYLE}</style>\n</head>\n<body>\n` +
    '<header>\n<h1>LOSAT Web report</h1>\n' +
    '<p class="statement"><strong>An application format of LOSAT Web.</strong> LOSAT Web made this page from the results ' +
    'of one run as the engine wrote them: the outfmt 6 rows, the outfmt 0 headings and sections, and the HSP records. ' +
    'It is not an NCBI BLAST report and was not compared with NCBI BLAST+. The page has no script and loads nothing.</p>\n' +
    `<p class="muted">Saved ${iso(head.exportedAt)}.</p>\n</header>\n` +
    '<section>\n<h2>Run</h2>\n<dl>\n' +
    item('Run', `${run.number}${run.group === undefined ? '' : `, group run ${run.group.position} of ${run.group.size}`}`) +
    (run.title === undefined ? '' : item('Job Title', escapeHtml(run.title))) +
    item('Program', `${escapeHtml(run.programLabel)}${run.task === undefined ? '' : ` (task ${escapeHtml(run.task)})`}`) +
    item('Arguments', `<code>${escapeHtml(run.argv.join(' '))}</code>`) +
    inputItem('Query', run.query) +
    inputItem('Subject', run.subject) +
    item('Threads requested', run.requestedThreads === 'auto' ? 'Auto' : String(run.requestedThreads)) +
    (engine === '' ? '' : item('Engine', escapeHtml(engine))) +
    (run.engineBuild === undefined ? '' : item('Build', escapeHtml(run.engineBuild))) +
    times +
    '</dl>\n</section>\n' +
    '<section>\n<h2>Verification</h2>\n' +
    `<p><strong>${escapeHtml(run.verification.label)}</strong></p>\n<ul>\n` +
    run.verification.details.map((line) => `<li>${escapeHtml(line)}</li>\n`).join('') +
    '</ul>\n' +
    run.verification.exceptions.map((line) => `<p class="notice">${escapeHtml(line)}</p>\n`).join('') +
    '<p class="muted">The badge speaks about the engine build and the run’s options; this page is a LOSAT Web format either way.</p>\n' +
    '</section>\n' +
    '<section>\n<h2>Command line</h2>\n' +
    '<p class="muted">The LOSAT command that writes each output from the same inputs (the names of the run’s -query and -subject).</p>\n' +
    `<ul>\n${commands}</ul>\n</section>\n` +
    '<section>\n<h2>Scope</h2>\n<dl>\n' +
    scopeItems(scope) +
    item('Alignments', head.alignments === false ? ALIGNMENTS_LEFT_OUT : ALIGNMENTS_INCLUDED) +
    '</dl>\n</section>\n'
  );
}

export interface ReportQuery {
  /** 0-based position of the query record in the run's input. */
  readonly position: number;
  readonly id: string;
  readonly length: number;
  readonly unit: Unit;
  /** The query's HSPs in the scope. */
  readonly hsps: number;
}

const TABLE_START =
  '<h3>HSPs (outfmt 6 fields, as written)</h3>\n<div class="table"><table>\n<thead><tr><th>HSP</th><th>Subject record</th>' +
  OUTFMT6_FIELDS.map((name) => `<th>${name}</th>`).join('') +
  '<th>outfmt 0</th></tr></thead>\n<tbody>\n';

/** A query's section, up to the rows of its HSP table. */
export function reportQueryStart(query: ReportQuery): string {
  const name = query.id === '' ? '' : `: ${escapeHtml(query.id)}`;
  return (
    `<section class="query">\n<h2>Query ${query.position + 1}${name}</h2>\n` +
    `<p class="muted">Length ${query.length} ${query.unit}; ${count(query.hsps, 'HSP')} in this report.</p>\n` +
    TABLE_START
  );
}

/** A row of the query's HSP table: the HSP, its subject record and its outfmt 6 fields as written. */
export function reportTableRow(hsp: ExportedHsp): string {
  let fields = '';
  for (const name of OUTFMT6_FIELDS) fields += `<td>${escapeHtml(hsp.outfmt6?.[name] ?? '')}</td>`;
  return `<tr><td>${hspLabel(hsp.qIdx, hsp.rank)}</td><td>${hsp.sIdx + 1}</td>${fields}<td>${hsp.inOutfmt0 ? 'shown' : 'not shown'}</td></tr>\n`;
}

export const REPORT_TABLE_END = '</tbody>\n</table></div>\n';
export const REPORT_ALIGNMENTS_START = '<h3>Alignments (outfmt 0, as written)</h3>\n';

/**
 * Text as written in a `<pre>`: the line end after the start tag is one that the HTML parser drops,
 * so a text that starts with a line end keeps it.
 */
const pre = (className: string, text: string): string => `<pre${className === '' ? '' : ` class="${className}"`}>\n${escapeHtml(text)}</pre>\n`;

/** A subject heading of outfmt 0 (defline and Length=), as written. */
export function reportHeading(heading: string): string {
  return pre('heading', heading);
}

/** An HSP's section of outfmt 0 (score lines and alignment), as written. */
export function reportSection(hsp: ExportedHsp, section: string): string {
  return `<p class="hsp-label">HSP ${hspLabel(hsp.qIdx, hsp.rank)}</p>\n${pre('section', section)}`;
}

/** The query's HSPs that outfmt 0 does not show, by label. */
export function reportNotInOutfmt0(labels: readonly string[]): string {
  return labels.length === 0 ? '' : `<p class="muted">outfmt 0 does not show ${labels.length === 1 ? 'HSP' : 'HSPs'} ${labels.join(', ')}.</p>\n`;
}

export const REPORT_QUERY_END = '</section>\n';

/** The page after the last query: the run's warnings, as the engine wrote them. */
export function reportTail(diagnostics: string): string {
  const warnings = diagnostics === '' ? '<p class="muted">The engine wrote no warnings.</p>\n' : pre('', diagnostics);
  return (
    `<section>\n<h2>Warnings</h2>\n${warnings}</section>\n` +
    '<footer class="muted"><p>LOSAT Web report: an application format of LOSAT Web, not an NCBI BLAST report.</p></footer>\n' +
    '</body>\n</html>\n'
  );
}
