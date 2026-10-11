// Writes the preview check's record: `<project>.json` for each project (http, chromium,
// firefox, webkit: every test with its status, duration, error and the record the test
// attached) and `summary.json` (one line per test). The records hold URLs, statuses and
// response headers of the site; no research data.
import { mkdirSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';
import type { FullResult, Reporter, TestCase, TestResult } from '@playwright/test/reporter';

interface Entry {
  readonly project: string;
  readonly title: string;
  readonly status: string;
  readonly durationMs: number;
  readonly skipReason?: string;
  readonly error?: string;
  readonly record?: unknown;
}

// eslint-disable-next-line no-control-regex
const ANSI = /\u001b\[[0-9;]*m/g;

export default class RecordReporter implements Reporter {
  private readonly entries: Entry[] = [];

  constructor(private readonly options: { out: string; origin: string }) {}

  onTestEnd(test: TestCase, result: TestResult): void {
    const attached = result.attachments.find((attachment) => attachment.name === 'record' && attachment.body !== undefined);
    const skip = test.annotations.find((annotation) => annotation.type === 'skip');
    const error = result.errors[0]?.message;
    this.entries.push({
      project: test.parent.project()?.name ?? 'unknown',
      title: test.title,
      status: result.status,
      durationMs: result.duration,
      ...(skip?.description === undefined || skip.description === '' ? {} : { skipReason: skip.description }),
      ...(error === undefined ? {} : { error: error.replace(ANSI, '').slice(0, 2000) }),
      ...(attached === undefined ? {} : { record: JSON.parse(attached.body!.toString('utf8')) as unknown }),
    });
  }

  onEnd(result: FullResult): void {
    mkdirSync(this.options.out, { recursive: true });
    const projects = [...new Set(this.entries.map((entry) => entry.project))];
    for (const project of projects) {
      const tests = this.entries.filter((entry) => entry.project === project);
      writeFileSync(join(this.options.out, `${project}.json`), `${JSON.stringify({ project, site: this.options.origin, tests }, null, 2)}\n`);
    }
    const summary = {
      site: this.options.origin,
      status: result.status,
      counts: Object.fromEntries(
        projects.map((project) => {
          const tests = this.entries.filter((entry) => entry.project === project);
          const count = (status: string) => tests.filter((entry) => entry.status === status).length;
          return [project, { passed: count('passed'), failed: count('failed') + count('timedOut'), skipped: count('skipped') }];
        }),
      ),
      tests: this.entries.map(({ project, title, status, skipReason, error }) => ({
        project,
        title,
        status,
        ...(skipReason === undefined ? {} : { skipReason }),
        ...(error === undefined ? {} : { error: error.split('\n')[0] }),
      })),
    };
    writeFileSync(join(this.options.out, 'summary.json'), `${JSON.stringify(summary, null, 2)}\n`);
  }
}
