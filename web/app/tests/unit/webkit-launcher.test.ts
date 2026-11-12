// The launcher through which the E2E tests start WebKit on Linux (tests/e2e/support/webkit-launcher.sh,
// playwright.config.ts): the program it starts has a core file limit of 1 byte, so that the kernel writes no
// core for a crashing WebKit process even where cores go to a program through a pipe (line WFL item PH).
import { execFileSync, spawnSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';

const LAUNCHER = fileURLToPath(new URL('../e2e/support/webkit-launcher.sh', import.meta.url));

describe.runIf(process.platform === 'linux')('the WebKit launcher of the E2E tests', () => {
  it('starts the launcher it is given, with its arguments, and a core file limit of 1 byte', () => {
    const output = execFileSync(LAUNCHER, ['-c', 'grep "^Max core file size" /proc/self/limits; printf "<%s>" "$@"', 'sh', 'a', 'b c'], {
      env: { ...process.env, LOSAT_WEB_WEBKIT_REAL: '/bin/sh' },
    }).toString();
    // The soft limit, which the kernel reads when a process dumps core.
    expect(output).toMatch(/^Max core file size\s+1\s+\S+\s+bytes\s*$/m);
    expect(output).toContain('<a><b c>');
  });

  it('starts the launcher without a limit where the hard limit is 0, which writes no core anyway', () => {
    const output = execFileSync('prlimit', ['--core=0:0', '--', LAUNCHER, '-c', 'grep "^Max core file size" /proc/self/limits'], {
      env: { ...process.env, LOSAT_WEB_WEBKIT_REAL: '/bin/sh' },
    }).toString();
    expect(output).toMatch(/^Max core file size\s+0\s+0\s+bytes\s*$/m);
  });

  it('says how to install WebKit when the launcher it is given does not exist', () => {
    const result = spawnSync(LAUNCHER, [], { env: { ...process.env, LOSAT_WEB_WEBKIT_REAL: '/nonexistent/pw_run.sh' } });
    expect(result.status).toBe(127);
    expect(result.stderr.toString()).toContain('npx playwright install webkit');
  });

  it('starts nothing when it is not given a launcher', () => {
    const env = { ...process.env };
    delete env['LOSAT_WEB_WEBKIT_REAL'];
    const result = spawnSync(LAUNCHER, [], { env });
    expect(result.status).not.toBe(0);
    expect(result.stderr.toString()).toContain('LOSAT_WEB_WEBKIT_REAL');
  });
});
