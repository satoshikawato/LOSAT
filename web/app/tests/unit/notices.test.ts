import { existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { NOTICE_FILES, npmLicensesText, npmPackages } from '../../build/notices';

const REPOSITORY = fileURLToPath(new URL('../../../..', import.meta.url));
const MODULES = fileURLToPath(new URL('../../../node_modules', import.meta.url));

describe('the license and the notices of the site', () => {
  it('copies files that the repository has', () => {
    for (const file of NOTICE_FILES) expect(existsSync(`${REPOSITORY}/${file.source}`), file.source).toBe(true);
  });

  it('names the packages of the modules in node_modules and of the helpers that Vite adds', () => {
    const packages = npmPackages([
      `${MODULES}/@vue/shared/dist/shared.esm-bundler.js`,
      `${MODULES}/@bjorn3/browser_wasi_shim/dist/wasi.js`,
      '\0vite/modulepreload-polyfill.js',
      '\0virtual:losat-engine',
      `${REPOSITORY}/web/app/src/main.ts`,
    ]);
    expect(packages.map((p) => p.name)).toEqual(['@bjorn3/browser_wasi_shim', '@vue/shared', 'rolldown', 'vite']);
  });

  it('stops at a virtual module whose owner it does not know', () => {
    expect(() => npmPackages(['\0some-plugin:helper'])).toThrow(/VIRTUAL_OWNERS/);
  });

  it('prints every license file of a package, and only the core part of Vite\'s', () => {
    const text = npmLicensesText(npmPackages([`${MODULES}/@bjorn3/browser_wasi_shim/dist/wasi.js`, '\0vite/modulepreload-polyfill.js']));
    expect(text).toContain('LICENSE-MIT:');
    expect(text).toContain('LICENSE-APACHE:');
    expect(text).toContain('Author: bjorn3');
    expect(text).toContain('Copyright (c) 2019-present, VoidZero Inc. and Vite contributors');
    expect(text).not.toContain('Licenses of bundled dependencies');
  });
});
