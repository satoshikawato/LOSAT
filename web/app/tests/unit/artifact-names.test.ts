// The names a test leaves for CI's artifact upload (tests/e2e/support/artifact-names.ts): the upload
// refuses a whole artifact for one name with " : < > | * ? CR or LF, so every e2e test fails on such a
// name, and evidence keeps its files under uploadable names.
import { mkdirSync, mkdtempSync, rmSync, writeFileSync } from 'node:fs';
import { tmpdir } from 'node:os';
import { join } from 'node:path';
import { afterEach, describe, expect, it } from 'vitest';
import { unuploadableNames, uploadableName } from '../e2e/support/artifact-names';
import { makeFixtures } from '../e2e/walk/fixtures';

const folders: string[] = [];
afterEach(() => {
  for (const folder of folders.splice(0)) rmSync(folder, { recursive: true, force: true });
});

function folder(): string {
  const made = mkdtempSync(join(tmpdir(), 'artifact-names-'));
  folders.push(made);
  return made;
}

describe('artifact names', () => {
  it('finds every name the upload refuses, in nested folders, and none of the others', () => {
    const root = folder();
    mkdirSync(join(root, 'native', 'run-1'), { recursive: true });
    for (const name of ['a"b.fa', 'c:d', 'e<f', 'g>h', 'i|j', 'k*l', 'm?n', 'o\rp', 'q\nr']) writeFileSync(join(root, 'native', 'run-1', name), '');
    for (const name of ["it's.faa", 'two genes.fna', 'a%3Cb.fa', 'trace.zip']) writeFileSync(join(root, 'native', 'run-1', name), '');
    mkdirSync(join(root, 'bad|folder'));
    expect(unuploadableNames(root).sort()).toEqual(
      ['bad|folder', ...['a"b.fa', 'c:d', 'e<f', 'g>h', 'i|j', 'k*l', 'm?n', 'o\rp', 'q\nr'].map((name) => join('native', 'run-1', name))].sort(),
    );
    expect(unuploadableNames(join(root, 'missing'))).toEqual([]);
  });

  it('gives every walk fixture a name the upload keeps, and two names two kept names', () => {
    const files = makeFixtures(20261119).files.map((file) => file.name);
    expect(files.some((name) => /[<>"]/.test(name))).toBe(true);
    const root = folder();
    const kept = files.map(uploadableName);
    for (const name of kept) writeFileSync(join(root, name), '');
    expect(unuploadableNames(root)).toEqual([]);
    expect(new Set(kept).size).toBe(new Set(files).size);
    expect(uploadableName('<svg onload="x">.fa')).toBe('%3Csvg onload=%22x%22%3E.fa');
    expect(uploadableName('100%3C.fa')).toBe('100%253C.fa');
    expect(uploadableName("it's.faa")).toBe("it's.faa");
  });
});
