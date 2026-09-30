// Contract suites (plan §2.1 L, §6.2 V-APP): lists of cases that every implementation of a
// port must pass. They do not depend on a test runner, so the same cases run in Vitest
// (Node), in browser workers under Playwright (tests/e2e/contracts.spec.ts), and in S09
// against the real reactor and Engine worker.

export interface ContractCase<Env> {
  readonly name: string;
  run(env: Env): Promise<void>;
}

export interface CaseResult {
  readonly name: string;
  readonly ok: boolean;
  readonly error?: string;
}

/**
 * Runs the cases one after another, each with a fresh environment. A case that does not
 * finish within `timeoutMs` fails, and the next case runs.
 */
export async function runCases<Env>(
  cases: readonly ContractCase<Env>[],
  makeEnv: () => Promise<Env> | Env,
  options: { readonly timeoutMs?: number; readonly onResult?: (result: CaseResult) => void } = {},
): Promise<CaseResult[]> {
  const timeoutMs = options.timeoutMs ?? 30_000;
  const results: CaseResult[] = [];
  for (const contractCase of cases) {
    let timer: ReturnType<typeof setTimeout> | undefined;
    const timeout = new Promise<never>((_, reject) => {
      timer = setTimeout(() => reject(new ContractError(`the case did not finish within ${timeoutMs} ms`)), timeoutMs);
    });
    try {
      await Promise.race([(async () => contractCase.run(await makeEnv()))(), timeout]);
      results.push({ name: contractCase.name, ok: true });
    } catch (error) {
      results.push({ name: contractCase.name, ok: false, error: describe(error) });
    } finally {
      clearTimeout(timer);
    }
    options.onResult?.(results[results.length - 1]!);
  }
  return results;
}

export class ContractError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'ContractError';
  }
}

export function check(condition: unknown, message: string): asserts condition {
  if (!condition) throw new ContractError(message);
}

export function same<T>(actual: T, expected: T, what: string): void {
  if (!Object.is(actual, expected)) {
    throw new ContractError(`${what}: expected ${show(expected)}, got ${show(actual)}`);
  }
}

/** Structural equality; object keys may be in any order. */
export function sameValue(actual: unknown, expected: unknown, what: string): void {
  if (!equalValues(actual, expected)) {
    throw new ContractError(`${what}: expected ${show(expected)}, got ${show(actual)}`);
  }
}

export function sameBytes(actual: Uint8Array, expected: Uint8Array, what: string): void {
  if (actual.length !== expected.length) {
    throw new ContractError(`${what}: expected ${expected.length} bytes, got ${actual.length}`);
  }
  for (let i = 0; i < actual.length; i++) {
    if (actual[i] !== expected[i]) {
      throw new ContractError(`${what}: byte ${i} is ${actual[i]}, expected ${expected[i]}`);
    }
  }
}

/** Expects `action` to throw an error whose "name: message" matches `pattern`. */
export function throws(action: () => unknown, pattern: RegExp, what: string): Error {
  try {
    action();
  } catch (error) {
    return matching(error, pattern, what);
  }
  throw new ContractError(`${what}: expected an error matching ${pattern}, but nothing was thrown`);
}

/** Expects `promise` to reject with an error whose "name: message" matches `pattern`. */
export async function rejects(promise: Promise<unknown>, pattern: RegExp, what: string): Promise<Error> {
  try {
    await promise;
  } catch (error) {
    return matching(error, pattern, what);
  }
  throw new ContractError(`${what}: expected a rejection matching ${pattern}, but it resolved`);
}

/** Resolves true if `promise` settles within `ms` milliseconds. */
export async function settlesWithin(promise: Promise<unknown>, ms: number): Promise<boolean> {
  let timer: ReturnType<typeof setTimeout> | undefined;
  const timeout = new Promise<false>((resolve) => {
    timer = setTimeout(() => resolve(false), ms);
  });
  const settled = promise.then(
    () => true,
    () => true,
  );
  try {
    return await Promise.race([settled, timeout]);
  } finally {
    clearTimeout(timer);
  }
}

/** Deterministic test bytes: byte `i` is `patternByte(i, seed)`. */
export function pattern(length: number, seed: number): Uint8Array {
  const bytes = new Uint8Array(length);
  for (let i = 0; i < length; i++) bytes[i] = patternByte(i, seed);
  return bytes;
}

export function patternByte(i: number, seed: number): number {
  return (i * 31 + seed * 17 + (i >>> 8)) & 0xff;
}

export function concatBytes(parts: readonly Uint8Array[]): Uint8Array {
  const bytes = new Uint8Array(parts.reduce((sum, part) => sum + part.length, 0));
  let offset = 0;
  for (const part of parts) {
    bytes.set(part, offset);
    offset += part.length;
  }
  return bytes;
}

export const MiB = 1024 * 1024;

function matching(error: unknown, pattern: RegExp, what: string): Error {
  const text = describe(error);
  if (!pattern.test(text)) throw new ContractError(`${what}: expected an error matching ${pattern}, got "${text}"`);
  return error instanceof Error ? error : new Error(text);
}

function describe(error: unknown): string {
  return error instanceof Error ? `${error.name}: ${error.message}` : String(error);
}

function show(value: unknown): string {
  const text = JSON.stringify(value);
  return text === undefined ? String(value) : text.length > 200 ? `${text.slice(0, 200)}...` : text;
}

function equalValues(a: unknown, b: unknown): boolean {
  if (Object.is(a, b)) return true;
  if (typeof a !== 'object' || typeof b !== 'object' || a === null || b === null) return false;
  if (Array.isArray(a) !== Array.isArray(b)) return false;
  const keysA = Object.keys(a);
  const keysB = Object.keys(b);
  if (keysA.length !== keysB.length) return false;
  return keysA.every(
    (key) =>
      Object.prototype.hasOwnProperty.call(b, key) &&
      equalValues((a as Record<string, unknown>)[key], (b as Record<string, unknown>)[key]),
  );
}
