// Prints the tables of RESULTS.md from results/<project>.json (node summarize.mjs).
import { readFileSync } from 'node:fs';

const MECHANICS = ['chromium', 'firefox', 'webkit'];
const STORAGE = ['chromium', 'chromium-stock-storage', 'firefox', 'firefox-stock-storage', 'webkit'];
const load = (name) => {
  try {
    return JSON.parse(readFileSync(new URL(`results/${name}.json`, import.meta.url), 'utf8'));
  } catch {
    return undefined;
  }
};
const results = Object.fromEntries([...new Set([...MECHANICS, ...STORAGE])].map((name) => [name, load(name)]));
const row = (project, id) => results[project]?.rows.find((r) => r.id === id);
const yes = (value) => (value ? 'yes' : 'no');

function mechanics(r) {
  if (r === undefined) return 'not run';
  if (r.id === 'b') {
    const exact = r.body !== null && /^\{"v":1,"release":"probe","client":"[0-9a-f-]{36}"\}$/.test(r.body);
    return `${r.visit?.status ?? r.visit?.error}; reached ${r.reached}; body exact ${yes(exact)}; Referer ${r.headers?.referer ?? 'none'}; Cookie ${r.headers?.cookie ?? 'none'}; CSP violations ${r.aCspViolations.length}`;
  }
  if (r.id === 'navigation') return `after b: ${r.b}; after a2: ${r.a2}`;
  if (!r.frameLoaded) {
    const csp = r.aCspViolations.map((v) => v.directive).join(',');
    if (r.hello === null && r.messages?.length) return 'no hello';
    return `frame blocked${csp ? ` (A's CSP: ${csp})` : r.mRequested ? ' (M answered; refused by COEP)' : ''}`;
  }
  const g = r.gScript?.ran ? 'gtag ran' : `gtag not run${r.mCspViolations?.length ? ` (M's CSP: ${r.mCspViolations.map((v) => v.directive).join(',')})` : ''}`;
  const c = r.collectReached.length ? `C got ${r.collectReached.join('+')}` : 'C got nothing';
  const checks = r.hello?.checksPass && r.resultChecksPass ? 'checks pass' : `checks fail (origin seen ${r.hello?.origin})`;
  return `frame loads; ${g}; ${c}; ${checks}`;
}

function opaque(r) {
  if (r === undefined) return 'not run';
  if (!r.frameLoaded && r.hello === null) return 'no hello reached A';
  const runs = (r.mRunMessages ?? []).map((m) => `${m.target}`).join(',') || 'none';
  return `hello origin "${r.hello.origin}" (check ${r.hello.checksPass ? 'passes' : 'refuses'}); M got run by: ${runs}; localStorage ${r.storage?.localError ?? 'ok'}`;
}

function identifier(r) {
  if (r === undefined) return 'not run';
  if (!r.frameLoaded) return 'frame blocked';
  const p = r.persistence;
  if (p === undefined) return 'not checked';
  const ls = r.storage?.localError ? `LS ${r.storage.localError}` : `LS kept ${yes(p.localKept)}${p.localKept ? (p.localSameAsTopLevel ? ' (unpartitioned)' : ' (partitioned)') : ''}`;
  const none = r.storage?.cookieAfter?.probe_none ? `cookie None kept ${yes(p.cookieNoneKept)}` : 'cookie None not set';
  const lax = r.storage?.cookieAfter?.probe_lax ? `cookie Lax kept ${yes(p.cookieLaxKept)}` : 'cookie Lax not set';
  return `${ls}; ${none}; ${lax}`;
}

const ids = ['b', 'a1', 'a2', 'a3', 'a3n', 'a4', 'a2-csp-unchanged', 'a2-no-sandbox', 'navigation'];
const what = (id) => MECHANICS.map((p) => row(p, id)?.what).find(Boolean) ?? id;
console.log(`| Variant | ${MECHANICS.map((p) => `${p} ${results[p]?.version ?? ''}`).join(' | ')} |`);
console.log(`|---|${MECHANICS.map(() => '---').join('|')}|`);
for (const id of ids) console.log(`| ${id}: ${what(id)} | ${MECHANICS.map((p) => mechanics(row(p, id))).join(' | ')} |`);
console.log();
console.log(`| Opaque sandbox (a2, sandbox="allow-scripts") | ${MECHANICS.map((p) => opaque(row(p, 'a2-sandbox-opaque'))).join(' | ')} |`);
console.log();
console.log(`| Identifier in M's frame, kept when A is opened again | ${STORAGE.join(' | ')} |`);
console.log(`|---|${STORAGE.map(() => '---').join('|')}|`);
for (const id of ['a2', 'a3', 'a4']) console.log(`| ${id} | ${STORAGE.map((p) => identifier(row(p, id))).join(' | ')} |`);
console.log();
for (const p of STORAGE) if (results[p]) console.log(`- ${p} ${results[p].version}: ${results[p].storageRules}`);
