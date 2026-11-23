// Writes src/domain/html-named-references.ts, the named character references of HTML that the
// viewer decodes in the outfmt 0 text of a run made with -html (domain/outfmt0-html.ts,
// htmlShownText), from the table of the WHATWG HTML Living Standard (section "Named character
// references"; https://html.spec.whatwg.org/entities.json). The standard says the table will
// never change. Run after a download of that file: node scripts/html-named-references.mjs entities.json
import { createHash } from 'node:crypto';
import { readFileSync, writeFileSync } from 'node:fs';
import { resolve } from 'node:path';

const [source] = process.argv.slice(2);
if (source === undefined) throw new Error('usage: node scripts/html-named-references.mjs <entities.json>');
const bytes = readFileSync(source);
const table = JSON.parse(bytes.toString('utf8'));
const OUT = resolve(import.meta.dirname, '../src/domain/html-named-references.ts');

/** A JavaScript string literal of `text` with every character outside printable ASCII escaped. */
function literal(text) {
  let out = "'";
  for (const char of text) {
    const code = char.codePointAt(0);
    if (char === "'" || char === '\\') out += `\\${char}`;
    else if (code >= 0x20 && code < 0x7f) out += char;
    else out += `\\u{${code.toString(16)}}`;
  }
  return `${out}'`;
}

const withSemicolon = new Map();
const legacy = [];
for (const [identifier, { characters }] of Object.entries(table)) {
  const match = /^&([A-Za-z][A-Za-z0-9]*)(;?)$/.exec(identifier);
  if (match === null) throw new Error(`unexpected identifier ${identifier}`);
  const [, name, semicolon] = match;
  if (semicolon === ';') withSemicolon.set(name, characters);
  else legacy.push([name, characters]);
}
// Every name without ";" is also in the table with ";" and the same characters: one map holds both.
for (const [name, characters] of legacy) {
  if (withSemicolon.get(name) !== characters) throw new Error(`&${name} differs from &${name};`);
}
const names = [...withSemicolon.keys()].sort();
const legacyNames = legacy.map(([name]) => name).sort();
const longest = Math.max(...legacyNames.map((name) => name.length));

const text = `// Written by scripts/html-named-references.mjs from the table of the named character references
// of the WHATWG HTML Living Standard (https://html.spec.whatwg.org/entities.json, SHA-256
// ${createHash('sha256').update(bytes).digest('hex')}; ${Object.keys(table).length} identifiers). Do not edit.
//
// Copyright © WHATWG (Apple, Google, Mozilla, Microsoft).
//
// The standard is licensed under the Creative Commons Attribution 4.0 International License; to
// the extent portions of it are incorporated into source code, such portions are licensed under
// the BSD 3-Clause License (THIRD_PARTY_NOTICES, section 12 (WHATWG)):
//
// Redistribution and use in source and binary forms, with or without modification, are permitted
// provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of
//    conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of
//    conditions and the following disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to
//    endorse or promote products derived from this software without specific prior written
//    permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
// IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
// FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
// WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

/**
 * The characters of each named character reference, by its name without "&" and ";" (the
 * ${names.length} names that the table lists with ";").
 */
export const NAMED_REFERENCES: ReadonlyMap<string, string> = new Map(
  Object.entries({
${names.map((name) => `    ${name}: ${literal(withSemicolon.get(name))},`).join('\n')}
  }),
);

/**
 * The ${legacyNames.length} names that the table also lists without ";": in text a browser decodes them
 * where no ";" follows (a name of NAMED_REFERENCES, with the same characters).
 */
export const LEGACY_NAMES: ReadonlySet<string> = new Set(
  ${literal(legacyNames.join(' '))}.split(' '),
);

/** The length of the longest name of LEGACY_NAMES. */
export const LEGACY_NAME_LENGTH = ${longest};
`;
writeFileSync(OUT, text);
console.log(`${OUT}: ${names.length} names, ${legacyNames.length} without ";"`);
