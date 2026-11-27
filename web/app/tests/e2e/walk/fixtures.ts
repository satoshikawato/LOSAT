// The inputs of the random walk: FASTA texts and files made from a seed, so that every search has
// hits in every mode. A protein P is drawn and written back as a coding sequence N (no stop
// codons), and the records are mutated copies and parts of P and N, a reverse complement, and
// unrelated sequences. Some deflines, a job title and file names carry HTML and script, which the
// page must show as text (walk/steps.ts counts `window.__walkXss`, which they would raise if run).
// A pure module (tests/unit/walk.test.ts).
import { createRandom, type Random } from './plan';

export type Alphabet = 'nucleotide' | 'protein';

export interface FixtureText {
  readonly id: string;
  readonly alphabet: Alphabet;
  readonly text: string;
}

export interface FixtureFile {
  readonly name: string;
  readonly alphabet: Alphabet;
  readonly text: string;
}

/** What a step pastes or opens: texts and files of each alphabet, and texts that are not plain FASTA. */
export interface FixtureLibrary {
  readonly texts: readonly FixtureText[];
  /** Texts that are not FASTA as NCBI reads it: no defline, CRLF line ends, a byte order mark. */
  readonly odd: readonly FixtureText[];
  readonly files: readonly FixtureFile[];
}

/** The script that a defline, a title or a file name carries: it counts on window.__walkXss if it ever runs. */
export const XSS_PAYLOAD = `<img src=x onerror="window.__walkXss=(window.__walkXss||0)+1">`;
export const XSS_SVG = `<svg onload="window.__walkXss=(window.__walkXss||0)+1">`;

const AMINO = 'ACDEFGHIKLMNPQRSTVWY';
const NUCLEOTIDES = 'ACGT';
/** The standard code's codons of each amino acid (no stop codons). */
const CODONS: Readonly<Record<string, readonly string[]>> = Object.freeze({
  A: ['GCT', 'GCC', 'GCA', 'GCG'],
  C: ['TGT', 'TGC'],
  D: ['GAT', 'GAC'],
  E: ['GAA', 'GAG'],
  F: ['TTT', 'TTC'],
  G: ['GGT', 'GGC', 'GGA', 'GGG'],
  H: ['CAT', 'CAC'],
  I: ['ATT', 'ATC', 'ATA'],
  K: ['AAA', 'AAG'],
  L: ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
  M: ['ATG'],
  N: ['AAT', 'AAC'],
  P: ['CCT', 'CCC', 'CCA', 'CCG'],
  Q: ['CAA', 'CAG'],
  R: ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
  S: ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
  T: ['ACT', 'ACC', 'ACA', 'ACG'],
  V: ['GTT', 'GTC', 'GTA', 'GTG'],
  W: ['TGG'],
  Y: ['TAT', 'TAC'],
});

function draw(random: Random, letters: string, length: number): string {
  let out = '';
  for (let i = 0; i < length; i++) out += letters[random.int(letters.length)];
  return out;
}

/** A copy of `sequence` with each letter replaced by another with probability `rate`. */
function mutate(random: Random, sequence: string, letters: string, rate: number): string {
  let out = '';
  for (const letter of sequence) {
    if (random.next() < rate) {
      const others = letters.replace(letter, '');
      out += others[random.int(others.length)];
    } else out += letter;
  }
  return out;
}

function reverseComplement(sequence: string): string {
  const pair: Record<string, string> = { A: 'T', C: 'G', G: 'C', T: 'A' };
  return [...sequence].reverse().map((letter) => pair[letter] ?? 'N').join('');
}

/** A record with its sequence in lines of `width`. */
function record(defline: string, sequence: string, width = 60): string {
  const lines = [];
  for (let i = 0; i < sequence.length; i += width) lines.push(sequence.slice(i, i + width));
  return `>${defline}\n${lines.join('\n')}\n`;
}

/** The walk's inputs for a seed. */
export function makeFixtures(seed: number): FixtureLibrary {
  const random = createRandom(seed ^ 0x0f17a5);
  const protein = draw(random, AMINO, 220);
  const coding = [...protein].map((amino) => random.pick(CODONS[amino]!)).join('');
  const nt = {
    whole: mutate(random, coding, NUCLEOTIDES, 0.02),
    part: mutate(random, coding.slice(90, 540), NUCLEOTIDES, 0.03),
    reverse: reverseComplement(mutate(random, coding.slice(150, 600), NUCLEOTIDES, 0.02)),
    other: draw(random, NUCLEOTIDES, 400),
    subject: coding,
  };
  const aa = {
    whole: mutate(random, protein, AMINO, 0.05),
    part: mutate(random, protein.slice(20, 170), AMINO, 0.08),
    other: draw(random, AMINO, 160),
    subject: protein,
  };
  const texts: FixtureText[] = [
    { id: 'nt-one', alphabet: 'nucleotide', text: record('nq1 coding sequence', nt.whole) },
    {
      id: 'nt-three',
      alphabet: 'nucleotide',
      text:
        record(`nq2 part ${XSS_PAYLOAD}`, nt.part) +
        record('nq3 reverse "strand" & <b>bold</b>', nt.reverse, 70) +
        record('nq4 unrelated 配列 é', nt.other),
    },
    { id: 'nt-subject', alphabet: 'nucleotide', text: record('ns1 subject', nt.subject) + record('ns2 lower', nt.part.toLowerCase()) },
    { id: 'aa-one', alphabet: 'protein', text: record('pq1 protein', aa.whole) },
    {
      id: 'aa-three',
      alphabet: 'protein',
      text: record(`pq2 part ${XSS_SVG}`, aa.part) + record('pq3 unrelated <i>x</i>', aa.other) + record('pq4 whole', aa.subject, 80),
    },
  ];
  const odd: FixtureText[] = [
    { id: 'nt-bare', alphabet: 'nucleotide', text: `${nt.whole.slice(0, 300)}\n` },
    { id: 'nt-crlf', alphabet: 'nucleotide', text: record('nq5 crlf', nt.part).replaceAll('\n', '\r\n') },
    { id: 'aa-bom', alphabet: 'protein', text: `﻿${record('pq5 bom', aa.part)}` },
    { id: 'aa-no-newline', alphabet: 'protein', text: record('pq6 no final newline', aa.whole).trimEnd() },
  ];
  const files: FixtureFile[] = [
    { name: 'query.fasta', alphabet: 'nucleotide', text: record('fq1 file query', nt.whole) },
    { name: 'two genes.fna', alphabet: 'nucleotide', text: record('fq2 gene a', nt.part) + record('fq3 gene b', nt.reverse) },
    { name: `${XSS_SVG}.fa`, alphabet: 'nucleotide', text: record(`fs1 ${XSS_PAYLOAD}`, nt.subject) },
    { name: 'proteins.faa', alphabet: 'protein', text: record('fp1 file protein', aa.whole) + record('fp2 file part', aa.part) },
    { name: "it's.faa", alphabet: 'protein', text: record('fp3 quote', aa.subject) },
  ];
  return { texts, odd, files };
}

/** The alphabet of a role's input in a mode. */
export function alphabetOf(program: string, role: 'query' | 'subject'): Alphabet {
  const query: Record<string, Alphabet> = { blastn: 'nucleotide', blastp: 'protein', psiblast: 'protein', blastx: 'nucleotide', tblastn: 'protein', tblastx: 'nucleotide' };
  const subject: Record<string, Alphabet> = { blastn: 'nucleotide', blastp: 'protein', psiblast: 'protein', blastx: 'protein', tblastn: 'nucleotide', tblastx: 'nucleotide' };
  const alphabet = (role === 'query' ? query : subject)[program];
  if (alphabet === undefined) throw new Error(`no mode ${program}`);
  return alphabet;
}

/**
 * The item that a draw in [0, 1) picks: mostly one of the mode's alphabet (`fitting`), sometimes
 * one of the other alphabet, sometimes an odd one.
 */
export function pickFor<T extends { alphabet: Alphabet }>(items: readonly T[], odd: readonly T[], alphabet: Alphabet, kind: number, which: number): T {
  const fitting = items.filter((item) => item.alphabet === alphabet);
  const other = items.filter((item) => item.alphabet !== alphabet);
  const pool = kind < 0.8 || odd.length === 0 ? (kind < 0.7 || other.length === 0 ? fitting : other) : odd;
  return pool[Math.min(pool.length - 1, Math.floor(which * pool.length))]!;
}
