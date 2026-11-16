// The default outfmt 6 row (BLAST+ "std" fields). The results screen shows the text of the
// fields as the engine wrote them (plan §4.4). It sorts with the engine values of the HSP
// record; the E value and bit score filters read the numbers that the rows show (shownScores),
// so that a row passes or not as it is shown (WA-RS-01). Nothing here formats a value.
import type { HspTable } from './hsp-table';

/** The twelve default fields, in their order (BLAST+ `-help`: "std"). */
export const OUTFMT6_FIELDS = Object.freeze([
  'qseqid',
  'sseqid',
  'pident',
  'length',
  'mismatch',
  'gapopen',
  'qstart',
  'qend',
  'sstart',
  'send',
  'evalue',
  'bitscore',
] as const);

export type Outfmt6Field = (typeof OUTFMT6_FIELDS)[number];
export type Outfmt6Row = Readonly<Record<Outfmt6Field, string>>;

/**
 * Splits one row (the bytes of an HSP's `out6` range, decoded) into its fields. Throws if
 * the row does not have the twelve fields, so a wrong range is never shown as values.
 */
export function splitOutfmt6Row(text: string): Outfmt6Row {
  const fields = text.replace(/\r?\n$/, '').split('\t');
  if (fields.length !== OUTFMT6_FIELDS.length) {
    throw new Error(`an outfmt 6 row has ${OUTFMT6_FIELDS.length} fields, this one ${fields.length}: ${JSON.stringify(text)}`);
  }
  return Object.freeze(Object.fromEntries(OUTFMT6_FIELDS.map((name, i) => [name, fields[i]!])) as Record<Outfmt6Field, string>);
}

const TAB = 0x09;
const EVALUE_FIELD = OUTFMT6_FIELDS.indexOf('evalue');
const BITSCORE_FIELD = OUTFMT6_FIELDS.indexOf('bitscore');

/**
 * The numbers that the outfmt 6 row of each HSP shows as its E value and bit score, as the
 * engine wrote them: NCBI writes an E value below 1e-180 as 0.0 and a bit score of 100 or more
 * without decimals, so these are not the record's values. The view filters compare these, so
 * that a row shown as 0.0 passes "E value ≤ 0". A record without an outfmt 6 row, or a field that
 * is not a number, keeps the record's value. Built once per run, without an object per row.
 */
export function shownScores(table: HspTable, out6: Uint8Array): Pick<HspTable, 'eValue' | 'bitScore'> {
  const eValue = Float64Array.from(table.eValue);
  const bitScore = Float64Array.from(table.bitScore);
  for (let row = 0; row < table.count; row++) {
    const start = table.out6Start[row]!;
    const end = Math.min(table.out6End[row]!, out6.length);
    if (start < 0) continue;
    let field = 0;
    let from = start;
    for (let at = start; at <= end && field <= BITSCORE_FIELD; at++) {
      const byte = at < end ? out6[at]! : TAB;
      if (byte !== TAB && byte !== 0x0a && byte !== 0x0d) continue;
      if (field === EVALUE_FIELD || field === BITSCORE_FIELD) {
        const text = String.fromCharCode(...out6.subarray(from, at)).trim();
        const value = text === '' ? NaN : Number(text);
        if (!Number.isNaN(value)) (field === EVALUE_FIELD ? eValue : bitScore)[row] = value;
      }
      field++;
      from = at + 1;
    }
  }
  return { eValue, bitScore };
}
