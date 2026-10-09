// The default outfmt 6 row (BLAST+ "std" fields). The results screen shows the text of the
// fields as the engine wrote them (plan §4.4); it never parses them into numbers to
// show, sort or filter (those use the engine values of the HSP record).

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
