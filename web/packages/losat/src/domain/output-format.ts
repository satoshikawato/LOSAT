// The outputs that every run of the engine writes (ABI v2 streams 0, 6 and 7), and the output of the
// argv's own -outfmt (stream 4: any other format that the CLI writes, and 6 or 7 with fields;
// docs/web/abi_v2.md §5, §7).
export type OutputFormat = 0 | 6 | 7;
export const OUTPUT_FORMATS: readonly OutputFormat[] = Object.freeze([0, 6, 7]);

/** The output of the argv's own -outfmt (ABI stream 4). */
export const CUSTOM_OUTPUT = 4;

/** One output of a run: a format that every run writes, or CUSTOM_OUTPUT. */
export type RunOutput = OutputFormat | typeof CUSTOM_OUTPUT;
