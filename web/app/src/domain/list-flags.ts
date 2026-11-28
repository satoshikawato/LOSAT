// The options of the tools' command lines whose value is a list that the command line splits (line WW, WW4b2 R3): one
// place for the forms (`listNameRefusal` in tool-argv.ts) and the Supported mark (`toolOptionKey` in verification.ts).
// No imports: the build's generator of the verification table (build/verification.ts) reads it too.

/**
 * makeblastdb's options whose value the command line splits at commas (`NStr::Split` of the value with the spaces at
 * its ends cut: NCBI makeblastdb.cpp `x_ProcessMaskData`; LOSAT writedb `mask_data/process.rs` `split_commas`).
 */
export const COMMA_LIST_FLAGS: ReadonlySet<string> = new Set(['-mask_data', '-mask_id', '-mask_desc', '-gi_mask_name']);

/** The number of entries of a comma list as the command line splits it: none for an empty value, else one more than its commas. */
export function commaListLength(value: string): number {
  const cut = value.trim();
  return cut === '' ? 0 : cut.split(',').length;
}
