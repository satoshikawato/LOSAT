// What the first bytes of an input say about it besides FASTA (WA-IN-04, WA-IN-05). The
// engine reads every input as FASTA, as the command line does (NCBI's reader): a GenBank,
// EMBL or FASTQ file is read as sequence lines without a defline, and a byte order mark
// makes the first line implausible data. The Web neither refuses nor changes such an input
// (the command line reads it the same way); it names what the input looks like, so that the
// source's card can say how to make it FASTA.

/** A text encoding mark at the start of the bytes, or another sequence file format. */
export type InputLead =
  | { readonly mark: 'UTF-8' | 'UTF-16' }
  | { readonly format: 'GenBank' | 'EMBL' | 'FASTQ' };

/** What the first bytes of an input (its preview) look like, or undefined for anything else. */
export function inputLead(bytes: Uint8Array): InputLead | undefined {
  if (bytes[0] === 0xef && bytes[1] === 0xbb && bytes[2] === 0xbf) return { mark: 'UTF-8' };
  if ((bytes[0] === 0xff && bytes[1] === 0xfe) || (bytes[0] === 0xfe && bytes[1] === 0xff)) return { mark: 'UTF-16' };
  const lines = new TextDecoder('latin1').decode(bytes).replace(/^\s+/, '').split(/\r?\n|\r/);
  const first = lines[0] ?? '';
  if (/^LOCUS\s/.test(first)) return { format: 'GenBank' };
  if (/^ID {3}\S/.test(first)) return { format: 'EMBL' };
  // A FASTQ record is four lines: `@` title, sequence, `+` (title again or nothing), qualities.
  // A long read's sequence line may run past the preview, so the third line is checked only when it is there.
  if (/^@\S/.test(first) && /^[A-Za-z*.-]*$/.test(lines[1] ?? '') && (lines.length < 3 || lines[2]!.startsWith('+'))) {
    return { format: 'FASTQ' };
  }
  return undefined;
}
