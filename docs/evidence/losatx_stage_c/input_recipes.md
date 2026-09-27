# Additional C inputs

All input paths and hashes used in acceptance are in inputs.sha256 and accepted_v2/results.json. A generated inputs remain unchanged and retain A's recipes.

- compressed_overflow.fna: header compressed_overflow, 6468 repetitions of TGG (19404nt), newline.
- compressed_overflow.faa: header protein, 8 W residues, newline. It activates long lookup chains and the SeqSrc minimum-length10 boundary.
- split_frame{1,2,3,-1,-2,-3}.fna: length19411, N background. A S06 seed coding block is placed across the overlap around9600nt, with frame-specific C padding/reverse complement. Exact recipe payload/offsets are recorded in input_recipes.json, extracted from these deterministic frozen synthetic inputs (not search expectations).
- skip_front/back/middle.fna: length29115, N background; the same seed coding block at24000,2400,both respectively. This yields a statistically invalid first/last/middle chunk. Metadata is retained and no subject search runs for invalid chunks.

input_recipes.json stores each FASTA header, sequence length and every non-N run with its zero-based offset, so these small synthetic files can be reconstructed byte-for-byte without biological input interpretation. These are fixture recipes only; expected search states still come exclusively from NCBI.
