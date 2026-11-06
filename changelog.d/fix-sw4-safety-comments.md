- Internal: `// SAFETY:` notes on 43 `unsafe` blocks that had none in the 5 lines above (SWEBOK 2026-11-04 audit,
  KA4). Each says what the operation needs and what guarantees it. No code change; the release binary's
  machine code is unchanged. Blocks whose guarantee rests on what callers of a public safe function pass are left
  without a note and listed for a later fix.
