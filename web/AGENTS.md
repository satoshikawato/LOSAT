# AGENTS.md — web/

Instructions for coding agents working under `web/`. The root [`AGENTS.md`](../AGENTS.md)
still governs everything under `LOSAT/`, including engine changes made for LOSAT Web.

## Authority

- [`PD-LOSAT-WEB-APP-BOUNDARY`](../docs/product_decisions/PD-LOSAT-WEB-APP-BOUNDARY.md)
  defines what application code may and may not do.
- [`docs/losat_web_gui_plan.md`](../docs/losat_web_gui_plan.md) is the implementation
  plan; [`docs/losat_web_gui_sessions/`](../docs/losat_web_gui_sessions/README.md) holds
  the session instructions.
- [`docs/web/losat_web_design_v0.1.md`](../docs/web/losat_web_design_v0.1.md) holds the
  application requirements, as amended by the plan's decisions (plan §0.4).
- [`docs/web/abi_v2.md`](../docs/web/abi_v2.md) is the contract between the engine
  adapter and the application.

## Rules

1. Do not compute or format BLAST-defined values (scores, E-values, identities,
   coverages, NCBI number formats, alignment text) under `web/`. Take them from engine
   output or engine records.
2. Store and export outfmt 0, 6 and 7 byte for byte. Never rewrite, filter or
   regenerate them. Filtered or selected exports use application formats (CSV, JSON).
3. Give the engine only an argv and FASTA bytes made of whole original records (the
   run snapshot). Never split one search into independent searches and join the
   results.
4. Never send research data (sequences, headers, file names, results, notes) over the
   network. Serve every runtime asset from the site itself; do not load scripts from
   other origins in documents that handle research data. The response headers,
   including the CSP, are defined once in `web/app/public/_headers`.
5. Keep the layers of `web/app/src` (plan §3.3): `domain` depends on nothing;
   `application` depends on `domain` and `ports`; `infra` implements `ports`; `ui`
   acts through `application` and may import only types and constants from `domain`
   and `ports`. Only `src/main.ts` and `src/composition.ts` may import every layer.
   ESLint enforces this, and `tests/unit/layers.test.ts` tests the rule itself.
6. User-facing text, help and errors are in English.
7. The FakeEngine is for development and tests only. Its outputs must say they are not
   search results, and the UI must show a banner while it is in use.
8. Pin exact dependency versions and commit `package-lock.json`. Justify any new
   runtime dependency in the commit message.
9. Keep application commits separate from engine (`LOSAT/`) commits.

## Commands

```bash
cd web/app
npm ci
npm run check   # typecheck, lint, unit tests, production build
npm run e2e     # Playwright against `vite preview` (production headers)
```

Run both before committing changes under `web/app`.
