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
   network. A fetch the user asks for (an accession, an Example set, the WindowMasker
   counts of a taxonomy ID; lines WF and WW, `docs/web/fetch_design.md`) sends only
   that identifier or the file's URL, to NCBI, EBI or CNCB (with "Fastest", the first fetch of a kind also asks each of them for one
   fixed public record, §10), and never on load or on Run. Serve every runtime asset from the site itself; do
   not load scripts from other origins in documents that handle research data. The
   response headers, including the CSP (`connect-src` names the fetch hosts), are
   defined once in `web/app/public/_headers`.
5. Keep the layers of `web/app/src` (plan §3.3): `domain` depends on nothing but the
   package's domain (`losat/domain`); `application` depends on `domain` and `ports` (the
   package's `losat/ports` included); `infra` implements `ports`; `ui` acts through
   `application` and may import only types and constants from `domain` and `ports`. The
   package's infra (`losat`) is imported only by `infra` and the composition root, and
   only `src/main.ts` and `src/composition.ts` may import every layer. The npm package
   `losat` (`web/packages/losat`: the host of the reactors, the Engine worker, the thread
   host and the WASI shim; plan LIB-4) keeps the same layers in its `src` (`domain`
   depends on nothing, `ports` on `domain`, `infra` implements `ports`) and imports
   nothing of the application or Vue. ESLint enforces this (the configurations of both),
   and `web/app/tests/unit/layers.test.ts` tests the rules themselves.
6. User-facing text, help and errors are in English.
7. The FakeEngine is for development and tests only. Its outputs must say they are not
   search results, and the UI must show a banner while it is in use.
8. Pin exact dependency versions and commit `web/package-lock.json`, the one lock file
   of the npm workspaces `web/app` and `web/packages/losat` (`web/package.json`).
   Justify any new runtime dependency in the commit message.
9. Keep application commits separate from engine (`LOSAT/`) commits.
10. The adapter (`web/adapter`) names the engine through `LOSAT::api` and the hidden
    façade `LOSAT::web_adapter` (`LOSAT/src/web_adapter.rs`) only. An engine item that the
    adapter needs is named there on its own (no glob, no module), and the façade keeps only
    items that the adapter uses; `web/adapter/tests/engine_boundary.rs` checks this. A
    function there must also be sound for every input a safe caller can give (review).

## Commands

```bash
cd web/app
npm ci          # installs the whole workspace into web/node_modules
npm run check   # the package's typecheck and lint, then the app's typecheck, lint, unit tests, production build
npm run e2e     # Playwright against `vite preview` (production headers)
```

Run both before committing changes under `web/app` or `web/packages/losat`. The package's
build (`dist/`, what `npm pack` and `npm publish` take) needs the reactors of
`web/adapter/tools/build_reactors.py`:

```bash
cd web/packages/losat
LOSAT_WEB_REACTORS=<output dir> npm run build
npm test        # the packed files, the manifest, the ABI table, the README's example
```

Publishing the package is the maintainer's act.
