# Product Decision: Bounded native-NCBI platform variance

- Decision ID: `PD-NCBI-PLATFORM-VARIANCE`
- Version: 1.1
- Date: 2026-09-02 (1.1: 2026-10-12, see Amendments)
- Status: Accepted

## Scope

This decision governs the six native-NCBI representative checks in the
Windows x64, macOS arm64, and macOS x64 release-certification jobs. It replaces
the former requirement that one LOSAT output simultaneously equal the frozen
Linux canonical bytes and every platform-local official NCBI binary output.

It does not change LOSAT production code, the 43-case LOSAT canonical manifest,
the existing BLASTN source-underdetermined class, or the accepted local-subject
TBLASTX `db_gencode` deviation.

## Source and executable authorities

Official NCBI BLAST+ source remains the sole behavioral authority. In NCBI
BLAST+ 2.17.0, `blast_format.cpp:770-832` selects the output formatter and
`tabular.cpp:59-91,1098-1108` expands and writes the requested tabular fields
and row terminators. Official NCBI executables remain validation oracles only;
they are never LOSAT runtime, build, feature, fallback, FFI, or subprocess
implementation dependencies.

Characterization run `33582928402` established that official NCBI 2.17.0
binaries can emit different raw output on the three supported native
platforms. That is an NCBI-side reference property. LOSAT is not authorized to
be platform-variable.

## Dual release contract

Release certification has two distinct gates.

### Gate A — `LOSAT_CANONICAL`

Each of the 43 LOSAT cases on each certified platform must satisfy:

```text
SHA256(raw LOSAT output bytes) == frozen PR 5 canonical SHA256
```

Gate A permits no normalization, sorting, field tolerance, platform baseline,
or platform key in production behavior. A native-NCBI fingerprint can never
compensate for a Gate A failure and never becomes expected LOSAT output.

#### Known frozen mismatches (version 1.1)

`LOSAT/tests/frozen_mismatch_allowlist.json` lists the Gate A cases whose
approved LOSAT output is no longer the frozen PR 5 bytes. For a listed case
the certifier requires, on each certified platform:

```text
SHA256(raw LOSAT output bytes) == allowed_actual_sha256 of the entry
```

and records the step as `ALLOWED_KNOWN_MISMATCH`, never as exact PR 5 bytes;
the platform summary then reports `all_exact_pr5_raw_bytes: false` and names
the case. Any other output of a listed case fails, including the frozen PR 5
bytes (a stale entry, which is to be removed). Every case that is not listed
keeps the frozen-hash rule above.

An entry is one byte string for all platforms. It is not a platform baseline,
a tolerance or a normalization, and the frozen hashes in
`LOSAT/tests/platform_native_v010_canonical.tsv` are not rewritten. The
certifier reads the allowlist from the candidate commit, requires each Gate A
entry to declare exactly the matrix and repeatability steps of its case
(`frozen_gate_a_steps`), and binds the file's SHA-256 into platform identity,
resume and the cross-platform aggregate.

Adding or changing an entry needs its own approved record of the new output.
The one entry, `blastn Sakai.MG1655.megablast` (Gate A steps 11, 52 and 53),
was approved on 2026-09-28
(`docs/evidence/losatx_stage_g_authority_v3/run-20260928T121516Z/`): the
frozen v0.1.0 bytes differ from NCBI BLAST+ 2.17.0 in 5 of 6482 output lines
(an equal-score greedy traceback), and the current LOSAT output is
byte-identical to NCBI for the same command.

### Gate B — `PLATFORM_NATIVE_NCBI_REFERENCE`

Each of the six existing native-NCBI representatives must resolve one exact
record by:

```text
authority version + platform ID + program + case ID
```

The selected record binds the NCBI release, official archive filename and
published checksum, executable hash and architecture, query and subject byte
identities, ordered argv, output format, controlled lexical paths, repository
working-directory contract, `shell=false`, empty certifier-controlled native
environment overrides, and the exact raw output fingerprint.

Gate B passes only when raw SHA-256, byte length, newline profile, and data-row
count all match. Raw SHA-256 is decisive. An unknown selector or fingerprint is
`NATIVE_NCBI_AUTHORITY_MISS` and hard-fails; there is no Linux, architecture,
OS-family, nearest-platform, or semantic fallback and no automatic whitelist
learning.

Windows CRLF bytes are normative fingerprint data. Platform-local TBLASTX row
order is normative fingerprint data. Neither may be normalized or sorted for
Gate B acceptance.

## Orthogonal classifications

The existing LOSAT parity class and the descriptive native-NCBI reference
class are separate axes. Native-versus-LOSAT structured analysis is diagnostic
and cannot rescue a Gate B raw mismatch.

`Sakai.MG1655.megablast` remains `SOURCE_UNDETERMINED_ACCEPTED` only on the
LOSAT parity axis. The TBLASTX d06 case remains
`approved_db_gencode_deviation / HSP_SET_DIFF`; the accepted TBLASTX deviation
ceiling remains exactly six. Native-NCBI platform variance creates no seventh
deviation and no new LOSAT parity exception.

## Evidence authority

The immutable authority is
`LOSAT/tests/ncbi_platform_variance_v010.json`. Its normative layer alone can
decide Gate B. Diagnostic metadata and characterization provenance never grant
acceptance.

The authority binds:

- characterization run `33582928402` at
  `be67253156ce10b852166709e9519ab54709fd80`;
- corrected analyzer commit
  `cfd213f57691607429c4e21edbf314899310d805`;
- corrected replay SHA-256
  `37e95cf9252ede022cc0d3f3f6a43f4776f78b3e31febd2a76bf2d2d667be2db`;
- PR 6 V3 `0e29e2201e2d1b03124c0b9d6698a81bfed8cec0` and run
  `33503928773`;
- the retained PR 5/Linux evidence lineage.

The R1 rich diagnostic evidence contains 11 exhaustive representation groups
(7 authoritative-seeded and 4 rich-only). Those hashes are frozen supporting
provenance with `runtime_enforced=false`. Rich fields are not searched during
normal V4 certification and are never an alternate acceptance path.

At runtime the certifier hashes the exact authority-file bytes. Authority
version and file SHA-256 bind platform identity, completion, resume, per-platform
summary, and cross-platform aggregate evidence. Evidence from another version
or file hash cannot resume or aggregate.

## Workflow claim

Each platform keeps the exact `43 + 6 + 12 = 61` execution plan. A matrix cell
may certify only its platform. The dependent hosted aggregator consumes exactly
the three expected artifacts, verifies their complete evidence manifests and
subgates, and alone may emit `CROSS_PLATFORM_NATIVE_CERTIFIED`.

## Change policy

An NCBI upgrade or newly observed official fingerprint requires a fresh
three-platform characterization, a new immutable authority version, review,
and explicit governance acceptance. Historical authority files are not edited
or auto-learned in place.

## Amendments

- 1.1 (2026-10-12, Owner): Gate A applies
  `LOSAT/tests/frozen_mismatch_allowlist.json` ("Known frozen mismatches"
  above), as the pull-request and Wasm regression checks already did. Under
  1.0 the certifier stopped at step 11 of 61 for every candidate that carries
  the approved Sakai output (run `38170565780` on `main`), so steps 12 to 61
  did not run. The Sakai and TBLASTX d06 ratchets, Gate B and the authority
  file are unchanged.

## Non-goals

- No LOSAT production, algorithm, comparator, traceback, scoring, pruning,
  filtering, ordering, statistics, or formatter change.
- No platform-specific LOSAT output path.
- No broad semantic tolerance or inferred NCBI equivalence class.
- No paired rich search in release certification.
- No expansion beyond the three platforms and six representatives recorded in
  the authority.
