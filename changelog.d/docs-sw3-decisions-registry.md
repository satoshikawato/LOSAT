- `docs/decisions/`: the line and stage IDs the repository cites (`LINE_IDS.md`, 151 IDs), a copy of the
  Owner decisions it cites by date or by decision-log entry (`OWNER_DECISIONS.md`), and the approval
  state of the 20 conditional switches added after the 2026-10-27 stop rule (`SWITCH_APPROVALS.md`).
  `scripts/line_ids.py list` lists the cited IDs; `check` names those cited in `CHANGELOG.md`,
  `changelog.d/` and `docs/` without a row (CI `rust lint` runs it; warnings only).
