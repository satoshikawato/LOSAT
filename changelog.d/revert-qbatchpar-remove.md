- `blastn`'s query-batch slots are removed: the switch `LOSAT_X_QBATCHPAR`, its tuning variables
  `LOSAT_X_QBATCHPARMEM`, `LOSAT_X_QBATCHPARSLOTS` and `LOSAT_X_QBATCHPARMISS`, and the switches that
  acted only on the slots, `LOSAT_X_DCTWOSLOTS`, `LOSAT_X_BNEIGHTSLOTS` and `LOSAT_X_BNQSLOTS`
  (Owner decision 2026-10-28: the slots' max RSS rise exceeded their time saving). A run ignores
  these variables, and the batch loop searches every query batch on the whole pool, as it did with
  the slots off. The `[X_STATS] QBATCHPAR` line is gone, the `[X_STATS] QSTREAM` line no longer has
  `slots=`, `slot_copies=`, `slot_records=` and `stale=`, and `LOSAT_X_QSTREAM`'s rule no longer has
  `off=slots`. Reports do not change.
