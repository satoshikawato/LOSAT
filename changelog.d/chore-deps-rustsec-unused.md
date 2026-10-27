- Dependencies: `anyhow` 1.0.102 -> 1.0.104 (RUSTSEC-2026-0190) and `crossbeam-epoch` 0.9.18 -> 0.9.21
  (RUSTSEC-2026-0204) in both lock files; the unused direct dependencies `lazy_static`, `bitvec` and
  `num_cpus` are removed (feature `parallel` is now `rayon` alone). No code or output change.
