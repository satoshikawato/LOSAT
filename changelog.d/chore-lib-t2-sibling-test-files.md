- The unit tests of `LOSAT/src` live in sibling `<stem>_tests.rs` files (`#[cfg(test)] #[path = ".."] mod tests;`),
  which the crate's `exclude` leaves out of the published `losat` crate: 6.28 -> 5.02 MiB, 62.8 % -> 50.2 % of the
  10 MiB that crates.io takes (dev of 2026-11-21). Test names, module paths and product code are unchanged.
  `python3 LOSAT/tests/check_inline_test_modules.py` (CI) fails on an inline test module and, with `--fix`, moves it.
