- Docs: the crate's rustdoc builds without a warning (47 broken intra-doc links, unclosed HTML tags
  and links to private items in doc comments fixed), and the `rust lint` check runs
  `cargo doc --no-deps --all-features` with `-D warnings`. Help text is unchanged. No engine change.
