- Docs: the crate's rustdoc builds without a warning. Of 47 warnings (broken intra-doc links,
  unclosed HTML tags, links to private items in doc comments), all are fixed. The
  `rust lint` check runs `cargo doc --no-deps --all-features` with `-D warnings`. No engine change.
