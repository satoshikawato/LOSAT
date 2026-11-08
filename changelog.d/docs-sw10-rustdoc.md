- Docs: the crate's rustdoc builds without a warning. Of 47 warnings (broken intra-doc links,
  unclosed HTML tags, links to private items in doc comments), 43 are fixed and 4 are allowed on
  the two items whose doc comment clap prints as help, so the help text is unchanged. The
  `rust lint` check runs `cargo doc --no-deps --all-features` with `-D warnings`. No engine change.
