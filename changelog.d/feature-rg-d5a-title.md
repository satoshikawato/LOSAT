- `makembindex -idmap` writes the `.map` ids that NCBI writes: the first word of the title as stored
  (HTML entities such as `&amp;` and a `,`, `;` or `~` at the end are kept; only the dots and spaces
  at the end are dropped; a third-party id's title starts `TPA:`), and for a database sequence
  without a title and with a tax id the first word of the organism's scientific name from the
  taxonomy database, with its first letter in upper case, or `No` without one. `windowmasker -ids`
  and `-exclude_ids` matched as text use the same titles, so a database sequence without a title
  and with a tax id is no longer refused. `makembindex` refuses a database sequence without letters
  (only a damaged database has one), where NCBI's stops with a null-pointer report.
