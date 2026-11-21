- `makembindex -idmap` writes the `.map` ids that NCBI writes: the first word of the title as stored
  (HTML entities such as `&amp;` and the punctuation at the end are kept, a third-party id's title
  starts `TPA:`), and for a database sequence without a title and with a tax id the scientific name
  from the taxonomy database, or `No`. `windowmasker -ids` and `-exclude_ids` matched as text use
  the same titles, so a database sequence without a title and with a tax id is no longer refused.
