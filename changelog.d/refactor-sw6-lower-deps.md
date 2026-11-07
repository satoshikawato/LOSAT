- Engine: `blastinput` no longer imports `algorithm`, `stats` no longer imports `blastinput`, and `utils` no longer imports
  `cli`: the protein query split sizes, the template type parser, the query reading type, the program entry names and the
  protein options check move to the layer that uses them. No output changes.
