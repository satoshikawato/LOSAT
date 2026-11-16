#!/bin/bash -e

cd "${SRC_DIR}/LOSAT"

cargo-bundle-licenses \
  --format yaml \
  --output "${SRC_DIR}/THIRDPARTY.yml"

cargo install \
  --no-track \
  --locked \
  --verbose \
  --root "${PREFIX}" \
  --path .

# NCBI's program names on PATH (relative links, so the prefix can move).
"${PREFIX}/bin/losat" --install-shims "${PREFIX}/bin"
