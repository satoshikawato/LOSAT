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
