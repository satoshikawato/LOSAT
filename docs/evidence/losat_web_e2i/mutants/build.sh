#!/usr/bin/env bash
# Build three mutants of LOSAT (SD fixture discrimination check), one after the other.
set -u
M=$HOME/.cache/losat-web-gui-target/sd-mutants
cd $M/wt/LOSAT
export RUSTUP_TOOLCHAIN=1.92.0
B=src/algorithm/blastn
build() {
  nice -n 10 cargo build --release --locked --target-dir $M/target > $M/build-$1.log 2>&1 && cp $M/target/release/LOSAT $M/LOSAT-$1 && echo "$1 built"
  git checkout -q -- src
}
git checkout -q -- src
# M1: the second template's lookup is skipped.
python3 - <<'PY'
from pathlib import Path
p = Path("src/algorithm/blastn/blast_engine/run.rs"); s = p.read_text()
old = "if let Some(second_template_index) = second_template_index {"
assert s.count(old) == 1; p.write_text(s.replace(old, "if let Some(second_template_index) = second_template_index.filter(|_| false) {"))
PY
build M1
# M2: dc-megablast's two-hit window 0 (one hit).
python3 - <<'PY'
from pathlib import Path
p = Path("src/algorithm/blastn/coordination.rs"); s = p.read_text()
old = "            window_size: 40,\n"
assert s.count(old) == 1; p.write_text(s.replace(old, "            window_size: 0,\n"))
PY
build M2
# M3: the 11/21 coding template scanned with the 11/18 scanner.
python3 - <<'PY'
from pathlib import Path
p = Path("src/algorithm/blastn/disc_lookup.rs"); s = p.read_text()
old = "    } else if template_type == DiscTemplateType::T11_21Coding {\n        DiscScanSubjectKind::Scan11_21\n"
assert s.count(old) == 1; p.write_text(s.replace(old, "    } else if template_type == DiscTemplateType::T11_21Coding {\n        DiscScanSubjectKind::Scan11_18\n"))
PY
build M3
# M4 (built afterwards, the same way): the 11/18 and 11/21 coding templates scanned by the general
# scanner, an equivalent mutant (NCBI's dedicated scanners return the same words in the same order).
python3 - <<'PY'
from pathlib import Path
p = Path("src/algorithm/blastn/disc_lookup.rs"); s = p.read_text()
for kind in ("Scan11_18", "Scan11_21"):
    old = f"        DiscScanSubjectKind::{kind}\n    }} else"
    assert s.count(old) == 1; s = s.replace(old, "        DiscScanSubjectKind::General\n    } else")
p.write_text(s)
PY
build M4
echo done
