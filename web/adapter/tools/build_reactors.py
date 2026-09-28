#!/usr/bin/env python3
"""Build the two LOSAT Web reactors (docs/web/abi_v2.md §2) and record their identity.

- losat-web-serial.wasm   wasm32-wasip1          no engine features
- losat-web-threads.wasm  wasm32-wasip1-threads  LOSAT `parallel` + `wasm-threads`

Each artifact gets a JSON record with the fields that LOSAT/tests/build_wasi_artifacts.py
records for the engine's own reactors, and artifacts.json adds the build identity check
(tools/check_build_identity.py, plan TD-6). Scripts that call cargo need
RUSTUP_TOOLCHAIN=1.92.0 in environments whose default toolchain differs.

Usage: build_reactors.py --target-dir DIR --output-dir DIR [--node node]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

ADAPTER = Path(__file__).resolve().parents[1]
ROOT = ADAPTER.parents[1]
ENGINE = ROOT / "LOSAT"
SPECS = {
    "serial": ("wasm32-wasip1", ["--lib"], [], "serial-reactor"),
    "threads": ("wasm32-wasip1-threads", ["--lib", "--features", "threads"], ["parallel", "wasm-threads"], "threaded-reactor"),
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--target-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--node", default="node")
    args = parser.parse_args()
    out, target = args.output_dir.resolve(), args.target_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    identity_check = subprocess.run([sys.executable, str(ADAPTER / "tools/check_build_identity.py")],
                                    capture_output=True, text=True)
    identity = json.loads(identity_check.stdout)
    if identity_check.returncode != 0:
        print(identity_check.stdout, end="")
        raise SystemExit(f"build identity check failed: {identity['failures']}")

    config = tomllib.loads((ADAPTER / ".cargo/config.toml").read_text())
    rust = subprocess.check_output(["rustc", "-vV"], text=True)
    sysroot = Path(subprocess.check_output(["rustc", "--print", "sysroot"], text=True).strip())
    node = json.loads(subprocess.check_output([args.node, "-p", "JSON.stringify(process.versions)"], text=True))
    records = {}
    for name, (triple, options, features, kind) in SPECS.items():
        argv = ["cargo", "build", "--release", "--locked", "--target", triple, *options,
                "--target-dir", str(target / name)]
        with (out / f"{name}.build.log").open("w") as log:
            subprocess.run(argv, cwd=ADAPTER, stdout=log, stderr=subprocess.STDOUT, check=True)
        built = target / name / triple / "release/losat_web_adapter.wasm"
        artifact = json.loads(subprocess.check_output(
            [args.node, str(ENGINE / "tests/wasi_artifact.js"), str(built), kind], text=True))
        destination = out / f"losat-web-{name}.wasm"
        shutil.copyfile(built, destination)
        artifact.update(
            artifact=destination.name, target=triple, engine_features=features, rust=rust, node=node,
            argv=argv, cwd=str(ADAPTER), target_rustflags=config["target"][triple]["rustflags"],
            reactor_link_args=["<rust-sysroot>/lib/rustlib/<target>/lib/self-contained/crt1-reactor.o",
                               "--entry=_initialize"],
            environment={k: v for k, v in os.environ.items()
                         if k.startswith(("RUSTFLAGS", "CARGO_ENCODED_RUSTFLAGS", "CARGO_BUILD_RUSTFLAGS"))},
            reactor_crt_sha256=digest(sysroot / "lib/rustlib" / triple / "lib/self-contained/crt1-reactor.o"),
            cargo_lock_sha256=digest(ADAPTER / "Cargo.lock"),
            engine_cargo_lock_sha256=digest(ENGINE / "Cargo.lock"),
            build_rs_sha256=digest(ENGINE / "build.rs"),
        )
        (out / f"losat-web-{name}.json").write_text(json.dumps(artifact, indent=2) + "\n")
        records[name] = artifact
        print(f"{destination.name}: {artifact['sha256']}", flush=True)
    (out / "artifacts.json").write_text(json.dumps({"artifacts": records, "build_identity": identity},
                                                   indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
