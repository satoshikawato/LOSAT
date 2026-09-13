#!/usr/bin/env python3
"""Build and inspect separate WASI commands/reactors; no oracle dependency."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tomllib

ROOT = Path(__file__).resolve().parents[2]
CRATE = ROOT / "LOSAT"
SPECS = {
    "serial-command": ("wasm32-wasip1", ["--bin", "LOSAT", "--no-default-features"], []),
    "threaded-command": ("wasm32-wasip1-threads", ["--bin", "LOSAT", "--features", "wasm-threads"], ["parallel", "wasm-threads"]),
    "threaded-reactor": ("wasm32-wasip1-threads", ["--lib", "--features", "wasm-threads"], ["parallel", "wasm-threads"]),
    "serial-reactor": ("wasm32-wasip1", ["--lib", "--no-default-features"], []),
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


# NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
# TBlastThreads the_threads(GetNumberOfThreads()); (*thread)->Run(); (*thread)->Join(&result);
# Artifact identity belongs to the host ABI, not to any NCBI implementation.
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--target-dir", type=Path, required=True)
    parser.add_argument("--node", default="node")
    parser.add_argument("--reverse-order", action="store_true")
    args = parser.parse_args()
    out, target = args.output_dir.resolve(), args.target_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    config = tomllib.loads((CRATE / ".cargo/config.toml").read_text())
    rust = subprocess.check_output(["rustc", "-vV"], text=True)
    sysroot = Path(subprocess.check_output(["rustc", "--print", "sysroot"], text=True).strip())
    node = subprocess.check_output([args.node, "-p", "JSON.stringify(process.versions)"], text=True)
    recorded, identities = {}, {}
    for kind in reversed(SPECS) if args.reverse_order else SPECS:
        triple, options, features = SPECS[kind]
        argv = ["cargo", "build", "--release", "--locked", "--target", triple, *options, "--target-dir", str(target / kind)]
        with (out / f"{kind}.build.log").open("w") as log:
            subprocess.run(argv, cwd=CRATE, stdout=log, stderr=subprocess.STDOUT, check=True)
        for previous, expected in recorded.items():
            if digest(previous) != expected:
                raise RuntimeError(f"build overwrote another artifact: {previous}")
        built = target / kind / triple / "release/LOSAT.wasm"
        identity = json.loads(subprocess.check_output([args.node, str(CRATE / "tests/wasi_artifact.js"), str(built), kind], text=True))
        destination = out / f"losat-{kind}.wasm"
        shutil.copyfile(built, destination)
        recorded[built] = digest(built)
        identity.update(target=triple, features=features, rust=rust, node=json.loads(node),
                        argv=argv, cwd=str(CRATE), target_rustflags=config["target"][triple]["rustflags"],
                        reactor_link_args=["<rust-sysroot>/lib/rustlib/<target>/lib/self-contained/crt1-reactor.o", "--entry=_initialize"] if kind.endswith("reactor") else [],
                        environment={k: v for k, v in os.environ.items() if k.startswith(("RUSTFLAGS", "CARGO_ENCODED_RUSTFLAGS", "CARGO_BUILD_RUSTFLAGS"))},
                        reactor_crt_sha256=digest(sysroot / "lib/rustlib" / triple / "lib/self-contained/crt1-reactor.o") if kind.endswith("reactor") else None,
                        cargo_lock_sha256=digest(CRATE / "Cargo.lock"), build_rs_sha256=digest(CRATE / "build.rs"))
        (out / f"losat-{kind}.json").write_text(json.dumps(identity, indent=2) + "\n")
        identities[kind] = identity
        print(f"{kind}: {identity['sha256']}", flush=True)
    runtime_files = {}
    for name in ["run_losat_wasi.js", "run_losat_wasi_threads.js", "wasi_thread_host.js", "wasi_artifact.js", "wasi_shared_memory.js"]:
        shutil.copyfile(CRATE / "tests" / name, out / name)
        runtime_files[name] = digest(out / name)
    (out / "runtime-files.json").write_text(json.dumps(runtime_files, indent=2) + "\n")
    (out / "artifacts.json").write_text(json.dumps(identities, indent=2) + "\n")


if __name__ == "__main__":
    main()
