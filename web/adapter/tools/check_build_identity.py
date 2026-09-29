#!/usr/bin/env python3
"""Build identity of the adapter and the engine (plan TD-6).

The adapter reactors must be built under the same conditions as the certified LOSAT
command-WASI builds. This check compares, and prints as JSON:
- `[profile.release]` of web/adapter/Cargo.toml and LOSAT/Cargo.toml;
- the target rustflags of web/adapter/.cargo/config.toml and LOSAT/.cargo/config.toml;
- the version of every dependency that both Cargo.lock files contain;
- the link arguments: the adapter has no build script of its own (no `build` key and no
  web/adapter/build.rs, which Cargo would find by itself), so its cdylib gets the
  reactor start-up (`crt1-reactor.o`, `--entry=_initialize`) from LOSAT/build.rs through
  the `LOSAT` dependency, exactly as LOSAT's own reactors do;
- that no other Cargo configuration applies: no `.cargo/config[.toml]` in the parent
  directories of either crate or in CARGO_HOME, and no profile, build or rustflags
  environment variable (Cargo reads all of these in addition to the compared files).

Usage: check_build_identity.py [--out FILE.json]   (exit status 1 on any difference)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
ENGINE = ROOT / "LOSAT"
ADAPTER = ROOT / "web" / "adapter"


def lock_versions(path: Path) -> dict[str, set[str]]:
    versions: dict[str, set[str]] = {}
    for package in tomllib.loads(path.read_text())["package"]:
        if "source" in package:  # registry packages; path packages have no source
            versions.setdefault(package["name"], set()).add(package["version"])
    return versions


ENVIRONMENT = ("RUSTFLAGS", "CARGO_ENCODED_RUSTFLAGS", "CARGO_BUILD_", "CARGO_PROFILE_", "CARGO_TARGET_")


def other_cargo_configs() -> list[str]:
    """Cargo configuration files that apply besides the two compared ones."""
    homes = [Path(os.environ.get("CARGO_HOME", Path.home() / ".cargo"))]
    directories = {*ENGINE.parents, *ADAPTER.parents}
    found = [directory / ".cargo" / name for directory in sorted(directories) for name in ("config", "config.toml")]
    found += [home / name for home in homes for name in ("config", "config.toml")]
    return [str(path) for path in found if path.exists()]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    engine_manifest = tomllib.loads((ENGINE / "Cargo.toml").read_text())
    adapter_manifest = tomllib.loads((ADAPTER / "Cargo.toml").read_text())
    engine_config = tomllib.loads((ENGINE / ".cargo/config.toml").read_text())
    adapter_config = tomllib.loads((ADAPTER / ".cargo/config.toml").read_text())
    engine_locks = lock_versions(ENGINE / "Cargo.lock")
    adapter_locks = lock_versions(ADAPTER / "Cargo.lock")

    shared = sorted(set(engine_locks) & set(adapter_locks))
    differing = {name: {"engine": sorted(engine_locks[name]), "adapter": sorted(adapter_locks[name])}
                 for name in shared if engine_locks[name] != adapter_locks[name]}
    checks = {
        "profile_release": {
            "engine": engine_manifest["profile"]["release"],
            "adapter": adapter_manifest["profile"]["release"],
        },
        "target_rustflags": {
            "engine": engine_config["target"],
            "adapter": adapter_config["target"],
        },
        "shared_dependencies": {"count": len(shared), "differing": differing},
        "link_arguments": {
            "adapter_build_script": adapter_manifest["package"].get("build"),
            "adapter_build_rs_exists": (ADAPTER / "build.rs").exists(),
            "engine_build_rs_sha256": hashlib.sha256((ENGINE / "build.rs").read_bytes()).hexdigest(),
        },
        "other_configuration": {
            "cargo_config_files": other_cargo_configs(),
            "environment": sorted(name for name in os.environ if name.startswith(ENVIRONMENT)),
        },
    }
    failures = []
    if checks["profile_release"]["engine"] != checks["profile_release"]["adapter"]:
        failures.append("profile_release")
    if checks["target_rustflags"]["engine"] != checks["target_rustflags"]["adapter"]:
        failures.append("target_rustflags")
    if differing:
        failures.append("shared_dependencies")
    if checks["link_arguments"]["adapter_build_script"] not in (None, False) or \
            checks["link_arguments"]["adapter_build_rs_exists"]:
        failures.append("link_arguments")
    if any(checks["other_configuration"].values()):
        failures.append("other_configuration")
    report = {"checks": checks, "failures": failures, "identical": not failures}
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.write_text(text)
    print(text, end="")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
