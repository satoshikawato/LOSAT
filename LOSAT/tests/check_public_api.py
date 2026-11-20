#!/usr/bin/env python3
"""Keep the public API of `losat::api` equal to its reviewed listing, `public_api.txt`.

No NCBI counterpart: this checks the library's own surface (docs/losat_library_surface_plan.md
§14), not ported behavior.

Only `api` is public: `lib.rs` hides the other modules with `#[doc(hidden)]` or keeps them
`pub(crate)`. cargo-public-api
lists that surface from rustdoc JSON, one line per item, and this check fails when the listing
differs from `public_api.txt`. A change of the surface is then a line in the pull request's diff
that a reviewer reads. After a deliberate change, run `update` and commit the file; when the check
fails in CI, the job also uploads the new listing (artifact `public-api`).

cargo-public-api builds rustdoc JSON with a nightly toolchain (it switches a stable one to
`nightly`) and reads one rustdoc JSON format: 0.52.0 reads format 57, which the nightlies from
2025-11-22 to 2026-06-23 write. Both are pinned below and in `.github/workflows/ci.yml` (job
`public-api`); `test_check_public_api.py` checks that the two agree. This script installs nothing:
a missing toolchain or tool is an error that names the command that installs it.

Usage:
  python3 LOSAT/tests/check_public_api.py check [--out FILE]
  python3 LOSAT/tests/check_public_api.py update
"""

from __future__ import annotations

import argparse
import difflib
import os
import subprocess
import sys
from pathlib import Path
from typing import Sequence

TOOL = "cargo-public-api"
TOOL_VERSION = "0.52.0"
TOOLCHAIN = "nightly-2026-05-24"
SNAPSHOT = Path("LOSAT/tests/public_api.txt")
INSTALL_TOOLCHAIN = f"rustup toolchain install --no-self-update --profile minimal {TOOLCHAIN}"
INSTALL_TOOL = f"cargo install --locked {TOOL}@{TOOL_VERSION}"


def listing_command(root: Path) -> list[str]:
    # `--simplified` once leaves out blanket impls (`impl<T> From<T> for T` ...). Auto traits
    # (`Send`, `Sync`, ...) and derived impls stay listed: losing one breaks a caller.
    manifest = str(root / "LOSAT" / "Cargo.toml")
    return ["cargo", f"+{TOOLCHAIN}", "public-api", "--manifest-path", manifest, "--simplified", "--color", "never"]


def toolchain_names(rustup_list: str) -> list[str]:
    """The names in `rustup toolchain list` output (`nightly-2026-05-24-x86_64-... (default)`)."""
    return [line.split()[0] for line in rustup_list.splitlines() if line.strip()]


def has_toolchain(names: Sequence[str], toolchain: str) -> bool:
    return any(name == toolchain or name.startswith(toolchain + "-") for name in names)


def tool_version(version_output: str) -> str | None:
    """`0.52.0` from `cargo-public-api 0.52.0`."""
    words = version_output.split()
    if len(words) >= 2 and words[0] == TOOL:
        return words[1]
    return None


def listing_diff(expected: str, actual: str) -> list[str]:
    """Unified diff lines from the snapshot to the current listing; empty when they are equal."""
    if expected == actual:
        return []
    return list(
        difflib.unified_diff(
            expected.splitlines(keepends=True),
            actual.splitlines(keepends=True),
            fromfile=str(SNAPSHOT),
            tofile="cargo public-api",
        )
    )


def environment_error() -> str | None:
    """Why the listing cannot be made here, or None. Never installs a toolchain or a tool."""
    try:
        listed = subprocess.run(
            ["rustup", "toolchain", "list"], capture_output=True, text=True, check=True
        ).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        return f"rustup toolchain list failed: {error}"
    if not has_toolchain(toolchain_names(listed), TOOLCHAIN):
        return f"toolchain {TOOLCHAIN} is not installed; install it with: {INSTALL_TOOLCHAIN}"
    try:
        found = subprocess.run([TOOL, "--version"], capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return f"{TOOL} is not on PATH; install it with: {INSTALL_TOOL}"
    version = tool_version(found)
    if version != TOOL_VERSION:
        return f"{TOOL} {version} found, {TOOL_VERSION} expected; install it with: {INSTALL_TOOL}"
    return None


def current_listing(root: Path) -> str:
    env = dict(os.environ)
    # rustup >= 1.28: never install the toolchain behind `+nightly-...` on the fly.
    env["RUSTUP_AUTO_INSTALL"] = "0"
    done = subprocess.run(listing_command(root), capture_output=True, text=True, env=env)
    if done.returncode != 0:
        sys.stderr.write(done.stderr)
        raise RuntimeError(f"cargo public-api exited with {done.returncode}")
    return done.stdout


def write_summary(diff: Sequence[str]) -> None:
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if not summary:
        return
    with open(summary, "a", encoding="utf-8") as out:
        out.write("### The public API of `losat::api` changed\n\n```diff\n")
        out.writelines(diff)
        out.write("```\n")


def parse_args(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("mode", choices=("check", "update"))
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--out", type=Path, help="also write the current listing to this file")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    problem = environment_error()
    if problem:
        print(f"check_public_api: {problem}", file=sys.stderr)
        return 2
    try:
        actual = current_listing(root)
    except RuntimeError as error:
        print(f"check_public_api: {error}", file=sys.stderr)
        return 2
    if args.out:
        args.out.write_text(actual, encoding="utf-8")
    snapshot = root / SNAPSHOT
    if args.mode == "update":
        snapshot.write_text(actual, encoding="utf-8")
        print(f"check_public_api: wrote {SNAPSHOT} ({len(actual.splitlines())} lines)")
        return 0
    expected = snapshot.read_text(encoding="utf-8") if snapshot.exists() else ""
    diff = listing_diff(expected, actual)
    if not diff:
        print(f"check_public_api: the public API equals {SNAPSHOT} ({len(actual.splitlines())} lines)")
        return 0
    sys.stdout.writelines(diff)
    write_summary(diff)
    print(
        f"check_public_api: the public API differs from {SNAPSHOT}. If the change is deliberate, run "
        "`python3 LOSAT/tests/check_public_api.py update` (or take the CI artifact `public-api`) and "
        "commit the file with the change.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
