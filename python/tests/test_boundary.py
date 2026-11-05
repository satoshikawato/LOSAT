"""The rules of python/AGENTS.md that a test can read in the sources."""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

PYTHON = Path(__file__).resolve().parents[1]
LOSAT = PYTHON.parent / "LOSAT"


def toml(path: Path) -> dict:
    return tomllib.loads(path.read_text())


def test_the_engine_is_reached_through_the_api_only():
    for path in (PYTHON / "src").glob("*.rs"):
        assert set(re.findall(r"\bLOSAT::(\w+)", path.read_text())) <= {"api"}, path


def test_the_process_stays_the_hosts():
    for path in (PYTHON / "src").glob("*.rs"):
        text = path.read_text()
        for word in ("build_global", "set_hook", "set_var", "set_current_dir", "Command::new"):
            assert word not in text, f"{word} in {path}"


def test_the_build_is_losats_but_for_panics():
    binding = toml(PYTHON / "Cargo.toml")["profile"]["release"]
    engine = toml(LOSAT / "Cargo.toml")["profile"]["release"]
    assert binding.pop("panic") == "unwind"
    assert engine.pop("panic") == "abort"
    assert binding == engine

    def rustflags(path: Path) -> list[str]:
        return toml(path)["target"]['cfg(target_arch = "x86_64")']["rustflags"]

    assert rustflags(PYTHON / ".cargo/config.toml") == rustflags(LOSAT / ".cargo/config.toml")
