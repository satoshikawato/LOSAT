#!/usr/bin/env python3
"""Enforce LOSAT's project-owned pure-Rust runtime boundary.

NCBI references:
- ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1379-1381
  qsort(hsp_list->hsp_array, hsp_list->hspcnt, sizeof(BlastHSP*),
        ScoreCompareHSPs);
- ncbi-blast/c++/src/algo/blast/core/link_hsps.c:483-486
  qsort(link_hsp_array, total_number_of_hsps, sizeof(LinkHSPStruct*),
        s_RevCompareHSPsTbx);

Those snippets define behavior that LOSAT ports from NCBI. This checker does
not implement that behavior; it prevents project production code from gaining
new imported/native implementation owners while the existing adapters are
removed in later, behavior-preserving pull requests.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence


INVENTORY_HEADER = ("rule", "path", "symbol", "classification")
TEMPORARY_DELEGATION = "temporary_production_algorithm_delegation"
REVIEWED_NON_ALGORITHM = "reviewed_non_algorithm_integration"

# NCBI references:
# - ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1379-1381
# - ncbi-blast/c++/src/algo/blast/core/link_hsps.c:483-486
# NCBI owns the comparator behavior at these call sites. LOSAT's completed
# migration owns sorting in Rust, so removed host-qsort adapter vocabulary is
# now a zero-state regression signal rather than an allowlist candidate.
OBSOLETE_QSORT_ADAPTER_NAMES = frozenset(
    {
        "native_hsp_qsort",
        "native_prelim_qsort",
        "native_qsort_blastn_hsps_by",
        "native_qsort_prelim_hits_by",
        "native_qsort_ungapped_hits_by",
        "output_qsort_native",
        "qsort_hit_indices",
        "qsort_hits_by",
        "qsort_subject_group_indices",
        "qsort_subject_groups_by",
        "qsort_ungapped_hits_by",
    }
)
OBSOLETE_QSORT_ADAPTER_PATHS = frozenset(
    {"LOSAT/src/algorithm/tblastx/ncbi_qsort.rs"}
)

IMPORT_BLOCK_RE = re.compile(
    r'\b(?:unsafe\s+)?extern\s+"(?P<abi>C|system)"\s*\{', re.MULTILINE
)
EXTERN_FUNCTION_RE = re.compile(
    r'\b(?:(?:pub(?:\([^)]*\))?)\s+)?(?:unsafe\s+)?extern\s+'
    r'"(?P<abi>C|system)"\s+fn\s+(?P<name>[A-Za-z_]\w*)\s*\(',
    re.MULTILINE,
)
FUNCTION_RE = re.compile(r"\bfn\s+([A-Za-z_]\w*)\s*(?:<[^>{}]*>)?\s*\(")
LINK_ATTRIBUTE_RE = re.compile(r"#\s*\[\s*(link(?:_name)?)\b([^\]]*)\]")
COMMAND_RE = re.compile(
    r"\b(?:[A-Za-z_]\w*\s*::\s*)*Command\s*::\s*new\s*\(",
    re.IGNORECASE,
)
COMMAND_ALIAS_RE = re.compile(
    r"\buse\s+std\s*::\s*process\s*::\s*(?:"
    r"Command\s+as\s+(?P<direct_alias>[A-Za-z_]\w*)|"
    r"\{[^}]*\bCommand\s+as\s+(?P<braced_alias>[A-Za-z_]\w*)[^}]*\}"
    r")\s*;"
)
DYNAMIC_LOADER_RE = re.compile(
    r"\b(?:libloading|dlopen|dlsym|LoadLibrary(?:A|W)?|GetProcAddress)\b"
)
LIBC_CALL_RE = re.compile(r"\blibc::(?P<symbol>[A-Za-z_]\w*)")
NATIVE_BUILD_RE = re.compile(
    r"\b(?:cc::Build|cmake::|bindgen::|Command\s*::\s*new\s*\(\s*"
    r'["\'](?:cc|gcc|clang|cl|cmake)["\'])'
)
NATIVE_DEPENDENCY_NAMES = {"bindgen", "cc", "cmake", "libc", "libloading"}


@dataclass(frozen=True, order=True)
class Finding:
    rule: str
    path: str
    symbol: str
    line: int
    classification: str
    detail: str

    @property
    def key(self) -> tuple[str, str, str]:
        return (self.rule, self.path, self.symbol)


@dataclass(frozen=True, order=True)
class Observation:
    kind: str
    path: str
    symbol: str
    line: int
    classification: str


@dataclass(frozen=True, order=True)
class AllowlistEntry:
    rule: str
    path: str
    symbol: str
    classification: str

    @property
    def key(self) -> tuple[str, str, str]:
        return (self.rule, self.path, self.symbol)


@dataclass
class AuditReport:
    production_files: int
    findings: list[Finding]
    observations: list[Observation]
    cargo_links_count: int


@dataclass
class Evaluation:
    allowed: list[tuple[Finding, AllowlistEntry]]
    unexpected: list[Finding]
    stale: list[AllowlistEntry]
    invalid_allowlist: list[str]

    @property
    def passed(self) -> bool:
        return not (self.unexpected or self.stale or self.invalid_allowlist)


def _line_number(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def _strip_rust_comments(text: str) -> str:
    """Replace Rust comments with spaces while preserving offsets/newlines."""

    chars = list(text)
    i = 0
    block_depth = 0
    in_line_comment = False
    in_string = False
    in_char = False
    escaped = False
    raw_terminator: str | None = None

    while i < len(chars):
        if in_line_comment:
            if chars[i] == "\n":
                in_line_comment = False
            else:
                chars[i] = " "
            i += 1
            continue

        if block_depth:
            if text.startswith("/*", i):
                chars[i] = chars[i + 1] = " "
                block_depth += 1
                i += 2
            elif text.startswith("*/", i):
                chars[i] = chars[i + 1] = " "
                block_depth -= 1
                i += 2
            else:
                if chars[i] != "\n":
                    chars[i] = " "
                i += 1
            continue

        if raw_terminator is not None:
            if text.startswith(raw_terminator, i):
                i += len(raw_terminator)
                raw_terminator = None
            else:
                i += 1
            continue

        if in_string:
            if escaped:
                escaped = False
            elif chars[i] == "\\":
                escaped = True
            elif chars[i] == '"':
                in_string = False
            i += 1
            continue

        if in_char:
            if escaped:
                escaped = False
            elif chars[i] == "\\":
                escaped = True
            elif chars[i] == "'":
                in_char = False
            i += 1
            continue

        raw_match = re.match(r'(?:br|rb|r)(?P<hashes>\#*)"', text[i:])
        if raw_match:
            hashes = raw_match.group("hashes")
            raw_terminator = '"' + hashes
            i += raw_match.end()
            continue
        if text.startswith("//", i):
            chars[i] = chars[i + 1] = " "
            in_line_comment = True
            i += 2
            continue
        if text.startswith("/*", i):
            chars[i] = chars[i + 1] = " "
            block_depth = 1
            i += 2
            continue
        if chars[i] == '"':
            in_string = True
        elif chars[i] == "'" and i + 2 < len(chars):
            # Lifetimes have no closing quote; only enter char mode when a
            # plausible closing quote is nearby.
            close = text.find("'", i + 1, min(i + 8, len(chars)))
            if close != -1:
                in_char = True
        i += 1

    return "".join(chars)


def _matching_brace(text: str, opening: int) -> int:
    depth = 0
    for index in range(opening, len(text)):
        if text[index] == "{":
            depth += 1
        elif text[index] == "}":
            depth -= 1
            if depth == 0:
                return index
    return len(text) - 1


def _enclosing_function(text: str, offset: int) -> str:
    matches = list(FUNCTION_RE.finditer(text, 0, offset))
    return matches[-1].group(1) if matches else "module_scope"


def _relative(path: Path, root: Path) -> str:
    return path.resolve().relative_to(root.resolve()).as_posix()


def _scan_rust_file(path: Path, root: Path) -> tuple[list[Finding], list[Observation]]:
    raw = path.read_text(encoding="utf-8")
    text = _strip_rust_comments(raw)
    rel = _relative(path, root)
    findings: list[Finding] = []
    observations: list[Observation] = []
    import_spans: list[tuple[int, int]] = []
    imported_symbols: list[tuple[str, int]] = []

    for match in IMPORT_BLOCK_RE.finditer(text):
        opening = text.find("{", match.start(), match.end())
        closing = _matching_brace(text, opening)
        import_spans.append((match.start(), closing + 1))
        block = text[opening + 1 : closing]
        for symbol_match in re.finditer(r"\bfn\s+([A-Za-z_]\w*)\s*\(", block):
            symbol = symbol_match.group(1)
            absolute = opening + 1 + symbol_match.start()
            imported_symbols.append((symbol, absolute))
            findings.append(
                Finding(
                    "rust.imported_abi",
                    rel,
                    symbol,
                    _line_number(text, absolute),
                    "production_native_boundary_candidate",
                    f'imported extern "{match.group("abi")}" function',
                )
            )

    for match in EXTERN_FUNCTION_RE.finditer(text):
        if any(start <= match.start() < end for start, end in import_spans):
            continue
        observations.append(
            Observation(
                "rust.implemented_abi",
                rel,
                match.group("name"),
                _line_number(text, match.start()),
                "rust_implemented_export_or_callback",
            )
        )

    without_imports = list(text)
    for start, end in import_spans:
        for index in range(start, end):
            if without_imports[index] != "\n":
                without_imports[index] = " "
    call_text = "".join(without_imports)
    for symbol, _ in imported_symbols:
        for call in re.finditer(rf"\b{re.escape(symbol)}\s*\(", call_text):
            owner = _enclosing_function(call_text, call.start())
            findings.append(
                Finding(
                    "rust.imported_abi_call",
                    rel,
                    f"{symbol}@{owner}",
                    _line_number(call_text, call.start()),
                    "production_native_call_candidate",
                    f"call to imported symbol {symbol} from {owner}",
                )
            )

    for match in LINK_ATTRIBUTE_RE.finditer(text):
        attribute = re.sub(r"\s+", " ", match.group(0)).strip()
        findings.append(
            Finding(
                "rust.link_attribute",
                rel,
                attribute,
                _line_number(text, match.start()),
                "production_native_link_candidate",
                "project-authored native link attribute",
            )
        )

    for match in DYNAMIC_LOADER_RE.finditer(text):
        findings.append(
            Finding(
                "rust.dynamic_loader",
                rel,
                match.group(0),
                _line_number(text, match.start()),
                "production_dynamic_load_candidate",
                "dynamic-loader token in production source",
            )
        )

    for match in LIBC_CALL_RE.finditer(text):
        findings.append(
            Finding(
                "rust.libc_call",
                rel,
                match.group("symbol"),
                _line_number(text, match.start()),
                "production_system_call_candidate",
                "direct libc call in project production source",
            )
        )

    command_matches = list(COMMAND_RE.finditer(text))
    for alias_match in COMMAND_ALIAS_RE.finditer(text):
        alias = alias_match.group("direct_alias") or alias_match.group("braced_alias")
        command_matches.extend(
            re.finditer(rf"\b{re.escape(alias)}\s*::\s*new\s*\(", text)
        )
    for match in sorted(command_matches, key=lambda item: item.start()):
        owner = _enclosing_function(text, match.start())
        findings.append(
            Finding(
                "rust.production_subprocess",
                rel,
                f"Command::new@{owner}",
                _line_number(text, match.start()),
                "production_subprocess_review_required",
                "production subprocess launch can delegate runtime behavior",
            )
        )

    for name in sorted(OBSOLETE_QSORT_ADAPTER_NAMES):
        match = re.search(rf"\b{re.escape(name)}\b", text)
        if match is not None:
            findings.append(
                Finding(
                    "rust.obsolete_qsort_adapter",
                    rel,
                    name,
                    _line_number(text, match.start()),
                    "obsolete_production_adapter",
                    "removed project-owned host-qsort adapter name",
                )
            )

    return findings, observations


def _scan_obsolete_adapter_paths(root: Path) -> list[Finding]:
    findings: list[Finding] = []
    for rel in sorted(OBSOLETE_QSORT_ADAPTER_PATHS):
        if (root / rel).exists():
            findings.append(
                Finding(
                    "path.obsolete_qsort_adapter",
                    rel,
                    Path(rel).name,
                    1,
                    "obsolete_production_adapter",
                    "removed project-owned host-qsort adapter path",
                )
            )
    return findings


def _dependency_tables(value: object, path: tuple[str, ...] = ()) -> Iterable[tuple[str, Mapping[str, object]]]:
    if not isinstance(value, dict):
        return
    for key, nested in value.items():
        next_path = path + (str(key),)
        if key in {"dependencies", "build-dependencies"} and isinstance(nested, dict):
            yield (".".join(next_path), nested)
        if isinstance(nested, dict):
            yield from _dependency_tables(nested, next_path)


def _scan_manifest(manifest: Path, root: Path) -> list[Finding]:
    if not manifest.exists():
        return []
    data = tomllib.loads(manifest.read_text(encoding="utf-8"))
    rel = _relative(manifest, root)
    findings: list[Finding] = []
    for table_name, dependencies in _dependency_tables(data):
        for name, specification in sorted(dependencies.items()):
            package_name = name
            if isinstance(specification, dict):
                declared_package = specification.get("package")
                if isinstance(declared_package, str):
                    package_name = declared_package
            if package_name in NATIVE_DEPENDENCY_NAMES:
                dependency_identity = (
                    name
                    if name == package_name
                    else f"{name}->package={package_name}"
                )
                findings.append(
                    Finding(
                        "cargo.native_dependency",
                        rel,
                        f"{dependency_identity}@{table_name}",
                        1,
                        "native_dependency_review_required",
                        "direct normal/build dependency can expose native code or ABI",
                    )
                )
    return findings


def _scan_build_scripts(root: Path) -> list[Finding]:
    findings: list[Finding] = []
    excluded = {".git", ".agents", ".tmp", "target", "tests"}
    build_scripts: list[Path] = []
    for current, directories, files in os.walk(root):
        directories[:] = sorted(directory for directory in directories if directory not in excluded)
        if "build.rs" in files:
            build_scripts.append(Path(current) / "build.rs")
    for path in sorted(build_scripts):
        text = _strip_rust_comments(path.read_text(encoding="utf-8"))
        rel = _relative(path, root)
        findings.append(
            Finding(
                "build.project_script",
                rel,
                "build.rs",
                1,
                "project_build_script_review_required",
                "project build script can introduce a native build route",
            )
        )
        for match in NATIVE_BUILD_RE.finditer(text):
            findings.append(
                Finding(
                    "build.native_code",
                    rel,
                    re.sub(r"\s+", " ", match.group(0)),
                    _line_number(text, match.start()),
                    "build_time_native_code_candidate",
                    "native compiler/binding tool invoked by project build script",
                )
            )
    return findings


def load_cargo_metadata(root: Path) -> Mapping[str, object]:
    manifest = root / "LOSAT" / "Cargo.toml"
    command = [
        "cargo",
        "metadata",
        "--format-version",
        "1",
        "--locked",
        "--manifest-path",
        str(manifest),
    ]
    completed = subprocess.run(
        command,
        cwd=root,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"cargo metadata failed ({completed.returncode}): {detail}")
    return json.loads(completed.stdout)


def _scan_cargo_links(metadata: Mapping[str, object]) -> list[Finding]:
    packages = {package["id"]: package for package in metadata.get("packages", [])}
    resolve = metadata.get("resolve") or {}
    root_id = resolve.get("root")
    nodes = {node["id"]: node for node in resolve.get("nodes", [])}
    if root_id is None:
        return []

    reachable_kinds: dict[str, set[str]] = {root_id: {"root"}}
    queue = [root_id]
    while queue:
        package_id = queue.pop(0)
        node = nodes.get(package_id, {})
        for dependency in node.get("deps", []):
            kinds = {
                (entry.get("kind") or "normal")
                for entry in dependency.get("dep_kinds", [])
            }
            kinds &= {"normal", "build"}
            if not kinds:
                continue
            dependency_id = dependency["pkg"]
            current = reachable_kinds.setdefault(dependency_id, set())
            before = len(current)
            current.update(kinds)
            if len(current) != before:
                queue.append(dependency_id)

    findings: list[Finding] = []
    for package_id, kinds in sorted(reachable_kinds.items()):
        package = packages.get(package_id)
        if not package or not package.get("links"):
            continue
        name = package["name"]
        version = package["version"]
        links = package["links"]
        kind_text = ",".join(sorted(kinds - {"root"})) or "normal"
        findings.append(
            Finding(
                "cargo.links_metadata",
                "LOSAT/Cargo.lock",
                f"{name}@{version}|links={links}|kinds={kind_text}",
                1,
                "review_required_metadata_not_delegation",
                "resolved normal/build dependency declares Cargo links metadata",
            )
        )
    return findings


def audit_repository(
    root: Path, *, metadata: Mapping[str, object] | None = None
) -> AuditReport:
    root = root.resolve()
    source_root = root / "LOSAT" / "src"
    rust_files = sorted(source_root.rglob("*.rs")) if source_root.exists() else []
    findings: list[Finding] = []
    observations: list[Observation] = []
    for path in rust_files:
        file_findings, file_observations = _scan_rust_file(path, root)
        findings.extend(file_findings)
        observations.extend(file_observations)

    findings.extend(_scan_manifest(root / "LOSAT" / "Cargo.toml", root))
    findings.extend(_scan_build_scripts(root))
    findings.extend(_scan_obsolete_adapter_paths(root))
    if metadata is None:
        metadata = load_cargo_metadata(root)
    cargo_links = _scan_cargo_links(metadata)
    findings.extend(cargo_links)

    unique_findings: dict[tuple[str, str, str], Finding] = {}
    for finding in sorted(findings):
        if finding.key in unique_findings:
            first = unique_findings[finding.key]
            raise RuntimeError(
                "boundary finding key is not exact: "
                f"{finding.key} occurs at lines {first.line} and {finding.line}"
            )
        unique_findings[finding.key] = finding

    return AuditReport(
        production_files=len(rust_files),
        findings=sorted(unique_findings.values()),
        observations=sorted(set(observations)),
        cargo_links_count=len(cargo_links),
    )


def load_allowlist(path: Path) -> list[AllowlistEntry]:
    lines = path.read_text(encoding="utf-8").splitlines()
    rows = [line for line in lines if line.strip() and not line.lstrip().startswith("#")]
    if not rows:
        raise ValueError(f"allowlist is empty: {path}")
    header = tuple(rows[0].split("\t"))
    if header != INVENTORY_HEADER:
        raise ValueError(
            f"inventory header must be {INVENTORY_HEADER}, found {header}: {path}"
        )
    entries: list[AllowlistEntry] = []
    seen: set[tuple[str, str, str]] = set()
    for line_number, row in enumerate(rows[1:], start=2):
        fields = row.split("\t")
        if len(fields) != len(INVENTORY_HEADER):
            raise ValueError(f"invalid inventory row {line_number}: expected 4 fields")
        entry = AllowlistEntry(*fields)
        if entry.key in seen:
            raise ValueError(f"duplicate allowlist key at row {line_number}: {entry.key}")
        seen.add(entry.key)
        entries.append(entry)
    return sorted(entries)


def evaluate_findings(
    findings: Sequence[Finding],
    allowlist: Sequence[AllowlistEntry],
    *,
    expected_classification: str = TEMPORARY_DELEGATION,
    require_empty_allowlist: bool = False,
) -> Evaluation:
    finding_by_key = {finding.key: finding for finding in findings}
    invalid: list[str] = []
    allow_by_key: dict[tuple[str, str, str], AllowlistEntry] = {}
    for entry in allowlist:
        if require_empty_allowlist:
            invalid.append(
                f"{entry.rule} {entry.path} {entry.symbol}: finalized production "
                "delegation allowlist must contain zero entries"
            )
            continue
        allow_by_key[entry.key] = entry
    unexpected = sorted(
        finding for key, finding in finding_by_key.items() if key not in allow_by_key
    )
    stale = sorted(entry for key, entry in allow_by_key.items() if key not in finding_by_key)
    allowed: list[tuple[Finding, AllowlistEntry]] = []

    for key in sorted(finding_by_key.keys() & allow_by_key.keys()):
        finding = finding_by_key[key]
        entry = allow_by_key[key]
        if entry.classification != expected_classification:
            invalid.append(
                f"{entry.rule} {entry.path} {entry.symbol}: expected classification "
                f"{expected_classification}, found {entry.classification}"
            )
            continue
        is_dependency_signal = finding.rule == "cargo.links_metadata"
        expects_dependency_signal = expected_classification == REVIEWED_NON_ALGORITHM
        if is_dependency_signal != expects_dependency_signal:
            invalid.append(
                f"{entry.rule} {entry.path} {entry.symbol}: dependency metadata "
                "signals and project-owned production findings must use separate inventories"
            )
            continue
        allowed.append((finding, entry))

    return Evaluation(allowed, unexpected, stale, sorted(invalid))


def _print_evaluation(category: str, accepted_label: str, evaluation: Evaluation) -> None:
    for finding, entry in evaluation.allowed:
        print(
            f"{accepted_label}\t"
            f"{finding.rule}\t{finding.path}:{finding.line}\t{finding.symbol}\t"
            f"{entry.classification}"
        )
    for finding in evaluation.unexpected:
        print(
            f"UNEXPECTED_{category}\t"
            f"{finding.rule}\t{finding.path}:{finding.line}\t{finding.symbol}\t"
            f"{finding.classification}\t{finding.detail}"
        )
    for entry in evaluation.stale:
        print(
            f"STALE_{category}\t"
            f"{entry.rule}\t{entry.path}\t{entry.symbol}\t{entry.classification}"
        )
    for error in evaluation.invalid_allowlist:
        print(f"INVALID_{category}\t{error}")


def _print_report(
    report: AuditReport,
    production_evaluation: Evaluation,
    dependency_evaluation: Evaluation,
) -> None:
    print("LOSAT pure-Rust runtime boundary audit")
    print(f"production_rust_files\t{report.production_files}")
    print(
        "temporary_production_findings\t"
        f"{len([finding for finding in report.findings if finding.rule != 'cargo.links_metadata'])}"
    )
    print(f"cargo_links_findings\t{report.cargo_links_count}")
    print(f"rust_implemented_abi\t{len(report.observations)}")
    for observation in report.observations:
        print(
            "INFO\t"
            f"{observation.kind}\t{observation.path}:{observation.line}\t"
            f"{observation.symbol}\t{observation.classification}"
        )
    _print_evaluation("TEMPORARY", "ALLOW_TEMPORARY", production_evaluation)
    _print_evaluation("DEPENDENCY", "REVIEWED_DEPENDENCY", dependency_evaluation)
    passed = production_evaluation.passed and dependency_evaluation.passed
    print(
        "summary\t"
        f"temporary_allowed={len(production_evaluation.allowed)}\t"
        f"temporary_unexpected={len(production_evaluation.unexpected)}\t"
        f"temporary_stale={len(production_evaluation.stale)}\t"
        f"temporary_invalid={len(production_evaluation.invalid_allowlist)}\t"
        f"dependencies_reviewed={len(dependency_evaluation.allowed)}\t"
        f"dependencies_unexpected={len(dependency_evaluation.unexpected)}\t"
        f"dependencies_stale={len(dependency_evaluation.stale)}\t"
        f"dependencies_invalid={len(dependency_evaluation.invalid_allowlist)}"
    )
    print(f"result\t{'PASS' if passed else 'FAIL'}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument("--allowlist", type=Path)
    parser.add_argument("--dependency-review", type=Path)
    parser.add_argument(
        "--metadata-file",
        type=Path,
        help="read saved cargo metadata JSON instead of invoking cargo metadata",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    root = args.root.resolve()
    allowlist_path = args.allowlist or root / "LOSAT" / "tests" / "pure_rust_runtime_allowlist.tsv"
    dependency_review_path = (
        args.dependency_review
        or root / "LOSAT" / "tests" / "pure_rust_runtime_dependency_review.tsv"
    )
    try:
        metadata = None
        if args.metadata_file:
            metadata = json.loads(args.metadata_file.read_text(encoding="utf-8"))
        report = audit_repository(root, metadata=metadata)
        allowlist = load_allowlist(allowlist_path)
        dependency_review = load_allowlist(dependency_review_path)
        production_findings = [
            finding for finding in report.findings if finding.rule != "cargo.links_metadata"
        ]
        dependency_signals = [
            finding for finding in report.findings if finding.rule == "cargo.links_metadata"
        ]
        production_evaluation = evaluate_findings(
            production_findings,
            allowlist,
            expected_classification=TEMPORARY_DELEGATION,
            require_empty_allowlist=True,
        )
        dependency_evaluation = evaluate_findings(
            dependency_signals,
            dependency_review,
            expected_classification=REVIEWED_NON_ALGORITHM,
        )
    except (OSError, RuntimeError, ValueError, json.JSONDecodeError) as error:
        print(f"boundary checker error: {error}", file=sys.stderr)
        return 2
    _print_report(report, production_evaluation, dependency_evaluation)
    return 0 if production_evaluation.passed and dependency_evaluation.passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
