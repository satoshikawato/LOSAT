#!/usr/bin/env python3
"""Count LOSAT's ports of NCBI BLAST+ unit-test cases against a ledger.

NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/blasthits_unit_test.cpp:1094-1099
```c++
BOOST_AUTO_TEST_CASE(testHSPListSort)
{
    const int kHspCnt = 10;
    setupHSPList(kHspCnt);
    Blast_HSPListSortByScore(m_HspList);
    BOOST_REQUIRE(Blast_HSPListIsSortedByScore(m_HspList));
```

NCBI keeps one Boost.Test case per checked behavior. LOSAT ports a case as a
Rust test next to the implementation, with the line

    // NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/blasthits_unit_test.cpp:1094-1099 testHSPListSort

above the test function, and a row in docs/evidence/ncbi_unit_cases/LEDGER.tsv.

Commands:
- extract: list the cases of the NCBI unit tests in scope (needs the NCBI
  source; run locally) into NCBI_CASES.tsv.
- scan: list the citation lines in LOSAT/src and LOSAT/tests and the test
  function under each.
- check: fail when a citation, a ledger row and the code disagree (CI; needs
  no NCBI source).
- report: Markdown table of ledger statuses per NCBI module.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence


CASES_HEADER = ("module", "file", "case", "line_start", "line_end", "kind", "note")
LEDGER_HEADER = ("module", "case", "status", "losat_test", "note", "since")
KINDS = ("AUTO", "FIXTURE", "DISABLED", "MACRO")
STATUSES = ("ported", "partial", "e2e", "to-port", "n-a", "superseded")
LINKED_STATUSES = ("ported", "partial")
COMMIT_KEY = "ncbi_commit"
DEFAULT_CASES = Path("docs/evidence/ncbi_unit_cases/NCBI_CASES.tsv")
DEFAULT_LEDGER = Path("docs/evidence/ncbi_unit_cases/LEDGER.tsv")
SCAN_DIRS = (Path("LOSAT/src"), Path("LOSAT/tests"))

# NCBI unit tests in scope (PLAN 2026-10-10 §2.3): the api modules that cover
# algo/blast/core and the input/format layers of the five LOSAT programs.
API_MODULES = (
    "aalookup", "aascan", "blast", "blastdiag", "blastengine", "blastextend",
    "blastfilter", "blasthits", "blastoptions", "blastsetup", "gapinfo",
    "gencode_singleton", "hspfilter_besthit", "hspfilter_culling", "hspstream",
    "linkhsp", "ntlookup", "ntscan", "nuclwordfinder", "optionshandle",
    "querydata", "queryinfo", "redoalignment", "scoreblk", "split_query",
    "stat", "subj_ranges", "traceback", "tracebacksearch", "bl2seq",
)
# (module, directory under c++/src, file glob). An api file is its own module;
# the other directories are one module each.
SOURCES = tuple(
    (name, "algo/blast/unit_tests/api", f"{name}_unit_test.cpp") for name in API_MODULES
) + (
    ("blast_format", "algo/blast/unit_tests/blast_format", "*.cpp"),
    ("blast_mt", "algo/blast/unit_tests/blast_mt", "*.cpp"),
    ("blastinput", "algo/blast/blastinput/unit_test", "blastinput_unit_test.cpp"),
    ("align_format", "objtools/align_format/unit_test", "*.cpp"),
)
# `#if` conditions that are false in NCBI's build: dead code, and
# SEQLOC_MIX_QUERY_OK, which no build defines.
FALSE_CONDITIONS = frozenset({"0", "SEQLOC_MIX_QUERY_OK"})
TRUE_CONDITIONS = frozenset({"1"})

CASE_RE = re.compile(r"\bBOOST_(AUTO|FIXTURE)_TEST_CASE\s*\(\s*([A-Za-z_]\w*)")
DIRECTIVE_RE = re.compile(r"^\s*#\s*(if|ifdef|ifndef|elif|else|endif)\b(.*)$")
DEFINE_RE = re.compile(r"^\s*#\s*define\b")
MACRO_CASE_RE = re.compile(
    r"^\s*#\s*define\s+(?P<name>[A-Za-z_]\w*)\((?P<params>[^)]*)\).*?"
    r"\bBOOST_(?P<kind>AUTO|FIXTURE)_TEST_CASE\s*\(\s*(?P<expr>[^,)]*?)\s*[,)]"
)
CITATION_MARK_RE = re.compile(r"NCBI unit test\s*\(")
CITATION_RE = re.compile(
    r"^\s*//+\s*NCBI unit test \((?P<commit>[0-9a-f]{7,40})\): "
    r"(?P<file>c\+\+/\S+?):(?P<lines>\d+(?:-\d+)?(?:,\d+(?:-\d+)?)*) "
    r"(?P<case>[A-Za-z_]\w*)\s*$"
)
FN_RE = re.compile(
    r"^\s*(?:pub(?:\([^)]*\))?\s+)?(?:const\s+)?(?:async\s+)?(?:unsafe\s+)?"
    r"(?:extern\s+\"[^\"]*\"\s+)?fn\s+([A-Za-z_]\w*)"
)
ATTRIBUTE_RE = re.compile(r"^\s*#\s*!?\s*\[")
TEST_ATTRIBUTE_RE = re.compile(r"^\s*#\s*\[\s*(?:tokio\s*::\s*)?test\b")
IGNORE_ATTRIBUTE_RE = re.compile(r"^\s*#\s*\[\s*ignore\b")
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")


class LedgerError(Exception):
    """A file that cannot be read as a case list or ledger."""


@dataclass(frozen=True)
class Case:
    module: str
    file: str
    case: str
    line_start: int
    line_end: int
    kind: str
    note: str = ""


@dataclass(frozen=True)
class Citation:
    path: str
    line: int
    commit: str
    file: str
    ranges: tuple[tuple[int, int], ...]
    case: str
    function: str | None
    is_test: bool
    is_ignored: bool


@dataclass(frozen=True)
class LedgerRow:
    module: str
    case: str
    status: str
    losat_test: str
    note: str
    since: str
    row: int

    @property
    def tests(self) -> list[str]:
        return [item.strip() for item in self.losat_test.split(";") if item.strip()]


@dataclass(frozen=True)
class Failure:
    module: str
    case: str
    reason: str


@dataclass
class CaseList:
    commit: str
    cases: list[Case]
    by_file_case: dict[tuple[str, str], Case] = field(default_factory=dict)
    by_module_case: dict[tuple[str, str], Case] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for case in self.cases:
            self.by_file_case[(case.file, case.case)] = case
            self.by_module_case[(case.module, case.case)] = case


# --- extract ---------------------------------------------------------------


def blank_cpp_comments_and_literals(text: str) -> str:
    """Blank C++ comments and the inside of literals, keeping offsets and newlines."""

    out = list(text)
    i = 0
    n = len(text)

    def blank(start: int, end: int) -> None:
        for index in range(start, end):
            if out[index] != "\n":
                out[index] = " "

    while i < n:
        if text.startswith("//", i):
            end = text.find("\n", i)
            end = n if end == -1 else end
            blank(i, end)
            i = end
        elif text.startswith("/*", i):
            end = text.find("*/", i + 2)
            end = n if end == -1 else end + 2
            blank(i, end)
            i = end
        elif text[i] == "R" and text.startswith('R"', i) and (i == 0 or not (text[i - 1].isalnum() or text[i - 1] == "_")):
            open_paren = text.find("(", i + 2)
            delimiter = text[i + 2 : open_paren] if open_paren != -1 else ""
            close = text.find(")" + delimiter + '"', open_paren + 1) if open_paren != -1 else -1
            end = n if close == -1 else close + len(delimiter) + 2
            blank(i + 2, end - 1)
            i = end
        elif text[i] in "\"'":
            quote = text[i]
            j = i + 1
            while j < n and text[j] != quote and text[j] != "\n":
                j += 2 if text[j] == "\\" else 1
            blank(i + 1, min(j, n))
            i = j + 1
        else:
            i += 1
    return "".join(out)


def _disabled_lines(lines: Sequence[str]) -> list[bool]:
    """Mark the lines that NCBI's build compiles out (`#if 0`, SEQLOC_MIX_QUERY_OK)."""

    # Each frame is "off", "on" or "unknown"; a line is disabled under any "off".
    stack: list[str] = []
    disabled: list[bool] = []
    for line in lines:
        match = DIRECTIVE_RE.match(line)
        if match:
            directive, rest = match.group(1), match.group(2).strip()
            if directive == "if":
                state = "off" if rest in FALSE_CONDITIONS else "on" if rest in TRUE_CONDITIONS else "unknown"
                stack.append(state)
            elif directive in ("ifdef", "ifndef"):
                stack.append("unknown")
            elif directive == "elif" and stack:
                stack[-1] = "off" if stack[-1] == "on" else "unknown"
            elif directive == "else" and stack:
                stack[-1] = {"off": "on", "on": "off"}.get(stack[-1], "unknown")
            elif directive == "endif" and stack:
                stack.pop()
        disabled.append("off" in stack)
    return disabled


def _matching_brace(text: str, opening: int) -> int:
    depth = 0
    for index in range(opening, len(text)):
        if text[index] == "{":
            depth += 1
        elif text[index] == "}":
            depth -= 1
            if depth == 0:
                return index
    raise LedgerError("unbalanced braces")


def extract_cases_from_text(module: str, file: str, text: str) -> list[Case]:
    blanked = blank_cpp_comments_and_literals(text)
    line_starts = [0]
    for index, char in enumerate(blanked):
        if char == "\n":
            line_starts.append(index + 1)

    def line_of(offset: int) -> int:
        low, high = 0, len(line_starts) - 1
        while low < high:
            mid = (low + high + 1) // 2
            if line_starts[mid] <= offset:
                low = mid
            else:
                high = mid - 1
        return low + 1

    lines = blanked.split("\n")
    disabled = _disabled_lines(lines)
    defines = _define_spans(lines)

    def in_define(offset: int) -> bool:
        line = line_of(offset)
        return any(start <= line <= end for start, end, _ in defines)

    def kind_of(line: int, macro_kind: str) -> str:
        if disabled[line - 1]:
            return "DISABLED"
        return "FIXTURE" if macro_kind == "FIXTURE" else "AUTO"

    cases: list[Case] = []
    for match in CASE_RE.finditer(blanked):
        if in_define(match.start()):
            continue
        line_start = line_of(match.start())
        opening = blanked.find("{", match.end())
        if opening == -1:
            raise LedgerError(f"{file}:{line_start}: no body for {match.group(2)}")
        try:
            line_end = line_of(_matching_brace(blanked, opening))
        except LedgerError as error:
            raise LedgerError(f"{file}:{line_start}: {error}") from None
        cases.append(Case(module, file, match.group(2), line_start, line_end, kind_of(line_start, match.group(1))))

    # A macro that declares a case (ntscan's DECLARE_TEST) is one MACRO row as
    # written in the source, named after the macro, from the `#define` to its
    # last invocation; the note gives how many cases it declares at run time.
    for define_start, _, definition in defines:
        macro = MACRO_CASE_RE.match(definition)
        if not macro:
            continue
        name = macro.group("name")
        params = [param.strip() for param in macro.group("params").split(",")]
        live = dead = 0
        first = last = 0
        for call in re.finditer(rf"\b{re.escape(name)}\s*\(", blanked):
            if in_define(call.start()):
                continue
            line = line_of(call.start())
            close = _matching_paren(blanked, call.end() - 1)
            args = [arg.strip() for arg in blanked[call.end() : close].split(",")]
            if len(args) != len(params):
                raise LedgerError(f"{file}:{line}: {name} takes {len(params)} arguments")
            if disabled[line - 1]:
                dead += 1
            else:
                live += 1
            first = first or line
            last = line_of(close)
        if not live + dead:
            continue
        parts = [part.strip() for part in macro.group("expr").split("##")]
        pattern = "".join(f"<{part}>" if part in params else part for part in parts)
        note = f"{live} cases at run time from {name}(...) at lines {first}-{last}, named {pattern}"
        if dead:
            note += f"; {dead} more compiled out"
        kind = "DISABLED" if disabled[define_start - 1] or not live else "MACRO"
        cases.append(Case(module, file, name, define_start, last, kind, note))
    cases.sort(key=lambda case: case.line_start)
    return cases


def _define_spans(lines: Sequence[str]) -> list[tuple[int, int, str]]:
    """Return (first line, last line, joined text) of each `#define`, with continuations."""

    spans: list[tuple[int, int, str]] = []
    index = 0
    while index < len(lines):
        if DEFINE_RE.match(lines[index]):
            start = index
            parts = [lines[index].rstrip()]
            while parts[-1].endswith("\\") and index + 1 < len(lines):
                parts[-1] = parts[-1][:-1]
                index += 1
                parts.append(lines[index].rstrip())
            spans.append((start + 1, index + 1, " ".join(parts)))
        index += 1
    return spans


def _matching_paren(text: str, opening: int) -> int:
    depth = 0
    for index in range(opening, len(text)):
        if text[index] == "(":
            depth += 1
        elif text[index] == ")":
            depth -= 1
            if depth == 0:
                return index
    raise LedgerError("unbalanced parentheses")


def extract_cases(ncbi_src: Path) -> list[Case]:
    src = ncbi_src / "c++" / "src"
    if not src.is_dir():
        raise LedgerError(f"not an NCBI BLAST+ source tree (no c++/src): {ncbi_src}")
    cases: list[Case] = []
    for module, directory, pattern in SOURCES:
        paths = sorted((src / directory).glob(pattern))
        if not paths:
            raise LedgerError(f"no NCBI unit-test source: c++/src/{directory}/{pattern}")
        for path in paths:
            file = f"c++/src/{directory}/{path.name}"
            text = path.read_text(encoding="utf-8", errors="replace")
            cases.extend(extract_cases_from_text(module, file, text))
    seen: dict[tuple[str, str], Case] = {}
    for case in cases:
        key = (case.module, case.case)
        if key in seen:
            first = seen[key]
            raise LedgerError(
                f"case name is not unique in module {case.module}: {case.case} at "
                f"{first.file}:{first.line_start} and {case.file}:{case.line_start}"
            )
        seen[key] = case
    return cases


def ncbi_commit(ncbi_src: Path) -> str:
    completed = subprocess.run(
        ["git", "-C", str(ncbi_src), "rev-parse", "HEAD"],
        check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    commit = completed.stdout.strip()
    if completed.returncode != 0 or not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise LedgerError(f"cannot read the NCBI source commit (git rev-parse HEAD): {ncbi_src}")
    return commit


def write_cases(path: Path, commit: str, cases: Iterable[Case]) -> None:
    lines = [
        "# NCBI BLAST+ unit-test cases in scope for LOSAT. Generated; do not edit.",
        "# python3 LOSAT/tests/ncbi_unit_case_ledger.py extract --ncbi-src \"$NCBI_SRC\" "
        f"--out {DEFAULT_CASES.as_posix()}",
        f"# {COMMIT_KEY}\t{commit}",
        "\t".join(CASES_HEADER),
    ]
    for case in cases:
        lines.append(
            "\t".join(
                (case.module, case.file, case.case, str(case.line_start), str(case.line_end), case.kind, case.note)
            )
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# --- reading the TSV files -------------------------------------------------


def _data_rows(path: Path, header: tuple[str, ...]) -> tuple[list[str], list[tuple[int, list[str]]]]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise LedgerError(f"cannot read {path}: {error}") from None
    comments = [line for line in lines if line.startswith("#")]
    rows = [(number, line) for number, line in enumerate(lines, start=1) if line.strip() and not line.startswith("#")]
    if not rows or tuple(rows[0][1].split("\t")) != header:
        raise LedgerError(f"{path}: header must be {chr(9).join(header)}")
    parsed: list[tuple[int, list[str]]] = []
    for number, line in rows[1:]:
        fields = line.split("\t")
        if len(fields) != len(header):
            raise LedgerError(f"{path}:{number}: expected {len(header)} tab-separated fields, found {len(fields)}")
        parsed.append((number, fields))
    return comments, parsed


def load_cases(path: Path) -> CaseList:
    comments, rows = _data_rows(path, CASES_HEADER)
    commit = ""
    for comment in comments:
        parts = comment.lstrip("# ").split("\t")
        if len(parts) == 2 and parts[0] == COMMIT_KEY:
            commit = parts[1].strip()
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise LedgerError(f"{path}: missing '# {COMMIT_KEY}<TAB><40-hex>' line")
    cases: list[Case] = []
    for number, (module, file, name, start, end, kind, note) in rows:
        if kind not in KINDS or not start.isdigit() or not end.isdigit():
            raise LedgerError(f"{path}:{number}: invalid case row")
        cases.append(Case(module, file, name, int(start), int(end), kind, note))
    return CaseList(commit, cases)


def load_ledger(path: Path) -> list[LedgerRow]:
    _, rows = _data_rows(path, LEDGER_HEADER)
    return [LedgerRow(*fields, row=number) for number, fields in rows]


# --- scan ------------------------------------------------------------------


def _line_kinds(lines: Sequence[str]) -> list[str]:
    """Classify Rust lines as comment, attr (first line of an attribute), attr+ (its
    continuation lines), blank, or code."""

    kinds: list[str] = []
    depth = 0
    for line in lines:
        if depth > 0:
            kinds.append("attr+")
            depth = max(0, depth + line.count("[") - line.count("]"))
        elif line.lstrip().startswith("//"):
            kinds.append("comment")
        elif ATTRIBUTE_RE.match(line):
            kinds.append("attr")
            depth = max(0, line.count("[") - line.count("]"))
        elif not line.strip():
            kinds.append("blank")
        else:
            kinds.append("code")
    return kinds


SKIPPED_KINDS = ("comment", "attr", "attr+")


@dataclass(frozen=True)
class FunctionInfo:
    name: str | None
    is_test: bool
    is_ignored: bool


def _function_below(lines: Sequence[str], kinds: Sequence[str], index: int) -> FunctionInfo:
    """The fn on the first code line at or after `index` (past comments and attributes),
    with the test and ignore attributes of the comment/attribute block around it."""

    attributes: list[str] = []
    back = index - 1
    while back >= 0 and kinds[back] in SKIPPED_KINDS:
        if kinds[back] == "attr":
            attributes.append(lines[back])
        back -= 1
    forward = index
    while forward < len(lines) and kinds[forward] in SKIPPED_KINDS:
        if kinds[forward] == "attr":
            attributes.append(lines[forward])
        forward += 1
    match = FN_RE.match(lines[forward]) if forward < len(lines) and kinds[forward] == "code" else None
    if not match:
        return FunctionInfo(None, False, False)
    return FunctionInfo(
        match.group(1),
        any(TEST_ATTRIBUTE_RE.match(attribute) for attribute in attributes),
        any(IGNORE_ATTRIBUTE_RE.match(attribute) for attribute in attributes),
    )


def _rust_files(root: Path) -> list[Path]:
    files: list[Path] = []
    for directory in SCAN_DIRS:
        base = root / directory
        if base.is_dir():
            files.extend(path for path in base.rglob("*.rs") if "target" not in path.relative_to(root).parts)
    return sorted(files)


def _parse_ranges(text: str) -> tuple[tuple[int, int], ...]:
    ranges = []
    for part in text.split(","):
        start, _, end = part.partition("-")
        ranges.append((int(start), int(end or start)))
    return tuple(ranges)


def scan_citations(root: Path) -> tuple[list[Citation], list[tuple[str, int, str]]]:
    """Return the well-formed citations and the malformed citation lines."""

    citations: list[Citation] = []
    malformed: list[tuple[str, int, str]] = []
    for path in _rust_files(root):
        lines = path.read_text(encoding="utf-8").splitlines()
        kinds: list[str] | None = None
        rel = path.relative_to(root).as_posix()
        for index, line in enumerate(lines):
            if not CITATION_MARK_RE.search(line):
                continue
            match = CITATION_RE.match(line)
            if not match:
                malformed.append((rel, index + 1, line.strip()))
                continue
            kinds = kinds or _line_kinds(lines)
            function = _function_below(lines, kinds, index + 1)
            citations.append(
                Citation(
                    rel, index + 1, match.group("commit"), match.group("file"),
                    _parse_ranges(match.group("lines")), match.group("case"),
                    function.name, function.is_test, function.is_ignored,
                )
            )
    return citations, malformed


class TestIndex:
    """The functions of each Rust file and their test attributes, read on demand."""

    def __init__(self, root: Path) -> None:
        self.root = root
        self._files: dict[str, dict[str, FunctionInfo] | None] = {}

    def functions(self, rel: str) -> dict[str, FunctionInfo] | None:
        if rel not in self._files:
            path = self.root / rel
            if not path.is_file():
                self._files[rel] = None
            else:
                lines = path.read_text(encoding="utf-8").splitlines()
                kinds = _line_kinds(lines)
                found: dict[str, FunctionInfo] = {}
                for index, line in enumerate(lines):
                    if kinds[index] != "code" or not FN_RE.match(line):
                        continue
                    info = _function_below(lines, kinds, index)
                    previous = found.get(info.name)
                    # A name defined twice (two test modules) counts as a running test if either is one.
                    if previous is None or (info.is_test and not info.is_ignored):
                        found[info.name] = info
                self._files[rel] = found
        return self._files[rel]

    def resolve(self, reference: str) -> str | None:
        """Return None when `reference` names an existing test or fixture, else the reason."""

        if "::" in reference:
            rel, _, name = reference.partition("::")
            functions = self.functions(rel)
            if functions is None:
                return f"losat_test file not found: {rel}"
            if name not in functions:
                return f"losat_test fn not found: {reference}"
            if not functions[name].is_test:
                return f"losat_test fn is not a #[test]: {reference}"
            if functions[name].is_ignored:
                return f"losat_test fn is #[ignore]d: {reference}"
            return None
        rel, _, anchor = reference.partition("#")
        if rel.endswith(".rs"):
            return f"losat_test names a Rust file without ::fn: {reference}"
        path = self.root / rel
        if not path.is_file():
            return f"losat_test fixture file not found: {rel}"
        if anchor and anchor not in path.read_text(encoding="utf-8", errors="replace"):
            return f"losat_test fixture id not found: {reference}"
        return None


# --- check -----------------------------------------------------------------


def check(root: Path, cases: CaseList, ledger: Sequence[LedgerRow]) -> tuple[list[Failure], dict[str, int]]:
    failures: list[Failure] = []
    citations, malformed = scan_citations(root)
    index = TestIndex(root)

    for rel, line, text in malformed:
        failures.append(Failure("-", "-", f"malformed citation at {rel}:{line}: {text}"))

    # (1) Each citation names an NCBI case, within its lines, above a #[test] fn.
    cited: dict[tuple[str, str], list[Citation]] = {}
    for citation in citations:
        where = f"{citation.path}:{citation.line}"
        case = cases.by_file_case.get((citation.file, citation.case))
        if case is None:
            failures.append(Failure("-", citation.case, f"cited case not in NCBI_CASES.tsv for {citation.file} ({where})"))
            continue
        if not (cases.commit.startswith(citation.commit) or citation.commit.startswith(cases.commit)):
            failures.append(Failure(case.module, case.case, f"citation commit {citation.commit} is not NCBI_CASES.tsv's {cases.commit[:8]} ({where})"))
        outside = [f"{a}-{b}" for a, b in citation.ranges if a > b or a < case.line_start or b > case.line_end]
        if outside:
            failures.append(Failure(case.module, case.case, f"cited lines {','.join(outside)} outside the case ({case.line_start}-{case.line_end}) ({where})"))
        if citation.function is None:
            failures.append(Failure(case.module, case.case, f"citation is not directly above a fn ({where})"))
        elif not citation.is_test:
            failures.append(Failure(case.module, case.case, f"citation is above fn {citation.function}, which is not a #[test] ({where})"))
        elif citation.is_ignored:
            failures.append(Failure(case.module, case.case, f"citation is above fn {citation.function}, which is #[ignore]d ({where})"))
        else:
            cited.setdefault((case.module, case.case), []).append(citation)

    # (4) Ledger rows are well formed, name NCBI cases, and are unique.
    rows: dict[tuple[str, str], LedgerRow] = {}
    for row in ledger:
        key = (row.module, row.case)
        if key in rows:
            failures.append(Failure(row.module, row.case, f"duplicate ledger row (rows {rows[key].row} and {row.row})"))
            continue
        rows[key] = row
        if key not in cases.by_module_case:
            failures.append(Failure(row.module, row.case, f"ledger row {row.row} names a case not in NCBI_CASES.tsv"))
        if row.status not in STATUSES:
            failures.append(Failure(row.module, row.case, f"status '{row.status}' is not one of {', '.join(STATUSES)}"))
        if row.status == "superseded" and not row.note.strip():
            failures.append(Failure(row.module, row.case, "superseded row needs the reason in note"))
        if row.status == "partial" and not row.tests and not row.note.strip() and key not in cited:
            failures.append(Failure(row.module, row.case, "partial row has no test, citation or evidence note"))
        if not DATE_RE.match(row.since):
            failures.append(Failure(row.module, row.case, f"since '{row.since}' is not YYYY-MM-DD"))

    # (2) Every losat_test named by a row exists; ported rows are linked to a test.
    linked_in_ledger: dict[str, int] = {}
    linked_in_code: dict[str, int] = {}
    for key, row in rows.items():
        reasons = [reason for reason in (index.resolve(test) for test in row.tests) if reason]
        for reason in reasons:
            failures.append(Failure(row.module, row.case, reason))
        if row.status not in LINKED_STATUSES:
            continue
        # A ported or partial case is linked by a citation or a `path::fn` test; a
        # fixture alone is e2e coverage.
        tests = [test for test in row.tests if "::" in test]
        if not (tests or key in cited):
            if row.tests:
                failures.append(Failure(row.module, row.case, f"{row.status} row names no #[test] fn and no citation (fixtures alone are e2e)"))
            elif row.status == "ported":
                failures.append(Failure(row.module, row.case, "ported row has no citation line and no losat_test"))
            continue
        linked_in_ledger[row.module] = linked_in_ledger.get(row.module, 0) + 1
        if key in cited or (tests and not reasons):
            linked_in_code[row.module] = linked_in_code.get(row.module, 0) + 1

    # (3) Each cited case has a ledger row that counts it as ported or partial.
    for (module, name), found in sorted(cited.items()):
        row = rows.get((module, name))
        where = ", ".join(f"{c.path}:{c.line}" for c in found)
        if row is None:
            failures.append(Failure(module, name, f"cited in code but has no ledger row ({where})"))
        elif row.status not in LINKED_STATUSES:
            failures.append(Failure(module, name, f"cited in code but ledger status is {row.status} ({where})"))

    # (5) Ratchet: the ledger may not count more linked ported/partial cases
    # than the code backs (a test was removed without fixing the ledger).
    for module in sorted(linked_in_ledger):
        ledger_count = linked_in_ledger[module]
        code_count = linked_in_code.get(module, 0)
        if ledger_count > code_count:
            failures.append(Failure(module, "-", f"ledger counts {ledger_count} linked ported/partial cases, code backs {code_count}"))

    summary = {
        "citations": len(citations),
        "cited_cases": len(cited),
        "ledger_rows": len(rows),
        "linked_ported_partial": sum(linked_in_ledger.values()),
        "backed_by_code": sum(linked_in_code.values()),
    }
    return failures, summary


# --- report ----------------------------------------------------------------


def report(cases: CaseList, ledger: Sequence[LedgerRow]) -> str:
    modules: list[str] = []
    totals: dict[str, int] = {}
    for case in cases.cases:
        if case.module not in modules:
            modules.append(case.module)
        totals[case.module] = totals.get(case.module, 0) + 1
    counts: dict[str, dict[str, int]] = {module: {} for module in modules}
    unlinked = 0
    for row in ledger:
        if (row.module, row.case) not in cases.by_module_case or row.status not in STATUSES:
            continue
        counts[row.module][row.status] = counts[row.module].get(row.status, 0) + 1
        if row.status == "partial" and not row.tests:
            unlinked += 1
    columns = (*STATUSES, "unlisted")
    lines = [
        "| NCBI module | cases | " + " | ".join(columns) + " |",
        "|---|---:|" + "---:|" * len(columns),
    ]
    grand = {column: 0 for column in columns}
    for module in modules:
        listed = sum(counts[module].values())
        values = [counts[module].get(status, 0) for status in STATUSES] + [totals[module] - listed]
        for column, value in zip(columns, values):
            grand[column] += value
        lines.append(f"| {module} | {totals[module]} | " + " | ".join(str(v) for v in values) + " |")
    lines.append(
        f"| **total** | **{len(cases.cases)}** | " + " | ".join(f"**{grand[c]}**" for c in columns) + " |"
    )
    lines.append("")
    lines.append(
        f"NCBI {cases.commit[:8]}. Partial rows linked only by an evidence note (no citation, no losat_test): {unlinked}."
    )
    return "\n".join(lines)


# --- command line ----------------------------------------------------------


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    extract_parser = commands.add_parser("extract", help="list the NCBI unit-test cases (needs the NCBI source)")
    extract_parser.add_argument("--ncbi-src", type=Path, required=True)
    extract_parser.add_argument("--out", type=Path, required=True)
    for name, text in (
        ("scan", "list the citation lines in LOSAT/src and LOSAT/tests"),
        ("check", "check citations, ledger and code against each other"),
        ("report", "Markdown table of ledger statuses per NCBI module"),
    ):
        command = commands.add_parser(name, help=text)
        command.add_argument("--root", type=Path, default=default_root)
        if name != "scan":
            command.add_argument("--cases", type=Path)
            command.add_argument("--ledger", type=Path)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        if args.command == "extract":
            cases = extract_cases(args.ncbi_src)
            write_cases(args.out, ncbi_commit(args.ncbi_src), cases)
            kinds = {kind: sum(1 for case in cases if case.kind == kind) for kind in KINDS}
            print(f"wrote {len(cases)} cases ({', '.join(f'{k} {v}' for k, v in kinds.items())}) to {args.out}")
            return 0
        root = args.root.resolve()
        if args.command == "scan":
            citations, malformed = scan_citations(root)
            for citation in citations:
                lines = ",".join(f"{a}-{b}" if a != b else str(a) for a, b in citation.ranges)
                function = citation.function or "-"
                test = "test" if citation.is_test else "not-a-test"
                print(f"{citation.path}:{citation.line}\t{citation.file}:{lines}\t{citation.case}\t{function}\t{test}")
            for rel, line, text in malformed:
                print(f"MALFORMED\t{rel}:{line}\t{text}")
            return 1 if malformed else 0
        cases = load_cases(args.cases or root / DEFAULT_CASES)
        ledger = load_ledger(args.ledger or root / DEFAULT_LEDGER)
        if args.command == "report":
            print(report(cases, ledger))
            return 0
        failures, summary = check(root, cases, ledger)
    except LedgerError as error:
        print(f"ncbi unit-case ledger error: {error}", file=sys.stderr)
        return 2
    for failure in failures:
        print(f"FAIL\t{failure.module}\t{failure.case}\t{failure.reason}")
    print("summary\t" + "\t".join(f"{key}={value}" for key, value in summary.items()) + f"\tfailures={len(failures)}")
    print(f"result\t{'FAIL' if failures else 'PASS'}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
