"""LOSAT's text of NCBI's frozen bytes: the one mapping of `LOSAT/src/wording.rs` (AGENTS.md rule 8).

LOSAT's messages say "Please refer to the LOSAT user manual." and leave out the internals of NCBI's
exception reports (`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exception 12). A freeze keeps NCBI's bytes;
a harness that compares LOSAT's standard error (or an archive's messages) with a freeze compares
it with `losat_bytes(freeze)`. The constants are read from `wording.rs`, and the rules below are
that file's, step for step; `LOSAT/tests/fixtures/wording/samples.tsv` holds the samples that both
check (`check_samples`, run on import).
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORDING_RS = ROOT / "LOSAT" / "src" / "wording.rs"
SAMPLES = ROOT / "LOSAT" / "tests" / "fixtures" / "wording" / "samples.tsv"


def _constants() -> dict[str, str]:
    text = WORDING_RS.read_text(encoding="utf-8")
    found = dict(re.findall(r'^pub const ([A-Z_]+): &str = "([^"\\]*)";$', text, re.M))
    serial = re.search(r"^pub const SERIAL_EXCEPTION_CLASSES: &\[&str\] = &\[([^\]]*)\];$", text, re.M)
    names = ["USER_MANUAL_FOOTER", "NCBI_USER_MANUAL_FOOTER", "NCBI_EXCEPTION_HEADER", "NCBI_POSITION", "POSITION"]
    missing = [name for name in names if name not in found]
    if missing or serial is None:
        raise RuntimeError(f"{WORDING_RS}: cannot read {missing or ['SERIAL_EXCEPTION_CLASSES']}")
    found["SERIAL_EXCEPTION_CLASSES"] = ",".join(re.findall(r'"([A-Za-z_]+)"', serial.group(1)))
    return found


_C = _constants()
USER_MANUAL_FOOTER = _C["USER_MANUAL_FOOTER"].encode()
NCBI_USER_MANUAL_FOOTER = _C["NCBI_USER_MANUAL_FOOTER"].encode()
_HEADER = _C["NCBI_EXCEPTION_HEADER"].encode()
_NCBI_POSITION = _C["NCBI_POSITION"].encode()
_POSITION = _C["POSITION"].encode()
_SERIAL = {name.encode() for name in _C["SERIAL_EXCEPTION_CLASSES"].split(",")}


def _digits(text: bytes, i: int) -> int | None:
    j = i
    while j < len(text) and text[j : j + 1].isdigit():
        j += 1
    return j if j > i else None


def _literal(text: bytes, i: int, expected: bytes) -> int | None:
    return i + len(expected) if text.startswith(expected, i) else None


def _tag(text: bytes, i: int) -> int | None:
    """The end of `(<class>::<code>) ` at `i` and of the function after it (`<function>() - `)."""
    j = _literal(text, i, b"(")
    if j is None:
        return None
    k = j
    while k < len(text) and (text[k : k + 1].isalnum() and text[k] < 128 or text[k : k + 1] == b"_"):
        k += 1
    name = text[j:k]
    if not (name.endswith(b"Exception") or name in _SERIAL):
        return None
    k = _literal(text, k, b"::")
    if k is None:
        return None
    close = [p for p in (text.find(b")", k), text.find(b"\n", k)) if p >= 0]
    if not close:
        return None
    k = _literal(text, min(close), b") ")
    if k is None:
        return None
    f = k
    while f < len(text) and (text[f : f + 1].isalnum() and text[f] < 128 or text[f : f + 1] in (b"_", b":", b"~")):
        f += 1
    if f > k:
        end = _literal(text, f, b"() - ")
        if end is not None:
            k = end
    return k


def _place(text: bytes, i: int) -> int | None:
    """The end of `T<thread> "<file>", line <n>: <severity>: ` at `i` and of its tag."""
    j = _digits(text, i + 1)
    j = j and _literal(text, j, b' "')
    if j is None:
        return None
    close = [p for p in (text.find(b'"', j), text.find(b"\n", j)) if p >= 0]
    if not close:
        return None
    j = _literal(text, min(close), b'", line ')
    j = j and _digits(text, j)
    j = j and _literal(text, j, b": ")
    if j is None:
        return None
    w = j
    while w < len(text) and text[w : w + 1].isalpha() and text[w] < 128:
        w += 1
    if w == j:
        return None
    j = _literal(text, w, b": ")
    return None if j is None else _tag(text, j)


def _glue(out: bytearray) -> None:
    """`wording::glue`: `: ` in place of a removed report or tag inside a line."""
    kept = out.rstrip(b" ")
    if kept and kept[-1:] not in (b":", b"\n"):
        if out.endswith(b" "):
            out.pop()
        out += b": "


def losat_bytes(text: bytes) -> bytes:
    """LOSAT's text of NCBI's `text` (`wording::losat_bytes`)."""
    if not any(needle in text for needle in (NCBI_USER_MANUAL_FOOTER, _HEADER, _NCBI_POSITION, b"::")):
        return text
    out = bytearray()
    i = 0
    while i < len(text):
        if text.startswith(NCBI_USER_MANUAL_FOOTER, i):
            out += USER_MANUAL_FOOTER
            i += len(NCBI_USER_MANUAL_FOOTER)
            continue
        if text.startswith(_HEADER, i):
            j = i + len(_HEADER)
            while j < len(text) and text[j] in b" \n":
                j += 1
            if j > i + len(_HEADER):
                end = _place(text, j) if text.startswith(b"T", j) else None
                if end is not None:
                    if out and out[-1:] not in (b" ", b"\n", b":"):
                        _glue(out)
                    i = end
                    continue
        if text.startswith(b"T", i) and out.endswith(b" "):
            end = _place(text, i)
            if end is not None:
                i = end
                continue
        if text.startswith(b"(", i):
            end = _tag(text, i)
            if end is not None:
                _glue(out)
                i = end
                continue
        if text.startswith(_NCBI_POSITION, i):
            out += _POSITION
            i += len(_NCBI_POSITION)
            continue
        out.append(text[i])
        i += 1
    return bytes(out)


_ASN_TEXT_RS = ROOT / "LOSAT" / "src" / "report" / "serial" / "asn_text.rs"
_STRING_WRAP = int(re.search(r"^const STRING_WRAP: usize = (\d+);$", _ASN_TEXT_RS.read_text(encoding="utf-8"), re.M).group(1))
_ARCHIVE_MESSAGE = b'message "'


def visible_string_line(line: bytes, value: bytes) -> bytes:
    """`asn_text::visible_string_line`: `line` and the VisibleString `value` after it, wrapped as
    LOSAT's ASN.1 text writer wraps them (`write_string`, `put_eol_at_word_end`)."""
    out = bytearray(b"\n" + line)
    line_length = len(line)

    def put(data: bytes) -> None:
        nonlocal line_length
        out.extend(data)
        line_length += len(data)

    def eol_at_word_end(limit: int) -> None:
        nonlocal line_length
        line_pos = line_length
        pos = len(out)
        good_place = False
        while pos > 0 and line_pos > 0:
            pos -= 1
            line_pos -= 1
            byte = out[pos]
            if line_pos <= limit and byte in b" \t\n\x0b\x0c\r'":
                good_place = True
                break
            if byte in b'\n"':
                break
        if pos > 0 and out[pos - 1] == 0x0A:
            good_place = False
        if not good_place:
            if line_pos < limit:
                pos += limit - line_pos
            while pos > 0 and out[pos - 1] == 0x22:
                pos -= 1
            if pos == 0:
                while pos < len(out) and out[pos] == 0x22:
                    pos += 1
        count = len(out) - pos
        out.insert(pos, 0x0A)
        line_length = count

    put(b'"')
    for byte in value:
        if line_length > _STRING_WRAP:
            eol_at_word_end(_STRING_WRAP)
        put(bytes([byte]))
        if byte == 0x22:
            put(b'"')
    put(b'"')
    return bytes(out[1:])


def _asn_string(text: bytes, quote: int) -> tuple[bytes, int] | None:
    value = bytearray()
    i = quote + 1
    while i < len(text):
        byte = text[i]
        if byte == 0x22 and text[i + 1 : i + 2] == b'"':
            value += b'"'
            i += 2
        elif byte == 0x22:
            return bytes(value), i + 1
        elif byte == 0x0A:
            i += 1
        else:
            value.append(byte)
            i += 1
    return None


def losat_output_bytes(output: bytes) -> bytes:
    """`wording::losat_output_bytes`: a whole output (standard output, an `-out` file, or both
    merged), whose ASN.1 text archive messages are mapped and wrapped again."""
    if _ARCHIVE_MESSAGE not in output:
        return losat_bytes(output)
    out = bytearray()
    done = 0
    search = 0
    while True:
        at = output.find(_ARCHIVE_MESSAGE, search)
        if at < 0:
            break
        line_start = output.rfind(b"\n", 0, at) + 1
        quote = at + len(_ARCHIVE_MESSAGE) - 1
        search = quote + 1
        if output[line_start:at].strip(b" "):
            continue
        parsed = _asn_string(output, quote)
        if parsed is None:
            continue
        value, end = parsed
        search = end
        mapped = losat_bytes(value)
        if mapped != value:
            out += losat_bytes(output[done:line_start])
            out += visible_string_line(output[line_start:quote], mapped)
            done = end
    if done == 0:
        return losat_bytes(output)
    out += losat_bytes(output[done:])
    return bytes(out)


def losat_text(text: str) -> str:
    return losat_bytes(text.encode()).decode()


def check_samples() -> int:
    count = 0
    for line in SAMPLES.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        ncbi, losat = (field.replace("\\n", "\n").replace("\\t", "\t") for field in line.split("\t"))
        if losat_text(ncbi) != losat or losat_text(losat) != losat:
            raise RuntimeError(f"{SAMPLES}: {ncbi!r} maps to {losat_text(ncbi)!r}, not {losat!r}")
        count += 1
    return count


def check_archive_wraps(paths) -> int:
    """`visible_string_line` writes each message of these frozen archives as NCBI wrote it."""
    count = 0
    for path in paths:
        data = Path(path).read_bytes()
        search = 0
        while (at := data.find(_ARCHIVE_MESSAGE, search)) >= 0:
            line_start = data.rfind(b"\n", 0, at) + 1
            quote = at + len(_ARCHIVE_MESSAGE) - 1
            search = quote + 1
            if data[line_start:at].strip(b" ") or (parsed := _asn_string(data, quote)) is None:
                continue
            value, end = parsed
            search = end
            if visible_string_line(data[line_start:quote], value) != data[line_start:end]:
                raise RuntimeError(f"{path}: the message at byte {at} is wrapped otherwise")
            count += 1
    return count


check_samples()
