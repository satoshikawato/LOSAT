"""LOSAT's text of NCBI's frozen bytes: the one mapping of `LOSAT/src/wording.rs` (AGENTS.md rule 8).

LOSAT states its own name and version (`BLASTP (LOSAT 0.2.0)` where NCBI writes `BLASTP 2.17.0+`),
writes `LOSAT engine error:` and the other prefixes where NCBI writes `BLAST engine error:`, says
"Please refer to the LOSAT user manual." and leaves out the internals of NCBI's exception reports
(`PD-LOSAT-CLI-NONSEARCH-DIFFERENCES` exceptions 12 and 14). A freeze keeps NCBI's bytes; a harness
that compares LOSAT's output or standard error with a freeze compares it with
`losat_output_bytes(freeze)` (`losat_bytes` for a message). The constants are read from
`wording.rs`, LOSAT's version from `LOSAT/Cargo.toml` and NCBI's from `LOSAT/src/lib.rs`; the rules
below are `wording.rs`'s, step for step; `LOSAT/tests/fixtures/wording/samples.tsv` holds the samples
that both check (`check_samples`, run on import).

On the command line (ad-hoc probes): `python3 -I LOSAT/tests/losat_wording.py map-output <file>`
writes LOSAT's bytes of NCBI's output file to standard output (`map-message` for a message).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
WORDING_RS = ROOT / "LOSAT" / "src" / "wording.rs"
LIB_RS = ROOT / "LOSAT" / "src" / "lib.rs"
CARGO_TOML = ROOT / "LOSAT" / "Cargo.toml"
SAMPLES = ROOT / "LOSAT" / "tests" / "fixtures" / "wording" / "samples.tsv"
OUTPUT_SAMPLES = SAMPLES.with_name("output_samples.tsv")


def _versions() -> tuple[str, str]:
    """LOSAT's version (`CARGO_PKG_VERSION`) and NCBI's release (`ncbi_blast_release!()`)."""
    cargo = re.search(r'^\[package\][^\[]*?^version = "([^"]+)"$', CARGO_TOML.read_text(encoding="utf-8"), re.M | re.S)
    release = re.search(r'macro_rules! ncbi_blast_release \{\s*\(\) => \{\s*"([^"]+)"', LIB_RS.read_text(encoding="utf-8"))
    if cargo is None or release is None:
        raise RuntimeError(f"cannot read the versions from {CARGO_TOML} and {LIB_RS}")
    return cargo.group(1), release.group(1)


_PREFIX_NAMES = ("QUERY_OPTIONS", "ENGINE", "DATABASE", "DATABASE_CREATION", "OPTIONS", "QUERY")


def _constants() -> dict[str, str]:
    text = WORDING_RS.read_text(encoding="utf-8")
    found = dict(re.findall(r'^pub const ([A-Z_]+): &str = "([^"\\]*)";$', text, re.M))
    serial = re.search(r"^pub const SERIAL_EXCEPTION_CLASSES: &\[&str\] = &\[([^\]]*)\];$", text, re.M)
    names = ["USER_MANUAL_FOOTER", "NCBI_USER_MANUAL_FOOTER", "NCBI_EXCEPTION_HEADER", "NCBI_POSITION", "POSITION",
             "RPSBPROC_TITLE", "NCBI_RPSBPROC_VERSION"]
    names += [f"{ncbi}{name}_ERROR" for name in _PREFIX_NAMES for ncbi in ("", "NCBI_")]
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
LOSAT_CRATE_VERSION, NCBI_RELEASE = _versions()
# `losat_version!()`, `report_version!()` and `NCBI_VERSION` (`crate::NCBI_BLAST_VERSION`) of wording.rs.
LOSAT_VERSION = f"LOSAT {LOSAT_CRATE_VERSION}".encode()
VERSION_NUMBER = LOSAT_CRATE_VERSION.encode()  # `wording::VERSION_NUMBER`
REPORT_VERSION = b"(" + LOSAT_VERSION + b")"
NCBI_VERSION = f"{NCBI_RELEASE}+".encode()
ERROR_PREFIXES = [(_C[f"NCBI_{name}_ERROR"].encode(), _C[f"{name}_ERROR"].encode()) for name in _PREFIX_NAMES]
QUERY_ERROR = _C["QUERY_ERROR"].encode()  # `wording::QUERY_ERROR`, the query reader's prefix
_RPSBPROC_TITLE = _C["RPSBPROC_TITLE"].encode()
_NCBI_RPSBPROC_BANNER = _RPSBPROC_TITLE + b" v" + _C["NCBI_RPSBPROC_VERSION"].encode()


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


_BER_VISIBLE_STRING = 0x1A  # `wording::BER_VISIBLE_STRING`: a binary archive's VisibleString tag
IDENT_GROWTH = len(LOSAT_VERSION) - len(NCBI_VERSION)  # `wording::IDENT_GROWTH`
_BYTE_POSITION = b" in string at byte "  # `wording::BYTE_POSITION`
# `wording::archive_byte_post`: a byte's post of a binary archive, whose trace names the archive.
_ARCHIVE_BYTE_POST = re.compile(rb"(Blast4-(?:archive|request)[A-Za-z0-9._-]*\nBad char \[0x[0-9A-F]+\]) in string at byte (\d+)")


def _after_program_name(out: bytearray) -> bool:
    """`wording::after_program_name`: `out` ends with a program's name in capitals and a space."""
    if not out.endswith(b" "):
        return False
    j = len(out) - 1
    while j > 0 and (chr(out[j - 1]).isascii() and (chr(out[j - 1]).isalnum() or out[j - 1] in b"_-")):
        j -= 1
    word = bytes(out[j : len(out) - 1])
    return any(65 <= b <= 90 for b in word) and not any(97 <= b <= 122 for b in word)


_XMLHELP_PROGRAM = b'<program type="'  # `wording::XMLHELP_PROGRAM`
_XMLHELP_VERSION = re.compile(rb"<version>([0-9.]+)</version>")


_XMLHELP_HEAD = re.compile(rb'(?:\A|[\n ])<program type="[^\n<>]*">[\n ]<name>[^<\n]*</name>[\n ]\Z')


def _xmlhelp_version(out: bytearray, text: bytes, i: int) -> int | None:
    """`wording::xmlhelp_version`: the end of `<version>V</version>` at `i` when `out` ends with the
    lines `<program type="...">` and `<name>...</name>` (each line feed a space in a collapsed message)."""
    if not (out.endswith(b"</name>\n") or out.endswith(b"</name> ")) or not _XMLHELP_HEAD.search(bytes(out)):
        return None
    m = _XMLHELP_VERSION.match(text, i)
    return m.end() if m else None


_BAD_CHAR_LINE = re.compile(rb"Bad char \[0x([0-9A-F]{1,2})\] in string at (?:line|byte) [0-9]+[\n ]")


def _quoted_piece(text: bytes, i: int) -> int | None:
    """The end of a byte post's line and of the piece it quotes at `i`, when the piece (up to the next
    line end) holds the byte (`wording::quoted_piece`), or None."""
    m = _BAD_CHAR_LINE.match(text, i)
    if not m:
        return None
    end = text.find(b"\n", m.end())
    end = len(text) if end < 0 else end
    return end if int(m.group(1), 16) in text[m.end():end] else None


def losat_bytes(text: bytes) -> bytes:
    """LOSAT's text of NCBI's `text` (`wording::losat_bytes`)."""
    needles = (NCBI_USER_MANUAL_FOOTER, _HEADER, _NCBI_POSITION, b"::", NCBI_VERSION, b" error:", _RPSBPROC_TITLE,
               _XMLHELP_PROGRAM)
    if not any(needle in text for needle in needles):
        return text
    out = bytearray()
    i = 0
    while i < len(text):
        if text.startswith(b"B", i) and (not out or out[-1:] in (b"\n", b" ")):
            # The piece that a byte's post quotes (the input's bytes) is copied as it is.
            end = _quoted_piece(text, i)
            if end is not None:
                out += text[i:end]
                i = end
                continue
        if text.startswith(NCBI_VERSION, i) and out.endswith(bytes([_BER_VISIBLE_STRING, len(NCBI_VERSION)])):
            out[-1] = len(LOSAT_VERSION)
            out += LOSAT_VERSION
            i += len(NCBI_VERSION)
            continue
        if text.startswith(NCBI_VERSION, i) and (not out or out[-1:] in (b"\n", b"\t", b'"', b">", b":", b" ")):
            out += REPORT_VERSION if _after_program_name(out) else LOSAT_VERSION
            i += len(NCBI_VERSION)
            continue
        prefix = next((pair for pair in ERROR_PREFIXES if text.startswith(pair[0], i)), None)
        if prefix is not None:
            out += prefix[1]
            i += len(prefix[0])
            continue
        if text.startswith(_NCBI_RPSBPROC_BANNER, i):
            out += _RPSBPROC_TITLE + b" " + REPORT_VERSION
            i += len(_NCBI_RPSBPROC_BANNER)
            continue
        if text.startswith(b"<", i):
            end = _xmlhelp_version(out, text, i)
            if end is not None:
                out += b"<version>" + VERSION_NUMBER + b"</version>"
                i = end
                continue
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
    merged), whose ASN.1 text archive messages are mapped and wrapped again, and the position of each
    byte's post of a binary archive moved by IDENT_GROWTH (`wording::archive_positions`)."""
    words = _output_words(output)
    if _BYTE_POSITION not in words:
        return words
    return _ARCHIVE_BYTE_POST.sub(lambda m: m.group(1) + _BYTE_POSITION + str(int(m.group(2)) + IDENT_GROWTH).encode(), words)


def _output_words(output: bytes) -> bytes:
    """`wording::output_words`: `losat_output_bytes` without the positions of the byte posts."""
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


_XMLHELP_DIR = ROOT / "LOSAT" / "src" / "ncbi_xmlhelp"
_XMLHELP_LOSAT = re.compile(rb'((?:\A|(?<=[\n ]))<program type="[^\n<>]*">[\n ]<name>([^<\n]*)</name>[\n ]<version>)'
                            + re.escape(VERSION_NUMBER)
                            + rb"(</version>)")


def _ncbi_xmlhelp_version(name: bytes) -> bytes:
    """NCBI's version of `-xmlhelp`'s program `name`, from its frozen XML (LOSAT's when there is none)."""
    path = _XMLHELP_DIR / f"{name.decode(errors='replace')}.xml"
    m = _XMLHELP_VERSION.search(path.read_bytes()) if path.is_file() else None
    return m.group(1) if m else VERSION_NUMBER


def ncbi_hash_form(losat: bytes, positions: bool = True) -> bytes:
    """LOSAT's output read back in NCBI's words, for a harness that holds only the hash of NCBI's bytes
    (no freeze to map with `losat_output_bytes`): the version statement, the error prefixes, rpsbproc's
    banner, the user-manual footer and `-xmlhelp`'s program version (NCBI's, from the program's frozen
    `LOSAT/src/ncbi_xmlhelp/<name>.xml`) become NCBI's again. Exact for those rules (each LOSAT text stands
    for one NCBI text); an exception report cannot be read back, so a harness with its freeze compares
    `losat_output_bytes(freeze)` instead. Pair it with `names_ncbi` so that an output that still writes
    NCBI's words is not taken for LOSAT's. Without `positions`, the positions of a binary archive's byte
    posts stay LOSAT's (`_output_words`, a run whose shorter title takes up the bytes the version adds)."""
    out = losat.replace(_RPSBPROC_TITLE + b" " + REPORT_VERSION, _NCBI_RPSBPROC_BANNER)
    out = _XMLHELP_LOSAT.sub(lambda m: m.group(1) + _ncbi_xmlhelp_version(m.group(2)) + m.group(3), out)
    if positions:
        out = _ARCHIVE_BYTE_POST.sub(lambda m: m.group(1) + _BYTE_POSITION + str(int(m.group(2)) - IDENT_GROWTH).encode(), out)
    out = out.replace(bytes([_BER_VISIBLE_STRING, len(LOSAT_VERSION)]) + LOSAT_VERSION,
                      bytes([_BER_VISIBLE_STRING, len(NCBI_VERSION)]) + NCBI_VERSION)
    out = out.replace(REPORT_VERSION, NCBI_VERSION).replace(LOSAT_VERSION, NCBI_VERSION)
    for ncbi, ours in ERROR_PREFIXES:
        out = out.replace(ours, ncbi)
    return out.replace(USER_MANUAL_FOOTER, NCBI_USER_MANUAL_FOOTER)


def names_ncbi(output: bytes) -> bytes | None:
    """The first of NCBI's words that LOSAT no longer writes (its version, a `BLAST ... error:` prefix,
    the BLAST+ user-manual footer) found in LOSAT's `output`, or None."""
    for needle in (NCBI_VERSION, NCBI_USER_MANUAL_FOOTER, *(ncbi for ncbi, _ in ERROR_PREFIXES)):
        at = output.find(needle)
        # NCBI's package names (`ncbi-blast-2.17.0+-x64-linux.tar.gz`) are not version statements.
        if at >= 0 and not (needle == NCBI_VERSION and output[max(0, at - 11) : at] == b"ncbi-blast-"):
            return needle
    return None


def unescape_sample(field: str) -> str:
    """A field of samples.tsv: `\\n` and `\\t` are a line feed and a tab, `\\xNN` the byte NN (below
    0x80), `{version}` is LOSAT's crate version, `{verlen}` the length byte of `LOSAT <version>` and
    `{at<n>}` the number n + IDENT_GROWTH (`wording::tests::unescape`)."""
    field = field.replace("\\n", "\n").replace("\\t", "\t").replace("{version}", LOSAT_CRATE_VERSION)
    field = field.replace("{verlen}", chr(len(LOSAT_VERSION)))
    field = re.sub(r"\{at(\d+)\}", lambda m: str(int(m.group(1)) + IDENT_GROWTH), field)
    return re.sub(r"\\x([0-7][0-9a-fA-F])", lambda m: chr(int(m.group(1), 16)), field)


def check_samples() -> int:
    count = 0
    for line in SAMPLES.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        ncbi, losat = (unescape_sample(field) for field in line.split("\t"))
        if losat_text(ncbi) != losat or losat_text(losat) != losat:
            raise RuntimeError(f"{SAMPLES}: {ncbi!r} maps to {losat_text(ncbi)!r}, not {losat!r}")
        # Without an exception report, LOSAT's text reads back as NCBI's (`ncbi_hash_form`).
        plain = not any(needle in ncbi.encode() for needle in (_HEADER, _NCBI_POSITION, b"::"))
        if plain and ncbi_hash_form(losat.encode()) != ncbi.encode():
            raise RuntimeError(f"{SAMPLES}: {losat!r} reads back as {ncbi_hash_form(losat.encode())!r}, not {ncbi!r}")
        count += 1
    return count


def check_output_samples() -> int:
    """output_samples.tsv: NCBI's frozen output, then LOSAT's (`wording::tests::the_output_samples_map_to_losats_bytes`)."""
    count = 0
    for line in OUTPUT_SAMPLES.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        ncbi, losat = (unescape_sample(field).encode("latin-1") for field in line.split("\t"))
        if losat_output_bytes(ncbi) != losat:
            raise RuntimeError(f"{OUTPUT_SAMPLES}: {ncbi!r} maps to {losat_output_bytes(ncbi)!r}, not {losat!r}")
        if ncbi_hash_form(losat) != ncbi:
            raise RuntimeError(f"{OUTPUT_SAMPLES}: {losat!r} reads back as {ncbi_hash_form(losat)!r}, not {ncbi!r}")
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
check_output_samples()


def _main(argv: list[str]) -> int:
    usage = ("usage: losat_wording.py map-output|map-message <ncbi-file>|-  (LOSAT's bytes of NCBI's, to standard output)\n"
             "       losat_wording.py hash-form <losat-file>|-  (LOSAT's output read back in NCBI's words, `ncbi_hash_form`;\n"
             "       exit 1, with nothing written, when it still writes NCBI's words, `names_ncbi`)")
    if len(argv) != 2 or argv[0] not in ("map-output", "map-message", "hash-form"):
        print(usage, file=sys.stderr)
        return 2
    data = sys.stdin.buffer.read() if argv[1] == "-" else Path(argv[1]).read_bytes()
    if argv[0] == "hash-form":
        found = names_ncbi(data)
        if found is not None:
            print(f"the output still writes NCBI's {found.decode()!r}", file=sys.stderr)
            return 1
        sys.stdout.buffer.write(ncbi_hash_form(data))
        return 0
    mapped = losat_output_bytes(data) if argv[0] == "map-output" else losat_bytes(data)
    sys.stdout.buffer.write(mapped)
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
