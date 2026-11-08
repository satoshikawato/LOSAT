#!/usr/bin/env python3
"""The one path prefix in frozen outputs, and the run root that stands for it.

The TBLASTN Stage G matrix (docs/evidence/tlosan_stage_g/matrix_162.jsonl) froze NCBI's bytes
of runs whose -query and -subject were absolute paths under FROZEN_PREFIX, the developer
machine's checkout. Outfmt 7 prints the -subject string in its `# Database:` line, and outfmt 0
in the `Database: User specified sequence set (Input: <path>).` lines, wrapped at 68 columns:

NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:544
    x_WrapOutputLine(dbinfo->definition, line_length, out);

(line_length is kFormatLineLength, 68). The wrap is NStr::Wrap with fWrap_FlatFile
(align_format_util.cpp:283), which breaks only at white space, '-' and ',' and otherwise cuts
at the width:

NCBI reference (598d8ae6): c++/src/corelib/ncbistr.cpp:5228
    else if (do_flat ? c == '-' : ispunct((unsigned char)c)) {

So where the wrapped lines break depends on the path's length and on those characters, and
replacing a checkout path of another length after the run does not give the frozen bytes.

The Owner's answer (2026-11-05, SWEBOK line SW8 item 6, question Q5) was
"(b) 比較する時に、接頭辞を置き換える": replace the prefix when comparing, so that the comparisons pass
wherever the repository is checked out, without a link to FROZEN_PREFIX. How (SW8 DECISIONS 29, Owner-delegated): the runs
pass RUN_ROOT instead of FROZEN_PREFIX (to_run), a relative path of the same length with none of
the break characters, with the checkout root as the working directory, so every line breaks where
NCBI's did; before hashing, the output and the standard error of a run whose arguments to_run
changed get FROZEN_PREFIX back (to_frozen), the same for every output format. The frozen files
stay unchanged. This is the one replacement of a checkout-path prefix and adds no other rule; the
comparisons' existing rules (the SAM @PG CL tag, db_search_fixtures.strip_sam_cl; check_losat.py
normalize_database_lines; dropin_corpus.py's <RUNDIR> and report_without_time; capture_outputs.py's
{OUT} and {LOSAT} in stderr) are their tools' own.
"""
from __future__ import annotations

import re

FROZEN_PREFIX = "/mnt/c/Users/genom/GitHub/LOSAT"
# `.` and `/.` x 15: 31 characters, the checkout root when run from it.
RUN_ROOT = "." + "/." * 15
# Where NCBI's flat-file wrap breaks (white space, '-', ','), and the quote characters of the
# ASN.1 and JSON writers, which break or escape there.
BREAK_CHARACTERS = " -,'\""

for _path in (FROZEN_PREFIX, RUN_ROOT):
    if any(c.isspace() or c in BREAK_CHARACTERS for c in _path):
        raise RuntimeError(f"{_path!r} has a character where an output line can break")
if len(RUN_ROOT) != len(FROZEN_PREFIX):
    raise RuntimeError(f"RUN_ROOT has {len(RUN_ROOT)} characters, FROZEN_PREFIX {len(FROZEN_PREFIX)}")

# The run root as a path's start: not after `/` (inside another path) or `.` (inside `../`).
NOT_AFTER = "/."
_RUN_ROOT_IN_OUTPUT = re.compile(b"(?<![" + re.escape(NOT_AFTER.encode()) + b"])"
                                 + re.escape(RUN_ROOT.encode() + b"/"))


def to_run(arg: str) -> str:
    """The argument to run: a path under FROZEN_PREFIX starts at RUN_ROOT instead."""
    if arg.startswith(FROZEN_PREFIX + "/"):
        return RUN_ROOT + arg[len(FROZEN_PREFIX):]
    return arg


def to_frozen(data: bytes) -> bytes:
    """The bytes of a run of to_run's arguments as the frozen run wrote them: every RUN_ROOT + `/`
    that starts a path is FROZEN_PREFIX + `/`."""
    return _RUN_ROOT_IN_OUTPUT.sub(FROZEN_PREFIX.encode() + b"/", data)


def case_table_entry() -> dict[str, str]:
    """to_frozen's replacement as data, for a runner in another language (web/adapter/tests/v_abi.js)."""
    return {"run": RUN_ROOT + "/", "frozen": FROZEN_PREFIX + "/", "not_after": NOT_AFTER}
