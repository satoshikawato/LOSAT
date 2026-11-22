#!/usr/bin/env python3
"""glibc's local time for each form of `TZ` and `TZDIR` (LOSAT/tests/unit/local_time_glibc_vectors.tsv).

NCBI BLAST+ on Linux prints makeblastdb's times with glibc's `localtime`. LOSAT finds the zone
where glibc finds it and reads it with the crate jiff (LOSAT/src/writedb/local_offset.rs); its
unit tests compare that with these vectors on every platform but Windows, and CI's Linux job
compares the vectors with its own glibc (`--check`), so that LOSAT is compared with glibc there.

The values come from glibc itself: Python's `time.localtime` (glibc's `localtime_r`) in a process
whose environment holds the form. The local time is written as an offset: the calendar time that
glibc gives, counted in seconds as UTC is, minus the instant (for a zone file with leap seconds,
less the leap seconds; elsewhere glibc's `tm_gmtoff`).

    python3 scripts/local_time_glibc_vectors.py --make-zones   # the zone files, with zic
    python3 scripts/local_time_glibc_vectors.py --out LOSAT/tests/unit/local_time_glibc_vectors.tsv
    python3 scripts/local_time_glibc_vectors.py --check   # the vectors are this host's glibc's

`--make-zones` writes the zone files of the forms into LOSAT/tests/unit/local_time/ with the
host's zic from the host's tzdata.zi: the versions of the format (1: the 32-bit block of a fat
file alone; 2: fat and slim; 3: footers past 24 hours; 4: a version 2 file labelled 4), leap
seconds (tzdata's, and one taken away), a zone file whose footer disagrees with its last transition, and files that are no zone
file. Rows (tab-separated, `{dir}` = that directory; `-` = `TZDIR` unset; `:60` where glibc writes
the second as 60, a leap second):

    form <TZ> <TZDIR> <instant>:<offset>[:60] ...

The instants are fixed ones from 1970 to 2100 and, for a form whose offset changes, each change
in the windows of `WINDOWS` (the instant before it and the instant of it).
"""
from __future__ import annotations

import argparse
import calendar
import datetime
import os
import platform
import shutil
import struct
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ZONES = ROOT / "LOSAT" / "tests" / "unit" / "local_time"
VECTORS = ROOT / "LOSAT" / "tests" / "unit" / "local_time_glibc_vectors.tsv"

FIXED = [0, 31536000, 63072000, 1000000000, 1700000000, 1750000000, 1767225599, 1767225600, 1785000000,
         2147483647, 2147483648, 2208988800, 4102444800, 4118000000, 4133980799]
# Years whose changes of offset are written (the first years, years of changes of law, now, the
# end of the 32-bit instants, the last years).
WINDOWS = [(1970, 1973), (1995, 1997), (2010, 2012), (2016, 2017), (2025, 2028), (2037, 2039), (2099, 2101)]

# Zones built with `zic -b slim` (the footer gives every change after the last transition).
SLIM = ["America/New_York", "America/Ojinaga", "America/Nuuk", "America/Santiago", "America/Scoresbysund",
        "America/Sao_Paulo", "America/St_Johns", "Africa/Casablanca", "Africa/Monrovia", "Antarctica/Troll",
        "Asia/Gaza", "Asia/Jerusalem", "Asia/Kathmandu", "Asia/Kolkata", "Asia/Tehran", "Asia/Tokyo",
        "Australia/Eucla", "Australia/Lord_Howe", "Europe/Dublin", "Europe/London", "Europe/Moscow",
        "Pacific/Apia", "Pacific/Chatham", "Pacific/Easter", "Pacific/Kiritimati"]

POSIX = [
    # Rules that tzdata's footers and NCBI's users write.
    "JST-9", "EST5", "UTC0", "GMT0", "<+09>-9", "<-03>3", "IST-5:30", "<+0545>-5:45", "<+1245>-12:45", "<-0930>9:30",
    "XXX-24", "XXX24", "XXX-24:00", "XXX-23:59:59", "XXX1:02:03",
    "EST5EDT,M3.2.0,M11.1.0", "CET-1CEST,M3.5.0,M10.5.0/3", "NZST-12NZDT,M9.5.0,M4.1.0/3",
    "AEST-10AEDT,M10.1.0,M4.1.0/3", "GMT0BST,M3.5.0/1,M10.5.0", "IST-1GMT0,M10.5.0,M3.5.0/1",
    "<+01>-1<+00>,M3.5.0/3,M10.5.0/2", "XXX5YYY4,M3.2.0,M11.1.0", "XXX5YYY6,M3.2.0,M11.1.0", "XXX5YYY-10,J60,J300",
    "<-02>2<-01>,M3.5.0/-1,M10.5.0/0", "IST-2IDT,M3.4.4/26,M10.5.0", "EET-2EEST,M3.4.4/50,M10.4.4/50",
    "<-04>4<-03>,M9.1.6/24,M4.1.6/24", "<-06>6<-05>,M9.1.6/22,M4.1.6/22", "XXX5YYY,M3.2.0/-167,M11.1.0/167",
    "LHST-10:30LHDT-11,M10.1.0,M4.1.0", "<+1030>-10:30<+11>-11,M10.1.0,M4.1.0",
    "CHAST-12:45CHADT,M9.5.0/2:45,M4.1.0/3:45", "NST3:30NDT,M3.2.0,M11.1.0", "PST8PDT,M3.2.0,M11.1.0",
    "XXX0YYY,M3.4.0,M9.3.0", "XXX5YYY,J60,J300", "XXX5YYY,59,300", "XXX5YYY,0,365", "XXX5YYY,J1,J365",
    "XXX5YYY,0/0,J365/25", "XXX5YYY,J365/25,0/0", "XXX5YYY,M3.2.0/2:30,M11.1.0/2:30:30", "XXX5YYY,M1.1.0,M12.5.6",
    "XXX-14YYY,M3.2.0,M11.1.0", "<A1B2>-1<C3D4>,M3.2.0,M11.1.0", "ABCDEFGHIJKLMNOP-9",
    # Forms that glibc reads otherwise than POSIX writes them (docs/product_decisions/PD-LOSAT-DATABASE-SEARCH.md).
    "XXX+25", "XXX5:70", "XXX99999", "AAA--8", "JST--9", "JST+-9", "XXX-168",
    "EST5EDT,M3.2.0", "EST5EDT,", "EST5EDT,M3,M11.1.0", "EST5EDT,M3.2,M11.1.0", "EST5EDT,,M11.1.0",
    "EST5EDT4", "XXX5YYY", "XXX5YYY4", "XXX-9YYY",
    "XXX1:2", "XXX-5:3", "XXX5YYY4:3,M3.2.0,M11.1.0",
    "JST-9 ", "JST-9:", "JST-9,", "JST-9junk", "EST5EDT,M3.2.0,M11.1.0 ", "EST5EDT,M3.2.0,M11.1.0,M12.1.0",
    "EST5EDT,M13.2.0,M11.1.0", "EST5EDT,M3.6.0,M11.1.0", "EST5EDT,M3.2.7,M11.1.0", "EST5EDT,J0,J300", "EST5EDT,366,300",
    "EST5EDT,m3.2.0,m11.1.0", "EST5EDT,M3.2.0/168,M11.1.0",
    # Rules whose change falls in another year in UTC than in local time (jiff keeps it in its UTC year).
    "XXX-24YYY,0/0,J365/25", "AAA-9BBB,0,365", "AAA-9BBB,0/0,J365/25", "AAA-9BBB,J60,365", "AAA1BBB,J60,365/-3",
    # Forms that glibc reads as UTC too.
    "", ":", "::", " ", "AB", "ABC", "garbage", "Asia/Nonexistent", "日本-9", "Tokyo Standard Time", "<+09>", "JST",
]

# (TZ, TZDIR) beyond the POSIX rules above, which are read with TZDIR={dir} (no file of their name).
FILES = [
    ("Asia/Tokyo", "{dir}"), (":Asia/Tokyo", "{dir}"), ("./Asia/Tokyo", "{dir}"), ("Asia/./Tokyo", "{dir}"),
    ("Asia//Tokyo", "{dir}"), ("../local_time/Asia/Tokyo", "{dir}"), ("Asia/Tokyo/", "{dir}"), ("Asia", "{dir}"),
    ("Dir", "{dir}"), ("notzif", "{dir}"), ("trunc", "{dir}"), ("Asia/Tokyo", "{dir}/"), ("Asia/Tokyo", "{dir}//"),
    ("Tokyo", "{dir}/Asia"), ("Asia/Tokyo", "relative/dir"), ("Asia/Tokyo", "/nonexistent"),
    ("{dir}/Asia/Tokyo", "-"), (":{dir}/Asia/Tokyo", "-"), ("/{dir}/Asia/Tokyo", "-"), ("{dir}/notzif", "-"),
    ("{dir}/Asia", "-"), ("/nonexistent/zone", "-"), ("/", "-"),
    ("JST-9", "{dir}/order"), (":JST-9", "{dir}/order"), ("Asia/Tokyo", "{dir}/order"),
    ("v1/America/New_York", "{dir}"), ("v1/Australia/Lord_Howe", "{dir}"), ("fat/America/New_York", "{dir}"),
    ("v4/America/New_York", "{dir}"), ("v4/Pacific/Easter", "{dir}"),
    ("leap/Asia/Tokyo", "{dir}"), ("leap/America/New_York", "{dir}"), ("leap/Europe/London", "{dir}"),
    ("v1/leap/Asia/Tokyo", "{dir}"), ("leap/negative", "{dir}"), ("ver1/leap/Asia/Tokyo", "{dir}"),
    ("tail/Asia/Tokyo", "{dir}"),
] + [(f"slim/{zone}", "{dir}") for zone in SLIM]


def make_zones(zic: str, tzdata: Path, leapseconds: Path) -> None:
    """Writes the zone files of the forms into ZONES (replacing it)."""
    shutil.rmtree(ZONES, ignore_errors=True)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for kind, options in [("slim", ["-b", "slim"]), ("fat", ["-b", "fat"]),
                              ("leap", ["-b", "fat", "-L", str(leapseconds)]),
                              ("leapslim", ["-b", "slim", "-L", str(leapseconds)])]:
            subprocess.run([zic, "-d", str(tmp / kind), *options, str(tzdata)], check=True, stderr=subprocess.DEVNULL)

        def put(name: str, data: bytes) -> None:
            path = ZONES / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)

        def read(kind: str, zone: str) -> bytes:
            return (tmp / kind / zone).read_bytes()

        def v1_block(data: bytes) -> int:
            isut, isstd, leap, times, types, chars = struct.unpack(">6l", data[20:44])
            return 44 + times * 5 + types * 6 + chars + leap * 8 + isstd + isut

        for zone in SLIM:
            put(f"slim/{zone}", read("slim", zone))
        put("Asia/Tokyo", read("slim", "Asia/Tokyo"))
        put("Dir/Zone", read("slim", "Asia/Tokyo"))
        put("order/JST-9", read("slim", "America/New_York"))
        put("order/Asia/Tokyo", read("slim", "America/New_York"))
        put("notzif", b"not a zone file\n")
        put("trunc", read("slim", "America/New_York")[:30])
        put("fat/America/New_York", read("fat", "America/New_York"))
        for zone in ["America/New_York", "Australia/Lord_Howe"]:
            data = bytearray(read("fat", zone)[:v1_block(read("fat", zone))])
            data[4] = 0
            put(f"v1/{zone}", bytes(data))
        for zone, kind in [("America/New_York", "fat"), ("Pacific/Easter", "slim")]:
            data = bytearray(read(kind, zone))
            data[4] = data[v1_block(data) + 4] = ord("4")
            put(f"v4/{zone}", bytes(data))
        for zone in ["Asia/Tokyo", "America/New_York", "Europe/London"]:
            put(f"leap/{zone}", read("leap", zone))
        data = bytearray(read("leap", "Asia/Tokyo")[:v1_block(read("leap", "Asia/Tokyo"))])
        data[4] = 0
        put("v1/leap/Asia/Tokyo", bytes(data))
        # Version '1' (no version of RFC 8536): glibc and jiff read the second block, which alone has the
        # leap seconds of a slim file.
        data = bytearray(read("leapslim", "Asia/Tokyo"))
        data[4] = data[v1_block(data) + 4] = ord("1")
        put("ver1/leap/Asia/Tokyo", bytes(data))
        # A zone file followed by more than a MiB of zeros, which glibc does not read.
        put("tail/Asia/Tokyo", read("slim", "Asia/Tokyo") + bytes((1 << 20) + 1))
        # Leap seconds that tzdata has never had: one taken away (1996) between two added.
        (tmp / "negative.zi").write_text("Zone negative 9:00 - JST\n")
        (tmp / "negative.leap").write_text("Leap 1972 Jun 30 23:59:60 + S\nLeap 1996 Jun 30 23:59:58 - S\n"
                                           "Leap 2016 Dec 31 23:59:60 + S\n")
        subprocess.run([zic, "-d", str(tmp / "negative"), "-b", "fat", "-L", str(tmp / "negative.leap"),
                        str(tmp / "negative.zi")], check=True, stderr=subprocess.DEVNULL)
        put("leap/negative", read("negative", "negative"))


# glibc's offsets: on stdin, one job a line, "TZ<tab>TZDIR<tab>spans<tab>instants" (TZ or TZDIR "\0" =
# unset; spans "lo:hi,..." whose changes of offset are written; instants "t,..."); on stdout, one line a
# job, "instant:offset ..." (":60" added where glibc writes the second as 60).
SERVER = r'''
import calendar, os, sys, time
def offset(t):
    return calendar.timegm(time.localtime(t)) - t
def point(t):
    tm = time.localtime(t)
    return f"{t}:{calendar.timegm(tm) - t}" + (":60" if tm.tm_sec == 60 else "")
def setenv(name, value):
    if value == "\0": os.environ.pop(name, None)
    else: os.environ[name] = value
def changes(lo, hi, points):
    t, before = lo, offset(lo)
    while t < hi:
        u = min(t + 7 * 86400, hi); after = offset(u)
        if after != before:
            a, b = t, u
            while b - a > 1:
                m = (a + b) // 2
                if offset(m) == before: a = m
                else: b = m
            points.update((b - 1, b))
        t, before = u, after
for line in sys.stdin:
    tz, tzdir, spans, fixed = line.rstrip("\n").split("\t")
    os.environ["TZ"] = "UTC0"; time.tzset()
    setenv("TZDIR", tzdir); setenv("TZ", tz); time.tzset()
    points = set(int(t) for t in fixed.split(",") if t)
    for span in spans.split(","):
        if span:
            lo, hi = span.split(":")
            changes(int(lo), int(hi), points)
    print(" ".join(point(t) for t in sorted(points)), flush=True)
'''


def glibc_offsets(jobs: list[tuple[str | None, str | None, list[tuple[int, int]], list[int]]]) -> list[str]:
    """glibc's offsets of each job (TZ or None for unset, TZDIR or None, spans, fixed instants)."""
    lines = []
    for tz, tzdir, spans, fixed in jobs:
        assert tz is None or not any(c in tz for c in "\t\n\0"), tz
        lines.append("\t".join(["\0" if tz is None else tz, "\0" if tzdir is None else tzdir,
                                ",".join(f"{lo}:{hi}" for lo, hi in spans), ",".join(map(str, fixed))]))
    env = {k: v for k, v in os.environ.items() if k not in ("TZ", "TZDIR")}
    out = subprocess.run([sys.executable, "-I", "-c", SERVER], input="\n".join(lines) + "\n", env=env,
                         capture_output=True, text=True, check=True).stdout
    return out.splitlines()


def check(forms: list[tuple[str, str]], rows: list[str]) -> int:
    """1 (with the forms printed) where this host's glibc gives other times than the vectors."""
    vectors = {}
    for line in VECTORS.read_text(encoding="utf-8").splitlines():
        if not line.startswith("#"):
            _, tz, tzdir, row = line.split("\t")
            if (tz, tzdir) in vectors:
                print(f"a second row of TZ={tz!r} TZDIR={tzdir!r}")
                return 1
            vectors[(tz, tzdir)] = row
    unequal = [form for form, row in zip(forms, rows) if vectors.pop(form, None) != row]
    for tz, tzdir in unequal + list(vectors):
        print(f"differs from glibc {platform.libc_ver()[1]}: TZ={tz!r} TZDIR={tzdir!r}")
    print(f"{len(forms)} forms, {len(unequal) + len(vectors)} differ from this host's glibc")
    return 1 if unequal or vectors else 0


def year_start(year: int) -> int:
    return calendar.timegm((year, 1, 1, 0, 0, 0))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--make-zones", action="store_true", help="write the zone files with zic")
    parser.add_argument("--zic", default="zic")
    parser.add_argument("--tzdata", type=Path, default=Path("/usr/share/zoneinfo/tzdata.zi"))
    parser.add_argument("--leapseconds", type=Path, default=Path("/usr/share/zoneinfo/leapseconds"))
    parser.add_argument("--out", type=Path)
    parser.add_argument("--check", action="store_true", help="compare the vectors with this host's glibc")
    args = parser.parse_args()
    if platform.libc_ver()[0] != "glibc":
        sys.exit("needs glibc (Linux)")
    if args.make_zones:
        make_zones(args.zic, args.tzdata, args.leapseconds)
    if not (args.out or args.check):
        return
    spans = [(year_start(a), year_start(b)) for a, b in WINDOWS]
    forms = [(tz, "{dir}") for tz in POSIX] + FILES
    jobs = [(tz.replace("{dir}", str(ZONES)), None if d == "-" else d.replace("{dir}", str(ZONES)), spans, FIXED)
            for tz, d in forms]
    rows = glibc_offsets(jobs)
    if args.check:
        sys.exit(check(forms, rows))
    tzdata = args.tzdata.read_text().splitlines()[0].removeprefix("# version ").strip()
    with open(args.out, "w", encoding="utf-8", newline="\n") as out:
        out.write(f"# glibc's local time (scripts/local_time_glibc_vectors.py), {datetime.date.today()}, "
                  f"glibc {platform.libc_ver()[1]}, zone files from tzdata {tzdata}\n")
        for (tz, tzdir), row in zip(forms, rows):
            out.write(f"form\t{tz}\t{tzdir}\t{row}\n")


if __name__ == "__main__":
    main()
