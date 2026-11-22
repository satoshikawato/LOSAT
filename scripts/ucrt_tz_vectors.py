#!/usr/bin/env python3
# `parse` below follows the code of Microsoft's C runtime (UCRT) that Microsoft publishes under the
# MIT License in the NuGet package Microsoft.Windows.SDK.CRTSource 10.0.22621.3
# (ucrt/time/tzset.cpp; also in THIRD_PARTY_NOTICES, section 6 (UCRT)). The package's LICENSE:
#
# MIT License
#
# Copyright (c) Microsoft Corporation.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

"""Vectors of the Windows C runtime's local time with `TZ` set (LOSAT/tests/unit/ucrt_local_time_vectors.tsv).

NCBI BLAST+ on Windows prints makeblastdb's times with UCRT's `localtime`. LOSAT ports UCRT's
`_tzset` and `_localtime64_s` (LOSAT/src/writedb/windows_zone.rs); its unit
tests compare the port with these vectors on every platform, and the Windows test
LOSAT/tests/windows_local_time.rs compares the runner's UCRT with them again.

The values come from UCRT itself: Windows PowerShell calls `_wputenv_s`, `_tzset`,
`_get_timezone`, `_get_daylight`, `_get_dstbias` and `_localtime64_s` of `ucrtbase.dll`, and
`GetStringTypeW` of `kernel32.dll`, in one process. Run it on Windows or under WSL (where it starts
`powershell.exe`):

    python3 scripts/ucrt_tz_vectors.py --out LOSAT/tests/unit/ucrt_local_time_vectors.tsv

Rows (tab-separated; a `TZ` value is written with `\\` for a backslash and `\\u{XXXX}` for a
UTF-16 unit outside printable ASCII):
    space <hex unit> <0|1>          `GetStringTypeW` with `CT_CTYPE1` has `C1_SPACE` (UCRT's `iswspace`
                                    above U+00FF)
    zone  <TZ> <_timezone> <_daylight> <_dstbias>   after `_tzset`
    time  <TZ> <unix seconds> <offset|none>         `_localtime64_s`: local minus UTC seconds,
                                                    `none` when it fails

`TZ` values whose numbers overflow UCRT's 32-bit arithmetic (undefined in C++; LOSAT computes
the mathematical value) are left out, as is the system zone (`TZ` unset), which is the host's.
The instants avoid the first and last three days of `localtime`'s range (another branch, which
LOSAT does not port: the current time is past them).
"""
from __future__ import annotations

import argparse
import datetime
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# The TZ forms of UCRT (`tzn[+|-]hh[:mm[:ss]][dzn]`) and its quirks: no sign or `+`, both signs,
# white space and digits that `wcstol` takes but the scan after it does not, numbers past 24,
# empty parts, POSIX and Olson forms, short and non-ASCII names.
TZ_VALUES = [
    "JST-9", "PST8PDT", "EST5EDT", "CST6CDT", "MST7MDT", "AKST9AKDT", "HST10", "MST7", "GMT0",
    "GMT0BST", "UTC", "UTC0", "CET-1CEST", "IST-5:30", "NPT-5:45", "ACST-9:30ACDT",
    "NZST-12NZDT", "NST3:30NDT", "AAA-14", "AAA+12", "AAA12", "AAA-12:34:56", "AAA12:34:56BBB",
    "AAA1:2:3", "AAA24", "AAA25", "AAA99", "AAA-99BBB", "AAA1000", "AAA-1000BBB", "AAA+8",
    "AAA+8BBB", "AAA++8", "AAA-+8", "AAA--8", "AAA+-8", "AAA -8", "AAA 8", "AAA\t8", "AAA 8",
    "AAA\u00858", "AAA　8", "AAA 8", "AAA​8", "AAA８", "AAA٨",
    "AAA༨BBB", "AAA8:３０", "AAA8:", "AAA8::", "AAA8:30:", "AAA:30", "AAA-:30",
    "AAA8:30:15:45", "AAA8.5", "AAA8,M3.2.0,M11.1.0", "AAA8D", "AB", "A", "ABC", "ABCD", "123-9",
    "-9", "+9", ":Asia/Tokyo", "Asia/Tokyo", "Tokyo Standard Time", "JST-9JDT", "EST5EDT4",
    "<+09>-9", "日本-9", "AAA-596523", "AAA-596523BBB", "AAA596523", "AAA-0:0:1BBB",
    "AAA0:59:59BBB", "AAA-24:00:00BBB", "AAA-48BBB",
]

DAY = 86_400
MAX_TIME64 = 0x7_9358_2AFF
# Only instants of `localtime`'s first branch (ucrt/time/localtime.cpp:58: more than three days
# after the epoch and before max_time_t, 3001-01-19 07:59:59 UTC): LOSAT does not port the branch
# for the three days at either end (the current time is never there), so the vectors claim no
# UCRT version's behaviour there.
INSTANTS = [
    3 * DAY + 1, 946_684_799, 946_684_800, 1_000_000_000, 1_700_000_000, 1_795_000_000,
    1_798_761_599, 1_798_761_600, 2_147_483_647, 2_147_483_648, 4_102_444_799, 4_102_444_800,
    4_107_542_400, MAX_TIME64 - 3 * DAY - 1,
]

PS_SCRIPT = r"""
param([string]$InFile)
$ErrorActionPreference = 'Stop'
Add-Type -TypeDefinition @"
using System; using System.Runtime.InteropServices;
public static class Ucrt {
  [DllImport("ucrtbase.dll", CallingConvention=CallingConvention.Cdecl, CharSet=CharSet.Unicode)]
  public static extern int _wputenv_s(string name, string value);
  [DllImport("ucrtbase.dll", CallingConvention=CallingConvention.Cdecl)]
  public static extern void _tzset();
  [DllImport("ucrtbase.dll", CallingConvention=CallingConvention.Cdecl)]
  public static extern int _localtime64_s([Out] int[] tm, ref long t);
  [DllImport("ucrtbase.dll", CallingConvention=CallingConvention.Cdecl)]
  public static extern int _get_dstbias(out int v);
  [DllImport("ucrtbase.dll", CallingConvention=CallingConvention.Cdecl)]
  public static extern int _get_timezone(out int v);
  [DllImport("ucrtbase.dll", CallingConvention=CallingConvention.Cdecl)]
  public static extern int _get_daylight(out int v);
  [DllImport("kernel32.dll", CharSet=CharSet.Unicode)]
  public static extern int GetStringTypeW(uint info, [MarshalAs(UnmanagedType.LPWStr)] string s, int n, [Out] ushort[] types);
}
"@
$ucrt = Get-Item "$env:SystemRoot\System32\ucrtbase.dll"
[Console]::Out.WriteLine("HOST`t" + [Environment]::OSVersion.VersionString + "`tucrtbase.dll " + $ucrt.VersionInfo.FileVersion)
$b0 = 0; [void][Ucrt]::_get_dstbias([ref]$b0)
[Console]::Out.WriteLine("START`t$b0")
foreach ($line in [System.IO.File]::ReadAllLines($InFile)) {
  $f = $line.Split("`t")
  if ($f[0] -eq 'space') {
    $types = New-Object uint16[] 1
    $ok = [Ucrt]::GetStringTypeW(1, [string][char][Convert]::ToInt32($f[1], 16), 1, $types)
    $sp = 0; if ($ok -ne 0 -and ($types[0] -band 8) -ne 0) { $sp = 1 }
    [Console]::Out.WriteLine("space`t$($f[1])`t$sp")
    continue
  }
  $v = -join ($f[1].Split(' ') | ForEach-Object { [char][Convert]::ToInt32($_, 16) })
  [void][Ucrt]::_wputenv_s('TZ', $v)
  [Ucrt]::_tzset()
  $tz = 0; $dl = 0; $db = 0
  [void][Ucrt]::_get_timezone([ref]$tz); [void][Ucrt]::_get_daylight([ref]$dl); [void][Ucrt]::_get_dstbias([ref]$db)
  $sb = New-Object System.Text.StringBuilder
  [void]$sb.Append("zone`t$($f[1])`t$tz`t$dl`t$db")
  foreach ($t in $f[2].Split(' ')) {
    $tm = New-Object int[] 9
    $tt = [long]$t
    $rc = [Ucrt]::_localtime64_s($tm, [ref]$tt)
    if ($rc -ne 0) { [void]$sb.Append("`t$t=E$rc") } else {
      [void]$sb.Append("`t$t=$($tm[5]),$($tm[7]),$($tm[2]),$($tm[1]),$($tm[0])")
    }
  }
  [Console]::Out.WriteLine($sb.ToString())
}
"""


def units(s: str) -> list[int]:
    b = s.encode("utf-16-le", "surrogatepass")
    return [int.from_bytes(b[i : i + 2], "little") for i in range(0, len(b), 2)]


def escape(s: str) -> str:
    out = []
    for u in units(s):
        if u == 0x5C:
            out.append("\\\\")
        elif 0x20 <= u < 0x7F:
            out.append(chr(u))
        else:
            out.append("\\u{%04X}" % u)
    return "".join(out)


# The parse of `tzset_from_environment_nolock` in exact integers (ucrt/time/tzset.cpp:264-335),
# here only to leave out the values whose 32-bit arithmetic overflows and to place the DST
# changes; the vectors themselves are UCRT's. White space above U+00FF is the set `spaces`.
def parse(tz: str, spaces: set[int]) -> tuple[int, bool, bool]:
    u = units(tz)
    at = lambda i: u[i] if i < len(u) else 0
    p = min(len(u), 3)
    neg = at(p) == ord("-")
    p += neg
    overflow = False

    def wcstol(i: int) -> int:
        while at(i) in (9, 10, 11, 12, 13, 0x20, 0x85, 0xA0) or at(i) in spaces:
            i += 1
        sign = -1 if at(i) == ord("-") else 1
        if at(i) in (ord("-"), ord("+")):
            i += 1
        n, read = 0, False
        while digit(at(i)) is not None:
            n, read, i = n * 10 + digit(at(i)), True, i + 1
        return max(-(2**31), min(2**31 - 1, sign * n)) if read else 0

    def fits(v: int) -> bool:
        return -(2**31) <= v < 2**31

    h = wcstol(p) * 3600
    overflow |= not fits(h)
    tzv = h
    while at(p) == ord("+") or ord("0") <= at(p) <= ord("9"):
        p += 1
    if at(p) == ord(":"):
        p += 1
        m = wcstol(p) * 60
        tzv += m
        overflow |= not fits(m) or not fits(tzv)
        while ord("0") <= at(p) <= ord("9"):
            p += 1
        if at(p) == ord(":"):
            p += 1
            tzv += wcstol(p)
            overflow |= not fits(tzv)
            while ord("0") <= at(p) <= ord("9"):
                p += 1
    if neg:
        tzv = -tzv
        overflow |= not fits(tzv)
    return tzv, at(p) != 0, overflow


def digit(c: int) -> int | None:
    for zero in (0x30, 0x660, 0x6F0, 0x966, 0x9E6, 0xA66, 0xAE6, 0xB66, 0xC66, 0xCE6, 0xD66,
                 0xE50, 0xED0, 0xF20, 0x1040, 0x17E0, 0x1810, 0xFF10):
        if zero <= c < zero + 10:
            return c - zero
    return None


def nth_sunday(year: int, month: int, n: int) -> datetime.date:
    first = datetime.date(year, month, 1)
    day = first + datetime.timedelta(days=(6 - first.weekday()) % 7 + 7 * (n - 1))
    if n == 5:
        while day.month != month:
            day -= datetime.timedelta(days=7)
        nxt = day + datetime.timedelta(days=7)
        if nxt.month == month:
            day = nxt
    return day


def changes(timezone: int) -> list[int]:
    """UTC instants around the United States DST changes of 2006, 2007 and 2026 for a zone of
    `_timezone`: the start at 02:00 standard time, the end at 02:00 daylight time (01:00
    standard time)."""
    out = []
    for year in (2006, 2007, 2026):
        rules = ((4, 1), (10, 5)) if year < 2007 else ((3, 2), (11, 1))
        for (month, week), hour in zip(rules, (2, 1)):
            day = nth_sunday(year, month, week)
            local = (day - datetime.date(1970, 1, 1)).days * DAY + hour * 3600
            out += [local + timezone + d for d in (-1, 0, 1)]
    return out


def offset(field: str, t: int) -> str:
    if field.startswith("E"):
        return "none"
    year, yday, hour, minute, sec = map(int, field.split(","))
    days = (datetime.date(1900 + year, 1, 1) - datetime.date(1970, 1, 1)).days + yday
    return str(days * DAY + hour * 3600 + minute * 60 + sec - t)


def powershell() -> str:
    for name in ("powershell.exe", "powershell"):
        if shutil.which(name):
            return name
    wsl = "/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe"
    if Path(wsl).exists():
        return wsl
    sys.exit("no Windows PowerShell: run this on Windows or under WSL")


def windows_path(path: Path) -> str:
    if platform.system() == "Windows":
        return str(path)
    return subprocess.run(["wslpath", "-w", str(path)], check=True, capture_output=True, text=True).stdout.strip()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    high = sorted({u for tz in TZ_VALUES for u in units(tz) if u > 0xFF})
    with tempfile.TemporaryDirectory(dir=args.out.parent) as tmp:
        tmp = Path(tmp)
        (tmp / "oracle.ps1").write_text(PS_SCRIPT, encoding="utf-8")
        # The spaces first: the instants of a zone depend on its offset.
        (tmp / "in1.txt").write_text("".join(f"space\t{u:04X}\n" for u in high), encoding="utf-8")
        ps = [powershell(), "-NoProfile", "-NonInteractive", "-ExecutionPolicy", "Bypass", "-File",
              windows_path(tmp / "oracle.ps1"), "-InFile"]
        cwd = "/mnt/c/Windows" if platform.system() != "Windows" else None
        first = subprocess.run(ps + [windows_path(tmp / "in1.txt")], check=True, capture_output=True,
                               text=True, cwd=cwd).stdout.replace("\r", "").splitlines()
        spaces = {int(f[1], 16) for f in (l.split("\t") for l in first) if f[0] == "space" and f[2] == "1"}
        rows, lines = [], []
        for tz in TZ_VALUES:
            timezone, daylight, overflow = parse(tz, spaces)
            if overflow:
                continue
            instants = sorted(set(INSTANTS + (changes(timezone) if daylight else [])))
            instants = [t for t in instants if 3 * DAY < t < MAX_TIME64 - 3 * DAY]
            hexed = " ".join(f"{u:04X}" for u in units(tz))
            lines.append(f"tz\t{hexed}\t{' '.join(map(str, instants))}\n")
            rows.append(tz)
        (tmp / "in2.txt").write_text("".join(lines), encoding="utf-8")
        second = subprocess.run(ps + [windows_path(tmp / "in2.txt")], check=True, capture_output=True,
                                text=True, cwd=cwd).stdout.replace("\r", "").splitlines()
    host = next(l for l in first if l.startswith("HOST\t")).split("\t", 1)[1]
    for out in (first, second):
        start = next(l for l in out if l.startswith("START\t"))
        if start != "START\t-3600":
            sys.exit(f"the process had run _tzset before: {start}")
    zones = [l.split("\t") for l in second if l.startswith("zone\t")]
    if len(zones) != len(rows):
        sys.exit(f"{len(zones)} zones for {len(rows)} TZ values")
    body = [f"space\t{u:04X}\t{int(u in spaces)}\n" for u in high]
    for tz, f in zip(rows, zones):
        body.append(f"zone\t{escape(tz)}\t{f[2]}\t{f[3]}\t{f[4]}\n")
        for item in f[5:]:
            t, value = item.split("=")
            body.append(f"time\t{escape(tz)}\t{t}\t{offset(value, int(t))}\n")
    head = (
        f"# UCRT's _tzset and _localtime64_s with TZ set (scripts/ucrt_tz_vectors.py), "
        f"{datetime.date.today().isoformat()}, {host}\n"
    )
    args.out.write_text(head + "".join(body), encoding="utf-8", newline="\n")
    print(f"{args.out}: {len(high)} characters, {len(rows)} TZ values, "
          f"{sum(1 for b in body if b.startswith('time'))} times")
    return 0


if __name__ == "__main__":
    sys.exit(main())
