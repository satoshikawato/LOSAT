#!/usr/bin/env python3
"""Line BDC, BDC5 probe O: the line ends of the batch files of blastdbcmd
(-entry_batch, -taxidlist, -ipg_batch) on one platform.

NcbiGetlineEOL (corelib/ncbistre.cpp:196-208) splits at '\\n' on Linux, at '\\n' and
then drops one '\\r' on Windows, at '\\r' or '\\n' on macOS; CArg_InputFile opens the
file (and stdin, '-') in text mode on Windows (corelib/ncbiargs.cpp:695-737). The
inputs are written here, byte for byte, so that no checkout converts them.

usage: bdc5_eol_probe.py --run-dir LOSAT/tests/blastdbcmd --out DIR -- PROGRAM...
PROGRAM is NCBI's blastdbcmd, or `losat blastdbcmd`. Every run: cwd = --run-dir, no
BLASTDB / NCBI_* / LOSAT_* variable, timeout 60 s. Writes DIR/in/*.txt, DIR/<id>.stdout,
DIR/<id>.stderr and DIR/manifest.tsv (id, rc, stdout bytes and SHA-256, stderr bytes and
SHA-256, stdin, argv).
"""

import argparse
import hashlib
import os
import platform
import subprocess
import sys
from pathlib import Path

# P00019.3 and AAA00018.1 hold no "_" or "|": a version 5 database reads them with NStr::ToUpper
# alone (s_PreProcessAccessionsForDBv5), which keeps a '\r'; CSeq_id's parser drops it from
# XP_000001.1 (e_us).
E1, E2 = "P00019.3", "AAA00018.1"

# name -> bytes. e_: -entry_batch of db/p5; t_: -taxidlist of j/lst/p5; p_: -ipg_batch of
# db/m/pmg (PIGs 1001 and 1003; "21" in the second field is a mask algorithm of pmg).
INPUTS = {
    "e_lf": f"{E1}\n{E2}\n",
    "e_crlf": f"{E1}\r\n{E2}\r\n",
    "e_cr": f"{E1}\r{E2}\r",
    "e_crcrlf": f"{E1}\r\r\n{E2}\r\n",
    "e_lfcr": f"{E1}\n\r{E2}\n",
    "e_mid": f"{E1}\r{E2}\n",
    "e_crlast": f"{E1}\r\n{E2}\r",
    "e_blank": f"\r\n\r\n{E1}\r\n \r\n\t\r\n{E2}\n",
    "e_sub": f"{E1}\n\x1a{E2}\n",
    "e_subin": f"{E1}\x1a\n{E2}\n",
    "e_crsub": f"{E1}\r\x1a\n{E2}\n",
    "e_subend": f"{E1}\n{E2}\n\x1a",
    "e_rng": f"{E2} 2-10\r\n{E1} 3-12\n",
    "e_rngcr": f"{E2} 2-10\r{E1} 3-12\r\n",
    "e_tab": f"{E2}\t2-10\r\n{E1}\t \r\n",
    "e_cr2": f"{E1}\r\r{E2}\r\n",
    "e_crlfcr": f"{E1}\r\n\r{E2}\r\n",
    "e_us": "XP_000001.1\r\nXP_000020.1\r\n",
    "t_lf": "9606\n562\n",
    "t_crlf": "9606\r\n562\r\n",
    "t_cr": "9606\r562\r",
    "t_crcrlf": "9606\r\r\n562\n",
    "t_blank": "\r\n9606\r\n",
    "t_sub": "9606\n\x1a562\n",
    "t_subend": "9606\n562\n\x1a",
    "p_lf": "1001 21\n1003\n",
    "p_crlf": "1001 21\r\n1003\r\n",
    "p_cr": "1001 21\r1003\r",
    "p_crcrlf": "1001 21\r\r\n1003\n",
    "p_sub": "1001 21\n\x1a1003\n",
}

ENTRY = ["-db", "db/p5", "-entry_batch"]
TAX = ["-db", "j/lst/p5", "-no_taxid_expansion", "-outfmt", "%a %T", "-taxidlist"]
PIG = ["-db", "db/m/pmg", "-ipg_batch"]

# (id, argv before the file, input, the file on the standard input instead, argv after)
CASES = [
    ("O01", ENTRY, "e_lf", False, ["-outfmt", "%a"]),
    ("O02", ENTRY, "e_crlf", False, ["-outfmt", "%a"]),
    ("O03", ENTRY, "e_cr", False, ["-outfmt", "%a"]),
    ("O04", ENTRY, "e_crcrlf", False, ["-outfmt", "%a"]),
    ("O05", ENTRY, "e_lfcr", False, ["-outfmt", "%a"]),
    ("O06", ENTRY, "e_mid", False, ["-outfmt", "%a"]),
    ("O07", ENTRY, "e_crlast", False, ["-outfmt", "%a"]),
    ("O08", ENTRY, "e_blank", False, ["-outfmt", "%a"]),
    ("O09", ENTRY, "e_sub", False, ["-outfmt", "%a"]),
    ("O10", ENTRY, "e_subin", False, ["-outfmt", "%a"]),
    ("O11", ENTRY, "e_crsub", False, ["-outfmt", "%a"]),
    ("O12", ENTRY, "e_subend", False, ["-outfmt", "%a"]),
    ("O13", ENTRY, "e_rng", False, []),
    ("O14", ENTRY, "e_rngcr", False, []),
    ("O15", ENTRY, "e_tab", False, []),
    ("O16", ENTRY, "e_cr2", False, ["-outfmt", "%a"]),
    ("O17", ENTRY, "e_crlfcr", False, ["-outfmt", "%a"]),
    ("O18", ENTRY, "e_crlf", False, ["-outfmt", "%a", "-get_dups"]),
    ("O19", ENTRY, "e_cr", False, ["-outfmt", "%a", "-get_dups"]),
    ("O20", ENTRY, "e_sub", False, ["-outfmt", "%a", "-get_dups"]),
    ("O21", ENTRY, "e_crlf", True, ["-outfmt", "%a"]),
    ("O22", ENTRY, "e_cr", True, ["-outfmt", "%a"]),
    ("O23", ENTRY, "e_sub", True, ["-outfmt", "%a"]),
    ("O24", ENTRY, "e_crcrlf", True, ["-outfmt", "%a"]),
    ("O25", ENTRY, "e_crlf", False, ["-outfmt", "%a", "-target_only"]),
    ("O26", ENTRY, "e_us", False, ["-outfmt", "%a"]),
    ("O31", TAX, "t_lf", False, []),
    ("O32", TAX, "t_crlf", False, []),
    ("O33", TAX, "t_cr", False, []),
    ("O34", TAX, "t_crcrlf", False, []),
    ("O35", TAX, "t_blank", False, []),
    ("O36", TAX, "t_sub", False, []),
    ("O37", TAX, "t_subend", False, []),
    ("O38", TAX, "t_crlf", True, []),
    ("O39", TAX, "t_cr", True, []),
    ("O41", PIG, "p_lf", False, []),
    ("O42", PIG, "p_crlf", False, []),
    ("O43", PIG, "p_cr", False, []),
    ("O44", PIG, "p_crcrlf", False, []),
    ("O45", PIG, "p_sub", False, []),
    ("O46", PIG, "p_crlf", True, []),
    ("O47", PIG, "p_cr", True, ["-outfmt", "%o %a"]),
]


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("program", nargs="+")
    args = parser.parse_args()
    out = args.out.resolve()
    (out / "in").mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, text in INPUTS.items():
        path = out / "in" / f"{name}.txt"
        path.write_bytes(text.encode("ascii"))
        paths[name] = path
    env = {
        k: v
        for k, v in os.environ.items()
        if not (k.startswith("NCBI_") or k.startswith("LOSAT_") or k.startswith("BLASTDB"))
    }
    rows = ["id\trc\tout_bytes\tout_sha256\terr_bytes\terr_sha256\tstdin\targv"]
    for cid, head, name, on_stdin, tail in CASES:
        argv = head + (["-"] if on_stdin else [str(paths[name])]) + tail
        stdin = open(paths[name], "rb") if on_stdin else subprocess.DEVNULL
        try:
            done = subprocess.run(
                args.program + argv,
                cwd=args.run_dir,
                env=env,
                stdin=stdin,
                capture_output=True,
                timeout=60,
            )
            rc, so, se = done.returncode, done.stdout, done.stderr
        except subprocess.TimeoutExpired as expired:
            rc, so, se = "timeout", expired.stdout or b"", expired.stderr or b""
        finally:
            if on_stdin:
                stdin.close()
        (out / f"{cid}.stdout").write_bytes(so)
        (out / f"{cid}.stderr").write_bytes(se)
        shown = " ".join(head + ([f"- < in/{name}.txt"] if on_stdin else [f"in/{name}.txt"]) + tail)
        rows.append(
            f"{cid}\t{rc}\t{len(so)}\t{sha(so)}\t{len(se)}\t{sha(se)}\t{name if on_stdin else '-'}\t{shown}"
        )
    rows.append("RUNDONE")
    (out / "manifest.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8", newline="\n")
    version = subprocess.run(args.program + ["-version"], capture_output=True, env=env)
    (out / "platform.txt").write_text(
        f"platform={platform.platform()}\nmachine={platform.machine()}\npython={sys.version.split()[0]}\n"
        f"program={' '.join(args.program)}\nversion_rc={version.returncode}\n",
        encoding="utf-8",
        newline="\n",
    )
    (out / "version.stdout").write_bytes(version.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
