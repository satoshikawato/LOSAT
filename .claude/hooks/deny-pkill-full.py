#!/usr/bin/env python3
"""PreToolUse hook (Bash): deny `pkill -f` and `pgrep -f`.

A full-command-line pattern also matches the shell that runs it whenever the
pattern appears in that shell's own command line (a wait loop, a `bash -c`
wrapper, a gate script that names its stages). LOSAT sessions killed their own
shells and gate runs this way, and wait loops never ended. Track processes by
PID instead.
"""
from __future__ import annotations

import json
import re
import sys

FULL = re.compile(r'\bp(?:kill|grep)\b[^;|&\n]*?\s(?:-[A-Za-z]*f[A-Za-z]*|--full)(?=\s|$)')


def main():
    try:
        payload = json.load(sys.stdin)
    except (json.JSONDecodeError, ValueError):
        return
    command = (payload.get('tool_input') or {}).get('command') or ''
    if not FULL.search(command):
        return
    print(json.dumps({
        'hookSpecificOutput': {
            'hookEventName': 'PreToolUse',
            'permissionDecision': 'deny',
            'permissionDecisionReason': (
                '`pkill -f` / `pgrep -f` match their own shell when the pattern is in its command line. '
                'Keep the PID (`cmd & echo $! > "$TASK_DIR/<name>.pid"`) and use `kill -0 "$PID"`, '
                '`kill "$PID"`, or `pkill -P "$PID"`; to find a program by name use `pgrep -x <name>`.'
            ),
        }
    }))


if __name__ == '__main__':
    main()
