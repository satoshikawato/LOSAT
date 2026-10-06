#!/usr/bin/env python3
"""Keep the verify_refs.py problems that are on lines added since the session start (a92fa902f)."""
import re, subprocess, sys
log = sys.argv[1]
diff = subprocess.run(['git', 'diff', '-U0', 'a92fa902f..HEAD', '--', '*.rs'], capture_output=True, text=True).stdout
added, cur = {}, None
for line in diff.split('\n'):
    if line.startswith('+++ b/'):
        cur = line[6:]; added.setdefault(cur, set())
    m = re.match(r'@@ -\S+ \+(\d+)(?:,(\d+))? @@', line)
    if m and cur:
        st = int(m.group(1)); n = int(m.group(2) or 1); added[cur].update(range(st, st + n))
bad = [l.rstrip() for l in open(log) if (m := re.match(r'(\S+?):(\d+) ', l)) and int(m.group(2)) in added.get(m.group(1), ())]
print('\n'.join(bad)); print(f'# problems on lines added in the session: {len(bad)}')
sys.exit(1 if bad else 0)
