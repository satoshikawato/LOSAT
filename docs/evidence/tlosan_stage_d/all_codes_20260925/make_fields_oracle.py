#!/usr/bin/env python3
"""Derive a diagnostic NCBI API formatter with unchanged search call state."""
from pathlib import Path
import sys
source = Path(sys.argv[1]).read_text()
old = '1, code, options->GetOptions().GetSumStatisticsMode());'
fields = 'qseqid sseqid score bitscore evalue nident positive length mismatch gaps gapopen qstart qend sstart send sframe'
new = f'1, code, options->GetOptions().GetSumStatisticsMode(), false, -1, "{fields}");'
assert source.count(old) == 1
Path(sys.argv[2]).write_text(source.replace(old, new))
