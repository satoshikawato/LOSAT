#!/usr/bin/env python3
"""
Analyze and compare BLAST/LOSAT output files.

This script compares LOSAT and NCBI BLAST+ outputs to identify:
- Hit count differences
- Score distribution differences
- Alignment length distribution differences
- Potential self-comparison issues

Usage:
    python analyze_blast_results.py <losat_output> <blast_output>
    python analyze_blast_results.py --dir <results_directory>
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict
import statistics

def parse_blast6(filepath):
    """Parse BLAST tabular format (outfmt 6) output."""
    hits = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 12:
                hit = {
                    'query': fields[0],
                    'subject': fields[1],
                    'identity': float(fields[2]),
                    'length': int(fields[3]),
                    'mismatches': int(fields[4]),
                    'gaps': int(fields[5]),
                    'qstart': int(fields[6]),
                    'qend': int(fields[7]),
                    'sstart': int(fields[8]),
                    'send': int(fields[9]),
                    'evalue': float(fields[10]),
                    'bitscore': float(fields[11]),
                }
                hits.append(hit)
    return hits

def analyze_hits(hits, name):
    """Analyze hit statistics."""
    if not hits:
        print(f"\n{name}: No hits found")
        return
    
    lengths = [h['length'] for h in hits]
    scores = [h['bitscore'] for h in hits]
    identities = [h['identity'] for h in hits]
    evalues = [h['evalue'] for h in hits]
    
    print(f"\n{name}:")
    print(f"  Total hits: {len(hits)}")
    
    print(f"  Alignment length:")
    print(f"    Min: {min(lengths)}, Max: {max(lengths)}")
    print(f"    Mean: {statistics.mean(lengths):.1f}, Median: {statistics.median(lengths)}")
    if len(lengths) > 1:
        print(f"    Std: {statistics.stdev(lengths):.1f}")
    
    print(f"  Bit score:")
    print(f"    Min: {min(scores):.1f}, Max: {max(scores):.1f}")
    print(f"    Mean: {statistics.mean(scores):.1f}, Median: {statistics.median(scores):.1f}")
    
    print(f"  Identity:")
    print(f"    Min: {min(identities):.1f}%, Max: {max(identities):.1f}%")
    print(f"    Mean: {statistics.mean(identities):.1f}%")
    
    # Check for very long alignments (potential self-comparison issue)
    very_long = [h for h in hits if h['length'] > 10000]
    if very_long:
        print(f"  WARNING: {len(very_long)} alignments > 10kb detected!")
        for h in very_long[:5]:  # Show first 5
            print(f"    {h['query']} vs {h['subject']}: {h['length']} bp, {h['bitscore']:.1f} bits")
    
    # Check for self-hits
    self_hits = [h for h in hits if h['query'] == h['subject']]
    if self_hits:
        print(f"  Self-hits: {len(self_hits)}")

def length_histogram(hits, name, bins=[0, 100, 500, 1000, 5000, 10000, 50000, float('inf')]):
    """Create a length histogram."""
    if not hits:
        return
    
    counts = defaultdict(int)
    for h in hits:
        for i in range(len(bins) - 1):
            if bins[i] <= h['length'] < bins[i+1]:
                label = f"{bins[i]}-{bins[i+1]-1 if bins[i+1] != float('inf') else '+'}"
                counts[label] += 1
                break
    
    print(f"\n{name} - Length distribution:")
    for i in range(len(bins) - 1):
        label = f"{bins[i]}-{bins[i+1]-1 if bins[i+1] != float('inf') else '+'}"
        count = counts[label]
        bar = '#' * min(count, 50)
        print(f"  {label:>12}: {count:>6} {bar}")

def compare_results(losat_hits, blast_hits):
    """Compare LOSAT and BLAST results."""
    print("\n" + "=" * 60)
    print("Comparison Summary")
    print("=" * 60)
    
    losat_count = len(losat_hits)
    blast_count = len(blast_hits)
    
    print(f"\nHit counts:")
    print(f"  LOSAT: {losat_count}")
    print(f"  BLAST: {blast_count}")
    
    if blast_count > 0:
        ratio = losat_count / blast_count * 100
        print(f"  Ratio: {ratio:.1f}%")
        
        if ratio > 120:
            print("  WARNING: LOSAT produces significantly more hits")
        elif ratio < 80:
            print("  WARNING: LOSAT produces significantly fewer hits")
    
    # Compare score ranges
    if losat_hits and blast_hits:
        losat_max_score = max(h['bitscore'] for h in losat_hits)
        blast_max_score = max(h['bitscore'] for h in blast_hits)
        losat_max_len = max(h['length'] for h in losat_hits)
        blast_max_len = max(h['length'] for h in blast_hits)
        
        print(f"\nMax bit score:")
        print(f"  LOSAT: {losat_max_score:.1f}")
        print(f"  BLAST: {blast_max_score:.1f}")
        
        print(f"\nMax alignment length:")
        print(f"  LOSAT: {losat_max_len}")
        print(f"  BLAST: {blast_max_len}")
        
        if losat_max_len > blast_max_len * 10:
            print("  WARNING: LOSAT max length is >10x BLAST max length")
            print("           This may indicate a self-comparison issue")

def analyze_directory(dirpath):
    """Analyze all result files in a directory."""
    dirpath = Path(dirpath)
    
    tsv_files = list(dirpath.glob("*.tsv"))
    if not tsv_files:
        print(f"No .tsv files found in {dirpath}")
        return
    
    results = {}
    for f in tsv_files:
        name = f.stem
        hits = parse_blast6(f)
        results[name] = hits
        analyze_hits(hits, name)
        length_histogram(hits, name)
    
    # Compare if we have both losat and blast results
    losat_results = {k: v for k, v in results.items() if 'losat' in k.lower()}
    blast_results = {k: v for k, v in results.items() if 'blast' in k.lower() and 'losat' not in k.lower()}
    
    if losat_results and blast_results:
        losat_key = list(losat_results.keys())[0]
        blast_key = list(blast_results.keys())[0]
        compare_results(losat_results[losat_key], blast_results[blast_key])

def main():
    parser = argparse.ArgumentParser(description='Analyze BLAST/LOSAT results')
    parser.add_argument('--dir', help='Directory containing result files')
    parser.add_argument('files', nargs='*', help='Individual result files to analyze')
    
    args = parser.parse_args()
    
    if args.dir:
        analyze_directory(args.dir)
    elif args.files:
        for f in args.files:
            hits = parse_blast6(f)
            analyze_hits(hits, Path(f).stem)
            length_histogram(hits, Path(f).stem)
        
        if len(args.files) == 2:
            losat_hits = parse_blast6(args.files[0])
            blast_hits = parse_blast6(args.files[1])
            compare_results(losat_hits, blast_hits)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()

