#!/usr/bin/env python3
"""
Compare LOSAT and NCBI BLAST+ integration test results.
Analyzes hit counts, distributions (length, e-value, bitscore, identity).
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

def analyze_distribution(values, name):
    """Calculate distribution statistics."""
    if not values:
        return {}
    
    sorted_vals = sorted(values)
    return {
        'count': len(values),
        'min': min(values),
        'max': max(values),
        'mean': statistics.mean(values),
        'median': statistics.median(values),
        'std': statistics.stdev(values) if len(values) > 1 else 0.0,
        'q25': sorted_vals[len(sorted_vals) // 4] if len(sorted_vals) > 0 else 0,
        'q75': sorted_vals[3 * len(sorted_vals) // 4] if len(sorted_vals) > 0 else 0,
    }

def compare_results(losat_file, blast_file):
    """Compare LOSAT and BLAST results."""
    losat_hits = parse_blast6(losat_file)
    blast_hits = parse_blast6(blast_file)
    
    print(f"\n{'='*80}")
    print(f"Comparison: {Path(losat_file).stem} vs {Path(blast_file).stem}")
    print(f"{'='*80}")
    
    # Hit counts
    print(f"\n## Hit Counts")
    print(f"| Tool | Count |")
    print(f"|------|-------|")
    print(f"| LOSAT | {len(losat_hits):,} |")
    print(f"| NCBI BLAST+ | {len(blast_hits):,} |")
    if len(blast_hits) > 0:
        ratio = len(losat_hits) / len(blast_hits) * 100
        print(f"| Ratio | {ratio:.1f}% |")
    
    # Length distribution
    if losat_hits and blast_hits:
        losat_lengths = [h['length'] for h in losat_hits]
        blast_lengths = [h['length'] for h in blast_hits]
        
        losat_len_stats = analyze_distribution(losat_lengths, "Length")
        blast_len_stats = analyze_distribution(blast_lengths, "Length")
        
        print(f"\n## Alignment Length Distribution")
        print(f"| Statistic | LOSAT | NCBI BLAST+ | Ratio |")
        print(f"|-----------|-------|-------------|-------|")
        print(f"| Count | {losat_len_stats['count']:,} | {blast_len_stats['count']:,} | {losat_len_stats['count']/blast_len_stats['count']*100 if blast_len_stats['count'] > 0 else 0:.1f}% |")
        print(f"| Min | {losat_len_stats['min']} | {blast_len_stats['min']} | - |")
        print(f"| Max | {losat_len_stats['max']} | {blast_len_stats['max']} | - |")
        print(f"| Mean | {losat_len_stats['mean']:.1f} | {blast_len_stats['mean']:.1f} | {losat_len_stats['mean']/blast_len_stats['mean']*100 if blast_len_stats['mean'] > 0 else 0:.1f}% |")
        print(f"| Median | {losat_len_stats['median']:.0f} | {blast_len_stats['median']:.0f} | - |")
        print(f"| Q25 | {losat_len_stats['q25']:.0f} | {blast_len_stats['q25']:.0f} | - |")
        print(f"| Q75 | {losat_len_stats['q75']:.0f} | {blast_len_stats['q75']:.0f} | - |")
    
    # Bit score distribution
    if losat_hits and blast_hits:
        losat_scores = [h['bitscore'] for h in losat_hits]
        blast_scores = [h['bitscore'] for h in blast_hits]
        
        losat_score_stats = analyze_distribution(losat_scores, "Bit Score")
        blast_score_stats = analyze_distribution(blast_scores, "Bit Score")
        
        print(f"\n## Bit Score Distribution")
        print(f"| Statistic | LOSAT | NCBI BLAST+ | Ratio |")
        print(f"|-----------|-------|-------------|-------|")
        print(f"| Min | {losat_score_stats['min']:.1f} | {blast_score_stats['min']:.1f} | - |")
        print(f"| Max | {losat_score_stats['max']:.1f} | {blast_score_stats['max']:.1f} | - |")
        print(f"| Mean | {losat_score_stats['mean']:.1f} | {blast_score_stats['mean']:.1f} | {losat_score_stats['mean']/blast_score_stats['mean']*100 if blast_score_stats['mean'] > 0 else 0:.1f}% |")
        print(f"| Median | {losat_score_stats['median']:.1f} | {blast_score_stats['median']:.1f} | - |")
    
    # E-value distribution
    if losat_hits and blast_hits:
        losat_evalues = [h['evalue'] for h in losat_hits]
        blast_evalues = [h['evalue'] for h in blast_hits]
        
        # Log scale for e-values
        losat_log_eval = [max(-300, min(0, h['evalue'])) if h['evalue'] > 0 else -300 for h in losat_hits]
        blast_log_eval = [max(-300, min(0, h['evalue'])) if h['evalue'] > 0 else -300 for h in blast_hits]
        
        losat_eval_stats = analyze_distribution(losat_evalues, "E-value")
        blast_eval_stats = analyze_distribution(blast_evalues, "E-value")
        
        print(f"\n## E-value Distribution")
        print(f"| Statistic | LOSAT | NCBI BLAST+ |")
        print(f"|-----------|-------|-------------|")
        print(f"| Min | {losat_eval_stats['min']:.2e} | {blast_eval_stats['min']:.2e} |")
        print(f"| Max | {losat_eval_stats['max']:.2e} | {blast_eval_stats['max']:.2e} |")
        print(f"| Mean | {losat_eval_stats['mean']:.2e} | {blast_eval_stats['mean']:.2e} |")
        print(f"| Median | {losat_eval_stats['median']:.2e} | {blast_eval_stats['median']:.2e} |")
    
    # Identity distribution
    if losat_hits and blast_hits:
        losat_identities = [h['identity'] for h in losat_hits]
        blast_identities = [h['identity'] for h in blast_hits]
        
        losat_id_stats = analyze_distribution(losat_identities, "Identity")
        blast_id_stats = analyze_distribution(blast_identities, "Identity")
        
        print(f"\n## Identity Distribution (%)")
        print(f"| Statistic | LOSAT | NCBI BLAST+ |")
        print(f"|-----------|-------|-------------|")
        print(f"| Min | {losat_id_stats['min']:.1f}% | {blast_id_stats['min']:.1f}% |")
        print(f"| Max | {losat_id_stats['max']:.1f}% | {blast_id_stats['max']:.1f}% |")
        print(f"| Mean | {losat_id_stats['mean']:.1f}% | {blast_id_stats['mean']:.1f}% |")
        print(f"| Median | {losat_id_stats['median']:.1f}% | {blast_id_stats['median']:.1f}% |")
    
    # Score ranges
    if losat_hits and blast_hits:
        print(f"\n## Score Range Distribution")
        ranges = [
            (0, 30, "< 30 bits"),
            (30, 50, "30-50 bits"),
            (50, 100, "50-100 bits"),
            (100, float('inf'), ">= 100 bits"),
        ]
        print(f"| Range | LOSAT | NCBI BLAST+ | Ratio |")
        print(f"|-------|-------|-------------|-------|")
        for min_score, max_score, label in ranges:
            losat_count = sum(1 for h in losat_hits if min_score <= h['bitscore'] < max_score)
            blast_count = sum(1 for h in blast_hits if min_score <= h['bitscore'] < max_score)
            ratio = losat_count / blast_count * 100 if blast_count > 0 else 0
            print(f"| {label} | {losat_count:,} | {blast_count:,} | {ratio:.1f}% |")

def main():
    parser = argparse.ArgumentParser(description='Compare LOSAT and NCBI BLAST+ integration test results')
    parser.add_argument('--losat-dir', default='./losat_out', help='Directory containing LOSAT output files')
    parser.add_argument('--blast-dir', default='./blast_out', help='Directory containing NCBI BLAST+ output files')
    
    args = parser.parse_args()
    
    losat_dir = Path(args.losat_dir)
    blast_dir = Path(args.blast_dir)
    
    # Test cases
    test_cases = [
        ('AP027280.AP027280.tlosatx.n1.out', 'AP027280.AP027280.tblastx.n1.out'),
        ('AP027280.AP027280.tlosatx.n8.out', 'AP027280.AP027280.tblastx.n8.out'),
        ('MjeNMV.MelaMJNV.tlosatx.n8.out', 'MjeNMV.MelaMJNV.tblastx.n8.out'),
        ('MelaMJNV.PemoMJNVA.tlosatx.n8.out', 'MelaMJNV.PemoMJNVA.tblastx.n8.out'),
        ('PemoMJNVA.PeseMJNV.tlosatx.n8.out', 'PemoMJNVA.PeseMJNV.tblastx.n8.out'),
        ('PeseMJNV.PemoMJNVB.tlosatx.n8.out', 'PeseMJNV.PemoMJNVB.tblastx.n8.out'),
        ('AP027131.AP027133.tlosatx.n8.out', 'AP027131.AP027133.tblastx.n8.out'),
    ]
    
    all_results = []
    
    for losat_file, blast_file in test_cases:
        losat_path = losat_dir / losat_file
        blast_path = blast_dir / blast_file
        
        if not losat_path.exists():
            print(f"Warning: LOSAT file not found: {losat_path}", file=sys.stderr)
            continue
        if not blast_path.exists():
            print(f"Warning: BLAST file not found: {blast_path}", file=sys.stderr)
            continue
        
        compare_results(losat_path, blast_path)
        all_results.append((losat_path, blast_path))
    
    if not all_results:
        print("No matching files found for comparison.", file=sys.stderr)
        return 1
    
    return 0

if __name__ == '__main__':
    sys.exit(main())

