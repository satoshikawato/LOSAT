#!/usr/bin/env python3
"""
Analyze debug output from LOSAT and NCBI BLAST+ to identify root cause of excessive hits.

This script parses debug logs and compares:
1. Cutoff values (cutoff_score_max, final_cutoff, linking cutoffs)
2. Effective search space (eff_searchsp)
3. HSP saving statistics
4. Linking filter statistics
5. Score distributions

Usage:
    python analyze_debug_output.py <losat_debug.log> [ncbi_debug.log]
"""

import sys
import re
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

def parse_losat_debug(log_file: str) -> Dict:
    """Parse LOSAT debug output."""
    data = {
        'cutoff_calc': [],
        'cutoff_update': [],
        'linking_cutoff': [],
        'linking_filter': [],
        'hsp_saving': [],
    }
    
    with open(log_file, 'r') as f:
        for line in f:
            # Parse cutoff calculation
            if '[DEBUG CUTOFF_CALC]' in line:
                match = re.search(r'eff_searchsp=(\d+)', line)
                if match:
                    data['cutoff_calc'].append({
                        'eff_searchsp': int(match.group(1)),
                        'line': line.strip()
                    })
            
            # Parse cutoff update
            if '[DEBUG CUTOFF_UPDATE]' in line:
                match = re.search(r'final_cutoff=(\d+)', line)
                if match:
                    data['cutoff_update'].append({
                        'final_cutoff': int(match.group(1)),
                        'line': line.strip()
                    })
            
            # Parse linking cutoff
            if '[DEBUG LINKING_CUTOFF]' in line:
                match = re.search(r'final_cutoff_small=(\d+).*final_cutoff_big=(\d+)', line)
                if match:
                    data['linking_cutoff'].append({
                        'cutoff_small': int(match.group(1)),
                        'cutoff_big': int(match.group(2)),
                        'line': line.strip()
                    })
            
            # Parse linking filter statistics
            if '[DEBUG LINKING_FILTER]' in line:
                match = re.search(r'INDEX 0.*filter_rate=([\d.]+)%', line)
                if match:
                    data['linking_filter'].append({
                        'index0_filter_rate': float(match.group(1)),
                        'line': line.strip()
                    })
            
            # Parse HSP saving statistics
            if '[DEBUG HSP_SAVING]' in line:
                match = re.search(r'saved=(\d+).*filtered_by_cutoff=(\d+)', line)
                if match:
                    data['hsp_saving'].append({
                        'saved': int(match.group(1)),
                        'filtered_by_cutoff': int(match.group(2)),
                        'line': line.strip()
                    })
    
    return data

def parse_ncbi_debug(log_file: str) -> Dict:
    """Parse NCBI BLAST+ debug output."""
    data = {
        'effective_search_space': None,
        'length_adjustment': None,
    }
    
    with open(log_file, 'r') as f:
        for line in f:
            # Extract effective search space
            match = re.search(r'Effective search space:\s*([\d.e+-]+)', line, re.IGNORECASE)
            if match:
                data['effective_search_space'] = float(match.group(1))
            
            # Extract length adjustment
            match = re.search(r'Length adjustment:\s*([\d.e+-]+)', line, re.IGNORECASE)
            if match:
                data['length_adjustment'] = float(match.group(1))
    
    return data

def compare_cutoffs(losat_data: Dict, ncbi_data: Dict):
    """Compare cutoff values between LOSAT and NCBI."""
    print("=" * 60)
    print("Cutoff Comparison")
    print("=" * 60)
    
    if losat_data['cutoff_calc']:
        print("\nLOSAT Cutoff Calculation:")
        for item in losat_data['cutoff_calc']:
            print(f"  eff_searchsp: {item['eff_searchsp']:,}")
    
    if losat_data['cutoff_update']:
        print("\nLOSAT Final Cutoff:")
        for item in losat_data['cutoff_update']:
            print(f"  final_cutoff: {item['final_cutoff']}")
    
    if ncbi_data.get('effective_search_space'):
        print(f"\nNCBI Effective Search Space: {ncbi_data['effective_search_space']:,.0f}")
    
    # Compare eff_searchsp
    if losat_data['cutoff_calc'] and ncbi_data.get('effective_search_space'):
        losat_eff = losat_data['cutoff_calc'][0]['eff_searchsp']
        ncbi_eff = ncbi_data['effective_search_space']
        diff = abs(losat_eff - ncbi_eff)
        ratio = losat_eff / ncbi_eff if ncbi_eff > 0 else 0
        print(f"\nComparison:")
        print(f"  LOSAT eff_searchsp: {losat_eff:,}")
        print(f"  NCBI eff_searchsp: {ncbi_eff:,.0f}")
        print(f"  Difference: {diff:,.0f} ({ratio:.2f}x)")

def compare_hsp_statistics(losat_data: Dict):
    """Analyze HSP saving and filtering statistics."""
    print("\n" + "=" * 60)
    print("HSP Saving Statistics")
    print("=" * 60)
    
    if losat_data['hsp_saving']:
        for item in losat_data['hsp_saving']:
            saved = item['saved']
            filtered = item['filtered_by_cutoff']
            total = saved + filtered
            if total > 0:
                filter_rate = (filtered / total) * 100
                print(f"\n  Total HSPs: {total:,}")
                print(f"  Saved: {saved:,} ({100-filter_rate:.1f}%)")
                print(f"  Filtered by cutoff: {filtered:,} ({filter_rate:.1f}%)")
    
    if losat_data['linking_filter']:
        print("\n" + "=" * 60)
        print("Linking Filter Statistics")
        print("=" * 60)
        for item in losat_data['linking_filter']:
            print(f"  INDEX 0 filter rate: {item['index0_filter_rate']:.2f}%")

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_debug_output.py <losat_debug.log> [ncbi_debug.log]")
        sys.exit(1)
    
    losat_log = sys.argv[1]
    ncbi_log = sys.argv[2] if len(sys.argv) > 2 else None
    
    print("Analyzing LOSAT debug output...")
    losat_data = parse_losat_debug(losat_log)
    
    ncbi_data = {}
    if ncbi_log:
        print("Analyzing NCBI debug output...")
        ncbi_data = parse_ncbi_debug(ncbi_log)
    
    # Compare cutoffs
    compare_cutoffs(losat_data, ncbi_data)
    
    # Analyze HSP statistics
    compare_hsp_statistics(losat_data)
    
    print("\n" + "=" * 60)
    print("Analysis complete")
    print("=" * 60)

if __name__ == '__main__':
    main()

