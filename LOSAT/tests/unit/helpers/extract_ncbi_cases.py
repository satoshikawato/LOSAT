#!/usr/bin/env python3
"""NCBI BLAST出力からテストケースを抽出してRustコードを生成するスクリプト

使用方法:
    python extract_ncbi_cases.py <ncbi_blast_output.txt> <query_length> <db_length> <db_num_seqs> [blastn|tblastx]
"""

import sys
import re
from typing import List, Dict

def parse_blast_output(filename: str) -> List[Dict]:
    """NCBI BLAST出力をパース"""
    cases = []
    
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            
            try:
                qseqid, sseqid, pident, length, mismatch, gapopen, \
                qstart, qend, sstart, send, evalue_str, bitscore_str = fields[:12]
                
                # 数値に変換
                length = int(length)
                bitscore = float(bitscore_str)
                
                # E-valueの処理（0.0や科学記法に対応）
                if evalue_str == '0.0':
                    evalue = 0.0
                else:
                    evalue = float(evalue_str)
                
                case = {
                    'qseqid': qseqid,
                    'sseqid': sseqid,
                    'alignment_length': length,
                    'expected_bit_score': bitscore,
                    'expected_evalue': evalue,
                    'identity': float(pident),
                    'mismatch': int(mismatch),
                    'gapopen': int(gapopen),
                }
                cases.append(case)
            except (ValueError, IndexError) as e:
                print(f"Warning: Skipping line due to error: {e}", file=sys.stderr)
                continue
    
    return cases

def calculate_raw_score_from_bit_score(bit_score: float, lambda_val: float, k: float) -> int:
    """Bit scoreからraw scoreを逆算"""
    import math
    ln2 = math.log(2.0)
    raw_score = (bit_score * ln2 + math.log(k)) / lambda_val
    return int(raw_score)

def generate_rust_code(cases: List[Dict], mode: str, q_len: int, db_len: int, db_num_seqs: int):
    """Rustテストケースコードを生成"""
    
    if mode == 'blastn':
        # BLASTN用のパラメータ（megablast default: reward=1, penalty=-2）
        params = {
            'lambda': 1.28,
            'k': 0.46,
            'h': 0.85,
            'alpha': 1.5,
            'beta': -2.0,
        }
    else:  # tblastx
        # TBLASTX用のパラメータ（BLOSUM62 default）
        params = {
            'lambda': 0.267,
            'k': 0.041,
            'h': 0.14,
            'alpha': 1.9,
            'beta': -30.0,
        }
    
    print("// Auto-generated test cases from NCBI BLAST output")
    print("// Add these to get_ncbi_blastn_evalue_cases() or get_ncbi_tblastx_evalue_cases()")
    print()
    
    for i, case in enumerate(cases[:10]):  # 最初の10ケースのみ
        raw_score = calculate_raw_score_from_bit_score(
            case['expected_bit_score'],
            params['lambda'],
            params['k']
        )
        
        if mode == 'blastn':
            print(f"        // Test case {i+1}: {case['qseqid']} vs {case['sseqid']}")
            print(f"        // Alignment length: {case['alignment_length']}, Identity: {case['identity']:.2f}%")
            print(f"        NcbiEvalueTestCase {{")
            print(f"            score: {raw_score},")
            print(f"            q_len: {q_len},")
            print(f"            db_len: {db_len},")
            print(f"            db_num_seqs: {db_num_seqs},")
            print(f"            params: KarlinParams {{")
            print(f"                lambda: {params['lambda']},")
            print(f"                k: {params['k']},")
            print(f"                h: {params['h']},")
            print(f"                alpha: {params['alpha']},")
            print(f"                beta: {params['beta']},")
            print(f"            }},")
            print(f"            expected_bit_score: {case['expected_bit_score']:.1},")
            if case['expected_evalue'] == 0.0:
                print(f"            expected_evalue: 0.0, // NCBI BLAST reports 0.0 for very small values")
            else:
                print(f"            expected_evalue: {case['expected_evalue']:.6e},")
            print(f"            tolerance: 0.1,")
            print(f"        }},")
        else:  # tblastx
            print(f"        // Test case {i+1}: {case['qseqid']} vs {case['sseqid']}")
            print(f"        // Alignment length: {case['alignment_length']} aa, Identity: {case['identity']:.2f}%")
            print(f"        NcbiEvalueTestCase {{")
            print(f"            score: {raw_score},")
            print(f"            aln_len: {case['alignment_length']},")
            print(f"            params: KarlinParams {{")
            print(f"                lambda: {params['lambda']},")
            print(f"                k: {params['k']},")
            print(f"                h: {params['h']},")
            print(f"                alpha: {params['alpha']},")
            print(f"                beta: {params['beta']},")
            print(f"            }},")
            print(f"            expected_bit_score: {case['expected_bit_score']:.1},")
            if case['expected_evalue'] == 0.0:
                print(f"            expected_evalue: 0.0, // NCBI BLAST reports 0.0 for very small values")
            else:
                print(f"            expected_evalue: {case['expected_evalue']:.6e},")
            print(f"            tolerance: 0.1,")
            print(f"        }},")
        print()

def main():
    if len(sys.argv) < 6:
        print("Usage: python extract_ncbi_cases.py <ncbi_output.txt> <q_len> <db_len> <db_num_seqs> <blastn|tblastx>")
        sys.exit(1)
    
    filename = sys.argv[1]
    q_len = int(sys.argv[2])
    db_len = int(sys.argv[3])
    db_num_seqs = int(sys.argv[4])
    mode = sys.argv[5]
    
    if mode not in ['blastn', 'tblastx']:
        print("Error: mode must be 'blastn' or 'tblastx'", file=sys.stderr)
        sys.exit(1)
    
    cases = parse_blast_output(filename)
    print(f"// Extracted {len(cases)} test cases from {filename}", file=sys.stderr)
    
    generate_rust_code(cases, mode, q_len, db_len, db_num_seqs)

if __name__ == '__main__':
    main()


