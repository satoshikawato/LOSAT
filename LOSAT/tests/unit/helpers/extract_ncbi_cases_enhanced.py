#!/usr/bin/env python3
"""NCBI BLAST出力からテストケースを抽出（改訂版）

重要な改善点:
1. ヘッダーから実効探索空間を抽出
2. -comp_based_stats 0 の使用を推奨
3. 複数のパラメータセットに対応
4. 様々な長さの配列に対応

使用方法:
    python extract_ncbi_cases_enhanced.py <ncbi_blast_output.txt> <query_length> <db_length> <db_num_seqs> [blastn|tblastx]
"""

import sys
import re
from typing import List, Dict, Optional

def extract_effective_search_space(header_lines: List[str]) -> Optional[float]:
    """NCBI BLASTヘッダーから実効探索空間を抽出"""
    for line in header_lines:
        # "Effective search space used: 1234567890" の形式を探す
        match = re.search(r'[Ee]ffective\s+search\s+space[:\s]+([\d.]+)', line)
        if match:
            try:
                return float(match.group(1))
            except ValueError:
                pass
    return None

def parse_blast_output(filename: str) -> tuple[List[Dict], Optional[float]]:
    """NCBI BLAST出力をパース（ヘッダー情報も含む）"""
    cases = []
    header_lines = []
    effective_space = None
    
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                header_lines.append(line)
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            
            try:
                qseqid, sseqid, pident, length, mismatch, gapopen, \
                qstart, qend, sstart, send, evalue_str, bitscore_str = fields[:12]
                
                length = int(length)
                bitscore = float(bitscore_str)
                
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
    
    # ヘッダーから実効探索空間を抽出
    effective_space = extract_effective_search_space(header_lines)
    if effective_space:
        print(f"Found effective search space in header: {effective_space}", file=sys.stderr)
    else:
        print("Warning: Could not find effective search space in header", file=sys.stderr)
    
    return cases, effective_space

def calculate_raw_score_from_bit_score(bit_score: float, lambda_val: float, k: float) -> int:
    """Bit scoreからraw scoreを逆算"""
    import math
    ln2 = math.log(2.0)
    raw_score = (bit_score * ln2 + math.log(k)) / lambda_val
    return int(raw_score)

def generate_rust_code(cases: List[Dict], mode: str, q_len: int, db_len: int, 
                       db_num_seqs: int, effective_space: Optional[float]):
    """Rustテストケースコードを生成（実効探索空間を含む）"""
    
    if mode == 'blastn':
        params = {
            'lambda': 1.28,
            'k': 0.46,
            'h': 0.85,
            'alpha': 1.5,
            'beta': -2.0,
        }
    else:  # tblastx
        params = {
            'lambda': 0.267,
            'k': 0.041,
            'h': 0.14,
            'alpha': 1.9,
            'beta': -30.0,
        }
    
    print("// Auto-generated test cases from NCBI BLAST output")
    print("// ⚠️ 重要: 以下のオプションで実行したデータを使用してください")
    print("//")
    if mode == 'blastn':
        print("// BLASTN: -dust no -reward 1 -penalty -2 -gapopen 0 -gapextend 0")
    else:
        print("// TBLASTX: -comp_based_stats 0 -seg no -gapopen 11 -gapextend 1 -matrix BLOSUM62")
    print("//")
    print("// 注意:")
    print("// - 実効探索空間のテストは1%の許容誤差を使用します（完全一致は不可能なため）")
    print("// - E-valueテストでは固定アライメント座標を使用してください（X-dropの影響を避けるため）")
    print()
    
    if effective_space:
        print(f"// Expected effective search space from NCBI BLAST header: {effective_space}")
        print()
    
    for i, case in enumerate(cases[:10]):  # 最初の10ケース
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
            print(f"            aln_len: 0,")
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
            if effective_space and i == 0:  # 最初のケースにのみ実効探索空間を追加
                print(f"            expected_effective_space: Some({effective_space:.1}),")
                print(f"            // ⚠️ 注意: 実効探索空間のテストは1%の許容誤差を使用します")
                print(f"            // NCBI BLASTの実効長計算は非常に複雑なため、完全一致は不可能です")
                print(f"            expected_length_adjustment: None, // Calculate or extract from NCBI BLAST")
            else:
                print(f"            expected_effective_space: None,")
                print(f"            expected_length_adjustment: None,")
            print(f"            tolerance: 0.1,")
            print(f"            test_name: \"{mode}_case_{i+1}\".to_string(),")
            print(f"        }},")
        else:  # tblastx
            print(f"        // Test case {i+1}: {case['qseqid']} vs {case['sseqid']}")
            print(f"        // Alignment length: {case['alignment_length']} aa, Identity: {case['identity']:.2f}%")
            print(f"        // ⚠️ 確認: -comp_based_stats 0 で実行しましたか？")
            print(f"        NcbiEvalueTestCase {{")
            print(f"            score: {raw_score},")
            print(f"            q_len: 0,")
            print(f"            db_len: 0,")
            print(f"            db_num_seqs: 0,")
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
            print(f"            expected_effective_space: None,")
            print(f"            expected_length_adjustment: None,")
            print(f"            tolerance: 0.1,")
            print(f"            test_name: \"{mode}_case_{i+1}\".to_string(),")
            print(f"        }},")
        print()

def main():
    if len(sys.argv) < 6:
        print("Usage: python extract_ncbi_cases_enhanced.py <ncbi_output.txt> <q_len> <db_len> <db_num_seqs> <blastn|tblastx>")
        print()
        print("⚠️  重要: 正解データ生成時の必須オプション")
        print("  BLASTN:")
        print("    -dust no  # フィルタリングOFF（必須）")
        print("    -reward 1 -penalty -2  # パラメータ明示（必須）")
        print("    -gapopen 0 -gapextend 0  # ギャップコスト明示（必須）")
        print()
        print("  TBLASTX:")
        print("    -comp_based_stats 0  # 組成補正OFF（必須）")
        print("    -seg no  # フィルタリングOFF（必須）")
        print("    -gapopen 11 -gapextend 1  # ギャップコスト明示（必須）")
        print("    -matrix BLOSUM62  # スコアリング行列明示（必須）")
        print()
        print("  - ヘッダーから実効探索空間を自動抽出します")
        print("  - 実効探索空間のテストは1%の許容誤差を使用します")
        sys.exit(1)
    
    filename = sys.argv[1]
    q_len = int(sys.argv[2])
    db_len = int(sys.argv[3])
    db_num_seqs = int(sys.argv[4])
    mode = sys.argv[5]
    
    if mode not in ['blastn', 'tblastx']:
        print("Error: mode must be 'blastn' or 'tblastx'", file=sys.stderr)
        sys.exit(1)
    
    cases, effective_space = parse_blast_output(filename)
    print(f"// Extracted {len(cases)} test cases from {filename}", file=sys.stderr)
    
    if mode == 'tblastx' and not any('comp_based_stats' in line or 'composition' in line.lower() 
                                      for line in open(filename) if line.startswith('#')):
        print("⚠️  Warning: TBLASTXデータですが、-comp_based_stats 0 の使用を確認できませんでした", 
              file=sys.stderr)
        print("   純粋な統計値と比較するため、-comp_based_stats 0 で再実行することを推奨します", 
              file=sys.stderr)
    
    generate_rust_code(cases, mode, q_len, db_len, db_num_seqs, effective_space)

if __name__ == '__main__':
    main()

