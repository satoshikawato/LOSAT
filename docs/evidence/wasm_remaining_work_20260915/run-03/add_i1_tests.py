# NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863,876-885
# H2_index--; if (b0) H2_index=next_larger; if (!(b0|b1|b2)) { ... }
# Boundary cases protect original predecessor selection and exact copied bits.
from pathlib import Path
import subprocess,difflib,json,hashlib
b=Path(__file__).resolve().parent;dest=b/'work/I1-source';p=dest/'LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs';raw=p.read_bytes();s=raw.decode().replace('\r\n','\n')
s=s.replace('A real call per HSP lets Wasm enter optimized scan code during a long group\n// invocation. Native retains inlining.', 'A call per HSP provides a repeated real function-entry boundary within a long\n// group invocation. Native retains inlining.')
test='''
    // NCBI reference: c++/src/algo/blast/core/link_hsps.c:827-863,876-885
    // for (H2_index=H_index-1; H2_index>1;) { ... if(b0) H2_index=next_larger; }
    // The end sentinels 0/1 cannot select a predecessor, even with tracing enabled.
    #[test]
    fn large_gap_scan_preserves_initial_state_and_zero_one_jumps() {
        let initial = (5, 3, -0.0f64, 9usize);
        let empty = scan_large_gap_predecessors(
            2, &[], &[], &[], 10, 10, false,
            initial.0, initial.1, initial.2, initial.3,
        );
        assert_eq!((empty.0, empty.1, empty.2.to_bits(), empty.3),
                   (initial.0, initial.1, initial.2.to_bits(), initial.3));
        for jump in [0, 1] {
            let helper = LhHelper {
                hsp_idx: SENTINEL_IDX, q_off_trim: 20, s_off_trim: 20,
                sum: [0, 5], next_larger: jump, maxsum1: 0,
            };
            let helpers = vec![helper; 3];
            for trace in [false, true] {
                let result = scan_large_gap_predecessors(
                    3, &helpers, &[], &[], 10, 10, trace,
                    initial.0, initial.1, initial.2, initial.3,
                );
                assert_eq!((result.0, result.1, result.2.to_bits(), result.3),
                           (initial.0, initial.1, initial.2.to_bits(), initial.3));
            }
        }
    }

    // NCBI reference: c++/src/algo/blast/core/link_hsps.c:838-859
    // b0=sum<=H_hsp_sum; b1=q_off_t<=H_query_etrim; b2=s_off_t<=H_sub_etrim;
    // if (!(b0|b1|b2)) { H2=H2_helper->ptr; H_hsp_xsum=H2->hsp_link.xsum[index]; }
    // Strict ties retain the first visited HSP, using the helper's mapped index.
    #[test]
    fn large_gap_scan_preserves_strict_ties_mapping_and_xsum_bits() {
        let (mut links, _) = tree_links(&[[0, 10], [0, 10]], &[true, true]);
        links[0].num[1] = 7;
        links[1].num[1] = 8;
        links[0].xsum[1] = f64::from_bits(1.0f64.to_bits() + 1);
        links[1].xsum[1] = f64::from_bits(2.0f64.to_bits() + 1);
        let hits = vec![mock_hit(0, 1, 0, 1, 0), mock_hit(2, 3, 2, 3, 1)];
        let helper = LhHelper {
            hsp_idx: 1, q_off_trim: 20, s_off_trim: 20,
            sum: [0, 10], next_larger: 0, maxsum1: 0,
        };
        let mut helpers = vec![helper; 4];
        helpers[3].hsp_idx = 0;
        helpers[3].next_larger = 2;
        for (query_offset, subject_offset, expected) in [(20, 20, 0), (10, 20, 1), (20, 10, 1)] {
            helpers[3].q_off_trim = query_offset;
            helpers[3].s_off_trim = subject_offset;
            for trace in [false, true] {
                let result = scan_large_gap_predecessors(
                    4, &helpers, &links, &hits, 10, 10, trace, 0, 0, 0.0, SENTINEL_IDX,
                );
                assert_eq!((result.0, result.1, result.2.to_bits(), result.3),
                    (links[expected].sum[1], links[expected].num[1],
                     links[expected].xsum[1].to_bits(), expected));
            }
        }
    }
'''
pos=s.rfind('\n}');s=s[:pos]+test+s[pos:];p.write_text(s);subprocess.run(['rustfmt',str(p)],check=True)
zz=p.read_text().splitlines(True);aa=raw.decode().splitlines(True);out=[]
for tag,i,j,k,l in difflib.SequenceMatcher(None,[x.rstrip('\r\n') for x in aa],[x.rstrip('\r\n') for x in zz],autojunk=False).get_opcodes():out.extend(aa[i:j] if tag=='equal' else zz[k:l] if tag in ['replace','insert'] else [])
p.write_bytes(''.join(out).encode())
(b/'I1-source-manifest.json').write_text(json.dumps({str(p.relative_to(dest)):hashlib.sha256(p.read_bytes()).hexdigest() for p in dest.rglob('*') if p.is_file()},indent=2)+'\n')
base=b/'work/baseline-source/LOSAT/src/algorithm/tblastx/sum_stats_linking/linking.rs';(b/'I1-vs-baseline.diff').write_text(''.join(difflib.unified_diff(base.read_text().splitlines(True),p.read_text().splitlines(True),fromfile='baseline/linking.rs',tofile='I1/linking.rs')))
