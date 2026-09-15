// NCBI reference: c++/src/algo/blast/core/greedy_align.c:71-77,523-526,571-576,666-678
// space->space_used = 0;
// last_seq2_off[d - 1][diag_lower-1] = kInvalidOffset;
// if (fence_hit && *fence_hit) { return 0; }
// last_seq2_off[d + 1] = (Int4*) s_GetMBSpace(mem_pool, ...);
// Diagnostic-only transcript compares every retained scratch cell after each call.
#[cfg(test)]
mod four_scratch_transcript {
    use super::*;
    #[test]
    fn four_scratch_transcript() {
        let mut mem = GreedyNonAffineMem::new();
        let mut random = 0x20260915u32;
        for iteration in 0..240 {
            let qlen = [0,1,2,3,4,7,8,15,16,31,32,63,64,65][iteration % 14];
            let slen = [65,64,63,32,31,16,15,8,7,4,3,2,1,0][iteration % 14];
            let query: Vec<u8> = (0..qlen).map(|i| (i % 4) as u8).collect();
            let mut subject: Vec<u8> = (0..slen).map(|i| {
                random = random.wrapping_mul(1664525).wrapping_add(1013904223);
                if random % 7 == 0 { ((random >> 8) % 16) as u8 } else {(i % 4) as u8}
            }).collect();
            if iteration % 13 == 0 && slen > 2 {subject[slen / 2] = FENCE_SENTRY;}
            for reverse in [false,true] {
                for traceback in [false,true] {
                    for max_dist in [1,2,4,16,64] {
                        let mut edit = GapPrelimEditBlock::new();
                        let mut seed = GreedySeed::default();
                        let (mut q, mut s, mut fence) = (0,0,false);
                        let result=blast_greedy_align(&query,qlen as i32,&subject,slen as i32,reverse,20,2,4,&mut q,&mut s,max_dist,&mut mem,traceback.then_some(&mut edit),4,&mut fence,&mut seed);
                        let mut hash=0xcbf29ce484222325u64;
                        let mut absorb=|v:u64| {hash=(hash ^ v).wrapping_mul(0x100000001b3);};
                        for row in [&mem.last_seq2_off[0],&mem.last_seq2_off[1],&mem.max_score,&mem.traceback_pool] {
                            absorb(row.len() as u64);for &v in row {absorb(v as u64);}
                        }
                        absorb(mem.traceback_used as u64);
                        absorb(edit.num_ops as u64);
                        for op in &edit.edit_ops {absorb(op.num as u64);absorb(op.op_type as u64);}
                        println!("SCRATCH {iteration} {reverse} {traceback} {max_dist} {result} {q} {s} {fence} {} {} {} {hash}",seed.start_q,seed.start_s,seed.match_length);
                    }
                }
            }
        }
    }
}
