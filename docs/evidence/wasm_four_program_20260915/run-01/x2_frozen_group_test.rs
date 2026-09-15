
    // NCBI c++/src/algo/blast/core/link_hsps.c:827-895,955-980:
    // b0 = sum <= H_hsp_sum; H2_index = next_larger;
    // H_hsp_link = H2; H->hsp->evalue = evalue;
    #[test]
    fn x2_frozen_group_transcript() {
        let params = KarlinParams::default();
        let contexts: Vec<_> = (0..3).map(|index| QueryContext {
            q_idx: 0, f_idx:index as u8, frame:index as i8 + 1,
            aa_seq: Vec::new(), aa_seq_nomask:None, aa_len:2000,
            orig_len:6000, frame_base:(index*2001) as i32,
            is_valid:true, karlin_params:params,
        }).collect();
        let mut helpers=Vec::new(); let mut links=Vec::new();
        for case in 0..128usize {
            let count = [8,16,17,32][case%4];
            let mut hits:Vec<_>=(0..count).map(|index| {
                let q=(count-index)*[7,10,13,30][(case/4)%4];
                let s=if case%3==0 {(count-index)*[7,10,13][case%3]} else {q+(index%3)};
                let len=[8,12,20][(case+index)%3];
                let score=[19,20,21,35,35,40,41][(index+case/4)%7];
                let mut h=mock_hit(q,q+len,s,s+len,score);
                h.ctx_idx=(index+case)%3; h.q_frame=contexts[h.ctx_idx].frame;
                h.q_orig_len=6000;h.s_orig_len=6000;h.hsp_list_order=index; h
            }).collect();
            hits.sort_by(|a,b| b.q_aa_start.cmp(&a.q_aa_start).then_with(||b.q_aa_end.cmp(&a.q_aa_end)).then_with(||b.s_aa_start.cmp(&a.s_aa_start)).then_with(||b.s_aa_end.cmp(&a.s_aa_end)));
            let cutoff=LinkHspCutoffs {cutoff_small_gap:if case%2==0 {20}else{0},cutoff_big_gap:20,gap_prob:0.5,ignore_small_gaps:case%2!=0};
            println!("STATE CASE {case}");
            let result=link_hsp_group_ncbi(hits,&params,&cutoff,0.5,false,6000,&contexts,&[0;6],&[0;3],&[4_000_000;3],&[params.k.ln();3],&mut helpers,&mut links);
            for (index,h) in result.iter().enumerate(){println!("STATE RESULT {index} {h:?} {}",h.e_value.to_bits());}
            for (index,h) in links.iter().enumerate(){println!("STATE LINK {index} {:?} {:?} {:?} {:?} {} {}",h.sum,h.num,h.link,h.xsum.map(f64::to_bits),h.linked_set,h.start_of_chain);}
        }
    }
