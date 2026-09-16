// NCBI reference: c++/src/algo/blast/core/link_hsps.c:603-688,827-895,955-980
// H2_index = next_larger; H_hsp_link = H2; H->hsp->evalue = evalue;
// Diagnostic-only transcript: preserve all floating-point bits and typed indices.
fn state_index(value: usize) -> Option<usize> {
    if value == SENTINEL_IDX { None } else { Some(value) }
}

fn state_snapshot(event: &str, helpers: &[LhHelper], links: &[HspLink], hits: &[UngappedHit], tree: &DualMaximumTree) {
    println!("STATE SNAP {event}");
    for (i,h) in helpers.iter().enumerate() {
        println!("STATE HELP {i} {:?}", (state_index(h.hsp_idx),h.q_off_trim,h.s_off_trim,h.sum,h.next_larger,h.maxsum1));
    }
    for (i,h) in links.iter().enumerate() {
        println!("STATE LINK {i} {:?}", ((h.score,h.ctx_idx,h.q_off_trim,h.s_off_trim,h.q_end_trim,h.s_end_trim),h.sum,h.xsum.map(f64::to_bits),h.num,h.link.map(state_index),(h.changed,h.linked_to,h.start_of_chain,h.linked_set),state_index(h.next_active),state_index(h.prev_active)));
    }
    for (i,h) in hits.iter().enumerate() {
        println!("STATE HIT {i} {h:?} E_BITS {}",h.e_value.to_bits());
    }
    println!("STATE TREE {} {:?}",tree.capacity,tree.nodes.iter().map(|n|(n.sum,n.index.map(state_index))).collect::<Vec<_>>());
}

// NCBI reference: c++/src/algo/blast/core/link_hsps.c:478-486,553-562,990-994
// link_hsp_array[index]->hsp = hsp_array[index];
// Unique stable identities retain chain ownership across reordered equal HSPs.
pub fn remaining_state_transcript() {
    let params=KarlinParams::default();
    let contexts:Vec<_>=(0..6).map(|index|QueryContext {
        q_idx:0,f_idx:index as u8,frame:if index<3 {index as i8+1}else{2-index as i8},
        aa_seq:Vec::new(),aa_seq_nomask:None,aa_len:2000,orig_len:6000,
        frame_base:(index*2001) as i32,is_valid:true,karlin_params:params,
    }).collect();
    let mut helpers=Vec::new(); let mut links=Vec::new();
    for case in 0..352usize {
        let count=[0,1,2,3,15,16,17,31,32,33,3][case%11];
        let subject_len=if case%16==0 {600001}else{6000};
        let mut hits:Vec<_>=(0..count).map(|index| {
            let j=if case>=176 {index/2*2}else{index};
            let q=(count-j)*[7,10,13,30][(case/4)%4];
            let s=if case%3==0 {(count-j)*[7,10,13][case%3]}else{q+(j%3)};
            let len=[8,12,20][(case+j)%3];
            let score=[19,20,21,35,35,40,41][(j+case/4)%7];
            let mut h=tests::mock_hit(q,q+len,s,s+len,score);
            h.ctx_idx=(index+case)%3+if case%2==0 {0}else{3};
            h.q_frame=contexts[h.ctx_idx].frame;
            h.s_frame=if case%3==0 {-1}else{1};
            h.q_orig_len=6000; h.s_orig_len=subject_len as usize;
            h.hsp_list_order=index;h.link_id=case*64+index+1;h
        }).collect();
        hits.sort_by(|a,b|b.q_aa_start.cmp(&a.q_aa_start).then_with(||b.q_aa_end.cmp(&a.q_aa_end)).then_with(||b.s_aa_start.cmp(&a.s_aa_start)).then_with(||b.s_aa_end.cmp(&a.s_aa_end)));
        let cutoff=LinkHspCutoffs {cutoff_small_gap:if case%2==0 {20}else{0},cutoff_big_gap:20,gap_prob:0.5,ignore_small_gaps:case%2!=0};
        println!("STATE CASE {case} COUNT {count} SUBJECT {subject_len}");
        let result=link_hsp_group_ncbi(hits,&params,&cutoff,0.5,false,subject_len,&contexts,&[0;6],&[0;6],&[4_000_000;6],&[params.k.ln();6],&mut helpers,&mut links);
        for (index,h) in result.iter().enumerate(){println!("STATE RESULT {index} {h:?} E_BITS {}",h.e_value.to_bits());}
    }
}
