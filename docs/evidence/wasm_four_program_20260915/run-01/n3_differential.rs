
// NCBI reference: c++/src/algo/blast/core/na_ungapped.c:262-349
// Uint1 q_byte = (q[0] << 6) | (q[1] << 4) | (q[2] << 2) | q[3];
// if (score >= reduced_cutoff) { s_NuclUngappedExtendExact(...); }
// Differential-only test: frozen pre-optimization loop vs candidate byte lookup.
#[cfg(test)] mod n3_differential {
use super::*;
fn baseline_approx(
    q_seq: &[u8],        // BLASTNA (1 byte/base)
    s_seq_packed: &[u8], // ncbi2na packed (4 bases/byte)
    q_off: usize,
    s_off: usize,
    s_match_end: usize,
    subject_len: usize,
    x_dropoff: i32,
    score_table: &[i32; 256],
    reduced_cutoff: i32,
    matrix: &[i32; BLASTNA_SIZE * BLASTNA_SIZE],
) -> UngappedData {
    let x = -x_dropoff;
    let q_len = q_seq.len();

    let len = (COMPRESSION_RATIO - (s_off % COMPRESSION_RATIO)) % COMPRESSION_RATIO;
    let q_ext = q_off + len;
    let s_ext = s_off + len;

    let mut q_idx = q_ext as isize;
    let mut s_idx = (s_ext / COMPRESSION_RATIO) as isize;
    let left_len = (q_ext.min(s_ext) / COMPRESSION_RATIO) as usize;

    let mut score = 0i32;
    let mut sum = 0i32;
    let mut new_q = q_idx;

    // Left extension in 4-base blocks.
    for _ in 0..left_len {
        if q_idx < 4 || s_idx == 0 {
            break;
        }
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:292-304
        // ```c
        // Uint1 s_byte = s[-1];
        // Uint1 q_byte = (q[-4] << 6) | (q[-3] << 4) | (q[-2] << 2) | q[-1];
        // ```
        // SAFETY: q_idx >= 4 and s_idx > 0 are enforced by the loop guard.
        let s_byte = unsafe { *s_seq_packed.get_unchecked((s_idx - 1) as usize) };
        let q_byte = unsafe {
            (*q_seq.get_unchecked((q_idx - 4) as usize) << 6)
                | (*q_seq.get_unchecked((q_idx - 3) as usize) << 4)
                | (*q_seq.get_unchecked((q_idx - 2) as usize) << 2)
                | *q_seq.get_unchecked((q_idx - 1) as usize)
        };
        sum += score_table[(q_byte ^ s_byte) as usize];
        if sum > 0 {
            new_q = q_idx - 4;
            score += sum;
            sum = 0;
        }
        if sum < x {
            break;
        }
        q_idx -= 4;
        s_idx -= 1;
    }

    let q_start = new_q.max(0) as usize;
    let s_start = s_ext - (q_ext - q_start);

    // Right extension in 4-base blocks.
    let mut q_idx = q_ext as isize;
    let mut s_idx = (s_ext / COMPRESSION_RATIO) as isize;
    let right_len = ((q_len - q_ext).min(subject_len - s_ext) / COMPRESSION_RATIO) as usize;
    sum = 0;
    new_q = q_idx;

    for _ in 0..right_len {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/na_ungapped.c:322-334
        // ```c
        // Uint1 s_byte = s[0];
        // Uint1 q_byte = (q[0] << 6) | (q[1] << 4) | (q[2] << 2) | q[3];
        // ```
        // SAFETY: q_idx + 3 and s_idx stay within bounds by the loop length
        // derived from query/subject lengths (NCBI uses the same limits).
        let s_byte = unsafe { *s_seq_packed.get_unchecked(s_idx as usize) };
        let q_byte = unsafe {
            (*q_seq.get_unchecked(q_idx as usize) << 6)
                | (*q_seq.get_unchecked((q_idx + 1) as usize) << 4)
                | (*q_seq.get_unchecked((q_idx + 2) as usize) << 2)
                | *q_seq.get_unchecked((q_idx + 3) as usize)
        };
        sum += score_table[(q_byte ^ s_byte) as usize];
        if sum > 0 {
            new_q = q_idx + 3;
            score += sum;
            sum = 0;
        }
        if sum < x {
            break;
        }
        q_idx += 4;
        s_idx += 1;
    }

    if score >= reduced_cutoff {
        return extend_hit_ungapped_exact_ncbi(
            q_seq,
            s_seq_packed,
            q_off,
            s_off,
            subject_len,
            x_dropoff,
            matrix,
        );
    }

    let s_match_len = s_match_end.saturating_sub(s_start);
    let right_len = if new_q >= q_start as isize {
        (new_q - q_start as isize + 1) as usize
    } else {
        0
    };
    let length = s_match_len.max(right_len);

    UngappedData {
        q_start,
        s_start,
        length,
        score,
    }
}
#[test] fn four_base_extension_matches_frozen_loop_at_cutoff_and_context_boundaries() {
 let table=build_nucl_score_table(2,-3);
 let matrix=crate::algorithm::blastn::alignment::build_blastna_matrix(2,-3);
 let mut random=0x20260915u32;
 let mut cases=0;
 for len in [16,17,31,32,33,63,64,65] {
  for masked in [false,true] {
   let query:Vec<u8>=(0..len).map(|i| {random=random.wrapping_mul(1664525).wrapping_add(1013904223);if masked && i%7==0 {((random>>8)%16) as u8}else{(i%4) as u8}}).collect();
   let subject:Vec<u8>=(0..len).map(|i| if i%13==0 {2}else{(i%4) as u8}).collect();
   let mut packed=vec![0u8; (len+3)/4];for (i,&value) in subject.iter().enumerate(){packed[i/4]|=value << (6-2*(i%4));}
   for phase in 0..4 {
    let mut concat=vec![15;phase];concat.extend_from_slice(&query);concat.extend_from_slice(&[15;4]);
    let all_bytes=build_query_four_base_bytes(&concat);
    let local_bytes=&all_bytes[phase..phase+len-3];
    for qoff in 0..len-11 {for soff in 0..len-11 {
     let approx=baseline_approx(&query,&packed,qoff,soff,soff+11,len,20,&table,i32::MAX,&matrix);
     for cutoff in [0,approx.score-1,approx.score,approx.score+1,i32::MAX] {
      for (q,bytes,offset) in [(&query[..],local_bytes,qoff),(&concat[..],&all_bytes[..],phase+qoff)] {
       let old=baseline_approx(q,&packed,offset,soff,soff+11,len,20,&table,cutoff,&matrix);
       let new=extend_hit_ungapped_approx_ncbi(q,bytes,&packed,offset,soff,soff+11,len,20,&table,cutoff,&matrix);
       assert_eq!((old.q_start,old.s_start,old.length,old.score),(new.q_start,new.s_start,new.length,new.score));cases+=1;
      }
     }
    }}
   }
  }
 }
 println!("N3_EXTENSION_COMPARISONS={cases}");
}
}
