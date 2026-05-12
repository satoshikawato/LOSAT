use LOSAT::algorithm::tblastx::lookup::compressed::{
    build_blosum62_compressed_lookup, build_blosum62_compressed_score_matrix, compressed_pv_test,
    compute_compressed_index, compute_encoded_compressed_index, BlastCompressedAaLookupTable,
};
use LOSAT::utils::matrix::ncbistdaa;

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:782-849
// ```c
// case 5:
//     new_cell->query_offsets[0] = backbone_cell->payload.query_offsets[2];
//     new_cell->query_offsets[1] = backbone_cell->payload.query_offsets[3];
//     new_cell->query_offsets[2] = query_offset;
// ...
// default:
//     cell_index = (num_entries -3) & COMPRESSED_HITS_CELL_MASK ;
//     if (cell_index == 0 ) {
//         new_cell->next = backbone_cell->payload.overflow_list.head;
//         backbone_cell->payload.overflow_list.head = new_cell;
//     }
// ```
#[test]
fn compressed_word_hit_overflow_matches_ncbi_scan_order() {
    let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
    for offset in 10..18 {
        lookup.add_word_hit(42, offset);
    }

    assert_eq!(lookup.backbone[42].num_used, 8);
    assert_eq!(lookup.overflow_cells.len(), 2);
    assert_eq!(
        lookup.query_offsets_in_scan_order(42),
        vec![10, 11, 12, 17, 13, 14, 15, 16]
    );
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1182-1242
// ```c
// PV_SET(pv, i, pv_array_bts);
// longest_chain = MAX(count, longest_chain);
// ```
#[test]
fn compressed_finalize_sets_pv_and_longest_chain() {
    let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
    lookup.add_word_hit(7, 100);
    lookup.add_word_hit(33, 200);
    lookup.add_word_hit(33, 201);
    lookup.finalize();

    assert_eq!(lookup.longest_chain, 2);
    assert!(compressed_pv_test(
        &lookup.pv,
        7,
        lookup.pv_array_bts as usize
    ));
    assert!(compressed_pv_test(
        &lookup.pv,
        33,
        lookup.pv_array_bts as usize
    ));
    assert!(!compressed_pv_test(
        &lookup.pv,
        34,
        lookup.pv_array_bts as usize
    ));
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:859-929
// ```c
// s_CompressedLookupAddEncoded(lookup, subject_word, query_offset);
// ...
// if (skip == 0)
//     s_CompressedLookupAddWordHit(lookup, index, query_offset);
// ```
#[test]
fn compressed_encoded_and_unencoded_indices_match() {
    let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
    let word = [
        ncbistdaa::S,
        ncbistdaa::T,
        ncbistdaa::A,
        ncbistdaa::G,
        ncbistdaa::W,
    ];
    let encoded_word: Vec<u8> = word
        .iter()
        .map(|&residue| lookup.compress_table[residue as usize])
        .collect();
    let encoded_index = compute_encoded_compressed_index(5, &encoded_word, 15);
    let (unencoded_index, skip) =
        compute_compressed_index(5, &word, 15, &lookup.scaled_compress_table);

    assert_eq!(skip, 0);
    assert_eq!(encoded_index as i32, unencoded_index);
    lookup.add_unencoded(&word, 321);
    assert_eq!(lookup.query_offsets_in_scan_order(encoded_index), vec![321]);
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1071-1176
// ```c
// s_CompressedAddWordHits(&info, query->sequence, offset);
// ...
// s_CompressedAddWordHitsCore(info, score, 0);
// ```
#[test]
fn compressed_neighbor_lookup_builds_representative_blosum62_word() {
    let query = [
        ncbistdaa::S,
        ncbistdaa::T,
        ncbistdaa::A,
        ncbistdaa::G,
        ncbistdaa::W,
    ];
    let lookup =
        build_blosum62_compressed_lookup(5, 19.3, &query, &[(0, 4)]).expect("compressed lookup");
    let matrix = build_blosum62_compressed_score_matrix(15).expect("compressed matrix");
    let total_hits: i32 = lookup.backbone.iter().map(|cell| cell.num_used).sum();

    assert_eq!(lookup.neighbor_matches, 285);
    assert_eq!(lookup.exact_matches, 1);
    assert_eq!(lookup.longest_chain, 1);
    assert_eq!(total_hits, 286);
    assert_eq!(matrix[ncbistdaa::W as usize][14], 1146);
}
