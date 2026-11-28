//! NCBI's writer of binary GI lists (`CBinaryListBuilder` of `eGi`,
//! c++/src/objtools/blast/seqdb_writer/writedb.cpp), which `blastdb_aliastool -gi_file_in` and
//! `-gilist` with a text list write, and which NCBI's readers of GI lists read
//! (`seqdb::id_list`, whose first words it shares).

use crate::seqdb::id_list::{GI_LIST_MAGIC, LONG_GI_LIST_MAGIC};

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/writedb.cpp:207-255
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_writer/writedb_convert.hpp:157-200
// An id is checked to be positive only in NCBI's debug builds (`_ASSERT`): a negative or zero
// id is written, and a negative one makes the list one of 8-byte ids (its upper 32 bits are not
// 0). The count is NCBI's `int` of the number of ids. Each word is one `write` to NCBI's stream
// (`s_WriteInt4`, `s_WriteInt8BE`), which decides what its file buffer writes before a closed pipe
// ends the process (blastdb_aliastool's GI lists, `writedb::ofstream::OfStream`).
/// The writes of a binary GI list of `ids` (in any order, repeats kept), one `write` per word:
/// the first word -1, or -2 when an id does not fit 32 bits; the number of ids; the ids sorted,
/// as 4-byte or 8-byte big-endian words.
pub(crate) fn write_gi_list(ids: &[i64], write: &mut dyn FnMut(&[u8])) {
    let eight = ids.iter().any(|&id| id >> 32 != 0);
    let magic = if eight {
        LONG_GI_LIST_MAGIC
    } else {
        GI_LIST_MAGIC
    };
    let mut sorted = ids.to_vec();
    sorted.sort_unstable();
    write(&magic.to_be_bytes());
    write(&(sorted.len() as i32).to_be_bytes());
    for id in sorted {
        if eight {
            write(&id.to_be_bytes());
        } else {
            write(&(id as i32).to_be_bytes());
        }
    }
}

/// The bytes of a binary GI list of `ids` (`write_gi_list`).
#[cfg(test)]
pub(crate) fn gi_list_bytes(ids: &[i64]) -> Vec<u8> {
    let mut out = Vec::new();
    write_gi_list(ids, &mut |word| out.extend_from_slice(word));
    out
}

#[cfg(test)]
#[path = "binary_list_tests.rs"]
mod tests;
