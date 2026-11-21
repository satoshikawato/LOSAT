//! NCBI's engine error of a query batch whose `-matrix` the score block cannot load, shared by the
//! ungapped searches of BLASTX and BLASTP (one owner of the text: `unknown_matrix_message_of`).

/// NCBI's error of a query batch whose score block cannot load the `-matrix` of an ungapped search
/// (BLASTX's `blastx::query_setup::unknown_matrix_message`, BLASTP's `unknown_matrix_message`), for
/// queries given by the messages that each has besides the error (BLASTP's queries hold the O
/// warning of a protein reader as well as `Sequence contains no data`, and its records are
/// another type): each query's error, then its messages as warnings,
/// sorted and without repeats as NCBI keeps them (`TSearchMessages::RemoveDuplicates`,
/// `report::query_warnings::sorted_query_messages`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:158-182 (cited at
/// `blastx::query_setup::unknown_matrix_message`: the text is `search_messages.ToString()`, every
/// message of each query)
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:919-926
pub(crate) fn unknown_matrix_message_of(messages: impl IntoIterator<Item = Vec<String>>) -> String {
    let mut message = String::new();
    for query in messages {
        message.push_str("Error: Unknown error code -1 ");
        for warning in crate::report::query_warnings::sorted_query_messages(&query) {
            message.push_str("Warning: ");
            message.push_str(warning);
            message.push(' ');
        }
    }
    message
}
