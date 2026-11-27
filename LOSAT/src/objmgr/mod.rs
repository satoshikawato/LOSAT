//! NCBI's object manager (`c++/src/objmgr/`) as LOSAT's programs meet it when they add a Seq-entry
//! of their input to a scope: psiblast and tblastn add the query Seq-entry of `-in_pssm`. The index
//! of the entry's annotations as the entry joins the scope, with its messages, exceptions and
//! LOSAT's refusals, is `annot_index`.

pub(crate) mod annot_index;
