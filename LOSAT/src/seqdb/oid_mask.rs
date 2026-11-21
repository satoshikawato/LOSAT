//! OID masks (`EOidMaskType`): the one owner of what makeblastdb's `-oid_masks` writes and the
//! alias key `OID_MASK_TYPE` reads — the mask types, the file of a volume and its name, and the
//! test of a predicted RefSeq accession that both sides apply to a sequence's ids.
//!
//! The file `<volume>.pxm` / `.nxm` of `fExcludeModel` is an OID list file of the volume, whose
//! layout `seqdb::oid_list` owns, with the bit 0 for an OID that is masked (`writedb::oid_list`
//! writes it).

// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdbcommon.hpp:147-150

use crate::seqid::accguide::{identify_seq_id, parse_flags, SeqIdWarning, F_PREDICTED};
use crate::seqid::ParsedSeqId;

/// `EOidMaskType::fExcludeModel`: sequences of predicted RefSeq accessions (`XP_`, `XM_`, `XR_`,
/// ...), the one type of makeblastdb 2.17.0.
pub(crate) const EXCLUDE_MODEL: u64 = 0x01;

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbcommon.cpp:2707-2717
/// `SeqDB_GetOidMaskFileExt` of `fExcludeModel`: the extension of a volume's mask file.
pub(crate) fn exclude_model_extension(protein: bool) -> &'static str {
    if protein {
        "pxm"
    } else {
        "nxm"
    }
}

/// `IdentifyAccession() & CSeq_id::fAcc_predicted` of `id`, with the default flags of the member
/// function; `warn` gets the accession guide's warnings, which NCBI posts at every call.
pub(crate) fn is_predicted(id: &ParsedSeqId, warn: &mut dyn FnMut(SeqIdWarning)) -> bool {
    identify_seq_id(id, parse_flags::ANY_RAW | parse_flags::ANY_LOCAL, warn) & F_PREDICTED != 0
}
