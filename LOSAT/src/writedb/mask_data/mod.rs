//! The mask data that `makeblastdb -mask_data` stores in a database: the
//! registry of masking algorithms, the masked ranges of each sequence, and
//! the column files (`.?aa`, `.?ab`, `.?ac`) titled `BlastDb/MaskData`.
//!
//! Line NP, M4 preparation: nothing calls this yet (`LOSAT makedb` still
//! refuses `-mask_data`). The mask files themselves are read by
//! `writedb::mask_info`.
