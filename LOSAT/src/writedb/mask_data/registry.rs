//! The registry of the masking algorithms of a database being written
//! (`CMaskInfoRegistry`,
//! `c++/src/objtools/blast/seqdb_writer/mask_info_registry.cpp`): it gives
//! each algorithm the id under which its masks are stored.

use std::collections::BTreeSet;

use crate::writedb::WriteDbError;

// NCBI reference (598d8ae6): c++/src/objects/blastdb/blastdb.asn:51-59
// ```asn
// Blast-filter-program ::= INTEGER {
//     not-set             (0),
//     dust                (10),
//     seg                 (20),
//     windowmasker        (30),
//     repeat              (40),
//     other               (100),
//     max                 (255)
// }
// ```
/// `eBlast_filter_program_not_set`.
pub const FILTER_PROGRAM_NOT_SET: i32 = 0;
/// `eBlast_filter_program_dust`.
pub const FILTER_PROGRAM_DUST: i32 = 10;
/// `eBlast_filter_program_seg`.
pub const FILTER_PROGRAM_SEG: i32 = 20;
/// `eBlast_filter_program_windowmasker`.
pub const FILTER_PROGRAM_WINDOWMASKER: i32 = 30;
/// `eBlast_filter_program_repeat`.
pub const FILTER_PROGRAM_REPEAT: i32 = 40;
/// `eBlast_filter_program_other`.
pub const FILTER_PROGRAM_OTHER: i32 = 100;
/// `eBlast_filter_program_max`.
pub const FILTER_PROGRAM_MAX: i32 = 255;

/// `CMaskInfoRegistry`: the algorithms registered so far and the ids in use.
#[derive(Debug, Default)]
pub struct MaskInfoRegistry {
    /// `m_UsedIds`.
    used_ids: BTreeSet<i32>,
    /// `m_RegisteredAlgos`: the program, name and options of each algorithm
    /// registered by program, and the id string of each registered by id,
    /// all in one list (NCBI compares them as strings).
    registered: Vec<Vec<u8>>,
}

fn exception(message: String) -> WriteDbError {
    WriteDbError::Creation(message)
}

impl MaskInfoRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/mask_info_registry.hpp:66-68
    // ```c++
    //     bool IsRegistered(int algo_id) const {
    //         return m_UsedIds.find(algo_id) != m_UsedIds.end();
    //     }
    // ```
    /// `IsRegistered`: whether `algo_id` has been given to an algorithm.
    pub fn is_registered(&self, algo_id: i32) -> bool {
        self.used_ids.contains(&algo_id)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/mask_info_registry.cpp:62-74
    // ```c++
    // int
    // CMaskInfoRegistry::x_FindNextValidIdWithinRange(int start, int stop)
    // {
    //     for (int id = start; id < stop &&
    //                          id < (int)eBlast_filter_program_max; id++) {
    //         if (m_UsedIds.find(id) == m_UsedIds.end()) {
    //             return id;
    //         }
    //     }
    //     string msg("Too many IDs in range " + NStr::IntToString(start));
    //     msg += "-" + NStr::IntToString(stop);
    //     NCBI_THROW(CWriteDBException, eArgErr, msg);
    // }
    // ```
    fn find_next_valid_id_within_range(&self, start: i32, stop: i32) -> Result<i32, WriteDbError> {
        let mut id = start;
        while id < stop && id < FILTER_PROGRAM_MAX {
            if !self.used_ids.contains(&id) {
                return Ok(id);
            }
            id += 1;
        }
        Err(exception(format!("Too many IDs in range {start}-{stop}")))
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/mask_info_registry.cpp:41-60
    // ```c++
    // int
    // CMaskInfoRegistry::x_AssignId(int start, int end)
    // {
    //     return x_FindNextValidIdWithinRange(start, end);
    // }
    //
    // int
    // CMaskInfoRegistry::x_AssignId(int start, int end, bool use_start)
    // {
    //     if (use_start) {
    //         if (m_UsedIds.find(start) != m_UsedIds.end()) {
    //             string msg("Masking algorithm with default arguments "
    //                        "already provided");
    //             NCBI_THROW(CWriteDBException, eArgErr, msg);
    //         }
    //         return start;
    //     } else {
    //         return x_FindNextValidIdWithinRange(start+1, end);
    //     }
    // }
    // ```
    fn assign_id_with_default(
        &self,
        start: i32,
        end: i32,
        use_start: bool,
    ) -> Result<i32, WriteDbError> {
        if use_start {
            if self.used_ids.contains(&start) {
                return Err(exception(
                    "Masking algorithm with default arguments already provided".to_string(),
                ));
            }
            Ok(start)
        } else {
            self.find_next_valid_id_within_range(start + 1, end)
        }
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/mask_info_registry.cpp:76-89
    // ```c++
    // int
    // CMaskInfoRegistry::Add(const string& id)
    // {
    //     if (find(m_RegisteredAlgos.begin(), m_RegisteredAlgos.end(), id)
    //             != m_RegisteredAlgos.end()) {
    //         NCBI_THROW(CWriteDBException, eArgErr, "Duplicate masking algorithm found.");
    //     }
    //
    //     m_RegisteredAlgos.push_back(id);
    //
    //     int algo_id = x_AssignId(eBlast_filter_program_other, eBlast_filter_program_max);
    //     m_UsedIds.insert(algo_id);
    //     return algo_id;
    // }
    // ```
    /// `Add(const string & id)`: registers an algorithm named by the user
    /// (`makeblastdb -mask_id`); its id is the first free one from 100 on.
    pub fn add_id(&mut self, id: &[u8]) -> Result<i32, WriteDbError> {
        if self.registered.iter().any(|known| known == id) {
            return Err(exception("Duplicate masking algorithm found.".to_string()));
        }
        self.registered.push(id.to_vec());
        let algo_id =
            self.find_next_valid_id_within_range(FILTER_PROGRAM_OTHER, FILTER_PROGRAM_MAX)?;
        self.used_ids.insert(algo_id);
        Ok(algo_id)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/mask_info_registry.cpp:91-141
    // ```c++
    // int
    // CMaskInfoRegistry::Add(EBlast_filter_program program,
    //                        const string& options, /* = string() */
    //                        const string& progname /* = string() */)
    // {
    //
    //     int algo_id = (int)program;
    //     string algo_descr = NStr::IntToString(algo_id) + progname + options;
    //
    //     if (find(m_RegisteredAlgos.begin(), m_RegisteredAlgos.end(), algo_descr)
    //             != m_RegisteredAlgos.end()) {
    //         NCBI_THROW(CWriteDBException, eArgErr, "Duplicate masking algorithm found.");
    //     }
    //
    //     m_RegisteredAlgos.push_back(algo_descr);
    //
    //     switch (program) {
    //     case eBlast_filter_program_dust:
    //         algo_id = x_AssignId(program, eBlast_filter_program_seg,
    //                              options.empty());
    //         break;
    //
    //     case eBlast_filter_program_seg:
    //         algo_id = x_AssignId(program, eBlast_filter_program_windowmasker,
    //                              options.empty());
    //         break;
    //
    //     case eBlast_filter_program_windowmasker:
    //         algo_id = x_AssignId(program, eBlast_filter_program_repeat,
    //                              options.empty());
    //         break;
    //
    //     case eBlast_filter_program_repeat:
    //         algo_id = x_AssignId(program, eBlast_filter_program_other);
    //         break;
    //
    //     case eBlast_filter_program_other:
    //         algo_id = x_AssignId(program, eBlast_filter_program_max);
    //         break;
    //
    //     case eBlast_filter_program_not_set:
    //     case eBlast_filter_program_max:
    //     default:
    //         string msg("Invalid filtering program: ");
    //         msg += NStr::IntToString((int)program);
    //         NCBI_THROW(CWriteDBException, eArgErr, msg);
    //     }
    //
    //     m_UsedIds.insert(algo_id);
    //     return algo_id;
    // }
    // ```
    // `progname` is the name of the mask file for makeblastdb. The two-
    // argument `x_AssignId` (repeat, other) takes the first free id from the
    // program's own value on and ignores the options.
    /// `Add(program, options, progname)`: registers an algorithm of a
    /// filtering program. `program` is the `Blast-filter-program` number.
    pub fn add_program(
        &mut self,
        program: i32,
        options: &[u8],
        progname: &[u8],
    ) -> Result<i32, WriteDbError> {
        let mut algo_descr = program.to_string().into_bytes();
        algo_descr.extend_from_slice(progname);
        algo_descr.extend_from_slice(options);

        if self.registered.contains(&algo_descr) {
            return Err(exception("Duplicate masking algorithm found.".to_string()));
        }
        self.registered.push(algo_descr);

        let algo_id = match program {
            FILTER_PROGRAM_DUST => {
                self.assign_id_with_default(program, FILTER_PROGRAM_SEG, options.is_empty())?
            }
            FILTER_PROGRAM_SEG => self.assign_id_with_default(
                program,
                FILTER_PROGRAM_WINDOWMASKER,
                options.is_empty(),
            )?,
            FILTER_PROGRAM_WINDOWMASKER => {
                self.assign_id_with_default(program, FILTER_PROGRAM_REPEAT, options.is_empty())?
            }
            FILTER_PROGRAM_REPEAT => {
                self.find_next_valid_id_within_range(program, FILTER_PROGRAM_OTHER)?
            }
            FILTER_PROGRAM_OTHER => {
                self.find_next_valid_id_within_range(program, FILTER_PROGRAM_MAX)?
            }
            _ => {
                return Err(exception(format!("Invalid filtering program: {program}")));
            }
        };

        self.used_ids.insert(algo_id);
        Ok(algo_id)
    }
}
