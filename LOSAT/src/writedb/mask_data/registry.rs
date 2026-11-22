//! The registry of the masking algorithms of a database being written
//! (`CMaskInfoRegistry`,
//! `c++/src/objtools/blast/seqdb_writer/mask_info_registry.cpp`): it gives
//! each algorithm the id under which its masks are stored.

use std::collections::BTreeSet;

use crate::seqdb::mask::FilterProgram;
use crate::writedb::WriteDbError;

/// The `Blast-filter-program` values the registry assigns ids from
/// (`crate::seqdb::mask::FILTER_PROGRAMS`).
pub(crate) const FILTER_PROGRAM_DUST: i32 = FilterProgram::Dust.value();
pub(crate) const FILTER_PROGRAM_SEG: i32 = FilterProgram::Seg.value();
pub(crate) const FILTER_PROGRAM_WINDOWMASKER: i32 = FilterProgram::Windowmasker.value();
pub(crate) const FILTER_PROGRAM_REPEAT: i32 = FilterProgram::Repeat.value();
pub(crate) const FILTER_PROGRAM_OTHER: i32 = FilterProgram::Other.value();
pub(crate) const FILTER_PROGRAM_MAX: i32 = FilterProgram::Max.value();

/// `CMaskInfoRegistry`: the algorithms registered so far and the ids in use.
#[derive(Debug, Default)]
pub(crate) struct MaskInfoRegistry {
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
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/mask_info_registry.hpp:66-68
    /// `IsRegistered`: whether `algo_id` has been given to an algorithm.
    pub(crate) fn is_registered(&self, algo_id: i32) -> bool {
        self.used_ids.contains(&algo_id)
    }

    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_writer/mask_info_registry.cpp:62-74
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
    /// `Add` of the `string` id: registers an algorithm named by the user
    /// (`makeblastdb -mask_id`); its id is the first free one from 100 on.
    pub(crate) fn add_id(&mut self, id: &[u8]) -> Result<i32, WriteDbError> {
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
    // `progname` is the name of the mask file for makeblastdb. The two-
    // argument `x_AssignId` (repeat, other) takes the first free id from the
    // program's own value on and ignores the options.
    /// `Add` of a program, its options and a program name: registers an algorithm of a
    /// filtering program. `program` is the `Blast-filter-program` number.
    pub(crate) fn add_program(
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
