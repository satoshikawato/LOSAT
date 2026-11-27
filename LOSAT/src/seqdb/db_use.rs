//! What a run can open of a database's files, by its program and options (line WQA D26 round 1).
//!
//! No NCBI counterpart: `-logfile` must not name a file that the run reads
//! ([`SeqDb::read_paths`](super::SeqDb::read_paths)), and which of a database's files a run can
//! open depends on what it is asked to do: its sequence, header and index files always, the
//! files of an id list or a taxonomy filter only with an option that asks for one, the masks only
//! with a mask option, the RPS files only for the RPS searches. A name that the run cannot open
//! is not a file it reads, so it is not refused (it is the run's log).

/// What a run can open of the files of a database (see the module).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct DbUse {
    /// The run is asked for a database mask (`-db_soft_mask`, `-db_hard_mask`, or an imported
    /// search strategy, which may hold either): the column files of the volumes and the GI-mask
    /// files.
    masks: bool,
    /// The run is asked for an id list or a taxonomy filter (`-gilist`, `-seqidlist`, `-taxids`
    /// and the others of [`DbUse::note_option`]): the ISAM and LMDB files that look an id up.
    id_filters: bool,
    /// The run is an RPS search: the four RPS files of each volume.
    rps: bool,
    /// The run is a tool that is not a search (`blastdbcmd`, `windowmasker`, ...): every file of
    /// the database that a tool of its kind may open, whatever its options say.
    tool: bool,
}

impl DbUse {
    /// The use of a tool that is not a search: the database's files that any tool opens.
    pub(crate) const TOOL: Self = Self {
        masks: true,
        id_filters: true,
        rps: false,
        tool: true,
    };

    /// The use of the search `program` before its options are read: its sequence, header and
    /// index files, and the RPS files of `rpsblast` and `rpstblastn`.
    pub(crate) fn search(program: &str) -> Self {
        Self {
            masks: false,
            id_filters: false,
            rps: matches!(program, "rpsblast" | "rpstblastn"),
            tool: false,
        }
    }

    /// Takes in the option `name` (its canonical name) of the run's arguments, when it was given
    /// a value.
    pub(crate) fn note_option(&mut self, name: &str) {
        match name {
            "db_soft_mask" | "db_hard_mask" => self.masks = true,
            "gilist" | "negative_gilist" | "seqidlist" | "negative_seqidlist" | "taxids"
            | "negative_taxids" | "taxidlist" | "negative_taxidlist" | "ipglist"
            | "negative_ipglist" => self.id_filters = true,
            // The strategy file may hold any of the options above.
            "import_search_strategy" => {
                self.masks = true;
                self.id_filters = true;
            }
            _ => {}
        }
    }

    /// Whether the run can open the GI-mask files of the database's `MASKLIST` and the column
    /// files of its volumes (`column::read_paths`).
    pub(super) fn reads_masks(self) -> bool {
        self.masks
    }

    /// Whether the run can open the file `.<mol><ext>` of a volume or an alias file (`ext` of
    /// `VOLUME_EXTENSIONS` or `LMDB_EXTENSIONS`). The column files `aa`, `ab` and `ac` are the
    /// columns' producer's (`column::read_paths`).
    pub(super) fn reads_extension(self, ext: &str) -> bool {
        match ext {
            // Opened by every search of a database: the volume's index, header and sequence
            // files, the string index of a version 4 database and the LMDB file of a version 5
            // one, and the alias file.
            "al" | "in" | "hr" | "sq" | "si" | "sd" | "db" => true,
            // Looked up with an id list or a taxonomy filter.
            "ni" | "nd" | "pi" | "pd" | "ti" | "td" | "os" | "ot" | "tf" | "to" => self.id_filters,
            // The GI index of a volume, opened with its GI masks.
            "og" => self.masks,
            "aa" | "ab" | "ac" => false,
            // The hash index, the metadata, the lock files: no search opens them.
            _ => self.tool,
        }
    }

    /// Whether the run can open the four RPS files of a volume ((`rps::RPS_EXTENSIONS`)).
    pub(super) fn reads_rps(self) -> bool {
        self.rps
    }

    /// Whether the run can open the metadata file `.<mol>js` of an alias file.
    pub(super) fn reads_alias_metadata(self) -> bool {
        self.tool
    }
}
