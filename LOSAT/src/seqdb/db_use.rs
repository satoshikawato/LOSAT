//! What a run can open of a database's files, by its program and options (line WQA D26 round 1),
//! and of the taxonomy files found beside them on the search path (line WQA D35).
//!
//! No NCBI counterpart: `-logfile` must not name a file that the run reads
//! ([`SeqDb::read_paths`](super::SeqDb::read_paths), [`taxdb::read_paths`](super::taxdb::read_paths)),
//! and which of a database's files a run can open depends on what it is asked to do: its
//! sequence, header and index files always, the files of an id list or a taxonomy filter only
//! with an option that asks for one, the masks only with a mask option, the RPS files only for
//! the RPS searches. A name that the run cannot open is not a file it reads, so it is not
//! refused (it is the run's log).
//!
//! The taxonomy files are looked for along the search path, the working directory first: NCBI's
//! taxonomy name database (`taxdb.bti`, then the `taxdb.btd` beside the one found) by every run
//! that reads a database and by a search whose report looks every taxid up without one, and
//! `taxonomy4blast.sqlite3` by a run that widens a taxonomy filter of its database to the taxids
//! below (WQA DECISIONS 62: the files that the native open trace shows the run checking or
//! opening).

/// What a run can open of the files of a database (see the module).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct DbUse {
    /// The run is asked for a database mask (`-db_soft_mask`, `-db_hard_mask`, or an imported
    /// search strategy, which may hold either): the column files of the volumes, their GI index
    /// and the GI-mask files.
    masks: bool,
    /// `-gilist` or `-negative_gilist`: the numeric ISAM files that look a GI up.
    gi_lists: bool,
    /// `-ipglist` or `-negative_ipglist`: the numeric ISAM files that look an IPG up.
    ipg_lists: bool,
    /// `-seqidlist` or `-negative_seqidlist`: the LMDB file of the sequence ids and the trace-id
    /// ISAM files.
    seqid_lists: bool,
    /// `-taxids`, `-negative_taxids`, `-taxidlist` or `-negative_taxidlist`: the LMDB files of
    /// the taxonomy ids.
    taxids: bool,
    /// The run is an RPS search: the four RPS files of each volume.
    rps: bool,
    /// The run is a tool that is not a search (`blastdbcmd`, `windowmasker`, ...): every file of
    /// the database that a tool of its kind may open, whatever its options say.
    tool: bool,
    /// The run looks NCBI's taxonomy name database up even without a database: a report of the
    /// search looks every taxid of its subjects up, those of `-subject` too (the plural names and
    /// the organism report, `TabularSpec::looks_up_every_taxid`), or the run imports a search
    /// strategy, which may name a database.
    looks_up_taxonomy: bool,
    /// `-no_taxid_expansion`: a taxonomy filter keeps its own taxids, and the run does not look
    /// for `taxonomy4blast.sqlite3`.
    no_taxid_expansion: bool,
}

impl DbUse {
    /// The use of a tool that is not a search: the database's files that any tool opens.
    pub(crate) const TOOL: Self = Self {
        masks: true,
        gi_lists: true,
        ipg_lists: true,
        seqid_lists: true,
        taxids: true,
        rps: false,
        tool: true,
        looks_up_taxonomy: false,
        no_taxid_expansion: false,
    };

    /// The use of the search `program` before its options are read: its sequence, header and
    /// index files, and the RPS files of `rpsblast` and `rpstblastn`.
    pub(crate) fn search(program: &str) -> Self {
        Self {
            masks: false,
            gi_lists: false,
            ipg_lists: false,
            seqid_lists: false,
            taxids: false,
            rps: matches!(program, "rpsblast" | "rpstblastn"),
            tool: false,
            looks_up_taxonomy: false,
            no_taxid_expansion: false,
        }
    }

    /// Takes in the option `name` (its canonical name) of the run's arguments, when it was given
    /// a value.
    pub(crate) fn note_option(&mut self, name: &str) {
        match name {
            "db_soft_mask" | "db_hard_mask" => self.masks = true,
            "gilist" | "negative_gilist" => self.gi_lists = true,
            "ipglist" | "negative_ipglist" => self.ipg_lists = true,
            "seqidlist" | "negative_seqidlist" => self.seqid_lists = true,
            "taxids" | "negative_taxids" | "taxidlist" | "negative_taxidlist" => {
                self.taxids = true;
            }
            // The strategy file may hold any of the options above, and a database.
            "import_search_strategy" => {
                self.masks = true;
                self.gi_lists = true;
                self.ipg_lists = true;
                self.seqid_lists = true;
                self.taxids = true;
                self.looks_up_taxonomy = true;
            }
            _ => {}
        }
    }

    /// Takes in the flag `name` (its canonical name) of the run's arguments.
    pub(crate) fn note_flag(&mut self, name: &str) {
        if name == "no_taxid_expansion" {
            self.no_taxid_expansion = true;
        }
    }

    /// Takes in that a report of the search looks every taxid of its subjects up in NCBI's
    /// taxonomy name database (`TabularSpec::looks_up_every_taxid` of an `-outfmt`).
    pub(crate) fn note_every_taxid(&mut self) {
        self.looks_up_taxonomy = true;
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
            // files, the alias file, and the string index of a version 4 database or the LMDB
            // file of a version 5 one (which of the two the volume's index file says, so both
            // are counted).
            "al" | "in" | "hr" | "sq" | "si" | "sd" | "db" => true,
            // Looked up with an id list or a taxonomy filter, by kind.
            "ni" | "nd" => self.gi_lists,
            "pi" | "pd" => self.ipg_lists,
            "ti" | "td" | "os" => self.seqid_lists,
            "ot" | "tf" | "to" => self.taxids,
            // The GI index of a volume, opened with its GI masks.
            "og" => self.masks,
            "aa" | "ab" | "ac" => false,
            // The hash index, the metadata, the lock files: no search opens them.
            _ => self.tool,
        }
    }

    /// Whether the run can open the four RPS files of a volume (`rps::RPS_EXTENSIONS`).
    pub(super) fn reads_rps(self) -> bool {
        self.rps
    }

    /// Whether the run can open the metadata file `.<mol>js` of an alias file.
    pub(super) fn reads_alias_metadata(self) -> bool {
        self.tool
    }

    /// Whether the run looks for NCBI's taxonomy name database on its search path
    /// (`taxdb::read_paths`): every run that reads a database (`reads_database`: `-db`, BLASTN's
    /// `-filtering_db`, a tool's database; the database's reader looks for `taxdb.bti` whatever
    /// the report), and a search that looks it up without one (`looks_up_taxonomy`).
    pub(crate) fn reads_taxonomy_names(self, reads_database: bool) -> bool {
        reads_database || self.looks_up_taxonomy
    }

    /// Whether the run looks for `taxonomy4blast.sqlite3` on its search path to widen the
    /// taxonomy filter of the database it reads (`reads_database`) to the taxids below
    /// (`blastinput::db_filter`, `blastdbcmd -taxids`): a taxonomy filter without
    /// `-no_taxid_expansion`.
    pub(crate) fn reads_taxonomy_tree(self, reads_database: bool) -> bool {
        (reads_database || self.looks_up_taxonomy) && self.taxids && !self.no_taxid_expansion
    }
}
