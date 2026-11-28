//! What a run can open of a database's files, by its program and options (line WQA D26 round 1).
//!
//! No NCBI counterpart: `-logfile` must not name a file that the run reads
//! ([`SeqDb::read_paths`](super::SeqDb::read_paths)), and which of a database's files a run can
//! open depends on what it is asked to do: its sequence, header and index files always, the
//! string index of a version 4 volume for a search, the files of an id list or a taxonomy filter
//! only with an option that asks for one, the masks only with a mask option, the RPS files only
//! for the RPS searches; a tool that reads the sequences in OID order (dustmasker, segmasker,
//! windowmasker, makembindex) opens no more than the sequence, header and index files, and
//! makeblastdb's copy of a database adds its masks. A name that the run cannot open is not a file
//! it reads, so it is not refused (it is the run's log).

/// What a run can open of the files of a database (see the module).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct DbUse {
    /// The run looks ids up in the volumes (a search, a tool that may look an id up): the string
    /// index of a version 4 volume. A tool that reads the sequences in OID order takes their ids
    /// from the definition lines and never opens it.
    string_index: bool,
    /// The run is asked for a database mask (`-db_soft_mask`, `-db_hard_mask`, an imported
    /// search strategy, which may hold either, or makembindex's `-db_mask`), or copies the
    /// database's masks (makeblastdb's `-input_type blastdb`): the column files of the volumes,
    /// their GI index and the GI-mask files.
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
    /// The run is a tool that is not a search and looks the database up by its options
    /// (`blastdbcmd`, `blastdb_aliastool`, `blastdbcheck`, ...): every file of the database that a
    /// tool of its kind may open, whatever its options say.
    tool: bool,
}

impl DbUse {
    /// The use of a tool that is not a search and looks the database up by its options: the
    /// database's files that any such tool opens.
    pub(crate) const TOOL: Self = Self {
        string_index: true,
        masks: true,
        gi_lists: true,
        ipg_lists: true,
        seqid_lists: true,
        taxids: true,
        rps: false,
        tool: true,
    };

    /// The use of a tool that reads every sequence and its definition line in OID order
    /// (`seqmasks_io::bdb_reader::MaskBdbReader`: dustmasker, segmasker, windowmasker and
    /// makembindex with a database input) before its options are read: the alias files and each
    /// volume's index, header and sequence files and LMDB file.
    pub(crate) const SEQUENCES: Self = Self {
        string_index: false,
        masks: false,
        gi_lists: false,
        ipg_lists: false,
        seqid_lists: false,
        taxids: false,
        rps: false,
        tool: false,
    };

    /// The use of `program` before its options are read: a search opens its sequence, header and
    /// index files, the string index, and the RPS files for `rpsblast` and `rpstblastn`; a tool
    /// that reads the sequences in OID order [`DbUse::SEQUENCES`]; makeblastdb, which copies a
    /// database's sequences with their masks (`-input_type blastdb`), those and the masks; any
    /// other tool [`DbUse::TOOL`].
    pub(crate) fn of_program(program: &str) -> Self {
        let search = Self {
            string_index: true,
            ..Self::SEQUENCES
        };
        match program {
            "blastn" | "blastp" | "blastx" | "tblastn" | "tblastx" | "psiblast" => search,
            "rpsblast" | "rpstblastn" => Self {
                rps: true,
                ..search
            },
            "dustmasker" | "segmasker" | "windowmasker" | "makembindex" => Self::SEQUENCES,
            "makeblastdb" => Self {
                masks: true,
                ..Self::SEQUENCES
            },
            _ => Self::TOOL,
        }
    }

    /// Takes in the option `name` (its canonical name) of the run's arguments, when it was given
    /// a value.
    pub(crate) fn note_option(&mut self, name: &str) {
        match name {
            "db_soft_mask" | "db_hard_mask" | "db_mask" => self.masks = true,
            "gilist" | "negative_gilist" => self.gi_lists = true,
            "ipglist" | "negative_ipglist" => self.ipg_lists = true,
            "seqidlist" | "negative_seqidlist" => self.seqid_lists = true,
            "taxids" | "negative_taxids" | "taxidlist" | "negative_taxidlist" => {
                self.taxids = true;
            }
            // The strategy file may hold any of the options above.
            "import_search_strategy" => {
                self.masks = true;
                self.gi_lists = true;
                self.ipg_lists = true;
                self.seqid_lists = true;
                self.taxids = true;
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
            // Opened by every run that reads a database: the volume's index, header and sequence
            // files, the alias file, and the LMDB file of a version 5 volume (the volume's index
            // file says which version it is, so it is counted for both).
            "al" | "in" | "hr" | "sq" | "db" => true,
            // The string index of a version 4 volume, which a run that looks ids up opens.
            "si" | "sd" => self.string_index,
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
}
