//! What a run can open of a database's files, by its program, its options and the database's own
//! alias files (line WQA D26).
//!
//! No NCBI counterpart: `-logfile` must not name a file that the run reads
//! ([`SeqDb::read_paths`](super::SeqDb::read_paths)), and which of a database's files a run can
//! open depends on what it is asked to do and on what the database holds. A search opens the
//! sequence, header and index files, and the id indexes of each volume (the string index, the GI
//! index and the trace-id index), where it looks its hits and its id lists up: which of them it
//! opens depends on the ids that the database holds, which the check does not read, so all of them
//! are counted. The LMDB lookup files of a version 5 volume are counted with an option or an alias
//! file's list that asks for their kind, the masks with a mask option, the RPS files for the RPS
//! searches and `blastdbcheck -cdd_delta`. A tool that reads the sequences in OID order
//! (dustmasker, segmasker, windowmasker, makembindex) opens no more than the sequence, header and
//! index files and the id indexes that the database's alias files' lists ask for, and
//! makeblastdb's copy of a database adds its masks. A name that the run cannot open is not a file
//! it reads, so it is not refused (it is the run's log).

use super::id_list::IdListKind;

/// What a run can open of the files of a database (see the module).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct DbUse {
    /// The string index of a version 4 volume (`.nsi`/`.nsd`): a search's lookups, a string id
    /// of an id list.
    string_index: bool,
    /// The GI index (`.nni`/`.nnd`): a search's lookups of the hits of a database with GIs, a GI
    /// list.
    gi_index: bool,
    /// The trace-id index (`.nti`/`.ntd`): a search's lookups of a database with trace ids, a
    /// trace-id list.
    trace_index: bool,
    /// The IPG index (`.npi`/`.npd`): `-ipglist` or `-negative_ipglist`.
    ipg_index: bool,
    /// The LMDB file of the sequence ids of a version 5 volume that a negative seqid list reads
    /// (`.nos`): `-seqidlist` or `-negative_seqidlist`, an alias file's `SEQIDLIST`.
    seqid_lmdb: bool,
    /// The LMDB files of the taxonomy ids of a version 5 volume (`.not`, `.ntf`, `.nto`):
    /// `-taxids`, `-negative_taxids`, `-taxidlist`, `-negative_taxidlist`, an alias file's
    /// `TAXIDLIST`.
    taxid_lmdb: bool,
    /// The run is asked for a database mask (`-db_soft_mask`, `-db_hard_mask`, an imported
    /// search strategy, which may hold either, or makembindex's `-db_mask`), or copies the
    /// database's masks (makeblastdb's `-input_type blastdb`): the column files of the volumes,
    /// their GI index and the GI-mask files.
    masks: bool,
    /// The run is an RPS search, or `blastdbcheck -cdd_delta`: the four RPS files of each volume.
    rps: bool,
    /// `blastdbcheck -cdd_delta`: the two files of the CDD delta search of each volume
    /// (`.wcounts`, `.obsr`) besides the RPS files.
    cdd_delta: bool,
    /// The run is a tool that is not a search and looks the database up by its options
    /// (`blastdbcmd`, `blastdb_aliastool`, `blastdbcheck`, `blast_formatter`, ...): every file of
    /// the database that a tool of its kind may open, whatever its options say.
    tool: bool,
}

impl DbUse {
    /// The use of a tool that is not a search and looks the database up by its options: the
    /// database's files that any such tool opens.
    pub(crate) const TOOL: Self = Self {
        string_index: true,
        gi_index: true,
        trace_index: true,
        ipg_index: true,
        seqid_lmdb: true,
        taxid_lmdb: true,
        masks: true,
        rps: false,
        cdd_delta: false,
        tool: true,
    };

    /// The use of a tool that reads every sequence and its definition line in OID order
    /// (`seqmasks_io::bdb_reader::MaskBdbReader`: dustmasker, segmasker, windowmasker and
    /// makembindex with a database input) before its options and the database's alias files are
    /// read: the alias files and each volume's index, header and sequence files and LMDB file.
    pub(crate) const SEQUENCES: Self = Self {
        string_index: false,
        gi_index: false,
        trace_index: false,
        ipg_index: false,
        seqid_lmdb: false,
        taxid_lmdb: false,
        masks: false,
        rps: false,
        cdd_delta: false,
        tool: false,
    };

    /// The use of `program` before its options are read: a search opens its sequence, header and
    /// index files and the id indexes, and the RPS files for `rpsblast` and `rpstblastn`; a tool
    /// that reads the sequences in OID order [`DbUse::SEQUENCES`]; makeblastdb, which copies a
    /// database's sequences with their masks (`-input_type blastdb`), those and the masks; any
    /// other tool [`DbUse::TOOL`].
    pub(crate) fn of_program(program: &str) -> Self {
        let search = Self {
            string_index: true,
            gi_index: true,
            trace_index: true,
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
    /// a value or is a flag that was given.
    pub(crate) fn note_option(&mut self, name: &str) {
        match name {
            "db_soft_mask" | "db_hard_mask" | "db_mask" => self.masks = true,
            "gilist" | "negative_gilist" => self.note_list(IdListKind::Gi),
            "ipglist" | "negative_ipglist" => self.ipg_index = true,
            "seqidlist" | "negative_seqidlist" => self.note_list(IdListKind::Si),
            "taxids" | "negative_taxids" | "taxidlist" | "negative_taxidlist" => {
                self.note_list(IdListKind::TaxId);
            }
            // The strategy file may hold any of the options above.
            "import_search_strategy" => {
                self.masks = true;
                self.ipg_index = true;
                for kind in IdListKind::ALL {
                    self.note_list(kind);
                }
            }
            "cdd_delta" => {
                self.rps = true;
                self.cdd_delta = true;
            }
            _ => {}
        }
    }

    /// Takes in an id list of `kind` over the database (an option or an alias file's `GILIST`,
    /// `TILIST`, `SEQIDLIST` or `TAXIDLIST`): the files that look its ids up. A seqid list may
    /// hold ids of every kind: the string index, the GI and trace-id indexes and the LMDB file of
    /// the ids.
    pub(super) fn note_list(&mut self, kind: IdListKind) {
        match kind {
            IdListKind::Gi => self.gi_index = true,
            IdListKind::Ti => self.trace_index = true,
            IdListKind::Si => {
                self.string_index = true;
                self.gi_index = true;
                self.trace_index = true;
                self.seqid_lmdb = true;
            }
            IdListKind::TaxId => self.taxid_lmdb = true,
        }
    }

    /// Whether the run can open the GI-mask files of the database's `MASKLIST` and the column
    /// files of its volumes (`column::read_paths`).
    pub(super) fn reads_masks(self) -> bool {
        self.masks
    }

    /// Whether the run can open the file `.<mol><ext>` of a volume (`ext` of
    /// `VOLUME_EXTENSIONS`) or the file `ext` beside a version 5 volume's LMDB file (`ext` of
    /// `LMDB_EXTENSIONS`, `filter::lmdb_file_beside`). The column files `aa`, `ab` and `ac` are
    /// the columns' producer's (`column::read_paths`).
    pub(super) fn reads_extension(self, ext: &str) -> bool {
        match ext {
            // Opened by every run that reads a database: the volume's index, header and sequence
            // files, the alias file, and the LMDB file of a version 5 volume.
            "al" | "in" | "hr" | "sq" | "db" => true,
            "si" | "sd" => self.string_index,
            "ni" | "nd" => self.gi_index,
            "ti" | "td" => self.trace_index,
            "pi" | "pd" => self.ipg_index,
            "os" => self.seqid_lmdb,
            "ot" | "tf" | "to" => self.taxid_lmdb,
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

    /// Whether the run can open the two files of the CDD delta search of a volume
    /// (`rps::CddFreqsKind::Wcounts`, `Obsr`).
    pub(super) fn reads_cdd_delta(self) -> bool {
        self.cdd_delta
    }

    /// Whether the run can open the metadata file `.<mol>js` of an alias file.
    pub(super) fn reads_alias_metadata(self) -> bool {
        self.tool
    }
}
