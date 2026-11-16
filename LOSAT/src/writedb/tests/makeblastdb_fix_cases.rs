//! Runs of makeblastdb 2.17.0 (Linux x86-64) for the tests of the makedb-v5 PR2 review fixes: plain
//! FASTA ids with control bytes or bytes from 0x80 up (`CSeq_id_find_pred`), `-taxid_map -`, a
//! directory as `-taxid_map` and `-taxid_map ""`: the input written by the test, the standard
//! input, every word, the exit status, the two output streams and, for exit status 0, the files
//! of the output directory (bytes, or the size and FNV-1a digest of a file over 1024 bytes).
//!
//! Generated (`pr2/fix/scripts/gen_fix_cases.py` of the makedb-v5 line's task folder); do not
//! edit by hand. `<D>` is the output directory, `<T>` the build time and `<S>` the seconds of the
//! timed line; the input is `<D>/in.fa`.

use super::makeblastdb_v5_cases::Expected;

/// One run of makeblastdb on an input of its own.
pub(super) struct FixCase {
    pub name: &'static str,
    /// The FASTA text of `<D>/in.fa`.
    pub input: &'static [u8],
    /// The standard input.
    pub stdin: &'static [u8],
    /// Whether the directory `<D>/mapdir` is made.
    pub mapdir: bool,
    /// Every word after `makeblastdb`.
    pub words: &'static [&'static str],
    pub status: i32,
    pub stdout: &'static [u8],
    pub stderr: &'static [u8],
    /// The date stored in the database (empty when none was written).
    pub date: &'static str,
    /// The files of `<D>` but the inputs, sorted by name (exit status 0 only).
    pub files: &'static [(&'static str, Expected)],
}

pub(super) const FIX_CASES: &[FixCase] = &[
    FixCase {
        name: "L_bytes_prot",
        input: b"\
            >a\x02b t1\nMKTAYIAKQRQISFVKSHFSRQ\n>c\x1fd t2\nMKTAYIAKQRQISFVKSHFSRQ\n>ab\xc3\xa9 t3\n\
            MKTAYIAKQRQISFVKSHFSRQ\n>abc\xe9def title\nMKTAYIAKQRQISFVKSHFSRQ\n>MW\xc3\xa9_1 t5\nMKT\
            AYIAKQRQISFVKSHFSRQ\n>e\x7ff t6\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|x t\x01g\x03h u\nMKTAYIAKQ\
            RQISFVKSHFSRQ\n>lcl|q\x02r z\nMKTAYIAKQRQISFVKSHFSRQ\n>sp|P12345|ABC_HUMAN h\xe9llo\nMKT\
            AYIAKQRQISFVKSHFSRQ\n>ref|NP_000001.1| ok\x01i\x05j s\nMKTAYIAKQRQISFVKSHFSRQ\n>XP_01152\
            0000.1\x80 raw high\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|abcdefghij\x02! k1\nMKTAYIAKQRQISFVKSH\
            FSRQ\n>lcl|abcdefghij k2\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|abcdefghij\x02z k3\nMKTAYIAKQRQIS\
            FVKSHFSRQ\n>plain_ok t\nMKTAYIAKQRQISFVKSHFSRQ\n>k\xff w\nMKTAYIAKQRQISFVKSHFSRQ\n",
        stdin: b"",
        mapdir: false,
        words: &["-in", "{D}/in.fa", "-dbtype", "prot", "-out", "{D}/db", "-title", "t", "-parse_seqids", "-hash_index"],
        status: 0,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from \
            FASTA; added 16 sequences in <S> seconds.\n\n\n",
        stderr: b"\
            Error: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0x2] in string at byt\
            e 8\n\x02b t1\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0x1F] \
            in string at byte 8\n\x1fd t2\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\n\
            Bad char [0x2] in string at byte 8\n\x02b t1\nError: (803.7) [makeblastdb] Blast-def-lin\
            e-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xa9 t3\nError: (803.7) [makebla\
            stdb] Blast-def-line-set.E.title\nBad char [0xA9] in string at byte 9\n\xc3\xa9 t3\nErro\
            r: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0x1F] in string at byte 8\
            \n\x1fd t2\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xE9] in \
            string at byte 8\n\xe9def title\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\
            \nBad char [0xC3] in string at byte 8\n\xc3\xa9 t3\nError: (803.7) [makeblastdb] Blast-d\
            ef-line-set.E.title\nBad char [0xA9] in string at byte 9\n\xc3\xa9 t3\nError: (803.7) [m\
            akeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xa9_1 \
            t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xA9] in string a\
            t byte 9\n\xc3\xa9_1 t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad ch\
            ar [0xE9] in string at byte 8\n\xe9def title\nError: (803.7) [makeblastdb] Blast-def-lin\
            e-set.E.seqid.E.local.str\nBad char [0x7F] in string at byte 23\ne\x7ff\nError: (803.7) \
            [makeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xa9_\
            1 t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xA9] in string\
            \x20at byte 9\n\xc3\xa9_1 t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.\
            local.str\nBad char [0x7F] in string at byte 23\ne\x7ff\nError: (803.7) [makeblastdb] Bl\
            ast-def-line-set.E.seqid.E.local.str\nBad char [0x2] in string at byte 22\nq\x02r\nError\
            : (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xE9] in string at byte 9\
            \nh\xe9llo\nError: (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.local.str\nBad cha\
            r [0x2] in string at byte 22\nq\x02r\nError: (803.7) [makeblastdb] Blast-def-line-set.E.\
            title\nBad char [0xE9] in string at byte 9\nh\xe9llo\nError: (803.7) [makeblastdb] Blast\
            -def-line-set.E.title\nBad char [0x80] in string at byte 8\n\x80 raw high\nError: (803.7\
            ) [makeblastdb] Blast-def-line-set.E.seqid.E.local.str\nBad char [0x2] in string at byte\
            \x2032\nabcdefghij\x02!\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad ch\
            ar [0x80] in string at byte 8\n\x80 raw high\nError: (803.7) [makeblastdb] Blast-def-lin\
            e-set.E.seqid.E.local.str\nBad char [0x2] in string at byte 32\nabcdefghij\x02!\nError: \
            (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.local.str\nBad char [0x2] in string a\
            t byte 32\nabcdefghij\x02z\nError: (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.lo\
            cal.str\nBad char [0x2] in string at byte 32\nabcdefghij\x02z\nError: (803.7) [makeblast\
            db] Blast-def-line-set.E.title\nBad char [0xFF] in string at byte 8\n\xff w\nError: (803\
            .7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xFF] in string at byte 8\n\xff \
            w\n",
        date: "Oct 30, 2026  4:38 AM",
        files: &[
            ("db.pdb", Expected::Digest(32768, 0x396f35cc41df4dd6)),
            ("db.phd", Expected::Bytes(b"\
                2082731551\x020\n2082731551\x021\n2082731551\x0210\n2082731551\x0211\n2082731551\x0212\n\
                2082731551\x0213\n2082731551\x0214\n2082731551\x0215\n2082731551\x022\n2082731551\x023\n\
                2082731551\x024\n2082731551\x025\n2082731551\x026\n2082731551\x027\n2082731551\x028\n208\
                2731551\x029\n")),
            ("db.phi", Expected::Bytes(b"\
                \x00\x00\x00\x01\x00\x00\x00\x02\x00\x00\x00\xd6\x00\x00\x00\x10\x00\x00\x00\x01\x00\x00\
                \x00@\x00\x00\x10\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xd6\
                \x00\x00\x004\x00\x00\x00A2082731551\x020\x00")),
            ("db.phr", Expected::Bytes(b"\
                0\x800\x80\xa0\x80\x1a\x05#b t1\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01a\x00\x00\
                \x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\
                \x1a\x05#d t2\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01c\x00\x00\x00\x00\x00\x00\x00\
                \x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x05## t3\x00\x00\
                \xa1\x800\x80\xa0\x80\xa1\x80\x1a\x02ab\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\
                \x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\n#def title\x00\x00\xa1\x800\x80\xa0\
                \x80\xa1\x80\x1a\x03abc\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\
                \x00\x00\x000\x800\x80\xa0\x80\x1a\x07##_1 t5\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\
                \x02MW\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\
                \x80\xa0\x80\x1a\x02t6\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x03e#f\x00\x00\x00\x00\
                \x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01t\
                \x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01x\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\
                \x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01z\x00\x00\xa1\x800\x80\xa0\
                \x80\xa1\x80\x1a\x03q#r\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\
                \x00\x00\x000\x800\x80\xa0\x80\x1a\x05h#llo\x00\x00\xa1\x800\x80\xa7\x800\x80\xa0\x80\
                \x1a\tABC_HUMAN\x00\x00\xa1\x80\x1a\x06P12345\x00\x00\xa2\x80\x1a\x08reviewed\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\
                \xa0\x80\x1a\x02ok\x00\x00\xa1\x800\x80\xa9\x800\x80\xa1\x80\x1a\tNP_000001\x00\x00\xa3\
                \x80\x02\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\
                \x00\x00\x000\x800\x80\xa0\x80\x1a\n# raw high\x00\x00\xa1\x800\x80\xa9\x800\x80\xa1\x80\
                \x1a\x0cXP_011520000\x00\x00\xa3\x80\x02\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\
                \xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x02k1\x00\x00\xa1\x80\
                0\x80\xa0\x80\xa1\x80\x1a\x0cabcdefghij#!\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\
                \x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x02k2\x00\x00\xa1\x800\x80\xa0\
                \x80\xa1\x80\x1a\nabcdefghij\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\
                \x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x02k3\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\
                \x0cabcdefghij#z\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\
                \x000\x800\x80\xa0\x80\x1a\x01t\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x08plain_ok\x00\
                \x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\
                \x80\x1a\x03# w\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01k\x00\x00\x00\x00\x00\x00\
                \x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x00")),
            ("db.pin", Expected::Bytes(b"\
                \x00\x00\x00\x05\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x01t\x00\x00\x00\x06db.pdb\
                \x00\x00\x00\x19Oct 30, 2026  4:38 AM\x00\x00\x00\x00\x00\x00\x00\x10`\x01\x00\x00\x00\
                \x00\x00\x00\x00\x00\x00\x16\x00\x00\x00\x00\x00\x00\x00-\x00\x00\x00Z\x00\x00\x00\x88\
                \x00\x00\x00\xbc\x00\x00\x00\xec\x00\x00\x01\x18\x00\x00\x01A\x00\x00\x01l\x00\x00\x01\
                \xbf\x00\x00\x01\xfc\x00\x00\x02D\x00\x00\x02y\x00\x00\x02\xac\x00\x00\x02\xe1\x00\x00\
                \x03\x11\x00\x00\x03<\x00\x00\x00\x01\x00\x00\x00\x18\x00\x00\x00/\x00\x00\x00F\x00\x00\
                \x00]\x00\x00\x00t\x00\x00\x00\x8b\x00\x00\x00\xa2\x00\x00\x00\xb9\x00\x00\x00\xd0\x00\
                \x00\x00\xe7\x00\x00\x00\xfe\x00\x00\x01\x15\x00\x00\x01,\x00\x00\x01C\x00\x00\x01Z\x00\
                \x00\x01q")),
            ("db.pjs", Expected::Bytes(b"\
                {\n  \"version\": \"1.2\",\n  \"dbname\": \"db\",\n  \"dbtype\": \"Protein\",\n  \"db-ve\
                rsion\": 5,\n  \"description\": \"t\",\n  \"number-of-letters\": 352,\n  \"number-of-seq\
                uences\": 16,\n  \"last-updated\": \"2026-10-30T04:38:00\",\n  \"number-of-volumes\": 1,\
                \n  \"bytes-total\": 51442,\n  \"bytes-to-cache\": 577,\n  \"files\": [\n    \"db.pdb\",\
                \n    \"db.phd\",\n    \"db.phi\",\n    \"db.phr\",\n    \"db.pin\",\n    \"db.pog\",\n \
                \x20  \"db.pos\",\n    \"db.pot\",\n    \"db.psq\",\n    \"db.ptf\",\n    \"db.pto\"\n  \
                ]\n}\n")),
            ("db.pog", Expected::Bytes(b"\
                \x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x10\x00\x00\x00\x00\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\
                \xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\
                \xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\
                \xff\xff\xff\xff\xff\xff\xff\xff")),
            ("db.pos", Expected::Bytes(b"\
                \x10\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x00\x00\
                \x00\x00\x07\x00\x00\x00\x00\x00\x00\x00\x0b\x00\x00\x00\x00\x00\x00\x00\x0e\x00\x00\x00\
                \x00\x00\x00\x00\x12\x00\x00\x00\x00\x00\x00\x00\x14\x00\x00\x00\x00\x00\x00\x00\x18\x00\
                \x00\x00\x00\x00\x00\x00\x1f\x00\x00\x00\x00\x00\x00\x00+\x00\x00\x00\x00\x00\x00\x00:\
                \x00\x00\x00\x00\x00\x00\x00G\x00\x00\x00\x00\x00\x00\x00R\x00\x00\x00\x00\x00\x00\x00_\
                \x00\x00\x00\x00\x00\x00\x00h\x00\x00\x00\x00\x00\x00\x00j\x00\x00\x00\x00\x00\x00\x00\
                \x01a\x01c\x02ab\x03abc\x02MW\x03e\x7ff\x01x\x03q\x02r\x06P12345\x0bNP_000001.1\x0eXP_01\
                1520000.1\x0cabcdefghij\x02!\nabcdefghij\x0cabcdefghij\x02z\x08plain_ok\x01k")),
            ("db.pot", Expected::Bytes(b"\
                \x10\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x00\x00\x00\
                \x00\x00\x03\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x00\x00\x00\x00\x05\x00\x00\x00\
                \x00\x00\x00\x00\x06\x00\x00\x00\x00\x00\x00\x00\x07\x00\x00\x00\x00\x00\x00\x00\x08\x00\
                \x00\x00\x00\x00\x00\x00\t\x00\x00\x00\x00\x00\x00\x00\n\x00\x00\x00\x00\x00\x00\x00\x0b\
                \x00\x00\x00\x00\x00\x00\x00\x0c\x00\x00\x00\x00\x00\x00\x00\x0d\x00\x00\x00\x00\x00\x00\
                \x00\x0e\x00\x00\x00\x00\x00\x00\x00\x0f\x00\x00\x00\x00\x00\x00\x00\x10\x00\x00\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\
                \x00")),
            ("db.psq", Expected::Bytes(b"\
                \x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\
                \x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\
                \x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\
                \x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\
                \x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\
                \x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\
                \x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\
                \x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\
                \x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\
                \x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\
                \x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\
                \x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\
                \x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\
                \x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\
                \x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\
                \x11\x10\x0f\x00")),
            ("db.ptf", Expected::Digest(16384, 0xc29f655d5a3815f8)),
            ("db.pto", Expected::Bytes(b"\
                \x10\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x02\x00\x00\x00\x03\x00\x00\x00\x04\x00\
                \x00\x00\x05\x00\x00\x00\x06\x00\x00\x00\x07\x00\x00\x00\x08\x00\x00\x00\t\x00\x00\x00\n\
                \x00\x00\x00\x0b\x00\x00\x00\x0c\x00\x00\x00\x0d\x00\x00\x00\x0e\x00\x00\x00\x0f\x00\x00\
                \x00")),
        ],
    },
    FixCase {
        name: "L_bytes_prot_v4",
        input: b"\
            >a\x02b t1\nMKTAYIAKQRQISFVKSHFSRQ\n>c\x1fd t2\nMKTAYIAKQRQISFVKSHFSRQ\n>ab\xc3\xa9 t3\n\
            MKTAYIAKQRQISFVKSHFSRQ\n>abc\xe9def title\nMKTAYIAKQRQISFVKSHFSRQ\n>MW\xc3\xa9_1 t5\nMKT\
            AYIAKQRQISFVKSHFSRQ\n>e\x7ff t6\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|x t\x01g\x03h u\nMKTAYIAKQ\
            RQISFVKSHFSRQ\n>lcl|q\x02r z\nMKTAYIAKQRQISFVKSHFSRQ\n>sp|P12345|ABC_HUMAN h\xe9llo\nMKT\
            AYIAKQRQISFVKSHFSRQ\n>ref|NP_000001.1| ok\x01i\x05j s\nMKTAYIAKQRQISFVKSHFSRQ\n>XP_01152\
            0000.1\x80 raw high\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|abcdefghij\x02! k1\nMKTAYIAKQRQISFVKSH\
            FSRQ\n>lcl|abcdefghij k2\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|abcdefghij\x02z k3\nMKTAYIAKQRQIS\
            FVKSHFSRQ\n>plain_ok t\nMKTAYIAKQRQISFVKSHFSRQ\n>k\xff w\nMKTAYIAKQRQISFVKSHFSRQ\n",
        stdin: b"",
        mapdir: false,
        words: &["-in", "{D}/in.fa", "-dbtype", "prot", "-out", "{D}/db", "-title", "t", "-parse_seqids", "-hash_index", "-blastdb_version", "4"],
        status: 0,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from \
            FASTA; added 16 sequences in <S> seconds.\n\n\n",
        stderr: b"\
            Error: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0x2] in string at byt\
            e 8\n\x02b t1\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0x1F] \
            in string at byte 8\n\x1fd t2\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\n\
            Bad char [0x2] in string at byte 8\n\x02b t1\nError: (803.7) [makeblastdb] Blast-def-lin\
            e-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xa9 t3\nError: (803.7) [makebla\
            stdb] Blast-def-line-set.E.title\nBad char [0xA9] in string at byte 9\n\xc3\xa9 t3\nErro\
            r: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0x1F] in string at byte 8\
            \n\x1fd t2\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xE9] in \
            string at byte 8\n\xe9def title\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\
            \nBad char [0xC3] in string at byte 8\n\xc3\xa9 t3\nError: (803.7) [makeblastdb] Blast-d\
            ef-line-set.E.title\nBad char [0xA9] in string at byte 9\n\xc3\xa9 t3\nError: (803.7) [m\
            akeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xa9_1 \
            t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xA9] in string a\
            t byte 9\n\xc3\xa9_1 t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad ch\
            ar [0xE9] in string at byte 8\n\xe9def title\nError: (803.7) [makeblastdb] Blast-def-lin\
            e-set.E.seqid.E.local.str\nBad char [0x7F] in string at byte 23\ne\x7ff\nError: (803.7) \
            [makeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xa9_\
            1 t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xA9] in string\
            \x20at byte 9\n\xc3\xa9_1 t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.\
            local.str\nBad char [0x7F] in string at byte 23\ne\x7ff\nError: (803.7) [makeblastdb] Bl\
            ast-def-line-set.E.seqid.E.local.str\nBad char [0x2] in string at byte 22\nq\x02r\nError\
            : (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xE9] in string at byte 9\
            \nh\xe9llo\nError: (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.local.str\nBad cha\
            r [0x2] in string at byte 22\nq\x02r\nError: (803.7) [makeblastdb] Blast-def-line-set.E.\
            title\nBad char [0xE9] in string at byte 9\nh\xe9llo\nError: (803.7) [makeblastdb] Blast\
            -def-line-set.E.title\nBad char [0x80] in string at byte 8\n\x80 raw high\nError: (803.7\
            ) [makeblastdb] Blast-def-line-set.E.seqid.E.local.str\nBad char [0x2] in string at byte\
            \x2032\nabcdefghij\x02!\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad ch\
            ar [0x80] in string at byte 8\n\x80 raw high\nError: (803.7) [makeblastdb] Blast-def-lin\
            e-set.E.seqid.E.local.str\nBad char [0x2] in string at byte 32\nabcdefghij\x02!\nError: \
            (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.local.str\nBad char [0x2] in string a\
            t byte 32\nabcdefghij\x02z\nError: (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.lo\
            cal.str\nBad char [0x2] in string at byte 32\nabcdefghij\x02z\nError: (803.7) [makeblast\
            db] Blast-def-line-set.E.title\nBad char [0xFF] in string at byte 8\n\xff w\nError: (803\
            .7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xFF] in string at byte 8\n\xff \
            w\n",
        date: "Oct 30, 2026  4:38 AM",
        files: &[
            ("db.phd", Expected::Bytes(b"\
                2082731551\x020\n2082731551\x021\n2082731551\x0210\n2082731551\x0211\n2082731551\x0212\n\
                2082731551\x0213\n2082731551\x0214\n2082731551\x0215\n2082731551\x022\n2082731551\x023\n\
                2082731551\x024\n2082731551\x025\n2082731551\x026\n2082731551\x027\n2082731551\x028\n208\
                2731551\x029\n")),
            ("db.phi", Expected::Bytes(b"\
                \x00\x00\x00\x01\x00\x00\x00\x02\x00\x00\x00\xd6\x00\x00\x00\x10\x00\x00\x00\x01\x00\x00\
                \x00@\x00\x00\x10\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xd6\
                \x00\x00\x004\x00\x00\x00A2082731551\x020\x00")),
            ("db.phr", Expected::Bytes(b"\
                0\x800\x80\xa0\x80\x1a\x05#b t1\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01a\x00\x00\
                \x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\
                \x1a\x05#d t2\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01c\x00\x00\x00\x00\x00\x00\x00\
                \x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x05## t3\x00\x00\
                \xa1\x800\x80\xa0\x80\xa1\x80\x1a\x02ab\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\
                \x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\n#def title\x00\x00\xa1\x800\x80\xa0\
                \x80\xa1\x80\x1a\x03abc\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\
                \x00\x00\x000\x800\x80\xa0\x80\x1a\x07##_1 t5\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\
                \x02MW\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\
                \x80\xa0\x80\x1a\x02t6\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x03e#f\x00\x00\x00\x00\
                \x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01t\
                \x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01x\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\
                \x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01z\x00\x00\xa1\x800\x80\xa0\
                \x80\xa1\x80\x1a\x03q#r\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\
                \x00\x00\x000\x800\x80\xa0\x80\x1a\x05h#llo\x00\x00\xa1\x800\x80\xa7\x800\x80\xa0\x80\
                \x1a\tABC_HUMAN\x00\x00\xa1\x80\x1a\x06P12345\x00\x00\xa2\x80\x1a\x08reviewed\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\
                \xa0\x80\x1a\x02ok\x00\x00\xa1\x800\x80\xa9\x800\x80\xa1\x80\x1a\tNP_000001\x00\x00\xa3\
                \x80\x02\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\
                \x00\x00\x000\x800\x80\xa0\x80\x1a\n# raw high\x00\x00\xa1\x800\x80\xa9\x800\x80\xa1\x80\
                \x1a\x0cXP_011520000\x00\x00\xa3\x80\x02\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\
                \xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x02k1\x00\x00\xa1\x80\
                0\x80\xa0\x80\xa1\x80\x1a\x0cabcdefghij#!\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\
                \x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x02k2\x00\x00\xa1\x800\x80\xa0\
                \x80\xa1\x80\x1a\nabcdefghij\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\
                \x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x02k3\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\
                \x0cabcdefghij#z\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\
                \x000\x800\x80\xa0\x80\x1a\x01t\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x08plain_ok\x00\
                \x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\
                \x80\x1a\x03# w\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01k\x00\x00\x00\x00\x00\x00\
                \x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x00")),
            ("db.pin", Expected::Bytes(b"\
                \x00\x00\x00\x04\x00\x00\x00\x01\x00\x00\x00\x01t\x00\x00\x00\x17Oct 30, 2026  4:38 AM\
                \x00\x00\x00\x00\x00\x10`\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x16\x00\x00\x00\x00\
                \x00\x00\x00-\x00\x00\x00Z\x00\x00\x00\x88\x00\x00\x00\xbc\x00\x00\x00\xec\x00\x00\x01\
                \x18\x00\x00\x01A\x00\x00\x01l\x00\x00\x01\xbf\x00\x00\x01\xfc\x00\x00\x02D\x00\x00\x02y\
                \x00\x00\x02\xac\x00\x00\x02\xe1\x00\x00\x03\x11\x00\x00\x03<\x00\x00\x00\x01\x00\x00\
                \x00\x18\x00\x00\x00/\x00\x00\x00F\x00\x00\x00]\x00\x00\x00t\x00\x00\x00\x8b\x00\x00\x00\
                \xa2\x00\x00\x00\xb9\x00\x00\x00\xd0\x00\x00\x00\xe7\x00\x00\x00\xfe\x00\x00\x01\x15\x00\
                \x00\x01,\x00\x00\x01C\x00\x00\x01Z\x00\x00\x01q")),
            ("db.pjs", Expected::Bytes(b"\
                {\n  \"version\": \"1.2\",\n  \"dbname\": \"db\",\n  \"dbtype\": \"Protein\",\n  \"db-ve\
                rsion\": 4,\n  \"description\": \"t\",\n  \"number-of-letters\": 352,\n  \"number-of-seq\
                uences\": 16,\n  \"last-updated\": \"2026-10-30T04:38:00\",\n  \"number-of-volumes\": 1,\
                \n  \"bytes-total\": 2159,\n  \"bytes-to-cache\": 561,\n  \"files\": [\n    \"db.phd\",\
                \n    \"db.phi\",\n    \"db.phr\",\n    \"db.pin\",\n    \"db.pog\",\n    \"db.psd\",\n \
                \x20  \"db.psi\",\n    \"db.psq\"\n  ]\n}\n")),
            ("db.pog", Expected::Bytes(b"\
                \x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x10\x00\x00\x00\x00\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\
                \xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\
                \xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\
                \xff\xff\xff\xff\xff\xff\xff\xff")),
            ("db.psd", Expected::Bytes(b"\
                a\x020\nab\x022\nabc\x023\nabc_human\x028\nabcdefghij\x02!\x0211\nabcdefghij\x0212\nabcd\
                efghij\x02z\x0213\nc\x021\ne\x7ff\x025\nk\x0215\nlcl|a\x020\nlcl|ab\x022\nlcl|abc\x023\n\
                lcl|abcdefghij\x02!\x0211\nlcl|abcdefghij\x0212\nlcl|abcdefghij\x02z\x0213\nlcl|c\x021\n\
                lcl|e\x7ff\x025\nlcl|k\x0215\nlcl|mw\x024\nlcl|plain_ok\x0214\nlcl|q\x02r\x027\nlcl|x\
                \x026\nmw\x024\nnp_000001\x029\nnp_000001.1\x029\np12345\x028\nplain_ok\x0214\nq\x02r\
                \x027\nx\x026\nxp_011520000\x0210\nxp_011520000.1\x0210\n")),
            ("db.psi", Expected::Bytes(b"\x00\x00\x00\x01\x00\x00\x00\x02\x00\x00\x01S\x00\x00\x00 \x00\x00\x00\x01\x00\x00\x00@\x00\x00\x10\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x01S\x00\x00\x004\x00\x00\x008a\x020\x00")),
            ("db.psq", Expected::Bytes(b"\
                \x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\
                \x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\
                \x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\
                \x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\
                \x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\
                \x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\
                \x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\
                \x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\
                \x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\
                \x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\
                \x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\
                \x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\
                \x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\
                \x08\x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\
                \x06\x11\x10\x0f\x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\
                \x11\x10\x0f\x00")),
        ],
    },
    FixCase {
        name: "L_bytes_nucl",
        input: b"\
            >a\x02b t1\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>c\x1fd t2\nACGTTGCAACGTAGCTAGCTA\
            GGATCCATGCATGCAAGCTT\n>ab\xc3\xa9 t3\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>abc\
            \xe9def title\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>MW\xc3\xa9_1 t5\nACGTTGCAACGT\
            AGCTAGCTAGGATCCATGCATGCAAGCTT\n>e\x7ff t6\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>l\
            cl|x t\x01g\x03h u\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>lcl|q\x02r z\nACGTTGCAAC\
            GTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>sp|P12345|ABC_HUMAN h\xe9llo\nACGTTGCAACGTAGCTAGCTAGGA\
            TCCATGCATGCAAGCTT\n>ref|NP_000001.1| ok\x01i\x05j s\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGC\
            AAGCTT\n>XP_011520000.1\x80 raw high\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>lcl|ab\
            cdefghij\x02! k1\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>lcl|abcdefghij k2\nACGTTGC\
            AACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>lcl|abcdefghij\x02z k3\nACGTTGCAACGTAGCTAGCTAGGATCC\
            ATGCATGCAAGCTT\n>plain_ok t\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>k\xff w\nACGTTG\
            CAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n",
        stdin: b"",
        mapdir: false,
        words: &["-in", "{D}/in.fa", "-dbtype", "nucl", "-out", "{D}/db", "-title", "t", "-parse_seqids"],
        status: 0,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences fr\
            om FASTA; added 16 sequences in <S> seconds.\n\n\n",
        stderr: b"\
            Error: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0x2] in string at byt\
            e 8\n\x02b t1\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0x1F] \
            in string at byte 8\n\x1fd t2\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\n\
            Bad char [0x2] in string at byte 8\n\x02b t1\nError: (803.7) [makeblastdb] Blast-def-lin\
            e-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xa9 t3\nError: (803.7) [makebla\
            stdb] Blast-def-line-set.E.title\nBad char [0xA9] in string at byte 9\n\xc3\xa9 t3\nErro\
            r: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0x1F] in string at byte 8\
            \n\x1fd t2\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xE9] in \
            string at byte 8\n\xe9def title\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\
            \nBad char [0xC3] in string at byte 8\n\xc3\xa9 t3\nError: (803.7) [makeblastdb] Blast-d\
            ef-line-set.E.title\nBad char [0xA9] in string at byte 9\n\xc3\xa9 t3\nError: (803.7) [m\
            akeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xa9_1 \
            t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xA9] in string a\
            t byte 9\n\xc3\xa9_1 t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad ch\
            ar [0xE9] in string at byte 8\n\xe9def title\nError: (803.7) [makeblastdb] Blast-def-lin\
            e-set.E.seqid.E.local.str\nBad char [0x7F] in string at byte 23\ne\x7ff\nError: (803.7) \
            [makeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xa9_\
            1 t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xA9] in string\
            \x20at byte 9\n\xc3\xa9_1 t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.\
            local.str\nBad char [0x7F] in string at byte 23\ne\x7ff\nError: (803.7) [makeblastdb] Bl\
            ast-def-line-set.E.seqid.E.local.str\nBad char [0x2] in string at byte 22\nq\x02r\nError\
            : (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xE9] in string at byte 9\
            \nh\xe9llo\nError: (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.local.str\nBad cha\
            r [0x2] in string at byte 22\nq\x02r\nError: (803.7) [makeblastdb] Blast-def-line-set.E.\
            title\nBad char [0xE9] in string at byte 9\nh\xe9llo\nError: (803.7) [makeblastdb] Blast\
            -def-line-set.E.title\nBad char [0x80] in string at byte 8\n\x80 raw high\nError: (803.7\
            ) [makeblastdb] Blast-def-line-set.E.seqid.E.local.str\nBad char [0x2] in string at byte\
            \x2032\nabcdefghij\x02!\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad ch\
            ar [0x80] in string at byte 8\n\x80 raw high\nError: (803.7) [makeblastdb] Blast-def-lin\
            e-set.E.seqid.E.local.str\nBad char [0x2] in string at byte 32\nabcdefghij\x02!\nError: \
            (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.local.str\nBad char [0x2] in string a\
            t byte 32\nabcdefghij\x02z\nError: (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.lo\
            cal.str\nBad char [0x2] in string at byte 32\nabcdefghij\x02z\nError: (803.7) [makeblast\
            db] Blast-def-line-set.E.title\nBad char [0xFF] in string at byte 8\n\xff w\nError: (803\
            .7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xFF] in string at byte 8\n\xff \
            w\n",
        date: "Oct 30, 2026  4:38 AM",
        files: &[
            ("db.ndb", Expected::Digest(32768, 0x6d893ee487275ae0)),
            ("db.nhr", Expected::Bytes(b"\
                0\x800\x80\xa0\x80\x1a\x05#b t1\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01a\x00\x00\
                \x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\
                \x1a\x05#d t2\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01c\x00\x00\x00\x00\x00\x00\x00\
                \x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x05## t3\x00\x00\
                \xa1\x800\x80\xa0\x80\xa1\x80\x1a\x02ab\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\
                \x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\n#def title\x00\x00\xa1\x800\x80\xa0\
                \x80\xa1\x80\x1a\x03abc\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\
                \x00\x00\x000\x800\x80\xa0\x80\x1a\x07##_1 t5\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\
                \x02MW\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\
                \x80\xa0\x80\x1a\x02t6\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x03e#f\x00\x00\x00\x00\
                \x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01t\
                \x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01x\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\
                \x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01z\x00\x00\xa1\x800\x80\xa0\
                \x80\xa1\x80\x1a\x03q#r\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\
                \x00\x00\x000\x800\x80\xa0\x80\x1a\x05h#llo\x00\x00\xa1\x800\x80\xa7\x800\x80\xa0\x80\
                \x1a\tABC_HUMAN\x00\x00\xa1\x80\x1a\x06P12345\x00\x00\xa2\x80\x1a\x08reviewed\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\
                \xa0\x80\x1a\x02ok\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x0bNP_000001.1\x00\x00\x00\
                \x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\n\
                # raw high\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x0eXP_011520000.1\x00\x00\x00\x00\
                \x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x02k1\
                \x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x0cabcdefghij#!\x00\x00\x00\x00\x00\x00\x00\
                \x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x02k2\x00\x00\xa1\
                \x800\x80\xa0\x80\xa1\x80\x1a\nabcdefghij\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\
                \x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x02k3\x00\x00\xa1\x800\x80\xa0\
                \x80\xa1\x80\x1a\x0cabcdefghij#z\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\
                \x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01t\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\
                \x08plain_ok\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x00\
                0\x800\x80\xa0\x80\x1a\x03# w\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01k\x00\x00\x00\
                \x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x00")),
            ("db.nin", Expected::Bytes(b"\
                \x00\x00\x00\x05\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x01t\x00\x00\x00\x06db.ndb\
                \x00\x00\x00\x19Oct 30, 2026  4:38 AM\x00\x00\x00\x00\x00\x00\x00\x10\x90\x02\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00)\x00\x00\x00\x00\x00\x00\x00-\x00\x00\x00Z\x00\x00\x00\x88\
                \x00\x00\x00\xbc\x00\x00\x00\xec\x00\x00\x01\x18\x00\x00\x01A\x00\x00\x01l\x00\x00\x01\
                \xbf\x00\x00\x01\xf3\x00\x00\x022\x00\x00\x02g\x00\x00\x02\x9a\x00\x00\x02\xcf\x00\x00\
                \x02\xff\x00\x00\x03*\x00\x00\x00\x01\x00\x00\x00\x0c\x00\x00\x00\x17\x00\x00\x00\"\x00\
                \x00\x00-\x00\x00\x008\x00\x00\x00C\x00\x00\x00N\x00\x00\x00Y\x00\x00\x00d\x00\x00\x00o\
                \x00\x00\x00z\x00\x00\x00\x85\x00\x00\x00\x90\x00\x00\x00\x9b\x00\x00\x00\xa6\x00\x00\
                \x00\xb1\x00\x00\x00\x0c\x00\x00\x00\x17\x00\x00\x00\"\x00\x00\x00-\x00\x00\x008\x00\x00\
                \x00C\x00\x00\x00N\x00\x00\x00Y\x00\x00\x00d\x00\x00\x00o\x00\x00\x00z\x00\x00\x00\x85\
                \x00\x00\x00\x90\x00\x00\x00\x9b\x00\x00\x00\xa6\x00\x00\x00\xb1\x00\x00\x00\xb1")),
            ("db.njs", Expected::Bytes(b"\
                {\n  \"version\": \"1.2\",\n  \"dbname\": \"db\",\n  \"dbtype\": \"Nucleotide\",\n  \"db\
                -version\": 5,\n  \"description\": \"t\",\n  \"number-of-letters\": 656,\n  \"number-of-\
                sequences\": 16,\n  \"last-updated\": \"2026-10-30T04:38:00\",\n  \"number-of-volumes\":\
                \x201,\n  \"bytes-total\": 51021,\n  \"bytes-to-cache\": 453,\n  \"files\": [\n    \"db.\
                ndb\",\n    \"db.nhr\",\n    \"db.nin\",\n    \"db.nog\",\n    \"db.nos\",\n    \"db.not\
                \",\n    \"db.nsq\",\n    \"db.ntf\",\n    \"db.nto\"\n  ]\n}\n")),
            ("db.nog", Expected::Bytes(b"\
                \x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x10\x00\x00\x00\x00\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\
                \xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\
                \xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\
                \xff\xff\xff\xff\xff\xff\xff\xff")),
            ("db.nos", Expected::Bytes(b"\
                \x10\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x00\x00\
                \x00\x00\x07\x00\x00\x00\x00\x00\x00\x00\x0b\x00\x00\x00\x00\x00\x00\x00\x0e\x00\x00\x00\
                \x00\x00\x00\x00\x12\x00\x00\x00\x00\x00\x00\x00\x14\x00\x00\x00\x00\x00\x00\x00\x18\x00\
                \x00\x00\x00\x00\x00\x00\x1f\x00\x00\x00\x00\x00\x00\x00+\x00\x00\x00\x00\x00\x00\x00:\
                \x00\x00\x00\x00\x00\x00\x00G\x00\x00\x00\x00\x00\x00\x00R\x00\x00\x00\x00\x00\x00\x00_\
                \x00\x00\x00\x00\x00\x00\x00h\x00\x00\x00\x00\x00\x00\x00j\x00\x00\x00\x00\x00\x00\x00\
                \x01a\x01c\x02ab\x03abc\x02MW\x03e\x7ff\x01x\x03q\x02r\x06P12345\x0bNP_000001.1\x0eXP_01\
                1520000.1\x0cabcdefghij\x02!\nabcdefghij\x0cabcdefghij\x02z\x08plain_ok\x01k")),
            ("db.not", Expected::Bytes(b"\
                \x10\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x00\x00\x00\
                \x00\x00\x03\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x00\x00\x00\x00\x05\x00\x00\x00\
                \x00\x00\x00\x00\x06\x00\x00\x00\x00\x00\x00\x00\x07\x00\x00\x00\x00\x00\x00\x00\x08\x00\
                \x00\x00\x00\x00\x00\x00\t\x00\x00\x00\x00\x00\x00\x00\n\x00\x00\x00\x00\x00\x00\x00\x0b\
                \x00\x00\x00\x00\x00\x00\x00\x0c\x00\x00\x00\x00\x00\x00\x00\x0d\x00\x00\x00\x00\x00\x00\
                \x00\x0e\x00\x00\x00\x00\x00\x00\x00\x0f\x00\x00\x00\x00\x00\x00\x00\x10\x00\x00\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\
                \x00")),
            ("db.nsq", Expected::Bytes(b"\
                \x00\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\
                \xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\
                \x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\
                \x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\
                \xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\
                \xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1")),
            ("db.ntf", Expected::Digest(16384, 0xc29f655d5a3815f8)),
            ("db.nto", Expected::Bytes(b"\
                \x10\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x02\x00\x00\x00\x03\x00\x00\x00\x04\x00\
                \x00\x00\x05\x00\x00\x00\x06\x00\x00\x00\x07\x00\x00\x00\x08\x00\x00\x00\t\x00\x00\x00\n\
                \x00\x00\x00\x0b\x00\x00\x00\x0c\x00\x00\x00\x0d\x00\x00\x00\x0e\x00\x00\x00\x0f\x00\x00\
                \x00")),
        ],
    },
    FixCase {
        name: "L_bytes_nucl_v4",
        input: b"\
            >a\x02b t1\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>c\x1fd t2\nACGTTGCAACGTAGCTAGCTA\
            GGATCCATGCATGCAAGCTT\n>ab\xc3\xa9 t3\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>abc\
            \xe9def title\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>MW\xc3\xa9_1 t5\nACGTTGCAACGT\
            AGCTAGCTAGGATCCATGCATGCAAGCTT\n>e\x7ff t6\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>l\
            cl|x t\x01g\x03h u\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>lcl|q\x02r z\nACGTTGCAAC\
            GTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>sp|P12345|ABC_HUMAN h\xe9llo\nACGTTGCAACGTAGCTAGCTAGGA\
            TCCATGCATGCAAGCTT\n>ref|NP_000001.1| ok\x01i\x05j s\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGC\
            AAGCTT\n>XP_011520000.1\x80 raw high\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>lcl|ab\
            cdefghij\x02! k1\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>lcl|abcdefghij k2\nACGTTGC\
            AACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>lcl|abcdefghij\x02z k3\nACGTTGCAACGTAGCTAGCTAGGATCC\
            ATGCATGCAAGCTT\n>plain_ok t\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>k\xff w\nACGTTG\
            CAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n",
        stdin: b"",
        mapdir: false,
        words: &["-in", "{D}/in.fa", "-dbtype", "nucl", "-out", "{D}/db", "-title", "t", "-parse_seqids", "-blastdb_version", "4"],
        status: 0,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences fr\
            om FASTA; added 16 sequences in <S> seconds.\n\n\n",
        stderr: b"\
            Error: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0x2] in string at byt\
            e 8\n\x02b t1\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0x1F] \
            in string at byte 8\n\x1fd t2\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\n\
            Bad char [0x2] in string at byte 8\n\x02b t1\nError: (803.7) [makeblastdb] Blast-def-lin\
            e-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xa9 t3\nError: (803.7) [makebla\
            stdb] Blast-def-line-set.E.title\nBad char [0xA9] in string at byte 9\n\xc3\xa9 t3\nErro\
            r: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0x1F] in string at byte 8\
            \n\x1fd t2\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xE9] in \
            string at byte 8\n\xe9def title\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\
            \nBad char [0xC3] in string at byte 8\n\xc3\xa9 t3\nError: (803.7) [makeblastdb] Blast-d\
            ef-line-set.E.title\nBad char [0xA9] in string at byte 9\n\xc3\xa9 t3\nError: (803.7) [m\
            akeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xa9_1 \
            t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xA9] in string a\
            t byte 9\n\xc3\xa9_1 t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad ch\
            ar [0xE9] in string at byte 8\n\xe9def title\nError: (803.7) [makeblastdb] Blast-def-lin\
            e-set.E.seqid.E.local.str\nBad char [0x7F] in string at byte 23\ne\x7ff\nError: (803.7) \
            [makeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xa9_\
            1 t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xA9] in string\
            \x20at byte 9\n\xc3\xa9_1 t5\nError: (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.\
            local.str\nBad char [0x7F] in string at byte 23\ne\x7ff\nError: (803.7) [makeblastdb] Bl\
            ast-def-line-set.E.seqid.E.local.str\nBad char [0x2] in string at byte 22\nq\x02r\nError\
            : (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xE9] in string at byte 9\
            \nh\xe9llo\nError: (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.local.str\nBad cha\
            r [0x2] in string at byte 22\nq\x02r\nError: (803.7) [makeblastdb] Blast-def-line-set.E.\
            title\nBad char [0xE9] in string at byte 9\nh\xe9llo\nError: (803.7) [makeblastdb] Blast\
            -def-line-set.E.title\nBad char [0x80] in string at byte 8\n\x80 raw high\nError: (803.7\
            ) [makeblastdb] Blast-def-line-set.E.seqid.E.local.str\nBad char [0x2] in string at byte\
            \x2032\nabcdefghij\x02!\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad ch\
            ar [0x80] in string at byte 8\n\x80 raw high\nError: (803.7) [makeblastdb] Blast-def-lin\
            e-set.E.seqid.E.local.str\nBad char [0x2] in string at byte 32\nabcdefghij\x02!\nError: \
            (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.local.str\nBad char [0x2] in string a\
            t byte 32\nabcdefghij\x02z\nError: (803.7) [makeblastdb] Blast-def-line-set.E.seqid.E.lo\
            cal.str\nBad char [0x2] in string at byte 32\nabcdefghij\x02z\nError: (803.7) [makeblast\
            db] Blast-def-line-set.E.title\nBad char [0xFF] in string at byte 8\n\xff w\nError: (803\
            .7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xFF] in string at byte 8\n\xff \
            w\n",
        date: "Oct 30, 2026  4:38 AM",
        files: &[
            ("db.nhr", Expected::Bytes(b"\
                0\x800\x80\xa0\x80\x1a\x05#b t1\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01a\x00\x00\
                \x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\
                \x1a\x05#d t2\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01c\x00\x00\x00\x00\x00\x00\x00\
                \x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x05## t3\x00\x00\
                \xa1\x800\x80\xa0\x80\xa1\x80\x1a\x02ab\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\
                \x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\n#def title\x00\x00\xa1\x800\x80\xa0\
                \x80\xa1\x80\x1a\x03abc\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\
                \x00\x00\x000\x800\x80\xa0\x80\x1a\x07##_1 t5\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\
                \x02MW\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\
                \x80\xa0\x80\x1a\x02t6\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x03e#f\x00\x00\x00\x00\
                \x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01t\
                \x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01x\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\
                \x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01z\x00\x00\xa1\x800\x80\xa0\
                \x80\xa1\x80\x1a\x03q#r\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\
                \x00\x00\x000\x800\x80\xa0\x80\x1a\x05h#llo\x00\x00\xa1\x800\x80\xa7\x800\x80\xa0\x80\
                \x1a\tABC_HUMAN\x00\x00\xa1\x80\x1a\x06P12345\x00\x00\xa2\x80\x1a\x08reviewed\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\
                \xa0\x80\x1a\x02ok\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x0bNP_000001.1\x00\x00\x00\
                \x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\n\
                # raw high\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x0eXP_011520000.1\x00\x00\x00\x00\
                \x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x02k1\
                \x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x0cabcdefghij#!\x00\x00\x00\x00\x00\x00\x00\
                \x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x02k2\x00\x00\xa1\
                \x800\x80\xa0\x80\xa1\x80\x1a\nabcdefghij\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\
                \x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x02k3\x00\x00\xa1\x800\x80\xa0\
                \x80\xa1\x80\x1a\x0cabcdefghij#z\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\
                \x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01t\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\
                \x08plain_ok\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x00\
                0\x800\x80\xa0\x80\x1a\x03# w\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x01k\x00\x00\x00\
                \x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x00")),
            ("db.nin", Expected::Bytes(b"\
                \x00\x00\x00\x04\x00\x00\x00\x00\x00\x00\x00\x01t\x00\x00\x00\x17Oct 30, 2026  4:38 AM\
                \x00\x00\x00\x00\x00\x10\x90\x02\x00\x00\x00\x00\x00\x00\x00\x00\x00)\x00\x00\x00\x00\
                \x00\x00\x00-\x00\x00\x00Z\x00\x00\x00\x88\x00\x00\x00\xbc\x00\x00\x00\xec\x00\x00\x01\
                \x18\x00\x00\x01A\x00\x00\x01l\x00\x00\x01\xbf\x00\x00\x01\xf3\x00\x00\x022\x00\x00\x02g\
                \x00\x00\x02\x9a\x00\x00\x02\xcf\x00\x00\x02\xff\x00\x00\x03*\x00\x00\x00\x01\x00\x00\
                \x00\x0c\x00\x00\x00\x17\x00\x00\x00\"\x00\x00\x00-\x00\x00\x008\x00\x00\x00C\x00\x00\
                \x00N\x00\x00\x00Y\x00\x00\x00d\x00\x00\x00o\x00\x00\x00z\x00\x00\x00\x85\x00\x00\x00\
                \x90\x00\x00\x00\x9b\x00\x00\x00\xa6\x00\x00\x00\xb1\x00\x00\x00\x0c\x00\x00\x00\x17\x00\
                \x00\x00\"\x00\x00\x00-\x00\x00\x008\x00\x00\x00C\x00\x00\x00N\x00\x00\x00Y\x00\x00\x00d\
                \x00\x00\x00o\x00\x00\x00z\x00\x00\x00\x85\x00\x00\x00\x90\x00\x00\x00\x9b\x00\x00\x00\
                \xa6\x00\x00\x00\xb1\x00\x00\x00\xb1")),
            ("db.njs", Expected::Bytes(b"\
                {\n  \"version\": \"1.2\",\n  \"dbname\": \"db\",\n  \"dbtype\": \"Nucleotide\",\n  \"db\
                -version\": 4,\n  \"description\": \"t\",\n  \"number-of-letters\": 656,\n  \"number-of-\
                sequences\": 16,\n  \"last-updated\": \"2026-10-30T04:38:00\",\n  \"number-of-volumes\":\
                \x201,\n  \"bytes-total\": 1750,\n  \"bytes-to-cache\": 437,\n  \"files\": [\n    \"db.n\
                hr\",\n    \"db.nin\",\n    \"db.nog\",\n    \"db.nsd\",\n    \"db.nsi\",\n    \"db.nsq\
                \"\n  ]\n}\n")),
            ("db.nog", Expected::Bytes(b"\
                \x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x10\x00\x00\x00\x00\x00\x00\
                \x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\
                \xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\
                \xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\
                \xff\xff\xff\xff\xff\xff\xff\xff")),
            ("db.nsd", Expected::Bytes(b"\
                a\x020\nab\x022\nabc\x023\nabc_human\x028\nabcdefghij\x02!\x0211\nabcdefghij\x0212\nabcd\
                efghij\x02z\x0213\nc\x021\ne\x7ff\x025\nk\x0215\nlcl|a\x020\nlcl|ab\x022\nlcl|abc\x023\n\
                lcl|abcdefghij\x02!\x0211\nlcl|abcdefghij\x0212\nlcl|abcdefghij\x02z\x0213\nlcl|c\x021\n\
                lcl|e\x7ff\x025\nlcl|k\x0215\nlcl|mw\x024\nlcl|np_000001.1\x029\nlcl|plain_ok\x0214\nlcl\
                |q\x02r\x027\nlcl|x\x026\nlcl|xp_011520000.1\x0210\nmw\x024\nnp_000001.1\x029\np12345\
                \x028\nplain_ok\x0214\nq\x02r\x027\nx\x026\nxp_011520000.1\x0210\n")),
            ("db.nsi", Expected::Bytes(b"\x00\x00\x00\x01\x00\x00\x00\x02\x00\x00\x01_\x00\x00\x00 \x00\x00\x00\x01\x00\x00\x00@\x00\x00\x10\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x01_\x00\x00\x004\x00\x00\x008a\x020\x00")),
            ("db.nsq", Expected::Bytes(b"\
                \x00\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\
                \xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\
                \x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\
                \x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\
                \xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\
                \xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1\x1b\xe4\x1b''(\xd4\xe4\xe4'\xc1")),
        ],
    },
    FixCase {
        name: "L_dup_prot",
        input: b"\
            >M\xc3\xbcller_a strain a\nMKTAYIAKQRQISFVKSHFSRQ\n>M\xc3\xbcller_b strain b\nMKTAYIAKQR\
            QISFVKSHFSRQ\n",
        stdin: b"",
        mapdir: false,
        words: &["-in", "{D}/in.fa", "-dbtype", "prot", "-out", "{D}/db", "-title", "t", "-parse_seqids"],
        status: 1,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from \
            FASTA; added 2 sequences in <S> seconds.\n\n\n",
        stderr: b"\
            Error: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in string at by\
            te 8\n\xc3\xbcller_a strain a\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\n\
            Bad char [0xBC] in string at byte 9\n\xc3\xbcller_a strain a\nError: (803.7) [makeblastd\
            b] Blast-def-line-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xbcller_b strai\
            n b\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xBC] in string \
            at byte 9\n\xc3\xbcller_b strain b\nError: (803.7) [makeblastdb] Blast-def-line-set.E.ti\
            tle\nBad char [0xC3] in string at byte 8\n\xc3\xbcller_a strain a\nError: (803.7) [makeb\
            lastdb] Blast-def-line-set.E.title\nBad char [0xBC] in string at byte 9\n\xc3\xbcller_a \
            strain a\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in st\
            ring at byte 8\n\xc3\xbcller_b strain b\nError: (803.7) [makeblastdb] Blast-def-line-set\
            .E.title\nBad char [0xBC] in string at byte 9\n\xc3\xbcller_b strain b\nBLAST Database c\
            reation error: Error: Duplicate seq_ids are found: \nLCL|M\n\n",
        date: "",
        files: &[
        ],
    },
    FixCase {
        name: "L_dup_prot_v4",
        input: b"\
            >M\xc3\xbcller_a strain a\nMKTAYIAKQRQISFVKSHFSRQ\n>M\xc3\xbcller_b strain b\nMKTAYIAKQR\
            QISFVKSHFSRQ\n",
        stdin: b"",
        mapdir: false,
        words: &["-in", "{D}/in.fa", "-dbtype", "prot", "-out", "{D}/db", "-title", "t", "-parse_seqids", "-blastdb_version", "4"],
        status: 1,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from \
            FASTA; added 2 sequences in <S> seconds.\n\n\n",
        stderr: b"\
            Error: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in string at by\
            te 8\n\xc3\xbcller_a strain a\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\n\
            Bad char [0xBC] in string at byte 9\n\xc3\xbcller_a strain a\nError: (803.7) [makeblastd\
            b] Blast-def-line-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xbcller_b strai\
            n b\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xBC] in string \
            at byte 9\n\xc3\xbcller_b strain b\nError: (803.7) [makeblastdb] Blast-def-line-set.E.ti\
            tle\nBad char [0xC3] in string at byte 8\n\xc3\xbcller_a strain a\nError: (803.7) [makeb\
            lastdb] Blast-def-line-set.E.title\nBad char [0xBC] in string at byte 9\n\xc3\xbcller_a \
            strain a\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in st\
            ring at byte 8\n\xc3\xbcller_b strain b\nError: (803.7) [makeblastdb] Blast-def-line-set\
            .E.title\nBad char [0xBC] in string at byte 9\n\xc3\xbcller_b strain b\nBLAST Database c\
            reation error: Error: Duplicate seq_ids are found: \nLCL|M\n\n",
        date: "",
        files: &[
        ],
    },
    FixCase {
        name: "L_dup_nucl",
        input: b"\
            >M\xc3\xbcller_a strain a\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n>M\xc3\xbcller_b s\
            train b\nACGTTGCAACGTAGCTAGCTAGGATCCATGCATGCAAGCTT\n",
        stdin: b"",
        mapdir: false,
        words: &["-in", "{D}/in.fa", "-dbtype", "nucl", "-out", "{D}/db", "-title", "t", "-parse_seqids"],
        status: 1,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Nucleotide\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences fr\
            om FASTA; added 2 sequences in <S> seconds.\n\n\n",
        stderr: b"\
            Error: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in string at by\
            te 8\n\xc3\xbcller_a strain a\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\n\
            Bad char [0xBC] in string at byte 9\n\xc3\xbcller_a strain a\nError: (803.7) [makeblastd\
            b] Blast-def-line-set.E.title\nBad char [0xC3] in string at byte 8\n\xc3\xbcller_b strai\
            n b\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xBC] in string \
            at byte 9\n\xc3\xbcller_b strain b\nError: (803.7) [makeblastdb] Blast-def-line-set.E.ti\
            tle\nBad char [0xC3] in string at byte 8\n\xc3\xbcller_a strain a\nError: (803.7) [makeb\
            lastdb] Blast-def-line-set.E.title\nBad char [0xBC] in string at byte 9\n\xc3\xbcller_a \
            strain a\nError: (803.7) [makeblastdb] Blast-def-line-set.E.title\nBad char [0xC3] in st\
            ring at byte 8\n\xc3\xbcller_b strain b\nError: (803.7) [makeblastdb] Blast-def-line-set\
            .E.title\nBad char [0xBC] in string at byte 9\n\xc3\xbcller_b strain b\nBLAST Database c\
            reation error: Error: Duplicate seq_ids are found: \nLCL|M\n\n",
        date: "",
        files: &[
        ],
    },
    FixCase {
        name: "T_map_stdin",
        input: b"\
            >ref|NP_000001.1| a\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|lowerid b\nMKTAYIAKQRQISFVKSHFSRQ\n>sp\
            |P12345|X_Y d\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|nomap e\nMKTAYIAKQRQISFVKSHFSRQ\n",
        stdin: b"NP_000001.1 9606\nlowerid 10090\nP12345 562\n",
        mapdir: false,
        words: &["-in", "{D}/in.fa", "-dbtype", "prot", "-out", "{D}/db", "-title", "t", "-parse_seqids", "-taxid_map", "-"],
        status: 0,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from \
            FASTA; added 4 sequences in <S> seconds.\n\n\n",
        stderr: b"",
        date: "Oct 30, 2026  4:38 AM",
        files: &[
            ("db.pdb", Expected::Digest(32768, 0x13360964b0bcf3cc)),
            ("db.phr", Expected::Bytes(b"\
                0\x800\x80\xa0\x80\x1a\x01a\x00\x00\xa1\x800\x80\xa9\x800\x80\xa1\x80\x1a\tNP_000001\x00\
                \x00\xa3\x80\x02\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x02%\x86\
                \x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01b\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\
                \x1a\x07lowerid\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x02'j\x00\x00\x00\x00\x00\
                \x000\x800\x80\xa0\x80\x1a\x01d\x00\x00\xa1\x800\x80\xa7\x800\x80\xa0\x80\x1a\x03X_Y\x00\
                \x00\xa1\x80\x1a\x06P12345\x00\x00\xa2\x80\x1a\x08reviewed\x00\x00\x00\x00\x00\x00\x00\
                \x00\x00\x00\xa2\x80\x02\x02\x022\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01e\x00\
                \x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x05nomap\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\
                \x02\x01\x00\x00\x00\x00\x00\x00\x00")),
            ("db.pin", Expected::Bytes(b"\
                \x00\x00\x00\x05\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x01t\x00\x00\x00\x06db.pdb\
                \x00\x00\x00\x19Oct 30, 2026  4:38 AM\x00\x00\x00\x00\x00\x00\x00\x04X\x00\x00\x00\x00\
                \x00\x00\x00\x00\x00\x00\x16\x00\x00\x00\x00\x00\x00\x00=\x00\x00\x00m\x00\x00\x00\xb7\
                \x00\x00\x00\xe4\x00\x00\x00\x01\x00\x00\x00\x18\x00\x00\x00/\x00\x00\x00F\x00\x00\x00]")),
            ("db.pjs", Expected::Bytes(b"\
                {\n  \"version\": \"1.2\",\n  \"dbname\": \"db\",\n  \"dbtype\": \"Protein\",\n  \"db-ve\
                rsion\": 5,\n  \"description\": \"t\",\n  \"number-of-letters\": 88,\n  \"number-of-sequ\
                ences\": 4,\n  \"last-updated\": \"2026-10-30T04:38:00\",\n  \"number-of-volumes\": 1,\n\
                \x20 \"number-of-taxids\": 4,\n  \"bytes-total\": 49794,\n  \"bytes-to-cache\": 205,\n  \
                \"files\": [\n    \"db.pdb\",\n    \"db.phr\",\n    \"db.pin\",\n    \"db.pog\",\n    \"\
                db.pos\",\n    \"db.pot\",\n    \"db.psq\",\n    \"db.ptf\",\n    \"db.pto\"\n  ]\n}\n")),
            ("db.pog", Expected::Bytes(b"\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x04\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff")),
            ("db.pos", Expected::Bytes(b"\
                \x04\x00\x00\x00\x00\x00\x00\x00\x0c\x00\x00\x00\x00\x00\x00\x00\x14\x00\x00\x00\x00\x00\
                \x00\x00\x1b\x00\x00\x00\x00\x00\x00\x00!\x00\x00\x00\x00\x00\x00\x00\x0bNP_000001.1\x07\
                lowerid\x06P12345\x05nomap")),
            ("db.pot", Expected::Bytes(b"\x04\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x00\x00\x00\x00\x00\x03\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x00\x00\x00\x00\x86%\x00\x00j'\x00\x002\x02\x00\x00\x00\x00\x00\x00")),
            ("db.psq", Expected::Bytes(b"\
                \x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\
                \x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\
                \x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\
                \x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00")),
            ("db.ptf", Expected::Digest(16384, 0x3952fb8b976f9f38)),
            ("db.pto", Expected::Bytes(b"\x01\x00\x00\x00\x03\x00\x00\x00\x01\x00\x00\x00\x02\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x01\x00\x00\x00")),
        ],
    },
    FixCase {
        name: "T_map_stdin_v4",
        input: b"\
            >ref|NP_000001.1| a\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|lowerid b\nMKTAYIAKQRQISFVKSHFSRQ\n>sp\
            |P12345|X_Y d\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|nomap e\nMKTAYIAKQRQISFVKSHFSRQ\n",
        stdin: b"NP_000001.1 9606\nlowerid 10090\nP12345 562\n",
        mapdir: false,
        words: &["-in", "{D}/in.fa", "-dbtype", "prot", "-out", "{D}/db", "-title", "t", "-parse_seqids", "-taxid_map", "-", "-blastdb_version", "4"],
        status: 0,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from \
            FASTA; added 4 sequences in <S> seconds.\n\n\n",
        stderr: b"",
        date: "Oct 30, 2026  4:38 AM",
        files: &[
            ("db.phr", Expected::Bytes(b"\
                0\x800\x80\xa0\x80\x1a\x01a\x00\x00\xa1\x800\x80\xa9\x800\x80\xa1\x80\x1a\tNP_000001\x00\
                \x00\xa3\x80\x02\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x02%\x86\
                \x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01b\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\
                \x1a\x07lowerid\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x02'j\x00\x00\x00\x00\x00\
                \x000\x800\x80\xa0\x80\x1a\x01d\x00\x00\xa1\x800\x80\xa7\x800\x80\xa0\x80\x1a\x03X_Y\x00\
                \x00\xa1\x80\x1a\x06P12345\x00\x00\xa2\x80\x1a\x08reviewed\x00\x00\x00\x00\x00\x00\x00\
                \x00\x00\x00\xa2\x80\x02\x02\x022\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01e\x00\
                \x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\x05nomap\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\
                \x02\x01\x00\x00\x00\x00\x00\x00\x00")),
            ("db.pin", Expected::Bytes(b"\
                \x00\x00\x00\x04\x00\x00\x00\x01\x00\x00\x00\x01t\x00\x00\x00\x17Oct 30, 2026  4:38 AM\
                \x00\x00\x00\x00\x00\x04X\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x16\x00\x00\x00\x00\
                \x00\x00\x00=\x00\x00\x00m\x00\x00\x00\xb7\x00\x00\x00\xe4\x00\x00\x00\x01\x00\x00\x00\
                \x18\x00\x00\x00/\x00\x00\x00F\x00\x00\x00]")),
            ("db.pjs", Expected::Bytes(b"\
                {\n  \"version\": \"1.2\",\n  \"dbname\": \"db\",\n  \"dbtype\": \"Protein\",\n  \"db-ve\
                rsion\": 4,\n  \"description\": \"t\",\n  \"number-of-letters\": 88,\n  \"number-of-sequ\
                ences\": 4,\n  \"last-updated\": \"2026-10-30T04:38:00\",\n  \"number-of-volumes\": 1,\n\
                \x20 \"bytes-total\": 616,\n  \"bytes-to-cache\": 189,\n  \"files\": [\n    \"db.phr\",\
                \n    \"db.pin\",\n    \"db.pog\",\n    \"db.psd\",\n    \"db.psi\",\n    \"db.psq\"\n  \
                ]\n}\n")),
            ("db.pog", Expected::Bytes(b"\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x04\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff")),
            ("db.psd", Expected::Bytes(b"\
                lcl|lowerid\x021\nlcl|nomap\x023\nlowerid\x021\nnomap\x023\nnp_000001\x020\nnp_000001.1\
                \x020\np12345\x022\nx_y\x022\n")),
            ("db.psi", Expected::Bytes(b"\
                \x00\x00\x00\x01\x00\x00\x00\x02\x00\x00\x00U\x00\x00\x00\x08\x00\x00\x00\x01\x00\x00\
                \x00@\x00\x00\x10\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00U\x00\
                \x00\x004\x00\x00\x00Blcl|lowerid\x021\x00")),
            ("db.psq", Expected::Bytes(b"\
                \x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\
                \x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\
                \x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\
                \x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00")),
        ],
    },
    FixCase {
        name: "T_map_and_in_stdin",
        input: b"",
        stdin: b"\
            NP_000001.1 9606\nlowerid 10090\nP12345 562\n>ref|NP_000001.1| a\nMKTAYIAKQRQISFVKSHFSRQ\
            \n>lcl|lowerid b\nMKTAYIAKQRQISFVKSHFSRQ\n>sp|P12345|X_Y d\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl\
            |nomap e\nMKTAYIAKQRQISFVKSHFSRQ\n",
        mapdir: false,
        words: &["-in", "-", "-dbtype", "prot", "-out", "{D}/db", "-title", "t", "-parse_seqids", "-taxid_map", "-"],
        status: 255,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\n",
        stderr: b"\
            Error: NCBI C++ Exception:\n    T0 \"/opt/conda/conda-bld/blast_1754905306566/work/c++/s\
            rc/corelib/ncbistr.cpp\", line 861: Error: (CStringException::eConvert) ncbi::NStr::Stri\
            ngToInt8() - Cannot convert string 'a' to Int8 (m_Pos = 0)\n\n",
        date: "",
        files: &[
        ],
    },
    FixCase {
        name: "T_map_only_stdin",
        input: b"",
        stdin: b"NP_000001.1 9606\nlowerid 10090\nP12345 562\n",
        mapdir: false,
        words: &["-in", "-", "-dbtype", "prot", "-out", "{D}/db", "-title", "t", "-parse_seqids", "-taxid_map", "-"],
        status: 0,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\n\nNo volumes were crea\
            ted.\n\n",
        stderr: b"\
            Error: [makeblastdb] No sequences matched any of the taxids provided.\n",
        date: "",
        files: &[
        ],
    },
    FixCase {
        name: "T_map_only_stdin_noeol",
        input: b"",
        stdin: b"NP_000001.1 9606\nlowerid 10090\nP12345 562",
        mapdir: false,
        words: &["-in", "-", "-dbtype", "prot", "-out", "{D}/db", "-title", "t", "-parse_seqids", "-taxid_map", "-"],
        status: 0,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\n",
        stderr: b"\
            Error: [makeblastdb] No sequences matched any of the taxids provided.\n",
        date: "",
        files: &[
        ],
    },
    FixCase {
        name: "T_map_dir",
        input: b"\
            >ref|NP_000001.1| a\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|lowerid b\nMKTAYIAKQRQISFVKSHFSRQ\n>sp\
            |P12345|X_Y d\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|nomap e\nMKTAYIAKQRQISFVKSHFSRQ\n",
        stdin: b"",
        mapdir: true,
        words: &["-in", "{D}/in.fa", "-dbtype", "prot", "-out", "{D}/db", "-title", "t", "-parse_seqids", "-taxid_map", "{D}/mapdir"],
        status: 0,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\nAdding sequences from \
            FASTA; added 4 sequences in <S> seconds.\n\n\n",
        stderr: b"",
        date: "Oct 30, 2026  4:38 AM",
        files: &[
            ("db.pdb", Expected::Digest(32768, 0x13360964b0bcf3cc)),
            ("db.phr", Expected::Bytes(b"\
                0\x800\x80\xa0\x80\x1a\x01a\x00\x00\xa1\x800\x80\xa9\x800\x80\xa1\x80\x1a\tNP_000001\x00\
                \x00\xa3\x80\x02\x01\x01\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\
                \x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01b\x00\x00\xa1\x800\x80\xa0\x80\xa1\x80\x1a\
                \x07lowerid\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\
                \x800\x80\xa0\x80\x1a\x01d\x00\x00\xa1\x800\x80\xa7\x800\x80\xa0\x80\x1a\x03X_Y\x00\x00\
                \xa1\x80\x1a\x06P12345\x00\x00\xa2\x80\x1a\x08reviewed\x00\x00\x00\x00\x00\x00\x00\x00\
                \x00\x00\xa2\x80\x02\x01\x00\x00\x00\x00\x00\x00\x000\x800\x80\xa0\x80\x1a\x01e\x00\x00\
                \xa1\x800\x80\xa0\x80\xa1\x80\x1a\x05nomap\x00\x00\x00\x00\x00\x00\x00\x00\xa2\x80\x02\
                \x01\x00\x00\x00\x00\x00\x00\x00")),
            ("db.pin", Expected::Bytes(b"\
                \x00\x00\x00\x05\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x01t\x00\x00\x00\x06db.pdb\
                \x00\x00\x00\x19Oct 30, 2026  4:38 AM\x00\x00\x00\x00\x00\x00\x00\x04X\x00\x00\x00\x00\
                \x00\x00\x00\x00\x00\x00\x16\x00\x00\x00\x00\x00\x00\x00<\x00\x00\x00k\x00\x00\x00\xb4\
                \x00\x00\x00\xe1\x00\x00\x00\x01\x00\x00\x00\x18\x00\x00\x00/\x00\x00\x00F\x00\x00\x00]")),
            ("db.pjs", Expected::Bytes(b"\
                {\n  \"version\": \"1.2\",\n  \"dbname\": \"db\",\n  \"dbtype\": \"Protein\",\n  \"db-ve\
                rsion\": 5,\n  \"description\": \"t\",\n  \"number-of-letters\": 88,\n  \"number-of-sequ\
                ences\": 4,\n  \"last-updated\": \"2026-10-30T04:38:00\",\n  \"number-of-volumes\": 1,\n\
                \x20 \"bytes-total\": 49779,\n  \"bytes-to-cache\": 205,\n  \"files\": [\n    \"db.pdb\"\
                ,\n    \"db.phr\",\n    \"db.pin\",\n    \"db.pog\",\n    \"db.pos\",\n    \"db.pot\",\n\
                \x20   \"db.psq\",\n    \"db.ptf\",\n    \"db.pto\"\n  ]\n}\n")),
            ("db.pog", Expected::Bytes(b"\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x04\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff\xff")),
            ("db.pos", Expected::Bytes(b"\
                \x04\x00\x00\x00\x00\x00\x00\x00\x0c\x00\x00\x00\x00\x00\x00\x00\x14\x00\x00\x00\x00\x00\
                \x00\x00\x1b\x00\x00\x00\x00\x00\x00\x00!\x00\x00\x00\x00\x00\x00\x00\x0bNP_000001.1\x07\
                lowerid\x06P12345\x05nomap")),
            ("db.pot", Expected::Bytes(b"\x04\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x00\x00\x00\x00\x02\x00\x00\x00\x00\x00\x00\x00\x03\x00\x00\x00\x00\x00\x00\x00\x04\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00")),
            ("db.psq", Expected::Bytes(b"\
                \x00\x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\
                \x0c\n\x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\
                \x12\x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00\x0c\n\x12\
                \x01\x16\t\x01\n\x0f\x10\x0f\t\x11\x06\x13\n\x11\x08\x06\x11\x10\x0f\x00")),
            ("db.ptf", Expected::Digest(16384, 0xc29f655d5a3815f8)),
            ("db.pto", Expected::Bytes(b"\x04\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00\x02\x00\x00\x00\x03\x00\x00\x00")),
        ],
    },
    FixCase {
        name: "T_map_empty_name",
        input: b"\
            >ref|NP_000001.1| a\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|lowerid b\nMKTAYIAKQRQISFVKSHFSRQ\n>sp\
            |P12345|X_Y d\nMKTAYIAKQRQISFVKSHFSRQ\n>lcl|nomap e\nMKTAYIAKQRQISFVKSHFSRQ\n",
        stdin: b"",
        mapdir: false,
        words: &["-in", "{D}/in.fa", "-dbtype", "prot", "-out", "{D}/db", "-title", "t", "-parse_seqids", "-taxid_map", ""],
        status: 1,
        stdout: b"\
            \n\nBuilding a new DB, current time: <T>\nNew DB name:   <D>/db\nNew DB title:  t\nSeque\
            nce type: Protein\nKeep MBits: T\nMaximum file size: 3000000000B\n",
        stderr: b"\
            Command line argument error: Argument \"taxid_map\". File is not accessible\n",
        date: "",
        files: &[
        ],
    },
];
