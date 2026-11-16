//! NCBI's stack as its reader of a PSSM fills it: whether NCBI BLAST+ 2.17.0 reads a PSSM whose
//! values nest in themselves, or may overflow its stack and crash, in a run of this command line,
//! environment and stack limit.
//!
//! NCBI reads a PSSM with `CObjectIStreamAsn`, whose functions recurse once per open value, on the
//! main thread's stack. A type that nests in itself can fill it: the reader overflows the stack and
//! the run ends with SIGSEGV before any output (gdb: the deepest frames are `CObjectIStreamAsn`
//! functions under `CPsiBlastArgs::ExtractAlgorithmOptions`; line RG D3-1 fix round 5,
//! `r5/bisect/gdb_*.txt`). The main thread's stack holds `RLIMIT_STACK` bytes (whole pages) from
//! its end, the kernel's layout of the command line, the environment and the auxiliary vector
//! takes the first of them, and the start goes down by a random 0 to 8,191 bytes
//! (`arch_align_stack`): NCBI reads a nest of `n` values of a type when the bytes left,
//! `room - top`, hold `base + per_value * n`, and may crash from the depth at which they do not
//! with the lowest start. LOSAT refuses from that depth (no result of NCBI's can be checked there)
//! and reads every shallower one.
//!
//! Measured (line RG D3-1 fix round 6, `r6/cal/`): NCBI 2.17.0 under `setarch -R` with an 8 MiB
//! stack, for each type at three depths (12 levels apart), the least `room - top` at which it reads,
//! by an environment variable's length (16-byte steps; the third depth read at its predicted bytes
//! and crashed 16 bytes below, for every type and both programs). Values of several types open at
//! once add their bytes per value (`r6/cal/check.*.txt`: Bioseq-sets around a descriptor's
//! User-objects and around a feature's Seq-loc mixes, Seq-align `disc` levels around a `std` segment's
//! Seq-loc mixes, and the three together: each read at its predicted bytes and crashed 16 below).

use super::table::named;

/// NCBI's bytes for the values of one type that nests in itself: `per_value` for each value open,
/// and `base` for everything else when the nest sits where it was measured (one type nested in
/// itself in a query Bioseq: a Bioseq-set: the query entry; Seq-loc-mix and Seq-loc-equiv: a
/// feature's location; User-object, Pub, OrgName: a descriptor; Seq-align: an annotation, its
/// innermost value a `std` segment of two Seq-loc intervals; Variation-ref: a feature;
/// SeqTable-multi-data: a Seq-table column's `int-scaled` data).
#[derive(Clone, Copy, Debug)]
struct Frames {
    ty: u16,
    per_value: i64,
    base: i64,
}

/// The frames of one program's reader (`Frames` of each recursive type).
#[derive(Debug)]
pub(crate) struct NcbiReader {
    frames: [Frames; 9],
}

/// tblastn's frames below the reader take 704 bytes fewer than psiblast's (every type, `r6/cal`).
const TBLASTN_LESS: i64 = 704;

const fn reader(less: i64) -> NcbiReader {
    const fn f(ty: u16, per_value: i64, base: i64, less: i64) -> Frames {
        Frames {
            ty,
            per_value,
            base: base - less,
        }
    }
    NcbiReader {
        frames: [
            f(named::BIOSEQ_SET, 1648, 13072, less),
            f(named::SEQ_LOC_MIX, 1216, 11920, less),
            f(named::SEQ_LOC_EQUIV, 1216, 11920, less),
            f(named::USER_OBJECT, 2016, 9424, less),
            f(named::SEQ_ALIGN, 1456, 12912, less),
            f(named::PUB, 1104, 12000, less),
            f(named::VARIATION_REF, 2112, 9776, less),
            f(named::SEQTABLE_MULTI_DATA, 880, 11360, less),
            f(named::ORGNAME, 1856, 10240, less),
        ],
    }
}

/// psiblast's reader (`r6/cal/summary.tsv`).
const PSIBLAST: NcbiReader = reader(0);
/// tblastn's reader.
const TBLASTN: NcbiReader = reader(TBLASTN_LESS);

/// NCBI's stack in a run: one program's reader and the bytes it may take.
#[derive(Clone, Copy, Debug)]
pub(crate) struct NcbiStack {
    reader: &'static NcbiReader,
    room: Room,
}

#[derive(Clone, Copy, Debug)]
enum Room {
    /// The room of this run (`room_of_this_run`).
    ThisRun,
    /// A room of the tests.
    #[cfg(test)]
    Bytes(i64),
}

/// psiblast's stack in this run.
pub(crate) const PSIBLAST_STACK: NcbiStack = NcbiStack {
    reader: &PSIBLAST,
    room: Room::ThisRun,
};

/// tblastn's stack in this run.
pub(crate) const TBLASTN_STACK: NcbiStack = NcbiStack {
    reader: &TBLASTN,
    room: Room::ThisRun,
};

impl NcbiStack {
    /// psiblast's or tblastn's reader with `room` bytes (the tests).
    #[cfg(test)]
    pub(crate) const fn with_room(tblastn: bool, room: i64) -> NcbiStack {
        NcbiStack {
            reader: if tblastn { &TBLASTN } else { &PSIBLAST },
            room: Room::Bytes(room),
        }
    }

    fn room(&self) -> i64 {
        match self.room {
            Room::ThisRun => {
                static ROOM: std::sync::OnceLock<i64> = std::sync::OnceLock::new();
                *ROOM.get_or_init(|| Layout::of_this_run().room())
            }
            #[cfg(test)]
            Room::Bytes(bytes) => bytes,
        }
    }
}

/// The values of NCBI's recursive types open in the reader, and the bytes of NCBI's stack they
/// take; the room is found once, when the first such value opens (a PSSM that psiblast writes has
/// none).
pub(super) struct Depth {
    stack: NcbiStack,
    room: Option<i64>,
    /// The `Frames` of each recursive value open, outermost first.
    open: Vec<u8>,
    counts: [u32; 9],
    /// `per_value` of the values open.
    bytes: i64,
    /// The position in `open` of the outermost value that is not a Bioseq-set.
    outer: Option<usize>,
}

impl Depth {
    pub(super) fn new(stack: NcbiStack) -> Depth {
        Depth {
            stack,
            room: None,
            open: Vec::new(),
            counts: [0; 9],
            bytes: 0,
            outer: None,
        }
    }

    /// A value of the type `ty` opens: `Ok(true)` for a recursive type (`close` it), the refusal's
    /// text when NCBI may overflow its stack with it open. The bytes of the values open are their
    /// `per_value` and the `base` of the outermost one that is not a Bioseq-set (the nests inside
    /// Bioseq-sets and Seq-loc mixes inside Seq-aligns' segments add their bytes to the outer nest's
    /// as measured), or a Bioseq-set's when only those are open.
    pub(super) fn open(&mut self, ty: u16) -> Result<bool, String> {
        let Some(at) = self
            .stack
            .reader
            .frames
            .iter()
            .position(|frames| frames.ty == ty)
        else {
            return Ok(false);
        };
        let frames = self.stack.reader.frames[at];
        self.counts[at] += 1;
        self.bytes += frames.per_value;
        if self.outer.is_none() && ty != named::BIOSEQ_SET {
            self.outer = Some(self.open.len());
        }
        self.open.push(at as u8);
        let base = match self.outer {
            Some(outer) => self.stack.reader.frames[usize::from(self.open[outer])].base,
            None => frames.base,
        };
        let room = *self.room.get_or_insert_with(|| self.stack.room());
        if base + self.bytes > room {
            return Err(format!(
                "a PSSM whose values nest so deep that NCBI BLAST+ may overflow its stack and crash with this run's stack limit and environment ({} {} values open at once)",
                self.counts[at],
                super::resolved(ty).name()
            ));
        }
        Ok(true)
    }

    /// The innermost recursive value closes.
    pub(super) fn close(&mut self) {
        let at = usize::from(self.open.pop().expect("a recursive value open"));
        self.counts[at] -= 1;
        self.bytes -= self.stack.reader.frames[at].per_value;
        if self.outer == Some(self.open.len()) {
            self.outer = None;
        }
    }
}

/// The lowest start of the main thread's stack below its layout: `arch_align_stack` takes
/// `get_random_int() % 8192` bytes.
// Linux 5.15 (the kernel of the measurements), arch/x86/kernel/process.c, `arch_align_stack`:
// ```c
// unsigned long arch_align_stack(unsigned long sp)
// {
// 	if (!(current->personality & ADDR_NO_RANDOMIZE) && randomize_va_space)
// 		sp -= get_random_int() % 8192;
// 	return sp & ~0xf;
// }
// ```
const LOWEST_START: u64 = 8191;

/// The default stack limit (`ulimit -s` 8192), where the run cannot tell its own.
const DEFAULT_LIMIT: u64 = 8 << 20;

/// The main thread's stack with no limit (`RLIMIT_STACK` infinity): the gap the kernel leaves below
/// the stack's end, at most 5/6 of the address space (`mmap_base`, arch/x86/mm/mmap.c:93-109).
const UNLIMITED: u64 = ((1 << 47) - 4096) / 6 * 5;

/// The auxiliary vector of Linux 5.15 on x86_64 in words, with its `AT_NULL` (20 entries), where
/// the run cannot read its own.
const AUXV_WORDS: u64 = 42;

/// What the kernel lays out at the end of the main thread's stack when a program starts, and the
/// stack's limit: the strings of the command line and environment, their counts and the auxiliary
/// vector (NCBI's process of this run has the same; the program's path and words differ from
/// `losat`'s by tens of bytes, which LOSAT cannot know).
#[derive(Debug, PartialEq)]
struct Layout {
    /// `RLIMIT_STACK`'s soft limit.
    limit: u64,
    /// From the stack's end to the first byte of `argv[0]`: a pointer's 8 bytes, the program's
    /// path (`AT_EXECFN`), the environment's strings and the arguments'.
    strings: u64,
    argc: u64,
    envc: u64,
    /// The auxiliary vector's words, with the `AT_NULL` pair.
    auxv_words: u64,
}

impl Layout {
    /// The bytes from the stack's end to the program's start (`mm->start_stack`) for a start
    /// `random` bytes lower.
    // Linux 5.15, fs/binfmt_elf.c, `create_elf_tables`:
    // ```c
    // 	p = arch_align_stack(p);
    // 	...
    // 		u_platform = (elf_addr_t __user *)STACK_ALLOC(p, len);
    // 	...
    // 	u_rand_bytes = (elf_addr_t __user *)
    // 		       STACK_ALLOC(p, sizeof(k_rand_bytes));
    // 	...
    // 	ei_index = elf_info - (elf_addr_t *)mm->saved_auxv;
    // 	sp = STACK_ADD(p, ei_index);
    //
    // 	items = (argc + 1) + (envc + 1) + 1;
    // 	bprm->p = STACK_ROUND(sp, items);
    // ```
    // (`ELF_PLATFORM` is "x86_64" with its NUL, `k_rand_bytes` 16 bytes; `STACK_ADD(sp, items)`
    // is `sp - items` in 8-byte words, `STACK_ROUND(sp, items)` `(sp - items) & ~15`; fs/exec.c
    // copies the strings below a pointer at the stack's end, the path first.)
    fn top(&self, random: u64) -> u64 {
        const END: u64 = 1 << 47;
        let p = END.saturating_sub(self.strings).saturating_sub(random) & !15;
        let p = p.saturating_sub(7 + 16 + 8 * self.auxv_words);
        let items = (self.argc + 1) + (self.envc + 1) + 1;
        END - (p.saturating_sub(8 * items) & !15)
    }

    /// The bytes that NCBI's reader may take with the lowest start.
    fn room(&self) -> i64 {
        let limit = (self.limit.min(UNLIMITED) & !4095) as i64;
        limit - self.top(LOWEST_START) as i64
    }

    /// The layout of this run: on Linux as the kernel made it (`/proc/self`), else (or when that
    /// cannot be read) of the arguments and environment the program gets, the default limit and
    /// Linux's auxiliary vector (WASI, which has neither a stack limit nor a stack of this form).
    fn of_this_run() -> Layout {
        #[cfg(target_os = "linux")]
        if let Some(layout) = Self::of_proc() {
            return layout;
        }
        let argv: Vec<_> = std::env::args_os().collect();
        let envs: Vec<_> = std::env::vars_os().collect();
        let len = |s: &std::ffi::OsStr| s.as_encoded_bytes().len() as u64;
        let path = argv.first().map_or(0, |s| len(s));
        Layout {
            limit: DEFAULT_LIMIT,
            strings: 8
                + path
                + 1
                + envs.iter().map(|(k, v)| len(k) + len(v) + 2).sum::<u64>()
                + argv.iter().map(|s| len(s) + 1).sum::<u64>(),
            argc: argv.len() as u64,
            envc: envs.len() as u64,
            auxv_words: AUXV_WORDS,
        }
    }

    #[cfg(target_os = "linux")]
    fn of_proc() -> Option<Layout> {
        let limits = std::fs::read_to_string("/proc/self/limits").ok()?;
        let soft = limits
            .lines()
            .find_map(|line| line.strip_prefix("Max stack size"))?
            .split_whitespace()
            .next()?;
        let limit = if soft == "unlimited" {
            UNLIMITED
        } else {
            soft.parse().ok()?
        };
        let maps = std::fs::read_to_string("/proc/self/maps").ok()?;
        let end = maps
            .lines()
            .find(|line| line.ends_with("[stack]"))?
            .split(['-', ' '])
            .nth(1)?;
        let end = u64::from_str_radix(end, 16).ok()?;
        // `(48) arg_start` of proc(5), after the command's name in parentheses.
        let stat = std::fs::read_to_string("/proc/self/stat").ok()?;
        let arg_start: u64 = stat[stat.rfind(')')? + 2..]
            .split(' ')
            .nth(45)?
            .parse()
            .ok()?;
        let environ = std::fs::read("/proc/self/environ").ok()?;
        let auxv = std::fs::read("/proc/self/auxv").ok()?;
        let words: Vec<u64> = auxv
            .chunks_exact(8)
            .map(|word| u64::from_ne_bytes(word.try_into().expect("8 bytes")))
            .collect();
        let auxv_words = words.chunks_exact(2).position(|pair| pair[0] == 0)? as u64 * 2 + 2;
        Some(Layout {
            limit,
            strings: end.checked_sub(arg_start)?,
            argc: std::env::args_os().count() as u64,
            envc: environ.iter().filter(|&&byte| byte == 0).count() as u64,
            auxv_words,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The kernel's layout (`setarch -R python3 -I r6/tools/stacktop.py --self`: `vm_end -
    /// start_stack` of the process, 6,192 and, with 1,234 more bytes in one more variable, 7,440).
    #[test]
    fn the_layout_is_the_kernels() {
        let layout = |strings, envc| Layout {
            limit: DEFAULT_LIMIT,
            strings,
            argc: 4,
            envc,
            auxv_words: 42,
        };
        assert_eq!(layout(5283, 60).top(0), 6192);
        assert_eq!(layout(6522, 61).top(0), 7440);
        assert_eq!(layout(5283, 60).top(LOWEST_START), 14384);
        assert_eq!(
            layout(5283, 60).room(),
            8388608 - 14384,
            "the room with the lowest start"
        );
        let unlimited = Layout {
            limit: u64::MAX,
            ..layout(5283, 60)
        };
        assert_eq!(unlimited.room(), (UNLIMITED & !4095) as i64 - 14384);
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn this_runs_layout_is_read() {
        let layout = Layout::of_proc().expect("/proc/self");
        assert!(layout.strings > 8 && layout.argc >= 1 && layout.auxv_words >= 2);
    }

    /// NCBI's measured thresholds (`r6/cal/summary.tsv`: depth, least bytes at which NCBI read):
    /// each is the model's `base + per_value * n`, `n` the values of the type open at once.
    #[test]
    fn the_frames_give_ncbis_measured_thresholds() {
        // (type, values in the measured nest of depth d: d + offset, d, least bytes (psiblast),
        // least bytes (tblastn) at depth d + 1 for tblastn where its crash depth was one more)
        let measured: [(u16, i64, i64, i64, i64, i64); 9] = [
            (named::BIOSEQ_SET, 0, 5076, 8378320, 5077, 8379264),
            (named::SEQ_LOC_MIX, 0, 6881, 8379216, 6881, 8378512),
            (named::SEQ_LOC_EQUIV, 0, 6881, 8379216, 6881, 8378512),
            (named::USER_OBJECT, 1, 4150, 8377840, 4151, 8379152),
            (named::SEQ_ALIGN, 1, 5745, 8379088, 5745, 8378384),
            (named::PUB, 1, 7578, 8379216, 7578, 8378512),
            (named::VARIATION_REF, 1, 3961, 8377520, 3962, 8378928),
            (named::SEQTABLE_MULTI_DATA, 1, 9508, 8379280, 9508, 8378576),
            (named::ORGNAME, 1, 4508, 8378944, 4508, 8378240),
        ];
        for (ty, offset, d, least, d_tbn, least_tbn) in measured {
            for (tblastn, d, least) in [(false, d, least), (true, d_tbn, least_tbn)] {
                let n = d + offset;
                for (room, reads) in [(least, true), (least - 16, false)] {
                    let mut depth = Depth::new(NcbiStack::with_room(tblastn, room));
                    let opened = (0..n).map(|_| depth.open(ty)).collect::<Vec<_>>();
                    assert_eq!(
                        opened.iter().all(Result::is_ok),
                        reads,
                        "{} tblastn {tblastn} room {room}",
                        super::super::type_info(ty).name()
                    );
                    if !reads {
                        assert!(opened[..opened.len() - 1].iter().all(Result::is_ok));
                    }
                }
            }
        }
    }

    /// Values of several types open at once (`r6/cal/check.psiblast.txt`): 2,537 Bioseq-sets around
    /// 2,075 User-objects need 8,373,600 bytes, 2,001 Seq-aligns around 4,485 Seq-loc mixes
    /// 8,380,128 (the Seq-align's base: its innermost `std` segment holds the mixes).
    #[test]
    fn nests_of_several_types_add_their_bytes() {
        let opens = |room: i64, nests: &[(u16, u32)]| {
            let mut depth = Depth::new(NcbiStack::with_room(false, room));
            nests
                .iter()
                .all(|&(ty, n)| (0..n).all(|_| depth.open(ty).is_ok()))
        };
        let user = [(named::BIOSEQ_SET, 2537), (named::USER_OBJECT, 2075)];
        assert!(opens(8373600, &user));
        assert!(!opens(8373600 - 16, &user));
        let align = [(named::SEQ_ALIGN, 2001), (named::SEQ_LOC_MIX, 4485)];
        assert!(opens(8380128, &align));
        assert!(!opens(8380128 - 16, &align));
        let three = [
            (named::BIOSEQ_SET, 1500),
            (named::SEQ_ALIGN, 1501),
            (named::SEQ_LOC_MIX, 3051),
        ];
        assert!(opens(8380384, &three));
        assert!(!opens(8380384 - 16, &three));
        // A value closes: its bytes go.
        let mut depth = Depth::new(NcbiStack::with_room(false, 13072 + 1648));
        assert_eq!(depth.open(named::BIOSEQ_SET), Ok(true));
        assert!(depth.open(named::BIOSEQ_SET).is_err());
        depth.close();
        depth.close();
        assert_eq!(depth.open(named::BIOSEQ_SET), Ok(true));
        assert_eq!(depth.open(named::SEQ_ID), Ok(false));
    }
}
