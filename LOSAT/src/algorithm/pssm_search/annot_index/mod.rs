//! The annotations of a PSSM's query as NCBI's object manager takes them when psiblast and
//! tblastn add the PSSM's Bioseq to their scope (line RG D3-1, `om/SPEC.md`): the Seq-tables'
//! column checks when each annotation is added (phase C, `CSeqTableInfo`), then the index of every
//! annotation in order (phase I, `CSeq_annot_Info::x_InitAnnotKeys`), with NCBI's messages, its
//! exceptions, and the failed diagnostic stream after its "Failed to parse location of" message
//! (`ncbi_exception::main_diag_failbit`). What NCBI does not do the same way on every run (a
//! message with a memory address, a crash, BitMagic's decoding, a signed overflow, the order of
//! memory addresses) is LOSAT's refusal, each naming its condition (`Refusal`).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:535-537 (tblastn_app.cpp:222-226)
//! ```c++
//! 		_ASSERT(pssm->HasQuery());
//! 		_ASSERT(pssm->GetQuery().IsSeq());  // single query only!
//! 		scope->AddTopLevelSeqEntry(pssm->SetQuery());
//! ```
//! NCBI reference (598d8ae6): c++/src/objmgr/scope_impl.cpp:230-238 (`CScope_Impl::AddSeq_entry`:
//! the TSE is made, each annotation's `CSeq_annot_Info` with it, then `x_ClearCacheOnNewData`
//! asks for the TSE's annotation Seq-ids, `CTSE_Info::GetAnnotIds`, which indexes them all,
//! `UpdateAnnotIndex`, tse_info.cpp:709, 1035)
//! ```c++
//!     CTSE_Lock tse_lock = ds_lock->GetDataSource().AddStaticTSE(entry);
//!     x_ClearCacheOnNewData(*tse_lock);
//! ```

mod align;
mod range;
mod table;
mod view;

use range::{Handle, HandleRangeMap, Mode, Range};
use view::V;

use crate::ncbi_exception::{post_line_in, NcbiException, Severity};
use crate::report::serial::asn_types::{table::named, Value};

/// What NCBI's object manager made of the query's annotations: its messages (standard error) in
/// order, whether its last left the main thread's diagnostic stream failed, and the exception
/// that ends the run (`CATCH_ALL`'s report, already in the failed stream's form when it is).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct AnnotOutcome {
    pub(crate) posts: Vec<u8>,
    pub(crate) muted: bool,
    pub(crate) abort: Option<Abort>,
}

/// `CATCH_ALL`'s report and exit status of an exception of the object manager.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Abort {
    pub(crate) exit: i32,
    pub(crate) message: String,
}

/// What LOSAT refuses: an annotation that NCBI's object manager does not index the same way on
/// every run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Refusal {
    /// `Empty region in <address> ...` (1204.1-.4, .7, .8) before the stream failed.
    EmptyRegion,
    /// An empty region of a feature, graph, alignment, location set or table location after the
    /// stream failed: NCBI prints the annotation's ASN.1 text, cut at its 4096-byte buffer.
    MutedDump,
    /// A null Seq-id handle or location of a table row before the stream failed.
    NullPointer,
    /// A Dense-seg whose division by `numseg` is by zero (SIGFPE).
    DivideByZero,
    /// A bit-bvector or bit-set-bvector that the index reads (BitMagic).
    BitVector,
    /// A Packed-seg read past its `present`, `starts`, `ids` or `lens`.
    PackedOutOfBounds,
    /// An int-delta sum that overflows its accumulator.
    DeltaOverflow,
    /// An int-scaled value that overflows `Int8`.
    ScaledOverflow,
    /// A `field-id` whose difference with the column group's base overflows `int`.
    FieldIdOverflow,
    /// Two Seq-ids of one table row's location with different messages whose order is that of
    /// memory addresses (one `m_Packed`).
    IdOrder,
    /// Two `giim` Seq-ids of one number in one location (or one discontinuous alignment), which
    /// NCBI takes for one Seq-id or two by the `giim` Seq-ids it holds at the time.
    GiimIds,
    /// A feature table of more rows than `table::ROWS_LIMIT`.
    TableRows,
    /// A Std-seg `dim` above `align::STD_DIM_LIMIT`.
    StdDim,
}

impl Refusal {
    /// The refusal's text, after "a PSSM query whose annot holds ...", in LOSAT's form
    /// ("... is not supported by LOSAT's <PROGRAM>", `in_pssm::PssmReader::refuse`).
    pub(crate) fn text(self) -> &'static str {
        match self {
            Refusal::EmptyRegion => "a PSSM query whose annot has an empty location (NCBI's object manager posts 'Empty region in' with a memory address)",
            Refusal::MutedDump => "a PSSM query whose annot has an empty location after NCBI's 'Failed to parse location' message (NCBI then prints the annotation's ASN.1 text, which LOSAT's writer does not write yet)",
            Refusal::NullPointer => "a PSSM query whose annot has a Seq-table row without a Seq-id or a location (NCBI's object manager stops with a NULL pointer report of memory addresses)",
            Refusal::DivideByZero => "a PSSM query whose annot has a Dense-seg of no segments with starts or strands (NCBI's object manager divides by zero and crashes)",
            Refusal::BitVector => "a PSSM query whose annot has a Seq-table bit-bvector or bit-set-bvector that NCBI's object manager reads (BitMagic data decoded without checks)",
            Refusal::PackedOutOfBounds => "a PSSM query whose annot has a Packed-seg with fewer present, starts, ids or lens than its dim and numseg (NCBI's object manager reads past them)",
            Refusal::DeltaOverflow => "a PSSM query whose annot has a Seq-table int-delta column whose sums leave their integer type (a signed overflow in NCBI's object manager)",
            Refusal::ScaledOverflow => "a PSSM query whose annot has a Seq-table int-scaled column whose values leave Int8 (a signed overflow in NCBI's object manager)",
            Refusal::FieldIdOverflow => "a PSSM query whose annot has a Seq-table column field-id below -2147483638 (a signed overflow in NCBI's object manager)",
            Refusal::IdOrder => "a PSSM query whose annot has a Seq-table row location whose Seq-ids get different messages and have no different gi, accession or general id numbers (NCBI's object manager writes them in the order of memory addresses)",
            Refusal::GiimIds => "a PSSM query whose annot has a location with two giim Seq-ids of one id number (NCBI's object manager takes them for one Seq-id or two depending on the giim Seq-ids it holds at that time)",
            Refusal::TableRows => "a PSSM query whose annot has a feature Seq-table of more than 10000000 rows (NCBI's object manager runs out of memory depending on the machine)",
            Refusal::StdDim => "a PSSM query whose annot has a Std-seg dim above 1000000 (NCBI's object manager allocates its rows depending on the machine's memory)",
        }
    }
}

/// What stops the index: an exception of NCBI's (`CException`, a `std::exception`), the null
/// pointer of a table row (`Refusal::NullPointer` before the stream failed, a report of nothing
/// after it), or LOSAT's refusal.
#[derive(Debug)]
pub(super) enum Stop {
    Exception(NcbiException),
    Std(&'static str),
    NullPointer,
    Refuse(Refusal),
}

impl Stop {
    fn exception(e: NcbiException) -> Stop {
        Stop::Exception(e)
    }

    fn std_exception(what: &'static str) -> Stop {
        Stop::Std(what)
    }
}

/// The message of an empty overlapping range: the annotation's kind
/// (`seq_annot_info.cpp` `ERR_POST_X(1..8, "Empty region in " ...)`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Event {
    /// 1204.1, a feature (with its ASN.1 text).
    Feature,
    /// 1204.2, a graph.
    Graph,
    /// 1204.3, an alignment.
    Align,
    /// 1204.4, a location set.
    Locs,
    /// 1204.7, a table's location.
    TableLocation,
    /// 1204.8, a feature table's row (no ASN.1 text).
    Row,
}

/// The state of the index: the program's diagnostic prefix, the messages so far and whether the
/// stream has failed.
pub(super) struct Indexer<'p> {
    program: &'p str,
    posts: Vec<u8>,
    muted: bool,
    caches: table::Caches,
}

/// The object manager's index of the annotations `annots` (`Seq-annot` values) of the query of a
/// PSSM read by `program` (`psiblast`, `tblastn`: the diagnostic prefix), or LOSAT's refusal.
pub(crate) fn index_annots(program: &str, annots: &[Value]) -> Result<AnnotOutcome, Refusal> {
    let mut ix = Indexer {
        program,
        posts: Vec::new(),
        muted: false,
        caches: table::Caches::default(),
    };
    let annots: Vec<V<'_>> = annots
        .iter()
        .map(|annot| V::new(named::SEQ_ANNOT, annot))
        .collect();
    let abort = match ix.run(&annots) {
        Ok(()) => None,
        Err(Stop::Refuse(refusal)) => return Err(refusal),
        Err(Stop::Exception(e)) => Some(crate::cli::NativeError::exception(&e)),
        Err(Stop::Std(what)) => Some(crate::cli::NativeError::std_exception(what)),
        Err(Stop::NullPointer) if ix.muted => Some(crate::cli::NativeError {
            exit: 255,
            message: String::new(),
        }),
        Err(Stop::NullPointer) => return Err(Refusal::NullPointer),
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:248-251, 260-263
    // ```c++
    //     catch (const CException& e) {                                           \
    //         LOG_POST(Error << "Error: " << e.what());                           \
    //         exit_code = BLAST_UNKNOWN_ERROR;                                    \
    //     }                                                                       \
    // ...
    //     catch (const std::exception& e) {                                       \
    //         LOG_POST(Error << "Error: " << e.what());                           \
    //         exit_code = BLAST_UNKNOWN_ERROR;                                    \
    //     }                                                                       \
    // ```
    // After the stream failed, the `LOG_POST` writes its line feed alone (`cli::catch_all_text`,
    // which the run applies once the stream is failed; the text kept here is that line feed).
    let abort = abort.map(|error| Abort {
        exit: error.exit,
        message: if ix.muted {
            "\n".to_string()
        } else {
            error.message
        },
    });
    Ok(AnnotOutcome {
        posts: ix.posts,
        muted: ix.muted,
        abort,
    })
}

impl Indexer<'_> {
    /// A post of the object manager (`ERR_POST_X`, `ERR_POST`; `id` empty without an error code),
    /// in the stream's state.
    fn post(&mut self, severity: Severity, id: &str, text: &[u8]) {
        let id = (!id.is_empty()).then_some(id);
        self.posts.extend(post_line_in(
            severity,
            id,
            Some(self.program),
            text,
            self.muted,
        ));
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:403-427, 474-500 (phase C:
    // `x_InitAnnotList` of each annotation in order makes a Seq-table's `CSeqTableInfo`; the other
    // kinds only list their objects) and :639-668 (phase I, `x_InitAnnotKeys`)
    // ```c++
    //     C_Data& data = m_Object->SetData();
    //     switch ( data.Which() ) {
    //     case C_Data::e_Ftable:
    //         x_InitFeatKeys(tse);
    //         break;
    //     case C_Data::e_Align:
    //         x_InitAlignKeys(tse);
    //         break;
    //     case C_Data::e_Graph:
    //         x_InitGraphKeys(tse);
    //         break;
    //     case C_Data::e_Locs:
    //         x_InitLocsKeys(tse);
    //         break;
    //     case C_Data::e_Seq_table:
    //         x_InitFeatTableKeys(tse);
    //         break;
    //     default:
    //         break;
    //     }
    // ```
    fn run(&mut self, annots: &[V<'_>]) -> Result<(), Stop> {
        let data: Vec<Option<(&'static str, V<'_>)>> = annots
            .iter()
            .map(|annot| annot.get("data").map(V::choice))
            .collect();
        let mut tables = Vec::with_capacity(annots.len());
        for kind in &data {
            tables.push(match kind {
                Some(("seq-table", table)) => Some(self.table_info(*table)?),
                _ => None,
            });
        }
        for ((annot, kind), table) in annots.iter().zip(&data).zip(&tables) {
            match kind {
                Some(("ftable", feats)) => {
                    for feat in feats.items() {
                        self.index_feature(feat)?;
                    }
                }
                Some(("align", aligns)) => {
                    for align in aligns.items() {
                        self.index_align(align)?;
                    }
                }
                Some(("graph", graphs)) => {
                    for graph in graphs.items() {
                        let loc = graph.get("loc").expect("mandatory");
                        let map = HandleRangeMap::of_location(loc, Mode::None);
                        self.refuse_empty(std::slice::from_ref(&map), Event::Graph)?;
                    }
                }
                Some(("locs", locs)) => self.index_locs(*annot, *locs)?,
                Some(("seq-table", _)) => self.index_table(table.as_ref().expect("phase C"))?,
                // `ids` (and an annotation without data) are not indexed.
                _ => {}
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/annot_object.cpp:426-447 (`x_ProcessFeat`)
    // ```c++
    //     hrmaps.resize(feat.IsSetProduct()? 2: 1);
    //     ...
    //     CHandleRangeMap::ETransSplicing mode = CHandleRangeMap::eNoTransSplicing;
    //     if ( feat.IsSetExcept_text() ) {
    //         if ( feat.GetExcept_text().find("trans-splicing") != NPOS ) {
    //             mode = CHandleRangeMap::eTransSplicing;
    //         }
    //         else if ( feat.GetExcept_text().find("circular RNA") != NPOS ) {
    //             mode = CHandleRangeMap::eCircularRNA;
    //         }
    //     }
    //     hrmaps[0].AddLocation(feat.GetLocation(), mode);
    //     if ( feat.IsSetProduct() ) {
    //         ...
    //         hrmaps[1].AddLocation(feat.GetProduct(), mode);
    //     }
    // ```
    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:741-786 (`x_InitFeatKeys`)
    // ```c++
    //         ITERATE ( vector<CHandleRangeMap>, hrmit, hrmaps ) {
    //             ...
    //             ITERATE ( CHandleRangeMap, hrit, *hrmit ) {
    //                 const CHandleRange& hr = hrit->second;
    //                 key.m_Range = hr.GetOverlappingRange();
    //                 if ( key.m_Range.Empty() ) {
    //                     ERR_POST_X(1, "Empty region in "<<GetDescription()<<" "<<
    //                                MSerial_AsnText<<*info.GetFeatFast());
    //                     continue;
    //                 }
    //                 ...
    //                 if ( hr.HasGaps() ) {
    //                     ...
    //                     if ( hr.IsCircular() ) {
    //                         key.m_Range = hr.GetCircularRangeStart();
    //                         x_Map(mapper, key, index);
    //                         key.m_Range = hr.GetCircularRangeEnd();
    //                     }
    //                 }
    //                 ...
    //                 x_Map(mapper, key, index);
    //             }
    //         }
    // ```
    fn index_feature(&mut self, feat: V<'_>) -> Result<(), Stop> {
        let except = feat.get("except-text").map(V::str).unwrap_or_default();
        let mode = if contains(except, b"trans-splicing") {
            Mode::TransSplicing
        } else if contains(except, b"circular RNA") {
            Mode::CircularRna
        } else {
            Mode::None
        };
        for loc in [feat.get("location"), feat.get("product")]
            .into_iter()
            .flatten()
        {
            let map = HandleRangeMap::of_location(loc, mode);
            self.map_location(&map, Event::Feature)?;
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:845-915 (`x_AddAlignKeys`,
    // `x_AddAlignIndex`)
    // ```c++
    //     if ( align.GetSegs().IsDisc() ) {
    //         // discontinuous alignments should be indexed separately
    //         // add keys for each sub-alignment
    //         ...
    //         for ( auto& align_ref : sub_aligns ) {
    //             vector<CHandleRangeMap> hrmaps;
    //             CAnnotObject_Info::x_ProcessAlign(hrmaps, *align_ref, master);
    //             for ( auto& hrm : hrmaps ) {
    //                 ITERATE ( CHandleRangeMap, hrit, hrm ) {
    //                     total_ranges_map[hrit->first] += hrit->second.GetOverlappingRange();
    //                 }
    //             }
    //             all_hrmaps.push_back(hrmaps);
    //         }
    //         auto hrmap_iter = all_hrmaps.begin();
    //         for ( auto& align_ref : sub_aligns ) {
    //             x_AddAlignIndex(info, *align_ref, mapper, *hrmap_iter++, &total_ranges_map);
    //         }
    //     }
    //     else {
    //         vector<CHandleRangeMap> hrmaps;
    //         CAnnotObject_Info::x_ProcessAlign(hrmaps, align, master);
    //         x_AddAlignIndex(info, align, mapper, hrmaps);
    //     }
    // ...
    //     ITERATE ( vector<CHandleRangeMap>, hrmit, hrmaps ) {
    //         ITERATE ( CHandleRangeMap, hrit, *hrmit ) {
    //             const CHandleRange& hr = hrit->second;
    //             SAnnotObject_Key key;
    //             key.m_Range = total_ranges_map? total_ranges_map->find(hrit->first)->second: hr.GetOverlappingRange();
    //             if ( key.m_Range.Empty() ) {
    //                 ERR_POST_X(3, "Empty region in "<<GetDescription()<<
    //                            " "<<MSerial_AsnText<<info.GetAlign());
    //                 continue;
    //             }
    //             ...
    //             x_Map(mapper, key, index);
    //         }
    //         ...
    //     }
    // ```
    // The mapped range is never empty here (no "Failed to parse location"): an alignment posts
    // only its sizes' warnings, or its empty regions (refused), in any order of its Seq-ids. The
    // Seq-ids of one row, and for a discontinuous alignment of all its rows, are one handle as
    // NCBI's trees take them (`Handle`).
    fn index_align(&mut self, align: V<'_>) -> Result<(), Stop> {
        let segs = align.get("segs").map(V::choice);
        if let Some(("disc", subs)) = segs {
            let mut all = Vec::new();
            for sub in subs.items() {
                let mut rows = Vec::new();
                align::process_align(self, &mut rows, sub)?;
                all.push(rows);
            }
            if HandleRangeMap::giim_repeats(all.iter().flatten()) {
                return Err(Stop::Refuse(Refusal::GiimIds));
            }
            let mut total: Vec<(&Handle, Range)> = Vec::new();
            for map in all.iter().flatten() {
                for (handle, hr) in &map.map {
                    let range = hr.overlapping();
                    match total.iter_mut().find(|(h, _)| *h == handle) {
                        Some((_, sum)) => sum.combine(range),
                        None => total.push((handle, range)),
                    }
                }
            }
            if total.iter().any(|(_, range)| range.is_empty()) {
                return Err(self.empty_region());
            }
        } else {
            let mut rows = Vec::new();
            align::process_align(self, &mut rows, align)?;
            for map in &rows {
                self.refuse_empty(std::slice::from_ref(map), Event::Align)?;
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:940-975 (`x_InitLocsKeys`)
    // ```c++
    //     // Only one referenced location per annot is allowed
    //     if ( m_ObjectIndex.GetInfos().size() != 1) {
    //         return;
    //     }
    //     ...
    //     info.GetMaps(hrmaps, master);
    //     ...
    //             if ( key.m_Range.Empty() ) {
    //                 ERR_POST_X(4, "Empty region in "<<GetDescription()<<" "<<
    //                            MSerial_AsnText<<info.GetLocs());
    //                 continue;
    //             }
    // ```
    // NCBI reference (598d8ae6): c++/src/objmgr/annot_object.cpp:396-416 (`GetMaps` of `locs`: the
    // first `region` of the annotation's `desc`, not the location itself)
    // ```c++
    //         if ( !annot.IsSetDesc() ) {
    //             break;
    //         }
    //         CConstRef<CSeq_loc> region;
    //         ITERATE(CSeq_annot::TDesc::Tdata, desc_it, annot.GetDesc().Get()) {
    //             if ( (*desc_it)->IsRegion() ) {
    //                 region.Reset(&(*desc_it)->GetRegion());
    //                 break;
    //             }
    //         }
    //         if ( region ) {
    //             hrmaps.resize(1);
    //             ...
    //             hrmaps[0].AddLocation(*region);
    //         }
    // ```
    fn index_locs(&mut self, annot: V<'_>, locs: V<'_>) -> Result<(), Stop> {
        if locs.len() != 1 {
            return Ok(());
        }
        let region = annot.get("desc").and_then(|desc| {
            desc.items()
                .map(V::choice)
                .find(|(variant, _)| *variant == "region")
                .map(|(_, region)| region)
        });
        if let Some(region) = region {
            let map = HandleRangeMap::of_location(region, Mode::None);
            self.refuse_empty(std::slice::from_ref(&map), Event::Locs)?;
        }
        Ok(())
    }

    /// An empty region of an annotation whose index maps only non-empty ranges (a graph, an
    /// alignment's row, a location set; `event` names its message): refused, in any order of the
    /// Seq-ids.
    fn refuse_empty(&mut self, maps: &[HandleRangeMap], event: Event) -> Result<(), Stop> {
        debug_assert!(matches!(event, Event::Graph | Event::Align | Event::Locs));
        if HandleRangeMap::giim_repeats(maps) {
            return Err(Stop::Refuse(Refusal::GiimIds));
        }
        let empty =
            |map: &HandleRangeMap| map.map.iter().any(|(_, hr)| hr.overlapping().is_empty());
        if maps.iter().any(empty) {
            return Err(self.empty_region());
        }
        Ok(())
    }

    /// `ERR_POST_X(n, "Empty region in " << GetDescription() ...)` of an empty region: before the stream
    /// failed, a message with the TSE's memory address (`CTSE_Info::GetDescription`,
    /// tse_info.cpp:2049-2056); after it, NCBI's ASN.1 text of the object (refused both), or for a
    /// table row its header alone (`Indexer::map_location`).
    fn empty_region(&self) -> Stop {
        Stop::Refuse(if self.muted {
            Refusal::MutedDump
        } else {
            Refusal::EmptyRegion
        })
    }

    /// The null pointer of a table row (`Stop::NullPointer`).
    fn null_pointer(&self) -> Stop {
        Stop::NullPointer
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:741-786 (quoted at
    // `index_feature`), 1048-1080 and 1119-1150 (a table's location and a table row's location:
    // the same loop with `ERR_POST_X(7, ...)` and `ERR_POST_X(8, "Empty region in "
    // <<GetDescription()<< " Seq-table row "<<row)`)
    // Each Seq-id of the map in NCBI's order (`HandleRangeMap::ncbi_order`): its overlapping range,
    // or for a circular one its two ranges, mapped (`x_map`); an empty one posts "Empty region"
    // (refused: a memory address, or after the stream failed the object's ASN.1 text), except for
    // a table row after the stream failed, whose post is its header alone, `Error: (1204.8)`.
    // Where NCBI's order is that of memory addresses (Seq-ids of one `m_Packed`: `m_Info <
    // handle.m_Info`, seq_id_handle.hpp:206, quoted at `HandleRangeMap::ncbi_order`), only such a
    // row's order changes the messages: its empty Seq-id's post (`Error: (1204.8)` after the
    // stream failed, "Empty region in <address>" before) and another's "Failed to parse location"
    // (`Error: (1204.6)`) in one group are refused (`om/` cases `s_o09_tie_gi5_acc5_u` and
    // `s_o10_tie_acc5v1_acc5v2_u`: psiblast and tblastn wrote them in different orders); any other
    // messages of a group are the same in every order (the 1204.6 messages differ only by how many
    // came before, an empty region of another kind is refused in every order).
    fn map_location(&mut self, map: &HandleRangeMap, event: Event) -> Result<(), Stop> {
        if HandleRangeMap::giim_repeats([map]) {
            return Err(Stop::Refuse(Refusal::GiimIds));
        }
        let empty = |hr: &range::HandleRange| hr.overlapping().is_empty();
        let fails = |hr: &range::HandleRange| {
            !empty(hr)
                && hr.has_gaps()
                && hr.circular
                && (hr.circular_start().is_empty() || hr.circular_end().is_empty())
        };
        for group in map.ncbi_order() {
            if event == Event::Row
                && group.iter().any(|&at| empty(&map.map[at].1))
                && group.iter().any(|&at| fails(&map.map[at].1))
            {
                return Err(Stop::Refuse(Refusal::IdOrder));
            }
            for at in group {
                let (handle, hr) = &map.map[at];
                let range = hr.overlapping();
                if range.is_empty() {
                    if event != Event::Row || !self.muted {
                        return Err(self.empty_region());
                    }
                    // After the stream failed: the header alone.
                    self.post(Severity::Error, "1204.8", b"");
                    continue;
                }
                if hr.has_gaps() && hr.circular {
                    self.x_map(hr.circular_start(), Some(handle.clone()))?;
                    self.x_map(hr.circular_end(), Some(handle.clone()))?;
                } else {
                    self.x_map(range, Some(handle.clone()))?;
                }
            }
        }
        Ok(())
    }

    // NCBI reference (598d8ae6): c++/src/objmgr/seq_annot_info.cpp:672-700
    // ```c++
    // bool CSeq_annot_Info::x_Map(const CTSEAnnotObjectMapper& mapper,
    //                             const SAnnotObject_Key& key,
    //                             const SAnnotObject_Index& index)
    // {
    //     if ( key.m_Range.Empty() ) {
    //         ...
    //         ERR_POST_X(6, "Failed to parse location of "<<s.rdbuf()<<" in "<<GetDescription());
    //         return false;
    //     }
    //     bool new_id = mapper.Map(key, index);
    //     m_ObjectIndex.AddMap(key, index);
    //     return new_id;
    // }
    // ```
    // The message stops at `s.rdbuf()`, which fails the stream (`ncbi_exception::
    // main_diag_failbit`). `mapper.Map` of a null handle (a table row without a Seq-id) asks
    // `CSeq_id_Handle::HaveMatchingHandles()` of it: `CCoreException::eNullPtr` (`Stop::NullPointer`).
    /// `x_Map(range, handle)`: `None` is a null handle.
    fn x_map(&mut self, range: Range, handle: Option<Handle>) -> Result<(), Stop> {
        if range.is_empty() {
            self.post(Severity::Error, "1204.6", b"Failed to parse location of ");
            self.muted = true;
            return Ok(());
        }
        if handle.is_none() {
            return Err(self.null_pointer());
        }
        Ok(())
    }
}

/// `string::find(needle) != NPOS`.
fn contains(haystack: &[u8], needle: &[u8]) -> bool {
    haystack
        .windows(needle.len())
        .any(|window| window == needle)
}

#[cfg(test)]
mod tests;
