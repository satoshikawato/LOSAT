//! The mean of the unit counts in a window (`CSeqMaskerScoreMean`), the score of the "mean"
//! trigger that blastn and the windowmasker application use.

use super::istat::UnitCounts;
use super::window::Window;
use super::WinMaskError;

pub(super) struct ScoreMean {
    sum: u32,
    start: u32,
    num: u32,
    scores: Vec<u32>,
    scores_start: usize,
}

impl ScoreMean {
    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_score_mean.cpp:83-90
    pub(super) fn new(window: &Window<'_>, ustat: &UnitCounts) -> Result<Self, WinMaskError> {
        let num = window.num_units() as u32;
        let mut score = Self {
            sum: 0,
            start: window.start(),
            num,
            scores: vec![0; num as usize],
            scores_start: 0,
        };
        score.fill_scores(window, ustat)?;
        Ok(score)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_score_mean.cpp:48-49
    #[inline]
    pub(super) fn value(&self) -> u32 {
        self.sum / self.num
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_score_mean.cpp:62-80
    /// `Uint4` arithmetic: the sum wraps as NCBI's does (it never needs to with real counts).
    pub(super) fn post_advance(
        &mut self,
        step: u32,
        window: &Window<'_>,
        ustat: &UnitCounts,
    ) -> Result<(), WinMaskError> {
        if step == 1 && window.unit_step() == 1 && window.start().wrapping_sub(self.start) == 1 {
            let slot = self.scores_start;
            self.sum = self.sum.wrapping_sub(self.scores[slot]);
            self.scores[slot] = ustat.at(window.unit(self.num as usize - 1))?;
            self.sum = self.sum.wrapping_add(self.scores[slot]);
            self.scores_start = if slot == self.num as usize - 1 {
                0
            } else {
                slot + 1
            };
            self.start = window.start();
            Ok(())
        } else {
            self.fill_scores(window, ustat)
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/winmask/seq_masker_score_mean.cpp:93-109
    fn fill_scores(&mut self, window: &Window<'_>, ustat: &UnitCounts) -> Result<(), WinMaskError> {
        self.sum = 0;
        self.scores_start = 0;
        for i in 0..self.num as usize {
            self.scores[i] = ustat.at(window.unit(i))?;
            self.sum = self.sum.wrapping_add(self.scores[i]);
        }
        self.start = window.start();
        Ok(())
    }
}
