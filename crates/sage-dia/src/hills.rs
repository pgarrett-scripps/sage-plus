//! Compact, m/z-sorted hill storage for one acquisition channel.
//!
//! A channel is either the MS1 stream or one DIA isolation window. Hills keep
//! only what rescoring reads (m/z, cycle span, apex, intensity profile), and
//! every profile lives in one shared `f32` arena, so a channel costs about
//! 24 bytes per hill plus 4 bytes per profile point. Cycle indices are local to
//! the channel: cycle `i` is the channel's `i`-th scan, at `rts[i]`.

use koth_core::Hill;

/// A hill reduced to what rescoring needs.
#[derive(Clone, Copy, Debug)]
pub struct CompactHill {
    pub mz: f32,
    pub intensity_max: f32,
    /// Intensity-weighted ion mobility (1/K0); 0 when the run has none.
    pub im: f32,
    /// First cycle of the profile.
    pub start: u32,
    /// Cycle of the highest point.
    pub apex: u32,
    /// Offset of the profile in the channel arena.
    offset: u32,
    /// Number of cycles in the profile.
    pub len: u32,
}

impl CompactHill {
    /// Last cycle of the profile (inclusive).
    pub fn end(&self) -> u32 {
        self.start + self.len - 1
    }
}

/// Hills of one channel, sorted by m/z, with the channel's cycle times.
pub struct Channel {
    /// Isolation window bounds in m/z (`0..inf` for MS1).
    pub lower: f64,
    pub upper: f64,
    /// Ion-mobility bounds (1/K0) of a diaPASEF box; `0..inf` without IM.
    pub im_lower: f64,
    pub im_upper: f64,
    /// Retention time (minutes) of each cycle of this channel.
    pub rts: Vec<f32>,
    hills: Vec<CompactHill>,
    mzs: Vec<f32>,
    /// Hill indices sorted by apex cycle, for time-first lookups.
    by_apex: Vec<u32>,
    arena: Vec<f32>,
}

impl Channel {
    /// Build a channel from koth hills whose `scan_start` / `scan_apex` are
    /// channel-local cycle indices (true for both koth MS1 and per-window MS2
    /// detection).
    pub fn from_hills(lower: f64, upper: f64, rts: Vec<f32>, hills: Vec<Hill>) -> Self {
        let mut compact = Vec::with_capacity(hills.len());
        let mut arena = Vec::with_capacity(hills.iter().map(|h| h.intensity_profile.len()).sum());
        for hill in hills {
            if hill.intensity_profile.is_empty() {
                continue;
            }
            compact.push(CompactHill {
                mz: hill.mz as f32,
                intensity_max: hill.intensity_max as f32,
                im: hill.im as f32,
                start: hill.scan_start as u32,
                apex: hill.scan_apex as u32,
                offset: arena.len() as u32,
                len: hill.intensity_profile.len() as u32,
            });
            arena.extend_from_slice(&hill.intensity_profile);
        }
        compact.sort_by(|a, b| a.mz.total_cmp(&b.mz));
        let mzs = compact.iter().map(|h| h.mz).collect();
        let mut by_apex: Vec<u32> = (0..compact.len() as u32).collect();
        by_apex.sort_by_key(|&i| compact[i as usize].apex);
        Channel {
            lower,
            upper,
            im_lower: 0.0,
            im_upper: f64::INFINITY,
            rts,
            hills: compact,
            mzs,
            by_apex,
            arena,
        }
    }

    /// Restrict the channel to an ion-mobility range (a diaPASEF box).
    pub fn with_im(mut self, im_lower: f64, im_upper: f64) -> Self {
        self.im_lower = im_lower;
        self.im_upper = im_upper;
        self
    }

    /// Does this channel isolate a precursor at `mz` (and `im`, when known)?
    pub fn contains(&self, mz: f64, im: f64) -> bool {
        self.lower <= mz
            && mz <= self.upper
            && (im == 0.0 || (self.im_lower <= im && im <= self.im_upper))
    }

    pub fn len(&self) -> usize {
        self.hills.len()
    }

    pub fn is_empty(&self) -> bool {
        self.hills.is_empty()
    }

    /// Heap bytes held by this channel.
    pub fn heap_bytes(&self) -> usize {
        self.hills.capacity() * std::mem::size_of::<CompactHill>()
            + self.mzs.capacity() * 4
            + self.by_apex.capacity() * 4
            + self.arena.capacity() * 4
            + self.rts.capacity() * 4
    }

    /// Hills whose m/z is within `ppm` of `mz`.
    pub fn find(&self, mz: f32, ppm: f32) -> &[CompactHill] {
        let tol = mz * ppm * 1e-6;
        let lo = self.mzs.partition_point(|&m| m < mz - tol);
        let hi = self.mzs.partition_point(|&m| m <= mz + tol);
        &self.hills[lo..hi]
    }

    /// Hills whose apex cycle lies in `lo..=hi`.
    pub fn apex_between(&self, lo: i64, hi: i64) -> impl Iterator<Item = &CompactHill> {
        let apex = |i: &u32| self.hills[*i as usize].apex as i64;
        let a = self.by_apex.partition_point(|i| apex(i) < lo);
        let b = self.by_apex.partition_point(|i| apex(i) <= hi);
        self.by_apex[a..b].iter().map(|&i| &self.hills[i as usize])
    }

    /// Like [`Channel::apex_between`], with each hill's index in the channel.
    pub fn apex_between_indexed(
        &self,
        lo: i64,
        hi: i64,
    ) -> impl Iterator<Item = (usize, &CompactHill)> {
        let apex = |i: &u32| self.hills[*i as usize].apex as i64;
        let a = self.by_apex.partition_point(|i| apex(i) < lo);
        let b = self.by_apex.partition_point(|i| apex(i) <= hi);
        self.by_apex[a..b]
            .iter()
            .map(|&i| (i as usize, &self.hills[i as usize]))
    }

    /// All hills, in m/z order (indices match `apex_between_indexed`).
    pub fn hills(&self) -> &[CompactHill] {
        &self.hills
    }

    /// Intensity profile of `hill`, one value per cycle from `hill.start`.
    pub fn profile(&self, hill: &CompactHill) -> &[f32] {
        &self.arena[hill.offset as usize..(hill.offset + hill.len) as usize]
    }

    /// Intensity of `hill` at channel cycle `cycle` (0 outside the hill).
    pub fn intensity(&self, hill: &CompactHill, cycle: i64) -> f32 {
        if cycle < hill.start as i64 || cycle > hill.end() as i64 {
            return 0.0;
        }
        self.arena[(hill.offset as i64 + cycle - hill.start as i64) as usize]
    }

    /// Cycle whose retention time is closest to `rt`.
    pub fn cycle_at(&self, rt: f32) -> usize {
        let i = self.rts.partition_point(|&t| t < rt);
        if i == 0 {
            0
        } else if i >= self.rts.len() {
            self.rts.len() - 1
        } else if (self.rts[i] - rt) < (rt - self.rts[i - 1]) {
            i
        } else {
            i - 1
        }
    }
}
