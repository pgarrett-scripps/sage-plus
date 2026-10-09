//! MS2 spectrum assembly for Bruker TDF input.
//!
//! `sage-plus-tdf` only reads frames and metadata tables. Everything here
//! reproduces the spectrum reader of timsrust 0.6.6 (`timsrust-tdf`), which
//! Sage Plus used through Beta 16, so search results do not change with the
//! reader:
//!
//! - ddaPASEF: one spectrum per precursor, merging the precursor's
//!   `PasefFrameMsMsInfo` scan ranges across frames.
//! - diaPASEF: one spectrum per (frame, window split), with the window groups of
//!   `DiaFrameMsMsWindows` expanded by `bruker_config.ms2.frame_splitting_params`.
//!
//! Each spectrum's TOF indices are summed per index, smoothed and centroided in
//! integer TOF space, exactly as timsrust does.

use sage_plus_tdf::{
    DiaFrameMsMsInfo, DiaFrameMsMsWindow, LinearMobilityScale, PasefFrameMsMsInfo,
};

/// A quadrupole isolation window: timsrust's `IsolationWindow`, with the same
/// floating-point operations so derived centers and widths match bit for bit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsolationWindow {
    pub lower: f64,
    pub upper: f64,
    pub collision_energy: f64,
}

impl IsolationWindow {
    pub fn from_center(center: f64, width: f64, collision_energy: f64) -> Self {
        let half_width = width / 2.0;
        Self {
            lower: center - half_width,
            upper: center + half_width,
            collision_energy,
        }
    }

    pub fn from_bounds(lower: f64, upper: f64, collision_energy: f64) -> Self {
        Self {
            lower,
            upper,
            collision_energy,
        }
    }

    pub fn center(&self) -> f64 {
        (self.lower + self.upper) / 2.0
    }

    pub fn width(&self) -> f64 {
        self.upper - self.lower
    }
}

/// Sum intensities of equal TOF indices; the result is sorted by TOF index.
pub fn group_and_sum(tof_indices: Vec<u32>, intensities: Vec<u64>) -> (Vec<u32>, Vec<u64>) {
    let mut points: Vec<(u32, u64)> = tof_indices.into_iter().zip(intensities).collect();
    // Integer sums do not depend on the order within a TOF index.
    points.sort_unstable_by_key(|&(tof, _)| tof);
    let mut tofs: Vec<u32> = Vec::with_capacity(points.len());
    let mut sums: Vec<u64> = Vec::with_capacity(points.len());
    for (tof, intensity) in points {
        match tofs.last() {
            Some(&last) if last == tof => *sums.last_mut().expect("parallel vectors") += intensity,
            _ => {
                tofs.push(tof);
                sums.push(intensity);
            }
        }
    }
    (tofs, sums)
}

/// Add the intensity of every neighbor within `window` TOF bins, using the
/// unsmoothed intensities. `tof_indices` must be sorted.
pub fn smooth(tof_indices: &[u32], intensities: &[u64], window: u32) -> Vec<u64> {
    let mut smoothed = intensities.to_vec();
    for (current, &current_tof) in tof_indices.iter().enumerate() {
        for (next, &next_tof) in tof_indices.iter().enumerate().skip(current + 1) {
            if next_tof - current_tof > window {
                break;
            }
            smoothed[current] += intensities[next];
            smoothed[next] += intensities[current];
        }
    }
    smoothed
}

/// Keep local maxima: of two points within `window` TOF bins, the smaller is
/// dropped (the later one on ties). `tof_indices` must be sorted.
pub fn centroid(tof_indices: &[u32], intensities: &[u64], window: u32) -> (Vec<u32>, Vec<u64>) {
    let mut keep = vec![true; tof_indices.len()];
    for (current, &current_tof) in tof_indices.iter().enumerate() {
        for (next, &next_tof) in tof_indices.iter().enumerate().skip(current + 1) {
            if next_tof - current_tof > window {
                break;
            }
            if intensities[current] < intensities[next] {
                keep[current] = false;
            } else {
                keep[next] = false;
            }
        }
    }
    let tofs = tof_indices
        .iter()
        .zip(&keep)
        .filter_map(|(&tof, &keep)| keep.then_some(tof))
        .collect();
    let values = intensities
        .iter()
        .zip(&keep)
        .filter_map(|(&value, &keep)| keep.then_some(value))
        .collect();
    (tofs, values)
}

/// Sum, smooth and centroid the raw points of one spectrum.
pub fn process(
    tof_indices: Vec<u32>,
    intensities: Vec<u64>,
    smoothing_window: u32,
    centroiding_window: u32,
) -> (Vec<u32>, Vec<u64>) {
    let (tofs, intensities) = group_and_sum(tof_indices, intensities);
    let intensities = smooth(&tofs, &intensities, smoothing_window);
    centroid(&tofs, &intensities, centroiding_window)
}

/// ddaPASEF spectra: the `PasefFrameMsMsInfo` row indices of each precursor,
/// in precursor order and, within a precursor, in table order. Rows without a
/// precursor are left out.
pub fn dda_groups(rows: &[PasefFrameMsMsInfo]) -> Vec<(u64, Vec<usize>)> {
    let mut order: Vec<(u64, usize)> = rows
        .iter()
        .enumerate()
        .filter_map(|(index, row)| row.precursor.map(|precursor| (precursor, index)))
        .collect();
    order.sort_by_key(|&(precursor, _)| precursor);
    let mut groups: Vec<(u64, Vec<usize>)> = Vec::new();
    for (precursor, index) in order {
        match groups.last_mut() {
            Some((last, indices)) if *last == precursor => indices.push(index),
            _ => groups.push((precursor, vec![index])),
        }
    }
    groups
}

/// The isolation windows of one diaPASEF window group, sorted by first scan.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct QuadGroup {
    pub scan_starts: Vec<usize>,
    pub scan_ends: Vec<usize>,
    pub windows: Vec<IsolationWindow>,
}

/// Window groups `1..=max(WindowGroup)`, indexed by group - 1. Groups without
/// rows are empty.
pub fn quad_groups(rows: &[DiaFrameMsMsWindow]) -> Result<Vec<QuadGroup>, String> {
    let count = rows.iter().map(|row| row.window_group).max().unwrap_or(0) as usize;
    let mut groups = vec![QuadGroup::default(); count];
    for row in rows {
        let group = (row.window_group as usize)
            .checked_sub(1)
            .ok_or("DiaFrameMsMsWindows has window group 0")?;
        let group = &mut groups[group];
        group.scan_starts.push(row.scan_num_begin as usize);
        group.scan_ends.push(row.scan_num_end as usize);
        group.windows.push(IsolationWindow::from_center(
            row.isolation_mz,
            row.isolation_width,
            row.collision_energy,
        ));
    }
    for group in &mut groups {
        let mut order: Vec<usize> = (0..group.scan_starts.len()).collect();
        order.sort_by_key(|&index| group.scan_starts[index]);
        *group = QuadGroup {
            scan_starts: order.iter().map(|&i| group.scan_starts[i]).collect(),
            scan_ends: order.iter().map(|&i| group.scan_ends[i]).collect(),
            windows: order.iter().map(|&i| group.windows[i]).collect(),
        };
    }
    Ok(groups)
}

/// How a scan range is split into spectra (timsrust's `QuadWindowExpansionStrategy`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Expansion {
    None,
    /// `n` overlapping halves: pieces of two `(end - start) / (n + 1)` steps.
    Even(usize),
    /// Pieces `span` wide in 1/K0, `step` apart, on the uncalibrated linear scale.
    UniformMobility {
        span: f64,
        step: f64,
    },
    /// Pieces `span` scans wide, `step` scans apart.
    UniformScan {
        span: usize,
        step: usize,
    },
}

/// One diaPASEF spectrum: scans `scan_start..scan_end` of `frame`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DiaSlice {
    pub frame: u64,
    pub scan_start: usize,
    pub scan_end: usize,
    pub window: IsolationWindow,
}

/// Split `start..end` as timsrust's `scan_range_subsplit` does.
pub fn subsplit(
    start: usize,
    end: usize,
    expansion: Expansion,
    mobility: &LinearMobilityScale,
) -> Vec<(usize, usize)> {
    match expansion {
        Expansion::None => vec![(start, end)],
        Expansion::Even(splits) => {
            let width = end.saturating_sub(start) / (splits + 1);
            (0..splits)
                .map(|split| (start + width * split, start + width * (split + 2)))
                .collect()
        }
        Expansion::UniformMobility { span, step } => {
            // timsrust truncates the inverse conversion to an unsigned scan.
            let scan = |im: f64| mobility.scan_number(im) as u32 as usize;
            let mut start_offset = start;
            let mut start_im = mobility.one_over_k0(start as u32 as f64);
            let mut end_offset = scan(start_im - span);
            let mut out = Vec::new();
            // A step that does not move would never end; timsrust hangs here.
            if step > 0.0 {
                while end_offset < end {
                    out.push((start_offset, end_offset));
                    start_im -= step;
                    start_offset = scan(start_im);
                    end_offset = scan(start_im - span);
                }
            }
            if start_offset < end {
                out.push((start_offset, end));
            }
            out
        }
        Expansion::UniformScan { span, step } => {
            let mut start_offset = start;
            let mut end_offset = start + span;
            let mut out = Vec::new();
            if step > 0 {
                while end_offset < end {
                    out.push((start_offset, end_offset));
                    start_offset += step;
                    end_offset += step;
                }
            }
            if start_offset < end {
                out.push((start_offset, end));
            }
            out
        }
    }
}

/// Expand every diaPASEF frame into its spectra, in `DiaFrameMsMsInfo` order.
///
/// With `per_window`, each isolation window of the frame's group is split on
/// its own (timsrust's `Quadrupole` splitting). Otherwise the group's whole
/// scan range is split and each piece covers the windows it overlaps
/// (`Window` splitting).
pub fn dia_slices(
    frames: &[DiaFrameMsMsInfo],
    groups: &[QuadGroup],
    per_window: bool,
    expansion: Expansion,
    mobility: &LinearMobilityScale,
) -> Result<Vec<DiaSlice>, String> {
    let mut slices = Vec::new();
    for info in frames {
        let group = (info.window_group as usize)
            .checked_sub(1)
            .and_then(|index| groups.get(index))
            .ok_or_else(|| {
                format!(
                    "frame {} uses window group {}, which has no DiaFrameMsMsWindows rows",
                    info.frame, info.window_group
                )
            })?;
        if per_window {
            for (k, window) in group.windows.iter().enumerate() {
                for (scan_start, scan_end) in subsplit(
                    group.scan_starts[k],
                    group.scan_ends[k],
                    expansion,
                    mobility,
                ) {
                    slices.push(DiaSlice {
                        frame: info.frame,
                        scan_start,
                        scan_end,
                        window: *window,
                    });
                }
            }
            continue;
        }
        let (Some(&group_start), Some(&group_end)) =
            (group.scan_starts.iter().min(), group.scan_ends.iter().max())
        else {
            return Err(format!("window group {} has no windows", info.window_group));
        };
        for (sws, swe) in subsplit(group_start, group_end, expansion, mobility) {
            let mut mz_min = f64::MAX;
            let mut mz_max = f64::MIN;
            let mut nce_sum = 0.0;
            let mut total_scan_width = 0.0;
            for (k, window) in group.windows.iter().enumerate() {
                let (gss, gse) = (group.scan_starts[k], group.scan_ends[k]);
                if swe <= gse || gss <= sws {
                    continue;
                }
                let half_isolation_width = window.width() / 2.0;
                let isolation_mz = window.center();
                mz_min = mz_min.min(isolation_mz - half_isolation_width);
                mz_max = mz_max.max(isolation_mz + half_isolation_width);
                let scan_width = gse.min(swe).saturating_sub(gss.max(sws)) as f64;
                nce_sum += window.collision_energy * scan_width;
                total_scan_width += scan_width;
            }
            slices.push(DiaSlice {
                frame: info.frame,
                scan_start: sws,
                scan_end: swe,
                window: IsolationWindow::from_bounds(mz_min, mz_max, nce_sum / total_scan_width),
            });
        }
    }
    Ok(slices)
}

#[cfg(test)]
#[path = "../tests/unit/tdf_spectra.rs"]
mod tests;
