//! Throwaway RSS instrumentation for the 2026-11 reader-memory audit.
//!
//! NOT FOR LANDING. It prints a `MEMMARK` line at stage boundaries with the
//! process's resident set size, and, for one file's spectra, how many bytes
//! the peak arrays hold against how many bytes they have reserved, so the
//! difference between the Thermo and mzML read paths can be attributed.
use sage_core::spectrum::{ProcessedSpectrum, RawSpectrum};

/// Resident set size of this process, in GiB.
pub fn rss_gib() -> f64 {
    let statm = std::fs::read_to_string("/proc/self/statm").unwrap_or_default();
    let pages: f64 = statm
        .split_whitespace()
        .nth(1)
        .and_then(|field| field.parse().ok())
        .unwrap_or(0.0);
    pages * 4096.0 / (1024.0 * 1024.0 * 1024.0)
}

/// One stage boundary.
pub fn mark(label: &str) {
    log::info!("MEMMARK stage={label} rss={:.3}", rss_gib());
}

/// What one file's raw spectra hold, as read, before processing.
pub fn raw_stats(file_id: usize, spectra: &[RawSpectrum]) {
    let mut peaks = 0usize;
    let mut mz_len = 0usize;
    let mut mz_cap = 0usize;
    let mut int_cap = 0usize;
    let mut id_bytes = 0usize;
    let mut ms1 = 0usize;
    for spectrum in spectra {
        peaks += spectrum.mz.len();
        mz_len += spectrum.mz.len();
        mz_cap += spectrum.mz.capacity();
        int_cap += spectrum.intensity.capacity();
        id_bytes += spectrum.id.capacity();
        if spectrum.ms_level == 1 {
            ms1 += 1;
        }
    }
    let used = (mz_len + mz_len) * 4;
    let reserved = (mz_cap + int_cap) * 4;
    log::info!(
        "MEMMARK raw file={file_id} spectra={} ms1={ms1} peaks={peaks} \
         peak_bytes_used={used} peak_bytes_reserved={reserved} \
         slack_bytes={} id_bytes={id_bytes} struct_bytes={} rss={:.3}",
        spectra.len(),
        reserved.saturating_sub(used),
        std::mem::size_of_val(spectra),
        rss_gib(),
    );
}

/// What one file's processed spectra hold, after `SpectrumProcessor`.
pub fn processed_stats(file_id: usize, ms1: &[ProcessedSpectrum], msn: &[ProcessedSpectrum]) {
    let bytes = |spectra: &[ProcessedSpectrum]| -> (usize, usize, usize) {
        let mut peaks = 0usize;
        let mut used = 0usize;
        let mut reserved = 0usize;
        for spectrum in spectra {
            peaks += spectrum.masses.len();
            // masses + intensities (4 each), charges (1), charge_is_known (1)
            used += spectrum.masses.len() * 10;
            reserved += spectrum.masses.capacity() * 4
                + spectrum.intensities.capacity() * 4
                + spectrum.charges.capacity()
                + spectrum.charge_is_known.capacity();
        }
        (peaks, used, reserved)
    };
    let (ms1_peaks, ms1_used, ms1_reserved) = bytes(ms1);
    let (msn_peaks, msn_used, msn_reserved) = bytes(msn);
    log::info!(
        "MEMMARK processed file={file_id} ms1={} ms1_peaks={ms1_peaks} \
         ms1_bytes_used={ms1_used} ms1_bytes_reserved={ms1_reserved} \
         msn={} msn_peaks={msn_peaks} msn_bytes_used={msn_used} \
         msn_bytes_reserved={msn_reserved} rss={:.3}",
        ms1.len(),
        msn.len(),
        rss_gib(),
    );
}
