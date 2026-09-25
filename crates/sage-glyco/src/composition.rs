//! Glycan compositions, the composition library, and oxonium and core Y-ion
//! evidence:
//!
//! * parsing glycan compositions (`HexNAc(4)Hex(5)Fuc(1)NeuAc(2)`, and the
//!   pGlyco-style one-letter form `N(4)H(5)F(1)A(2)`),
//! * a mass-sorted composition library that answers "which compositions
//!   explain this precursor delta", including isotope-error shifts,
//! * oxonium-ion evidence for gating spectra as glyco or non-glyco,
//! * the N-glycan core Y-ion ladder (peptide plus partial glycan) and its
//!   evidence in a spectrum.
//!
//! Masses are monoisotopic residue masses (the monosaccharide minus water), so
//! a composition's mass is exactly the delta it adds to the peptide.

use std::fmt::{Display, Write as _};

use sage_core::mass::{Tolerance, NEUTRON, PROTON};
use sage_core::spectrum::ProcessedSpectrum;

/// Monosaccharide residues searched in N- and mucin-type O-glycans.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Monosaccharide {
    HexNAc,
    Hex,
    /// Deoxyhexose; in mammalian glycans this is fucose.
    Fuc,
    NeuAc,
    NeuGc,
}

impl Monosaccharide {
    pub const ALL: [Monosaccharide; 5] = [
        Monosaccharide::HexNAc,
        Monosaccharide::Hex,
        Monosaccharide::Fuc,
        Monosaccharide::NeuAc,
        Monosaccharide::NeuGc,
    ];

    /// Monoisotopic residue mass (C8H13NO5 for HexNAc, and so on).
    pub const fn mass(self) -> f64 {
        match self {
            Monosaccharide::HexNAc => 203.079_373,
            Monosaccharide::Hex => 162.052_824,
            Monosaccharide::Fuc => 146.057_909,
            Monosaccharide::NeuAc => 291.095_417,
            Monosaccharide::NeuGc => 307.090_331,
        }
    }

    pub const fn name(self) -> &'static str {
        match self {
            Monosaccharide::HexNAc => "HexNAc",
            Monosaccharide::Hex => "Hex",
            Monosaccharide::Fuc => "Fuc",
            Monosaccharide::NeuAc => "NeuAc",
            Monosaccharide::NeuGc => "NeuGc",
        }
    }

    fn parse(token: &str) -> Option<Self> {
        Some(match token {
            "HexNAc" | "N" => Monosaccharide::HexNAc,
            "Hex" | "H" => Monosaccharide::Hex,
            "Fuc" | "dHex" | "F" => Monosaccharide::Fuc,
            "NeuAc" | "A" => Monosaccharide::NeuAc,
            "NeuGc" | "G" => Monosaccharide::NeuGc,
            _ => return None,
        })
    }
}

/// Counts of each monosaccharide, indexed like [`Monosaccharide::ALL`].
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct GlycanComposition(pub [u8; 5]);

impl GlycanComposition {
    pub fn count(&self, residue: Monosaccharide) -> u8 {
        self.0[residue as usize]
    }

    pub fn mass(&self) -> f64 {
        Monosaccharide::ALL
            .iter()
            .map(|residue| residue.mass() * self.count(*residue) as f64)
            .sum()
    }

    /// Whether `other` could be a fragment of this composition.
    pub fn contains(&self, other: &GlycanComposition) -> bool {
        self.0.iter().zip(other.0.iter()).all(|(a, b)| a >= b)
    }

    /// Parse `HexNAc(4)Hex(5)Fuc(1)NeuAc(2)` or `N(4)H(5)F(1)A(2)`. Residues
    /// may appear in any order but only once; a count of zero is allowed.
    pub fn parse(text: &str) -> Result<Self, String> {
        let mut counts = [0u8; 5];
        let mut seen = [false; 5];
        let mut rest = text.trim();
        if rest.is_empty() {
            return Err("empty glycan composition".into());
        }
        while !rest.is_empty() {
            let open = rest
                .find('(')
                .ok_or_else(|| format!("expected `(` in glycan composition `{text}`"))?;
            let close = rest
                .find(')')
                .ok_or_else(|| format!("expected `)` in glycan composition `{text}`"))?;
            if close < open {
                return Err(format!("malformed glycan composition `{text}`"));
            }
            let name = &rest[..open];
            let residue = Monosaccharide::parse(name)
                .ok_or_else(|| format!("unknown monosaccharide `{name}` in `{text}`"))?;
            let count: u8 = rest[open + 1..close]
                .parse()
                .map_err(|_| format!("invalid count for `{name}` in `{text}`"))?;
            if std::mem::replace(&mut seen[residue as usize], true) {
                return Err(format!("`{name}` repeated in glycan composition `{text}`"));
            }
            counts[residue as usize] = count;
            rest = rest[close + 1..].trim_start();
        }
        Ok(Self(counts))
    }
}

impl Display for GlycanComposition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut out = String::new();
        for residue in Monosaccharide::ALL {
            let count = self.count(residue);
            if count > 0 {
                write!(out, "{}({count})", residue.name())?;
            }
        }
        f.write_str(&out)
    }
}

/// A composition library sorted by mass, for precursor delta lookup.
#[derive(Clone, Debug, Default)]
pub struct GlycanLibrary {
    compositions: Vec<(f64, GlycanComposition)>,
}

impl GlycanLibrary {
    pub fn new(compositions: impl IntoIterator<Item = GlycanComposition>) -> Self {
        let mut compositions: Vec<_> = compositions
            .into_iter()
            .map(|composition| (composition.mass(), composition))
            .collect();
        compositions.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
        compositions.dedup_by(|a, b| a.1 == b.1);
        Self { compositions }
    }

    /// One composition per non-empty, non-comment line.
    pub fn parse(text: &str) -> Result<Self, String> {
        let compositions = text
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .map(GlycanComposition::parse)
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self::new(compositions))
    }

    pub fn len(&self) -> usize {
        self.compositions.len()
    }

    pub fn is_empty(&self) -> bool {
        self.compositions.is_empty()
    }

    pub fn iter(&self) -> impl Iterator<Item = &(f64, GlycanComposition)> {
        self.compositions.iter()
    }

    /// Lightest and heaviest composition masses: the precursor window an
    /// open glyco search must cover below the observed mass.
    pub fn mass_range(&self) -> Option<(f64, f64)> {
        Some((self.compositions.first()?.0, self.compositions.last()?.0))
    }

    /// Compositions whose mass explains `delta` (observed precursor mass minus
    /// peptide mass) within `tolerance_da`, for each allowed isotope error.
    /// Returns `(library index, isotope error)`.
    pub fn explain(
        &self,
        delta: f64,
        tolerance_da: f64,
        isotope_errors: std::ops::RangeInclusive<i8>,
    ) -> Vec<(usize, i8)> {
        let mut hits = Vec::new();
        for isotope in isotope_errors {
            let target = delta - isotope as f64 * NEUTRON as f64;
            let lo = self
                .compositions
                .partition_point(|(mass, _)| *mass < target - tolerance_da);
            for (index, (mass, _)) in self.compositions.iter().enumerate().skip(lo) {
                if *mass > target + tolerance_da {
                    break;
                }
                hits.push((index, isotope));
            }
        }
        hits
    }

    pub fn get(&self, index: usize) -> Option<&GlycanComposition> {
        self.compositions
            .get(index)
            .map(|(_, composition)| composition)
    }

    /// Pairs of distinct compositions a search cannot separate by precursor
    /// mass alone: masses within `tolerance_da` after any isotope shift in
    /// `0..=max_isotope_shift` neutrons. Returns `(i, j, shift)` with `i < j`.
    pub fn indistinguishable_pairs(
        &self,
        tolerance_da: f64,
        max_isotope_shift: u8,
    ) -> Vec<(usize, usize, u8)> {
        let mut pairs = Vec::new();
        for (i, (mass_i, _)) in self.compositions.iter().enumerate() {
            for shift in 0..=max_isotope_shift {
                for (index, _) in
                    self.explain(mass_i + shift as f64 * NEUTRON as f64, tolerance_da, 0..=0)
                {
                    if index != i {
                        pairs.push((i.min(index), i.max(index), shift));
                    }
                }
            }
        }
        pairs.sort_unstable();
        pairs.dedup();
        pairs
    }
}

/// A singly charged glycan fragment used as glyco evidence.
#[derive(Copy, Clone, Debug)]
pub struct OxoniumIon {
    pub name: &'static str,
    pub mz: f32,
}

/// Oxonium ions checked for gating, most diagnostic first. Values are [M+H]+.
pub const OXONIUM_IONS: [OxoniumIon; 10] = [
    OxoniumIon {
        name: "HexNAc",
        mz: 204.086_6,
    },
    OxoniumIon {
        name: "HexNAc-C2H6O3",
        mz: 138.055,
    },
    OxoniumIon {
        name: "HexNAc-2H2O",
        mz: 168.065_5,
    },
    OxoniumIon {
        name: "HexNAc-H2O",
        mz: 186.076_1,
    },
    OxoniumIon {
        name: "HexNAc-C2H4O2",
        mz: 144.065_5,
    },
    OxoniumIon {
        name: "HexHexNAc",
        mz: 366.139_5,
    },
    OxoniumIon {
        name: "NeuAc",
        mz: 292.102_7,
    },
    OxoniumIon {
        name: "NeuAc-H2O",
        mz: 274.092_1,
    },
    OxoniumIon {
        name: "NeuGc",
        mz: 308.097_6,
    },
    OxoniumIon {
        name: "Hex",
        mz: 163.060_1,
    },
];

/// Oxonium evidence for one spectrum.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OxoniumEvidence {
    /// Which of [`OXONIUM_IONS`] were observed.
    pub matched: Vec<bool>,
    /// Summed oxonium intensity over total ion current.
    pub intensity_fraction: f32,
}

impl OxoniumEvidence {
    pub fn count(&self) -> usize {
        self.matched.iter().filter(|m| **m).count()
    }

    /// The gate used by MSFragger-Glyco and pGlyco-style filters: the HexNAc
    /// oxonium ion plus at least `min_other` other oxonium ions.
    pub fn is_glyco(&self, min_other: usize) -> bool {
        self.matched.first().copied().unwrap_or(false) && self.count() > min_other
    }
}

/// Most intense peak whose singly charged m/z matches `mz`. Oxonium ions
/// are 1+, so peaks the deisotoper assigned a higher charge are skipped.
pub(crate) fn find_singly_charged(
    query: &ProcessedSpectrum,
    mz: f32,
    tolerance: Tolerance,
) -> Option<usize> {
    let neutral = mz - PROTON;
    let (lo, hi) = tolerance.bounds(neutral);
    let start = query.masses.partition_point(|mass| *mass < lo);
    (start..query.masses.len())
        .take_while(|&idx| query.masses[idx] <= hi)
        .filter(|&idx| query.charges.get(idx).copied().unwrap_or(1) <= 1)
        .max_by(|&a, &b| query.intensities[a].total_cmp(&query.intensities[b]))
}

pub fn oxonium_evidence(query: &ProcessedSpectrum, tolerance: Tolerance) -> OxoniumEvidence {
    let mut summed = 0.0;
    let matched = OXONIUM_IONS
        .iter()
        .map(|ion| match find_singly_charged(query, ion.mz, tolerance) {
            Some(idx) => {
                summed += query.intensities[idx];
                true
            }
            None => false,
        })
        .collect();
    OxoniumEvidence {
        matched,
        intensity_fraction: if query.total_ion_current > 0.0 {
            summed / query.total_ion_current
        } else {
            0.0
        },
    }
}

/// Cross-ring 0,2X fragment of the innermost HexNAc (C4H5NO), which stays on
/// the peptide in HCD of N-glycopeptides.
pub const HEXNAC_CROSS_RING: f64 = 83.037_114;

/// One Y ion: the peptide carrying a partial glycan.
#[derive(Clone, Debug, PartialEq)]
pub struct YIon {
    pub label: String,
    /// Mass added to the bare peptide.
    pub delta: f64,
}

/// Core N-glycan Y ions consistent with `composition`: Y0 (bare peptide),
/// the 0,2X cross-ring ion, and the chitobiose/trimannosyl core ladder with
/// and without core fucose. Only fragments contained in the composition are
/// generated, so the ladder doubles as composition evidence.
pub fn n_glycan_core_y_ions(composition: &GlycanComposition) -> Vec<YIon> {
    const CORE: [(u8, u8); 6] = [(0, 0), (1, 0), (2, 0), (2, 1), (2, 2), (2, 3)];
    let mut ions = vec![
        YIon {
            label: "Y0".into(),
            delta: 0.0,
        },
        YIon {
            label: "Y0+0,2X".into(),
            delta: HEXNAC_CROSS_RING,
        },
    ];
    let fucosylated = composition.count(Monosaccharide::Fuc) > 0;
    for (hexnac, hex) in CORE.into_iter().skip(1) {
        for fuc in 0..=u8::from(fucosylated) {
            let partial = GlycanComposition([hexnac, hex, fuc, 0, 0]);
            if composition.contains(&partial) && partial != *composition {
                ions.push(YIon {
                    label: format!("Y[{partial}]"),
                    delta: partial.mass(),
                });
            }
        }
    }
    ions
}

/// Y-ion evidence for a candidate peptide mass and composition.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct YIonEvidence {
    pub generated: usize,
    pub matched: usize,
    /// Summed matched Y-ion intensity over total ion current.
    pub intensity_fraction: f32,
    /// Whether the bare peptide (Y0) or HexNAc-peptide (Y1) was observed:
    /// the ions that anchor the peptide mass independently of composition.
    pub anchored: bool,
}

pub fn y_ion_evidence(
    query: &ProcessedSpectrum,
    peptide_mass: f32,
    composition: &GlycanComposition,
    precursor_charge: u8,
    tolerance: Tolerance,
) -> YIonEvidence {
    let ions = n_glycan_core_y_ions(composition);
    let mut evidence = YIonEvidence {
        generated: ions.len(),
        ..Default::default()
    };
    let mut summed = 0.0;
    for ion in &ions {
        let neutral = peptide_mass + ion.delta as f32;
        let best = (1..=precursor_charge.max(1))
            .filter_map(|charge| {
                let mz = neutral / charge as f32 + PROTON;
                find_at_charge(query, mz, charge, tolerance)
            })
            .max_by(|&a, &b| query.intensities[a].total_cmp(&query.intensities[b]));
        if let Some(idx) = best {
            evidence.matched += 1;
            summed += query.intensities[idx];
            if ion.label == "Y0" || ion.label == "Y[HexNAc(1)]" {
                evidence.anchored = true;
            }
        }
    }
    if query.total_ion_current > 0.0 {
        evidence.intensity_fraction = summed / query.total_ion_current;
    }
    evidence
}

/// Peak at `mz` for `charge`: a deisotoped peak of that charge, or a peak of
/// unknown charge (stored as singly charged) at that m/z.
pub(crate) fn find_at_charge(
    query: &ProcessedSpectrum,
    mz: f32,
    charge: u8,
    tolerance: Tolerance,
) -> Option<usize> {
    let known = (mz - PROTON) * charge as f32;
    let (lo, hi) = tolerance.bounds(known);
    let start = query.masses.partition_point(|mass| *mass < lo);
    let deisotoped = (start..query.masses.len())
        .take_while(|&idx| query.masses[idx] <= hi)
        .filter(|&idx| query.has_known_charge(idx) && query.charges[idx] == charge);
    let unknown = mz - PROTON;
    let (lo, hi) = tolerance.bounds(unknown);
    let start = query.masses.partition_point(|mass| *mass < lo);
    let undetermined = (start..query.masses.len())
        .take_while(|&idx| query.masses[idx] <= hi)
        .filter(|&idx| !query.has_known_charge(idx));
    deisotoped
        .chain(undetermined)
        .max_by(|&a, &b| query.intensities[a].total_cmp(&query.intensities[b]))
}

/// A biosynthetically plausible mammalian N-glycan composition space, used
/// to size the ambiguity problem. It is a superset of curated lists such as
/// the 182 human N-glycans shipped with Byonic and MSFragger-Glyco, not a
/// substitute for them.
pub fn n_glycan_composition_space() -> Vec<GlycanComposition> {
    let mut out = Vec::new();
    for hexnac in 2..=7u8 {
        for hex in 3..=12u8 {
            // High-mannose and hybrid glycans carry few antennary HexNAc.
            if hex > 6 && hexnac > 3 {
                continue;
            }
            for fuc in 0..=3u8 {
                for neuac in 0..=4u8 {
                    for neugc in 0..=2u8 {
                        let antennae = hexnac - 2;
                        if neuac + neugc > antennae || fuc > antennae + 1 {
                            continue;
                        }
                        out.push(GlycanComposition([hexnac, hex, fuc, neuac, neugc]));
                    }
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spectrum(peaks: &[(f32, u8, bool, f32)]) -> ProcessedSpectrum {
        let mut peaks = peaks.to_vec();
        peaks.sort_by(|a, b| a.0.total_cmp(&b.0));
        ProcessedSpectrum {
            level: 2,
            masses: peaks.iter().map(|p| p.0).collect(),
            charges: peaks.iter().map(|p| p.1).collect(),
            charge_is_known: peaks.iter().map(|p| p.2).collect(),
            intensities: peaks.iter().map(|p| p.3).collect(),
            total_ion_current: peaks.iter().map(|p| p.3).sum(),
            ..Default::default()
        }
    }

    #[test]
    fn parses_both_notations() {
        let long = GlycanComposition::parse("HexNAc(4)Hex(5)Fuc(1)NeuAc(2)").unwrap();
        let short = GlycanComposition::parse("N(4)H(5)F(1)A(2)").unwrap();
        assert_eq!(long, short);
        assert_eq!(long.to_string(), "HexNAc(4)Hex(5)Fuc(1)NeuAc(2)");
        // A2G2FS2, the common core-fucosylated disialylated biantennary glycan.
        assert!((long.mass() - 2350.8304).abs() < 1e-3, "{}", long.mass());
        assert!(GlycanComposition::parse("HexNAc(2)HexNAc(1)").is_err());
        assert!(GlycanComposition::parse("Kdn(1)").is_err());
    }

    #[test]
    fn oxonium_ions_match_residue_masses() {
        let hexnac = Monosaccharide::HexNAc.mass() as f32 + PROTON;
        assert!((OXONIUM_IONS[0].mz - hexnac).abs() < 1e-3);
        let hexhexnac =
            (Monosaccharide::HexNAc.mass() + Monosaccharide::Hex.mass()) as f32 + PROTON;
        assert!((OXONIUM_IONS[5].mz - hexhexnac).abs() < 1e-3);
        let neuac = Monosaccharide::NeuAc.mass() as f32 + PROTON;
        assert!((OXONIUM_IONS[6].mz - neuac).abs() < 1e-3);
    }

    #[test]
    fn explains_delta_with_isotope_error() {
        let library = GlycanLibrary::parse(
            "# comment\nHexNAc(2)Hex(5)\nHexNAc(4)Hex(5)Fuc(1)\nHexNAc(4)Hex(5)NeuAc(1)\n",
        )
        .unwrap();
        assert_eq!(library.len(), 3);
        let target = GlycanComposition::parse("HexNAc(4)Hex(5)Fuc(1)").unwrap();
        let delta = target.mass() + NEUTRON as f64;
        let hits = library.explain(delta, 0.02, 0..=1);
        assert_eq!(hits.len(), 1);
        assert_eq!(library.get(hits[0].0), Some(&target));
        assert_eq!(hits[0].1, 1);
    }

    #[test]
    fn known_isomeric_and_near_isobaric_pairs() {
        // NeuAc + Hex and NeuGc + Fuc have the same elemental formula.
        let a = GlycanComposition::parse("HexNAc(4)Hex(5)NeuAc(1)").unwrap();
        let b = GlycanComposition::parse("HexNAc(4)Hex(4)Fuc(1)NeuGc(1)").unwrap();
        assert!((a.mass() - b.mass()).abs() < 1e-4);
        // Two Fuc sit 17 mDa from one NeuAc plus one neutron: only an
        // isotope-error assignment separates them.
        let c = GlycanComposition::parse("HexNAc(4)Hex(5)Fuc(2)").unwrap();
        let d = GlycanComposition::parse("HexNAc(4)Hex(5)NeuAc(1)").unwrap();
        let gap = c.mass() - (d.mass() + NEUTRON as f64);
        assert!(gap.abs() < 0.03, "{gap}");
    }

    #[test]
    fn oxonium_gate_and_y_ladder() {
        let composition = GlycanComposition::parse("HexNAc(4)Hex(5)Fuc(1)").unwrap();
        let peptide = 1500.7;
        let y1 = peptide + Monosaccharide::HexNAc.mass() as f32;
        let query = spectrum(&[
            (OXONIUM_IONS[0].mz - PROTON, 1, false, 100.0),
            (OXONIUM_IONS[1].mz - PROTON, 1, false, 40.0),
            (OXONIUM_IONS[5].mz - PROTON, 1, false, 20.0),
            // Y0 observed as a deisotoped 2+ peak, Y1 as an undetermined 2+ m/z.
            (peptide, 2, true, 30.0),
            (y1 / 2.0 + PROTON - PROTON, 1, false, 50.0),
            (500.0, 1, false, 10.0),
        ]);
        let tol = Tolerance::Ppm(-20.0, 20.0);
        let oxonium = oxonium_evidence(&query, tol);
        assert_eq!(oxonium.count(), 3);
        assert!(oxonium.is_glyco(1));

        let ladder = n_glycan_core_y_ions(&composition);
        // Y0, 0,2X, and ten core fragments (five core ions, with and without Fuc).
        assert_eq!(ladder.len(), 12, "{ladder:?}");
        let evidence = y_ion_evidence(&query, peptide, &composition, 3, tol);
        assert_eq!(evidence.matched, 2);
        assert!(evidence.anchored);

        let blank = spectrum(&[(500.0, 1, false, 10.0)]);
        assert!(!oxonium_evidence(&blank, tol).is_glyco(1));
    }

    /// Sizes the composition ambiguity problem for the design doc:
    /// `cargo test -p sage-glyco composition::tests::ambiguity_report -- --ignored --nocapture`
    #[test]
    #[ignore]
    fn ambiguity_report() {
        let space = n_glycan_composition_space();
        let human = space
            .iter()
            .copied()
            .filter(|c| c.count(Monosaccharide::NeuGc) == 0);
        for (label, library) in [
            ("human (no NeuGc)", GlycanLibrary::new(human)),
            ("mammalian (with NeuGc)", GlycanLibrary::new(space.clone())),
        ] {
            let (lo, hi) = library.mass_range().unwrap();
            println!(
                "{label}: {} compositions, mass range {lo:.1}..{hi:.1} Da",
                library.len()
            );
            for (tolerance, isotopes) in [(0.001, 0), (0.02, 0), (0.02, 1), (0.02, 2), (0.05, 2)] {
                let pairs = library.indistinguishable_pairs(tolerance, isotopes);
                let mut involved: Vec<usize> = pairs.iter().flat_map(|p| [p.0, p.1]).collect();
                involved.sort_unstable();
                involved.dedup();
                println!(
                "tolerance {tolerance} Da, isotope shifts 0..={isotopes}: {} pairs, {} of {} compositions ambiguous",
                pairs.len(),
                involved.len(),
                library.len()
            );
            }
            for (i, j, shift) in library.indistinguishable_pairs(0.02, 1).into_iter().take(6) {
                let (a, b) = (library.get(i).unwrap(), library.get(j).unwrap());
                println!("  {a} ~ {b} (+{shift} n, {:.4} Da)", b.mass() - a.mass());
            }
        }
    }
}
