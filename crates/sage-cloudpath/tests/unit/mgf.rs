use sage_core::{
    mass::Tolerance,
    spectrum::{RawSpectrum, Representation},
};

use super::{MgfError, MgfReader};

fn make_ions_section_spectrum_0() -> String {
    let s = r#"
        BEGIN IONS
        TITLE=spectrum 0, 1/K0=0.966
        RTINSECONDS=0.8963232289
        PEPMASS=367.069682741984 56700.5185546875
        CHARGE=2+ and 3+
        TOL=10
        TOLU=ppm
        148.2041016 
        169.5001831 4608.2421875
        226.0483246 5335.4907226563
        228.3407898 30918.244140625
        322.5945435 5311.5737304688
        1144.66272 6260.8315429688
        END IONS
        "#;
    return String::from(s);
}

fn run_asserts_for_spectrum_0(s: &RawSpectrum) {
    assert_eq!(s.id, "spectrum 0, 1/K0=0.966");
    assert_eq!(s.ms_level, 2);
    assert_eq!(s.representation, Representation::Centroid);
    assert_eq!(s.precursors.len(), 2);
    assert_eq!(s.precursors[0].charge, Some(2));
    assert_eq!(s.precursors[1].charge, Some(3));
    assert_eq!(s.precursors[0].inverse_ion_mobility, Some(0.966));
    assert_eq!(s.precursors[1].inverse_ion_mobility, Some(0.966));
    assert!((s.precursors[0].mz - 367.069682741984).abs() < 0.0001);
    assert_eq!(s.precursors[0].intensity, Some(56700.5185546875));
    assert_eq!(
        s.precursors[0].isolation_window,
        Some(Tolerance::Ppm(-10.0, 10.0))
    );
    assert!((s.precursors[1].mz - 367.069682741984).abs() < 0.0001);
    assert_eq!(s.precursors[1].intensity, Some(56700.5185546875));
    assert_eq!(
        s.precursors[1].isolation_window,
        Some(Tolerance::Ppm(-10.0, 10.0))
    );
    assert!((s.scan_start_time - 0.8963232289 / 60.0).abs() < 0.0001);
    assert_eq!(s.ion_injection_time, 0.0);
    assert_eq!(s.intensity.len(), s.mz.len());
    assert!((s.mz[3] - 228.3407898).abs() < 0.0001);
    assert!((s.intensity[0] - 1.0).abs() < 0.0001);
}

#[tokio::test]
async fn parse_spectrum() -> Result<(), MgfError> {
    let s = make_ions_section_spectrum_0();
    let mut spectra = MgfReader::with_file_id(0).parse(s)?;

    assert_eq!(spectra.len(), 1);
    let s = spectra.pop().unwrap();

    run_asserts_for_spectrum_0(&s);
    Ok(())
}

#[test]
fn empty_input_returns_an_error() {
    assert!(matches!(
        MgfReader::with_file_id(0).parse(String::new()),
        Err(MgfError::MissingBeginIons)
    ));
}

#[test]
fn missing_begin_marker_returns_an_error() {
    assert!(matches!(
        MgfReader::with_file_id(0).parse("TITLE=orphan\n100.0 20.0\n".into()),
        Err(MgfError::MissingBeginIons)
    ));
}

#[test]
fn malformed_spectrum_reports_its_line() {
    let error = MgfReader::with_file_id(0)
        .parse("BEGIN IONS\nTITLE=test\nPEPMASS=invalid\nEND IONS\n".into())
        .unwrap_err();
    assert!(matches!(
        error,
        MgfError::MalformedLine {
            line: 3,
            message: _
        }
    ));
}

#[test]
fn missing_end_marker_returns_an_error() {
    assert!(matches!(
        MgfReader::with_file_id(0).parse("BEGIN IONS\nTITLE=test\nPEPMASS=500\n100 20\n".into()),
        Err(MgfError::UnterminatedSpectrum)
    ));
}

#[test]
fn empty_ions_block_is_unterminated() {
    assert!(matches!(
        MgfReader::with_file_id(0).parse("BEGIN IONS\n".into()),
        Err(MgfError::UnterminatedSpectrum)
    ));
}

#[tokio::test]
async fn parse_two_spectra() -> Result<(), MgfError> {
    let mut content = "# a comment at the beginning of the file".to_string();
    content.push_str(&make_ions_section_spectrum_0());
    content.push_str("\n\n");
    content.push_str(&make_ions_section_spectrum_0());

    let spectra = MgfReader::with_file_id(0).parse(content)?;
    assert_eq!(spectra.len(), 2);
    spectra
        .iter()
        .for_each(|spec: &RawSpectrum| run_asserts_for_spectrum_0(spec));
    Ok(())
}

#[tokio::test]
/// Example taken from https://www.matrixscience.com/help/data_file_help.html
async fn parse_mgf_matrixscience_example_1() -> Result<(), MgfError> {
    let s = r#"
        COM=10 pmol digest of Sample X15
        ITOL=1
        ITOLU=Da
        MODS=Carbamidomethyl (C)
        IT_MODS=Oxidation (M)
        MASS=Monoisotopic
        USERNAME=Lou Scene
        USEREMAIL=leu@altered-state.edu
        CHARGE=2+ and 3+
        BEGIN IONS
        TITLE=Spectrum 1
        PEPMASS=983.6
        846.60 73
        846.80 44
        847.60 67
        1640.10 291
        1640.60 54
        1895.50 49
        END IONS

        BEGIN IONS
        TITLE=Spectrum 2
        PEPMASS=1084.9
        SCANS=3
        RTINSECONDS=25
        345.10 237
        370.20 128
        460.20 108
        1673.30 1007
        1674.00 974
        1675.30 79
        END IONS
        "#;
    let mut spectra = MgfReader::with_file_id(0).parse(s.to_string())?;
    assert_eq!(spectra.len(), 2);

    let s = spectra.pop().unwrap();
    assert_eq!(s.precursors.len(), 2);
    assert_eq!(s.precursors[0].charge, Some(2));
    assert_eq!(s.precursors[1].charge, Some(3));
    assert_eq!(s.precursors[0].isolation_window, None);
    Ok(())
}

#[tokio::test]
/// Example taken from https://www.matrixscience.com/help/data_file_help.html
async fn parse_mgf_matrixscience_example_2() -> Result<(), MgfError> {
    let s = r#"
        # following lines define parameters.
        # NB no spaces allowed on either side of the = symbol
        COM=My favourite protein has been eaten by an enzyme
        CLE=Trypsin
        CHARGE=2+
        # following line will be treated as a peptide mass
        1024.6
        # following line is a sequence query, which must
        # conform precisely to sequence query syntax rules
        2321 seq(n-ACTL) comp(2[C])
        # so is this
        1896 ions(345.6:24.7,347.8:45.4, ... ,1024.7:18.7)
        # An MS/MS ions query is delimited by the tags
        # BEGIN IONS and END IONS. Space(s)
        # are used to separate mass and intensity values
        BEGIN IONS
        TITLE=The first peptide - dodgy peak detection, so extra wide tolerance
        PEPMASS=896.05 25674.3
        CHARGE=3+
        TOL=3
        TOLU=Da
        SEQ=n-AC[DHK]
        COMP=2[H]0[M]3[DE]*[K]
        240.1 3
        242.1 12
        245.2 32
        1623.7 55
        1624.7 23
        END IONS
        "#;
    let mut spectra = MgfReader::with_file_id(0).parse(s.to_string())?;
    assert_eq!(spectra.len(), 1);

    let s = spectra.pop().unwrap();
    assert_eq!(s.precursors.len(), 1);
    assert_eq!(s.precursors[0].charge, Some(3));
    assert_eq!(
        s.precursors[0].isolation_window,
        Some(Tolerance::Da(-3.0, 3.0))
    );
    Ok(())
}
