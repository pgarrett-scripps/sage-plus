use sage_core::{mass::Tolerance, spectrum::Representation};

use super::{MzMLError, MzMLReader};

#[tokio::test]
async fn reads_fragment_charge_binary_array() -> Result<(), MzMLError> {
    let input = r#"
    <mzML><run><spectrumList count="1">
      <spectrum id="scan=charges">
        <cvParam accession="MS:1000511" value="2"/>
        <cvParam accession="MS:1000127"/>
        <cvParam accession="MS:1000285" value="60"/>
        <binaryDataArrayList count="3">
          <binaryDataArray>
            <cvParam accession="MS:1000514"/>
            <cvParam accession="MS:1000523"/>
            <cvParam accession="MS:1000576"/>
            <binary>AAAAAAAAeUAAAAAAAEB/QAAAAAAAwIJA</binary>
          </binaryDataArray>
          <binaryDataArray>
            <cvParam accession="MS:1000515"/>
            <cvParam accession="MS:1000521"/>
            <cvParam accession="MS:1000576"/>
            <binary>AAAgQQAAoEEAAPBB</binary>
          </binaryDataArray>
          <binaryDataArray>
            <cvParam accession="MS:1000516"/>
            <cvParam accession="MS:1000519"/>
            <cvParam accession="MS:1000576"/>
            <binary>AgAAAAAAAAADAAAA</binary>
          </binaryDataArray>
        </binaryDataArrayList>
      </spectrum>
    </spectrumList></run></mzML>
    "#;

    let spectra = MzMLReader::with_file_id(2).parse(input.as_bytes()).await?;

    assert_eq!(spectra.len(), 1);
    assert_eq!(spectra[0].mz, [400.0, 500.0, 600.0]);
    assert_eq!(spectra[0].intensity, [10.0, 20.0, 30.0]);
    assert_eq!(spectra[0].fragment_charges, Some(vec![2, 0, 3]));
    Ok(())
}

#[tokio::test]
async fn resolves_referenceable_param_groups_in_all_supported_scopes() -> Result<(), MzMLError> {
    let input = r#"
    <mzML>
      <referenceableParamGroupList count="5">
        <referenceableParamGroup id="spectrum-params">
          <cvParam accession="MS:1000511" value="2"/>
          <cvParam accession="MS:1000127"/>
          <cvParam accession="MS:1000285" value="100"/>
        </referenceableParamGroup>
        <referenceableParamGroup id="scan-params">
          <cvParam accession="MS:1000016" value="120" unitAccession="UO:0000010"/>
          <cvParam accession="MS:1000927" value="8.5"/>
        </referenceableParamGroup>
        <referenceableParamGroup id="precursor-params">
          <cvParam accession="MS:1000827" value="500.2"/>
          <cvParam accession="MS:1000828" value="1.0"/>
          <cvParam accession="MS:1000829" value="1.5"/>
        </referenceableParamGroup>
        <referenceableParamGroup id="ion-params">
          <cvParam accession="MS:1000744" value="500.1"/>
          <cvParam accession="MS:1000041" value="3"/>
          <cvParam accession="MS:1002815" value="1.2"/>
        </referenceableParamGroup>
        <referenceableParamGroup id="mz-array-params">
          <cvParam accession="MS:1000514"/>
          <cvParam accession="MS:1000523"/>
          <cvParam accession="MS:1000576"/>
        </referenceableParamGroup>
      </referenceableParamGroupList>
      <run><spectrumList count="1">
        <spectrum id="scan=1">
          <referenceableParamGroupRef ref="spectrum-params"/>
          <scanList count="1"><scan>
            <referenceableParamGroupRef ref="scan-params"/>
          </scan></scanList>
          <precursorList count="1"><precursor>
            <referenceableParamGroupRef ref="precursor-params"/>
            <selectedIonList count="1"><selectedIon>
              <referenceableParamGroupRef ref="ion-params"/>
            </selectedIon></selectedIonList>
          </precursor></precursorList>
          <binaryDataArrayList count="1"><binaryDataArray>
            <referenceableParamGroupRef ref="mz-array-params"/>
            <binary>AAAAAAAAWUAAAAAAAEBZQA==</binary>
          </binaryDataArray></binaryDataArrayList>
        </spectrum>
      </spectrumList></run>
    </mzML>
    "#;

    let spectra = MzMLReader::with_file_id(4).parse(input.as_bytes()).await?;
    assert_eq!(spectra.len(), 1);
    let spectrum = &spectra[0];
    assert_eq!(spectrum.file_id, 4);
    assert_eq!(spectrum.ms_level, 2);
    assert_eq!(
        spectrum.representation,
        sage_core::spectrum::Representation::Centroid
    );
    assert_eq!(spectrum.scan_start_time, 2.0);
    assert_eq!(spectrum.ion_injection_time, 8.5);
    assert_eq!(spectrum.mz, vec![100.0, 101.0]);
    assert_eq!(spectrum.precursors.len(), 1);
    assert_eq!(spectrum.precursors[0].mz, 500.1);
    assert_eq!(spectrum.precursors[0].charge, Some(3));
    assert_eq!(spectrum.precursors[0].inverse_ion_mobility, Some(1.2));
    assert_eq!(
        spectrum.precursors[0].isolation_window,
        Some(sage_core::mass::Tolerance::Da(-1.0, 1.5))
    );
    Ok(())
}

#[tokio::test]
async fn parse_spectrum_issue_78() -> Result<(), MzMLError> {
    let s = r#"
        <spectrum id="spectrum=2442" index="286" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum" />
            <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" />
            <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="638.352905273437955"/>
            <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="113.885513305664006"/>
            <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="793.395202636718977"/>
            <cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="147.290603637695"/>
            <cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="769.255798339843977"/>
            <userParam name="filter string" type="xsd:string" value="ITMS + c NSI d w Full ms2 457.72@cid35.00 [115.00-930.00]"/>
            <userParam name="preset scan configuration" type="xsd:string" value="2"/>
            <scanList count="1">
                <cvParam cvRef="MS" accession="MS:1000795" name="no combination" />
                <scan >
                    <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="1503.96166992188" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
                    <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:" type="xsd:double" value="457.723968505858977"/>
                    <scanWindowList count="1">
                        <scanWindow>
                            <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="115" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="930" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        </scanWindow>
                    </scanWindowList>
                </scan>
            </scanList>
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                            <cvParam cvRef="MS" accession="MS:1002815" name="inverse reduced ion mobility" value="1.078628" unitAccession="MS:1002814" unitName="volt-second per square centimeter"/>
                        </selectedIon>
                    </selectedIonList>
                    <activation>
                        <cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" />
                        <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="35.0"/>
                    </activation>
                </precursor>
            </precursorList>
            <binaryDataArrayList count="2">
                <binaryDataArray encodedLength="1088">
                    <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>AAAAoExpYkAAAACA3MpkQAAAAACph2VAAAAAAE4wZkAAAACAlMdmQAAAAECZAmdAAAAAwP9jaEAAAADgj4ZoQAAAAGC7HWlAAAAAAOXFaUAAAADg+4dqQAAAAMC1pmpAAAAA4IGFa0AAAACAaUZsQAAAACBzYW1AAAAAANCjbUAAAACAQ6duQAAAAIDsxG5AAAAAQKIlb0AAAACA5z9vQAAAAIDuw29AAAAAAJQicEAAAAAg9UZwQAAAAKCeVHBAAAAAIInEcEAAAACAcs5wQAAAAOA6BHFAAAAAADoOcUAAAAAgfcRxQAAAAOA68nFAAAAAoPExckAAAADATKVyQAAAAMC10nJAAAAAwBJHc0AAAAAA7FNzQAAAAIAYkXNAAAAAgJzRc0AAAABgE2R0QAAAAMCrc3RAAAAAgE+zdEAAAAAAhMR0QAAAAIC64XRAAAAA4Cf/dEAAAADgy3B1QAAAAMCVgnVAAAAAoDugdUAAAACAX/Z1QAAAAAAAB3ZAAAAAgO4XdkAAAABAqEJ2QAAAAIDp8nZAAAAAIAgRd0AAAACggzR3QAAAAODwT3dAAAAAIHJsd0AAAAAA4YJ3QAAAAGC91ndAAAAAAL3id0AAAADg0xZ4QAAAAOA5NXhAAAAAYDaPeEAAAACgK7p4QAAAACCm0XhAAAAA4GHkeEAAAADgyPJ4QAAAAOB5/3hAAAAAoFtNeUAAAADA8H15QAAAAGAHtXlAAAAAoD7HeUAAAAAAEtR5QAAAAGCx5XlAAAAA4NEJekAAAAAgtVN6QAAAACDCX3pAAAAAIAqmekAAAACg4OR6QAAAAGDymnxAAAAAICV/fUAAAAAgd6Z9QAAAAKDYA4BAAAAAoCoVgEAAAACA/kOAQAAAAKCpYoBAAAAA4MycgEAAAADA3DyBQAAAAKCbrIFAAAAAoPC6gUAAAADgV22CQAAAACABY4NAAAAAQE+qg0AAAADA0vKDQAAAAEDz+oNAAAAAoIxrhEAAAADg6euEQAAAAIAuDIVAAAAAoOwjhUAAAACgZUuFQAAAAADdm4VAAAAAoCzrh0AAAABgYvWHQAAAAOALCohA</binary>
                </binaryDataArray>
                <binaryDataArray encodedLength="544">
                    <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
                    <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>3FlbQDg/ZUB8w3FAV2fMQMiOnkCXfP4/T2I2QC6qskAnhOZA/NU2QCc2QEAI1UhAQcAbQRrziUBmHq5AXutSQWZDbkAZGWdAzt6lQYNptUDSFDNBoY4IQAYaQEDeT7Q/16HGP9GtXUCITrQ/Rxu0Pzhc6j9mpjZAX1X8P7tPQ0AqxS5BZTzZPye+m0B7Sa5AfPsPQRr/W0CYwBRBwDh3QMAmtD/nq6E/bJHGPxJ9UUDsy/dAoCYMQRM2a0BkAR9Boo5pQMV0VEArYu5A4kaMQAyTI0BQPRJAML3TQCKVCED85+tArObGP1BVP0EtJuVAdyKAQFjctkFQa2NBixMTQXyyjUFX8eo/IHelQTdFcEFo1zZAhagsQAO53EBIugRB0M+gQfhBgkH0MsJAbGlIQZXg+EHe6CZBsbA2QHMHOECtW6BAjE2oQUpZckBasZ1AtKl3QEZYIUHkip1AQX7TQPqF60GNuaE/USk2QGLF40Im65ZAmXqlQBGuSUC70KBAAneMQeK3aEB87MVA5NigQE/Wb0BO475A</binary>
                </binaryDataArray>
            </binaryDataArrayList>
        </spectrum>
        "#;
    let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

    assert_eq!(spectra.len(), 1);
    let s = spectra.pop().unwrap();

    assert_eq!(s.id, "spectrum=2442");
    assert_eq!(s.ms_level, 2);
    assert_eq!(s.representation, Representation::Centroid);
    assert_eq!(s.precursors.len(), 1);
    assert_eq!(s.precursors[0].charge, Some(2));
    assert!((s.precursors[0].mz - 457.723968) < 0.0001);
    assert!(match s.precursors[0].inverse_ion_mobility {
        Some(x) => (x - 1.0786) < 0.0001,
        None => false,
    });
    assert_eq!(
        s.precursors[0].isolation_window,
        Some(Tolerance::Da(-1.5, 0.75))
    );
    assert!((s.scan_start_time - 25.066).abs() < 0.0001);
    assert_eq!(s.ion_injection_time, 0.0);
    assert_eq!(s.intensity.len(), s.mz.len());
    Ok(())
}

#[tokio::test]
async fn parse_spectrum_issue_117() -> Result<(), MzMLError> {
    // The issue was that some converters write the ion mobility as part of the selected ion (as in the last test)
    // and some write it as part of the scan, as in this test. This test checks that it can be read
    // fom the scan section.
    let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum" />
            <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" />
            <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="638.352905273437955"/>
            <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="113.885513305664006"/>
            <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="793.395202636718977"/>
            <userParam name="filter string" type="xsd:string" value="ITMS + c NSI d w Full ms2 457.72@cid35.00 [115.00-930.00]"/>
            <scanList count="1">
                <cvParam cvRef="MS" accession="MS:1000795" name="no combination" />
                <scan >
                    <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="1503.96166992188" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
                    <cvParam cvRef="MS" accession="MS:1002815" name="inverse reduced ion mobility" value="1.078628" unitAccession="MS:1002814" unitName="volt-second per square centimeter"/>
                    <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:" type="xsd:double" value="457.723968505858977"/>
                    <scanWindowList count="1">
                        <scanWindow>
                            <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="115" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="930" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        </scanWindow>
                    </scanWindowList>
                </scan>
            </scanList>
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                    <activation>
                        <cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" />
                        <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="35.0"/>
                    </activation>
                </precursor>
            </precursorList>
            <binaryDataArrayList count="2">
                <binaryDataArray encodedLength="1088">
                    <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>AAAAoExpYkAAAACA3MpkQAAAAACph2VAAAAAAE4wZkAAAACAlMdmQAAAAECZAmdAAAAAwP9jaEAAAADgj4ZoQAAAAGC7HWlAAAAAAOXFaUAAAADg+4dqQAAAAMC1pmpAAAAA4IGFa0AAAACAaUZsQAAAACBzYW1AAAAAANCjbUAAAACAQ6duQAAAAIDsxG5AAAAAQKIlb0AAAACA5z9vQAAAAIDuw29AAAAAAJQicEAAAAAg9UZwQAAAAKCeVHBAAAAAIInEcEAAAACAcs5wQAAAAOA6BHFAAAAAADoOcUAAAAAgfcRxQAAAAOA68nFAAAAAoPExckAAAADATKVyQAAAAMC10nJAAAAAwBJHc0AAAAAA7FNzQAAAAIAYkXNAAAAAgJzRc0AAAABgE2R0QAAAAMCrc3RAAAAAgE+zdEAAAAAAhMR0QAAAAIC64XRAAAAA4Cf/dEAAAADgy3B1QAAAAMCVgnVAAAAAoDugdUAAAACAX/Z1QAAAAAAAB3ZAAAAAgO4XdkAAAABAqEJ2QAAAAIDp8nZAAAAAIAgRd0AAAACggzR3QAAAAODwT3dAAAAAIHJsd0AAAAAA4YJ3QAAAAGC91ndAAAAAAL3id0AAAADg0xZ4QAAAAOA5NXhAAAAAYDaPeEAAAACgK7p4QAAAACCm0XhAAAAA4GHkeEAAAADgyPJ4QAAAAOB5/3hAAAAAoFtNeUAAAADA8H15QAAAAGAHtXlAAAAAoD7HeUAAAAAAEtR5QAAAAGCx5XlAAAAA4NEJekAAAAAgtVN6QAAAACDCX3pAAAAAIAqmekAAAACg4OR6QAAAAGDymnxAAAAAICV/fUAAAAAgd6Z9QAAAAKDYA4BAAAAAoCoVgEAAAACA/kOAQAAAAKCpYoBAAAAA4MycgEAAAADA3DyBQAAAAKCbrIFAAAAAoPC6gUAAAADgV22CQAAAACABY4NAAAAAQE+qg0AAAADA0vKDQAAAAEDz+oNAAAAAoIxrhEAAAADg6euEQAAAAIAuDIVAAAAAoOwjhUAAAACgZUuFQAAAAADdm4VAAAAAoCzrh0AAAABgYvWHQAAAAOALCohA</binary>
                </binaryDataArray>
                <binaryDataArray encodedLength="544">
                    <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
                    <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>3FlbQDg/ZUB8w3FAV2fMQMiOnkCXfP4/T2I2QC6qskAnhOZA/NU2QCc2QEAI1UhAQcAbQRrziUBmHq5AXutSQWZDbkAZGWdAzt6lQYNptUDSFDNBoY4IQAYaQEDeT7Q/16HGP9GtXUCITrQ/Rxu0Pzhc6j9mpjZAX1X8P7tPQ0AqxS5BZTzZPye+m0B7Sa5AfPsPQRr/W0CYwBRBwDh3QMAmtD/nq6E/bJHGPxJ9UUDsy/dAoCYMQRM2a0BkAR9Boo5pQMV0VEArYu5A4kaMQAyTI0BQPRJAML3TQCKVCED85+tArObGP1BVP0EtJuVAdyKAQFjctkFQa2NBixMTQXyyjUFX8eo/IHelQTdFcEFo1zZAhagsQAO53EBIugRB0M+gQfhBgkH0MsJAbGlIQZXg+EHe6CZBsbA2QHMHOECtW6BAjE2oQUpZckBasZ1AtKl3QEZYIUHkip1AQX7TQPqF60GNuaE/USk2QGLF40Im65ZAmXqlQBGuSUC70KBAAneMQeK3aEB87MVA5NigQE/Wb0BO475A</binary>
                </binaryDataArray>
            </binaryDataArrayList>
        </spectrum>
        "#;
    let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

    assert_eq!(spectra.len(), 1);
    let s = spectra.pop().unwrap();
    assert!(match s.precursors[0].inverse_ion_mobility {
        Some(x) => (x - 1.0786) < 0.0001,
        None => false,
    });

    // The rest of these assertions just make sure the integrity of the spectrum is maintained
    assert_eq!(s.id, "spectrum=8678309");
    assert_eq!(s.ms_level, 2);
    assert_eq!(s.representation, Representation::Centroid);
    assert_eq!(s.precursors.len(), 1);
    assert_eq!(s.precursors[0].charge, Some(2));
    assert!((s.precursors[0].mz - 457.723968) < 0.0001);
    assert_eq!(
        s.precursors[0].isolation_window,
        Some(Tolerance::Da(-1.5, 0.75))
    );
    assert!((s.scan_start_time - 25.066).abs() < 0.0001);
    assert_eq!(s.ion_injection_time, 0.0);
    assert_eq!(s.intensity.len(), s.mz.len());
    Ok(())
}

#[tokio::test]
async fn parse_spectrum_issue_210() -> Result<(), MzMLError> {
    // Handle cases where both isolation window target m/z is set and different than selected ion m/z
    let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                </precursor>
            </precursorList>
        </spectrum>
        "#;
    let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

    assert_eq!(spectra.len(), 1);
    let s = spectra.pop().unwrap();
    assert!((s.precursors[0].mz - 457.723968) < 0.0001);
    assert_eq!(
        s.precursors[0].isolation_window,
        Some(Tolerance::Da(-1.5, 0.75))
    );

    // Check different ordering of fields in mzML
    let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1">
                <precursor>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                </precursor>
            </precursorList>
        </spectrum>
        "#;
    let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

    assert_eq!(spectra.len(), 1);
    let s = spectra.pop().unwrap();
    assert!((s.precursors[0].mz - 457.723968) < 0.0001);
    assert_eq!(
        s.precursors[0].isolation_window,
        Some(Tolerance::Da(-1.5, 0.75))
    );

    // Check fallback keeping iso window m/z of fields in mzML
    let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                </precursor>
            </precursorList>
        </spectrum>
        "#;
    let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

    assert_eq!(spectra.len(), 1);
    let s = spectra.pop().unwrap();
    assert!((s.precursors[0].mz - 457.75) < 0.0001);
    assert_eq!(
        s.precursors[0].isolation_window,
        Some(Tolerance::Da(-1.5, 0.75))
    );
    Ok(())
}
