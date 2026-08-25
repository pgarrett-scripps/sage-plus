# Raw-format test fixtures

These fixtures provide end-to-end coverage for the Thermo RAW and Bruker TDF readers.

They were copied from the MIT-licensed
[`tacular-omics/spxtacular`](https://github.com/tacular-omics/spxtacular) repository at commit
`6d5c6145095a1a8bceb036be0a7caffeb09dbfc7`.

- `thermo/Angiotensin_325-CID.raw` is the MIT-licensed fisher-py test file. It contains ten
  profile-mode FTMS CID MS2 scans of angiotensin from an Orbitrap Fusion Lumos.
- `bruker/example_dia.d` is a compact Bruker timsTOF DIA fixture containing the SQLite metadata
  and paired binary peak data required by the TDF reader.

The fixtures contain no patient or confidential data. Keep them unchanged so parser regressions
remain reproducible across platforms.
