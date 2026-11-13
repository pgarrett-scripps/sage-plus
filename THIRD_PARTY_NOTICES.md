# Third-party notices

Sage Plus reads Bruker TDF files with `sage-plus-tdf`, a fork of MannLabs TimsRust 0.6.6
licensed `Apache-2.0 AND MIT`. It reads Thermo RAW files with `sage-plus-raw`, a fork of
Sigilweaver OpenTFRaw licensed Apache-2.0. Their license texts, notices and source
attribution are shipped in `licenses/sage-plus-tdf/` and `licenses/sage-plus-raw/`.
The Sage Plus license remains in the root `LICENSE` file.

## Unimod

Sage Plus compiles modification names and masses from the Unimod database
(https://www.unimod.org) into the executable. The data, in
`crates/sage/data/unimodifications.json`, is a JSON extract of Unimod's
`unimod.xml` (data version 17:06:2025 11:31) and carries Unimod's notice:
"Copyright (C) 2002-2006 Unimod; this information may be copied, distributed
and/or modified under certain conditions, but it comes WITHOUT ANY WARRANTY;
see the accompanying Design Science License for more details." The extract
remains under the Design Science License. Its source data is that JSON file in
the Sage Plus source distribution.

The notice, a description of the conversion, and the full Design Science
License text are in `crates/sage/data/LICENSE-unimod.txt`. Executable archives
include it as `LICENSE-unimod.txt`, and containers include it at
`/app/licenses/unimod.txt`.

## UniProt (bundled contaminant FASTA)

`assets/contaminants/universal-contaminants.fasta` is an opt-in contaminant
database, searched only when `database.contaminant_fasta` names it. 379 of its
381 entries are protein sequences from UniProtKB release 2026_04 (4 November
2026), which The UniProt Consortium publishes under the Creative Commons
Attribution 4.0 International License (CC BY 4.0,
https://creativecommons.org/licenses/by/4.0/). They are reproduced with
UniProt's own headers and are unmodified. Attribution: The UniProt Consortium,
*UniProt: the Universal Protein Knowledgebase*, Nucleic Acids Research,
https://doi.org/10.1093/nar/gkae1010.

Which proteins to include came from the universal contaminant library of
Frankenfield AM, Ni J, Ahmed M, Hao L, *Protein Contaminants Matter: Building
Universal Protein Contaminant Libraries for DDA and DIA Proteomics*, Journal of
Proteome Research 2022, 21(9), 2104-2113, doi:10.1021/acs.jproteome.2c00145.
That library states no license, so none of its files are redistributed here:
only the list of accessions was used, and every sequence was downloaded from
UniProt. The remaining two entries are the FLAG and HA epitope tags, which have
no UniProt record.

`assets/contaminants/README.md` records the release, download date, checksum and
rebuild procedure. Executable archives include the directory as
`assets/contaminants/`, and containers include it at `/app/assets/contaminants/`.
