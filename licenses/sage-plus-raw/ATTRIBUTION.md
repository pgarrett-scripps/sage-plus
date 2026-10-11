# Credits

## Upstream OpenTFRaw

This independently maintained fork derives from [OpenTFRaw](https://github.com/Sigilweaver/OpenTFRaw), originally authored by Nathan Riley and contributors. Original copyright 2026 Sigilweaver Holdings LLC and Apache-2.0 licensing are retained. See NOTICE, LICENSE and CONTRIBUTORS.md. The fork does not own the upstream opentfraw registry identity or website.

Private-branch qualification also uses publicly distributed ProteoWizard RAW and reference mzML fixtures at the revision recorded in REVIEW.md and scripts/public-fixtures.json. Those files are downloaded and hash-checked during tests, not redistributed in the Rust package.

## Prior art

### unfinnigan

Gene Selkov, 2010-2012. Perl and Python reverse-engineering of the Thermo RAW binary format.
The most thorough prior independent analysis of the format, covering versions 57, 62, 63, 64,
and 66. Field names and layout notes from unfinnigan were cross-referenced when validating
field offsets.

Source: https://github.com/prvst/unfinnigan

## Rust dependencies

- [thiserror](https://github.com/dtolnay/thiserror) -- derive macro for Error impls (David Tolnay, MIT/Apache-2.0)
