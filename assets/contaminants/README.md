# Universal protein contaminant FASTA

`universal-contaminants.fasta` holds 381 entries: the proteins that turn up in
bottom-up proteomics samples whatever the experiment was, from skin and hair,
serum and milk in reagents, the trypsin and other enzymes added on purpose, resin
and antibody proteins from affinity steps, and the common affinity tags.

Searching it is **opt-in and off by default**. Point
`database.contaminant_fasta` at this file; nothing changes in a configuration
that does not name it. There is no keyword, no automatic default and no tag on
any accession: the proteins are ordinary targets, and Sage Plus knows they are
contaminants because they came from this file. See
[DOCS.md](../../DOCS.md#the-bundled-contaminant-fasta).

```json
{
  "database": {
    "fasta": "my-proteome.fasta",
    "contaminant_fasta": "assets/contaminants/universal-contaminants.fasta"
  }
}
```

## Where the list comes from

The membership of the list is the universal contaminant library of Frankenfield,
Ni, Ahmed and Hao, which identified these proteins across 3,667 proteomics
experiments:

> Frankenfield AM, Ni J, Ahmed M, Hao L. Protein Contaminants Matter: Building
> Universal Protein Contaminant Libraries for DDA and DIA Proteomics. *Journal of
> Proteome Research* 2022, 21(9), 2104-2113. doi:10.1021/acs.jproteome.2c00145.
> PMC10040255. Library: https://github.com/HaoGroup-ProtContLib

**The sequences here are not copied from that library.** As of 2026-11-13 that
repository states no license: it has no `LICENSE` file and the GitHub API reports
`"license": null`, so its FASTA carries no permission to redistribute. Only the
accession list was taken from it, and every sequence was downloaded again from
UniProt, which publishes its data under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

- Source of the sequences: UniProtKB release **2026_04** (4 November 2026),
  downloaded 2026-11-13 from `https://rest.uniprot.org/uniprotkb/accessions`,
  one request per accession.
- Licence of the sequences: CC BY 4.0, attributed to The UniProt Consortium.
  Cite: The UniProt Consortium, *UniProt: the Universal Protein Knowledgebase*,
  Nucleic Acids Research, https://doi.org/10.1093/nar/gkae1010.
- Headers are UniProt's own, unchanged. The source library prefixed every
  accession with `Cont_`; that prefix is **not** reproduced, because Sage Plus
  identifies a contaminant by the file it came from rather than by a keyword in
  its name.

### What differs from the 2022 library

- 379 of the 381 entries are UniProt records, rebuilt by accession. Every one of
  those accessions still resolves in release 2026_04.
- 367 of those sequences are identical to the 2022 library. 12 have been revised
  in UniProt since: `O76011`, `P60372`, `Q3LI67`, `Q9BYR0`, `A0A3Q1M032`,
  `E1BCW0`, `E1BJK2`, `G3N188`, `G5E513`, `P01030`, `P06868`, `A0A3Q1M3L6`. This
  file carries the current UniProt sequence for each, which is the point of
  rebuilding rather than copying.
- The remaining 2 entries are the FLAG and HA epitope tags. They have no UniProt
  record, so they are written with a `tag|` header that does not pretend to be
  one, from the published tag sequences.

## Checksum

```
sha256  f087348de96e5f2c316b3933e7d967b9c3f85a6f3f5bd5be5877ea5395bfe823  universal-contaminants.fasta
```

Verify with `sha256sum -c`, or:

```sh
sha256sum assets/contaminants/universal-contaminants.fasta
```

The same file ships in the release archives under `assets/contaminants/` and in the
container at `/app/assets/contaminants/universal-contaminants.fasta`.

## Rebuilding it

The list is the 381 accessions of the source library, with `Cont_` stripped. To
rebuild from a newer UniProt release, request each accession from
`https://rest.uniprot.org/uniprotkb/accessions?accessions=<accession>&format=fasta`,
keep UniProt's headers, add the two tag entries unchanged, and record the new
release, date and checksum here. An automated client should identify itself with
a generic name; no personal contact details belong in that header.
