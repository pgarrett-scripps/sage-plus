"""Build the CAMPI metaproteome search references.

Fecal samples: the sample-specific metagenome database (DB2MG) plus the human
reference, and the full catalog database written by prepare_databases.py.
SIHUMIx sample: the eight-species UniProt reference alone, and the same
reference combined with the whole gut catalog. Catalog and metagenome entries
are cleaned with the rule in prepare_databases.py.
"""

import argparse
import json
from pathlib import Path

from prepare_databases import clean, records, write


def append(handle, path, counts):
    for name, sequence in records(path):
        sequence = clean(sequence)
        if sequence is None:
            counts["dropped"] += 1
            continue
        counts["kept"] += 1
        counts["residues"] += len(sequence)
        write(handle, name, sequence)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campi", type=Path, required=True)
    parser.add_argument("--human", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    receipt = {}
    plans = {
        "gut-db2mg-human.fasta": [args.human, args.campi / "GUT_DB2MG.faa"],
        "sihumi-reference.fasta": [args.campi / "SIHUMI_DB1UNIPROT.faa"],
        "sihumi-igc.fasta": [args.campi / "SIHUMI_DB1UNIPROT.faa", args.campi / "GUT_DB1IGC.faa"],
    }
    for output, sources in plans.items():
        receipt[output] = {}
        with (args.out / output).open("w") as handle:
            for source in sources:
                counts = {"kept": 0, "dropped": 0, "residues": 0}
                append(handle, source, counts)
                receipt[output][source.name] = counts
    (args.out / "campi-databases.json").write_text(json.dumps(receipt, indent=1))


if __name__ == "__main__":
    main()
