#!/usr/bin/env python3
"""Write Bruker SDK scan-to-1/K0 reference values for the TimsCalibration unit tests.

Requires the Bruker `libtimsdata` library through the `tdfpy` package, which is not a Sage
Plus dependency. Values come from `tims_scannum_to_oneoverk0` for every calibration row of
each acquisition, at fractional scans, scans beyond the acquisition range, and scans beyond
the C8 and C9 model limits.

    python3 benchmarks/generate_tims_calibration_fixture.py \
        --tdfpy /path/to/site-packages --data /path/to/tdfpy/tests/data
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUTPUT = REPO / "crates/sage-cloudpath/tests/data/bruker/tims_calibration_sdk.json"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tdfpy", type=Path, required=True,
                        help="directory containing the tdfpy package")
    parser.add_argument("--data", type=Path, required=True,
                        help="directory containing example_dda.d, example_dia.d, example_prm.d")
    args = parser.parse_args()
    sys.path.insert(0, str(args.tdfpy))
    import numpy as np
    from tdfpy.timsdata import TimsData

    cases = []
    for name in ["example_dda", "example_dia", "example_prm"]:
        path = args.data / f"{name}.d"
        tdf = path / "analysis.tdf"
        db = sqlite3.connect(tdf)
        rows = db.execute(
            "SELECT Id, ModelType, C0, C1, C2, C3, C4, C5, C6, C7, C8, C9 FROM TimsCalibration"
        ).fetchall()
        frames = dict(db.execute(
            "SELECT TimsCalibration, min(Id) FROM Frames GROUP BY TimsCalibration"
        ).fetchall())
        scan_count = db.execute("SELECT max(NumScans) FROM Frames").fetchone()[0]
        sdk = TimsData(str(path))
        for calibration_id, model_type, *c in rows:
            slope = (c[2] - c[3]) / c[1]
            # Scans where the intermediate value reaches the C8 and C9 limits.
            at_c8 = c[4] + c[0] + (c[2] - c[8]) / slope
            at_c9 = c[4] + c[0] + (c[2] - c[9]) / slope
            scans = sorted({
                0.0, 0.5, 1.0, 17.25, scan_count / 2, scan_count - 1.0, float(scan_count),
                -25.0, scan_count + 25.0,
                at_c8 - 300, at_c8 - 1, at_c8 + 1, at_c8 + 300,
                at_c9 - 300, at_c9 - 1, at_c9 + 1, at_c9 + 300,
            })
            values = sdk.scanNumToOneOverK0(frames[calibration_id], np.array(scans, dtype=np.float64))
            cases.append({
                "source": f"{name}.d",
                "tdf_sha256": hashlib.sha256(tdf.read_bytes()).hexdigest(),
                "calibration_id": calibration_id,
                "model_type": model_type,
                "coefficients": c,
                "scans": scans,
                "one_over_k0": [float(value) for value in values],
            })
    OUTPUT.write_text(json.dumps({
        "description": "Bruker libtimsdata tims_scannum_to_oneoverk0 reference values. Scans "
                       "include fractional values, points beyond the acquisition range, and "
                       "points beyond the C8 and C9 model limits.",
        "generator": "benchmarks/generate_tims_calibration_fixture.py",
        "cases": cases,
    }, indent=1))
    print(f"wrote {len(cases)} cases to {OUTPUT.relative_to(REPO)}")


if __name__ == "__main__":
    main()
