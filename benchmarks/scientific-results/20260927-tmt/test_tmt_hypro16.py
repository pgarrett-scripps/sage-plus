"""Missing reporter channels: 0.0 (upstream, Sage Plus <= Beta 12), null (Beta 13) and NaN are all missing.

    UV_CACHE_DIR=/mnt/data1/build-cache/uv-cache uv run --no-project --with pyarrow --with numpy \
        python benchmarks/scientific-results/20260927-tmt/test_tmt_hypro16.py
"""
import sys
import tempfile
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).parent))
import tmt_hypro16 as t  # noqa: E402


def test_missing_channels():
    path = Path(tempfile.mkdtemp()) / "x.parquet"
    rep = [[None, 0.0, float("nan")] + [1.0] * 13, None, [5.0] * 16]
    n = len(rep)
    pq.write_table(pa.table({
        "filename": ["a.mzML"] * n, "scannr": ["controllerType=0 controllerNumber=1 scan=1", "scan=2", "3"],
        "peptide": ["PEPTIDEK"] * n, "stripped_peptide": ["PEPTIDEK"] * n, "proteins": ["sp|X|A_HUMAN"] * n,
        "is_decoy": [False] * n, "rank": [1] * n, "spectrum_q": [0.0] * n, "peptide_q": [0.0] * n,
        "protein_q": [0.0] * n, "reporter_ion_intensity": pa.array(rep, type=pa.list_(pa.float32())),
    }), path)
    rows = t.load_psms(path)
    assert np.isnan(rows["rep"][0, :3]).all() and (rows["rep"][0, 3:] == 1.0).all()
    assert list(rows["has_reporter"]) == [True, False, True]
    assert rows["scan"] == [1, 2, 3]


if __name__ == "__main__":
    test_missing_channels()
    print("ok")
