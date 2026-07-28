"""Shared fixtures and synthetic-data builders for the core-library tests."""

import sys
from pathlib import Path

import h5py
import numpy as np
import pytest
import tifffile

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture(scope="session")
def test_data(repo_root: Path) -> Path:
    d = repo_root / "test_data"
    if not d.is_dir():
        pytest.skip("test_data/ not present")
    return d


@pytest.fixture(scope="session")
def sample_dls(test_data: Path) -> Path:
    p = test_data / "sample_dls.xlsx"
    if not p.is_file():
        pytest.skip("test_data/sample_dls.xlsx not present")
    return p


@pytest.fixture(scope="session")
def real_tiffs(test_data: Path):
    tiffs = sorted((test_data / "tiff_images").glob("*.tif"))
    if not tiffs:
        pytest.skip("test_data/tiff_images/ has no TIFFs")
    return tiffs


# ── Synthetic builders ─────────────────────────────────────────────────

def build_fv3000_info(channels, n_lambda: int = 2) -> str:
    """
    Build a minimal Olympus-FV3000-shaped ``Info`` metadata blob.

    ``channels`` is a list of ``(laser_id, wavelength, transmissivity, voltage)``.
    A couple of lambda-phase entries are appended so tests exercise the
    skip-lambda branch.
    """
    lines = []
    seen = {}
    for laser_id, wl, trans, _v in channels:
        seen[laser_id] = (wl, trans)
    # Include an inactive laser so the transmissivity filter is exercised.
    lines.append("- Laser LD405 transmissivity = 0.0")
    for laser_id, (wl, trans) in seen.items():
        lines.append(f"- Laser {laser_id} transmissivity = {trans}")

    idx = 0
    for laser_id, _wl, _t, _v in channels:
        idx += 1
        lines.append(
            f" channel laserDataId #{idx:02d} = {laser_id}_65794Imaging_main_phase_1"
        )
    for i in range(n_lambda):
        idx += 1
        lines.append(
            f" channel laserDataId #{idx:02d} = "
            f"LD405_65792Imaging_main_phase_Lambda"
        )

    for i, (_l, _w, _t, volt) in enumerate(channels, start=1):
        lines.append(f" pmt voltage #{i:02d} = {volt}")
    for i in range(len(channels) + 1, idx + 1):
        lines.append(f" pmt voltage #{i:02d} = 0")

    return "\n".join(lines)


def write_synthetic_tiff(path: Path, channels, n_pages: int, size: int = 32):
    """Write a multi-page ImageJ TIFF carrying FV3000-shaped metadata."""
    rng = np.random.default_rng(0)
    stack = rng.integers(0, 4096, size=(n_pages, size, size), dtype=np.uint16)
    tifffile.imwrite(
        str(path), stack, imagej=True,
        metadata={"Info": build_fv3000_info(channels)},
    )
    return path


def write_detection_mat(path: Path, A: np.ndarray, c: np.ndarray, hval: np.ndarray):
    """
    Write a MATLAB-v7.3-shaped ``detection_v2.mat``.

    ``frameInfo`` fields are stored as HDF5 object references, matching what
    cmeAnalysis produces for a multi-frame struct array.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(str(path), "w") as f:
        refs = f.create_group("#refs#")
        fi = f.create_group("frameInfo")
        for name, arr in (("A", A), ("c", c), ("hval_Ar", hval)):
            ds = refs.create_dataset(f"{name}_0", data=np.asarray(arr, dtype=float))
            ref_ds = fi.create_dataset(name, shape=(1, 1), dtype=h5py.ref_dtype)
            ref_ds[0, 0] = ds.ref
    return path


# Two cells with hand-checked expected outcomes at k_std = 2.0.
CELL1 = {
    "A": np.array([[10.0, 100.0], [0.1, 1.0], [5.0, 50.0],
                   [4.0, 40.0], [0.2, 2.0], [20.0, 200.0]]),
    "c": np.array([[1.0, 1.0], [1.0, 1.0], [1.0, 1.0],
                   [1.0, 1.0], [1.0, 1.0], [5.0, 5.0]]),
    "h": np.array([[1.0, 1.0], [1.0, 1.0], [1.0, 1.0],
                   [0.0, 0.0], [1.0, 1.0], [1.0, 1.0]]),
    "expected_threshold": 4.648090,   # mean(c)+2*std(c) = 1.666667+2*1.490712
    "expected_kept": 3,               # rows 0, 2, 5
}
CELL2 = {
    "A": np.array([[3.0, 30.0], [1.0, 10.0], [5.0, 50.0], [2.0, 20.0]]),
    "c": np.array([[2.0, 2.0], [2.0, 2.0], [2.0, 2.0], [2.0, 2.0]]),
    "h": np.array([[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]),
    "expected_threshold": 2.0,        # std == 0
    "expected_kept": 2,               # rows 0, 2 (strictly greater than 2.0)
}


@pytest.fixture
def condition_folder(tmp_path: Path) -> Path:
    """A two-cell, two-channel condition folder with synthetic detections."""
    condition = tmp_path / "488nm_3.0pct_580V_561nm_3.2pct_500V"
    for name, spec in (("cell1", CELL1), ("cell2", CELL2)):
        write_detection_mat(
            condition / name / "ch1" / "Detection" / "detection_v2.mat",
            spec["A"], spec["c"], spec["h"],
        )
        (condition / name / "ch2").mkdir(parents=True, exist_ok=True)
    return condition
