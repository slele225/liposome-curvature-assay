"""
Tests for analyze_cpp.py -- the native (cme_detect) counterpart of
analyze_matlab.py.

Most tests use synthetic ``detection_cpp.tsv`` tables built by
``tests.conftest.write_detection_tsv`` from the same hand-checked cells
(CELL1/CELL2) as the analyze_matlab tests, so the two scripts can be
compared row for row. The end-to-end tests at the bottom run the real
``cme_detect.exe`` on synthetic images and are skipped when the native
backend has not been built.
"""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import tifffile

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import analyze_cpp  # noqa: E402
import analyze_matlab  # noqa: E402
from tests.conftest import (  # noqa: E402
    CELL1,
    CELL2,
    write_detection_mat,
    write_detection_tsv,
)

CME_DETECT_EXE = Path(
    os.environ.get(
        "CME_DETECT_EXE",
        REPO_ROOT / "native_detection" / "build" / "Release" / "cme_detect.exe",
    )
)


def run_cpp(monkeypatch, argv):
    """Run analyze_cpp.main() with the given CLI args, return exit code."""
    monkeypatch.setattr(sys, "argv", ["analyze_cpp.py"] + list(argv))
    with pytest.raises(SystemExit) as exc:
        analyze_cpp.main()
    return exc.value.code


def run_matlab(monkeypatch, argv):
    monkeypatch.setattr(sys, "argv", ["analyze_matlab.py"] + list(argv))
    with pytest.raises(SystemExit) as exc:
        analyze_matlab.main()
    return exc.value.code


def read_table(path):
    """Parse a written TXT into (comment_lines, header_cols, data_rows)."""
    comments, header, rows = [], None, []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            comments.append(line)
        elif header is None:
            header = line.split("\t")
        elif line.strip():
            rows.append(line.split("\t"))
    return comments, header, rows


BASE_ARGS = ["--channels", "ch1,ch2", "--lipid-channel", "ch1"]


# ── Parity with analyze_matlab.py ──────────────────────────────────────

def test_filtered_output_identical_to_analyze_matlab(tmp_path, monkeypatch):
    """Same puncta -> byte-identical filtered table from both backends."""
    condition = tmp_path / "condition"
    for name, spec in (("cell1", CELL1), ("cell2", CELL2)):
        write_detection_mat(
            condition / name / "ch1" / "Detection" / "detection_v2.mat",
            spec["A"], spec["c"], spec["h"],
        )
        write_detection_tsv(
            condition / name / "ch1" / "Detection" / "detection_cpp.tsv",
            spec["A"], spec["c"], spec["h"], channels=("ch1", "ch2"),
        )
        (condition / name / "ch2").mkdir(parents=True, exist_ok=True)

    assert run_matlab(monkeypatch, ["--input", str(condition)] + BASE_ARGS) == 0
    assert run_cpp(
        monkeypatch,
        ["--input", str(condition)] + BASE_ARGS
        + ["--output-name", "filtered_cpp.txt", "--raw-output-name", "raw_cpp.txt"],
    ) == 0

    mat_filtered = (condition / "filtered_puncta_A_values.txt").read_text()
    cpp_filtered = (condition / "filtered_cpp.txt").read_text()
    assert cpp_filtered == mat_filtered
    assert len(read_table(condition / "filtered_cpp.txt")[2]) == 5

    # Raw tables: identical rows and annotations; only the two provenance
    # comment lines name a different source file.
    mat_c, mat_h, mat_rows = read_table(condition / "raw_puncta_values.txt")
    cpp_c, cpp_h, cpp_rows = read_table(condition / "raw_cpp.txt")
    assert cpp_h == mat_h
    assert cpp_rows == mat_rows
    diff = [(a, b) for a, b in zip(mat_c, cpp_c) if a != b]
    assert len(diff) == 2
    assert all("detection_v2.mat" in a and "detection_cpp.tsv" in b for a, b in diff)


def test_default_filtered_output(cpp_condition_folder, monkeypatch):
    code = run_cpp(
        monkeypatch, ["--input", str(cpp_condition_folder)] + BASE_ARGS
    )
    assert code == 0

    comments, header, rows = read_table(
        cpp_condition_folder / "filtered_puncta_A_values.txt"
    )
    assert comments[0] == "# Filtered puncta A values"
    assert header == ["source_image", "A_ch1", "A_ch2"]
    assert len(rows) == 5
    assert sorted(float(r[1]) for r in rows) == [3.0, 5.0, 5.0, 10.0, 20.0]
    assert "cell1|row3" not in {r[0] for r in rows}      # hval_Ar == 0
    assert any("hval policy: master" in c for c in comments)
    assert any("ch1 only" in c for c in comments)


def test_raw_file_contains_every_punctum(cpp_condition_folder, monkeypatch):
    run_cpp(monkeypatch, ["--input", str(cpp_condition_folder)] + BASE_ARGS)

    comments, header, rows = read_table(
        cpp_condition_folder / "raw_puncta_values.txt"
    )
    assert comments[0] == "# RAW CMEanalysis puncta values"
    assert comments[3] == (
        "# Rows contain the raw A, c, and hval_Ar values loaded from "
        "detection_cpp.tsv."
    )
    assert any("Source: cell*/ch1/Detection/detection_cpp.tsv" in c for c in comments)
    assert any("# Rows: 10" in c for c in comments)
    assert "filtered" not in " ".join(comments).lower()
    assert header == [
        "source_image",
        "A_ch1", "c_ch1", "hval_ch1",
        "A_ch2", "c_ch2", "hval_ch2",
        "passes_A_threshold_master", "passes_hval_policy",
        "passes_final_filter",
    ]
    assert len(rows) == 10
    by_source = {r[0]: r for r in rows}
    assert by_source["cell1|row3"][3] == "0"
    assert [float(v) for v in by_source["cell1|row0"][1:7]] == [
        10.0, 1.0, 1.0, 100.0, 1.0, 1.0
    ]
    assert sum(int(r[-1]) for r in rows) == 5


def test_lipid_only_subset_of_channels(cpp_condition_folder, monkeypatch):
    """--channels may omit trailing slave channels of the cme_detect run."""
    code = run_cpp(
        monkeypatch,
        ["--input", str(cpp_condition_folder),
         "--channels", "ch1", "--lipid-channel", "ch1"],
    )
    assert code == 0
    _, header, rows = read_table(cpp_condition_folder / "raw_puncta_values.txt")
    assert header == [
        "source_image", "A_ch1", "c_ch1", "hval_ch1",
        "passes_A_threshold_master", "passes_hval_policy",
        "passes_final_filter",
    ]
    assert len(rows) == 10


# ── Channel-order validation ───────────────────────────────────────────

def test_channel_order_valid_ch2_master(tmp_path, monkeypatch):
    """cme_detect run with --channels ch2,ch1 --master ch2: the table lives
    under ch2/ and its columns are in ch2, ch1 order."""
    condition = tmp_path / "condition"
    write_detection_tsv(
        condition / "cell1" / "ch2" / "Detection" / "detection_cpp.tsv",
        CELL1["A"], CELL1["c"], CELL1["h"], channels=("ch2", "ch1"),
    )
    (condition / "cell1" / "ch1").mkdir(parents=True, exist_ok=True)

    code = run_cpp(
        monkeypatch,
        ["--input", str(condition),
         "--channels", "ch2,ch1", "--lipid-channel", "ch2"],
    )
    assert code == 0
    _, header, rows = read_table(condition / "filtered_puncta_A_values.txt")
    assert header == ["source_image", "A_ch2", "A_ch1"]
    assert len(rows) == CELL1["expected_kept"]


@pytest.mark.parametrize(
    "channels,lipid",
    [("ch1,ch2", "ch2"), ("ch2,ch1", "ch1")],
)
def test_cli_channel_order_mismatch_rejected(
    cpp_condition_folder, monkeypatch, capsys, channels, lipid
):
    """--lipid-channel not first in --channels fails loudly, writes nothing."""
    code = run_cpp(
        monkeypatch,
        ["--input", str(cpp_condition_folder),
         "--channels", channels, "--lipid-channel", lipid],
    )
    assert code == 1
    out = capsys.readouterr().out
    assert "FIRST entry" in out
    assert "cme_detect --channels" in out
    assert "--channels ch2,ch1 --lipid-channel ch2" in out
    assert not (cpp_condition_folder / "raw_puncta_values.txt").exists()
    assert not (cpp_condition_folder / "filtered_puncta_A_values.txt").exists()


def test_tsv_channel_order_mismatch_rejected(tmp_path, monkeypatch, capsys):
    """The TSV header records cme_detect's order; a different --channels
    list (even if internally consistent) is rejected."""
    condition = tmp_path / "condition"
    # Table written with ch2 as master, but the caller claims ch1 was.
    write_detection_tsv(
        condition / "cell1" / "ch1" / "Detection" / "detection_cpp.tsv",
        CELL1["A"], CELL1["c"], CELL1["h"], channels=("ch2", "ch1"),
    )
    code = run_cpp(
        monkeypatch,
        ["--input", str(condition),
         "--channels", "ch1,ch2", "--lipid-channel", "ch1"],
    )
    assert code == 1
    out = capsys.readouterr().out
    assert "channel order mismatch" in out
    assert "--channels ch2,ch1" in out
    assert not (condition / "raw_puncta_values.txt").exists()
    assert not (condition / "filtered_puncta_A_values.txt").exists()


def test_missing_channel_column_reported(tmp_path, monkeypatch, capsys):
    condition = tmp_path / "condition"
    write_detection_tsv(
        condition / "cell1" / "ch1" / "Detection" / "detection_cpp.tsv",
        CELL1["A"][:, :1], CELL1["c"][:, :1], CELL1["h"][:, :1],
        channels=("ch1",),
    )
    code = run_cpp(
        monkeypatch,
        ["--input", str(condition),
         "--channels", "ch1,ch2", "--lipid-channel", "ch1"],
    )
    assert code == 1
    out = capsys.readouterr().out
    assert "channel order mismatch" in out


# ── hval policies (same hand-built cell as the analyze_matlab tests) ───

POLICY_A = np.array([[10.0, 1.0], [10.0, 1.0], [-0.5, 0.0]])
POLICY_C = np.array([[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]])
POLICY_H = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])


@pytest.fixture
def policy_folder(tmp_path):
    condition = tmp_path / "condition"
    write_detection_tsv(
        condition / "cell1" / "ch1" / "Detection" / "detection_cpp.tsv",
        POLICY_A, POLICY_C, POLICY_H, channels=("ch1", "ch2"),
    )
    return condition


def _run_policy(monkeypatch, folder, policy):
    return run_cpp(
        monkeypatch,
        ["--input", str(folder)] + BASE_ARGS + ["--hval-filter", policy],
    )


def test_hval_filter_master_is_default(policy_folder, monkeypatch):
    assert run_cpp(monkeypatch, ["--input", str(policy_folder)] + BASE_ARGS) == 0
    _, _, rows = read_table(policy_folder / "filtered_puncta_A_values.txt")
    assert [r[0] for r in rows] == ["cell1|row0"]


def test_hval_filter_none_keeps_hval_zero_rows(policy_folder, monkeypatch):
    assert _run_policy(monkeypatch, policy_folder, "none") == 0
    comments, _, rows = read_table(policy_folder / "filtered_puncta_A_values.txt")
    assert [r[0] for r in rows] == ["cell1|row0", "cell1|row1"]
    assert any("hval policy: none" in c for c in comments)


def test_hval_filter_all_requires_every_channel(policy_folder, monkeypatch):
    assert _run_policy(monkeypatch, policy_folder, "all") == 1
    assert not (policy_folder / "filtered_puncta_A_values.txt").exists()
    comments, _, rows = read_table(policy_folder / "raw_puncta_values.txt")
    assert len(rows) == 3
    assert any("hval policy = all" in c for c in comments)
    assert all(r[-1] == "0" for r in rows)


def test_raw_retains_negative_amplitudes(policy_folder, monkeypatch):
    _run_policy(monkeypatch, policy_folder, "master")
    _, _, rows = read_table(policy_folder / "raw_puncta_values.txt")
    assert float(rows[2][1]) == -0.5
    assert len(rows) == 3


# ── End-to-end with the real native executable ─────────────────────────

needs_exe = pytest.mark.skipif(
    not CME_DETECT_EXE.is_file(),
    reason=f"native backend not built ({CME_DETECT_EXE}); "
           "see native_detection/README.md",
)


def _synthetic_condition(root: Path, n_cells: int = 2, size: int = 128):
    """Two-channel condition folder with bright Gaussian spots on a flat
    background (ch1 = master, twice as bright as ch2)."""
    rng = np.random.default_rng(0)
    yy, xx = np.mgrid[0:size, 0:size]
    for cell in range(1, n_cells + 1):
        imgs = {"ch1": np.full((size, size), 200.0), "ch2": np.full((size, size), 150.0)}
        for cy in range(16, size - 8, 24):
            for cx in range(16, size - 8, 24):
                x0 = cx + rng.uniform(-3, 3)
                y0 = cy + rng.uniform(-3, 3)
                g = np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * 1.5 ** 2))
                imgs["ch1"] += 1000.0 * g
                imgs["ch2"] += 500.0 * g
        for ch, img in imgs.items():
            noisy = np.clip(img + rng.normal(0, 8, img.shape), 0, 65535)
            out = root / f"cell{cell}" / ch
            out.mkdir(parents=True)
            tifffile.imwrite(str(out / "img.tif"), noisy.astype(np.uint16))
    return root


def _cme_detect(*args):
    return subprocess.run(
        [str(CME_DETECT_EXE), *args],
        capture_output=True, text=True, timeout=600,
    )


@needs_exe
def test_cme_detect_help():
    r = _cme_detect("--help")
    assert r.returncode == 0
    for flag in ("--input", "--channels", "--master", "--seed", "--sigma",
                 "--threads", "--output", "--no-matlab-layout"):
        assert flag in r.stdout


@needs_exe
def test_cme_detect_rejects_master_not_first(tmp_path):
    cond = _synthetic_condition(tmp_path / "cond", n_cells=1, size=64)
    r = _cme_detect("--input", str(cond), "--channels", "ch1,ch2",
                    "--master", "ch2", "--sigma", "1.5,1.5")
    assert r.returncode == 2
    assert "must equal the first entry of --channels" in r.stderr
    assert not (cond / "cme_detect_output").exists()


@needs_exe
def test_cme_detect_then_analyze_cpp(tmp_path, monkeypatch):
    cond = _synthetic_condition(tmp_path / "cond")
    r = _cme_detect("--input", str(cond), "--channels", "ch1,ch2",
                    "--master", "ch1", "--sigma", "1.5,1.5",
                    "--seed", "1", "--threads", "1")
    assert r.returncode == 0, r.stderr
    assert "Detections:" in r.stdout

    out_dir = cond / "cme_detect_output"
    for name in ("detections_all.tsv", "summary.tsv", "sigma.tsv"):
        assert (out_dir / name).is_file()
    sigma = (out_dir / "sigma.tsv").read_text().splitlines()
    assert sigma[0] == "channel\tsigma"
    assert sigma[1].split("\t")[0] == "ch1"

    # MATLAB layout: only the master channel gets Detection/.
    for cell in ("cell1", "cell2"):
        det = cond / cell / "ch1" / "Detection"
        assert (det / "detection_cpp.tsv").is_file()
        assert (det / "dmasks.tif").is_file()
        assert not (cond / cell / "ch2" / "Detection").exists()
        header, cols, n = analyze_cpp.read_detection_tsv(det / "detection_cpp.tsv")
        assert analyze_cpp.tsv_channel_order(header) == ["ch1", "ch2"]
        assert n > 0
        A1 = np.array([float(v) for v in cols["A_ch1"]])
        A2 = np.array([float(v) for v in cols["A_ch2"]])
        # ch1 is the brighter channel in the synthetic data.
        assert np.median(A1) > np.median(A2) > 0

    code = run_cpp(monkeypatch, ["--input", str(cond)] + BASE_ARGS)
    assert code == 0
    _, header, rows = read_table(cond / "filtered_puncta_A_values.txt")
    assert header == ["source_image", "A_ch1", "A_ch2"]
    assert rows
    _, _, raw_rows = read_table(cond / "raw_puncta_values.txt")
    assert len(rows) <= len(raw_rows)
    assert {r[0].split("|")[0] for r in raw_rows} == {"cell1", "cell2"}


@needs_exe
def test_cme_detect_ch2_master_layout(tmp_path, monkeypatch):
    """Folder names carry no master/slave meaning: ch2 may be the master."""
    cond = _synthetic_condition(tmp_path / "cond", n_cells=1)
    r = _cme_detect("--input", str(cond), "--channels", "ch2,ch1",
                    "--master", "ch2", "--sigma", "1.5,1.5", "--no-masks")
    assert r.returncode == 0, r.stderr
    det = cond / "cell1" / "ch2" / "Detection" / "detection_cpp.tsv"
    assert det.is_file()
    assert not (cond / "cell1" / "ch1" / "Detection").exists()
    header, _, _ = analyze_cpp.read_detection_tsv(det)
    assert analyze_cpp.tsv_channel_order(header) == ["ch2", "ch1"]

    assert run_cpp(
        monkeypatch,
        ["--input", str(cond), "--channels", "ch2,ch1", "--lipid-channel", "ch2"],
    ) == 0
    _, header, rows = read_table(cond / "filtered_puncta_A_values.txt")
    assert header == ["source_image", "A_ch2", "A_ch1"]
    assert rows
