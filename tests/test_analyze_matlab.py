"""
Tests for analyze_matlab.py -- hval filter policies and the always-written
raw diagnostic output.

Uses the synthetic detection_v2.mat builders from conftest. The
``condition_folder`` fixture provides two cells with hand-checked expected
outcomes at k_std = 2.0 (5 of 10 puncta kept under the default policy).
"""

import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import analyze_matlab  # noqa: E402
from tests.conftest import CELL1, write_detection_mat  # noqa: E402


def run_analyze(monkeypatch, argv):
    """Run analyze_matlab.main() with the given CLI args, return exit code."""
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


# ── Default behaviour is unchanged ─────────────────────────────────────

def test_default_filtered_output_unchanged(condition_folder, monkeypatch):
    code = run_analyze(
        monkeypatch, ["--input", str(condition_folder)] + BASE_ARGS
    )
    assert code == 0

    comments, header, rows = read_table(
        condition_folder / "filtered_puncta_A_values.txt"
    )
    assert comments[0] == "# Filtered puncta A values"
    assert header == ["source_image", "A_ch1", "A_ch2"]
    assert len(rows) == 5
    assert sorted(float(r[1]) for r in rows) == [3.0, 5.0, 5.0, 10.0, 20.0]
    # hval_Ar == 0 still excludes a punctum under the default policy.
    assert "cell1|row3" not in {r[0] for r in rows}
    # The comments state exactly which hval policy was used.
    assert any("hval policy: master" in c for c in comments)
    assert any("ch1 only" in c for c in comments)


# ── Raw diagnostic file ────────────────────────────────────────────────

def test_raw_file_contains_every_punctum(condition_folder, monkeypatch):
    run_analyze(monkeypatch, ["--input", str(condition_folder)] + BASE_ARGS)

    comments, header, rows = read_table(
        condition_folder / "raw_puncta_values.txt"
    )
    assert comments[0] == "# RAW CMEanalysis puncta values"
    assert comments[1] == (
        "# NO Python-side background/amplitude threshold has been applied."
    )
    assert comments[2] == "# NO hval filtering has been applied."
    assert comments[3] == (
        "# Rows contain the raw A, c, and hval_Ar values loaded from "
        "detection_v2.mat."
    )
    assert any("Lipid/master channel: ch1" in c for c in comments)
    assert any("Channels in order: ch1, ch2" in c for c in comments)
    assert any("# Rows: 10" in c for c in comments)
    assert "filtered" not in " ".join(comments).lower()

    assert header == [
        "source_image",
        "A_ch1", "c_ch1", "hval_ch1",
        "A_ch2", "c_ch2", "hval_ch2",
        "passes_A_threshold_master", "passes_hval_policy",
        "passes_final_filter",
    ]
    # Every punctum from both cells, unfiltered: 6 + 4.
    assert len(rows) == 10
    # hval == 0 rows are retained (cell1 row 3 has h = 0).
    by_source = {r[0]: r for r in rows}
    assert by_source["cell1|row3"][3] == "0"
    # A/c/hval round-trip against the synthetic input for one row.
    assert [float(v) for v in by_source["cell1|row0"][1:7]] == [
        10.0, 1.0, 1.0, 100.0, 1.0, 1.0
    ]
    # Annotations agree with the filtered result, but remove nothing.
    assert sum(int(r[-1]) for r in rows) == 5


def test_raw_columns_scale_with_channel_count(condition_folder, monkeypatch):
    """Single-channel run: the raw header is generated programmatically."""
    run_analyze(
        monkeypatch,
        ["--input", str(condition_folder),
         "--channels", "ch1", "--lipid-channel", "ch1"],
    )
    _, header, rows = read_table(condition_folder / "raw_puncta_values.txt")
    assert header == [
        "source_image", "A_ch1", "c_ch1", "hval_ch1",
        "passes_A_threshold_master", "passes_hval_policy",
        "passes_final_filter",
    ]
    assert len(rows) == 10


# ── Channel-order validation ───────────────────────────────────────────
#
# --channels must be given in MATLAB loadConditionData / frameInfo order,
# whose first entry is the master/source channel — so the first entry
# must equal --lipid-channel. Either physical folder may be the master.

def test_channel_order_valid_ch1_master(condition_folder, monkeypatch):
    """--channels ch1,ch2 --lipid-channel ch1 is accepted."""
    code = run_analyze(
        monkeypatch,
        ["--input", str(condition_folder),
         "--channels", "ch1,ch2", "--lipid-channel", "ch1"],
    )
    assert code == 0


def test_channel_order_valid_ch2_master(tmp_path, monkeypatch):
    """--channels ch2,ch1 --lipid-channel ch2 is accepted (ch2 selected
    first in MATLAB, so the detection lives under ch2/)."""
    condition = tmp_path / "condition"
    write_detection_mat(
        condition / "cell1" / "ch2" / "Detection" / "detection_v2.mat",
        CELL1["A"], CELL1["c"], CELL1["h"],
    )
    (condition / "cell1" / "ch1").mkdir(parents=True, exist_ok=True)

    code = run_analyze(
        monkeypatch,
        ["--input", str(condition),
         "--channels", "ch2,ch1", "--lipid-channel", "ch2"],
    )
    assert code == 0
    _, header, rows = read_table(condition / "filtered_puncta_A_values.txt")
    # Columns follow the given (MATLAB) order: master ch2 first.
    assert header == ["source_image", "A_ch2", "A_ch1"]
    assert len(rows) == CELL1["expected_kept"]


@pytest.mark.parametrize(
    "channels,lipid",
    [("ch1,ch2", "ch2"), ("ch2,ch1", "ch1")],
)
def test_channel_order_mismatch_rejected(
    condition_folder, monkeypatch, capsys, channels, lipid
):
    """--lipid-channel not first in --channels fails loudly, writes nothing."""
    code = run_analyze(
        monkeypatch,
        ["--input", str(condition_folder),
         "--channels", channels, "--lipid-channel", lipid],
    )
    assert code == 1
    out = capsys.readouterr().out
    # The error explains the required MATLAB selection / frameInfo order.
    assert "FIRST entry" in out
    assert "loadConditionData" in out
    assert "frameInfo" in out
    assert "master/source channel" in out
    assert "--channels ch2,ch1 --lipid-channel ch2" in out
    # Rejected before any output file is written.
    assert not (condition_folder / "raw_puncta_values.txt").exists()
    assert not (condition_folder / "filtered_puncta_A_values.txt").exists()


# ── hval policies ──────────────────────────────────────────────────────
#
# Hand-built cell: threshold = mean(c)+2*std(c) = 1.0 (std == 0).
#   row |  A_ch1 | h_ch1 | h_ch2 | A > thr
#   ----+--------+-------+-------+--------
#    0  |  10.0  |   1   |   0   |  yes     <- kept by master, none
#    1  |  10.0  |   0   |   1   |  yes     <- kept by none only
#    2  |  -0.5  |   1   |   1   |  no      <- kept by nothing; raw keeps it

POLICY_A = np.array([[10.0, 1.0], [10.0, 1.0], [-0.5, 0.0]])
POLICY_C = np.array([[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]])
POLICY_H = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])


@pytest.fixture
def policy_folder(tmp_path):
    condition = tmp_path / "condition"
    write_detection_mat(
        condition / "cell1" / "ch1" / "Detection" / "detection_v2.mat",
        POLICY_A, POLICY_C, POLICY_H,
    )
    return condition


def _run_policy(monkeypatch, folder, policy):
    return run_analyze(
        monkeypatch,
        ["--input", str(folder)] + BASE_ARGS + ["--hval-filter", policy],
    )


def test_hval_filter_master_is_default(policy_folder, monkeypatch):
    code = run_analyze(
        monkeypatch, ["--input", str(policy_folder)] + BASE_ARGS
    )
    assert code == 0
    _, _, rows = read_table(policy_folder / "filtered_puncta_A_values.txt")
    assert [r[0] for r in rows] == ["cell1|row0"]


def test_hval_filter_none_keeps_hval_zero_rows(policy_folder, monkeypatch):
    code = _run_policy(monkeypatch, policy_folder, "none")
    assert code == 0
    comments, _, rows = read_table(
        policy_folder / "filtered_puncta_A_values.txt"
    )
    assert [r[0] for r in rows] == ["cell1|row0", "cell1|row1"]
    assert any("hval policy: none" in c for c in comments)


def test_hval_filter_all_requires_every_channel(policy_folder, monkeypatch):
    """No row has hval == 1 in both channels AND A above threshold."""
    code = _run_policy(monkeypatch, policy_folder, "all")
    assert code == 1
    assert not (policy_folder / "filtered_puncta_A_values.txt").exists()
    # ... but the raw diagnostic file is still written, with every row.
    comments, _, rows = read_table(policy_folder / "raw_puncta_values.txt")
    assert len(rows) == 3
    assert any("hval policy = all" in c for c in comments)
    assert all(r[-1] == "0" for r in rows)


def test_raw_retains_negative_amplitudes(policy_folder, monkeypatch):
    _run_policy(monkeypatch, policy_folder, "master")
    _, _, rows = read_table(policy_folder / "raw_puncta_values.txt")
    assert float(rows[2][1]) == -0.5
    assert len(rows) == 3


def test_raw_annotations_follow_policy(policy_folder, monkeypatch):
    _run_policy(monkeypatch, policy_folder, "none")
    _, _, rows = read_table(policy_folder / "raw_puncta_values.txt")
    # With policy 'none' every row passes the hval annotation ...
    assert [r[-2] for r in rows] == ["1", "1", "1"]
    # ... and the final-filter annotation is the A threshold alone.
    assert [r[-1] for r in rows] == ["1", "1", "0"]
