"""Tests for core.detection_io -- .mat loading, thresholding, TSV round trip."""

import numpy as np
import pytest

from core.detection_io import (
    apply_threshold,
    build_puncta_frame,
    cell_thresholds,
    find_cell_dirs,
    load_condition_detections,
    parse_channel_list,
    read_puncta_table,
    write_puncta_file,
)
from tests.conftest import CELL1, CELL2

CHANNELS = ["ch1", "ch2"]


def test_parse_channel_list():
    assert parse_channel_list("ch1,ch2") == ["ch1", "ch2"]
    assert parse_channel_list(" ch1 ") == ["ch1"]
    with pytest.raises(ValueError):
        parse_channel_list(" , ")


def test_find_cell_dirs_numeric_order(tmp_path):
    for n in (10, 2, 1):
        (tmp_path / f"cell{n}").mkdir()
    assert [p.name for p in find_cell_dirs(tmp_path)] == ["cell1", "cell2", "cell10"]


def test_load_condition_detections(condition_folder):
    det = load_condition_detections(condition_folder, CHANNELS, "ch1")

    assert det.n_cells == 2
    assert not det.errors
    assert not det.missing
    assert det.master_index == 0
    assert det.total_puncta == 10

    cell1 = det.cells[0]
    assert cell1.cell_name == "cell1"
    assert cell1.A.shape == (6, 2)
    np.testing.assert_allclose(cell1.A, CELL1["A"])
    np.testing.assert_allclose(cell1.c, CELL1["c"])


def test_load_rejects_bad_master():
    with pytest.raises(ValueError):
        load_condition_detections(".", ["ch1"], "ch2")


def test_missing_detection_recorded(condition_folder):
    target = condition_folder / "cell2" / "ch1" / "Detection" / "detection_v2.mat"
    target.unlink()
    det = load_condition_detections(condition_folder, CHANNELS, "ch1")
    assert det.n_cells == 1
    assert len(det.missing) == 1
    assert det.missing[0] == target


def test_build_puncta_frame_columns(condition_folder):
    det = load_condition_detections(condition_folder, CHANNELS, "ch1")
    frame = build_puncta_frame(det)

    assert len(frame) == 10
    for col in ("cell", "row", "source_image", "A_ch1", "A_ch2", "c", "hval_Ar"):
        assert col in frame.columns
    assert frame["source_image"].iloc[0] == "cell1|row0"
    # 'c' and 'hval_Ar' come from the master channel only.
    np.testing.assert_allclose(frame["c"].to_numpy()[:6], CELL1["c"][:, 0])


def test_cell_thresholds_match_hand_computation(condition_folder):
    det = load_condition_detections(condition_folder, CHANNELS, "ch1")
    frame = build_puncta_frame(det)
    thr = cell_thresholds(frame, k_std=2.0)

    assert thr["cell1"] == pytest.approx(CELL1["expected_threshold"], abs=1e-6)
    assert thr["cell2"] == pytest.approx(CELL2["expected_threshold"], abs=1e-12)


def test_apply_threshold_counts(condition_folder):
    det = load_condition_detections(condition_folder, CHANNELS, "ch1")
    frame = build_puncta_frame(det)
    result = apply_threshold(frame, "ch1", k_std=2.0)

    assert result.total_seen == 10
    assert result.total_kept == CELL1["expected_kept"] + CELL2["expected_kept"] == 5
    assert result.fraction_kept == pytest.approx(0.5)

    by_cell = {s.cell_name: s for s in result.per_cell}
    assert by_cell["cell1"].kept == CELL1["expected_kept"]
    assert by_cell["cell2"].kept == CELL2["expected_kept"]
    assert by_cell["cell1"].seen == 6

    # hval_Ar == 0 must exclude a punctum even when it is bright enough.
    assert "cell1|row3" not in set(result.frame["source_image"])


def test_threshold_is_monotonic_in_k(condition_folder):
    """Raising k-std can only keep fewer puncta -- the GUI slider relies on this."""
    det = load_condition_detections(condition_folder, CHANNELS, "ch1")
    frame = build_puncta_frame(det)

    counts = [apply_threshold(frame, "ch1", k).total_kept
              for k in (0.0, 1.0, 2.0, 5.0, 20.0)]
    assert counts == sorted(counts, reverse=True)
    assert counts[0] >= 5
    assert counts[-1] == 0


def test_apply_threshold_unknown_column(condition_folder):
    det = load_condition_detections(condition_folder, CHANNELS, "ch1")
    frame = build_puncta_frame(det)
    with pytest.raises(KeyError):
        apply_threshold(frame, "ch9", k_std=2.0)


def test_write_and_read_round_trip(condition_folder, tmp_path):
    det = load_condition_detections(condition_folder, CHANNELS, "ch1")
    frame = build_puncta_frame(det)
    result = apply_threshold(frame, "ch1", k_std=2.0)

    out = write_puncta_file(
        tmp_path / "filtered_puncta_A_values.txt",
        result.frame, CHANNELS, "ch1", 2.0,
    )
    text = out.read_text(encoding="utf-8")
    lines = text.splitlines()

    assert lines[0] == "# Filtered puncta A values"
    assert lines[1] == "# Lipid/master channel: ch1"
    assert lines[5] == "source_image\tA_ch1\tA_ch2"
    assert len(lines) == 6 + 5  # 5 header lines + column header + 5 rows

    back = read_puncta_table(out)
    assert list(back.columns) == ["source_image", "A_ch1", "A_ch2"]
    assert len(back) == 5
    np.testing.assert_allclose(
        sorted(back["A_ch1"].to_numpy()), [3.0, 5.0, 5.0, 10.0, 20.0]
    )
