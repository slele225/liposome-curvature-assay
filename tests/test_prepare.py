"""Tests for core.prepare -- metadata parsing, planning and execution."""

from pathlib import Path

import numpy as np
import pytest

from core.prepare import (
    ChannelInfo,
    center_crop_2d,
    channels_summary,
    execute_plan,
    find_tiffs,
    frame_range_warnings,
    group_name,
    parse_int_list,
    parse_metadata,
    plan_preparation,
)
from tests.conftest import write_synthetic_tiff


# ── Pure helpers ───────────────────────────────────────────────────────

def test_parse_int_list():
    assert parse_int_list("0,1,2") == [0, 1, 2]
    assert parse_int_list(" 3 , 4 ") == [3, 4]
    assert parse_int_list("") == []
    for bad in ("-1", "a", "1,-2"):
        with pytest.raises(ValueError):
            parse_int_list(bad)


def test_center_crop_2d():
    img = np.arange(64, dtype=np.uint16).reshape(8, 8)
    assert center_crop_2d(img, 1).shape == (8, 8)
    assert center_crop_2d(img, 1) is img
    cropped = center_crop_2d(img, 2)
    assert cropped.shape == (4, 4)
    # Centre quarter of an 8x8 starts at (2, 2).
    assert cropped[0, 0] == img[2, 2]


def test_group_name_and_summary():
    chans = [
        ChannelInfo("LD488", 488, 5.0, 580),
        ChannelInfo("LD561", 561, 3.2, 500),
    ]
    assert group_name(chans) == "488nm_5.0pct_580V_561nm_3.2pct_500V"
    assert channels_summary(chans) == "[488 (5.0%, 580V), 561 (3.2%, 500V)]"


# ── Metadata parsing against the real microscope files ─────────────────

def test_parse_metadata_real_tiffs(real_tiffs):
    """
    Pins current behaviour on the bundled FV3000 files.

    NOTE: these TIFFs have 3 image pages (SizeC = 3) but the parser reports 6
    channels, because the FV3000 ``Info`` block lists the microscope's whole
    channel configuration table rather than only the acquired channels. This
    test documents the existing behaviour rather than endorsing it -- see
    test_frame_range_warning_flags_channel_page_mismatch below.
    """
    channels = parse_metadata(real_tiffs[0])
    assert len(channels) == 6
    assert [c.wavelength for c in channels] == sorted(c.wavelength for c in channels)
    assert {c.laser_id for c in channels} == {"LD488", "LD561"}
    # Inactive lasers (transmissivity == 0) must be excluded.
    assert all(c.transmissivity > 0 for c in channels)


def test_plan_is_read_only(real_tiffs, tmp_path):
    """Planning must never create the output tree."""
    out = tmp_path / "never_created"
    plan = plan_preparation(real_tiffs[0].parent, out, None, 1)
    assert not out.exists()
    assert len(plan.items) == len(real_tiffs)
    assert plan.output_paths, "plan should list files it would write"
    assert all(str(p).startswith(str(out)) for p in plan.output_paths)


def test_frame_range_warning_flags_channel_page_mismatch(real_tiffs, tmp_path):
    """The 6-channels-vs-3-pages mismatch is surfaced before any write."""
    plan = plan_preparation(real_tiffs[0].parent, tmp_path / "out", None, 1)
    warnings = frame_range_warnings(plan)
    assert warnings, "expected out-of-range frames to be reported"
    assert "frames out of range" in warnings[0]


def test_execute_records_error_instead_of_raising(real_tiffs, tmp_path):
    """A file whose frames exceed its pages is reported, not fatal."""
    out = tmp_path / "out"
    plan = plan_preparation(real_tiffs[0].parent, out, None, 1)
    result = execute_plan(plan)
    assert result.errors
    assert result.processed == 0
    assert any("out of range" in msg for _p, msg in result.errors)


# ── Round trip on a synthetic two-channel TIFF ─────────────────────────

def test_full_round_trip(tmp_path):
    src = tmp_path / "raw"
    src.mkdir()
    channels = [("LD488", 488, 5.0, 580), ("LD561", 561, 3.2, 500)]
    for i in range(3):
        write_synthetic_tiff(src / f"img{i:03d}.tif", channels, n_pages=2)

    assert len(find_tiffs(src)) == 3

    out = tmp_path / "prepared"
    plan = plan_preparation(src, out, None, crop_divisor=1)

    assert not frame_range_warnings(plan)
    assert len(plan.planned) == 3
    group = "488nm_5.0pct_580V_561nm_3.2pct_500V"
    assert plan.group_counts == {group: 3}
    assert [it.cell_index for it in plan.planned] == [1, 2, 3]
    assert len(plan.output_paths) == 6  # 3 files x 2 channels

    result = execute_plan(plan)
    assert result.processed == 3
    assert not result.errors
    assert len(result.written) == 6

    for p in result.written:
        assert p.is_file()
    assert (out / group / "cell1" / "ch1" / "img000_ch1.tif").is_file()
    assert (out / group / "cell1" / "ch2" / "img000_ch2.tif").is_file()
    assert (out / group / "cell3" / "ch2" / "img002_ch2.tif").is_file()


def test_crop_divisor_applied(tmp_path):
    src = tmp_path / "raw"
    src.mkdir()
    channels = [("LD488", 488, 5.0, 580)]
    write_synthetic_tiff(src / "a.tif", channels, n_pages=1, size=32)

    out = tmp_path / "prepared"
    plan = plan_preparation(src, out, None, crop_divisor=2)
    result = execute_plan(plan)

    assert result.processed == 1
    import tifffile
    with tifffile.TiffFile(str(result.written[0])) as tif:
        assert tif.pages[0].shape == (16, 16)


def test_frames_length_mismatch_is_recorded(tmp_path):
    src = tmp_path / "raw"
    src.mkdir()
    channels = [("LD488", 488, 5.0, 580), ("LD561", 561, 3.2, 500)]
    write_synthetic_tiff(src / "a.tif", channels, n_pages=2)

    plan = plan_preparation(src, tmp_path / "out", [0], crop_divisor=1)
    assert plan.failed
    assert "--frames has 1 entries" in plan.failed[0].reason


def test_no_tiffs_raises(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(FileNotFoundError):
        plan_preparation(empty, tmp_path / "out", None, 1)
