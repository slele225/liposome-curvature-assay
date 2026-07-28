"""
Regression tests for ID-driven channel resolution in prepare_input.py.

The metadata below is the verbatim (relevant) portion of the ``Info`` block
from ``data/olivia_data/dataHis eGFP SLiC post exess removal.tif``. All 21
files in that folder carry an identical channel/detector configuration.

Ground truth for this acquisition:

    Page | Channel | Laser | Detector      | Voltage | Dye       | Role
    -----+---------+-------+---------------+---------+-----------+---------
     0   | CH1     | 488   | FV30-SD_D_1   | 652     | EGFP      | protein
     1   | CH2     | 561   | FV30-SD_D_2   | 876     | Texas Red | lipid
     2   | CH3     | 561   | FV31-LETD     | 213     | --        | DIC

Detector #1 links channel ID 716c811e-..., which is not CH1/CH2/CH3 -- it is
a configured-but-unused detector and must never surface in the output.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pytest
import tifffile

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from prepare_input import (  # noqa: E402
    MetadataError,
    assign_slots,
    build_channel_map,
    channel_order_description,
    group_name,
    override_warnings,
    process_images,
    resolve_channels,
)


REAL_INFO = """\
 BitsPerPixel = 12
 DimensionOrder = XYCZT
 PixelType = uint16
 SizeC = 3
 SizeT = 1
 SizeX = 512
 SizeY = 512
 SizeZ = 1
- Channel CH1 ID = caa3e599-a28b-4656-95e3-9838e3db2af7
- Channel CH1 linked laser index = 1
- Channel CH2 ID = 7db1983d-bdd3-4b26-8313-0421c2b486fa
- Channel CH2 linked laser index = 2
- Channel CH3 ID = 3e0774a9-a564-4bc8-b3f3-5deebf943ecf
- Channel CH3 linked laser index = 2
- Detector ID #1 = __FV30-SD_D_1__
- Detector ID #2 = __FV31-LETD_D_1__
- Detector ID #3 = __FV30-SD_D_1__
- Detector ID #4 = __FV30-SD_D_2__
- Detector linked channel ID #1 = 716c811e-e24e-4e3b-958d-482d101dddd5
- Detector linked channel ID #2 = 3e0774a9-a564-4bc8-b3f3-5deebf943ecf
- Detector linked channel ID #3 = caa3e599-a28b-4656-95e3-9838e3db2af7
- Detector linked channel ID #4 = 7db1983d-bdd3-4b26-8313-0421c2b486fa
- Detector voltage #1 = 0.0
- Detector voltage #2 = 213.0
- Detector voltage #3 = 652.0
- Detector voltage #4 = 876.0
- Laser LD405 ID = LD405_65792
- Laser LD405 data ID = LD405_65792Imaging_main_phase_1
- Laser LD405 transmissivity = 1.0
- Laser LD488 ID = LD488_65794
- Laser LD488 data ID = LD488_65794Imaging_main_phase_1
- Laser LD488 transmissivity = 17.7
- Laser LD488 wavelength = 488.0
- Laser LD561 ID = LD561_65796
- Laser LD561 data ID = LD561_65796Imaging_main_phase_1
- Laser LD561 transmissivity = 15.3
- Laser LD561 wavelength = 561.0
- Laser LD640 ID = LD640_65798
- Laser LD640 data ID = LD640_65798Imaging_main_phase_1
- Laser LD640 transmissivity = 1.0
channel deviceName #1 = SD1
channel deviceName #2 = SD1
channel deviceName #3 = SD1
channel deviceName #4 = TD
channel dyeId #1 = EGFP_DYE
channel dyeId #2 = Texas Red_DYE
channel dyeName #1 = EGFP
channel dyeName #2 = Texas Red
channel laserDataId #01 = LD488_65794Imaging_main_phase_1
channel laserDataId #02 = LD561_65796Imaging_main_phase_1
channel laserDataId #03 = LD561_65796Imaging_main_phase_1
channel laserDataId #04 = LD405_65792Imaging_main_phase_Lambda
channel laserDataId #05 = LD488_65794Imaging_main_phase_Lambda
channel name #1 = CH1
channel name #2 = CH2
channel name #3 = CH3
dyeData excitationWavelength #1 = 488
dyeData excitationWavelength #2 = 595
dyeData name #1 = EGFP
dyeData name #2 = Texas Red
laser name #1 = LD405
laser name #2 = LD488
laser name #3 = LD561
laser name #4 = LD640
"""

PHANTOM_DETECTOR_UUID = "716c811e-e24e-4e3b-958d-482d101dddd5"


# Three genuine fluorescence channels, no DIC, so every resolved channel is a
# candidate slot. Used to exercise --frames naming a subset.
THREE_CHANNEL_INFO = """\
 SizeC = 3
- Channel CH1 ID = aaaaaaaa-0001
- Channel CH1 linked laser index = 0
- Channel CH2 ID = bbbbbbbb-0002
- Channel CH2 linked laser index = 1
- Channel CH3 ID = cccccccc-0003
- Channel CH3 linked laser index = 2
- Detector ID #1 = __FV30-SD_D_3__
- Detector ID #2 = __FV30-SD_D_1__
- Detector ID #3 = __FV30-SD_D_2__
- Detector linked channel ID #1 = cccccccc-0003
- Detector linked channel ID #2 = aaaaaaaa-0001
- Detector linked channel ID #3 = bbbbbbbb-0002
- Detector voltage #1 = 700.0
- Detector voltage #2 = 500.0
- Detector voltage #3 = 600.0
- Laser LD405 transmissivity = 2.0
- Laser LD405 wavelength = 405.0
- Laser LD488 transmissivity = 5.0
- Laser LD488 wavelength = 488.0
- Laser LD561 transmissivity = 3.2
- Laser LD561 wavelength = 561.0
laser name #1 = LD405
laser name #2 = LD488
laser name #3 = LD561
"""

THREE_CHANNEL_FULL_NAME = "405nm_2.0pct_500V_488nm_5.0pct_600V_561nm_3.2pct_700V"
THREE_CHANNEL_SUBSET_NAME = "488nm_5.0pct_600V_405nm_2.0pct_500V"


# ── The resolved mapping ───────────────────────────────────────────────

def test_resolves_ground_truth_table():
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)

    assert len(channels) == 2, "DIC channel must not become a ch<i> folder"

    ch1, ch2 = channels
    assert (ch1["page"], ch1["channel"], ch1["laser_id"]) == (0, "CH1", "LD488")
    assert (ch1["wavelength"], ch1["transmissivity"]) == (488, 17.7)
    assert (ch1["detector_id"], ch1["voltage"]) == ("FV30-SD_D_1", 652)

    assert (ch2["page"], ch2["channel"], ch2["laser_id"]) == (1, "CH2", "LD561")
    assert (ch2["wavelength"], ch2["transmissivity"]) == (561, 15.3)
    assert (ch2["detector_id"], ch2["voltage"]) == ("FV30-SD_D_2", 876)

    assert len(excluded) == 1
    assert excluded[0]["channel"] == "CH3"
    assert excluded[0]["page"] == 2
    assert excluded[0]["detector_id"] == "FV31-LETD_D_1"
    assert excluded[0]["voltage"] == 213
    assert "LETD" in excluded[0]["reason"]


def test_channels_ordered_by_wavelength_ascending():
    channels, _ = resolve_channels(REAL_INFO, n_pages=3)
    wavelengths = [c["wavelength"] for c in channels]
    assert wavelengths == sorted(wavelengths)
    assert wavelengths[0] == 488, "ch1 must be the lowest-wavelength channel"


def test_group_name_reflects_corrected_pairings():
    channels, _ = resolve_channels(REAL_INFO, n_pages=3)
    assert group_name(channels) == "488nm_17.7pct_652V_561nm_15.3pct_876V"


def test_phantom_detector_is_absent_everywhere():
    """Detector #1 links no acquired channel: it must not leak into output."""
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)

    # Its 0 V reading was what the old positional pairing surfaced.
    assert all(c["voltage"] != 0 for c in channels + excluded)

    # Only slots 2, 3 and 4 describe real channels; slot 1 is never read.
    assert {c["detector_id"] for c in channels + excluded} == {
        "FV30-SD_D_1", "FV30-SD_D_2", "FV31-LETD_D_1"
    }

    blob = repr(build_channel_map(channels, excluded, "g"))
    assert PHANTOM_DETECTOR_UUID not in blob


def test_dye_names_resolve_from_the_gated_channel_block():
    """
    Both dyes resolve: the `channel dyeName` block has one entry per
    fluorescence channel, and its EGFP entry agrees with the exact
    488 nm excitation match that can be checked independently.
    """
    channels, _ = resolve_channels(REAL_INFO, n_pages=3)
    assert [c["dye"] for c in channels] == ["EGFP", "Texas Red"]
    assert all("dyeName" in c["dye_source"] for c in channels)


def test_dye_block_discarded_when_it_contradicts_an_exact_match():
    """One disagreement throws the whole positional block away."""
    info = REAL_INFO.replace(
        "channel dyeName #1 = EGFP", "channel dyeName #1 = Texas Red"
    ).replace("channel dyeName #2 = Texas Red", "channel dyeName #2 = EGFP")
    channels, _ = resolve_channels(info, n_pages=3)
    # Falls back to the exact excitation match, which only covers 488.
    assert channels[0]["dye"] == "EGFP"
    assert channels[1]["dye"] is None


def test_dye_block_ignored_when_it_does_not_align():
    """Wrong entry count -> the block is ragged and cannot be trusted."""
    info = REAL_INFO.replace(
        "channel dyeName #2 = Texas Red",
        "channel dyeName #2 = Texas Red\nchannel dyeName #3 = Cy5",
    )
    channels, _ = resolve_channels(info, n_pages=3)
    assert channels[0]["dye"] == "EGFP"
    assert channels[1]["dye"] is None


# ── --frames must not desync the manifest from the pixels ──────────────

def test_frames_override_carries_channel_metadata_into_the_slot():
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)
    slots, slot_excluded = assign_slots(channels, excluded, [1, 0])
    manifest = build_channel_map(slots, excluded, group_name(slots), "1,0")

    ch1 = manifest["channels"]["ch1"]
    assert ch1["source_page"] == 1
    assert ch1["wavelength_nm"] == 561
    assert ch1["pmt_voltage_v"] == 876
    assert ch1["channel_name"] == "CH2"
    assert ch1["detector_id"] == "FV30-SD_D_2"
    assert ch1["dye"] == "Texas Red"

    ch2 = manifest["channels"]["ch2"]
    assert ch2["source_page"] == 0
    assert ch2["wavelength_nm"] == 488
    assert ch2["pmt_voltage_v"] == 652
    assert ch2["channel_name"] == "CH1"
    assert ch2["detector_id"] == "FV30-SD_D_1"
    assert ch2["dye"] == "EGFP"


def test_channel_order_does_not_claim_ascending_under_override():
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)
    slots, slot_excluded = assign_slots(channels, excluded, [1, 0])
    manifest = build_channel_map(slots, excluded, group_name(slots), "1,0")

    assert "NOT wavelength-ascending" in manifest["channel_order"]
    assert "1,0" in manifest["channel_order"]
    assert manifest["frames_argument"] == "1,0"


def test_default_order_is_unchanged():
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)
    slots, slot_excluded = assign_slots(channels, excluded, None)
    manifest = build_channel_map(slots, excluded, group_name(slots), None)

    assert slots == channels
    assert manifest["channel_order"] == (
        "excitation wavelength ascending (ch1 = lowest)"
    )
    assert manifest["frames_argument"] is None
    assert manifest["channels"]["ch1"]["source_page"] == 0
    assert manifest["channels"]["ch1"]["wavelength_nm"] == 488
    assert manifest["channels"]["ch2"]["source_page"] == 1
    assert manifest["channels"]["ch2"]["wavelength_nm"] == 561


def test_frames_matching_natural_order_is_not_flagged():
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)
    slots, slot_excluded = assign_slots(channels, excluded, [0, 1])
    description = channel_order_description(slots, "0,1")

    assert "NOT wavelength-ascending" not in description
    assert "coincides with wavelength-ascending" in description
    assert not override_warnings(slots, "0,1", "g")


def test_non_ascending_override_warns():
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)
    slots, slot_excluded = assign_slots(channels, excluded, [1, 0])
    warnings = override_warnings(slots, "1,0", "g")

    assert warnings
    assert "non-wavelength-ascending" in warnings[0]
    assert not override_warnings(channels, None, "g")


def test_frames_forcing_an_excluded_channel_is_flagged():
    """Writing the DIC page into a slot is allowed but never silent."""
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)
    slots, slot_excluded = assign_slots(channels, excluded, [0, 2])
    manifest = build_channel_map(slots, slot_excluded, group_name(slots), "0,2")

    assert manifest["channels"]["ch2"]["channel_name"] == "CH3"
    assert "LETD" in manifest["channels"]["ch2"]["excluded_by_rule"]
    # CH3 is written; CH2 (page 1) was not named, so it takes its place.
    assert [c["channel_name"] for c in manifest["excluded_channels"]] == ["CH2"]
    assert any("non-fluorescence" in w for w in override_warnings(slots, "0,2", "g"))


def test_frames_arity_and_range_are_validated():
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)
    with pytest.raises(ValueError, match="only 2 channels"):
        assign_slots(channels, excluded, [0, 1, 2])
    with pytest.raises(IndexError, match="page 9"):
        assign_slots(channels, excluded, [0, 9])
    with pytest.raises(ValueError, match="repeats"):
        assign_slots(channels, excluded, [1, 1])
    with pytest.raises(ValueError, match="at least one"):
        assign_slots(channels, excluded, [])


def test_group_name_follows_slot_order():
    """ch1's laser comes first in the folder name, so --frames reorders it."""
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)

    default = group_name(assign_slots(channels, excluded, None)[0])
    remapped = group_name(assign_slots(channels, excluded, [1, 0])[0])

    assert default == "488nm_17.7pct_652V_561nm_15.3pct_876V"
    assert remapped == "561nm_15.3pct_876V_488nm_17.7pct_652V"
    assert default != remapped


def test_group_name_omits_excluded_channels():
    """The DIC channel is not a slot, so it never reaches the folder name."""
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)
    name = group_name(assign_slots(channels, excluded, None)[0])

    assert excluded[0]["voltage"] == 213
    assert "213V" not in name
    assert name.count("nm_") == 2


def test_group_name_describes_advertises_slot_order():
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)
    slots, slot_excluded = assign_slots(channels, excluded, [1, 0])
    manifest = build_channel_map(slots, excluded, group_name(slots), "1,0")

    assert manifest["group"] == "561nm_15.3pct_876V_488nm_17.7pct_652V"
    assert "slot" in manifest["group_name_describes"]
    assert "NOT the slot order" not in manifest["group_name_describes"]


# ── --frames naming a subset of the resolved channels ──────────────────

def test_three_channel_baseline():
    """All three channels resolve and are wavelength-ascending by default."""
    channels, excluded = resolve_channels(THREE_CHANNEL_INFO, n_pages=3)
    assert not excluded
    assert [c["channel"] for c in channels] == ["CH1", "CH2", "CH3"]
    assert [c["wavelength"] for c in channels] == [405, 488, 561]
    assert [c["voltage"] for c in channels] == [500, 600, 700]
    assert group_name(channels) == THREE_CHANNEL_FULL_NAME


def test_frames_subset_keeps_only_the_named_pages():
    channels, excluded = resolve_channels(THREE_CHANNEL_INFO, n_pages=3)
    slots, slot_excluded = assign_slots(channels, excluded, [1, 0])

    assert len(slots) == 2
    assert [c["page"] for c in slots] == [1, 0]
    assert [c["channel"] for c in slots] == ["CH2", "CH1"]
    assert [c["wavelength"] for c in slots] == [488, 405]
    assert [c["voltage"] for c in slots] == [600, 500]


def test_frames_subset_group_name_drops_the_unnamed_channel():
    channels, excluded = resolve_channels(THREE_CHANNEL_INFO, n_pages=3)
    slots, _ = assign_slots(channels, excluded, [1, 0])

    assert group_name(slots) == THREE_CHANNEL_SUBSET_NAME
    # The dropped channel's 561 nm / 700 V never reach the folder name.
    assert "561nm" not in group_name(slots)
    assert "700V" not in group_name(slots)


def test_frames_subset_records_the_unnamed_channel_as_excluded():
    channels, excluded = resolve_channels(THREE_CHANNEL_INFO, n_pages=3)
    slots, slot_excluded = assign_slots(channels, excluded, [1, 0])
    manifest = build_channel_map(slots, slot_excluded, group_name(slots), "1,0")

    assert set(manifest["channels"]) == {"ch1", "ch2"}
    assert len(manifest["excluded_channels"]) == 1

    dropped = manifest["excluded_channels"][0]
    assert dropped["channel_name"] == "CH3"
    assert dropped["source_page"] == 2
    assert dropped["wavelength_nm"] == 561
    assert dropped["pmt_voltage_v"] == 700
    assert dropped["reason"] == "excluded by explicit --frames selection"


def test_frames_subset_of_one():
    channels, excluded = resolve_channels(THREE_CHANNEL_INFO, n_pages=3)
    slots, slot_excluded = assign_slots(channels, excluded, [2])

    assert [c["channel"] for c in slots] == ["CH3"]
    assert group_name(slots) == "561nm_3.2pct_700V"
    assert [c["channel"] for c in slot_excluded] == ["CH1", "CH2"]
    assert all(
        c["reason"] == "excluded by explicit --frames selection"
        for c in slot_excluded
    )


def test_full_arity_frames_still_valid_on_three_channels():
    """The pre-existing one-entry-per-channel form is unchanged."""
    channels, excluded = resolve_channels(THREE_CHANNEL_INFO, n_pages=3)
    slots, slot_excluded = assign_slots(channels, excluded, [2, 1, 0])

    assert [c["channel"] for c in slots] == ["CH3", "CH2", "CH1"]
    assert slot_excluded == []
    assert group_name(slots) == "561nm_3.2pct_700V_488nm_5.0pct_600V_405nm_2.0pct_500V"


def test_frames_subset_rejects_duplicates_and_out_of_range():
    channels, excluded = resolve_channels(THREE_CHANNEL_INFO, n_pages=3)
    with pytest.raises(ValueError, match="repeats"):
        assign_slots(channels, excluded, [0, 0])
    with pytest.raises(ValueError, match="repeats"):
        assign_slots(channels, excluded, [1, 2, 1])
    with pytest.raises(IndexError, match="page 3"):
        assign_slots(channels, excluded, [0, 3])
    with pytest.raises(ValueError, match="only 3 channels"):
        assign_slots(channels, excluded, [0, 1, 2, 0])


def test_frames_subset_end_to_end(tmp_path):
    """Three-channel file, --frames 1,0: two slots on disk, third recorded."""
    import json

    src = tmp_path / "raw"
    src.mkdir()
    stack = np.stack([
        np.full((16, 16), fill, dtype=np.uint16) for fill in (11, 22, 33)
    ])
    tifffile.imwrite(
        str(src / "three.tif"), stack, imagej=True,
        metadata={"Info": THREE_CHANNEL_INFO},
    )

    out = tmp_path / "prepared"
    _run(src, out, [1, 0], "1,0")

    group_dir = next(p for p in out.iterdir() if p.is_dir())
    assert group_dir.name == THREE_CHANNEL_SUBSET_NAME

    cell = group_dir / "cell1"
    assert sorted(p.name for p in cell.iterdir()) == ["ch1", "ch2"]

    manifest = json.loads((group_dir / "channel_map.json").read_text())
    assert manifest["group"] == THREE_CHANNEL_SUBSET_NAME

    for slot, info in manifest["channels"].items():
        written = tifffile.imread(str(next((cell / slot).glob("*.tif"))))
        assert np.array_equal(written, stack[info["source_page"]])

    assert manifest["channels"]["ch1"]["source_page"] == 1
    assert manifest["channels"]["ch2"]["source_page"] == 0
    dropped = manifest["excluded_channels"]
    assert [c["channel_name"] for c in dropped] == ["CH3"]
    assert dropped[0]["reason"] == "excluded by explicit --frames selection"


# ── End-to-end: the manifest must match the bytes on disk ──────────────

def _write_source_tiff(path):
    """A 3-page TIFF carrying the real FV3000 metadata, one constant per page."""
    stack = np.stack([
        np.full((16, 16), fill, dtype=np.uint16) for fill in (111, 222, 333)
    ])
    tifffile.imwrite(
        str(path), stack, imagej=True, metadata={"Info": REAL_INFO}
    )
    return stack


def _run(src, out, frame_order, frames_arg):
    with pytest.raises(SystemExit) as exc:
        process_images(str(src), str(out), frame_order, 1, False, frames_arg=frames_arg)
    assert exc.value.code == 0, "processing reported errors"


@pytest.mark.parametrize(
    "frame_order, frames_arg", [(None, None), ([1, 0], "1,0"), ([0, 2], "0,2")]
)
def test_written_pixels_match_the_manifest(tmp_path, frame_order, frames_arg):
    """
    Read each ch<i> TIFF back and assert it is byte-identical to the raw page
    the manifest names in that slot's source_page.

    This holds regardless of how resolution is refactored, and is the
    assertion that catches a manifest drifting away from the pixels.
    """
    import json

    src = tmp_path / "raw"
    src.mkdir()
    raw = _write_source_tiff(src / "sample.tif")

    out = tmp_path / "prepared"
    _run(src, out, frame_order, frames_arg)

    group_dir = next(p for p in out.iterdir() if p.is_dir())
    manifest = json.loads((group_dir / "channel_map.json").read_text())

    assert manifest["channels"], "manifest recorded no channels"
    for slot, info in manifest["channels"].items():
        written = tifffile.imread(
            str(next((group_dir / "cell1" / slot).glob("*.tif")))
        )
        expected = raw[info["source_page"]]
        assert np.array_equal(written, expected), (
            f"{slot} claims source_page {info['source_page']} but its pixels "
            f"do not match that raw page"
        )


def test_end_to_end_frames_override_slot_metadata(tmp_path):
    """The reported reproduction, checked against pixels rather than order."""
    import json

    src = tmp_path / "raw"
    src.mkdir()
    _write_source_tiff(src / "sample.tif")

    out = tmp_path / "prepared"
    _run(src, out, [1, 0], "1,0")

    group_dir = next(p for p in out.iterdir() if p.is_dir())
    assert group_dir.name == "561nm_15.3pct_876V_488nm_17.7pct_652V"
    manifest = json.loads((group_dir / "channel_map.json").read_text())

    assert manifest["group"] == group_dir.name
    assert manifest["channels"]["ch1"]["source_page"] == 1
    assert manifest["channels"]["ch1"]["wavelength_nm"] == 561
    assert manifest["channels"]["ch1"]["channel_name"] == "CH2"
    assert manifest["channels"]["ch2"]["source_page"] == 0
    assert manifest["channels"]["ch2"]["wavelength_nm"] == 488
    assert "NOT wavelength-ascending" in manifest["channel_order"]


def test_end_to_end_default_is_wavelength_ascending(tmp_path):
    import json

    src = tmp_path / "raw"
    src.mkdir()
    _write_source_tiff(src / "sample.tif")

    out = tmp_path / "prepared"
    _run(src, out, None, None)

    group_dir = next(p for p in out.iterdir() if p.is_dir())
    assert group_dir.name == "488nm_17.7pct_652V_561nm_15.3pct_876V"
    manifest = json.loads((group_dir / "channel_map.json").read_text())

    assert manifest["channel_order"] == (
        "excitation wavelength ascending (ch1 = lowest)"
    )
    assert manifest["channels"]["ch1"]["wavelength_nm"] == 488
    assert manifest["channels"]["ch1"]["source_page"] == 0
    assert manifest["channels"]["ch2"]["wavelength_nm"] == 561
    assert len(manifest["excluded_channels"]) == 1


# ── The bug itself: no positional pairing ──────────────────────────────

def _permute_detector_block(info, order):
    """Rewrite the detector block so its slots appear in a different order."""
    fields = ("Detector ID", "Detector linked channel ID", "Detector voltage")
    values = {
        field: dict(
            re.findall(rf"{field} #(\d+) = (\S+)", info)
        )
        for field in fields
    }
    lines = []
    for line in info.splitlines():
        matched = next(
            (f for f in fields if re.search(rf"{f} #\d+ = ", line)), None
        )
        if matched is None:
            lines.append(line)
            continue
        slot = re.search(r"#(\d+)", line).group(1)
        new_slot = order[int(slot) - 1]
        lines.append(f"- {matched} #{slot} = {values[matched][str(new_slot)]}")
    return "\n".join(lines) + "\n"


def test_detector_order_does_not_change_the_mapping():
    """
    The detector block's order is arbitrary relative to the channel block.

    Reversing it must not change a single resolved value -- if resolution
    ever regresses to zip/enumerate pairing, this is what catches it.
    """
    baseline, baseline_excluded = resolve_channels(REAL_INFO, n_pages=3)
    shuffled_info = _permute_detector_block(REAL_INFO, [4, 3, 2, 1])
    shuffled, shuffled_excluded = resolve_channels(shuffled_info, n_pages=3)

    assert shuffled == baseline
    assert shuffled_excluded == baseline_excluded


def test_laser_index_is_not_read_as_the_channel_number():
    """CH1 links laser index 1, which is LD488 (0-based) -- not LD405."""
    channels, _ = resolve_channels(REAL_INFO, n_pages=3)
    assert channels[0]["laser_id"] == "LD488"
    assert all(c["laser_id"] != "LD405" for c in channels)


def test_channel_count_ignores_laser_transmissivity():
    """
    All four lasers have transmissivity > 0 here; only 3 channels exist.

    Counting active lasers was the original over-count (4 channels for a
    3-page TIFF).
    """
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)
    assert len(channels) + len(excluded) == 3


# ── Failing loudly ─────────────────────────────────────────────────────

def test_sizec_page_count_mismatch_names_both_numbers():
    with pytest.raises(MetadataError) as exc:
        resolve_channels(REAL_INFO, n_pages=4)
    message = str(exc.value)
    assert "3" in message and "4" in message


def test_missing_channel_id_block_raises():
    stripped = "\n".join(
        line for line in REAL_INFO.splitlines() if "Channel CH" not in line
    )
    with pytest.raises(MetadataError):
        resolve_channels(stripped, n_pages=3)


def test_missing_sizec_raises():
    stripped = "\n".join(
        line for line in REAL_INFO.splitlines() if "SizeC" not in line
    )
    with pytest.raises(MetadataError):
        resolve_channels(stripped, n_pages=3)


def test_channel_without_a_detector_raises():
    stripped = "\n".join(
        line for line in REAL_INFO.splitlines()
        if "Detector linked channel ID #3" not in line
    )
    with pytest.raises(MetadataError, match="CH1"):
        resolve_channels(stripped, n_pages=3)


def test_empty_metadata_raises():
    with pytest.raises(MetadataError):
        resolve_channels("", n_pages=3)


# ── Manifest ───────────────────────────────────────────────────────────

def test_channel_map_records_the_audit_trail():
    channels, excluded = resolve_channels(REAL_INFO, n_pages=3)
    group = group_name(channels)
    manifest = build_channel_map(channels, excluded, group)

    assert manifest["group"] == group
    assert set(manifest["channels"]) == {"ch1", "ch2"}

    ch1 = manifest["channels"]["ch1"]
    assert ch1["source_page"] == 0
    assert ch1["channel_name"] == "CH1"
    assert ch1["wavelength_nm"] == 488
    assert ch1["transmissivity_pct"] == 17.7
    assert ch1["detector_id"] == "FV30-SD_D_1"
    assert ch1["pmt_voltage_v"] == 652
    assert ch1["dye"] == "EGFP"

    assert manifest["channels"]["ch2"]["source_page"] == 1
    assert manifest["channels"]["ch2"]["pmt_voltage_v"] == 876

    assert len(manifest["excluded_channels"]) == 1
    dropped = manifest["excluded_channels"][0]
    assert dropped["channel_name"] == "CH3"
    assert dropped["source_page"] == 2
    assert "LETD" in dropped["reason"]


# ── Non-regression for well-formed two-channel files ───────────────────

def test_plain_two_channel_file_is_unaffected():
    """No DIC channel, counts already lined up -> same structure as before."""
    info = """\
 SizeC = 2
- Channel CH1 ID = aaaa-1
- Channel CH1 linked laser index = 1
- Channel CH2 ID = bbbb-2
- Channel CH2 linked laser index = 2
- Detector ID #1 = __FV30-SD_D_1__
- Detector ID #2 = __FV30-SD_D_2__
- Detector linked channel ID #1 = aaaa-1
- Detector linked channel ID #2 = bbbb-2
- Detector voltage #1 = 580.0
- Detector voltage #2 = 500.0
- Laser LD488 data ID = LD488_1Imaging_main_phase_1
- Laser LD488 transmissivity = 5.0
- Laser LD488 wavelength = 488.0
- Laser LD561 data ID = LD561_1Imaging_main_phase_1
- Laser LD561 transmissivity = 3.2
- Laser LD561 wavelength = 561.0
laser name #1 = LD405
laser name #2 = LD488
laser name #3 = LD561
"""
    channels, excluded = resolve_channels(info, n_pages=2)
    assert not excluded
    assert [c["page"] for c in channels] == [0, 1]
    assert group_name(channels) == "488nm_5.0pct_580V_561nm_3.2pct_500V"


def test_lambda_phase_channel_is_excluded():
    """The lambda-phase exclusion rule survives the rewrite."""
    info = REAL_INFO.replace(
        "- Laser LD488 data ID = LD488_65794Imaging_main_phase_1",
        "- Laser LD488 data ID = LD488_65794Imaging_main_phase_Lambda",
    )
    channels, excluded = resolve_channels(info, n_pages=3)
    assert [c["channel"] for c in channels] == ["CH2"]
    assert any("lambda" in c["reason"].lower() for c in excluded)


def test_device_name_td_rule_applies_when_block_aligns():
    """
    deviceName is only trusted when it has one entry per acquired channel.

    In the real file the block is ragged (4 entries, 3 channels), so it is
    ignored there and LETD does the work. Here it aligns, so it applies.
    """
    info = """\
 SizeC = 2
- Channel CH1 ID = aaaa-1
- Channel CH1 linked laser index = 1
- Channel CH2 ID = bbbb-2
- Channel CH2 linked laser index = 2
- Detector ID #1 = __FV30-SD_D_1__
- Detector ID #2 = __FV30-SD_D_2__
- Detector linked channel ID #1 = aaaa-1
- Detector linked channel ID #2 = bbbb-2
- Detector voltage #1 = 580.0
- Detector voltage #2 = 500.0
- Laser LD488 transmissivity = 5.0
- Laser LD488 wavelength = 488.0
- Laser LD561 transmissivity = 3.2
- Laser LD561 wavelength = 561.0
channel deviceName #1 = SD1
channel deviceName #2 = TD
laser name #1 = LD405
laser name #2 = LD488
laser name #3 = LD561
"""
    channels, excluded = resolve_channels(info, n_pages=2)
    assert [c["channel"] for c in channels] == ["CH1"]
    assert excluded[0]["channel"] == "CH2"
    assert "TD" in excluded[0]["reason"]
