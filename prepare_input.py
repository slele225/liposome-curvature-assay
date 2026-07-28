"""
Split multi-frame TIFFs into per-channel folders for MATLAB detection.

Channel resolution is driven by the image data and by the explicit ID
links in the Olympus FV3000 metadata -- never by laser transmissivity
and never by positional pairing of the metadata blocks.

How a channel is resolved
-------------------------
1. The channel count comes from ``SizeC``, cross-checked against the
   actual page count of the TIFF. A disagreement is a hard error.
2. Each acquired channel ``CH<n>`` is resolved by ID:

       Channel CH<n> ID              = <uuid>
       Channel CH<n> linked laser index = <k>     (0-based into `laser name #<k+1>`)
       Detector linked channel ID #<m> = <uuid>   (find the m matching the channel)
       Detector voltage #<m>, Detector ID #<m>

   The detector block is NOT in the same order as the channel block, so
   the two are joined on the channel UUID. Detectors whose linked channel
   ID matches no acquired channel are configured-but-unused and are
   ignored entirely.
3. Non-fluorescence channels (transmitted light / DIC) and lambda-phase
   channels are dropped -- see the exclusion rules below.
4. Survivors are ordered by excitation wavelength ascending, so ``ch1``
   is always the lowest-wavelength channel.

Folder name format, in ch1, ch2, ... slot order (so it tracks --frames):
    {wave1}nm_{trans1}pct_{volt1}V_{wave2}nm_{trans2}pct_{volt2}V_...

Each output group folder also gets a ``channel_map.json`` recording the
full resolved mapping and every excluded channel with its reason, so the
ch1/ch2 assignment stays auditable without reopening the raw files.
"""

import os
import re
import sys
import json
import argparse
import tifffile


class MetadataError(ValueError):
    """FV3000 metadata is missing, malformed, or cannot be resolved unambiguously."""


# ── Exclusion rules ────────────────────────────────────────────────────
#
# Named rules rather than magic strings buried in conditionals. Each maps
# to a human-readable reason recorded in channel_map.json.

# Olympus transmitted-light detector (e.g. "FV31-LETD_D_1"). This is the
# DIC / brightfield channel: it carries no fluorescence signal and must
# never become a ch<i> folder.
TRANSMITTED_LIGHT_DETECTOR_TOKEN = "LETD"

# The channel-configuration block labels the transmitted-light path "TD"
# (as opposed to "SD1"/"SD2" for the spectral fluorescence detectors).
# See _device_names_by_channel for why this is only conditionally usable.
TRANSMITTED_LIGHT_DEVICE_NAME = "TD"

# Spectral (lambda) scan phases are configuration entries, not acquired
# imaging channels. Only "phase_1" main imaging phases are kept.
LAMBDA_PHASE_TOKEN = "Lambda"


# ── Metadata field patterns ────────────────────────────────────────────
#
# The FV3000 Info block mixes two line styles: dash-prefixed summary
# lines ("- Detector ID #1 = ...") and bare table lines ("laser name #1 =
# ..."). The leading "-" is therefore optional throughout.

SIZE_C_RE = re.compile(r"^\s*SizeC\s*=\s*(\d+)\s*$", re.M)

CHANNEL_ID_RE = re.compile(
    r"^\s*-?\s*Channel\s+CH(\d+)\s+ID\s*=\s*(\S+)\s*$", re.M
)
CHANNEL_LASER_INDEX_RE = re.compile(
    r"^\s*-?\s*Channel\s+CH(\d+)\s+linked\s+laser\s+index\s*=\s*(\d+)\s*$", re.M
)

LASER_NAME_RE = re.compile(r"^\s*laser\s+name\s+#0*(\d+)\s*=\s*(\S+)\s*$", re.M)
LASER_WAVELENGTH_RE = re.compile(
    r"^\s*-?\s*Laser\s+(\S+)\s+wavelength\s*=\s*([0-9.]+)\s*$", re.M
)
LASER_TRANS_RE = re.compile(
    r"^\s*-?\s*Laser\s+(\S+)\s+transmissivity\s*=\s*([0-9.]+)\s*$", re.M
)
LASER_DATA_ID_RE = re.compile(
    r"^\s*-?\s*Laser\s+(\S+)\s+data\s+ID\s*=\s*(\S+)\s*$", re.M
)

DETECTOR_LINKED_CHANNEL_RE = re.compile(
    r"^\s*-?\s*Detector\s+linked\s+channel\s+ID\s+#0*(\d+)\s*=\s*(\S+)\s*$", re.M
)
DETECTOR_ID_RE = re.compile(
    r"^\s*-?\s*Detector\s+ID\s+#0*(\d+)\s*=\s*(\S+)\s*$", re.M
)
DETECTOR_VOLTAGE_RE = re.compile(
    r"^\s*-?\s*Detector\s+voltage\s+#0*(\d+)\s*=\s*([0-9.]+)\s*$", re.M
)

DEVICE_NAME_RE = re.compile(
    r"^\s*channel\s+deviceName\s+#0*(\d+)\s*=\s*(\S+)\s*$", re.M
)

DYE_NAME_RE = re.compile(r"^\s*dyeData\s+name\s+#0*(\d+)\s*=\s*(.+?)\s*$", re.M)
DYE_EXCITATION_RE = re.compile(
    r"^\s*dyeData\s+excitationWavelength\s+#0*(\d+)\s*=\s*([0-9.]+)\s*$", re.M
)
CHANNEL_DYE_NAME_RE = re.compile(
    r"^\s*channel\s+dyeName\s+#0*(\d+)\s*=\s*(.+?)\s*$", re.M
)

LASER_DIGITS_RE = re.compile(r"(\d+)")


def _indexed(pattern, info, cast=str):
    """Parse a '#<i> = <value>' block into {index: value}."""
    return {int(i): cast(v) for i, v in pattern.findall(info)}


def _named(pattern, info, cast=str):
    """Parse a '<name> = <value>' block into {name: value}."""
    return {n: cast(v) for n, v in pattern.findall(info)}


def _strip_detector_id(raw):
    """'__FV30-SD_D_1__' -> 'FV30-SD_D_1'."""
    return raw.strip("_")


def _device_names_by_channel(info, size_c):
    """
    Return {channel_number: deviceName}, or {} if it cannot be trusted.

    The ``channel <field> #<i>`` sub-blocks in the FV3000 Info dump are
    ragged -- on the bundled files ``channel name`` has 7 entries,
    ``channel laserDataId`` has 10 and ``channel deviceName`` has 4, all
    describing different slices of the instrument configuration. They are
    therefore NOT co-indexed with each other or with CH1..CH<SizeC>, and
    reading deviceName #<n> as "the device for CH<n>" silently mislabels
    channels.

    We only use this block when it has exactly one entry per acquired
    channel, which is the one case where alignment with CH1..CH<SizeC> is
    defensible. Otherwise the transmitted-light detector is caught by
    TRANSMITTED_LIGHT_DETECTOR_TOKEN, which is joined by UUID and always
    reliable.
    """
    devices = _indexed(DEVICE_NAME_RE, info)
    if sorted(devices) != list(range(1, size_c + 1)):
        return {}
    return devices


def _dyes_by_excitation(info):
    """
    Return {excitation_wavelength: dye_name} for unambiguous exact matches.

    ``dyeData`` entries are matched to a channel only when the dye's
    excitation wavelength equals the channel's laser wavelength exactly,
    and only when that wavelength maps to a single dye. Anything fuzzier
    (nearest-wavelength, positional) risks mislabelling, and a wrong dye
    label is worse than none.
    """
    names = _indexed(DYE_NAME_RE, info)
    excitations = _indexed(DYE_EXCITATION_RE, info, lambda v: int(float(v)))

    by_wavelength = {}
    for idx, name in names.items():
        wavelength = excitations.get(idx)
        if wavelength is None:
            continue
        by_wavelength.setdefault(wavelength, set()).add(name)

    return {
        wavelength: next(iter(dye_names))
        for wavelength, dye_names in by_wavelength.items()
        if len(dye_names) == 1
    }


# How a channel's dye label was established, recorded in channel_map.json so
# the confidence behind each label is visible.
DYE_SOURCE_EXCITATION = "exact dyeData excitationWavelength match"
DYE_SOURCE_CHANNEL_BLOCK = (
    "channel dyeName block, cross-checked against dyeData excitation wavelengths"
)


def _assign_dyes(info, fluorescence_channels):
    """
    Return {channel_number: (dye_name_or_None, source_or_None)}.

    Two independent sources, neither sufficient alone:

    * ``dyeData excitationWavelength #<i>`` matched exactly to the channel's
      laser wavelength. Unambiguous, but incomplete -- Texas Red records
      595 nm, so it never matches the 561 nm channel exciting it.
    * the ``channel dyeName #<i>`` block. Complete, but one of the ragged
      ``channel <field> #<i>`` blocks carrying no UUIDs, so it can only be
      read positionally.

    The dyeName block is used only when it has exactly one entry per
    fluorescence channel (the same gating rule as deviceName -- dyes exist
    for fluorescence channels, not for the transmitted-light one) AND every
    entry that can be independently cross-checked against an exact
    excitation match agrees with it. A single disagreement discards the
    whole block: a wrong dye label is worse than a missing one.
    """
    by_excitation = _dyes_by_excitation(info)
    ordered = sorted(fluorescence_channels, key=lambda c: c["channel_number"])
    exact = {
        c["channel_number"]: by_excitation.get(c["wavelength"]) for c in ordered
    }

    block = _indexed(CHANNEL_DYE_NAME_RE, info)
    positional = {}
    if ordered and sorted(block) == list(range(1, len(ordered) + 1)):
        candidate = {
            c["channel_number"]: block[i]
            for i, c in enumerate(ordered, start=1)
        }
        if all(
            exact[n] is None or exact[n] == candidate[n] for n in candidate
        ):
            positional = candidate

    assignments = {}
    for c in ordered:
        n = c["channel_number"]
        if positional:
            assignments[n] = (positional[n], DYE_SOURCE_CHANNEL_BLOCK)
        elif exact[n] is not None:
            assignments[n] = (exact[n], DYE_SOURCE_EXCITATION)
        else:
            assignments[n] = (None, None)
    return assignments


def _laser_wavelength(laser_id, explicit_wavelengths):
    """
    Excitation wavelength for a laser, preferring the explicit metadata field.

    Not every laser gets a ``Laser <id> wavelength`` line (the bundled
    files omit it for LD405 and LD640), so fall back to the digits in the
    laser name.
    """
    if laser_id in explicit_wavelengths:
        return int(float(explicit_wavelengths[laser_id]))
    digits = LASER_DIGITS_RE.search(laser_id)
    if not digits:
        raise MetadataError(
            f"cannot determine excitation wavelength for laser '{laser_id}': "
            f"no 'Laser {laser_id} wavelength' entry and no digits in its name"
        )
    return int(digits.group(1))


def resolve_channels(info, n_pages=None):
    """
    Resolve acquired channels from an FV3000 ``Info`` block.

    Returns ``(channels, excluded)``. ``channels`` are the surviving
    fluorescence channels ordered by excitation wavelength ascending;
    ``excluded`` records every dropped channel with a reason.

    Raises MetadataError on anything ambiguous or unparseable -- a
    wrong-but-silent channel assignment produces plausible-looking output
    and quietly corrupts the analysis, so this fails loudly instead.
    """
    if not info:
        raise MetadataError("no ImageJ 'Info' metadata block found")

    # 1. Channel count from the data, cross-checked against the pages.
    size_c_match = SIZE_C_RE.search(info)
    if not size_c_match:
        raise MetadataError("metadata has no 'SizeC' entry")
    size_c = int(size_c_match.group(1))
    if size_c < 1:
        raise MetadataError(f"metadata reports SizeC = {size_c}")

    if n_pages is not None and n_pages != size_c:
        raise MetadataError(
            f"channel count mismatch: metadata reports SizeC = {size_c} but "
            f"the TIFF has {n_pages} pages"
        )

    # 2. Acquired channels, by ID.
    channel_ids = {int(n): uuid for n, uuid in CHANNEL_ID_RE.findall(info)}
    if not channel_ids:
        raise MetadataError(
            "metadata has no 'Channel CH<n> ID' entries; channels cannot be "
            "resolved by ID"
        )
    expected = list(range(1, size_c + 1))
    if sorted(channel_ids) != expected:
        raise MetadataError(
            f"expected channels CH1..CH{size_c} (from SizeC), found "
            f"{sorted('CH%d' % n for n in channel_ids)}"
        )

    laser_indices = {
        int(n): int(k) for n, k in CHANNEL_LASER_INDEX_RE.findall(info)
    }
    laser_names = _indexed(LASER_NAME_RE, info)
    if not laser_names:
        raise MetadataError("metadata has no 'laser name #<i>' entries")

    explicit_wavelengths = _named(LASER_WAVELENGTH_RE, info)
    transmissivities = _named(LASER_TRANS_RE, info, float)
    laser_data_ids = _named(LASER_DATA_ID_RE, info)

    # 3. Detector block, joined to channels on the channel UUID. Detectors
    #    linked to a UUID that is not an acquired channel are unused
    #    configuration and are dropped here.
    detector_ids = _indexed(DETECTOR_ID_RE, info)
    detector_voltages = _indexed(DETECTOR_VOLTAGE_RE, info, lambda v: int(float(v)))
    detector_by_channel_uuid = {}
    for m, uuid in DETECTOR_LINKED_CHANNEL_RE.findall(info):
        m = int(m)
        if uuid not in channel_ids.values():
            continue
        if uuid in detector_by_channel_uuid:
            raise MetadataError(
                f"channel ID {uuid} is linked by more than one detector "
                f"(#{detector_by_channel_uuid[uuid]} and #{m})"
            )
        detector_by_channel_uuid[uuid] = m

    device_names = _device_names_by_channel(info, size_c)

    resolved, excluded = [], []

    for number in expected:
        name = f"CH{number}"
        uuid = channel_ids[number]

        # -- excitation laser, by index into the laser table (0-based) --
        if number not in laser_indices:
            raise MetadataError(f"{name} has no 'linked laser index' entry")
        laser_index = laser_indices[number]
        laser_slot = laser_index + 1
        if laser_slot not in laser_names:
            raise MetadataError(
                f"{name} links laser index {laser_index} (slot #{laser_slot}) "
                f"but only slots {sorted(laser_names)} exist"
            )
        laser_id = laser_names[laser_slot]

        wavelength = _laser_wavelength(laser_id, explicit_wavelengths)
        if laser_id not in transmissivities:
            raise MetadataError(
                f"{name} uses laser {laser_id} but no "
                f"'Laser {laser_id} transmissivity' entry was found"
            )
        transmissivity = transmissivities[laser_id]

        # -- detector, by UUID join --
        detector_slot = detector_by_channel_uuid.get(uuid)
        if detector_slot is None:
            raise MetadataError(
                f"{name} (ID {uuid}) is not linked by any detector"
            )
        if detector_slot not in detector_ids:
            raise MetadataError(f"no 'Detector ID #{detector_slot}' entry")
        if detector_slot not in detector_voltages:
            raise MetadataError(f"no 'Detector voltage #{detector_slot}' entry")
        detector_id = _strip_detector_id(detector_ids[detector_slot])
        voltage = detector_voltages[detector_slot]

        record = {
            # DimensionOrder is XYCZT with SizeZ = SizeT = 1, so page order
            # follows channel order: CH1 -> page 0, CH2 -> page 1, ...
            "page": number - 1,
            "channel": name,
            "channel_number": number,
            "laser_id": laser_id,
            "wavelength": wavelength,
            "transmissivity": transmissivity,
            "detector_id": detector_id,
            "voltage": voltage,
            "dye": None,
            "dye_source": None,
        }

        reason = _exclusion_reason(record, device_names, laser_data_ids)
        if reason:
            excluded.append(dict(record, reason=reason))
        else:
            resolved.append(record)

    # 4. Dyes, once the fluorescence channels are known -- the dyeName block
    #    is gated on their count, so this cannot run inside the loop.
    for number, (dye, source) in _assign_dyes(info, resolved).items():
        for record in resolved:
            if record["channel_number"] == number:
                record["dye"], record["dye_source"] = dye, source

    # 5. Wavelength ascending; channel number breaks ties so the order is
    #    deterministic when two channels share an excitation laser.
    resolved.sort(key=lambda c: (c["wavelength"], c["channel_number"]))
    return resolved, excluded


def _exclusion_reason(record, device_names, laser_data_ids):
    """Return why this channel is not a fluorescence channel, or None to keep it."""
    if TRANSMITTED_LIGHT_DETECTOR_TOKEN in record["detector_id"]:
        return (
            f"transmitted-light/DIC channel: detector "
            f"{record['detector_id']} is an "
            f"{TRANSMITTED_LIGHT_DETECTOR_TOKEN} detector"
        )

    device = device_names.get(record["channel_number"])
    if device == TRANSMITTED_LIGHT_DEVICE_NAME:
        return (
            f"transmitted-light/DIC channel: channel deviceName is "
            f"'{TRANSMITTED_LIGHT_DEVICE_NAME}'"
        )

    data_id = laser_data_ids.get(record["laser_id"], "")
    if LAMBDA_PHASE_TOKEN in data_id:
        return f"lambda-phase channel: laserDataId '{data_id}'"

    return None


def parse_metadata(tiff_path):
    """Resolve the channels of one TIFF. See resolve_channels."""
    with tifffile.TiffFile(tiff_path) as tif:
        info = (tif.imagej_metadata or {}).get("Info", "")
        n_pages = len(tif.pages)
    return resolve_channels(info, n_pages=n_pages)


# ── Presentation ───────────────────────────────────────────────────────

def group_name(slots):
    """
    Folder name like 488nm_17.7pct_652V_561nm_15.3pct_876V.

    Built from the slot list, so the name reads in ch1, ch2, ... order and
    changes when --frames reorders the slots. Channels dropped by the
    exclusion rules are not slots and so never appear.
    """
    parts = [
        f"{c['wavelength']}nm_{c['transmissivity']:.1f}pct_{c['voltage']}V"
        for c in slots
    ]
    return "_".join(parts)


def channels_summary(channels):
    """Summary like '[488 (17.7%, 652V), 561 (15.3%, 876V)]'."""
    parts = [
        f"{c['wavelength']} ({c['transmissivity']:.1f}%, {c['voltage']}V)"
        for c in channels
    ]
    return "[" + ", ".join(parts) + "]"


def format_channel_table(slots, excluded, indent="      "):
    """Human-readable slot table plus exclusions with reasons."""
    lines = []
    for i, c in enumerate(slots, start=1):
        dye = f"  {c['dye']}" if c["dye"] else ""
        forced = "  [rule-excluded, forced by --frames]" if "reason" in c else ""
        lines.append(
            f"{indent}ch{i}  page {c['page']}  {c['channel']}  "
            f"{c['wavelength']}nm  {c['transmissivity']:.1f}%  "
            f"{c['voltage']}V  {c['detector_id']}{dye}{forced}"
        )
    written_pages = {c["page"] for c in slots}
    dropped = [c for c in excluded if c["page"] not in written_pages]
    if dropped:
        lines.append(f"{indent}excluded:")
        for c in dropped:
            lines.append(
                f"{indent}  {c['channel']} (page {c['page']}, "
                f"{c['wavelength']}nm, {c['detector_id']}): {c['reason']}"
            )
    return lines


def override_warnings(slots, frames_arg, group):
    """
    Warnings for a --frames override that departs from documented behaviour.

    ch1-is-lowest-wavelength is a contract the README and the downstream
    steps rely on. Overriding it is legitimate, but must never be silent.
    """
    if frames_arg is None:
        return []

    warnings = []
    if not is_wavelength_ascending(slots):
        order = ", ".join(
            f"ch{i}={c['wavelength']}nm" for i, c in enumerate(slots, start=1)
        )
        warnings.append(
            f"WARNING: --frames {frames_arg} puts group '{group}' in a "
            f"non-wavelength-ascending order ({order}). The README and the "
            f"downstream steps document ch1 as the lowest wavelength; check "
            f"channel_map.json before setting --lipid-col/--protein-col."
        )
    for i, c in enumerate(slots, start=1):
        if "reason" in c:
            warnings.append(
                f"WARNING: --frames {frames_arg} writes page {c['page']} "
                f"({c['channel']}) to ch{i}, but that channel was excluded as "
                f"non-fluorescence: {c['reason']}"
            )
    return warnings


ASCENDING_ORDER = "excitation wavelength ascending (ch1 = lowest)"

# The group folder name is built from the same slot list as the manifest, so
# it reads in ch1, ch2, ... order and tracks --frames. Stated in the manifest
# so the name's meaning is explicit wherever it is read.
GROUP_NAME_MEANING = (
    "wavelength/transmissivity/PMT voltage of the channels written to each "
    "slot, in ch1, ch2, ... order. Reflects --frames, so it is only "
    "wavelength-ascending when the slot order is -- see `channel_order`. "
    "Channels dropped by the exclusion rules do not appear."
)


# A resolved fluorescence channel that --frames simply did not name. Not a
# metadata problem -- the operator chose a subset -- but it is recorded so the
# manifest accounts for every channel in the file.
FRAMES_DESELECTED_REASON = "excluded by explicit --frames selection"


def assign_slots(channels, excluded, frame_order):
    """
    Partition the resolved channels into ``(slots, excluded)``.

    ``slots`` is the ordered ``[ch1, ch2, ...]`` assignment and is the single
    source of truth: the written TIFFs, the group folder name and
    channel_map.json all come from it, so a slot's metadata can never drift
    from the page actually written to it.

    Without ``--frames`` the slots are the resolved fluorescence channels in
    wavelength-ascending order. With ``--frames`` each slot takes the
    requested page *together with the metadata of the channel that page
    belongs to*, which is what keeps the manifest honest under an override.

    ``--frames`` may name fewer pages than there are resolved channels,
    meaning "keep only these, in this order". The channels left unnamed move
    into the returned ``excluded`` list, so they stay out of both the slots
    and the group folder name.
    """
    if frame_order is None:
        return list(channels), list(excluded)

    if not frame_order:
        raise ValueError("--frames must name at least one page")

    if len(frame_order) > len(channels):
        raise ValueError(
            f"--frames has {len(frame_order)} entries but this TIFF resolved "
            f"only {len(channels)} channels"
        )

    repeated = sorted({p for p in frame_order if frame_order.count(p) > 1})
    if repeated:
        raise ValueError(
            f"--frames repeats page(s) {repeated}; each page may be used at "
            f"most once"
        )

    by_page = {c["page"]: c for c in list(channels) + list(excluded)}
    slots = []
    for page in frame_order:
        channel = by_page.get(page)
        if channel is None:
            raise IndexError(
                f"--frames requests page {page}, which is not a channel page "
                f"(valid pages: 0..{len(by_page) - 1})"
            )
        slots.append(channel)

    written = {c["page"] for c in slots}
    remaining = [
        dict(c, reason=FRAMES_DESELECTED_REASON)
        for c in channels
        if c["page"] not in written
    ]
    remaining += [c for c in excluded if c["page"] not in written]
    remaining.sort(key=lambda c: c["page"])
    return slots, remaining


def is_wavelength_ascending(slots):
    wavelengths = [c["wavelength"] for c in slots]
    return wavelengths == sorted(wavelengths)


def channel_order_description(slots, frames_arg):
    """Describe the actual slot order -- including when --frames overrode it."""
    if frames_arg is None:
        return ASCENDING_ORDER
    if is_wavelength_ascending(slots):
        return (
            f"custom order via --frames {frames_arg} "
            f"(coincides with wavelength-ascending)"
        )
    return f"custom order via --frames {frames_arg} (NOT wavelength-ascending)"


def build_channel_map(slots, excluded, group, frames_arg=None):
    """
    The channel_map.json payload for one output group.

    ``slots`` must be the same ordered list used to write the TIFFs. Slot
    order is never re-derived here.
    """
    written_pages = {c["page"] for c in slots}

    def slot_entry(c):
        entry = {
            "source_page": c["page"],
            "channel_name": c["channel"],
            "laser_id": c["laser_id"],
            "wavelength_nm": c["wavelength"],
            "transmissivity_pct": c["transmissivity"],
            "detector_id": c["detector_id"],
            "pmt_voltage_v": c["voltage"],
            "dye": c["dye"],
            "dye_source": c["dye_source"],
        }
        # A rule-excluded channel forced into a slot by --frames keeps its
        # reason, so the manifest never presents DIC data as fluorescence.
        if "reason" in c:
            entry["excluded_by_rule"] = c["reason"]
        return entry

    return {
        "group": group,
        "group_name_describes": GROUP_NAME_MEANING,
        "channel_order": channel_order_description(slots, frames_arg),
        "frames_argument": frames_arg,
        "channels": {
            f"ch{i}": slot_entry(c) for i, c in enumerate(slots, start=1)
        },
        "excluded_channels": [
            {
                "source_page": c["page"],
                "channel_name": c["channel"],
                "laser_id": c["laser_id"],
                "wavelength_nm": c["wavelength"],
                "detector_id": c["detector_id"],
                "pmt_voltage_v": c["voltage"],
                "reason": c["reason"],
            }
            for c in excluded
            if c["page"] not in written_pages
        ],
    }


# ── Frames argument ────────────────────────────────────────────────────

def parse_int_list(s):
    if not s.strip():
        return []
    parts = [p.strip() for p in s.split(",") if p.strip()]
    out = []
    for p in parts:
        if p.startswith("-") or not p.isdigit():
            raise ValueError(f"Invalid index: '{p}'")
        out.append(int(p))
    return out


# ── Image processing ───────────────────────────────────────────────────

def center_crop_2d(img2d, divisor):
    if divisor <= 1:
        return img2d
    h, w = img2d.shape[:2]
    ch, cw = h // divisor, w // divisor
    y0, x0 = (h - ch) // 2, (w - cw) // 2
    return img2d[y0 : y0 + ch, x0 : x0 + cw]


def process_images(
    input_folder, output_base, frame_order, crop_divisor, dry_run, frames_arg=None
):
    tiff_files = sorted(
        f for f in os.listdir(input_folder) if f.lower().endswith((".tif", ".tiff"))
    )
    if not tiff_files:
        print("Error: No .tif/.tiff files found.")
        sys.exit(1)

    group_stats = {}
    group_maps = {}
    warned_groups = set()
    processed = 0
    skipped = 0
    errors = 0

    print(f"\nFound {len(tiff_files)} TIFF files to process...")
    if dry_run:
        print("(DRY RUN — no files will be written)")
    print("=" * 60)

    for tiff_file in tiff_files:
        input_path = os.path.join(input_folder, tiff_file)
        base_name = os.path.splitext(tiff_file)[0]

        try:
            channels, excluded = parse_metadata(input_path)
            if not channels:
                print(
                    f"  ⚠ {tiff_file}: no fluorescence channels remain after "
                    f"exclusions, skipping"
                )
                skipped += 1
                continue

            # One source of truth: slots drive the written pixels, the folder
            # name and the manifest. slot_excluded accounts for every channel
            # that did not get a slot, whether by rule or by --frames.
            slots, slot_excluded = assign_slots(channels, excluded, frame_order)
            fo = [c["page"] for c in slots]

            # Built from the same slots as the manifest, so the folder name
            # reads in ch1, ch2, ... order.
            v_group = group_name(slots)
            channel_map = build_channel_map(
                slots, slot_excluded, v_group, frames_arg
            )

            # Every cell in a group must share one mapping, otherwise the
            # single channel_map.json would misdescribe some of them.
            if v_group in group_maps and group_maps[v_group] != channel_map:
                raise MetadataError(
                    f"channel mapping differs from earlier files in group "
                    f"'{v_group}'; refusing to write an inconsistent "
                    f"channel_map.json"
                )
            group_maps[v_group] = channel_map

            if v_group not in warned_groups:
                warned_groups.add(v_group)
                for warning in override_warnings(slots, frames_arg, v_group):
                    print(warning, file=sys.stderr)

            group_stats[v_group] = group_stats.get(v_group, 0) + 1
            cell_idx = group_stats[v_group]

            print(f"  {tiff_file}: {len(slots)} channels -> {v_group}/cell{cell_idx}/")
            for line in format_channel_table(slots, slot_excluded):
                print(line)

            if dry_run:
                processed += 1
                continue

            cell_dir = os.path.join(output_base, v_group, f"cell{cell_idx}")
            os.makedirs(cell_dir, exist_ok=True)

            with tifffile.TiffFile(input_path) as tif:
                n_pages = len(tif.pages)
                bad_frames = [fi for fi in fo if fi < 0 or fi >= n_pages]
                if bad_frames:
                    raise IndexError(
                        f"Frames out of range: {bad_frames} (n_pages={n_pages})"
                    )

                for save_idx, frame_idx in enumerate(fo):
                    channel_num = save_idx + 1
                    channel_dir = os.path.join(cell_dir, f"ch{channel_num}")
                    os.makedirs(channel_dir, exist_ok=True)

                    frame = tif.pages[frame_idx].asarray()
                    img_to_save = center_crop_2d(frame, crop_divisor)

                    file_out = f"{base_name}_ch{channel_num}.tif"
                    tifffile.imwrite(
                        os.path.join(channel_dir, file_out), img_to_save
                    )

            map_path = os.path.join(output_base, v_group, "channel_map.json")
            with open(map_path, "w", encoding="utf-8") as fh:
                json.dump(channel_map, fh, indent=2, ensure_ascii=False)
                fh.write("\n")

            processed += 1

        except Exception as e:
            print(f"  ✗ Error processing {tiff_file}: {e}")
            errors += 1

    print("=" * 60)
    if dry_run:
        print("\nDRY RUN COMPLETE — no files were written.")
    else:
        print("\nPROCESSING COMPLETE!")
    print(f"  ✓ Processed: {processed} files")
    if skipped:
        print(f"  ⚠ Skipped (no fluorescence channels): {skipped} files")
    print(f"  ✗ Errors: {errors} files")
    if not dry_run:
        print(f"  Output saved to: {output_base}")

    if group_stats:
        print("\nGroups found:")
        for group, count in sorted(group_stats.items()):
            print(f"  • {group}: {count} cells")

    sys.exit(1 if errors > 0 else 0)


# ── CLI ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Split TIFF channels by laser/detector metadata "
                    "and prepare MATLAB input."
    )
    parser.add_argument("--input", required=True, help="Input folder containing TIFFs")
    parser.add_argument("--output", required=True, help="Base output folder")
    parser.add_argument(
        "--frames",
        default=None,
        help="Comma-separated source page indices naming the pages to keep, "
             "in ch1, ch2, ... order (e.g. '0,1'). May name fewer pages than "
             "there are resolved channels, in which case the unnamed ones are "
             "dropped. Indices must be unique and in range. Optional — "
             "channel resolution already yields the correct pages, so this is "
             "rarely needed.",
    )
    parser.add_argument(
        "--crop", type=int, default=1, help="Center crop divisor (1 = no crop)"
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Resolve channels and print the per-channel table plus any "
             "excluded channels, without writing files.",
    )
    args = parser.parse_args()

    # The summary uses ✓/✗/⚠, which the Windows cp1252 default kills as soon
    # as stdout is a pipe or file rather than a console.
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")

    if not os.path.isdir(args.input):
        print(f"Error: Input folder does not exist: {args.input}")
        sys.exit(1)

    if args.crop < 1:
        print("Error: --crop must be >= 1")
        sys.exit(1)

    if args.frames is None:
        frame_order = None
    else:
        try:
            frame_order = parse_int_list(args.frames)
            if not frame_order:
                print("Error: --frames must contain at least one index.")
                sys.exit(1)
            # File-independent, so report it once here rather than per TIFF.
            repeated = sorted(
                {p for p in frame_order if frame_order.count(p) > 1}
            )
            if repeated:
                print(
                    f"Error: --frames repeats page(s) "
                    f"{', '.join(str(p) for p in repeated)}; each page may be "
                    f"used at most once."
                )
                sys.exit(1)
        except ValueError as e:
            print(f"Error parsing --frames: {e}")
            sys.exit(1)

    if not args.dry_run:
        os.makedirs(args.output, exist_ok=True)

    print("=" * 60)
    print("IMAGE PROCESSOR — Generic Multi-Channel Splitter")
    print("=" * 60)
    print(f"Input folder:  {args.input}")
    print(f"Output folder: {args.output}")
    print(f"Frame order:   {frame_order if frame_order is not None else '(auto: resolved from metadata)'}")
    print(f"Crop divisor:  {args.crop}")
    if args.dry_run:
        print("Mode:          DRY RUN")

    process_images(
        args.input, args.output, frame_order, args.crop, args.dry_run,
        frames_arg=args.frames,
    )


if __name__ == "__main__":
    main()
