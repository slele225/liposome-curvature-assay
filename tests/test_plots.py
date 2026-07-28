"""
Tests for core.plots.

Every builder must return a bare matplotlib Figure -- never a pyplot-managed
one -- so the same object can be embedded in Qt or saved by the caller.
"""

import numpy as np
import pytest
from matplotlib.figure import Figure

from core.calibration import calibrate, load_dls_number, run_bootstrap
from core.curvature import curvature_curve, normalize
from core.plots import (
    apply_transform,
    figure_amplitude_histograms,
    figure_bootstrap,
    figure_curvature,
    figure_diameter_distribution,
    figure_dls_comparison,
    figure_dls_distribution,
    figure_dls_overlay,
    figure_overlay,
    figure_protein_amplitude,
    figure_threshold_preview,
    save_figure,
    zoom_limits_raw,
)


@pytest.fixture(scope="module")
def puncta():
    rng = np.random.default_rng(3)
    lipid = rng.lognormal(mean=3.0, sigma=0.4, size=800)
    protein = lipid * rng.uniform(0.5, 1.5, size=800)
    return lipid, protein


@pytest.fixture(scope="module")
def curve(puncta):
    lipid, protein = puncta
    return curvature_curve(lipid, protein, bin_width=2.0, conversion_factor=2.5)


@pytest.fixture(scope="module")
def calibration(sample_dls):
    d, w = load_dls_number(sample_dls)
    rng = np.random.default_rng(11)
    sqrt_A = rng.choice(d, size=3000, p=w / w.sum()) / 2.5
    return calibrate(d, w, sqrt_A)


def _assert_bare_figure(fig):
    """A pyplot figure carries a number; a bare Figure does not."""
    assert isinstance(fig, Figure)
    assert getattr(fig, "number", None) is None
    assert fig.axes, "figure should have at least one axes"


# ── Helpers ────────────────────────────────────────────────────────────

def test_apply_transform():
    vals = np.array([4.0, 16.0])
    raw, lbl = apply_transform(vals, "raw")
    np.testing.assert_allclose(raw, vals)
    assert lbl == "A"

    sq, lbl = apply_transform(vals, "sqrt")
    np.testing.assert_allclose(sq, [2.0, 4.0])
    assert lbl == "√A"

    ls, lbl = apply_transform(vals, "log_sqrt")
    np.testing.assert_allclose(ls, np.log([2.0, 4.0]))
    assert lbl == "log(√A)"


def test_apply_transform_handles_nonpositive():
    out, _ = apply_transform(np.array([0.0, -1.0]), "log_sqrt")
    assert np.all(np.isfinite(out))


def test_zoom_limits_raw_full_range():
    vals = np.arange(100.0)
    lo, hi = zoom_limits_raw(vals, 100)
    assert lo < 0 and hi > 99


# ── Figure builders ────────────────────────────────────────────────────

def test_figure_curvature(curve):
    fig = figure_curvature(curve, title="t", bin_width=2.0)
    _assert_bare_figure(fig)
    ax = fig.axes[0]
    assert ax.get_xlabel() == "Liposome diameter (nm)"
    assert len(ax.collections) == 2  # scatter of puncta + scatter of bin means


def test_figure_curvature_ypad_tightens_axis(curve):
    wide = figure_curvature(curve, "t", 2.0, y_pad=1.0).axes[0].get_ylim()
    tight = figure_curvature(curve, "t", 2.0, y_pad=0.1).axes[0].get_ylim()
    assert (tight[1] - tight[0]) < (wide[1] - wide[0])


def test_figure_overlay(curve):
    norm = normalize(curve.bin_centres, curve.bin_means, "rightmost")
    fig = figure_overlay(
        [(curve.bin_centres, norm, "cond A", 2.5, curve.n_points)],
        normalize_to="rightmost",
    )
    _assert_bare_figure(fig)
    assert "fold enrichment vs flat" in fig.axes[0].get_ylabel()


def test_figure_amplitude_histograms_one_and_two_channels(puncta):
    lipid, protein = puncta
    one = figure_amplitude_histograms({"A_ch1": lipid})
    _assert_bare_figure(one)
    assert len(one.axes) == 1

    two = figure_amplitude_histograms({"A_ch1": lipid, "A_ch2": protein})
    _assert_bare_figure(two)
    assert len(two.axes) == 2


@pytest.mark.parametrize("transform", ["raw", "sqrt", "log_sqrt"])
def test_figure_amplitude_histograms_transforms(puncta, transform):
    lipid, _ = puncta
    fig = figure_amplitude_histograms({"A_ch1": lipid}, transform=transform)
    _assert_bare_figure(fig)


def test_figure_diameter_distribution(puncta):
    lipid, _ = puncta
    fig = figure_diameter_distribution(lipid, conversion_factor=2.5)
    _assert_bare_figure(fig)
    assert "diameter" in fig.axes[0].get_xlabel().lower()


def test_figure_protein_amplitude(puncta):
    _, protein = puncta
    _assert_bare_figure(figure_protein_amplitude(protein))


def test_figure_dls_overlay(calibration):
    fig = figure_dls_overlay(calibration)
    _assert_bare_figure(fig)
    assert len(fig.axes) == 3          # two panels plus the nm twin axis
    assert fig.axes[0].get_xscale() == "log"


def test_figure_bootstrap(calibration):
    boot = run_bootstrap(calibration, n_boot=8, k_per=300)
    _assert_bare_figure(figure_bootstrap(boot))


def test_figure_dls_distribution(sample_dls):
    from core.calibration import load_dls_section
    d, w = load_dls_section(sample_dls, "number")
    _assert_bare_figure(figure_dls_distribution(d, w))
    _assert_bare_figure(figure_dls_distribution(d, w, log=True, zoom_pct=100))


def test_figure_dls_comparison(sample_dls, puncta):
    from core.calibration import load_dls_section
    d, w = load_dls_section(sample_dls, "number")
    lipid, protein = puncta
    fig = figure_dls_comparison(
        d, w,
        channels=[("A_ch1", "Lipid"), ("A_ch2", "EGFP")],
        fluor_sqrt={"A_ch1": np.sqrt(lipid), "A_ch2": np.sqrt(protein)},
    )
    _assert_bare_figure(fig)
    assert len(fig.axes) == 3          # DLS panel + two channel panels


def test_figure_threshold_preview(puncta):
    lipid, _ = puncta
    fig = figure_threshold_preview(lipid, threshold=20.0, kept=300, total=800)
    _assert_bare_figure(fig)
    assert "300" in fig.axes[0].get_title()
    assert "800" in fig.axes[0].get_title()


# ── Saving ─────────────────────────────────────────────────────────────

def test_save_figure_creates_parents(curve, tmp_path):
    out = tmp_path / "nested" / "deeper" / "fig.png"
    returned = save_figure(figure_curvature(curve, "t", 2.0), out, dpi=72)
    assert returned == out
    assert out.is_file()
    assert out.stat().st_size > 0


def test_builders_do_not_write(curve, tmp_path, monkeypatch):
    """Building a figure must never touch the filesystem."""
    before = set(tmp_path.rglob("*"))
    monkeypatch.chdir(tmp_path)
    figure_curvature(curve, "t", 2.0)
    figure_amplitude_histograms({"A_ch1": curve.diameter})
    assert set(tmp_path.rglob("*")) == before
