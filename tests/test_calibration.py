"""Tests for core.calibration against the bundled Zetasizer export."""

import numpy as np
import pytest

from core.calibration import (
    calibrate,
    dls_bin_edges,
    load_dls_number,
    load_dls_section,
    load_sqrt_A,
    overlay_cost,
    ratio_of_means_conversion,
    run_bootstrap,
)

K_TRUE = 2.5  # nm per sqrt(A) used to synthesise the fluorescence sample


@pytest.fixture(scope="module")
def dls(sample_dls):
    return load_dls_number(sample_dls)


@pytest.fixture(scope="module")
def synthetic_sqrt_A(dls):
    """
    Draw diameters from the DLS number distribution and divide by K_TRUE.

    A correct calibration must recover K_TRUE from this.
    """
    d, w = dls
    rng = np.random.default_rng(12345)
    sampled = rng.choice(d, size=40000, p=w / w.sum())
    return sampled / K_TRUE


def test_load_dls_number_shape_and_mean(dls):
    d, w = dls
    assert len(d) == len(w) == 69
    assert d[0] == pytest.approx(0.4)
    assert np.all(np.diff(d) > 0), "diameter bins must be ascending"
    weighted_mean = float(np.sum(d * w) / np.sum(w))
    assert weighted_mean == pytest.approx(80.1166, abs=1e-3)


def test_load_dls_section_drops_empty_bins(sample_dls, dls):
    d_all, _w_all = dls
    d_plot, w_plot = load_dls_section(sample_dls, "number")
    assert len(d_plot) < len(d_all), "plotting loader should drop zero-weight bins"
    assert np.all(w_plot > 0)


def test_load_dls_section_missing_raises(sample_dls):
    with pytest.raises(ValueError, match="Could not find"):
        load_dls_section(sample_dls, "nonsense")


def test_dls_bin_edges_bracket_centres(dls):
    d, _w = dls
    edges = dls_bin_edges(d)
    assert len(edges) == len(d) + 1
    assert np.all(np.diff(edges) > 0)
    assert np.all(edges[:-1] < d)
    assert np.all(d < edges[1:])


def test_ratio_of_means(dls, synthetic_sqrt_A):
    d, w = dls
    factor, dls_mean, fluor_mean = ratio_of_means_conversion(d, w, synthetic_sqrt_A)
    assert dls_mean == pytest.approx(80.1166, abs=1e-3)
    assert factor == pytest.approx(dls_mean / fluor_mean)
    assert factor == pytest.approx(K_TRUE, rel=0.05)


def test_calibrate_recovers_known_factor(dls, synthetic_sqrt_A):
    d, w = dls
    cal = calibrate(d, w, synthetic_sqrt_A)

    assert cal.conversion_factor == pytest.approx(K_TRUE, rel=0.05)
    assert cal.implied_mean_diameter == pytest.approx(
        float(np.mean(synthetic_sqrt_A)) * cal.conversion_factor
    )
    assert cal.n_puncta == len(synthetic_sqrt_A)
    assert cal.chi2 >= 0

    # Shapes needed to rebuild the overlay figure without re-reading the file.
    assert cal.bin_edges.shape == (len(d) + 1,)
    assert cal.bin_widths.shape == (len(d),)
    assert cal.dls_density_norm.shape == (len(d),)
    assert cal.fluor_density_norm.shape == (len(d),)
    assert cal.dls_density_norm.max() == pytest.approx(1.0)
    assert cal.fluor_density_norm.max() == pytest.approx(1.0)


def test_methods_agree_on_clean_data(dls, synthetic_sqrt_A):
    d, w = dls
    cal = calibrate(d, w, synthetic_sqrt_A)
    assert cal.percent_difference < 5.0
    assert not cal.methods_disagree(tolerance_pct=5.0)
    assert cal.methods_disagree(tolerance_pct=0.0)


def test_overlay_cost_minimised_at_truth(dls, synthetic_sqrt_A):
    d, w = dls
    cal = calibrate(d, w, synthetic_sqrt_A)
    best = overlay_cost(cal.conversion_factor, synthetic_sqrt_A, cal.bin_edges,
                        cal.dls_density_norm, cal.bin_widths)
    for offset in (0.5, 1.5):
        worse = overlay_cost(cal.conversion_factor * offset, synthetic_sqrt_A,
                             cal.bin_edges, cal.dls_density_norm, cal.bin_widths)
        assert worse > best


def test_overlay_cost_infinite_when_empty(dls):
    d, w = dls
    edges = dls_bin_edges(d)
    widths = np.diff(edges)
    density = (w / widths) / (w / widths).max()
    # A tiny k pushes every punctum below the first bin edge.
    assert overlay_cost(1e-12, np.array([1.0, 2.0]), edges, density, widths) == np.inf


def test_bootstrap(dls):
    d, w = dls
    rng = np.random.default_rng(7)
    sqrt_A = rng.choice(d, size=2000, p=w / w.sum()) / K_TRUE
    cal = calibrate(d, w, sqrt_A)

    seen = []
    boot = run_bootstrap(cal, n_boot=12, k_per=500,
                         progress=lambda i, n: seen.append(i))

    assert len(boot.overlay_factors) == 12
    assert len(boot.rom_factors) == 12
    assert seen == list(range(1, 13))
    assert boot.overlay_cv > 0
    assert boot.rom_cv > 0
    assert isinstance(boot.overlay_is_tighter, bool)
    assert cal.bootstrap is boot
    assert np.mean(boot.rom_factors) == pytest.approx(K_TRUE, rel=0.10)


def test_load_sqrt_A_is_positive_only(tmp_path):
    p = tmp_path / "filtered.txt"
    p.write_text(
        "# comment\nsource_image\tA_ch1\tA_ch2\n"
        "c|r0\t4\t1\nc|r1\t0\t1\nc|r2\t9\t1\nc|r3\t-3\t1\n",
        encoding="utf-8",
    )
    vals = load_sqrt_A(p, "A_ch1")
    np.testing.assert_allclose(sorted(vals), [2.0, 3.0])


def test_load_sqrt_A_bad_column(tmp_path):
    p = tmp_path / "filtered.txt"
    p.write_text("source_image\tA_ch1\nc|r0\t4\n", encoding="utf-8")
    with pytest.raises(ValueError, match="not found"):
        load_sqrt_A(p, "A_ch7")
