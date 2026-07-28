"""Tests for core.curvature -- amplitude to diameter to protein density."""

import numpy as np
import pytest

from core.curvature import (
    amplitude_to_diameter,
    bin_by_diameter,
    compute_protein_density,
    curvature_curve,
    normalize,
)


def test_conversion_factor_mode():
    lipid = np.array([4.0, 9.0, 16.0])
    d, scale = amplitude_to_diameter(lipid, conversion_factor=10.0)
    assert scale == 10.0
    np.testing.assert_allclose(d, [20.0, 30.0, 40.0])


def test_dls_mean_diameter_mode():
    lipid = np.array([4.0, 9.0, 16.0])  # sqrt -> 2, 3, 4; mean 3
    d, scale = amplitude_to_diameter(lipid, dls_mean_diameter_nm=90.0)
    assert scale == pytest.approx(30.0)
    np.testing.assert_allclose(d, [60.0, 90.0, 120.0])
    assert np.mean(d) == pytest.approx(90.0)


def test_negative_amplitudes_clipped():
    d, _ = amplitude_to_diameter(np.array([-4.0, 9.0]), conversion_factor=1.0)
    np.testing.assert_allclose(d, [0.0, 3.0])


def test_requires_exactly_one_mode():
    with pytest.raises(ValueError):
        amplitude_to_diameter(np.array([1.0]))


def test_conversion_factor_wins_over_mean():
    """Matches the CLI, which rejects both but prefers the factor internally."""
    d, scale = amplitude_to_diameter(
        np.array([4.0]), dls_mean_diameter_nm=100.0, conversion_factor=3.0
    )
    assert scale == 3.0
    np.testing.assert_allclose(d, [6.0])


def test_protein_density():
    diameter = np.array([10.0, 20.0])
    protein = np.array([np.pi * 100.0, np.pi * 400.0])
    density, valid = compute_protein_density(protein, diameter)
    assert valid.all()
    np.testing.assert_allclose(density, [1.0, 1.0], rtol=1e-9)


def test_protein_density_zero_diameter_invalid():
    density, valid = compute_protein_density(np.array([1.0]), np.array([0.0]))
    assert not valid[0]


def test_bin_by_diameter():
    x = np.array([0.0, 0.4, 1.2, 1.6, 5.0])
    y = np.array([1.0, 3.0, 10.0, 20.0, 7.0])
    centres, means = bin_by_diameter(x, y, bin_width=1.0)

    assert centres[0] == pytest.approx(0.5)
    assert means[0] == pytest.approx(2.0)     # (1 + 3) / 2
    assert means[1] == pytest.approx(15.0)    # (10 + 20) / 2
    assert len(centres) == len(means)
    # Empty bins are dropped entirely.
    assert len(centres) == 3


def test_bin_by_diameter_empty():
    centres, means = bin_by_diameter(np.array([]), np.array([]), 1.0)
    assert centres.size == 0 and means.size == 0


def test_normalize_modes():
    centres = np.array([1.0, 2.0, 3.0])
    means = np.array([4.0, 2.0, 8.0])

    np.testing.assert_allclose(normalize(centres, means, "none"), means)
    np.testing.assert_allclose(normalize(centres, means, "rightmost"), means / 8.0)
    np.testing.assert_allclose(normalize(centres, means, "leftmost"), means / 4.0)
    np.testing.assert_allclose(normalize(centres, means, "minimum"), means / 2.0)

    with pytest.raises(ValueError):
        normalize(centres, means, "sideways")


def test_normalize_rejects_nonpositive_reference():
    with pytest.raises(ValueError):
        normalize(np.array([1.0, 2.0]), np.array([1.0, 0.0]), "rightmost")


def test_curvature_curve_with_cutoff():
    lipid = np.array([1.0, 4.0, 9.0, 100.0])      # sqrt -> 1, 2, 3, 10
    protein = np.array([1.0, 1.0, 1.0, 1.0])

    curve = curvature_curve(lipid, protein, bin_width=1.0, conversion_factor=10.0)
    assert curve.n_points == 4
    assert curve.n_points_before_cutoff == 4

    cut = curvature_curve(lipid, protein, bin_width=1.0,
                          conversion_factor=10.0, diameter_cutoff=50.0)
    assert cut.n_points_before_cutoff == 4
    assert cut.n_points == 3                      # the 100 nm punctum is dropped
    assert cut.diameter.max() <= 50.0
    assert cut.scale == 10.0
