"""Tests for nebulosa._kde — core KDE algorithms."""

import numpy as np
import pytest

from nebulosa._kde import _bandwidth, _get_dens, calculate_density, wkde2d


class TestBandwidth:
    def test_returns_positive_float(self):
        rng = np.random.default_rng(0)
        x = rng.normal(size=100)
        h = _bandwidth(x)
        assert isinstance(h, float)
        assert h > 0

    def test_larger_for_wider_data(self):
        rng = np.random.default_rng(0)
        narrow = rng.normal(scale=1, size=200)
        wide = rng.normal(scale=10, size=200)
        assert _bandwidth(wide) > _bandwidth(narrow)


class TestWkde2d:
    def test_output_shape(self):
        rng = np.random.default_rng(1)
        n_pts = 50
        x = rng.normal(size=n_pts)
        y = rng.normal(size=n_pts)
        w = rng.uniform(size=n_pts)
        result = wkde2d(x, y, w, n=80)
        assert result["x"].shape == (80,)
        assert result["y"].shape == (80,)
        assert result["z"].shape == (80, 80)

    def test_default_grid_size(self):
        rng = np.random.default_rng(2)
        n_pts = 30
        result = wkde2d(
            rng.normal(size=n_pts),
            rng.normal(size=n_pts),
            rng.uniform(size=n_pts),
        )
        assert result["z"].shape == (100, 100)

    def test_density_nonnegative(self):
        rng = np.random.default_rng(3)
        n_pts = 100
        result = wkde2d(
            rng.normal(size=n_pts),
            rng.normal(size=n_pts),
            rng.uniform(size=n_pts),
        )
        assert np.all(result["z"] >= 0)

    def test_mismatched_lengths_raises(self):
        with pytest.raises(ValueError, match="same length"):
            wkde2d(np.array([1, 2]), np.array([1, 2, 3]), np.array([1, 2]))

    def test_infinite_values_raises(self):
        with pytest.raises(ValueError, match="infinite"):
            wkde2d(
                np.array([1.0, np.inf]),
                np.array([1.0, 2.0]),
                np.array([1.0, 1.0]),
            )

    def test_zero_weights_produce_zero_density(self):
        rng = np.random.default_rng(4)
        n_pts = 50
        result = wkde2d(
            rng.normal(size=n_pts),
            rng.normal(size=n_pts),
            np.zeros(n_pts),
        )
        # With all-zero weights, density should be zero everywhere
        # (0 * norm.pdf = 0 for all grid points)
        assert np.allclose(result["z"], 0)


class TestGetDens:
    def test_returns_correct_length(self):
        rng = np.random.default_rng(5)
        n_pts = 60
        coords = rng.normal(size=(n_pts, 2))
        result = wkde2d(coords[:, 0], coords[:, 1], rng.uniform(size=n_pts))
        dens = _get_dens(coords, result["x"], result["y"], result["z"])
        assert dens.shape == (n_pts,)

    def test_values_nonnegative(self):
        rng = np.random.default_rng(6)
        n_pts = 40
        coords = rng.normal(size=(n_pts, 2))
        result = wkde2d(coords[:, 0], coords[:, 1], rng.uniform(size=n_pts))
        dens = _get_dens(coords, result["x"], result["y"], result["z"])
        assert np.all(dens >= 0)


class TestCalculateDensity:
    def test_wkde_returns_vector(self, adata_simple):
        coords = adata_simple.obsm["X_umap"]
        w = np.asarray(adata_simple.X.toarray()[:, 0]).flatten()
        result = calculate_density(w, coords, method="wkde")
        assert result.shape == (adata_simple.n_obs,)

    def test_ks_returns_vector(self, adata_simple):
        coords = adata_simple.obsm["X_umap"]
        w = np.asarray(adata_simple.X.toarray()[:, 0]).flatten()
        result = calculate_density(w, coords, method="ks")
        assert result.shape == (adata_simple.n_obs,)

    def test_wkde_raw_returns_dict(self, adata_simple):
        coords = adata_simple.obsm["X_umap"]
        w = np.asarray(adata_simple.X.toarray()[:, 0]).flatten()
        result = calculate_density(w, coords, method="wkde", map=False)
        assert isinstance(result, dict)
        assert "x" in result and "y" in result and "z" in result

    def test_zero_weights_returns_zeros(self, adata_simple):
        coords = adata_simple.obsm["X_umap"]
        w = np.zeros(adata_simple.n_obs)
        with pytest.warns(UserWarning, match="All weights are zero"):
            result = calculate_density(w, coords, method="wkde")
        assert np.all(result == 0)

    def test_invalid_method_raises(self, adata_simple):
        coords = adata_simple.obsm["X_umap"]
        w = np.ones(adata_simple.n_obs)
        with pytest.raises(ValueError, match="Unknown method"):
            calculate_density(w, coords, method="invalid")

    def test_density_higher_where_expressed(self, adata_simple):
        """Cells with high marker1 expression (cluster 1, top-right)
        should have higher density than non-expressing cells."""
        coords = adata_simple.obsm["X_umap"]
        w = np.asarray(adata_simple.X.toarray()[:, 0]).flatten()
        dens = calculate_density(w, coords, method="wkde")
        # Mean density for cluster 1 (cells 0-49, high marker1)
        # should exceed mean density for cluster 3 (cells 100-149, low marker1)
        assert dens[:50].mean() > dens[100:150].mean()
