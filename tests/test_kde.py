"""Tests for pynebulosa._kde — core KDE algorithms."""

import numpy as np
import pytest
from scipy.stats import gaussian_kde

from pynebulosa._kde import _bandwidth, _get_dens, calculate_density, wkde2d


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

    def test_zero_iqr_falls_back_to_std(self):
        """When all values are identical, IQR=0 and std is used."""
        x = np.ones(50) + np.array([0.0] * 49 + [0.001])
        h = _bandwidth(x)
        assert h > 0

    def test_constant_array(self):
        """Constant array has IQR=0 and std~=0, bandwidth should still work."""
        x = np.full(100, 5.0)
        h = _bandwidth(x)
        assert np.isfinite(h)


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

    def test_infinite_values_x_raises(self):
        with pytest.raises(ValueError, match="infinite"):
            wkde2d(
                np.array([1.0, np.inf]),
                np.array([1.0, 2.0]),
                np.array([1.0, 1.0]),
            )

    def test_infinite_values_y_raises(self):
        with pytest.raises(ValueError, match="infinite"):
            wkde2d(
                np.array([1.0, 2.0]),
                np.array([1.0, np.nan]),
                np.array([1.0, 1.0]),
            )

    def test_non_finite_lims_raises(self):
        """Covers line 85: lims validation."""
        with pytest.raises(ValueError, match="finite"):
            wkde2d(
                np.array([1.0, 2.0]),
                np.array([1.0, 2.0]),
                np.array([1.0, 1.0]),
                lims=(0.0, np.inf, 0.0, 1.0),
            )

    def test_zero_weights_produce_zero_density(self):
        rng = np.random.default_rng(4)
        n_pts = 50
        result = wkde2d(
            rng.normal(size=n_pts),
            rng.normal(size=n_pts),
            np.zeros(n_pts),
        )
        assert np.allclose(result["z"], 0)

    def test_custom_bandwidth(self):
        rng = np.random.default_rng(7)
        n_pts = 50
        x, y, w = rng.normal(size=n_pts), rng.normal(size=n_pts), rng.uniform(size=n_pts)
        result = wkde2d(x, y, w, h=(0.5, 0.5))
        assert result["z"].shape == (100, 100)
        assert np.all(result["z"] >= 0)

    def test_custom_adjust(self):
        rng = np.random.default_rng(8)
        n_pts = 50
        x, y, w = rng.normal(size=n_pts), rng.normal(size=n_pts), rng.uniform(size=n_pts)
        narrow = wkde2d(x, y, w, adjust=0.5)
        wide = wkde2d(x, y, w, adjust=2.0)
        # Wider bandwidth should produce a smoother (lower max) density
        assert narrow["z"].max() > wide["z"].max()

    def test_custom_lims(self):
        rng = np.random.default_rng(9)
        n_pts = 50
        x, y, w = rng.normal(size=n_pts), rng.normal(size=n_pts), rng.uniform(size=n_pts)
        result = wkde2d(x, y, w, lims=(-5.0, 5.0, -5.0, 5.0))
        assert result["x"][0] == pytest.approx(-5.0)
        assert result["x"][-1] == pytest.approx(5.0)


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

    def test_coords_outside_grid_clipped(self):
        """Coordinates beyond grid range should be clipped, not error."""
        grid_x = np.linspace(0, 1, 10)
        grid_y = np.linspace(0, 1, 10)
        z = np.ones((10, 10))
        # Coords well outside grid
        coords = np.array([[-10.0, -10.0], [10.0, 10.0]])
        dens = _get_dens(coords, grid_x, grid_y, z)
        assert dens.shape == (2,)
        assert np.all(np.isfinite(dens))


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
        result = calculate_density(w, coords, method="wkde", map_to_cells=False)
        assert isinstance(result, dict)
        assert "x" in result and "y" in result and "z" in result

    def test_ks_raw_returns_kde_object(self, adata_simple):
        """Covers line 220: ks method with map_to_cells=False."""
        coords = adata_simple.obsm["X_umap"]
        w = np.asarray(adata_simple.X.toarray()[:, 0]).flatten()
        result = calculate_density(w, coords, method="ks", map_to_cells=False)
        assert isinstance(result, gaussian_kde)

    def test_zero_weights_returns_zeros(self, adata_simple):
        coords = adata_simple.obsm["X_umap"]
        w = np.zeros(adata_simple.n_obs)
        with pytest.warns(UserWarning, match="All weights are zero"):
            result = calculate_density(w, coords, method="wkde")
        assert np.all(result == 0)

    def test_zero_weights_raw_returns_empty(self, adata_simple):
        """Covers line 197: zero weights with map_to_cells=False."""
        coords = adata_simple.obsm["X_umap"]
        w = np.zeros(adata_simple.n_obs)
        with pytest.warns(UserWarning, match="All weights are zero"):
            result = calculate_density(w, coords, method="wkde", map_to_cells=False)
        assert isinstance(result, dict)
        assert len(result["x"]) == 0

    def test_invalid_method_raises(self, adata_simple):
        coords = adata_simple.obsm["X_umap"]
        w = np.ones(adata_simple.n_obs)
        with pytest.raises(ValueError, match="Unknown method"):
            calculate_density(w, coords, method="invalid")

    def test_density_higher_where_expressed(self, adata_simple):
        coords = adata_simple.obsm["X_umap"]
        w = np.asarray(adata_simple.X.toarray()[:, 0]).flatten()
        dens = calculate_density(w, coords, method="wkde")
        assert dens[:50].mean() > dens[100:150].mean()

    def test_negative_weights_clipped(self, adata_simple):
        """Negative expression values (z-scored data) should be clipped to 0."""
        coords = adata_simple.obsm["X_umap"]
        w = np.array([-1.0, -2.0] * 100, dtype=float)
        with pytest.warns(UserWarning, match="All weights are zero"):
            result = calculate_density(w, coords, method="wkde")
        assert np.all(result == 0)

    def test_mixed_negative_positive_weights(self, adata_simple):
        """Negative weights should be clipped; positive ones kept."""
        coords = adata_simple.obsm["X_umap"]
        w = np.full(adata_simple.n_obs, -1.0)
        w[:50] = 5.0  # Only cluster 1 has positive weights
        dens = calculate_density(w, coords, method="wkde")
        assert dens[:50].mean() > dens[100:150].mean()

    def test_accepts_list_input(self):
        """Input arrays can be plain Python lists."""
        coords = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]]
        w = [1.0, 2.0, 3.0]
        result = calculate_density(w, coords)
        assert len(result) == 3
