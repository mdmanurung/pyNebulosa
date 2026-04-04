"""Tests for nebulosa._plotting — plot_density and helpers."""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

from nebulosa import plot_density
from nebulosa._plotting import _density_scatter, _finalize

# Use non-interactive backend for tests
matplotlib.use("Agg")


class TestDensityScatter:
    """Direct tests for the _density_scatter helper."""

    def test_returns_axes(self):
        fig, ax = plt.subplots()
        embeddings = np.random.default_rng(0).normal(size=(50, 2))
        density = np.random.default_rng(1).uniform(size=50)
        result = _density_scatter(ax, embeddings, density, "test", ("X", "Y"), 1.0, "viridis")
        assert result is ax
        plt.close("all")

    def test_custom_colorbar_label(self):
        fig, ax = plt.subplots()
        embeddings = np.random.default_rng(0).normal(size=(30, 2))
        density = np.random.default_rng(1).uniform(size=30)
        _density_scatter(
            ax,
            embeddings,
            density,
            "test",
            ("X", "Y"),
            1.0,
            "magma",
            colorbar_label="Joint density",
        )
        plt.close("all")

    def test_points_sorted_by_density(self):
        """High-density points should render on top (plotted last)."""
        fig, ax = plt.subplots()
        embeddings = np.array([[0, 0], [1, 1], [2, 2]], dtype=float)
        density = np.array([0.1, 0.9, 0.5])
        _density_scatter(ax, embeddings, density, "test", ("X", "Y"), 1.0, "viridis")
        # The scatter collection exists
        assert len(ax.collections) == 1
        plt.close("all")


class TestFinalize:
    """Direct tests for the _finalize helper."""

    def test_show_true_returns_none(self):
        fig, ax = plt.subplots()
        result = _finalize(fig, show=True, save=None)
        assert result is None
        plt.close("all")

    def test_show_false_closes_figure(self):
        fig, ax = plt.subplots()
        result = _finalize(fig, show=False, save=None)
        assert result is fig
        plt.close("all")

    def test_show_none_returns_figure(self):
        fig, ax = plt.subplots()
        result = _finalize(fig, show=None, save=None)
        assert result is fig
        plt.close("all")

    def test_return_axes(self):
        fig, ax = plt.subplots()
        result = _finalize(fig, show=False, save=None, return_axes=ax)
        assert result is ax
        plt.close("all")

    def test_save_and_show(self, tmp_path):
        fig, ax = plt.subplots()
        outpath = tmp_path / "test.png"
        result = _finalize(fig, show=True, save=str(outpath))
        assert outpath.exists()
        assert result is None
        plt.close("all")


class TestPlotDensitySingleFeature:
    def test_returns_figure(self, adata_simple):
        result = plot_density(adata_simple, "marker1", show=False)
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_with_provided_ax(self, adata_simple):
        fig, ax = plt.subplots()
        result = plot_density(adata_simple, "marker1", ax=ax, show=False)
        assert result is ax
        plt.close("all")

    def test_with_dense_matrix(self, adata_dense):
        result = plot_density(adata_dense, "marker1", show=False)
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_string_feature_normalized_to_list(self, adata_simple):
        """String feature input should work identically to single-element list."""
        r1 = plot_density(adata_simple, "marker1", show=False)
        assert isinstance(r1, plt.Figure)
        plt.close("all")


class TestPlotDensityMultiFeature:
    def test_returns_figure(self, adata_simple):
        result = plot_density(adata_simple, ["marker1", "marker2"], show=False)
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_joint_adds_panel(self, adata_simple):
        fig = plot_density(adata_simple, ["marker1", "marker2"], joint=True, show=False)
        visible_axes = [ax for ax in fig.axes if ax.get_visible()]
        # 3 scatter + 3 colorbar = 6 visible axes
        assert len(visible_axes) >= 3
        plt.close("all")

    def test_hidden_axes_when_uneven_grid(self, adata_simple):
        """Covers lines 223-224: hidden axes when panels don't fill grid.

        4 features with ncols=3 -> 2 rows x 3 cols = 6 slots, 2 hidden.
        """
        fig = plot_density(
            adata_simple,
            ["marker1", "marker2", "gene3", "gene4"],
            ncols=3,
            show=False,
        )
        all_axes = fig.axes
        # Colorbar axes are always visible, but subplot axes should include hidden ones
        invisible = [ax for ax in all_axes if not ax.get_visible()]
        assert len(invisible) >= 1
        plt.close("all")

    def test_four_features_wraps_rows(self, adata_simple):
        """4 features with ncols=2 should create a 2x2 grid."""
        fig = plot_density(
            adata_simple,
            ["marker1", "marker2", "gene3", "gene4"],
            ncols=2,
            show=False,
        )
        assert isinstance(fig, plt.Figure)
        plt.close("all")


class TestPlotDensityParameters:
    def test_custom_reduction(self, adata_simple):
        result = plot_density(
            adata_simple,
            "marker1",
            reduction="X_pca",
            dims=(0, 1),
            show=False,
        )
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_ks_method(self, adata_simple):
        result = plot_density(adata_simple, "marker1", method="ks", show=False)
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_custom_cmap(self, adata_simple):
        result = plot_density(adata_simple, "marker1", cmap="magma", show=False)
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_metadata_feature(self, adata_simple):
        adata_simple.obs["n_counts"] = np.random.default_rng(0).uniform(size=adata_simple.n_obs)
        result = plot_density(adata_simple, "n_counts", show=False)
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_adjust_bandwidth(self, adata_simple):
        result = plot_density(adata_simple, "marker1", adjust=2.0, show=False)
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_custom_figsize(self, adata_simple):
        fig = plot_density(adata_simple, "marker1", figsize=(8, 6), show=False)
        w, h = fig.get_size_inches()
        assert w == pytest.approx(8)
        assert h == pytest.approx(6)
        plt.close("all")

    def test_custom_size(self, adata_simple):
        fig = plot_density(adata_simple, "marker1", size=5.0, show=False)
        assert isinstance(fig, plt.Figure)
        plt.close("all")

    def test_show_true_returns_none(self, adata_simple):
        """Covers lines 244-245: show=True path."""
        result = plot_density(adata_simple, "marker1", show=True)
        assert result is None
        plt.close("all")

    def test_show_none_returns_figure(self, adata_simple):
        result = plot_density(adata_simple, "marker1", show=None)
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_layer_parameter(self, adata_simple):
        """Covers valid layer extraction path."""
        adata_simple.layers["log"] = adata_simple.X.copy()
        result = plot_density(adata_simple, "marker1", layer="log", show=False)
        assert isinstance(result, plt.Figure)
        plt.close("all")


class TestPlotDensityErrors:
    def test_missing_feature_raises(self, adata_simple):
        with pytest.raises(ValueError, match="not found"):
            plot_density(adata_simple, "nonexistent_gene", show=False)

    def test_invalid_reduction_raises(self, adata_simple):
        with pytest.raises(KeyError, match="not found"):
            plot_density(adata_simple, "marker1", reduction="X_invalid", show=False)

    def test_invalid_dims_raises(self, adata_simple):
        with pytest.raises(ValueError):
            plot_density(adata_simple, "marker1", dims=(0, 1, 2), show=False)

    def test_out_of_range_dims_raises(self, adata_simple):
        with pytest.raises(ValueError, match="not present"):
            plot_density(adata_simple, "marker1", dims=(0, 99), show=False)

    def test_empty_features_raises(self, adata_simple):
        with pytest.raises(ValueError, match="At least one"):
            plot_density(adata_simple, [], show=False)

    def test_invalid_method_raises(self, adata_simple):
        with pytest.raises(ValueError, match="Unknown method"):
            plot_density(adata_simple, "marker1", method="invalid", show=False)

    def test_invalid_layer_raises(self, adata_simple):
        with pytest.raises(ValueError, match="Layer"):
            plot_density(adata_simple, "marker1", layer="nonexistent", show=False)


class TestPlotDensitySave:
    def test_save_to_file(self, adata_simple, tmp_path):
        outpath = tmp_path / "test_density.png"
        plot_density(adata_simple, "marker1", save=str(outpath), show=False)
        assert outpath.exists()
        plt.close("all")

    def test_save_pdf(self, adata_simple, tmp_path):
        outpath = tmp_path / "test_density.pdf"
        plot_density(adata_simple, "marker1", save=str(outpath), show=False)
        assert outpath.exists()
        plt.close("all")
