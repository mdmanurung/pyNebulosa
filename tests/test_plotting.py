"""Tests for nebulosa._plotting — plot_density and helpers."""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

from nebulosa import plot_density


# Use non-interactive backend for tests
matplotlib.use("Agg")


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


class TestPlotDensityMultiFeature:
    def test_returns_figure(self, adata_simple):
        result = plot_density(
            adata_simple, ["marker1", "marker2"], show=False
        )
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_joint_adds_panel(self, adata_simple):
        fig = plot_density(
            adata_simple, ["marker1", "marker2"], joint=True, show=False
        )
        # Should have 3 axes: marker1, marker2, joint
        visible_axes = [ax for ax in fig.axes if ax.get_visible()]
        # Each scatter has an associated colorbar axes, so 6 total (3 scatter + 3 cbar)
        assert len(visible_axes) >= 3
        plt.close("all")


class TestPlotDensityParameters:
    def test_custom_reduction(self, adata_simple):
        result = plot_density(
            adata_simple, "marker1", reduction="X_pca",
            dims=(0, 1), show=False,
        )
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_ks_method(self, adata_simple):
        result = plot_density(
            adata_simple, "marker1", method="ks", show=False
        )
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_custom_cmap(self, adata_simple):
        result = plot_density(
            adata_simple, "marker1", cmap="magma", show=False
        )
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_metadata_feature(self, adata_simple):
        """Plotting a numeric obs column should work."""
        # Add a numeric metadata column
        adata_simple.obs["n_counts"] = np.random.default_rng(0).uniform(
            size=adata_simple.n_obs
        )
        result = plot_density(adata_simple, "n_counts", show=False)
        assert isinstance(result, plt.Figure)
        plt.close("all")

    def test_adjust_bandwidth(self, adata_simple):
        result = plot_density(
            adata_simple, "marker1", adjust=2.0, show=False
        )
        assert isinstance(result, plt.Figure)
        plt.close("all")


class TestPlotDensityErrors:
    def test_missing_feature_raises(self, adata_simple):
        with pytest.raises(ValueError, match="not found"):
            plot_density(adata_simple, "nonexistent_gene", show=False)

    def test_invalid_reduction_raises(self, adata_simple):
        with pytest.raises(KeyError, match="not found"):
            plot_density(
                adata_simple, "marker1", reduction="X_invalid", show=False
            )

    def test_invalid_dims_raises(self, adata_simple):
        with pytest.raises(ValueError):
            plot_density(
                adata_simple, "marker1", dims=(0, 1, 2), show=False
            )

    def test_out_of_range_dims_raises(self, adata_simple):
        with pytest.raises(ValueError, match="not present"):
            plot_density(
                adata_simple, "marker1", dims=(0, 99), show=False
            )

    def test_empty_features_raises(self, adata_simple):
        with pytest.raises(ValueError, match="At least one"):
            plot_density(adata_simple, [], show=False)

    def test_invalid_method_raises(self, adata_simple):
        with pytest.raises(ValueError, match="Unknown method"):
            plot_density(
                adata_simple, "marker1", method="invalid", show=False
            )

    def test_invalid_layer_raises(self, adata_simple):
        with pytest.raises(ValueError, match="Layer"):
            plot_density(
                adata_simple, "marker1", layer="nonexistent", show=False
            )


class TestPlotDensitySave:
    def test_save_to_file(self, adata_simple, tmp_path):
        outpath = tmp_path / "test_density.png"
        plot_density(adata_simple, "marker1", save=str(outpath), show=False)
        assert outpath.exists()
        plt.close("all")
