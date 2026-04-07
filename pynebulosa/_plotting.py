"""Plotting functions for gene-weighted density visualization.

Port of R/plotting.R (``plot_density_``) and the orchestration
logic from R/helpers.R (``.plot_final_density``).
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from pynebulosa._kde import calculate_density
from pynebulosa._utils import _get_embeddings, _get_feature_data

if TYPE_CHECKING:
    import matplotlib.axes
    import matplotlib.figure
    from anndata import AnnData


def _density_scatter(
    ax: matplotlib.axes.Axes,
    embeddings: np.ndarray,
    density: np.ndarray,
    title: str,
    dim_names: tuple[str, str],
    size: float,
    cmap: str,
    colorbar_label: str = "Density",
) -> matplotlib.axes.Axes:
    """Create a single density scatter panel.

    Port of R's ``plot_density_`` function.

    Parameters
    ----------
    ax
        Matplotlib axes to draw on.
    embeddings
        Cell coordinates, shape ``(N, 2)``.
    density
        Density values per cell, length N.
    title
        Plot title (typically the feature name).
    dim_names
        Axis labels for x and y.
    size
        Point size.
    cmap
        Matplotlib colormap name (e.g. ``"viridis"``).
    colorbar_label
        Label for the colorbar.

    Returns
    -------
    The axes with the scatter plot drawn.
    """
    # Sort by density so high-density points render on top
    order = np.argsort(density)
    x = embeddings[order, 0]
    y = embeddings[order, 1]
    c = density[order]

    sc = ax.scatter(
        x,
        y,
        c=c,
        cmap=cmap,
        s=size,
        edgecolors="none",
        rasterized=True,
    )
    ax.set_xlabel(dim_names[0])
    ax.set_ylabel(dim_names[1])
    ax.set_title(title, fontsize=14)

    # Style to match R version: clean background
    ax.set_facecolor("white")
    ax.tick_params(colors="black")
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
        spine.set_color("black")

    plt.colorbar(sc, ax=ax, label=colorbar_label, shrink=0.8)

    return ax


def plot_density(
    adata: AnnData,
    features: str | list[str],
    *,
    joint: bool = False,
    reduction: str | None = None,
    layer: str | None = None,
    dims: tuple[int, int] = (0, 1),
    method: str = "wkde",
    adjust: float = 1.0,
    size: float = 1.0,
    cmap: str = "viridis",
    figsize: tuple[float, float] | None = None,
    ncols: int = 3,
    show: bool | None = None,
    save: str | None = None,
    ax: matplotlib.axes.Axes | None = None,
) -> matplotlib.figure.Figure | matplotlib.axes.Axes | None:
    """Plot gene-weighted 2D kernel density.

    Main user-facing function. Port of R's ``plot_density()`` generic.

    Creates scatter plots of cell embeddings colored by gene-weighted
    kernel density estimates. For multiple features, creates a grid of
    subplots. When ``joint=True``, appends a panel showing the product
    of individual densities.

    Parameters
    ----------
    adata
        AnnData object with computed dimensionality reduction.
    features
        Gene name(s) or obs column name(s) to visualize.
    joint
        If ``True`` and multiple features are given, append a joint
        density plot (product of individual densities).
    reduction
        Key in ``adata.obsm`` (e.g. ``"X_umap"``). Auto-detected
        if ``None``.
    layer
        Layer in ``adata.layers`` for expression data. Uses
        ``adata.X`` if ``None``.
    dims
        0-indexed dimension pair to plot. Default: ``(0, 1)``.
    method
        KDE method: ``"wkde"`` (custom weighted KDE) or ``"ks"``
        (scipy gaussian_kde with weights).
    adjust
        Bandwidth adjustment factor.
    size
        Point size for scatter plot.
    cmap
        Matplotlib colormap name. Default: ``"viridis"``.
    figsize
        Figure size ``(width, height)`` in inches.
    ncols
        Number of columns for multi-panel layout.
    show
        Show the plot. If ``None``, uses matplotlib's default.
    save
        Path to save the figure.
    ax
        Pre-existing axes (only for single-feature plots).

    Returns
    -------
    For a single feature with *ax* provided: the ``Axes``.
    For multiple features or when creating a new figure: the ``Figure``.
    Returns ``None`` when ``show=True``.

    Examples
    --------
    >>> import scanpy as sc
    >>> import pynebulosa as nb
    >>> adata = sc.datasets.pbmc3k_processed()
    >>> nb.plot_density(adata, "CD4")
    >>> nb.plot_density(adata, ["CD8A", "CCR7"], joint=True)
    """
    # Normalize features to list
    if isinstance(features, str):
        features = [features]
    if not features:
        raise ValueError("At least one feature must be specified")

    # Extract data (dims validated inside _get_embeddings)
    embeddings, dim_names = _get_embeddings(adata, reduction, dims)
    feature_data = _get_feature_data(adata, features, layer)

    # Calculate densities for each feature
    densities = {}
    for i, feat in enumerate(features):
        w = feature_data[:, i]
        densities[feat] = calculate_density(w, embeddings, method=method, adjust=adjust)

    panel_densities = [densities[f] for f in features]
    panel_labels = list(features)
    colorbar_labels = ["Density"] * len(features)

    if joint and len(features) > 1:
        joint_dens = np.prod(panel_densities, axis=0)
        joint_label = " ".join(f"{f}+" for f in features)
        panel_densities.append(joint_dens)
        panel_labels.append(joint_label)
        colorbar_labels.append("Joint density")

    n_panels = len(panel_labels)

    # Single feature with provided axes
    if n_panels == 1 and ax is not None:
        _density_scatter(
            ax,
            embeddings,
            panel_densities[0],
            panel_labels[0],
            dim_names,
            size,
            cmap,
            colorbar_labels[0],
        )
        return _finalize(ax.figure, show, save, return_axes=ax)

    # Create figure with subplots
    nrows = math.ceil(n_panels / ncols)
    actual_ncols = min(n_panels, ncols)
    if figsize is None:
        figsize = (5 * actual_ncols, 4.5 * nrows)

    fig, axes = plt.subplots(nrows, actual_ncols, figsize=figsize, squeeze=False)

    for i in range(n_panels):
        row, col = divmod(i, actual_ncols)
        _density_scatter(
            axes[row, col],
            embeddings,
            panel_densities[i],
            panel_labels[i],
            dim_names,
            size,
            cmap,
            colorbar_labels[i],
        )

    # Hide unused axes
    for i in range(n_panels, nrows * actual_ncols):
        row, col = divmod(i, actual_ncols)
        axes[row, col].set_visible(False)

    fig.tight_layout()

    return _finalize(fig, show, save)


def _is_inline_backend() -> bool:
    """Check if matplotlib is using a Jupyter inline backend."""
    return "inline" in matplotlib.get_backend()


def _finalize(
    fig_or_ax,
    show: bool | None,
    save: str | None,
    return_axes=None,
):
    """Handle show/save/return logic."""
    fig = fig_or_ax if isinstance(fig_or_ax, plt.Figure) else fig_or_ax.figure

    if save is not None:
        fig.savefig(save, dpi=300, bbox_inches="tight")

    if show is True:
        plt.show()
        return None

    if show is False:
        plt.close(fig)
        return return_axes if return_axes is not None else fig

    # show is None — auto-detect: in Jupyter inline backend, show and
    # return None to prevent the figure from being displayed twice
    # (once by the backend and once by Jupyter's repr of the return value).
    if _is_inline_backend():
        plt.show()
        return None

    return return_axes if return_axes is not None else fig
