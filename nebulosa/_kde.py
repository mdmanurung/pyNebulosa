"""Core kernel density estimation algorithms.

Port of the R package Nebulosa's KDE functions (R/kde.R).
Implements weighted 2D kernel density estimation for single-cell
gene expression visualization.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.stats import gaussian_kde, norm


def _bandwidth(x: np.ndarray) -> float:
    """Compute bandwidth using Silverman's rule of thumb.

    Approximates R's ``ks::hpi()`` plug-in bandwidth selector.

    Parameters
    ----------
    x
        1D array of data values.

    Returns
    -------
    Scalar bandwidth estimate.
    """
    n = len(x)
    std = np.std(x, ddof=1)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    # Silverman's rule: 0.9 * min(std, IQR/1.34) * n^(-1/5)
    a = min(std, iqr / 1.34) if iqr > 0 else std
    return 0.9 * a * n ** (-0.2)


def wkde2d(
    x: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    h: tuple[float, float] | None = None,
    adjust: float = 1.0,
    n: int = 100,
    lims: tuple[float, float, float, float] | None = None,
) -> dict:
    """Weighted 2D kernel density estimation.

    Direct port of the R function ``Nebulosa:::wkde2d``.

    Parameters
    ----------
    x
        1D array of coordinates for dimension 1 (N cells).
    y
        1D array of coordinates for dimension 2 (N cells).
    w
        1D array of weights (e.g. gene expression values).
    h
        Bandwidths for x and y directions. If ``None``, computed
        automatically using Silverman's rule.
    adjust
        Bandwidth adjustment factor. Default: 1.0.
    n
        Number of grid points in each direction. Default: 100.
    lims
        Grid limits as ``(x_min, x_max, y_min, y_max)``.
        If ``None``, uses data range.

    Returns
    -------
    Dictionary with keys ``"x"``, ``"y"`` (grid coordinates, each
    length *n*) and ``"z"`` (density matrix, shape ``(n, n)``).
    """
    nx = len(x)
    if not (nx == len(y) == len(w)):
        raise ValueError("x, y, and w must have the same length")
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise ValueError("Missing or infinite values in coordinates are not allowed")

    if lims is None:
        lims = (x.min(), x.max(), y.min(), y.max())
    if not all(np.isfinite(v) for v in lims):
        raise ValueError("Only finite values are allowed in lims")

    # Bandwidth selection
    if h is None:
        h = (_bandwidth(x), _bandwidth(y))
    h = (h[0] * adjust, h[1] * adjust)

    # Create grid
    gx = np.linspace(lims[0], lims[1], n)
    gy = np.linspace(lims[2], lims[3], n)

    # Standardized distances: (n, N)
    ax = (gx[:, None] - x[None, :]) / h[0]
    ay = (gy[:, None] - y[None, :]) / h[1]

    # Weighted Gaussian kernel values: (n, N)
    wx = norm.pdf(ax) * w[None, :]
    wy = norm.pdf(ay) * w[None, :]

    # Density matrix: (n, n)
    w_total = w.sum() * h[0] * h[1]
    if w_total == 0:
        z = np.zeros((n, n))
    else:
        z = wx @ wy.T / w_total

    return {"x": gx, "y": gy, "z": np.asarray(z)}


def _get_dens(
    coords: np.ndarray,
    grid_x: np.ndarray,
    grid_y: np.ndarray,
    z: np.ndarray,
) -> np.ndarray:
    """Map grid-based density estimates back to individual data points.

    Port of R's ``get_dens`` using ``np.searchsorted`` (equivalent
    to R's ``findInterval``).

    Parameters
    ----------
    coords
        Cell coordinates, shape ``(N, 2)``.
    grid_x
        Grid x-coordinates from :func:`wkde2d`.
    grid_y
        Grid y-coordinates from :func:`wkde2d`.
    z
        Density matrix from :func:`wkde2d`, shape ``(n, n)``.

    Returns
    -------
    1D array of density values, one per cell.
    """
    ix = np.searchsorted(grid_x, coords[:, 0]).clip(0, len(grid_x) - 1)
    iy = np.searchsorted(grid_y, coords[:, 1]).clip(0, len(grid_y) - 1)
    return z[ix, iy]


def _get_dens_ks(
    coords: np.ndarray,
    kde_result: gaussian_kde,
) -> np.ndarray:
    """Map scipy gaussian_kde density estimates to individual data points.

    Parameters
    ----------
    coords
        Cell coordinates, shape ``(N, 2)``.
    kde_result
        Fitted ``scipy.stats.gaussian_kde`` object.

    Returns
    -------
    1D array of density values, one per cell.
    """
    return kde_result(coords.T)


def calculate_density(
    w: np.ndarray,
    coords: np.ndarray,
    method: str = "wkde",
    adjust: float = 1.0,
    map: bool = True,
) -> np.ndarray | dict:
    """Estimate weighted kernel density.

    Port of the R function ``Nebulosa:::calculate_density``.

    Parameters
    ----------
    w
        Weight vector (e.g. gene expression), length N.
    coords
        2D coordinates (e.g. UMAP embedding), shape ``(N, 2)``.
    method
        KDE method: ``"wkde"`` (custom weighted KDE) or ``"ks"``
        (scipy's ``gaussian_kde`` with weights).
    adjust
        Bandwidth adjustment factor (only for ``"wkde"``).
    map
        If ``True``, return per-cell density values. If ``False``,
        return the raw density estimation result.

    Returns
    -------
    If *map* is ``True``, a 1D array of density values (length N).
    If *map* is ``False``, a dict (for ``"wkde"``) or
    ``gaussian_kde`` object (for ``"ks"``).
    """
    w = np.asarray(w, dtype=float)
    coords = np.asarray(coords, dtype=float)

    if method not in ("wkde", "ks"):
        raise ValueError(f"Unknown method '{method}'. Use 'wkde' or 'ks'.")

    # Handle zero-expression case
    w_sum = w.sum()
    if w_sum == 0:
        warnings.warn(
            "All weights are zero; returning zero densities.",
            stacklevel=2,
        )
        if map:
            return np.zeros(len(w))
        return {"x": np.array([]), "y": np.array([]), "z": np.array([])}

    # Normalize weights (matching R: w / sum(w) * length(w))
    w_norm = w / w_sum * len(w)

    if method == "wkde":
        dens = wkde2d(
            x=coords[:, 0],
            y=coords[:, 1],
            w=w_norm,
            adjust=adjust,
        )
        if map:
            return _get_dens(coords, dens["x"], dens["y"], dens["z"])
        return dens

    else:  # method == "ks"
        dens = gaussian_kde(
            coords.T,
            weights=w_norm,
        )
        if map:
            return _get_dens_ks(coords, dens)
        return dens
