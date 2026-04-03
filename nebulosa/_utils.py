"""Utility helpers for data extraction and validation.

Port of R/helpers.R and data extraction logic from R/methods.R.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from anndata import AnnData


def _validate_dims(dims: tuple[int, int]) -> None:
    """Ensure exactly 2 dimensions are specified.

    Parameters
    ----------
    dims
        Tuple of dimension indices.

    Raises
    ------
    ValueError
        If *dims* does not have exactly 2 elements.
    """
    if len(dims) != 2:
        raise ValueError(
            f"Exactly 2 dimensions are required, got {len(dims)}"
        )


def _get_embeddings(
    adata: AnnData,
    reduction: str | None,
    dims: tuple[int, int],
) -> tuple[np.ndarray, tuple[str, str]]:
    """Extract 2D cell embeddings from an AnnData object.

    Parameters
    ----------
    adata
        AnnData object with computed reductions in ``.obsm``.
    reduction
        Key in ``adata.obsm`` (e.g. ``"X_umap"``). If ``None``,
        auto-detects: tries ``"X_umap"``, then the last key.
    dims
        0-indexed dimension indices to extract.

    Returns
    -------
    Tuple of ``(embeddings, dim_labels)`` where *embeddings* is
    an ``(N, 2)`` array and *dim_labels* are display names.

    Raises
    ------
    KeyError
        If the specified or detected reduction is not found.
    ValueError
        If the requested dimensions are out of range.
    """
    _validate_dims(dims)

    if not adata.obsm:
        raise KeyError("No reductions found in adata.obsm")

    obsm_keys = list(adata.obsm.keys())

    if reduction is None:
        if "X_umap" in obsm_keys:
            reduction = "X_umap"
        else:
            reduction = obsm_keys[-1]

    if reduction not in adata.obsm:
        raise KeyError(
            f"Reduction '{reduction}' not found in adata.obsm. "
            f"Available: {obsm_keys}"
        )

    emb = np.asarray(adata.obsm[reduction])

    # Validate dimension indices
    for d in dims:
        if d < 0 or d >= emb.shape[1]:
            raise ValueError(
                f"Dimension {d} not present in '{reduction}' "
                f"(has {emb.shape[1]} dimensions)"
            )

    cell_embeddings = emb[:, list(dims)]

    # Create readable labels (e.g. "UMAP 1", "UMAP 2")
    base = reduction.replace("X_", "").upper()
    dim_labels = (f"{base} {dims[0] + 1}", f"{base} {dims[1] + 1}")

    return cell_embeddings, dim_labels


def _get_feature_data(
    adata: AnnData,
    features: list[str],
    layer: str | None,
) -> np.ndarray:
    """Extract feature data (gene expression or metadata) from AnnData.

    Checks ``adata.obs`` first (metadata), then gene expression
    from ``adata.X`` or ``adata.layers``.

    Parameters
    ----------
    adata
        AnnData object.
    features
        List of feature names (gene names or obs column names).
    layer
        Layer name for expression data. ``None`` uses ``adata.X``.

    Returns
    -------
    Array of shape ``(N, n_features)``.

    Raises
    ------
    ValueError
        If any feature is not found.
    """
    # Check if all features are in metadata
    if all(f in adata.obs.columns for f in features):
        return adata.obs[features].to_numpy(dtype=float)

    # Check gene expression
    missing = [f for f in features if f not in adata.var_names]
    if missing:
        # Check if missing are in obs
        obs_missing = [f for f in missing if f not in adata.obs.columns]
        if obs_missing:
            raise ValueError(
                f"Feature(s) not found in var_names or obs: "
                f"{', '.join(obs_missing)}"
            )

    # Get expression matrix
    if layer is None:
        mat = adata.X
    else:
        if layer not in adata.layers:
            raise ValueError(
                f"Layer '{layer}' not found. "
                f"Available: {list(adata.layers.keys())}"
            )
        mat = adata.layers[layer]

    # Extract columns for each feature
    var_index = {name: i for i, name in enumerate(adata.var_names)}
    result = np.empty((adata.n_obs, len(features)))
    for i, feat in enumerate(features):
        if feat in var_index:
            col = mat[:, var_index[feat]]
            if sp.issparse(col):
                col = col.toarray()
            result[:, i] = np.asarray(col).flatten()
        else:
            result[:, i] = adata.obs[feat].to_numpy(dtype=float)

    return result
